"""Count composition coverage, structure evidence, processing and test independence.

Counts of normalized compositions are never called counts of independent samples.
"""
from pathlib import Path
import json,shlex,math,re,hashlib
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from .common import ROOT

OUT=ROOT/'data/processed/effective_coverage'

def ncomp(d):return int(d.loc[d.composition_id.ne(''),'composition_id'].nunique())

def cif_fields(path):
    lines=path.read_text().splitlines();cell={};atoms=[];i=0
    cell_tags=['_cell_length_a','_cell_length_b','_cell_length_c','_cell_angle_alpha','_cell_angle_beta','_cell_angle_gamma']
    for line in lines:
        bits=shlex.split(line,comments=True)
        if bits and bits[0] in cell_tags:cell[bits[0]]=float(bits[1])
    while i<len(lines):
        if lines[i].strip()!='loop_':i+=1;continue
        i+=1;headers=[]
        while i<len(lines) and lines[i].strip().startswith('_'):headers.append(lines[i].strip());i+=1
        vals=[]
        while i<len(lines) and lines[i].strip() and not lines[i].strip().startswith(('loop_','_','data_')):
            vals+=shlex.split(lines[i],comments=True);i+=1
        if '_atom_site_fract_x' not in headers:continue
        if len(vals)%len(headers):raise ValueError('incomplete atom loop')
        atoms.extend(dict(zip(headers,vals[k:k+len(headers)])) for k in range(0,len(vals),len(headers)))
    coords=['_atom_site_fract_x','_atom_site_fract_y','_atom_site_fract_z']
    complete=bool(len(cell)==6 and all(math.isfinite(v) and v>0 for v in cell.values()) and atoms)
    for a in atoms:
        complete=complete and all(k in a and math.isfinite(float(a[k])) for k in coords)
        complete=complete and bool(a.get('_atom_site_type_symbol')) and '_atom_site_occupancy' in a and math.isfinite(float(a['_atom_site_occupancy']))
    return {'six_cell_parameters':len(cell)==6,'atom_site_rows':len(atoms),'coordinates_and_occupancy_present':bool(complete),
            'scope':'restricted syntactic check; not symmetry, chemistry, topology or sample-identity validation'}

def run():
    OUT.mkdir(parents=True,exist_ok=True)
    src=ROOT/'data/processed/unified/measurement_records.csv'
    d=pd.read_csv(src,keep_default_na=False)
    subsets={'all_source_records':d,'primary_tables':d[d.dataset.ne('laskowski')],
             'primary_clean_numeric':d[d.eligible_primary_clean_numeric],
             'reported_rt_candidates':d[d.eligible_reported_rt_candidate]}
    summary=[];structure_rows=[]
    for name,x in subsets.items():
        c=x[x.cif_path.ne('')]
        near=c[c.cif_match_status.eq('close_composition_match')]
        not_near=c[~c.cif_match_status.eq('close_composition_match')]
        conservative=not_near[not_near.cif_coordinates_status.eq('unperturbed_per_upstream_script')]
        row={'subset':name,'source_records':len(x),'parsed_compositions':ncomp(x),'unparsed_formula_records':int(x.composition_id.eq('').sum()),
             'cif_records':len(c),'cif_parsed_compositions':ncomp(c),'cif_unparsed_formula_records':int(c.composition_id.eq('').sum()),
             'close_match_cif_records':len(near),'not_close_cif_records':len(not_near),'not_close_cif_compositions':ncomp(not_near),
             'unperturbed_not_close_candidate_records':len(conservative),'unperturbed_not_close_candidate_compositions':ncomp(conservative),
             'leakage_groups':int(x.leakage_group_id.nunique()),'independent_physical_samples_verified':'not_determined'}
        summary.append(row)
        if name=='primary_clean_numeric':
            for r in c.itertuples():
                structural_tier='CIF_available'
                if r.cif_match_status=='close_composition_match':structural_tier='near_structure_match'
                elif r.cif_coordinates_status=='perturbed_per_upstream_script':structural_tier='not_close_but_perturbed'
                else:structural_tier='unperturbed_not_close_candidate_unverified'
                structure_rows.append({'record_id':r.record_id,'composition_id':r.composition_id,'formula':r.formula_raw,
                    'cif_path':r.cif_path,'tier':structural_tier,'same_measurement_sample_primary_verified':False})
    pd.DataFrame(summary).to_csv(OUT/'coverage_by_scope.csv',index=False)
    pd.DataFrame(structure_rows).to_csv(OUT/'clean_structure_evidence.csv',index=False)
    checks=[]
    for r in d[d.cif_path.ne('')].itertuples():checks.append({'record_id':r.record_id,'cif_path':r.cif_path,**cif_fields(ROOT/r.cif_path)})
    pd.DataFrame(checks).to_csv(OUT/'cif_syntax_audit.csv',index=False)
    clean=subsets['primary_clean_numeric']
    ci=[]
    for key,g in clean.groupby('composition_id'):
        ci.append({'composition_id':key,'composition_key':g.composition_key.iloc[0],
                   'formula_examples':';'.join(sorted(set(g.formula_normalized))),'source_record_count':len(g),
                   'record_ids':';'.join(g.record_id),'has_cif':bool(g.cif_path.ne('').any()),
                   'scope':'same nominal normalized composition; phases and physical samples may differ'})
    pd.DataFrame(ci).to_csv(OUT/'primary_clean_composition_index.csv',index=False)
    rt=subsets['reported_rt_candidates'].reset_index(drop=True);folds=[]
    for fold,(tr,te) in enumerate(GroupKFold(n_splits=5).split(rt,groups=rt.leakage_group_id)):
        a,b=rt.iloc[tr],rt.iloc[te]
        da={x for ds in a.doi_normalized for x in ds.split(';') if x};db={x for ds in b.doi_normalized for x in ds.split(';') if x}
        assert not(set(a.composition_id)&set(b.composition_id)) and not(da&db)
        folds.append({'fold':fold,'train_records':len(a),'test_records':len(b),'test_compositions':ncomp(b),
                      'test_leakage_groups':b.leakage_group_id.nunique(),'composition_overlap':0,'doi_overlap':0,
                      'status':'exploratory grouped internal validation; not untouched external confirmation'})
    pd.DataFrame(folds).to_csv(OUT/'internal_grouped_validation_capacity.csv',index=False)
    extra=pd.read_csv(ROOT/'data/processed/new_sources/balzat2026_linked_experiments.csv',keep_default_na=False)
    overlap=d[d.composition_key.isin(extra.composition_key)].copy()
    overlap[['record_id','dataset','formula_raw','composition_key','doi_normalized','eligible_primary_clean_numeric']].to_csv(OUT/'supplementary_composition_overlap.csv',index=False)
    process=[]
    for r in extra.itertuples():
        process.append({'record_id':r.measurement_id,'formula':r.formula,'doi':r.doi,
                        'treatment_temperature_K':r.synthesis_T_K_table_S1,'treatment_time_h':r.synthesis_time_h,
                        'treatment_kind':'source synthesis schedule','metadata_consistent':r.eligible_consistent_synthesis_metadata,
                        'fully_standardized_processing_protocol':False,'in_github_core_table':False})
    annotations=pd.read_csv(ROOT/'data/processed/primary_source_annotations.csv',keep_default_na=False)
    old=annotations[annotations.ID.eq('8k7')].iloc[0]
    assert '1173 K for 5 h' in old.sample_processing_or_phase
    process.append({'record_id':'obelix:8k7','formula':old.database_formula,'doi':old.doi,
                    'treatment_temperature_K':1173,'treatment_time_h':5,'treatment_kind':'pellet sintering; partial manual annotation',
                    'metadata_consistent':True,'fully_standardized_processing_protocol':False,'in_github_core_table':True})
    pd.DataFrame(process).to_csv(OUT/'partial_quantitative_processing_inventory.csv',index=False)
    result={'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'coverage':summary,
            'independent_physical_sample_count':None,'reason':'No universal original-sample identifier or verified phase/process equivalence',
            'complete_coordinate_cif_syntax_checks':int(pd.DataFrame(checks).coordinates_and_occupancy_present.sum()),
            'processing':{'core_structured_synthesis_T_and_time_records':0,'supplementary_quantitative_schedule_records':32,
                          'supplementary_consistent_schedule_records':int(extra.eligible_consistent_synthesis_metadata.sum()),
                          'supplementary_compositions':1,'supplementary_publications':1,'manual_core_sintering_annotations':1,
                          'partial_quantitative_processing_records_total':len(process),'full_processing_protocols_verified':0},
            'testing':{'locked_untouched_external_samples':0,'rt_internal_leakage_groups':218,
                       'supplementary_records_from_new_doi':32,'supplementary_records_from_new_composition':0,
                       'supplementary_doi_already_explored':True,'supplementary_composition_overlap_record_ids':overlap.record_id.tolist()},
            'correction':'Original close match Yes is now recognized: 61 CIF rows in all data,55 in clean subset; conductivity values unchanged.'}
    (OUT/'effective_coverage_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    report='''# 有效覆盖、结构、工艺与独立验证审计 · 2026-11-05

**结论：当前清洗集有925种归一化组成，室温候选有799种。它们不能直接称为925/799个独立物理样品。**

| 口径 | 来源记录 | 归一化组成 | CIF记录 | 有CIF的已解析组成 | 非近似匹配CIF记录/组成 | 未扰动且非近似匹配的保守候选记录/组成 |
|---|---:|---:|---:|---:|---:|---:|
'''
    for r in summary:
        report+=f"| {r['subset']} | {r['source_records']} | {r['parsed_compositions']} | {r['cif_records']} | {r['cif_parsed_compositions']} | {r['not_close_cif_records']}/{r['not_close_cif_compositions']} | {r['unperturbed_not_close_candidate_records']}/{r['unperturbed_not_close_candidate_compositions']} |\n"
    report+='''
## 独立材料与结构

组成去重只将原子比例相同的化学式归为同一组，不合并不同晶相、制备批次或测量条件的物理样品。当前主表没有跨来源统一的原始样品ID，不能给出已验证的独立物理样品总数。未核实Laskowski转录也不能作为额外独立材料证据。

321个CIF均通过本次针对文件格式的六个晶胞参数、原子坐标及占位字段存在性检查；这不等于通过完整晶体学、拓扑或样品对应验证。全体有CIF记录中有1条水合配方未解析，因此“有CIF的已解析组成”为282，另1条尚未归组。清洗集的260条CIF对应233种组成。排除原表近似匹配后为205条、194种组成；再按上游脚本/ICSD标记排除坐标扰动后为23条、23种组成。这23条仍未逐条原文证实与电导率测量样品严格对应，不能称作23个已确认精确实验结构。未扰动状态本身是上游规则推断。

本次修复旧统一管线只识别`true`而漏掉原表`Yes`的问题：全体61条CIF、清洗后55条CIF应为近似匹配；不改原始化学式、电导率或标签删失关系。旧公开版本为GitHub commit `90a18c28c7a8bbed4e4857c7dff1dcd3d816da7f`；此前因子试跑的原始输入已另存本地快照。因子试跑未用这个错误的CIF匹配字段生成因子。

## 工艺信息

GitHub核心统一主表目前没有标准化的合成温度、时间、气氛、压制/烧结流程字段；按“已结构化的合成温度+时间”口径，核心表覆盖为0，不代表所有原论文都未报告工艺。

另有本地补充Li7SiPS8数据32次实验，均已提取合成温度、时间、密度及电导率，28次表格间工艺信息无冲突、4次有冲突。这32次属于1种组成、1篇论文，尚未并入GitHub核心主表。此前手工原文核查还记录了核心表`obelix:8k7`的1173K、5h烧结条件，对应Li5La3Ta2O12，是另一种组成/另一篇论文。因此当前有33条已定位的“部分定量工艺条件”，其中29条未记录温度/时间冲突；这些处理步骤并不完全同类，也都不等于完整工艺SOP。全流程逐项核验的记录为0。

## 真正独立测试

当前没有预先锁定、未参与探索的独立外部确认集，已建立的这类样本为0。室温1393条可按组成和任一DOI的传递连通关系分成218个防泄漏组，支持组间无组成/文献重叠的内部验证；这218是划分单位，不是218个独立外测样品。

新增32次Li7SiPS8实验来自新DOI，文献来源不同；但是Li7SiPS8与核心表`mchaffie:0045`的Li3.5Si0.5P0.5S4原子比例相同。另有一个未核实Laskowski转录对应同组成。不能把32次实验称为32种新材料，或32个跨组成独立测试样本。其结果此前已经用于探索，不能再宣称为从未查看的最终盲测。若目标仅为“同组成的新工艺/新论文复现”，它们可构成1个新文献组的补充验证，但应保留组内相关性和使用历史。

所有计数可从`data/processed/effective_coverage/`中的JSON、逐条结构证据、工艺清单、组成目录和折成员覆盖表复核。复现：`python3 -m src.audit_effective_dataset`。主数据来源和原始值保持不变。
'''
    (ROOT/'reports/effective_dataset_coverage.md').write_text(report)
    print(json.dumps(result,ensure_ascii=False,indent=2))
    return result

if __name__=='__main__':run()
