"""Lossless, offline integration of the requested experimental data sources.

Never equate a common composition/publication with a confirmed common sample.
Laskowski's damaged secondary transcription is provisional, not primary-verified.
"""
from pathlib import Path
from io import StringIO
from collections import defaultdict
import sys, json, re, hashlib, math, unicodedata, zipfile
from functools import lru_cache
import numpy as np
import pandas as pd
from .common import ROOT, RAW, SHA, composition, composition_key, chemistry
from .quality_rules import warning_reasons

DEST=ROOT/'data/processed/unified'
SOURCES=RAW/'unified_sources'
UPSTREAM=RAW/f'upstream/OBELiX-{SHA}'
ELEMENTS=set('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og'.split())

@lru_cache(maxsize=64)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def short_hash(value):return hashlib.sha256(value.encode()).hexdigest()[:20]
def number(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else np.nan
    except (ValueError,TypeError):return np.nan
def clean_text(v):return '' if v is None or str(v).lower()=='nan' else str(v).strip()

def parse_sigma(v):
    raw=clean_text(v)
    normalized=unicodedata.normalize('NFKC',raw).replace('−','-').replace('–','-')
    normalized=''.join(c for c in normalized if not unicodedata.combining(c))
    normalized=re.sub(r'\s+','',normalized).replace('≤','<=').replace('≥','>=')
    m=re.fullmatch(r'([<>]=?|~|≈)?([+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)',normalized)
    if not m:return (np.nan,'unparsed' if raw else 'missing',np.nan,normalized)
    x=float(m[2]);rel=m[1] or '='
    if x<=0 or not math.isfinite(x):return(np.nan,'invalid_nonpositive',np.nan,normalized)
    return x,rel,math.log10(x) if rel=='=' else np.nan,normalized

def formula_syntax(raw):
    s=unicodedata.normalize('NFKC',clean_text(raw))
    s=s.translate(str.maketrans('₀₁₂₃₄₅₆₇₈₉','0123456789'))
    if len(s)>=2 and s[0]==s[-1] and s[0] in ['"',"'"]:s=s[1:-1].strip()
    return s

def parse_formula(raw):
    s=formula_syntax(raw)
    # Whitespace may encode detached PDF subscripts: do not silently move digits.
    if re.search(r'\d\s+\d|[A-Za-z]\s+[A-Z]',s):return '',{},'ambiguous_detached_subscripts'
    s=''.join(s.split())
    if not s:return '',{},'missing'
    if any(t in s.lower() for t in ['wt%','mol%','glass','ceramic']):return '',{},'mixture_or_additive'
    if '*' in s or '·' in s:return '',{},'adduct_requires_review'
    if '-' in s or '+' in s:return '',{},'mixture_phase_or_symbolic_formula'
    try:
        if s.count('(')!=s.count(')') or s.count('[')!=s.count(']'):raise ValueError('unbalanced')
        c=composition(s)
        if set(c)-ELEMENTS:raise ValueError('unknown element')
        if any(n<0 or not math.isfinite(n) for n in c.values()) or sum(c.values())<=0:raise ValueError('bad count')
        c={k:v for k,v in c.items() if v>0}
        if 'Li' not in c:raise ValueError('no Li')
        return composition_key(s),c,'parsed'
    except (ValueError,IndexError,ZeroDivisionError):return '',{},'unparsed'

def dois(v):
    s=unicodedata.normalize('NFKC',clean_text(v)).lower()
    # Whitespace/semicolon separate source fields; retain multiple DOIs individually.
    return sorted(set(m.rstrip('.,;') for m in re.findall(r'10\.\d{4,9}/[^\s;]+',s)))

def base(source,sid,row,path,payload):
    return {'record_id':source+':'+str(sid),'dataset':source,'source_record_id':str(sid),
            'source_row':row,'source_file':str(path.relative_to(ROOT)), 'source_sha256':digest(path),
            'source_payload_json':json.dumps(payload,ensure_ascii=False,sort_keys=True),
            'formula_raw':'','doi_raw':'','sigma_raw':'','temperature_raw':'','temperature_unit_raw':'',
            'measurement_temperature_K':np.nan,'temperature_status':'unknown','reported_rt':False,
            'conductivity_type':'unknown','measurement_method_status':'direct_or_extrapolated_unknown',
            'lowest_extrapolation_temperature_K':np.nan,'structural_family_raw':'','space_group_raw':'',
            'space_group_number':np.nan,'icsd_id':'','cif_path':'','cif_match_status':'not_available',
            'cif_coordinates_status':'not_available','activation_energy_eV':np.nan,
            'bulk_sigma_raw':'','total_sigma_raw':'','source_note':'','quality_flags_initial':'',
            'provenance_status':'primary_public_table','original_obelix_split':'',
            'lattice_a_A':np.nan,'lattice_b_A':np.nan,'lattice_c_A':np.nan,
            'lattice_alpha_deg':np.nan,'lattice_beta_deg':np.nan,'lattice_gamma_deg':np.nan}

def load_obelix():
    path=RAW/'obelix/all.csv';d=pd.read_csv(path,keep_default_na=False,dtype=str)
    train=set(pd.read_csv(RAW/'obelix/train.csv').ID.astype(str))
    records=[];cifs=[];aux=[]
    archive=RAW/'obelix/all_cifs.zip'
    with zipfile.ZipFile(archive) as z:
        names={Path(n).stem:n for n in z.namelist() if n.endswith('.cif')}
        for i,v in d.iterrows():
            r=base('obelix',v['ID'],i+2,path,v.to_dict())
            r.update(formula_raw=v['Reduced Composition'],doi_raw=v['DOI'],sigma_raw=v['Ionic conductivity (S cm-1)'],
                     temperature_status='reported_room_temperature_exact_T_unknown',reported_rt=True,
                     structural_family_raw=v['Family'],space_group_raw=v['Space group'],
                     space_group_number=number(v['Space group #']),icsd_id=v['ICSD ID'],source_note=v['note'],
                     quality_flags_initial=warning_reasons(v['note']),bulk_sigma_raw=v['IC (Bulk)'],total_sigma_raw=v['IC (Total)'],
                     original_obelix_split='train' if v['ID'] in train else 'test')
            for out,key in [('lattice_a_A','a'),('lattice_b_A','b'),('lattice_c_A','c'),('lattice_alpha_deg','alpha'),('lattice_beta_deg','beta'),('lattice_gamma_deg','gamma')]:r[out]=number(v[key])
            for col,kind in [('IC (Total)','total'),('IC (Bulk)','bulk')]:
                if v[col]:
                    val,rel,log,_=parse_sigma(v[col]);main,mrel,_,_=parse_sigma(r['sigma_raw'])
                    if rel==mrel and math.isfinite(val) and math.isclose(val,main,rel_tol=1e-6):r['conductivity_type']=kind
                    aux.append({'parent_record_id':r['record_id'],'kind':kind,'sigma_raw':v[col],'sigma_value_S_cm':val,
                                'sigma_relation':rel,'log10_sigma_S_cm':log,'source_field':col,'source_file':r['source_file']})
            if v['ID'] in names:
                n=names[v['ID']];raw=z.read(n)
                target=DEST/'structures'/f"{v['ID']}.cif";target.parent.mkdir(exist_ok=True)
                target.write_bytes(raw)
                r.update(cif_path=str(target.relative_to(ROOT)),cif_match_status='close_composition_match' if v['close match'].lower()=='true' else 'source_matched_not_independently_verified',
                         cif_coordinates_status='perturbed_per_upstream_script' if v['ICSD ID'] else 'unperturbed_per_upstream_script')
                cifs.append({'record_id':r['record_id'],'cif_path':r['cif_path'],'archive_member':n,
                             'cif_sha256':hashlib.sha256(raw).hexdigest(),'archive_sha256':digest(archive),
                             'coordinate_status':r['cif_coordinates_status'],'match_status':r['cif_match_status'],
                             'license':'CC-BY-4.0; released CIFs; no new ICSD access granted'})
            records.append(r)
    return records,cifs,aux

def load_liion():
    # The full live official HTML table includes the 2024 correction to record 703.
    path=SOURCES/'liion_original.html'
    d=pd.read_html(StringIO(path.read_text()),keep_default_na=False)[0]
    d=d[pd.to_numeric(d[0],errors='coerce').notna()].copy()
    d.columns=['ID','composition','source','temperature','target','log_target','family','ChemicalFamily']
    assert len(d)==820 and d.ID.nunique()==820
    d.to_csv(DEST/'liion_official_table.csv',index=False)
    records=[]
    for ordinal,(_,v) in enumerate(d.iterrows(),1):
        v=v.map(clean_text);r=base('liion',v.ID,ordinal,path,v.to_dict())
        t=number(v.temperature)
        r.update(formula_raw=v.composition,doi_raw=v.source,sigma_raw=v.target,temperature_raw=v.temperature,
                 temperature_unit_raw='C',measurement_temperature_K=t+273.15,
                 temperature_status='numeric_temperature_in_official_table',reported_rt=bool(15<=t<=35),
                 structural_family_raw=v.family,source_note='Official source table ordinal; not an HTML line number.')
        if v.ID=='703':r['source_note']+=' Official 2024-10-18 correction used: Li6NBr3, 3.10e-7 S/cm, DOI 10.1002/anie.198000741.'
        records.append(r)
    old=pd.read_csv(UPSTREAM/'data/misc/LiIonDatabase.csv',dtype=str,keep_default_na=False)
    changes=[]
    # The old bundled CSV has two ID922 entries; the live table corrects one to929.
    # Compare table ordinals rather than treating that corrupted ID as unique.
    assert len(old)==len(d)
    for ordinal,(_,v) in enumerate(d.iterrows()):
        for field in ['ID','composition','source','temperature','target','family','ChemicalFamily']:
            a=clean_text(old.iloc[ordinal][field]);b=clean_text(v[field])
            missing={'','#N/A','nan','N/A'}
            same=(a==b) or (a in missing and b in missing) or (field in ['target','temperature'] and math.isclose(number(a),number(b),rel_tol=1e-7,abs_tol=1e-13))
            if not same:changes.append({'source_record_id':v.ID,'table_ordinal':ordinal+1,'field':field,'obelix_bundled_value':a,'official_live_value':b})
    pd.DataFrame(changes).to_csv(DEST/'liion_version_differences.csv',index=False)
    return records

def load_caltech():
    path=RAW/'mchaffie_conductivity.csv';d=pd.read_csv(path,dtype=str,keep_default_na=False);records=[]
    for i,v in d.iterrows():
        r=base('mchaffie',f'{i+1:04d}',i+2,path,v.to_dict());t=number(v.lowest_extrapolation_temperature_K)
        r.update(formula_raw=v.compound,doi_raw=v.conductivity_doi,sigma_raw=v.conductivity_siemens_per_cm,
                 icsd_id=v.icsd_collectioncode,lowest_extrapolation_temperature_K=t,reported_rt=True,
                 temperature_status='reported_rt_target_exact_T_unknown',
                 measurement_method_status='rt_extrapolated_from_higher_T' if math.isfinite(t) and t>308.15 else 'direct_or_extrapolated_unknown')
        records.append(r)
    return records

def load_laskowski_transcription():
    """One candidate per nonempty conductivity cell; all remain secondary/unverified.

    Only join adjacent formula fragments when the conductivity row has no formula
    and both neighbors have no conductivity. Ambiguities remain quarantined.
    """
    path=UPSTREAM/'data/misc/laskowski_semi-fromatted.csv'
    d=pd.read_csv(path,header=None,keep_default_na=False,dtype=str)
    lookup=pd.read_csv(UPSTREAM/'data/misc/doi_lookup_table.csv',dtype=str).set_index('ref. number').DOI.to_dict()
    records=[];transform=[]
    for i,v in d[d[1].ne('')].iterrows():
        r=base('laskowski',f'csvline{i+1:04d}',i+1,path,dict(zip([str(n) for n in d.columns],v)))
        indices=[i];formula=v[0];flags=['secondary_transcription_not_primary_verified']
        # Do not infer a complete formula from one side of a detached label row.
        if not formula:
            neighbors=[j for j in [i-1,i+1] if 0<=j<len(d) and d.at[j,1]=='']
            indices=sorted([i]+neighbors)
            fragments=[d.at[j,0] for j in indices if d.at[j,0]]
            formula=''.join(fragments)
            flags.append('adjacent_fragments_joined_requires_primary_check')
        cells=d.loc[indices]
        def cell(col):
            vals=list(dict.fromkeys(x for x in cells[col] if x))
            if len(vals)>1:flags.append(f'conflicting_transcription_column_{col}')
            return vals[0] if len(vals)==1 else ''
        ref=cell(7)
        r.update(formula_raw=formula,sigma_raw=v[1],doi_raw=lookup.get(ref,''),reported_rt=True,
                 temperature_status='reported_rt_secondary_transcription',provenance_status='secondary_damaged_transcription_unverified',
                 structural_family_raw=cell(5),space_group_raw=cell(3),space_group_number=number(cell(4)),
                 icsd_id=cell(6),activation_energy_eV=number(cell(2).replace('\u0305','')),
                 source_note=f'Laskowski reference number: {ref}; transcription CSV lines: '+','.join(str(j+1) for j in indices),
                 quality_flags_initial=';'.join(flags))
        r['source_payload_json']=json.dumps({'transcription_rows':[{'csv_line':j+1,'cells':d.loc[j].tolist()} for j in indices],
                                             'reference_number':ref,'doi_lookup_file':str((UPSTREAM/'data/misc/doi_lookup_table.csv').relative_to(ROOT))},ensure_ascii=False)
        if len(indices)>1:transform.append({'record_id':r['record_id'],'input_lines':','.join(str(j+1) for j in indices),
                                          'operation':'concatenate_adjacent_formula_fragments','output_formula':formula,'status':'unverified'})
        records.append(r)
    pd.DataFrame(transform).to_csv(DEST/'laskowski_fragment_recovery_log.csv',index=False)
    # Preserve every damaged raw line, including fragments not confidently assigned.
    q=d.copy();q.columns=['formula_fragment','sigma_fragment','Ea_fragment','SG_fragment','SG_number_fragment','phase_fragment','ICSD_fragment','reference_fragment']
    q.insert(0,'csv_line',np.arange(1,len(q)+1));q.to_csv(DEST/'laskowski_all_transcription_fragments.csv',index=False)
    return records

def normalize(records):
    for r in records:
        key,c,status=parse_formula(r['formula_raw']);val,rel,log,norm=parse_sigma(r['sigma_raw'])
        flags=[x for x in r.pop('quality_flags_initial').split(';') if x]
        if status!='parsed':flags.append('formula_'+status)
        if rel in ['unparsed','missing','invalid_nonpositive']:flags.append('conductivity_'+rel)
        ds=dois(r['doi_raw'])
        if not ds:flags.append('doi_missing_or_unparsed')
        # A preprint-based mapping concern, not a confirmed correction.
        if '10.1021/acs.chemmater.0c04150' in ds and set(c)=={'Li','Sn','S','Se'}:flags.append('Sn_S_Se_mapping_requires_final_paper_review')
        formula_normalized=''.join(formula_syntax(r['formula_raw']).split()) if status=='parsed' else ''
        r.update(formula_normalized=formula_normalized,
                 formula_parse_status=status,composition_key=key,composition_id='cmp_'+short_hash(key) if key else '',
                 formula_normalization_actions='strip_matching_outer_quotes;unicode_nfkc;remove_whitespace' if formula_normalized!=clean_text(r['formula_raw']) else 'none',
                 elements=';'.join(sorted(c)),n_elements=len(c) if c else np.nan,
                 stoichiometry_json=json.dumps(c,sort_keys=True) if c else '',
                 chemistry_class=chemistry(formula_normalized) if key else '',
                 doi_normalized=';'.join(ds),sigma_value_S_cm=val,sigma_relation=rel,
                 sigma_unit='S/cm',sigma_text_normalized=norm,log10_sigma_S_cm=log,
                 quality_flags=';'.join(sorted(set(flags))))
        r['publication_composition_group']='pc_'+short_hash(key+'|'+ds[0]) if key and ds else ''
        r['license_status']={'obelix':'CC-BY-4.0','liion':'University of Liverpool academic attribution; commercial licence required',
                             'mchaffie':'CC0-1.0','laskowski':'OBELiX CC-BY-4.0 transcription; original RSC SI rights not established'}[r['dataset']]
        r['eligible_numeric_label']=bool(key and rel=='=' and math.isfinite(log))
        r['eligible_primary_clean_numeric']=bool(r['eligible_numeric_label'] and not flags and r['dataset']!='laskowski')
        r['eligible_reported_rt_candidate']=bool(r['eligible_primary_clean_numeric'] and r['reported_rt'] and r['measurement_method_status']!='rt_extrapolated_from_higher_T')
    return pd.DataFrame(records)

def link_records(d):
    """Pairwise links are evidence tiers, never destructive sample deduplication."""
    groups=defaultdict(list)
    for idx,r in d.iterrows():
        if r.composition_key:
            for doi in r.doi_normalized.split(';'):
                if doi:groups[(r.composition_key,doi)].append(idx)
    links=[];seen=set()
    for (comp,doi),indices in groups.items():
        for a_pos,a in enumerate(indices):
            for b in indices[a_pos+1:]:
                x,y=d.loc[a],d.loc[b]
                if x.dataset==y.dataset:continue
                k=(x.record_id,y.record_id,doi)
                if k in seen:continue
                seen.add(k)
                t1,t2=x.measurement_temperature_K,y.measurement_temperature_K
                t_known=bool(pd.notna(t1) and pd.notna(t2));t_match=bool(t_known and abs(t1-t2)<=0.5)
                t_conflict=bool(t_known and abs(t1-t2)>0.5)
                sg_conflict=bool(pd.notna(x.space_group_number) and pd.notna(y.space_group_number) and x.space_group_number!=y.space_group_number)
                type_conflict=x.conductivity_type!='unknown' and y.conductivity_type!='unknown' and x.conductivity_type!=y.conductivity_type
                v1,v2=x.sigma_value_S_cm,y.sigma_value_S_cm
                delta=abs(math.log10(v1/v2)) if pd.notna(v1) and pd.notna(v2) else np.nan
                agree=bool(x.sigma_relation==y.sigma_relation and pd.notna(delta) and delta<=0.01)
                tier='same_composition_publication_only'
                if t_conflict or sg_conflict or type_conflict:tier='condition_or_phase_conflict'
                elif agree and t_match and x.conductivity_type==y.conductivity_type!='unknown':tier='strong_metadata_agreement_sample_identity_unverified'
                elif agree:tier='label_agreement_conditions_incomplete'
                elif pd.notna(delta) and delta>1:tier='large_label_discrepancy_needs_review'
                links.append({'record_id_a':x.record_id,'record_id_b':y.record_id,'dataset_a':x.dataset,'dataset_b':y.dataset,
                              'composition_id':x.composition_id,'shared_doi':doi,'link_tier':tier,'abs_delta_log10_sigma':delta,
                              'same_relation_and_label_within_0.01dex':agree,'both_numeric_T':t_known,
                              'temperature_conflict':t_conflict,'space_group_conflict':sg_conflict,'conductivity_type_conflict':type_conflict,
                              'same_physical_sample_confirmed':False})
    cols=['record_id_a','record_id_b','dataset_a','dataset_b','composition_id','shared_doi','link_tier','abs_delta_log10_sigma','same_relation_and_label_within_0.01dex','both_numeric_T','temperature_conflict','space_group_conflict','conductivity_type_conflict','same_physical_sample_confirmed']
    out=pd.DataFrame(links,columns=cols)
    out.to_csv(DEST/'cross_source_links.csv',index=False)
    out[out.link_tier.isin(['condition_or_phase_conflict','large_label_discrepancy_needs_review'])].to_csv(DEST/'cross_source_conflicts.csv',index=False)
    return out

def split_groups(d):
    # Union all shared compositions and all individual DOIs to prevent either form
    # of overlap in downstream group splits. This is not sample deduplication.
    parent=list(range(len(d)))
    def find(i):
        while i!=parent[i]:parent[i]=parent[parent[i]];i=parent[i]
        return i
    seen={}
    for i,r in d.iterrows():
        keys=(['c:'+r.composition_key] if r.composition_key else [])+['d:'+x for x in r.doi_normalized.split(';') if x]
        for key in keys:
            if key in seen:parent[find(i)]=find(seen[key])
            else:seen[key]=i
    members=defaultdict(list)
    for i in range(len(d)):members[find(i)].append(d.at[i,'record_id'])
    ids={k:'split_'+short_hash('|'.join(sorted(v))) for k,v in members.items()}
    d['leakage_group_id']=[ids[find(i)] for i in range(len(d))]

def write_outputs(d,links,cifs,aux):
    split_groups(d)
    d.to_csv(DEST/'measurement_records.csv',index=False)
    # Conservative cleaned view; ambiguous source transcriptions stay in master.
    clean=d[d.eligible_primary_clean_numeric].copy()
    clean.to_csv(DEST/'clean_numeric_records.csv',index=False)
    rt=d[d.eligible_reported_rt_candidate].copy();rt.to_csv(DEST/'reported_rt_candidates.csv',index=False)
    d[d.quality_flags.ne('')].to_csv(DEST/'review_queue.csv',index=False)
    d[d.sigma_relation.isin(['<','<=','>','>='])].to_csv(DEST/'censored_records.csv',index=False)
    pd.DataFrame(cifs).to_csv(DEST/'structure_manifest.csv',index=False)
    pd.DataFrame(aux).to_csv(DEST/'obelix_auxiliary_bulk_total.csv',index=False)
    cleaning=[]
    for r in d.itertuples():
        for field,before,after,rule in [
            ('formula',r.formula_raw,r.formula_normalized,'unicode/outer_quotes/whitespace; no inferred stoichiometry'),
            ('conductivity',r.sigma_raw,str(r.sigma_value_S_cm),'parse S/cm value and preserve relation separately'),
            ('DOI',r.doi_raw,r.doi_normalized,'extract and lowercase individual DOI strings')]:
            if before!=after:cleaning.append({'record_id':r.record_id,'field':field,'input':before,'output':after,'rule':rule})
        if r.dataset=='liion':cleaning.append({'record_id':r.record_id,'field':'temperature','input':r.temperature_raw+' C','output':str(r.measurement_temperature_K)+' K','rule':'add273.15'})
        if r.sigma_relation!='=':cleaning.append({'record_id':r.record_id,'field':'regression_label','input':r.sigma_raw,'output':'null','rule':'nonexact, censored or invalid label excluded from log target'})
        if r.quality_flags:cleaning.append({'record_id':r.record_id,'field':'quality_flags','input':r.source_note,'output':r.quality_flags,'rule':'retain record in master and exclude from primary clean view'})
    pd.DataFrame(cleaning).to_csv(DEST/'cleaning_log.csv',index=False)
    comps=[]
    for key,g in d[d.composition_id.ne('')].groupby('composition_id',sort=True):
        comps.append({'composition_id':key,'composition_key':g.composition_key.iloc[0],
                      'formula_examples':';'.join(sorted(set(g.formula_normalized))), 'n_source_records':len(g),
                      'datasets':';'.join(sorted(set(g.dataset))),'n_publication_dois':len(set(x for ds in g.doi_normalized for x in ds.split(';') if x)),
                      'record_ids':';'.join(g.record_id),'scope':'composition only; may include multiple phases/samples/temperatures'})
    pd.DataFrame(comps).to_csv(DEST/'composition_index.csv',index=False)
    coverage=[]
    for name,g in d.groupby('dataset',sort=False):
        coverage.append({'dataset':name,'imported_source_records':len(g),'parsed_formulas':int(g.composition_key.ne('').sum()),
                         'exact_numeric_labels':int(g.sigma_relation.eq('=').sum()),'numeric_temperature':int(g.measurement_temperature_K.notna().sum()),
                         'clean_primary_numeric':int(g.eligible_primary_clean_numeric.sum()),'reported_rt_candidates':int(g.eligible_reported_rt_candidate.sum()),
                         'cif_files':int(g.cif_path.ne('').sum()),'import_status':'provisional_secondary_transcription' if name=='laskowski' else 'complete_public_table'})
    coverage.append({'dataset':'ddse','imported_source_records':0,'import_status':'bulk_export_not_public; collaborator_access_required'})
    pd.DataFrame(coverage).to_csv(DEST/'source_coverage.csv',index=False)
    descriptions={
        'record_id':'稳定来源记录标识；不是独立样品标识。', 'source_record_id':'原数据 ID；Laskowski 为原转录 CSV 电导率所在行。',
        'dataset':'来源：obelix/liion/mchaffie/laskowski；DDSE 无整库数据，不制造占位测量行。',
        'source_file':'项目相对路径，指向可追溯的原始文件。', 'source_sha256':'来源文件字节内容的 SHA256。',
        'source_row':'CSV 行号（包含表头）；LiIon 为官网表内数据行序号；Laskowski 无表头。',
        'source_payload_json':'完整来源字段或转录碎片，未自动覆盖原数值。',
        'formula_raw':'来源化学式或明确标记的邻行拼接结果。', 'formula_normalized':'规范 Unicode、成对外层引号和空白；不猜补元素、下标、相比例。',
        'formula_parse_status':'parsed 或明确的拒绝解析原因；parsed 仅证明语法通过，不证明实验化学式正确。',
        'formula_normalization_actions':'仅记录语法规范化操作；原式保留在 formula_raw。',
        'doi_raw':'原来源的文献标识字符串。','doi_normalized':'提取并小写的 DOI，多个用分号分隔；未在线核验每个 DOI 是否有效。',
        'sigma_raw':'来源原始电导率文本，单位由来源表指定为 S/cm。','sigma_unit':'本次四个来源统一为 S/cm；没有把 mS/cm 当成 S/cm。',
        'sigma_text_normalized':'仅去除排版组合字符、空白并规范比较符号后的文本。',
        'temperature_raw':'来源明确给出的温度文本；无数值不填。','temperature_unit_raw':'来源温度单位；LiIon 为摄氏度 C。',
        'temperature_status':'区分明确数值温度、来源室温声明和二手转录的室温声明。',
        'composition_key':'元素原子分数，10 位小数；仅作组成匹配，不确认晶体或样品相同。',
        'composition_id':'组成键的 SHA256 前20位；不是结构 ID。','publication_composition_group':'组成+首个 DOI 的检索组；多 DOI 的完整配对见 links。',
        'sigma_value_S_cm':'S/cm；删失记录为界限值，不能直接作回归标签。','sigma_relation':'=、<、<=、>、>=、~、≈；保留删失/近似关系。',
        'log10_sigma_S_cm':'仅精确正数标签取 log10；删失/近似/解析失败为空。',
        'measurement_temperature_K':'仅来源明确数值温度时填入；LiIon 摄氏度+273.15。未把 RT 强填298.15K。',
        'reported_rt':'来源声称 RT，或 LiIon 数值温度15–35°C的操作性窗口。',
        'structural_family_raw':'来源结构家族/相标记；不是统一验证的晶体结构。',
        'space_group_raw':'来源空间群原始文本；PDF转录残损未猜补。','space_group_number':'来源空间群编号，缺失为空。',
        'icsd_id':'ICSD 标识而非结构文件；不据此宣称已取得 CIF 或许可。',
        'cif_path':'仅实际随本数据交付的 CIF 路径。','cif_match_status':'来源提供的结构对应程度，含 close match；未独立核验全部结构。',
        'activation_energy_eV':'仅 Laskowski 二手转录给出的活化能，单位 eV；未核实值不进主清洗子集。',
        'bulk_sigma_raw':'OBELiX 原始 bulk 电导率字段，S/cm；另表保存为辅助标签。',
        'total_sigma_raw':'OBELiX 原始 total 电导率字段，S/cm；另表保存为辅助标签。',
        'source_note':'原注释、来源更正说明或转录定位说明。',
        'original_obelix_split':'OBELiX 原始 train/test；不能直接沿用为合并后独立验证划分。',
        'conductivity_type':'bulk/total 只在对应原字段与主标签一致时填写；否则 unknown。',
        'measurement_method_status':'是否外推；McHaffie 有高温最低测量标记则标为 RT 外推；na 不代表直接测量。',
        'lowest_extrapolation_temperature_K':'McHaffie 原始最低测量/外推温度；不是该 RT 标签的实际测量温度。',
        'provenance_status':'primary_public_table 或 secondary_damaged_transcription_unverified。',
        'quality_flags':'所有质量问题；不因数值看起来合理就补正。',
        'eligible_primary_clean_numeric':'可解析锂组成、精确正数标签、无质量标记、非未核实二手转录。不是已验证同质训练集。',
        'eligible_reported_rt_candidate':'上述条件+来源RT/15–35°C且未明确高温外推；未知bulk/total或直接性仍可能存在。',
        'leakage_group_id':'共享任一组成或任一 DOI 的传递连通组；用于分组划分，不能视为重复样品。',
        'cif_coordinates_status':'按上游 randomize_cifs.py 与 ICSD ID 推断；不是逐个结构恢复原始坐标。',
        'license_status':'按来源保留不同条款；合并文件不统一改为开放商业许可。',
        'elements':'解析组成中非零元素列表，分号分隔。','n_elements':'解析组成中非零元素数。',
        'stoichiometry_json':'解析的名义元素计数；名义混合配方不等于单一晶相。',
        'chemistry_class':'按元素集合启发式分类，见 common.chemistry；未作为实验相鉴定。',
        'eligible_numeric_label':'组成可解析且精确正数标签；仍可能有来源/质量问题，不等于可直接训练。',
    }
    for axis in ['a','b','c']:descriptions[f'lattice_{axis}_A']=f'OBELiX 原始晶格常数 {axis}，单位 Å；缺失不填。'
    for angle in ['alpha','beta','gamma']:descriptions[f'lattice_{angle}_deg']=f'OBELiX 原始晶格角 {angle}，单位度。'
    dictionary=[]
    for col in d.columns:
        desc=descriptions[col]
        dictionary.append({'field':col,'dtype':str(d[col].dtype),'meaning':desc,'non_null_count':int((d[col].notna()&d[col].astype(str).ne('')).sum())})
    pd.DataFrame(dictionary).to_csv(DEST/'data_dictionary.csv',index=False)
    summary={'source_counts':d.dataset.value_counts().to_dict(),'total_source_records':len(d),
             'clean_primary_numeric_records':len(clean),'reported_rt_candidates':len(rt),
             'parsed_composition_groups':len(comps),'cross_source_candidate_links':len(links),
             'confirmed_sample_duplicates_removed':0,'cif_files':len(cifs),
             'censored_records':int(d.sigma_relation.isin(['<','<=','>','>=']).sum()),
             'formula_review_records':int(d.formula_parse_status.ne('parsed').sum()),
             'leakage_groups':int(d.leakage_group_id.nunique()),
             'ddse_status':'not_obtained_bulk_export_restricted',
             'laskowski_status':'1345 conductivity cells from secondary transcription; reported target1346, not fully recovered or primary-verified'}
    (DEST/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
    # Excel strings must never be interpreted as formulae from external sources.
    sys.path.insert(0,str(ROOT/'experiments/vendor'))
    import openpyxl
    from openpyxl.cell.cell import ILLEGAL_CHARACTERS_RE
    from openpyxl.cell import WriteOnlyCell
    wb=openpyxl.Workbook(write_only=True)
    sheets=[('all_source_records',d),('clean_numeric',clean),('reported_RT_candidates',rt),
            ('cross_source_links',links),('structures',pd.DataFrame(cifs)),('aux_bulk_total',pd.DataFrame(aux)),
            ('composition_index',pd.DataFrame(comps)),('source_coverage',pd.DataFrame(coverage)),('data_dictionary',pd.DataFrame(dictionary))]
    for name,frame in sheets:
        ws=wb.create_sheet(name);ws.append(list(frame.columns))
        for row in frame.itertuples(index=False,name=None):
            cells=[]
            for value in row:
                if pd.isna(value):value=None
                elif isinstance(value,str):
                    value=ILLEGAL_CHARACTERS_RE.sub('',value)
                    cell=WriteOnlyCell(ws,value=value);cell.data_type='s';value=cell
                cells.append(value)
            ws.append(cells)
    wb.save(DEST/'unified_electrolytes.xlsx')
    return summary

def manifest():
    entries=[]
    jobs=[('obelix',RAW/'obelix/all.csv',f'https://github.com/NRC-Mila/OBELiX/tree/{SHA}','CC-BY-4.0','complete_public_table'),
          ('obelix_cifs',RAW/'obelix/all_cifs.zip',f'https://github.com/NRC-Mila/OBELiX/tree/{SHA}','CC-BY-4.0','321_released_cifs'),
          ('liion',SOURCES/'liion_original.html','https://pcwww.liverpool.ac.uk/~msd30/lmds/LiIonDatabase.html','University of Liverpool: academic attribution; commercial permission required','complete_public_table_corrected_2024'),
          ('mchaffie',RAW/'mchaffie_conductivity.csv','https://data.caltech.edu/records/23mvv-6gk43','CC0-1.0','complete_public_table'),
          ('laskowski',UPSTREAM/'data/misc/laskowski_semi-fromatted.csv',f'https://github.com/NRC-Mila/OBELiX/blob/{SHA}/data/misc/laskowski_semi-fromatted.csv','CC-BY-4.0 repository transcription; original SI rights unresolved','damaged_secondary_transcription_1345_label_cells'),
          ('laskowski_doi_lookup',UPSTREAM/'data/misc/doi_lookup_table.csv',f'https://github.com/NRC-Mila/OBELiX/blob/{SHA}/data/misc/doi_lookup_table.csv','CC-BY-4.0 repository','reference_lookup_not_independently_verified'),
          ('ddse_access_evidence',SOURCES/'ddse_browser_body.txt','https://www.digbat.org/','bulk data not publicly exported','public_UI_says_collaborator_only_bulk_copy')]
    for name,path,url,license,status in jobs:
        entries.append({'dataset':name,'local_path':str(path.relative_to(ROOT)),'url':url,'sha256':digest(path),
                        'bytes':path.stat().st_size,'license':license,'status':status,'retrieval_date':'2026-11-05'})
    entries.append({'dataset':'laskowski_primary_SI','url':'https://authors.library.caltech.edu/records/64j4e-7fx68',
                    'status':'not_downloaded: institution HTTP521; publisher SI HTTP404 / article HTTP403','reported_records':1346})
    entries.append({'dataset':'ddse','url':'https://www.digbat.org/','status':'not_downloaded: public UI explicitly restricts bulk copy to collaborators',
                    'reported_snapshot_materials_user':2448,'current_UI_experimental_materials_all_ions':3859,
                    'current_UI_experimental_entries_all_ions':28074,'current_UI_last_update':'2026.07.30','current_UI_is_not_Li_only':True})
    (DEST/'source_manifest.json').write_text(json.dumps(entries,ensure_ascii=False,indent=2))

def run():
    DEST.mkdir(parents=True,exist_ok=True)
    ob,cifs,aux=load_obelix();records=ob+load_liion()+load_caltech()+load_laskowski_transcription()
    d=normalize(records);links=link_records(d);summary=write_outputs(d,links,cifs,aux);manifest()
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    return summary

if __name__=='__main__':run()
