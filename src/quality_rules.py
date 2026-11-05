"""Versioned flags for unresolved source-table notes, independent of conductivity."""
import re

QUALITY_RULE_VERSION='2.0'
QUALITY_RULES={
 'explicit_incorrect':r'incorrect|not correct|wrong composition',
 'composition_pending':r'composition should be corrected|check the composition|composition check|true formula for lattice',
 'structure_match_pending':r'check cif file matching|lattice parameters not matching cif|lattice parameter should be checked|lattice parameter double check|structure may be different',
 'measurement_or_phase_pending':r'no conductivity|could be deleted|should be deleted|high T structure|non linear extrapolation',
 'duplicate_pending':r'duplicat',
 'label_pending':r'ic # check|IC number should be read',
}

def warning_reasons(note):
    note='' if note is None else str(note)
    return ';'.join(k for k,p in QUALITY_RULES.items() if re.search(p,note,re.I))
