#!/usr/bin/env python3
"""Pure-stdlib source data verification; run inside repository candidate or beside its data/ directory.
Does not perform DFT; does not alter frozen source values.
"""
from pathlib import Path
import csv, json, statistics, hashlib, tarfile, math, re, sys
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'
if not DATA.is_dir(): DATA=ROOT/'source_tables'
assert DATA.is_dir(),f'Missing data directory {DATA}'
EH_TO_EV=27.211386245988

def rows(n):
 with (DATA/n).open(encoding='utf-8-sig', newline='') as f:return list(csv.DictReader(f))
def near(a,b,t=1e-8):
 if not math.isclose(float(a),float(b),rel_tol=0,abs_tol=t):
  raise AssertionError(f'expected {b}, got {a}; tolerance {t}')

m1e=rows('M1_ENTITIES.csv');m1=rows('M1_SHIFTS.csv')
assert len(m1e)==65 and len(m1)==38
entities={r['entity_id']:r for r in m1e}
for r in m1e:
 near((float(r['reduced_Eh'])-float(r['neutral_Eh']))*EH_TO_EV,r['reduction_energy_eV'],1e-8)
for r in m1:
 par=entities[r['parent_id']];add=entities[r['adduct_id']]
 near(float(add['reduction_energy_eV'])-float(par['reduction_energy_eV']),r['delta_delta_Ered_eV'],1e-9)
assert (sum(float(r['delta_delta_Ered_eV'])<0 for r in m1),sum(float(r['delta_delta_Ered_eV'])>0 for r in m1))==(27,11)
exceptions={r['adduct_id'] for r in rows('M1_FOUR_RAW_OUTPUT_EXCEPTIONS.csv')}
assert exceptions=={'Q0059_R0_SME1','Q0059_R0_SME2','Q0297_R3_SME1','Q0471_R0_SME1'}
assert len([r for r in m1 if r['adduct_id'] not in exceptions])==34

m2=rows('M2R_CORE.csv');m3=rows('M3R_SHIFTS.csv');vintage=rows('M2R_10SITE_VINTAGE_SENSITIVITY.csv')
assert len(m2)==10 and len(m3)==4 and len(vintage)==10
median_original=statistics.median(abs(float(r['difference_vs_sme_eV'])) for r in m2)
median_mixed=statistics.median(abs(float(r['corrected_diff'])) for r in vintage)
near(median_original,0.17589955941670565,1e-12)
near(median_mixed,0.13786928312282415,1e-12)
corrected={r['site_id']:r for r in m3}
for r in vintage:
 if r['site'] in corrected:
  near(r['corrected_accys'],corrected[r['site']]['m3r_final_shift_eV'],1e-12)

m5=rows('M5_DECOMPOSITION.csv');m56=rows('M5_M6_COMPARISON.csv');m6=rows('M6_24_ENERGIES.csv');geo=rows('M6_CYSTEINE_AMIDE_CONTACT_GEOMETRY.csv')
assert len(m5)==len(m56)==4 and len(m6)==24 and len(geo)==8
for r in m5:
 s=float(r['sme_relaxed_shift_eV']);c=float(r['counterfactual_sme_on_accysome_geometry_shift_eV']); a=float(r['accysome_m3r_shift_eV'])
 near(c-s,r['geometry_mediated_component_eV'],1e-10)
 near(a-c,r['sidechain_chemical_component_eV'],1e-10)
 near(a-s,r['total_accysome_minus_sme_eV'],1e-10)
 near(float(r['decomposition_residual_eV']),0,1e-9)
energies={(r['site'],r['kind'],r['state']):r for r in m6}
assert len(energies)==24
for r in m6:
 assert r['status']=='PASS' and r['state'] in ('neutral','reduced')
 assert (int(r['charge']),int(r['multiplicity']))==((0,1) if r['state']=='neutral' else (-1,2))
for q in m56:
 site=q['site']
 red={kind:(float(energies[(site,kind,'reduced')]['energy_hartree'])-float(energies[(site,kind,'neutral')]['energy_hartree']))*EH_TO_EV for kind in ('parent','sme','accysome')}
 sme=red['sme']-red['parent'];acc=red['accysome']-red['parent']
 near(sme,q['m6_sme_shift_eV'],3e-9)
 near(acc,q['m6_accysome_shift_eV'],3e-9)
 near(acc-sme,q['m6_full_minus_sme_eV'],3e-9)
 assert (sme*acc<0)==(q['m6_sign_reversal'].lower()=='true')

# Full M6 raw-file checks, including exact input geometry SHA-256 references.
archive=ROOT/'raw_m6/QUINONE_M6_COMPLETE_REPRODUCIBILITY_PACKAGE.tar.gz'
raw_count=0
if archive.exists():
 with tarfile.open(archive,'r:gz') as tar:
  names=set(tar.getnames())
  for r in m6:
   stem=f"quinone_m6_validation/jobs/{r['site']}__{r['kind']}__{r['state']}"
   for ending in ('input.xyz','job.inp','job.out'):
    assert stem+'/'+ending in names,f'Missing M6 raw {stem}/{ending}'
   xyz=tar.extractfile(stem+'/input.xyz').read()
   assert hashlib.sha256(xyz).hexdigest()==r['sha256'],stem+' geometry hash mismatch'
   out=tar.extractfile(stem+'/job.out').read()
   assert b'ORCA TERMINATED NORMALLY' in out,stem+' not terminated normally'
   raw_count+=1

print(json.dumps({'status':'PASS','m1_entities':len(m1e),'m1_site_shifts':len(m1),'m1_raw_exception_sites':len(exceptions),'m2r_sites':len(m2),'m2r_original_median_abs_eV':median_original,'m2r_mixed_vintage_median_abs_eV':median_mixed,'m5_cases':len(m5),'m6_energies':len(m6),'m6_raw_triples_verified':raw_count,'experimental_endpoint_match':False,'notice':'Downstream arithmetic and M6 raw-file verification only; not a rerun or upstream M1-M5 raw recovery.'},indent=2))
