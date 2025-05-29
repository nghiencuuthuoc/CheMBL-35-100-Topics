-- File: T7_1_get_egfr_ligands.sql
SELECT 
  mol.molregno, 
  cs.canonical_smiles, 
  act.standard_value AS ic50_nm
FROM activities AS act
JOIN assays AS ass ON act.assay_id = ass.assay_id
JOIN target_dictionary AS td ON ass.tid = td.tid
JOIN molecule_dictionary AS mol ON act.molregno = mol.molregno
JOIN compound_structures AS cs ON mol.molregno = cs.molregno
WHERE td.pref_name ILIKE '%EGFR%'
  AND act.standard_type = 'IC50'
  AND act.standard_units = 'nM'
  AND act.standard_value::text ~ '^[0-9\\.]+$'
LIMIT 100;
