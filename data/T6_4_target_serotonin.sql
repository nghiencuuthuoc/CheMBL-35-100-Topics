-- T6_1_target_egfr.sql
SELECT DISTINCT mol.chembl_id, mol.molregno, cs.canonical_smiles, act.standard_value, act.standard_type
FROM compound_structures AS cs
JOIN molecule_dictionary AS mol ON cs.molregno = mol.molregno
JOIN activities AS act ON act.molregno = mol.molregno
JOIN assays AS a ON act.assay_id = a.assay_id
JOIN target_dictionary AS td ON a.tid = td.tid
WHERE td.pref_name ILIKE '%serotonin%'
  AND act.standard_type = 'IC50'
  AND act.standard_value::text ~ '^[0-9\\.]+$'
LIMIT 100;
