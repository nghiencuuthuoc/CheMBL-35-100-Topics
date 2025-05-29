-- T5_1_extract_drug_like_data.sql
SELECT md.chembl_id, cs.canonical_smiles, act.standard_value::float AS ic50_nm
FROM activities act
JOIN compound_structures cs ON act.molregno = cs.molregno
JOIN molecule_dictionary md ON md.molregno = act.molregno
WHERE act.standard_type = 'IC50'
  AND act.standard_value IS NOT NULL
  AND act.standard_value::text ~ '^[0-9\.]+$'
LIMIT 100;
