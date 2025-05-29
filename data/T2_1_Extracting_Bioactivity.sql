SELECT md.chembl_id,
       cs.canonical_smiles,
       act.standard_value::float AS ic50_nM,   -- Tên cột rõ ràng
       act.standard_type
FROM activities act
JOIN compound_structures cs ON act.molregno = cs.molregno
JOIN molecule_dictionary md ON md.molregno = cs.molregno
WHERE act.standard_type = 'IC50'
  AND act.standard_value::text ~ '^[0-9\.]+$'
  AND act.standard_value::float > 0
LIMIT 100;
