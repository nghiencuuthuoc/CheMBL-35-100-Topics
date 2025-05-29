-- T1: Trích xuất dữ liệu IC50 từ ChEMBL (phiên bản sửa lỗi)
SELECT 
    md.chembl_id, 
    cs.canonical_smiles, 
    act.standard_value AS ic50_nM,
    act.standard_type,
    act.standard_relation
FROM activities act
JOIN compound_structures cs ON act.molregno = cs.molregno
JOIN molecule_dictionary md ON md.molregno = cs.molregno
WHERE act.standard_type = 'IC50'
  AND act.standard_units = 'nM'
  AND act.standard_value::text ~ '^[0-9\.]+$'
  AND act.standard_value::float BETWEEN 1 AND 10000;
