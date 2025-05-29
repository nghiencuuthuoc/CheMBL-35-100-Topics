-- T10_1_extract_smiles.sql
SELECT mol.molregno, mol.canonical_smiles
FROM compound_structures AS mol
JOIN compound_records AS cr ON mol.molregno = cr.molregno
WHERE mol.canonical_smiles IS NOT NULL
LIMIT 100;

