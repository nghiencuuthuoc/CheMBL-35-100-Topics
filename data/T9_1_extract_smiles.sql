-- file: T9_1_extract_smiles.sql
SELECT mol.molregno, mol.canonical_smiles
FROM compound_structures AS mol
JOIN compound_properties AS prop ON mol.molregno = prop.molregno
WHERE mol.canonical_smiles IS NOT NULL
LIMIT 100;
