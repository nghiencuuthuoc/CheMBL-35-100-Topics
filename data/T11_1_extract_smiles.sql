-- File: T11_1_extract_smiles.sql
SELECT 
    m.molregno, 
    cs.canonical_smiles
FROM compound_structures cs
JOIN molecule_dictionary m ON cs.molregno = m.molregno
LIMIT 100;
