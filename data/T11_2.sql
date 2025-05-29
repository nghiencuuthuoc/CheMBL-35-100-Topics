-- Only extract compounds from approved drugs
SELECT 
    m.molregno, cs.canonical_smiles
FROM compound_structures cs
JOIN molecule_dictionary m ON cs.molregno = m.molregno
WHERE m.max_phase = 4
LIMIT 1000;
