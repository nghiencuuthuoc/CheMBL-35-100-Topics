-- File: T3_leadlike_extract.sql
SELECT md.chembl_id, cs.canonical_smiles
FROM compound_structures cs
JOIN molecule_dictionary md ON cs.molregno = md.molregno
WHERE cs.canonical_smiles IS NOT NULL
LIMIT 100;
