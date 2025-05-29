-- File: T4_1_extract_smiles.sql
SELECT md.chembl_id, cs.canonical_smiles
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
WHERE cs.canonical_smiles IS NOT NULL
LIMIT 100;