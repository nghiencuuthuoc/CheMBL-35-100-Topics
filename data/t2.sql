SELECT md.chembl_id, cs.canonical_smiles, act.pchembl_value
FROM activities act
JOIN molecule_dictionary md ON md.molregno = act.molregno
JOIN compound_structures cs ON cs.molregno = md.molregno
WHERE act.pchembl_value IS NOT NULL
LIMIT 100;
