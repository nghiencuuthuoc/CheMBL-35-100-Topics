-- Extract 100 IC50 entries with molfile
SELECT
    act.molregno,
    act.standard_value,
    act.standard_units,
    act.assay_id,
    cs.molfile
FROM
    activities act
JOIN
    compound_structures cs ON act.molregno = cs.molregno
WHERE
    act.standard_type = 'IC50'
    AND act.standard_units = 'nM'
    AND act.standard_value IS NOT NULL
    AND act.standard_value > 0
    AND act.standard_value < 10000
LIMIT 100;
