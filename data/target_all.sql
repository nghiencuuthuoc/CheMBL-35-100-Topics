SELECT pref_name, target_type, organism, COUNT(DISTINCT act.activity_id) AS n_activities
FROM target_dictionary td
JOIN assays a ON td.tid = a.tid
JOIN activities act ON a.assay_id = act.assay_id
GROUP BY pref_name, target_type, organism
ORDER BY n_activities DESC
-- LIMIT 20
