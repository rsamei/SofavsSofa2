"""Paired DeLong AUROC contrasts in BigQuery, with mid-rank placement values.

For a contrast sum_k c_k * AUC_k, the per-subject combined placement value sum_k c_k * V_k gives
Var = Var(V10)/m + Var(V01)/n, the DeLong covariance applied to the contrast.
"""
from bq_client import PROJECT, query


def delong_sql(source, sets, sys_expr, contrasts, y="hospital_expire_flag"):
    """source: FROM expression (a backquoted table or a parenthesised subquery with an alias);
    sets: {name: SQL condition}; sys_expr: {system: SQL expression}; contrasts: {name: {system: coefficient}}.
    Every system also gets its own AUROC row (contrast 'auc_<system>')."""
    systems = list(sys_expr)
    sets_sql = ",\n      ".join(f"STRUCT('{k}' AS an, ({v}) AS inc)" for k, v in sets.items())
    sys_sql = ",\n      ".join(f"STRUCT('{s}' AS sys, {sys_expr[s]} AS x)" for s in systems)
    pv_cols = ",\n    ".join(f"MAX(IF(sys = '{s}', v, NULL)) AS v_{s}" for s in systems)
    defs = []
    for nm, coef in [(f"auc_{s}", {s: 1}) for s in systems] + list(contrasts.items()):
        coefs = ", ".join(f"{float(coef.get(s, 0))} AS k_{s}" for s in systems)
        defs.append(f"STRUCT('{nm}' AS contrast, {coefs})")
    combo = " + ".join(f"d.k_{s} * p.v_{s}" for s in systems)
    return f"""
WITH long AS (
  SELECT st.an, sy.sys, stay_id, {y} AS y, sy.x
  FROM {source}, UNNEST([
      {sets_sql}
    ]) AS st, UNNEST([
      {sys_sql}
    ]) AS sy
  WHERE st.inc
),
cc AS (SELECT an, stay_id FROM long GROUP BY an, stay_id HAVING COUNTIF(x IS NULL) = 0),
r AS (
  SELECT l.an, l.sys, l.stay_id, l.y, l.x,
    RANK() OVER (PARTITION BY l.an, l.sys ORDER BY l.x) + (COUNT(*) OVER (PARTITION BY l.an, l.sys, l.x) - 1) / 2.0 AS pr,
    RANK() OVER (PARTITION BY l.an, l.sys, l.y ORDER BY l.x) + (COUNT(*) OVER (PARTITION BY l.an, l.sys, l.y, l.x) - 1) / 2.0 AS cr,
    COUNTIF(l.y = 1) OVER (PARTITION BY l.an, l.sys) AS m, COUNTIF(l.y = 0) OVER (PARTITION BY l.an, l.sys) AS n
  FROM long l JOIN cc USING (an, stay_id)
),
pl AS (SELECT an, sys, stay_id, y, IF(y = 1, (pr - cr) / n, 1 - (pr - cr) / m) AS v FROM r),
pv AS (SELECT an, stay_id, ANY_VALUE(y) AS y, {pv_cols} FROM pl GROUP BY an, stay_id),
defs AS (SELECT * FROM UNNEST([
      {",".join(defs)}
    ])),
cv AS (SELECT p.an, d.contrast, p.y, {combo} AS c FROM pv p CROSS JOIN defs d)
SELECT an, contrast, COUNTIF(y = 1) AS n_pos, COUNTIF(y = 0) AS n_neg,
  AVG(IF(y = 1, c, NULL)) AS est,
  SQRT(VAR_SAMP(IF(y = 1, c, NULL)) / COUNTIF(y = 1) + VAR_SAMP(IF(y = 0, c, NULL)) / COUNTIF(y = 0)) AS se
FROM cv GROUP BY an, contrast ORDER BY an, contrast
"""


def run_delong(source, sets, sys_expr, contrasts, y="hospital_expire_flag"):
    """Returns rows {set, contrast, n_pos, n_neg, est, se, lo, hi}."""
    f, rows, _ = query(PROJECT, delong_sql(source, sets, sys_expr, contrasts, y))
    out = []
    for r in rows:
        d = dict(zip(f, r))
        est = float(d["est"])
        se = float(d["se"]) if d["se"] is not None else 0.0
        out.append({"set": d["an"], "contrast": d["contrast"], "n_pos": int(d["n_pos"]), "n_neg": int(d["n_neg"]),
                    "est": est, "se": se, "lo": est - 1.96 * se, "hi": est + 1.96 * se})
    return out
