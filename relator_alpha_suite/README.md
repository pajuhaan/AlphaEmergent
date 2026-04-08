# Relator Alpha Suite

This directory contains a modular Python implementation of the load-bearing
formulas used in the manuscript **"Alpha — The Emergent Fine-Structure Constant"**.

## Design goals

- Modular and independently executable scripts.
- English comments and referee-friendly notation.
- Article-standard symbols, including Unicode where helpful.
- No external CODATA `alpha` value is inserted into the default production path.
- The universal-QED route is implemented strictly as a one-way cross-check.
- Default pipeline output reproduces the article-style numerical tables.
- Terminal tables are rendered with `rich` using protected number formatting,
  so fixed-point and scientific notation are not visually truncated.

## File layout

- `compute_lambda_geom.py`  
  Computes the vector-channel representative `Λ_geom`.

- `compute_dc_article.py`  
  Computes `D_C(α)` from the article scalar closure. Supports the
  `current_*`, `refined_*`, and `realized_*` evaluator families.

- `compute_dc_qed.py`  
  Reconstructs the QED-induced scalar branch from universal pure-photonic
  `A_1^(2n)` coefficients.

- `solve_alpha_alp.py`  
  Solves the ALP closure for `α` using either the article branch or the
  QED-induced cross-check branch.

- `run_relator_alpha_pipeline.py`  
  Runs the publication-style numerical pipeline and prints the article tables in
  order: scalar block, vector block, internal `α` lock, external comparisons,
  pure-photonic `g-2` diagnostics, and appendix sensitivity / stability audits.

## Important note on the article conventions

The manuscript uses two closely related local conventions for the visible-source
renormalization factor associated with `ρ_dyn`.

To reproduce the published tables *numerically* without mixing conventions, the
code separates three scalar model families:

- `current_*`  
  visible baseline, with the current-text source normalization;
- `refined_*`  
  physical refined branch used for the physical-root and `g-2` tables;
- `realized_*`  
  appendix / realized-lock convention used where the manuscript's local series
  and sensitivity tables follow the alternate normalization path.

This split is deliberate and is what allows the default pipeline to match the
article tables cleanly instead of forcing a single convention onto all tables.

## Important honesty note about `Λ_geom`

The manuscript does **not** reduce the full operator-exact `Λ_geom` to a
one-line closed formula. Therefore the default code computes the article's
**diagnostic representative**

```text
Λ_geom^(mean) = Λ_ind + ΔΛ_UV→IR^(mean)(ℓ₀) + ΔΛ_out(η₀) R_χ^(mean).
```

Optional hooks are exposed for custom `P_IR^(χ)` and `R_χ` values if you later
replace the representative inputs by a stronger operator evaluation.

## Quick start

```bash
pip install -r requirements.txt
python run_relator_alpha_pipeline.py
```

The default pipeline prints the numerical tables corresponding to:

- fixed scalar geometry and static data,
- visible dynamic amplitudes and dynamic block,
- hidden Jacobi diagnostics,
- scalar physical roots and scalar series coefficients,
- vector exact-status and source ledgers,
- internal `α` lock and external benchmark comparisons,
- pure-photonic Relator `g-2` tables,
- appendix exact sensitivities, numerical sensitivities,
  rank stability, and precision stability.

You may also execute each module independently:

```bash
python compute_lambda_geom.py
python compute_dc_article.py --model current_rank5
python compute_dc_qed.py --order 5
python solve_alpha_alp.py --branch article --model current_rank5
python solve_alpha_alp.py --branch qed --order 5
```

## Dependencies

See `requirements.txt`.
