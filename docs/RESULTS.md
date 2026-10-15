# Manuscript-to-code map

This map targets the **3 October 2026 integrated manuscript**, `template.tex`, identified by its hash in [provenance.json](../provenance.json). That manuscript contains one numerical table and no figures. All finite computational claims are in Section 6, **Exact examples and validated computation** (`sec:computation`). The following are expected results from the supplied exact records; a run's `report.json` records its validation outcome.

## Computational results

| Manuscript location or claim | Reproduction source | Expected result |
|---|---|---|
| Section 6: regular rational references and exact first-jet ranks | `certificates/reference_certificate.py`; `reference_certificates.json` | First `n = 3, 4, 6, 8` points of the stated stereographic family; ranks `3, 5, 9, 13` |
| Section 6: centered control | Same script, `centered_negative_control()` | Six rational nonantipodal points, sum exactly zero, selected-tree rank `8 < 9`; repeated pair distances exclude the strict chamber |
| Section 6: near-centered positive reference | Same script, `near_centered_points()` | Perturb appended stereographic coordinates at vertices `i = 3, 4, 5` by `((i+1)ε, (i+2)ε)`, with `ε = 1/10000`; exact rank `9` |
| Appendix A (`app:baseline`): instantiated analytic baseline | Same script, `quantitative_bounds()` | Exact Gram inverse, certified positive singular-value lower bound, and baseline radius/delay values for all five positive references |
| Theorem 5.1 (`loc:thm:local`), instantiated in Section 6 | `certificates/local_certificate.py`; `local_certificates.json` | Five accepted interval certificates with exact left inverses, positive tube/MST/matching margins, and total derivative error `θ ≤ 1/2` |
| Table 1 (`tab:certificates`) | `build_certificate_table.py`; `certificate_table.tex`; `certificates/printed_bounds_verification.json` | Five rows with conservative rounding; inverse constant recomputed at the printed delay endpoint |
| Section 6: improvement over Appendix A | `radius_gain`, `delay_gain` in `local_certificates.json` | Simultaneous gains at the selected radius/delay pairs, listed below |
| Section 6: interval-kernel checks | `certificates/verify_local.py`; `verification.json` | 300 seeded rational cases and one exact inverse-square-root derivative check |
| Section 6: reference Jacobian and temporal-derivative checks | Same verification script | Five exact Jacobian-enclosure comparisons; five comparisons of `DF₂` with an independently written Hessian/product-rule formula |
| Section 6: accepted and oversized candidate checks | Same verification script | Five final certificates recomputed; five candidates with both radius and delay multiplied by 100 rejected |
| Section 6: witness diagnostics supporting Theorem 3.1 | `proof_checks/universal_witness.py`; `results.json` | Twenty exact cases; twelve direct dense determinant comparisons; details below |

The regular stereographic parameters are `(0,0), (1/2,0), (1/5,3/5), (-3/5,2/7), (-2/5,-4/5), (4/5,-3/5), (5/4,2/5), (-3/2,-1/4)`, with `p = (2a, 2b, 1-a²-b²)/(1+a²+b²)`. All positive references have `p₁ = (4/5, 0, 3/5)`. Exact coordinates, tangent blocks, selected MST edges, Jacobians, and inverse Gram matrices are retained in `certificates/reference_certificates.json`.

## Table 1

| Reference | Printed formation radius η | Printed delay τ | θ upper | K at printed τ, upper |
|---|---:|---:|---:|---:|
| `n = 3` | `1.872e-4` | `2.104e-3` | `0.496` | `1.04e4` |
| `n = 4` | `5.773e-5` | `5.519e-4` | `0.497` | `1.98e5` |
| `n = 6` | `3.027e-5` | `4.322e-4` | `0.498` | `1.73e5` |
| `n = 8` | `5.993e-6` | `1.204e-4` | `0.499` | `2.66e6` |
| Near-centered `n = 6` | `6.671e-9` | `3.443e-4` | `0.499` | `1.90e9` |

Radii are in radians and delays in the manuscript's consensus time units. The guarantee covers every formation in the open ball of printed radius and every `0 < |h| ≤ τ`. Radius and delay are rounded downward; the two upper bounds are rounded upward. The exact records contain the larger unrounded radius and delay values.

Equation `loc:eq:inverse` gives

```text
K_h = 2/(1 - θ) × (c₀ + 2c₁/|h|).
```

The table generator evaluates this expression at the **printed** delay before rounding its result upward. The endpoint constant must not be reused for smaller delays.

## Baseline comparison

These values are derived from exact rational records; decimals here are for inspection. Gains refer to the unrounded certificates, not ratios of the displayed table entries.

| Reference | Baseline η₀, approximately | Baseline τ₀, approximately | Radius gain | Delay gain |
|---|---:|---:|---:|---:|
| `n = 3` | `4.2550689620e-6` | `3.5175236753e-5` | `43.996875` | `59.8383279419` |
| `n = 4` | `4.2874607090e-7` | `2.6582256396e-6` | `134.6625` | `207.6353723145` |
| `n = 6` | `4.7122424631e-7` | `2.3187224819e-6` | `64.24375` | `186.4323510742` |
| `n = 8` | `7.6584543582e-8` | `2.8263343465e-7` | `78.25625` | `426.3131494141` |
| Near-centered `n = 6` | `1.2108630696e-10` | `5.9582151045e-10` | `55.1` | `577888.8` |

These reproduce the manuscript's regular-case claims of radius gains of about **44–135** and delay gains of about **60–426**. The near-centered delay gain is about **5.78 × 10⁵**, while its exact certified formation radius is about **6.6718555136 × 10⁻⁹** radians, below **7 × 10⁻⁹**. A large delay gain does not imply spatial robustness.

The baseline uses a 32-bit dyadic lower bound for the relevant reciprocal square root. The local implementation uses outward-rounded intervals with endpoints in `2^-96 Z`. The radius/delay search expands both parameters, then holds the radius fixed while expanding the delay. Its selected pair is sufficient, not proven optimal.

## Witness and arithmetic diagnostics

| Witness size n | Tree cases | Expected rank 2n−3 | Direct dense determinant comparison |
|---:|---|---:|:---:|
| 3 | Path, star, binary, seeded random | 3 | Yes |
| 4 | Path, star, binary, seeded random | 5 | Yes |
| 8 | Path, star, binary, seeded random | 13 | Yes |
| 16 | Path, star, binary, seeded random | 29 | No |
| 32 | Path, star, binary, seeded random | 61 | No |

All twenty cases check the exact block identity and nonzero factorized witness minor. Twelve also compare the complete dense rational determinant with the factorized value. The random trees use a seeded Prüfer construction with seed `20261002`. The witness need not have the tested tree as its MST.

The 300 elementary interval cases independently use seed `20261002`. They exercise endpoint containment, addition, multiplication, squaring, division where the denominator excludes zero, and square roots. The inverse-square-root derivative check encloses the exact value `4/5` and derivative `-48/125` at the stated test point. These are finite software checks, not an exhaustive proof of the arithmetic implementation.

## Proofs and finite evidence

| Manuscript result | What establishes it | Role of computation |
|---|---|---|
| Theorem 3.1 (`thm:rank`): generic rank `2n−3` for every labeled tree and every `n ≥ 3` | Pole–equator witness, explicit symbolic minor, and genericity argument | Twenty finite diagnostics of the witness construction |
| First-jet identifiability (`cor:jet`) and two sufficiently close snapshots (`cor:two`) | Observation changes, inverse function theorem, and uniform perturbation argument | Exact rank checks at five chosen references |
| Theorem 4.1 (`thm:rate`): optimal inverse-conditioning order `Θ(1/|h|)` | Analytic upper and lower bounds | No empirical scaling fit is claimed or required |
| Theorem 5.1 (`loc:thm:local`): uniform local recovery and full formation-ball coverage | Analytic implication from validated inequalities | Five exact/interval instantiations of those inequalities |
| Proposition 7.1 (`prop:centered`): centered rank at most `n+2` | Factorization of the derivative rows through the three-dimensional sum variation | One exact centered example with rank `8 < 9` |
| Appendix A (`app:baseline`) | Closed-form analytic sufficient bounds | Exact baseline evaluation at five references |

The geometric and forward-stability results in Section 2, including the corrected two-point formulas, are analytic background. They are not numerical experiments.

## Reproduction boundaries

- Full mode regenerates references, reruns candidate selection, replays the accepted certificates, runs witness diagnostics, and regenerates the table.
- Verify mode regenerates references and the centered control, then replays the supplied final positive certificates without rerunning candidate selection. It also runs witness diagnostics and table generation.
- `elapsed_seconds` is machine-dependent and excluded from local-certificate comparisons. Reference `display_only` floating values are excluded; exact rational data and substantive checks are compared. The LaTeX table is compared byte-for-byte.
- The original scripts and records are immutable baselines for the runner and are identified by [provenance.json](../provenance.json). A local run report records the packaging validation; historical supplied reports describe the original internal checks.
- No current-manuscript claim requires an omitted exploratory rank survey, numerical conditioning fit, application-performance experiment, or noisy-diagram estimator. The superseded phase-circle argument is outside this supplement's scope.

Recovery is local, under the specified known dynamics and known sampling times. Reflection ambiguity remains; centered rank deficiency alone does not exclude nonlinear local injectivity. Certificate failure does not establish nonidentifiability, and no global or unknown-dynamics reconstruction result is claimed. Software checks are additional internal verification, not external peer review or proof-assistant formalization. Rocq sources and reviewer reports are excluded.
