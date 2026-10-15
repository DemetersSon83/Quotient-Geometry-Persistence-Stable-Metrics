# Quotient geometry: reproducibility supplement

Computational supplement to **Quotient Geometry and Persistence-Stable Metrics for Swarm Configurations**, subtitle **Generic Local Recovery from Two Temporal H₀ Snapshots**, by Mark M. Bailey, PhD. This package targets the integrated manuscript revision dated 3 October 2026.

It reproduces the five validated certificates and Table 1, the centered rank-deficient control, interval and derivative checks, and twenty tree-witness diagnostics. See [the manuscript-to-code map](docs/RESULTS.md) for expected results and the distinction between proofs and finite checks.

## Run

Use Python 3.12; the reference local interpreter is CPython 3.12.14. The computation uses only the Python standard library. No package installation, network access, or LaTeX installation is needed.

From the repository root:

```sh
# Full regeneration, including the candidate search.
python3 reproduce.py

# Regenerate references and replay accepted certificates without rerunning the search.
python3 reproduce.py --mode verify

# Choose a new output directory; it must not already exist.
python3 reproduce.py --output /path/to/new-run
```

**Do not use `python -O` or `PYTHONOPTIMIZE`: assertions implement verification conditions.** The runner rejects optimized Python.

Each run writes to a fresh `outputs/<mode>-<timestamp>-<unique>/` directory by default. It preserves the supplied scripts and baseline results. Output includes:

- `report.json`: stage outcomes and comparisons with the supplied baselines;
- `summary.csv`: certificate results for inspection;
- logs and regenerated files in the original `certificates/`, `proof_checks/`, and root table layout.

JSON comparisons use the exact structural data, excluding `elapsed_seconds` in local certificates and `display_only` reference values. The generated LaTeX table is compared byte-for-byte. The report is the record of whether a particular run passed; successful execution alone is not a claim about unrun environments.

## What runs

| Stage | Script | Full | Verify |
|---|---|:---:|:---:|
| Exact references, analytic baselines, centered control | `certificates/reference_certificate.py` | ✓ | ✓ |
| Candidate radius/delay search | `certificates/local_certificate.py` | ✓ | — |
| Accepted-certificate replay and independent checks | `certificates/verify_local.py` | ✓ | ✓ |
| Tree-witness diagnostics | `proof_checks/universal_witness.py` | ✓ | ✓ |
| Conservatively rounded Table 1 | `build_certificate_table.py` | ✓ | ✓ |

`certificates/interval_core.py` supplies the shared outward-rounded interval arithmetic. Verify mode uses the supplied accepted certificate records as its replay inputs. Both modes regenerate the exact reference data, including the centered negative control.

The original scripts also remain directly runnable at their original paths. They write beside their source files by default; use `reproduce.py` to keep baseline files unchanged.

## Expected results

- Five accepted certificates: regular references with `n = 3, 4, 6, 8`, and a near-centered `n = 6` reference.
- Exact positive-reference ranks `3, 5, 9, 13, 9`; centered control rank `8` in a quotient of dimension `9`.
- 300 elementary interval cases, one inverse-square-root derivative check, five reference Jacobian checks, and five independent temporal-derivative checks.
- Five accepted-certificate replays and five rejected oversized candidates.
- Twenty witness cases across paths, stars, binary trees, and seeded random trees at `n = 3, 4, 8, 16, 32`; twelve also receive a direct dense determinant check.
- Exact regeneration of the five-row [printed certificate table](certificate_table.tex). Values, gains, and manuscript labels are listed in [docs/RESULTS.md](docs/RESULTS.md).

Full regeneration and exact baseline comparisons passed locally on 5 October 2026 with CPython 3.12.14. See [validation/report.json](validation/report.json), [the replay report](validation/verify_report.json), and [runner checks](validation/packaging_controls.json). The five-row LaTeX table matched byte-for-byte. These reports and the accompanying logs record actual local executions.

GitHub Actions is configured for full reproduction. The workflow has not been run remotely. Individual local runs also produce their own `report.json`.

## Interpretation and scope

The all-tree rank theorem and the sharp `Θ(1/|h|)` conditioning law are proved in the manuscript. The finite software checks support reproducibility; they do not establish those universal claims by sampling.

The certificates give sufficient local guarantees under known homogeneous complete spherical consensus and known sampling times. They cover full formation-metric balls and all nonzero delays up to each certified endpoint. Reflection ambiguity remains. The bounds are not optimal condition numbers or maximal operating limits, and a rejected candidate does not establish nonidentifiability. No estimator for arbitrary noisy diagrams, global reconstruction, or recovery under unknown dynamics is claimed.

Candidate selection is heuristic. Final acceptance uses exact rational calculations and outward-rounded intervals with 96-bit dyadic endpoints. Radii and delay endpoints in Table 1 are rounded downward; derivative errors and inverse bounds are rounded upward. The displayed inverse constant applies at the printed delay endpoint; use the manuscript's `K_h` formula at smaller delays.

## Provenance

[provenance.json](provenance.json) identifies the source archive and manuscript revision and records SHA-256 hashes of the original study scripts and baseline outputs. Those files are preserved byte-for-byte. The isolated runner, documentation, and CI configuration are packaging additions.

This repository contains the computational supplement. The manuscript remains the source of the proofs. Rocq sources and reviewer reports are outside this package's scope; these software checks are not proof-assistant formalization or external peer review.

No public repository or archival identifier has been assigned, and no license grant has been selected. See [GitHub setup](docs/GITHUB.md) for publishing this initialized repository to a destination you control.
