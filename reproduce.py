#!/usr/bin/env python3
"""Reproduce the manuscript computations in a fresh, isolated directory.

No external packages are required. Original study files remain unchanged.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    "certificates/interval_core.py",
    "certificates/reference_certificate.py",
    "certificates/local_certificate.py",
    "certificates/verify_local.py",
    "proof_checks/universal_witness.py",
    "build_certificate_table.py",
]
OUTPUTS = [
    "certificates/reference_certificates.json",
    "certificates/local_certificates.json",
    "certificates/verification.json",
    "proof_checks/results.json",
    "certificates/printed_bounds_verification.json",
    "certificate_table.tex",
]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def check_source_integrity():
    manifest = read_json(ROOT / "provenance.json")
    for item in manifest["files"]:
        path = ROOT / item["path"]
        if not path.is_file() or sha256(path) != item["sha256"]:
            raise RuntimeError("Preserved study file differs from provenance: " + item["path"])
    return len(manifest["files"])


def mathematical_content(path, value):
    """Remove only documented display/runtime fields, never certificate data."""
    if path == "certificates/local_certificates.json":
        value = {key: item for key, item in value.items() if key != "elapsed_seconds"}
    if path == "certificates/reference_certificates.json":
        value = dict(value)
        value["certificates"] = [
            {key: item for key, item in row.items() if key != "display_only"}
            for row in value["certificates"]
        ]
    return value


def compare_output(work, relative_path):
    expected, actual = ROOT / relative_path, work / relative_path
    if not actual.is_file():
        raise RuntimeError("Expected output was not generated: " + relative_path)
    if relative_path.endswith(".json"):
        left = mathematical_content(relative_path, read_json(expected))
        right = mathematical_content(relative_path, read_json(actual))
        equal = left == right
        method = "exact parsed JSON comparison with documented exclusions"
    else:
        equal = expected.read_bytes() == actual.read_bytes()
        method = "byte-for-byte comparison"
    if not equal:
        raise RuntimeError("Reproduction differs from manuscript baseline: " + relative_path)
    return {"path": relative_path, "status": "PASS", "method": method,
            "baseline_sha256": sha256(expected), "generated_sha256": sha256(actual)}


def run_step(work, script):
    # -E also protects children when a caller ignored PYTHONOPTIMIZE only
    # in the parent interpreter. The study programs depend on assertions.
    command = [sys.executable, "-E", "-B", script]
    log = work / "logs" / (Path(script).stem + ".log")
    print("Running " + script, flush=True)
    started = time.perf_counter()
    with log.open("w", encoding="utf-8") as stream:
        result = subprocess.run(command, cwd=work, stdout=stream,
                                stderr=subprocess.STDOUT, check=False)
    record = {"script": script, "command": ["python3", "-E", "-B", script],
              "exit_code": result.returncode,
              "elapsed_seconds": round(time.perf_counter() - started, 6),
              "log": str(log.relative_to(work))}
    if result.returncode:
        print(log.read_text(encoding="utf-8"), file=sys.stderr)
        raise RuntimeError(f"{script} failed with exit code {result.returncode}; see {log}")
    return record


def summarize(work):
    references = read_json(work / OUTPUTS[0])
    local = read_json(work / OUTPUTS[1])
    verification = read_json(work / OUTPUTS[2])
    witnesses = read_json(work / OUTPUTS[3])
    printed = {row["label"]: row for row in read_json(work / OUTPUTS[4])}
    rows = []
    for row in local["results"]:
        bounds = printed[row["label"]]
        rows.append({"label": row["label"], "n": row["n"],
                     "eta_print_exact": bounds["eta_down"],
                     "tau_print_exact": bounds["tau_down"],
                     "theta_print_upper_exact": bounds["theta_up"],
                     "K_at_printed_tau_upper_exact": bounds["K_at_printed_tau_up"],
                     "radius_gain_exact": row["radius_gain"],
                     "delay_gain_exact": row["delay_gain"],
                     "radius_gain_display": float(Fraction(row["radius_gain"])),
                     "delay_gain_display": float(Fraction(row["delay_gain"]))})
    with (work / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return {
        "accepted_certificates": len(local["results"]),
        "reference_ranks": {row["label"]: row["checks"]["rank"]
                            for row in references["certificates"]},
        "centered_control_rank": references["negative_controls"][0]["jet_rank_exact"],
        "centered_control_quotient_dimension": references["negative_controls"][0]["quotient_dimension"],
        "elementary_interval_cases": verification["elementary_rational_interval_cases"],
        "independent_temporal_derivative_checks": sum(
            row["temporal_derivative_hessian_crosscheck"] for row in verification["results"]),
        "oversized_candidates_rejected": sum(
            row["enlarged_candidate_not_certified"] for row in verification["results"]),
        "tree_witness_cases": len(witnesses["results"]),
        "dense_determinant_comparisons": sum(
            row["complete_minor_checked_directly"] for row in witnesses["results"]),
        "printed_table_rows": len(printed),
    }


def main():
    if sys.flags.optimize or not __debug__:
        print("ERROR: Assertions must be enabled. Remove -O/-OO and PYTHONOPTIMIZE.", file=sys.stderr)
        return 2
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["full", "verify"], default="full",
                        help="full: regenerate candidate search; verify: replay stored final candidates")
    parser.add_argument("--output", type=Path, help="new output directory; must not already exist")
    args = parser.parse_args()
    now = datetime.now(timezone.utc)
    name = args.mode + "-" + now.strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]
    work = (args.output or ROOT / "outputs" / name).resolve()
    if work.exists():
        parser.error("Output directory already exists; choose a new directory: " + str(work))
    try:
        original_count = check_source_integrity()
    except (OSError, ValueError, RuntimeError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        return 1
    work.mkdir(parents=True)
    (work / "logs").mkdir()
    report = {"status": "RUNNING", "mode": args.mode, "started_utc": now.isoformat(),
              "python_version": platform.python_version(), "python_implementation": platform.python_implementation(),
              "platform": platform.platform(), "assertions_enabled": __debug__,
              "preserved_study_files_verified": original_count,
              "candidate_search_regenerated": args.mode == "full", "steps": [], "comparisons": [],
              "comparison_exclusions": ["local_certificates.json: top-level elapsed_seconds",
                                         "reference_certificates.json: per-certificate display_only"],
              "scope": "Finite exact-arithmetic manuscript computations; not formal verification of the universal theorems."}
    started = time.perf_counter()
    try:
        for relative_path in SCRIPTS:
            target = work / relative_path
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / relative_path, target)
        if args.mode == "verify":
            shutil.copyfile(ROOT / OUTPUTS[1], work / OUTPUTS[1])
        sequence = ["certificates/reference_certificate.py"]
        if args.mode == "full":
            sequence.append("certificates/local_certificate.py")
        sequence += ["certificates/verify_local.py", "proof_checks/universal_witness.py",
                     "build_certificate_table.py"]
        for script in sequence:
            report["steps"].append(run_step(work, script))
        for relative_path in OUTPUTS:
            comparison = compare_output(work, relative_path)
            if args.mode == "verify" and relative_path == OUTPUTS[1]:
                comparison["role"] = "supplied candidate input; scalar bounds replayed by verify_local.py; search not regenerated"
            report["comparisons"].append(comparison)
        report["summary"] = summarize(work)
        report["status"] = "PASS"
    except Exception as error:
        report["status"] = "FAIL"
        report["error"] = str(error)
        print("ERROR: " + str(error), file=sys.stderr)
    finally:
        report["elapsed_seconds"] = round(time.perf_counter() - started, 6)
        write_json(work / "report.json", report)
    print(report["status"] + ": " + str(work / "report.json"), flush=True)
    if report["status"] == "PASS":
        print("5 certificates; 300 elementary interval cases; 20 tree witnesses; printed table reproduced.")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
