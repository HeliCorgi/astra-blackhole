"""Compile physics-claim labels from machine-readable audit obligations.

This checker does not determine which quantization is physically correct.
It validates declared evidence and prevents claim promotion around unresolved
or failed blocking gates.

Copyright 2026 HeliCorgi
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[1]
ALLOWED_STATUS = {"PASS", "FAIL", "PARTIAL", "PENDING", "N/A"}
ALLOWED_GATES = {
    "GR_REDUCTION", "LEAN_SEMANTIC", "CLOCK", "INNER_PRODUCT_DOMAIN",
    "FACTOR_ORDERING", "REGULATOR", "KNOWN_LIMIT", "NUMERICAL",
    "CLAIM_COMPILER",
}

def dotted_get(obj, path: str):
    cur = obj
    for part in path.split("."):
        if isinstance(cur, list):
            cur = cur[int(part)]
        else:
            cur = cur[part]
    return cur

def validate_evidence(item: dict) -> None:
    path = ROOT / item["path"]
    if not path.exists():
        raise ValueError(f"Missing evidence path: {item['path']}")
    expected = item.get("expect", {})
    if expected:
        if path.suffix != ".json":
            raise ValueError(f"JSON expectations require JSON evidence: {item['path']}")
        data = json.loads(path.read_text())
        for key, value in expected.items():
            actual = dotted_get(data, key)
            if actual != value:
                raise ValueError(
                    f"Evidence mismatch {item['path']}:{key}: expected {value!r}, got {actual!r}"
                )

def load_and_validate(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != 1:
        raise ValueError("Unsupported physics audit schema")
    obligations = data.get("obligations")
    if not isinstance(obligations, list) or not obligations:
        raise ValueError("No obligations declared")
    ids = [item.get("id") for item in obligations]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate obligation id")
    for item in obligations:
        if item.get("status") not in ALLOWED_STATUS:
            raise ValueError(f"Invalid status for {item.get('id')}: {item.get('status')}")
        if item.get("gate") not in ALLOWED_GATES:
            raise ValueError(f"Invalid gate for {item.get('id')}: {item.get('gate')}")
        evidence = item.get("evidence", [])
        if item["status"] in {"PASS", "FAIL", "PARTIAL"} and not evidence:
            raise ValueError(f"{item['id']} needs evidence for status {item['status']}")
        for ev in evidence:
            validate_evidence(ev)
    return data

def aggregate_gate_status(obligations: list[dict]) -> dict[str, str]:
    grouped: dict[str, list[str]] = {}
    for item in obligations:
        grouped.setdefault(item["gate"], []).append(item["status"])
    result = {}
    for gate, statuses in grouped.items():
        s = set(statuses)
        if "FAIL" in s:
            result[gate] = "FAIL"
        elif "PENDING" in s:
            result[gate] = "PENDING"
        elif "PARTIAL" in s:
            result[gate] = "PARTIAL"
        elif s == {"N/A"}:
            result[gate] = "N/A"
        else:
            result[gate] = "PASS"
    return result

def compile_claims(data: dict) -> dict:
    obligations = data["obligations"]
    labels: list[str] = []
    lean = [x for x in obligations if x["gate"] == "LEAN_SEMANTIC"]
    if lean and all(x["status"] == "PASS" for x in lean):
        labels.append("ALGEBRAICALLY VERIFIED")
    if any(x["gate"] == "NUMERICAL" and x["status"] == "PASS" for x in obligations):
        labels.append("MODEL-INTERNAL NUMERICAL RESULT")
    if any(x["gate"] == "CLOCK" and x["status"] == "FAIL" and x.get("finding") == "dependent" for x in obligations):
        labels.append("CLOCK-DEPENDENT")
    if any(x["gate"] == "FACTOR_ORDERING" and x["status"] == "FAIL" and x.get("finding") == "sensitive" for x in obligations):
        labels.append("ORDERING-SENSITIVE")
    if any(x["gate"] == "REGULATOR" and x["status"] == "FAIL" for x in obligations):
        labels.append("REGULATOR-UNSTABLE")
    semi = [x for x in obligations if x["id"] == "known_limits.semiclassical"]
    if semi and semi[0]["status"] == "PASS":
        labels.append("SEMICLASSICAL LIMIT PASSED")
    blockers = [x["id"] for x in obligations if x.get("blocking", False) and x["status"] not in {"PASS", "N/A"}]
    if blockers:
        labels.append("PHYSICAL INTERPRETATION NOT IDENTIFIED")
    order = [
        "ALGEBRAICALLY VERIFIED", "MODEL-INTERNAL NUMERICAL RESULT",
        "CLOCK-DEPENDENT", "ORDERING-SENSITIVE", "REGULATOR-UNSTABLE",
        "SEMICLASSICAL LIMIT PASSED", "PHYSICAL INTERPRETATION NOT IDENTIFIED",
    ]
    labels = [label for label in order if label in labels]
    return {
        "schema": 1,
        "target": data["target"],
        "labels": labels,
        "gate_status": aggregate_gate_status(obligations),
        "promotion_blockers": blockers,
        "failed_obligations": [x["id"] for x in obligations if x["status"] == "FAIL"],
        "partial_obligations": [x["id"] for x in obligations if x["status"] == "PARTIAL"],
        "pending_obligations": [x["id"] for x in obligations if x["status"] == "PENDING"],
        "interpretation_rule": "Labels describe only passed scoped gates. A blocking PENDING/PARTIAL/FAIL forces PHYSICAL INTERPRETATION NOT IDENTIFIED.",
    }

def self_test() -> None:
    base = {
        "target": {"id": "synthetic"},
        "obligations": [
            {"id": "l1", "gate": "LEAN_SEMANTIC", "status": "PASS", "blocking": False},
            {"id": "n1", "gate": "NUMERICAL", "status": "PASS", "blocking": False},
            {"id": "r1", "gate": "REGULATOR", "status": "FAIL", "blocking": True},
            {"id": "s1", "gate": "KNOWN_LIMIT", "status": "PENDING", "blocking": True},
        ],
    }
    claims = compile_claims(base)
    expected = {
        "ALGEBRAICALLY VERIFIED", "MODEL-INTERNAL NUMERICAL RESULT",
        "REGULATOR-UNSTABLE", "PHYSICAL INTERPRETATION NOT IDENTIFIED",
    }
    if set(claims["labels"]) != expected:
        raise AssertionError(claims)
    base["obligations"][2]["status"] = "PASS"
    base["obligations"][3]["status"] = "PASS"
    claims = compile_claims(base)
    if "PHYSICAL INTERPRETATION NOT IDENTIFIED" in claims["labels"]:
        raise AssertionError("cleared blockers still emitted physical-interpretation block")
    with tempfile.TemporaryDirectory() as directory:
        bad = Path(directory) / "bad.json"
        bad.write_text(json.dumps({"schema": 1, "obligations": []}))
        try:
            load_and_validate(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("empty obligations unexpectedly validated")
    print("Physics audit compiler self-test passed.")

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--obligations", type=Path)
    parser.add_argument("--claims-out", type=Path)
    parser.add_argument("--check-against", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    if args.obligations is None:
        if args.self_test:
            return
        parser.error("--obligations is required unless only --self-test is used")
    data = load_and_validate(args.obligations)
    claims = compile_claims(data)
    text = json.dumps(claims, ensure_ascii=False, indent=2) + "\n"
    if args.claims_out:
        args.claims_out.parent.mkdir(parents=True, exist_ok=True)
        args.claims_out.write_text(text)
    else:
        print(text, end="")
    if args.check_against:
        expected = json.loads(args.check_against.read_text())
        if claims != expected:
            raise RuntimeError("physics-audit claim compilation differs from saved baseline")

if __name__ == "__main__":
    main()
