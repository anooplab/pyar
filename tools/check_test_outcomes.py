"""Reject unreviewed xfails and unexpected skipped-test causes in CI JUnit."""

import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET


def check(report, inventory):
    approved = json.loads(Path(inventory).read_text())
    skipped = []
    executed = 0
    for case in ET.parse(report).iter("testcase"):
        executed += 1
        outcome = case.find("skipped")
        if outcome is None:
            continue
        key = case.attrib["classname"] + "::" + case.attrib["name"]
        reason = outcome.attrib.get("message", "")
        if outcome.attrib.get("type") == "pytest.xfail":
            if approved.get(key) != reason:
                raise ValueError(f"Unreviewed xfail: {key}: {reason}")
        else:
            # Optional tests are allowed only at reviewed locations. Reasons and
            # counts remain visible in the retained JUnit report.
            module = case.attrib["classname"].split(".")[1]
            allowed = json.loads(Path("tests/optional_test_modules.json").read_text())
            if module not in allowed:
                raise ValueError(f"Unreviewed skip: {key}: {reason}")
            skipped.append({"test": key, "reason": reason})
    if not executed:
        raise ValueError("No tests executed")
    print(json.dumps({"tests": executed, "optional_skips": skipped}, indent=2))


if __name__ == "__main__":
    check(sys.argv[1], sys.argv[2])
