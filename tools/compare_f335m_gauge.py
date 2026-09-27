"""Compare the new gauge JSON against the committed ad hoc artifact.

Kept in the repo as a script rather than a test because it needs the real
40 GB cal set; run it after changing gauge.py to prove the subcommand still
reproduces the recorded F335M numbers.

The reference is read out of git, so this keeps working after the ad hoc
artifact has been replaced by the subcommand's output.
"""

import json
import subprocess
import sys
from pathlib import Path

FIELDS = ("n_matched", "dx_median", "dy_median", "offset_median", "offset_mad")
REFERENCE = "out/f335m_step2_detectors.json"


def reference() -> dict:
    """The pre-subcommand artifact, from git or from a local copy."""
    local = Path("out/f335m_step2_detectors.committed.json")
    if local.exists():
        return json.loads(local.read_text(encoding="utf-8"))
    blob = subprocess.run(
        ["git", "show", f"HEAD:{REFERENCE}"],
        capture_output=True, text=True, check=True,
    ).stdout
    return json.loads(blob)


old = reference()
old = old.get("groups", old)
new = json.loads(Path(REFERENCE).read_text(encoding="utf-8"))

if set(old) != set(new):
    print(f"FAIL group sets differ: {set(old) ^ set(new)}")
    sys.exit(1)

worst, worst_field = 0.0, ""
for key in sorted(old):
    for field in FIELDS:
        delta = abs(float(old[key][field]) - float(new[key][field]))
        if delta > worst:
            worst, worst_field = delta, f"{key}.{field}"

print(f"groups compared        : {len(new)}")
print(f"fields per group       : {len(FIELDS)}")
print(f"max |difference|       : {worst}")
print(f"worst field            : {worst_field or 'none (bit-identical)'}")
print(f"match rates (new only) : "
      + ", ".join(f"{k} {v['n_matched']}/{v['n_frame_stars']}"
                  for k, v in list(new.items())[:2]))

sys.exit(0 if worst == 0.0 else 1)
