"""Export the completed022 summary without changing or selecting results."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
summary = json.loads((ROOT / 'results/naibbe-unknown/summary.json').read_text())
rows = summary['summary']
assert len(rows) == 48 and sum(r['primary'] for r in rows) == 24
assert sum(r['primary'] and r['gate_pass'] for r in rows) == summary['primary_gates_passed']
with (ROOT / 'results/naibbe-unknown-summary.csv').open('x', newline='') as output:
    writer = csv.DictWriter(output, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps({'status': 'EXPORTED', 'rows': len(rows)}))
