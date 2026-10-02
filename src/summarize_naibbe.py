"""Export already-frozen021 scores; does not generate predictions or alter gates."""
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    freeze=json.loads((ROOT/'results/naibbe-oracle/score-freeze.json').read_text())
    for path,digest in freeze['hashes'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    result=json.loads((ROOT/'results/naibbe-oracle/summary.json').read_text())
    rows=result['summary']
    assert len(rows)==64 and sum(row['primary'] for row in rows)==24
    output=ROOT/'results/naibbe-oracle-summary.csv'
    with output.open('x',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    print(json.dumps({'status':'EXPORTED_FROZEN_SCORES','rows':len(rows),
                      'sha256':hashlib.sha256(output.read_bytes()).hexdigest()}))


if __name__=='__main__':main()
