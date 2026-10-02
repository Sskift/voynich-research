"""Download known public mirror images, without altering their pixels."""
import concurrent.futures
import datetime
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAMES = ['088r','088v','089r','089v','099r','099v','100r','100v','101r','101v','102r','102v']


def fetch(name):
    url = f'https://www.voynich.com/folios/color/{name}.jpg'
    path = ROOT/'data/images'/f'{name}.jpg'
    if not path.exists():
        subprocess.run(['curl','-fLsS','--max-time','60',url,'-o',str(path)],check=True)
    data = path.read_bytes()
    assert data[:2] == b'\xff\xd8', (name, 'not JPEG')
    result = dict(file=str(path.relative_to(ROOT)),url=url,bytes=len(data),
                  sha256=hashlib.sha256(data).hexdigest(),
                  retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  source='Public mirror, linked by https://www.voynich.com/folios/',
                  note='Mirror filename is NOT an IVTFF folio ID; consult pharma-image-coverage.json for pharmaceutical scans.')
    print(name, len(data), flush=True)
    return result


if __name__ == '__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(fetch, NAMES))
    (ROOT/'data/pharma-image-sources.json').write_text(json.dumps(records,indent=2)+'\n')
