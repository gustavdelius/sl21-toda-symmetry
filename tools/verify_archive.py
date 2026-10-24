"""Verify imported files without modifying them. Standard library only."""
from pathlib import Path
import hashlib,json,sys
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'provenance/import-manifest.json').read_text())
bad=[]
for row in m['files']:
    p=root/row['path']
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
        bad.append(row['path'])
print(f"{'FAIL' if bad else 'PASS'} {len(m['files'])} imported files checked")
for name in bad:print('CHANGED OR MISSING',name)
sys.exit(bool(bad))
