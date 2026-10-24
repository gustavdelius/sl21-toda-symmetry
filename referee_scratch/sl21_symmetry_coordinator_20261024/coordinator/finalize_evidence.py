"""Read-only preservation audit; writes only this investigation's new records."""
from pathlib import Path
import datetime
import hashlib
import json
import re
import subprocess

root = Path('/home/gustav/Git/affine_toda')
scratch = root/'referee_scratch/sl21_symmetry_coordinator_20261024'
report_names = [
    'referee_sl21_symmetry_coordinator_report.md',
    'referee_sl21_coord_stage0_audit_report.md',
    'referee_sl21_coord_stage1_subtractions_report.md',
    'referee_sl21_coord_stage1_kernels_report.md',
    'referee_sl21_coord_stage2_companion_report.md',
    'referee_sl21_coord_stage1_contact_bridge_report.md',
    'referee_sl21_coord_final_adversarial_report.md',
]
reports = [root/'briefs'/name for name in report_names]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def git(*args):
    return subprocess.check_output(['git', *args], cwd=root, text=True)

baseline = json.loads((scratch/'preexisting_hashes.json').read_text())
external = json.loads((scratch/'coordinator/external_source_hash.json').read_text())
missing, changed = [], []
for name, expected in {**baseline, **external}.items():
    path = root/name
    if not path.is_file():
        missing.append(name)
    elif digest(path) != expected:
        changed.append(name)

initial = json.loads((scratch/'initial_state.json').read_text())
record = {
    'date_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'preexisting_repository_files_checked': len(baseline),
    'external_source_files_checked': len(external),
    'missing': missing,
    'changed': changed,
    'commit': git('rev-parse','HEAD').strip(),
    'branch': git('branch','--show-current').strip(),
    'tracked_diff': git('diff','--stat'),
    'staged_diff': git('diff','--cached','--stat'),
    'status': git('status','--porcelain=v1'),
    'fresh_reports': [str(p.relative_to(root)) for p in reports],
    'scope': 'Fresh reports and scratch only; no book integration, git mutation or render.',
}
assert not missing and not changed, record
assert record['commit'] == initial['commit'].strip()
assert record['branch'] == initial['branch'].strip()
assert not record['tracked_diff'] and not record['staged_diff']
(scratch/'final_preservation.json').write_text(json.dumps(record, indent=2)+'\n')

manifest_path = scratch/'output_manifest.json'
files = sorted([p for p in scratch.rglob('*') if p.is_file() and p != manifest_path] + reports)
manifest = {str(p.relative_to(root)): digest(p) for p in files}
manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')

# Local evidence links in the synthesis must resolve, including final records.
main_report = reports[0]
broken = []
for target in re.findall(r'\]\(([^)]+)\)', main_report.read_text()):
    if '://' in target or target.startswith('#'):
        continue
    if not (main_report.parent/target.split('#')[0]).is_file():
        broken.append(target)
assert not broken, broken

print(f'PASS {len(baseline)} pre-existing repository files unchanged')
print(f'PASS {len(external)} external primary source unchanged')
print('PASS commit and branch unchanged; tracked and staged diffs empty')
print(f'PASS all local synthesis links resolve; {len(reports)} fresh reports')
print(f'PASS output manifest records {len(manifest)} fresh files, excluding itself')
