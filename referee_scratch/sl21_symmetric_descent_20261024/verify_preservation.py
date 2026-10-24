"""Read-only comparison against this referee pass's initial snapshot."""
import hashlib,json,subprocess
from pathlib import Path
base=Path(__file__).resolve().parent
root=base.parents[1]
baseline=json.loads((base/'preexisting_hashes.json').read_text())
initial=json.loads((base/'initial_state.json').read_text())
changed=[]
for name,digest in baseline.items():
    path=root/name
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
        changed.append(name)
def git(*args):
    return subprocess.check_output(['git',*args],cwd=root,text=True)
state={'commit':git('rev-parse','HEAD'),'branch':git('branch','--show-current'),
       'diff':git('diff'),'staged':git('diff','--cached')}
result={'baseline_file_count':len(baseline),'changed_preexisting_files':changed,
        'git_state_matches_initial':state==initial,'current_state':state}
(base/'preservation.json').write_text(json.dumps(result,indent=2)+'\n')
print(('PASS' if not changed else 'FAIL')+' preservation of '+str(len(baseline))+' pre-existing files')
print(('PASS' if state==initial else 'FAIL')+' HEAD, branch, tracked diff and index match initial state')
if changed or state!=initial:raise SystemExit(1)
