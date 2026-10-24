"""Read-only preservation check for the two-barred referee pass."""
from pathlib import Path
import json,hashlib,subprocess
base=Path(__file__).resolve().parent;root=base.parents[1]
old=json.loads((base/'preexisting_hashes.json').read_text())
initial=json.loads((base/'initial_state.json').read_text())
hashfile=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
changed=[p for p,h in old.items() if not (root/p).is_file() or hashfile(root/p)!=h]
git=lambda *args:subprocess.check_output(['git',*args],cwd=root,text=True)
state={'commit':git('rev-parse','HEAD'),'branch':git('branch','--show-current'),'diff':git('diff'),'staged':git('diff','--cached')}
result={'baseline_file_count':len(old),'changed_preexisting_files':changed,'git_state_matches_initial':state==initial,'current_state':state}
(base/'preservation.json').write_text(json.dumps(result,indent=2)+'\n')
print(('PASS ' if not changed else 'FAIL ')+'preservation of '+str(len(old))+' pre-existing files')
print(('PASS ' if state==initial else 'FAIL ')+'HEAD, branch, tracked diff and index match initial state')
if changed or state!=initial:raise SystemExit(1)
