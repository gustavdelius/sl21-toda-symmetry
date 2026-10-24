"""Run a small completed recent suite in a disposable copy of its evidence tree."""
from pathlib import Path
import os,shutil,subprocess,sys,tempfile
root=Path(__file__).resolve().parents[1]
relative=Path('referee_scratch/sl21_two_barred_20261024')
scripts={'coordinator/check_two_barred_v2.py':56,'coordinator/check_contact_diamond_v2.py':7,'coordinator/check_balancing.py':5}
failures=[]
with tempfile.TemporaryDirectory(prefix='sl21-checks-') as tmp:
    work=Path(tmp);shutil.copytree(root/relative,work/relative)
    for name,expected in scripts.items():
        run=subprocess.run([sys.executable,'-u',str(work/relative/name)],cwd=work,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)
        count=sum(line.startswith('PASS ') for line in run.stdout.splitlines())
        ok=run.returncode==0 and count==expected and not any(line.startswith('FAIL ') for line in run.stdout.splitlines())
        print(f"{'PASS' if ok else 'FAIL'} {name}: {count}/{expected} checks",flush=True)
        if not ok:failures.append(name);print(run.stdout)
print('This smoke test checks selected formulas, not every historical script or the analytic proofs.')
sys.exit(bool(failures))
