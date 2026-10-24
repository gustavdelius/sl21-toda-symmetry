import itertools
import sympy as S
checks=[]
def test(name,b):
    b=bool(b);checks.append(b);print(('PASS' if b else 'FAIL')+' '+name)
# Separate labels from physical momentum. All fusion additions use both, never net charge alone.
vertices={'A':(-1,0,0,1,0),'B':(1,0,0,0,1),'K':(1,0,1,0,1),'F':(-1,-1,-1,1,0),'L':(1,1,1,0,1)}
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def cocycle_exp(x,y):return x[0]*(y[1]+y[2])
def like(x,y):return x[3]*y[3]+x[4]*y[4]
for names in itertools.product(vertices,repeat=3):
    x,y,z=[vertices[n] for n in names]
    test('cocycle and label factorization '+''.join(names),
      (cocycle_exp(x,y)+cocycle_exp(add(x,y),z)-cocycle_exp(y,z)-cocycle_exp(x,add(y,z)))%2==0
      and like(x,y)+like(add(x,y),z)==like(y,z)+like(x,add(y,z)))
# Every matching of four derivative legs is recovered by cross contractions along any binary tree.
def matchings(ids):
    if not ids:return [frozenset()]
    i,*rest=ids
    ans=matchings(rest)
    for j in rest:
        ans += [m|{tuple(sorted((i,j)))} for m in matchings([k for k in rest if k!=j])]
    return ans
all_match=set(matchings([0,1,2,3]))
def leaves(tree):return [tree] if isinstance(tree,int) else leaves(tree[0])+leaves(tree[1])
def tree_match(tree):
    if isinstance(tree,int):return {frozenset()}
    left,right=tree;L=leaves(left);R=leaves(right);ans=set()
    for lm in tree_match(left):
      for rm in tree_match(right):
        used={i for pair in lm|rm for i in pair}
        lf=[i for i in L if i not in used];rf=[i for i in R if i not in used]
        for n in range(min(len(lf),len(rf))+1):
          for ll in itertools.combinations(lf,n):
            for rr in itertools.permutations(rf,n):
              ans.add(lm|rm|frozenset(tuple(sorted(p)) for p in zip(ll,rr)))
    return ans
for tree in [(((0,1),2),3),((0,(1,2)),3),((0,1),(2,3)),(0,((1,2),3)),(0,(1,(2,3)))]:
    test('Wick matching factorization '+str(tree),tree_match(tree)==all_match)
# Regular normal-field derivative terms correspond to unmatched legs and remain arbitrary.
test('all 10 derivative matchings retained',len(all_match)==10)
# Disjoint B slot tests retain only one labelled B permutation in a middle A chamber.
s=S.symbols('s',nonzero=True);c=1+s+1/s
restricted=sum(coef*s**(-abs(pos-2)) for pos,coef in zip([3,2,1,0],[1,-c,c,-1]))
test('distinct B slot tests do not obey simple Serre polynomial',S.factor(restricted-(1/s-s))==0)
print('DIAGNOSTIC one-permutation middle-chamber coefficient =',S.factor(restricted))
# Tail-only lower pairs differ from action-dressed old ones, even though quartics agree.
t=S.symbols('t');eta=t-2;gamma=t/(1-t)
test('BB tail-only pair versus action-dressed pair',S.expand(2*t*eta-2*(t-1)*eta-2*eta)==0)
test('AA tail-only pair versus action-dressed pair',S.factor(-2*gamma-2+2/(1-t))==0)
print('SUMMARY',sum(checks),'PASS',len(checks)-sum(checks),'FAIL')
if not all(checks):raise SystemExit(1)
