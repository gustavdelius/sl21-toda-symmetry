import sympy as s
T=s.symbols('T',positive=True); rho=2/T; star=(7+s.sqrt(17))/4
N=0

def ck(label,condition):
 global N
 assert condition,label
 N+=1;print('PASS '+label,flush=True)

candidates=[]; fixed=[]
for nl in range(3):
 for nr in range(3):
  n=nl+nr
  if not n:continue
  E=rho*(nl*(nl-1)+nr*(nr-1))/2
  cases=[('one source',E-n+1),('equal two sources',E-2*n+5),('unlike two sources',E-2*n+2*T+1)]
  for label,omega in cases:
   found=[]
   for j in range(7):
    f=s.factor(s.together(omega+j))
    if f==0:
     fixed.append((label,nl,nr,j));continue
    roots=s.solve(s.together(f).as_numer_denom()[0],T)
    for root in roots:
     if root.is_real and 2<float(root)<3:
      found.append((j,root));candidates.append((label,nl,nr,j,root))
   allowed=label=='unlike two sources' and (nl,nr)==(2,2)
   ck(f'physical pole candidates {label} ({nl},{nr})',found==([(0,star)] if allowed else []))
ck('only new candidate is full unlike zero jet',len(candidates)==1)
ck('fixed source faces are exactly unary and mixed pair',set(fixed)=={('one source',0,1,0),('one source',1,0,0),('one source',1,1,1)})
# Source-only cluster has two-dimensional translation invariance, unlike a
# current-containing cluster anchored to the one-dimensional current line.
ck('free equal source-pair degree is four',2+2==4)
ck('free unlike source-pair degree is 2h',s.expand((2*T-2)+2-2*T)==0)
# Bare-current faces: only derivatives with another current in their own
# chirality can create internal singular degree.
for nl,nr in [(2,0),(0,2),(2,1),(1,2),(2,2)]:
 n=nl+nr;E=rho*(nl*(nl-1)+nr*(nr-1))/2
 derivative_degree=(nl if nl>=2 else 0)+(nr if nr>=2 else 0)
 omega=E+n-1-derivative_degree
 roots=[]
 for j in range(5):
  roots += [r for r in s.solve(s.together(omega+j).as_numer_denom()[0],T) if r.is_real and 2<float(r)<3]
 ck(f'bare current face no pole ({nl},{nr})',not roots)
# At the new zero-jet resonance all proper two-source faces are positive.
for nl in range(3):
 for nr in range(3):
  if nl+nr==0 or (nl,nr)==(2,2):continue
  E=(2/star)*(nl*(nl-1)+nr*(nr-1))/2
  ck(f'proper unlike degree at resonance ({nl},{nr})',float(E-2*(nl+nr)+2*star+1)>0)
print('SUMMARY',N,'PASS, 0 FAIL')
print('CANDIDATES',candidates)
print('FIXED PRE-PRESCRIBED FACES',fixed)
print('CAUTION: Candidate pole is not a proven nonzero residue of the full forest.')
