"""Independent finite-family closure checks: leaf identities and unlike faces."""
import itertools
import sympy as s
u=s.symbols('u')
x=s.symbols('x0:3')
beta,C=s.symbols('beta C')

for n in (1,2,3):
    xx=x[:n]
    vand=s.prod(xx[i]-xx[j] for i in range(n) for j in range(i+1,n))
    F=vand/s.prod(xi-u for xi in xx)
    # Partial-fraction coefficient at pole xi, obtained without solving a
    # desired identity. It is polynomial because Vandermonde cancels poles.
    coeff=[s.cancel(vand/s.prod(xj-xi for xj in xx if xj!=xi)) for xi in xx]
    split=sum(ci/(xi-u) for ci,xi in zip(coeff,xx))
    assert s.factor(F-split)==0
    assert all(s.denom(ci)==1 for ci in coeff)
    scalar=F*(C+beta*sum(1/(u-xi) for xi in xx))
    assert s.factor(scalar-(C*split-beta*s.diff(split,u)))==0
    print('PASS',n,'equal-root constituents reduce to polynomial-weighted anti leaves, including derivative current')

t,g,rho,h=s.symbols('t g rho h',positive=True)
# Gram calculations in the d,a0,a2 basis, no imported matrix/table.
G=s.Matrix([[-rho,1,1],[1,-1,1-t],[1,1-t,-1]])
d=s.Matrix([1,0,0]);a2=s.Matrix([0,0,1]);a0=s.Matrix([0,1,0]);T=a0+a2
dot=lambda v,w:(v.T*G*w)[0]
assert s.expand(dot(d+a2,a2))==0
assert s.expand(dot(d+a2,d+a2)-(1-rho))==0
assert s.expand(dot(T,d)-2)==0
print('PASS K2/barK2 oscillator orthogonality; U length; even-root anti power')

# Enumerate all new unlike-source faces of B^3 A barB from the Gram data.
rows=[]
for na in range(2):
 for nb in range(4):
  for nd in range(2):
   N=na+nb+nd
   if not N:continue
   # Like-pair regulator, unlike hol pair, two source contractions,
   # source/source powers, derivative cost and anchored relative measure.
   curr=g*nb*(nb-1)/2-rho*na*nb
   source=2*(na-nb-nd)
   pair=2*h-2
   deriv=nb+nd
   measure=N+4-1
   degree=s.expand(curr+source+pair-deriv+measure)
   expected=curr+3*na-2*(nb+nd)+2*h+1
   assert s.expand(degree-expected)==0
   phys=s.factor(degree.subs({g:2/t,rho:2/t,h:t}))
   # Exact positivity on 2<t<3: expand t*degree at t=2+x. All
   # coefficients are nonnegative, with positive nonconstant or constant.
   q=s.symbols('q',positive=True)
   pol=s.Poly(s.cancel(t*phys).subs(t,q+2).expand(),q)
   assert all(c>=0 for c in pol.all_coeffs()) and any(c>0 for c in pol.all_coeffs())
   rows.append((na,nb,nd,degree,phys))
for na,nb,nd,degree,phys in rows:
 print('INFO unlike face',na,nb,nd,'degree',degree,'physical',phys)
print('PASS all',len(rows),'current-containing unlike-source faces positive for 2<t<3')
assert rows[-1][3]==3*g+2*h-3*rho-4
print('PASS full five-current unlike-source degree 3(g-rho)+2h-4')
print('INFO First-source proper faces and equal-source leaf proof are separate requirements, not consequences of this positivity count.')
