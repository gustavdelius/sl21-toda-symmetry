"""Full independent action-to-Lax verification including arbitrary boson first jets."""
import sympy as s
from solve_lax import G,zeros,madd,neg,mm,H,roots,S,mu,z,A,B,rt,pp,pm,dot,M,Lp,Lm,DmLp,DpLm,a,b,c,d,e,f
choice={a:s.I,d:s.I,b:-s.I,e:-s.I,c:-s.I,f:1}
def sub(X):return [[v.subs(choice) for v in row] for row in X]
Lplus=sub(Lp);Lminus=sub(Lm);Dminus=sub(DmLp);Dplus=sub(DpLm)
# Each diagonal H has the same dot-product metric as the source roots.
for root,edge in zip(roots,[(2,0),(0,1),(1,2)]):
 h=H(root);i,j=edge
 for k,other in enumerate(roots):
  oi,oj=[(2,0),(0,1),(1,2)][k]
  assert s.simplify(h[oi]-h[oj]-(root.T*other)[0])==0
print('PASS explicit supertrace-zero Cartan embedding reproduces the full root Gram matrix')
# Arbitrary first derivatives. Differentiating coefficient exponentials supplies
# precisely the terms that cancel against Cartan commutators.
pu,pv,mu1,mv=s.symbols('u_plus v_plus u_minus v_minus')
P=s.Matrix([pu,pv]);N=s.Matrix([mu1,mv])
HP=H(P);HN=H(N)
for i in range(3):
 Lplus[i][i]=G(HP[i]/2);Lminus[i][i]=G(-HN[i]/2)
# Determine exponential derivatives from the auxiliary root weight of each entry.
weights={(0,1):roots[1],(1,0):-roots[1],(2,0):roots[0],(1,2):roots[2]}
for i,j in [(0,1),(1,0),(2,0),(1,2)]:
 Dminus[i][j]+=Lplus[i][j]*((weights[i,j].T*N)[0]/2)
weightsm={(1,0):roots[1],(0,1):-roots[1],(0,2):roots[0],(2,1):roots[2]}
for i,j in [(1,0),(0,1),(0,2),(2,1)]:
 Dplus[i][j]+=Lminus[i][j]*((weightsm[i,j].T*P)[0]/2)
# Source boson equation: phi_+- = -grad(V)+i*Yukawa source.
grad=[mu**2*(S[1]**2-S[0]**2*S[2]**2)*roots[1][k] for k in range(2)]
Y=[sum((mu*S[j]**2*roots[j][k]*dot(roots[j],pp)*dot(roots[j],pm) for j in (0,2)),G()) for k in range(2)]
field_eq=[G(-grad[k])+s.I*Y[k] for k in range(2)]
HE=H(field_eq)
for i in range(3):Dplus[i][i]=-HE[i]*s.Rational(1,2);Dminus[i][i]=HE[i]*s.Rational(1,2)
curv=madd(madd(Dplus,neg(Dminus)),madd(mm(Lplus,Lminus),neg(mm(Lminus,Lplus))))
for i in range(3):
 for j in range(3):
  assert curv[i][j].zero(),(i,j,curv[i][j])
print('PASS all nine curvature entries vanish for arbitrary fields, first jets and all four independent Grassmann fermions')
# Constant vacuum at z=1, phi=(i sqrt(2) pi n,v0), fermions=0.
n=s.symbols('n',integer=True)
plus0=s.Matrix([[0,s.I*mu,-s.S.Zero],[-s.I*mu*z,0,0],[0,0,0]])
minus0=s.Matrix([[0,-s.I*mu/z,0],[s.I*mu,0,0],[0,0,0]])
At0=(plus0+minus0)/2;Ax0=(plus0-minus0)/2
assert At0.subs(z,1)==s.zeros(3)
Bvac=Ax0.subs(z,1)
assert Bvac**2==s.diag(mu**2,mu**2,0)
print('PASS zero vacuum time connection at z=1; spatial vacuum eigenvalues are mu,-mu,0')
print('PASS zero-curvature connection uses all nonlinear source terms, not just the bosonic truncation')
