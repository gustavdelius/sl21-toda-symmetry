import sympy as S
from itertools import combinations
l0,l2,S0,S2,S02,S00,S22,beta=S.symbols('l0 l2 S0 S2 S02 S00 S22 beta');ncheck=0
for N in [1,2,3,4]:
 B=S.symbols('B:'+str(N));K0=S.symbols('K0:'+str(N));K2=S.symbols('K2:'+str(N))
 F=S.prod(B);T=S.prod(B[i]+l0*K0[i]+l2*beta*K2[i] for i in range(N))
 E=1+l0*S0+l2*S2+l0*l2*S02+l0**2*S00/2+l2**2*S22/2
 U0=sum(K0[i]*S.prod(B[j] for j in range(N) if j!=i) for i in range(N))
 U2=sum(beta*K2[i]*S.prod(B[j] for j in range(N) if j!=i) for i in range(N))
 U02=sum(K0[i]*beta*K2[j]*S.prod(B[k] for k in range(N) if k not in [i,j]) for i in range(N) for j in range(N) if i!=j)
 U22=sum(beta**2*K2[i]*K2[j]*S.prod(B[k] for k in range(N) if k not in [i,j]) for i in range(N) for j in range(N) if i!=j)
 actual=S.Poly(S.expand(E*T),l0,l2)
 for label,e in [('unlike coefficient',actual.coeff_monomial(l0*l2)-(S02*F+S0*U2+S2*U0+U02)),('equal coefficient',actual.coeff_monomial(l2**2)-(S22*F/2+S2*U2+U22/2))]:
  e=S.expand(e);print(('PASS' if e==0 else 'FAIL')+' N='+str(N)+' '+label+' residual='+str(e));assert e==0;ncheck+=1
 print('INFO N='+str(N)+' distinct source allocations: bulk 1; one unary '+str(2*N)+'; two unary '+str(N*(N-1)))
print('TOTAL',ncheck,'PASS')
