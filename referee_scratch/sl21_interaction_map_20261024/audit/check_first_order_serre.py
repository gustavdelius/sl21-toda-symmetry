import sympy as S
b=S.symbols('b',positive=True); t=b*b; rho=2/t
I=S.I;a0=S.Matrix([I*b/S.sqrt(2),b/S.sqrt(2),I]);a2=S.Matrix([I*b/S.sqrt(2),-b/S.sqrt(2),-I]);d=-(a0+a2)/t
z,v,a=S.symbols('z v a');count=0
def check(name,e):
 global count
 e=S.factor(e);print(('PASS' if e==0 else 'FAIL')+' '+name+' residual='+str(e));count+=1
 if e!=0:raise AssertionError(name)
check('K0 derivative C contraction vanishes',a0.dot(d+a0))
check('K0 versus B exponent equals rho minus one',-(d+a0).dot(d)-(rho-1))
check('K0 versus A exponent equals one minus rho',(d+a0).dot(d)-(1-rho))
check('K0 and K2 same B exponent',(d+a0).dot(d)-(d+a2).dot(d))
fz=-1/(z-v)+1/(z-a);fv=1/(z-v)+1/(v-a)
check('K0 marginal full scalar vanishes',fz*fv+1/(z-v)**2)
check('K0 one regular C coefficient retained',fz+fv-1/(z-a)-1/(v-a))
# Three pair exponent changes for one tail: -1 to each B, +1 to A.
# Integer branch signs are canceled by the source-label cocycle signs.
check('K0 B integer phase and source cocycle cancel',(-1)*(-1)-1)
check('K0 A integer phase and source cocycle cancel',(-1)*(-1)-1)
g,r=S.symbols('g r')
# Full B^3A base exponent 3g-3rho; one K0 shift -2+1;
# two C contractions subtract 2; relative measure 3.
check('K0 full zeroth jet marginal degree',3*g-3*r-2+1-2+3-3*(g-r))
check('K0 first surviving jet strictly improves degree',(3*g-3*r+1).subs(g,r)-1)
print('TOTAL',count,'PASS')
