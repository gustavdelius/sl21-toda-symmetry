import sympy as S
x,y,b=S.symbols('x y b')
for name,expr in [
 ('two-tail antichiral vertex splits into single edges',(x-y)/(x*y)-(1/y-1/x)),
 ('two-tail barC contraction splits into single double poles',-b*(x-y)/(x*y)*(1/x+1/y)-b*(1/x**2-1/y**2))]:
 print(('PASS ' if S.factor(expr)==0 else 'FAIL ')+name)
