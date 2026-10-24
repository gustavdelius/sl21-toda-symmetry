import sympy as S
sig,q,C,Cb,psi0,psiz,psiw,psib=S.symbols('sig q C Cb psi0 psiz psiw psib');I=S.I
# Independently integrate: outer 1/(z-w), P, zP; then sign/barw²,
# edge w/barw², edge 1/barw. No field-jet relation is imposed here.
Msign=S.pi*sig;Mw=I*S.pi*sig;Mone=I*S.pi*sig
before=q**2*psiz*(-I*S.pi*Msign+S.pi*Mw)+q**2*S.pi*(psiw*Mw+psib*Mone)-I*S.pi*q*C*psi0*Msign-S.pi*q*Cb*psi0*Mone
raw=-before/S.pi
expected=-I*S.pi*sig*(q**2*(psiw+psib)-q*(C+Cb)*psi0)
checks=[('full independent pair jet',S.expand(raw-expected)),('outer test derivative',S.expand(raw).coeff(psiz)),('node0 switching contact survives',S.expand(raw.subs({q:1,psi0:0,psiw:1,psib:0}))+I*S.pi*sig),('node2 switching coefficient',S.expand(raw.subs({q:S.Symbol('t')-1,psi0:0,psiw:1,psib:0}))+I*S.pi*sig*(S.Symbol('t')-1)**2)]
for n,e in checks:
 e=S.simplify(e);print(('PASS' if e==0 else 'FAIL')+' '+n+' residual='+str(e));assert e==0
wrong=S.simplify(raw.subs({q:S.Symbol('t')-1,psi0:0,psiw:1,psib:0})-0)
print('EXPECTED NONZERO constant-spectator shortcut residual='+str(wrong))
print('TOTAL',len(checks),'PASS; one explicit counterexample retained')
