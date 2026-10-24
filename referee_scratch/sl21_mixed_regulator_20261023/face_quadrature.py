#!/usr/bin/env python3
"""Unexpanded finite-cutoff check of the one-interaction current-pole face.

The inner exclusion disk and the two current-time contours share one cutoff.
The exponent and Gaussian infrared test function are retained in full.
"""
import math
import cmath
from scipy.integrate import quad

t=math.sqrt(2)
q=t-1
c=math.sqrt(2)
C,D,Cbar,Dbar=0.13,-0.17,-0.07,0.09
KL=(1-t)*D-C
KR=(1-t)*Dbar-Cbar
AR=-t*Dbar-Cbar
outer_radius=0.4
expected=math.pi*q*t*(KL+KR)
values=[]
checks=[]

def check(name,ok):
    checks.append(bool(ok))
    print(('PASS ' if ok else 'FAIL ')+name,flush=True)

for eps in [0.005,0.0025,0.00125]:
    delta=c*eps
    # Work in scaled coordinates. Add x and -x before quadrature to remove
    # the odd 1/epsilon piece without subtracting two large final integrals.
    def weight(z):
        return cmath.exp(eps*(KL*z+AR*z.conjugate())
                         -eps**2*abs(z)**2/outer_radius**2)
    def volume_kernel(xx,yy):
        z=complex(xx,yy)
        return weight(z)*(-eps*Cbar/z.conjugate()+q/z.conjugate()**2)
    def slice_volume(yy):
        start=math.sqrt(max(0,1-yy*yy))
        val,err=quad(lambda xx:(volume_kernel(xx,yy)+volume_kernel(-xx,yy)).real,
                     start,math.inf,epsabs=1e-10,epsrel=1e-10,limit=250)
        return val
    # The slice formula changes at y=1, where the exclusion disk ends.
    volume=2*(quad(slice_volume,0,1,epsabs=2e-9,epsrel=2e-9,limit=200)[0]
              +quad(slice_volume,1,c,epsabs=2e-9,epsrel=2e-9,limit=200)[0])
    def edge_kernel(xx):
        z=complex(xx,c)
        return weight(z)*(-Cbar/z.conjugate()+q/(eps*z.conjugate()**2))
    # integral (partial_w Theta_delta) F = - integral Im F(x+i delta) dx
    edge=-quad(lambda xx:(edge_kernel(xx)+edge_kernel(-xx)).imag,
               0,math.inf,epsabs=1e-10,epsrel=1e-10,limit=250)[0]
    result=KL*volume-q*edge
    values.append(result)
    relative_error=abs((result-expected)/expected)
    check('unexpanded finite face eps='+str(eps),relative_error<0.09)
    print('  value='+format(result,'.12g')+' expected='+format(expected,'.12g')
          +' relative_error='+format(relative_error,'.6g'),flush=True)

# The infrared test function gives a linear cutoff correction. Extrapolate
# the two smallest cutoffs without fitting the predicted limit.
extrapolated=2*values[2]-values[1]
check('unexpanded face converges to the calculated finite contact',
      abs((extrapolated-expected)/expected)<0.001)
check('raw errors decrease with the cutoff',
      abs(values[2]-expected)<abs(values[1]-expected)<abs(values[0]-expected))
print(str(sum(checks))+'/'+str(len(checks))+' checks passed',flush=True)
raise SystemExit(0 if all(checks) else 1)
