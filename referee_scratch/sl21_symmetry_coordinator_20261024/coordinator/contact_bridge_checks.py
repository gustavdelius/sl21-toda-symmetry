"""Coordinator's independent scalar/field-jet reconstruction of the final bridge."""
import sympy as s

t, h, v, vb = s.symbols('t h v vb', nonzero=True)
c, e, qw, alpha = s.symbols('c e qw alpha')
eta, beta = h-t, t-1
d = -(c+e)/t

def check(name, expression):
    assert s.factor(expression) == 0, (name, s.factor(expression))
    print('PASS '+name)

# These scalar simple coefficients come from expanding the rational current
# factor. Derivatives of its collapsed h-kernel have exponent h-2.
check('node2 scalar Laurent violation', (2-t)-(-(h-2))-eta)
check('node0 scalar Laurent violation', (t-2)-(h-2)+eta)

# Keep the full antichiral field jet and the remaining source contraction.
kbar = e-(h-1)/vb
db = c-1/vb
qbar = d+e-(h-2)/vb
check('node2 complete antichiral first jet', beta*kbar-db-t*qbar-eta/vb)

# Divide all contacts by 2 i nu lambda0 lambda2 Q. The strip constants
# were integrated independently in intake_checks.py; retain qw as a full jet.
action = beta*(t*qw+eta/v-db+beta*kbar)
tail1 = -t*(qw+eta/v)
tail2 = t*(db-beta*(qw+kbar))
check('node2 action including strip edge', action-beta*(t*(qw+qbar)+eta*(1/v+1/vb)))
check('both tails with all field jets', tail1+tail2+t*(t*(qw+qbar)+eta*(1/v+1/vb)))
check('raw mixed contact including regulator', action+tail1+tail2+t*(qw+qbar)+eta*(1/v+1/vb))

db0 = c-beta/vb
kbar0 = c-(h-1)/vb
check('node0 complete leftover', eta/v+db0-kbar0-eta*(1/v+1/vb))

raw = action+tail1+tail2
matched = raw-alpha*(qw+qbar+eta*(1/v+1/vb))
check('insertionwise finite matching', matched+(t+alpha)*(qw+qbar)+(1+alpha)*eta*(1/v+1/vb))
check('canonical physical endpoint coefficient', (t+alpha).subs(alpha,t*(t-1))-t*t)
whole_zero_form = raw-alpha*(qw+qbar)
check('off-diagonal convention difference', matched-whole_zero_form+alpha*eta*(1/v+1/vb))
print('INFO This checks separated-source coefficients, not continuation or a quantum domain.')
