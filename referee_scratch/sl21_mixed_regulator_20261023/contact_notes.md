# Working derivation checkpoint

This is scratch, not a completed verdict. Existing reports remain unchanged.

## Coordinate correction

Working z=x+i tau, t_M=-i tau gives z=x-t_M. Thus phi=(Phi-integral dotPhi)/2, not the plus split of the mixed-charges report. Source BL z=(tau+i x)/2, and working z=2i source bar-z. Correct [phi,barphi]=+pi i. Put nu=exp(-i pi/t); R diagonal nu^2, off-diagonal nu^-2. The compatible cocycle is the previous local sign times exp[-i pi/t(n_L m_R-n_R m_L)]. The rephasing f is also inverted. Local action/current/endpoint counterterms are unchanged. Topological weights for the previous fixed Htop expression reverse.

## Direct tail-current contacts

For B, let v=d+a2, q=t-1, K_L=v.partial X, K_R=v.barpartial X, C=a0.partial X, Cbar=a0.barpartial X. Factor nu lambda2 out of the two products. J bar-H has double-pole coefficient tq U and simple-pole coefficient -t partial U. H bar-J has the same double-pole coefficient tq U, and simple coefficient tq partial U+t^2 barpartial U (the holomorphic Taylor expansion of H is essential). The holomorphic Plemelj sign is negative, the antiholomorphic sign positive. Their combined density contact is 2pi i nu lambda2 t^2 partial_x U. AA similarly gives nu gamma^2 lambda1 partial_x F_A. These are the tail sectors; they do NOT yet include the J bar-J action-insertion contact.

## A finite current-product contact, independent of t convergence

Integrate the outer J_B current over x at tau=+delta minus tau=-delta. Normalize its simple and double poles as Cauchy distributions. Put the inner bar-J at the origin. With Theta_delta=1_{|Im w|<delta}, its current/Y2 kernel produces

    -2pi i exp(D w) [(C+qD)Theta_delta+q partial_w Theta_delta].

Here D=d.partial X and C+qD=-K_L. An interaction disk |w|<epsilon around the inner current is excluded. With c=delta/epsilon, the three elementary face integrals are

    I(c)  = integral_{|w|>epsilon, |Im w|<delta} d2w / bar-w^2,
    J1(c) = integral_{|w|>epsilon} d2w (partial_w Theta_delta)/bar-w,
    Jw(c) = integral_{|w|>epsilon} d2w (partial_w Theta_delta) w/bar-w^2.

For 0<c<1:

    I=2(arcsin c+c sqrt(1-c^2)),
    J1=-2 arcsin c,
    Jw=-2 arcsin c-4c sqrt(1-c^2).

For c>=1 they are respectively pi,-pi,-pi. The limit c->0 gives zero for all three. The derivative-of-strip term must be retained; omitting it gives a noncovariant incorrect answer even for c>=1.

Let A_R=a2.barpartial X, so q A_R-Cbar=t K_R. Expanding the fused fields to the order that survives epsilon->0 gives the bare-J/action-sector density contact

    -2i nu lambda2 q {K_L[I-q Jw]-t K_R J1} U.

For c>=1 this is -2pi i nu lambda2 t(t-1) partial_x U. Added to the tail-only contact, the unadjusted radial/Cauchy prescription gives coefficient t, not t^2. The analogous A calculation gives gamma^2/t rather than gamma^2. These are finite product-normalization differences, not an established correction to Bernard-LeClair. Need numerical checks of I,J1,Jw and careful report scope.

Canonical matching: impose the independent equal-time product condition [J_B,bar-J_B]_R=0, as realized by the anchored canonical fields and their spatial derivatives. This requires a finite product contact cancelling the action-sector expression. For c>=1 it adds +2pi i nu lambda2 t(t-1) partial_x U. It is an additional mixed-product normalization, not the action counterterm or a current-insertion counterterm. Need distinguish a defined formal product scheme from an established physical operator realization.

## Mass-squared forest / analytic continuation proposal

In the lambda0 background, the same collision subtraction acts on the descent chain:

    [Y2]_0 = Y2-lambda0 epsilon^(2t) P/t,
    [K2]_0 = K2-lambda0 epsilon^(2t) L/t,
    [bar-K2]_0 = bar-K2-lambda0 epsilon^(2t) bar-L/t,
    [U]_0 = U-lambda0 epsilon^(2t-2)F_B/(t-1)
               -lambda0 epsilon^(2t) C Cbar F_B/t.

The first expression produces exactly the action composite counterterm upon multiplication by lambda2. The next two produce exactly the previously verified current tails. The last gives the endpoint subtraction. J_B has zero fixed-current conservation coefficient with Y0 (double minus simple coefficients cancel). Hence after current-product contact normalization one can try to commute the contour descent past the Y0 insertion. The potential obstruction is the simultaneous collision boundary, not the separated OPE.

All primitive mass-squared joint corners after lower-order current-pole normalization vanish in Re t>1. Fresh kernel:

    K02=nu |w0-w2|^(2t-2) f(u;w0,w2) f(bar-v;bar-w0,bar-w2),
    f(u;w0,w2)=[1/(u-w0)+(t-1)/(u-w2)]/[(u-w0)(u-w2)].

Its kernel degree is 2t-8; two contour and two area measures give 2t-2. J bar-K Y0 and action-P contacts including epsilon^(2t) have the same worst power. A current-L contact has power 2t-1 (one derivative remains uncontracted). Single-current/two-Y and Y0/Y2 subfaces have power 2t. Same-node joint corners and K2/bar-K2 contacts have power 2. These facts are independently checked in regulator_checks.py. They do NOT alone prove cancellation of finite current-pole faces; the first-order contact above is a counterexample to using naive power counting alone.

Possible completion: define the mixed product independently by its Gaussian expansion, subtract all one-current pole faces with the canonical matching condition just stated, prove by contour deformation that the remaining second-order term is the one-Y0 dressing of the leading endpoint in Re t>1, and continue the identity meromorphically with no extra finite terms. Then the formal canonically matched bracket would be nu/(2pi i) Delta[t^2 lambda2 U]_ren through m^2. Must explain the subtraction of the finite first-order product contact and its dressing; otherwise this is not established. Do not claim an affine algebra or a global Hilbert-space charge.

## Files and verification so far

- New scratch regulator_checks.py and regulator_checks.out: 52/52 checks passed before the finite-face tests are added.
- Source PDF pages 105 and 106 were rendered and inspected in this directory.
- inputs_before.sha256 records AGENTS.md and all six preceding reports.
- Proposed new report: briefs/referee_sl21_mixed_regulator_report.md (not yet created).
