# Repeated-current mixed overlaps and the first mixed/Serre union

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed commit: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

This is an independent referee calculation for the common-product continuation. No subagents were spawned. No existing book file, report, earlier script, index, or git state was changed. New work is confined to this report and `referee_scratch/sl21_common_product_20261024/overlap/`. The old scripts were not imported or rerun as evidence.

## Result and verdicts

The lowest repeated-current overlap, **BBbarB at order lambda2**, passes its actual local mixed-face test with the canonical finite matching. An independent g-regulator does produce an extra term, and its coefficient must be calculated: it is t(g-rho) in the moving-current simple pole. After the action, both tails, and finite matching are combined, the local contact is the insertion-wise derivative of the canonically normalized endpoint product. Its regulator violation has no pole identically on g=rho in the BBbarB union. This is a finite-order projected overlap result, not a proof of associativity for all products.

For **B^3 A barB at order lambda2**, the proposed full-cluster degree 3(g-rho)-1 is correct. Its first jet is potentially marginal. I found a nonzero residue in one actual decorated mixed-face spatial chamber, so the earlier h-regulator argument cannot be transplanted by merely saying that two prescriptions differ by g-rho. However, that chamber residue is **not a fullword anomaly**. I give an explicit inherited-height, acyclic prescription for the collapsed contact in which the whole-word family is entire in g; consequently all its chamber residues cancel. The unresolved comparison is whether the full source/action/tail forest realizes that prescription, with the same lower products and cocycles.

| Item | Verdict | Independent evidence and scope |
|---|---|---|
| Repeated-B mixed face has the same regulator coefficient as the separated-Yukawa h face | wrong | its simple-pole discrepancy is t(g-rho), not g-rho |
| Actual BBbarB local action/tail/matching contact | correct as (6) below | full spectator derivative-current Laurent expansion and all first jets |
| BBbarB union has a pole identically on g=rho | wrong at this order | degree g-2, Taylor divisors g-2+j, inherited proper faces |
| Full B^3 A barB/Y2 degree is 3(g-rho)-1 | correct | direct Wick-power and integration-dimension count |
| A g-rho factor alone makes decorated finite-matching changes equivalent | wrong as a proof | explicit nonzero single-chamber first-jet residue |
| That residue proves a full five-current anomaly | unclear, and not inferred | fusion choices, chambers and their actual phases must be assembled |
| Acyclic inherited-height contact family has no diagonal pole | correct for the explicit construction in section 5 | proper positive-cone Fourier pushforward, independently checked orientation data |
| Full massive mixed/Serre product already equals that contact family | unclear | action/source forest comparison remains outstanding |

No theorem label or physical symmetry-algebra claim is recommended.

## 1. Definitions and actual repeated-spectator Laurent coefficients

Use X(z)X(w) approximately -log(z-w), t=b^2, rho=2/t,

- a0^2=a2^2=-1, a0.a2=1-t, a0.d=a2.d=1, d^2=-rho;
- C=a0.partial X, B=:C exp(dX):, A=:exp(-dX):;
- beta=t-1, eta=t-2;
- K2=:exp[(d+a2)X+a2 Xbar]:, U=:exp[(d+a2)(X+Xbar)]:;
- H_B=-t lambda2 K2, with the reflected tail on the barred current;
- canonical finite mixed matching alpha lambda2 U, alpha=t(t-1).

Keep a B spectator at a, a moving B at z, and Y2 at w. Multiply the physical BB factor by (z-a)^delta, delta=g-rho. The normal fields and all derivative contractions remain physical. Strip the spectator/barred factors that do not depend on z. The exact holomorphic Wick expression is

\[
\frac{(z-a)^g}{(z-w)(a-w)}
\left\{
\left[C(z)-\frac1{z-a}+\frac\beta{z-w}\right]
\left[C(a)-\frac1{a-z}+\frac\beta{a-w}\right]
+\frac1{(z-a)^2}
\right\}.
\tag{1}
\]

The derivative-pair term in braces is essential. Let K(w) be the actual regulated B(a)K2(w) normal product, including its scalar contractions. With E=(w-a)^g/(a-w), its prefactor is

\[
K/E=C(a)+\frac{\beta-1}{a-w}=:b(w).
\]

Expanding z=w+q, including exp[d(X(z)-X(w))]=1+q d.partial X(w)+..., gives

\[
D_2=\beta K,\qquad
\boxed{S_2=-\partial_wK+\frac{t\delta}{w-a}K.}
\tag{2}
\]

For transparency, put r=w-a and D=d.partial X(w). After stripping E, the simple coefficient from (1) is

\[
C(w)b+\beta Db+\frac{\beta g-1}{r}b+
\frac{1-\beta}{r^2}.
\]

On the other hand

\[
(\partial_wK)/E=(g-1)b/r+(\beta-1)/r^2-b[C(w)+\beta D].
\]

Their sum is (tg-2)b/r=t(g-rho)b/r. This uses a0+beta d=-(d+a2), not a scalar-background assumption. The extra factor t arises because the g-regulator is differentiated in the moving-current double pole as well as in the fused endpoint derivative. The source-h regulator does not have this same dependence.

An arbitrary moving test f(z) changes the simple residue by f'(w)D2. This term is retained before doing the strip integral.

## 2. Action, both tails, and the canonical mixed contact

Put the barred current at x and define Q_g=B(a)U(x) with its inherited BB factor (x-a)^delta. All products in this paragraph include their ordinary uncontracted normal fields. Write Q_x for the left derivative of Q_g and Q_bar for its right derivative, so the endpoint spatial derivative is Q_x+Q_bar. Let r=x-a.

Suppress the one common ordered cocycle and charge normalization, and normalize the contact by 2 i nu lambda2. The strip calculation with I=pi and Jw=J1=-pi, using (2), gives

\[
\begin{array}{c|l}
\text{one-Y2 action}&-\beta[t(Q_x+Q_{\rm bar})-t\delta Q_g/r]\\
\text{moving B / barred K2 tail}&t[Q_x-t\delta Q_g/r]\\
\text{moving K2 / barred B tail}&t[\beta Q_x+tQ_{\rm bar}].
\end{array}
\tag{3}
\]

The antiholomorphic identity used here is beta a2-a0=t(d+a2). There is no extra antiholomorphic g factor: both repeated labels are chiral. Adding all three rows gives

\[
\boxed{\mathcal C_{\rm raw}/(2i\nu\lambda_2)
=t[\partial_xQ_g-\delta Q_g/r].}
\tag{4}
\]

This equals t times the insertion-wise derivative R_g partial_x Q_physical. It includes the derivative of the spectator C contractions. Replacing the spectator by a bare exponential would not establish (4).

The canonical finite matching, transported insertion-wise, adds

\[
\alpha[\partial_xQ_g-\delta Q_g/r].
\tag{5}
\]

Since t+alpha=t^2, the matched face is

\[
\boxed{\mathcal C_{\rm can}/(2i\nu\lambda_2)
=t^2[\partial_xQ_g-\delta Q_g/(x-a)].}
\tag{6}
\]

A convention that differentiates the *whole* regulated endpoint instead differs by an explicit delta term. The distinction must be made before studying its integrated union.

For the first moving-test jet, the action contribution is proportional to beta^2(I+Jw)=0, and the two double-pole tail terms are opposite. The minimal finite matching has no independent derivative-of-test contact. Higher moving-test jets cannot enter this isolated double-pole face. Thus (6) is not restricted to a constant moving test at this order. This does not settle every independent many-variable contact in higher products.

### The BBbarB union

There are two relative current coordinates and two source coordinates. The BB vertex contributes g; the three current/source powers contribute -3; the three possible derivative contractions contribute -3. Hence the full relative degree is g-2. Its successive Taylor divisors are g-2+j, j a nonnegative integer. On g=rho these are not identically zero. The proper repeated-current/source faces have their inherited g-1+j divisors, and the isolated mixed face has the fixed Cauchy normalization just computed. Generic rho avoids their discrete exceptional values.

Consequently a delta times meromorphic remainder from the BBbarB union vanishes on a generic diagonal. This verifies the repeated-current version of the projected braided Leibniz contact through lambda2 in the continued family, with the canonical t^2 insertion on either mixed face. It does **not** prove a literal sharp-cutoff identity at arbitrary physical coupling, or arbitrary-order associativity. The lower BB subtractions must retain both original labels when fused: a pair produces two cross-regulator factors against a further B, not one. The fresh script checks this pullback explicitly.

## 3. The five-current scalar channel and its barred tail

For B positions z1,z2,z3 and minority A at a, set

\[
R_{g,\rho}=\prod_{i<j}(z_i-z_j)^g\prod_i(z_i-a)^{-\rho},
\quad G(w)=\frac{(a-w)^2}{\prod_i(w-z_i)^2}.
\]

The scalar Wick polynomial is reconstructed directly from

\[
f_i=-\sum_{j\ne i}(z_i-z_j)^{-1}+(z_i-a)^{-1}
+\beta(z_i-w)^{-1}
\]

as f1 f2 f3 plus the three derivative-pair terms. Its source vertex factor gives

\[
K_{B^3A/Y2}=\frac{\beta\eta}{\rho}\partial_a[R_{g,\rho}G(w)].
\tag{7}
\]

The fresh check uses two asymmetric B triples and symbolic a,w,t. Adding a barred B at u supplies beta/(ubar-wbar)^2 in the fully contracted action channel. Thus the scalar action channel has both an a derivative and an ubar derivative. The three chiral K2 tails retain the latter derivative too.

The barred K2 tail is different: its holomorphic kernel is (7) evaluated at w=u. The immediate a primitive contains G(u). Write its exact partial fractions as

\[
G(u)=\sum_i\left[\frac{S_i}{(u-z_i)^2}+\frac{R_i}{u-z_i}\right],
\quad S_i=\frac{(a-z_i)^2}{\prod_{j\ne i}(z_i-z_j)^2}.
\]

The R_i are individually nonzero at generic positions, while their sum is zero. Thus there is no rational u primitive; a primitive contains R_i log(u-z_i). The sum-zero condition removes the overall scale logarithm, but it does not remove branch choices or proper-face contributions. An argument that gives every tail two ordinary rational derivatives would therefore be incomplete.

The full count can be checked without assuming (7). There are four relative current coordinates and two source coordinates; the current/current exponent is 3g-3rho, current/source exponent -3, and four contracted C factors cost four powers. The degree is

\[
\boxed{3(g-\rho)-1.}
\tag{8}
\]

The first jet has degree 3(g-rho), so generic-coupling exclusions alone cannot remove its possible diagonal pole.

## 4. An actual nonzero chamber residue of the decorated matching change

Fuse one B and the barred B at x, leaving B(z1), B(z2), A(a). Up to its fixed ordered phase, the scalar endpoint product is

\[
Q_g=\frac\eta\rho\partial_a[R_{g,\rho}S_x],\qquad
S_x=\frac{(a-x)^2}{(x-z_1)^2(x-z_2)^2}.
\tag{9}
\]

This identity was checked afresh by contracting the two remaining C factors, including their pair contraction. It is not inferred from endpoint-weight algebra. Insertion-wise versus whole-endpoint differentiation differs by

\[
\delta\frac\eta\rho
\left[\frac1{x-z_1}+\frac1{x-z_2}\right]
\partial_a[R_{g,\rho}S_x].
\tag{10}
\]

A local smooth test with a nonzero first a jet transfers the derivative in (10) onto that jet. To inspect one actual spatial chamber set

\[
x=0,\quad z_1=uR,\quad z_2=R,\quad a=uyR,
\qquad 0<u,y<1.
\]

The radial contribution is delta times the integral of R^(3delta-1). Removing the common coefficient and branch phase, the positive angular magnitude at g=rho is

\[
I(\rho)=\int_0^1\!du\int_0^1\!dy\,
(1+u)u^{-\rho}(1-u)^\rho
 y^{2-\rho}(1-y)^{-\rho}(1-uy)^{-\rho}.
\tag{11}
\]

For -1<rho<0 every displayed endpoint is integrable and the integrand is strictly positive. Therefore this chamber has a nonzero normalized diagonal residue I(rho)/3. At rho=-sqrt(2)/5 the independent quadrature gives

\[
I=0.28038157816939,\qquad I/3=0.09346052605646.
\]

For delta=.04,.02,.01,.005,.0025 the regulated values are .08566075,.08945381,.09142947,.09243794,.09294745, approaching that nonzero limit.

This calculation rules out discarding (10) solely because it contains delta. It does **not** establish a difference between the complete physical word products. A spatial chamber indicator is not a smooth test through the full diagonal. The other chambers, the other two mixed-fusion choices, their cocycles, and the proper faces can cancel this residue. Equal identical tests can also conceal a first-jet contact that independent slot tests detect. No fullword anomaly is claimed.

## 5. A concrete common contact prescription that cancels the total residue

There is an explicit way to define the *collapsed contact sector* so that its complete residue vanishes. Preserve the constituent-label multiset at every fusion. Give every remaining chiral coordinate, including U's left coordinate, a single inherited position in one total ordering of imaginary heights. Use that same height for its physical Wick factors and every inherited g factor. Continue the complete boundary-value product, before partitioning it into spatial chambers.

The following direct Fourier argument proves regularity of this prescription, including arbitrary smooth slot tests. For ordered heights c1>...>cn, orient each pair from the earlier to the later height. In an initial half-plane, a factor (xi-xj+i0)^(-pij) is the Fourier transform, up to its elementary entire phase, of

\[
\frac{k_{ij,+}^{p_{ij}-1}}{\Gamma(p_{ij})}.
\]

This one-dimensional distribution extends to an entire family in p by integration by parts: multiplying by 1/Gamma cancels its Taylor-subtraction poles. For a graph of such factors define qi as outgoing minus incoming positive edge momenta. Then

\[
\sum_i c_iq_i=\sum_{i\to j}(c_i-c_j)k_{ij}
\ge \min_{i\to j}(c_i-c_j)\sum_{i\to j}k_{ij}.
\tag{12}
\]

The incidence map is therefore proper on the positive edge cone. Its pushforward of the tensor product of the entire one-dimensional distributions is well defined and entire in all pair exponents. Wick derivative contractions merely shift exponents by integers and sum finitely many terms. Normal-field Taylor coefficients and arbitrary compact smooth spatial tests preserve this conclusion.

For the decorated contact (9), all singular coordinate factors are chiral; the right U vertex is inert with respect to the remaining coordinates. Thus (12) applies if the inherited-height rule is obeyed. The complete family multiplying delta in (10) has no diagonal pole. Its fullword diagonal restriction is zero. Consequently the nonzero individual chamber residues in section 4 **must cancel** in this particular whole-word prescription; no equal-test restriction is needed for that conclusion.

This is a construction and derivation, not an appeal to an abstract algebraic witness. The script checks (12) for all 24 total orders of the four remaining labels. Its role is to verify the orientation data; the inequality proves the general statement.

The premise is substantive. Opposite i0 assignments on two factors of the same coordinate difference admit a positive circulation and destroy the properness argument. A fusion that keeps the original moving-B height in R_g but uses a conflicting barred-current height in the physical U factor is not covered. The source/action/tail problem has both chiralities and two-dimensional integrations; it cannot be declared acyclic by this purely chiral argument.

Thus this section supplies a viable, regular contact normalization and eliminates a purported chamber-only no-go. It does not yet show that the complete five-current source forest induces it. That comparison must include the barred-tail logarithms, moving source holes, all three mixed fusions, and the like-current lower counterterms. The phase of the full endpoint product must be inherited from the actual cocycle, not reconstructed from a desired quantum-group relation.

## 6. Remaining ceiling and evidence

The concrete advancement is the BBbarB mixed-face calculation (1)–(6), including its first spatial-test jet and full-union divisor check. It is the lowest repeated-current mixed overlap required by a common product. The five-current calculation identifies the precise first-jet danger and gives a regular, coherent contact prescription that can remove it. Completing the full five-current action/tail/forest comparison remains the next discriminating calculation. A valid completion may prove equality with section 5 or exhibit a higher local contact forced by the mismatch. Neither outcome is established here.

Fresh scripts and retained outputs:

- `referee_scratch/sl21_common_product_20261024/overlap/check_overlap.py` and `check_overlap.out`: 19 PASS, 0 FAIL. Full Laurent coefficients, action/tail/matching jets, scalar five-current action derivative, barred-tail residues, collision degrees, fused-label pullback.
- `referee_scratch/sl21_common_product_20261024/overlap/check_contact_residue.py` and `check_contact_residue.out`: 29 PASS, 0 FAIL. Actual endpoint Wick identity, chamber residue and continuation, all 24 acyclic orders, opposite-orientation diagnostic.
- `check_contact_residue_initial.out` retains one failed *symbolic comparison* before simplification of two equal rational expressions. The calculation was unchanged; replacing structural equality by simplification makes that identity pass. No discarded mathematical counterexample is hidden.

These checks are supporting evidence. The Fourier proof, exact Laurent derivation and collision analysis specify their scope; the PASS counts do not certify full associativity.

## Integration

None. Independently verified: the t(g-rho) repeated-spectator pole term, all local BBbarB action/tail/finite contacts, the five-current degree, and the concrete decorated-contact chamber residue. Constructed: an acyclic inherited-height extension of the collapsed chiral contact sector, with an entire-family argument. Rejected: treating the h-source coefficient as the g-current coefficient; discarding a regulator multiple before checking its union; and inferring a full anomaly from one spatial chamber. Unresolved: matching this extension to the complete five-current mixed/Serre source forest, arbitrary-order product coherence, and the physical quantum algebra. No existing report or book text was edited or relabelled.
