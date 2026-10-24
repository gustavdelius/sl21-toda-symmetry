# Root-system sl(2|1) Toda: boundary and product compatibility test

Date: 2026-10-24. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: independently test the algebraic compatibility of the proposed canonical mixed products, continued same-chirality quartics, and conditional boundary limits. `AGENTS.md` was read. New work is confined to this report and `referee_scratch/sl21_product_compatibility_20261024/`. No existing file was edited, no previous script was imported or rerun, no subagent was spawned, and no git mutation, commit, push, or book render was performed.

## Outcome

There is no obstruction at the algebraic overlap tested here. The proposed boundary differences have precisely the two powers of the boundary variable that pass a coefficientwise quartic/mixed compatibility test. Their two amplitudes remain independent. An explicit four-dimensional associative matrix model realizes the candidate relations with both amplitudes arbitrary. This is evidence of algebraic consistency, not a construction of this Toda theory's charges or an identification of its symmetry algebra.

There is also a concrete, lower-order condition on the missing three-contour product. Let $A,B,\bar A,\bar B$ denote the four candidate charges and define the same-chirality composite

\[
\mathcal D=BA-s^{-1}AB,\qquad s=e^{-2\pi i/t}.
\]

Write the diagonal mixed defects as $R_A-L_A$ and $R_B-L_B$, including their common factor $\nu/(2\pi i)$ in the endpoint definitions. Here $R$ and $L$ are right and left spatial endpoints, not chirality labels. Any common associative product must give

\[
\boxed{[\mathcal D,\bar A]=(1-s^{-2})BR_A,\qquad
[\mathcal D,\bar B]=(s^{-1}-s)AL_B.}
\tag{P1}
\]

These conditions already apply before replacing endpoints by vacuum expectation values. They specify the part of the three-contour contact that the finite canonical matching must produce. They also show why the separated same-chirality phase cannot be promoted to the global relation $BA=s^{-1}AB$ when the displayed right-hand sides act nontrivially.

| Item | Verdict | Independent evidence and scope |
|---|---|---|
| Endpoint weights needed in (P1) | **Correct on the stated canonical cover** | Direct integration of the canonical delta function and root products |
| Conditions (P1) | **Correct as necessary associative identities** | Direct expansion of each three-letter word |
| Compatibility of endpoint powers with the quartics | **Correct coefficientwise** | Free-word expansion with symbolic $s$, detailed below |
| A relation between the two endpoint amplitudes forced by this test | **Wrong** | Explicit matrix model permits arbitrary independent amplitudes |
| Compatibility of the actual mixed and four-fold field-theory prescriptions | **Unclear** | A common regulated three-contour extension has not been evaluated |
| Physical vacuum sector or nonzero endpoint limits | **Unclear** | Not supplied by an abstract matrix model |

The “wrong” entry rejects a possible inference, not a sentence in an existing report. No existing text requires a correction from this pass, and no theorem label is proposed.

## 1. Inputs and what is re-derived

The candidate binary relations are

\[
A\bar A-s\bar AA=f_A,\qquad B\bar B-s\bar BB=f_B,
\qquad A\bar B=s^{-1}\bar BA,\qquad B\bar A=s^{-1}\bar AB,
\tag{P2}
\]

with $f_a=R_a-L_a$. The mixed-regulator report, equations (19)–(20), specifies the finite changes of these endpoints. The mixed-audit report, line 153, correctly retains the freedom in density contacts not detected by a constant outer spatial test function. The mass-two Serre report, lines 14–22, specifies the two quartics as four-fold continued products rather than iterated associative products. Those restrictions matter here.

This pass does not recompute the radial finite-contact coefficients or the mass-two continuation integrals. Instead it takes their candidate identities as the object of an independent algebraic consistency test. No conclusion about the equality of these different product constructions is assumed in deriving the necessary conditions.

For the endpoint calculation use $X=(\varphi^u,\varphi^v,h)$ and

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\quad
a_2=(ib/\sqrt2,-b/\sqrt2,-i),\quad
T=a_0+a_2,\quad d=-T/t,\quad t=b^2.
\]

The explicit vectors give $Td=2$. With

\[
\varphi(x)=\tfrac12\left(\Phi(x)-\int_{-\infty}^x\dot\Phi\right),
\quad [\Phi(y),\dot\Phi(x)]=4\pi i\delta(x-y),
\]

one obtains $[\Phi(y),\varphi(x)]=-2\pi i\Theta(x-y)$. Thus a current of disorder momentum $qd$ has $[n_L,Q_q]=-2qQ_q$, $[n_R,Q_q]=0$ when $L<x<R$, where $n_{L,R}=T\Phi(L,R)/(2\pi i)$. The endpoint exponents $k=d+T$ and $d+a_2$ satisfy $kd=2-2/t$ and $(d+a_2)d=1-2/t$. Their integer difference has trivial exponential phase. Consequently

\[
L_AB=s^2BL_A,\quad L_AA=s^{-2}AL_A,
\quad L_BB=s^{-2}BL_B,\quad L_BA=s^2AL_B.
\tag{P3}
\]

Right endpoints commute with currents supported to their left. The corresponding antichiral weights are reversed. Counterterms with the same momenta inherit these phases; derivative contacts are absent at separated endpoints. These statements concern the same anchored canonical cover and scalar cocycle convention as the reports. A different zero-mode or boundary-domain realization would need its own check.

## 2. The first useful three-product condition

In any associative product, the elementary braided Leibniz identity is

\[
XY Z-r_Xr_Y ZXY
=X(YZ-r_YZY)+r_Y(XZ-r_XZX)Y.
\tag{P4}
\]

It follows by expanding the right side, so no Jacobi formula is imported. Apply (P4) to $BA$ and $AB$ using (P2). Their total phase with either barred generator is one. Direct subtraction gives

\[
[\mathcal D,\bar A]=Bf_A-s^{-2}f_AB,
\qquad
[\mathcal D,\bar B]=s^{-1}(f_BA-Af_B).
\tag{P5}
\]

Substitution of $f_a=R_a-L_a$ and (P3) gives (P1). This is a useful projected associativity condition: it is weaker than determining every component of the three-product associator, but it fixes a nonzero three-contour contact using only independently separated endpoint products.

If the right boundary is fixed in a common orbit, write $K=\nu^N$, $\nu^2=s$, so that $KA=s^{-1}AK$ and $KB=sBK$. The conditional scalar boundary formulas become

\[
R_A=a,\quad L_A=aK^2,\qquad R_B=b,\quad L_B=bK^{-2}.
\tag{P6}
\]

Here $a=\nu\mathcal A_\alpha s^{-r}/(2\pi i)$ and $b=\nu\mathcal B_\alpha s^r/(2\pi i)$; $r$ is the fixed right label. Equation (P1) then reads

\[
[\mathcal D,\bar A]=a(1-s^{-2})B,
\qquad [\mathcal D,\bar B]=b(s^{-1}-s)AK^{-2}.
\tag{P7}
\]

For generic $s$ and invertible $K$, setting $\mathcal D=0$ in such a representation forces $aB=bA=0$. No such pair relation is imposed by the reported quartics. They leave room for the required composite. At the level of endpoint operators, the corresponding conclusion is only that $BR_A$ and $AL_B$ must vanish if $\mathcal D$ vanishes; zero divisors or special domains cannot be excluded by this argument.

## 3. Crossing the quartics with the barred generators

Let

\[
\mathcal S_B=B^3A-cB^2AB+cBAB^2-AB^3,
\qquad c=1+s+s^{-1},
\]

and define $\mathcal S_A$ by interchanging $A$ and $B$. To test whether (P2) preserves the proposed quartic ideal, move one barred generator through each word using (P2). For an occurrence of a matching letter $X$ in $W=UXV$, its contribution is

\[
(\text{phase of }V)\ U f_X(K)V.
\tag{P8}
\]

For $f_X=K^m$, moving the endpoint to the right multiplies it by $s^{m(\#B(V)-\#A(V))}$. This gives a short, completely explicit calculation in free words.

For the minority contact in $\mathcal S_B\bar A-s^{-2}\bar A\mathcal S_B$, the coefficient of $B^3K^m$ is

\[
p(s^{m-1}),\qquad
p(x)=1-cx+cx^2-x^3=(1-x)(1-sx)(1-s^{-1}x).
\tag{P9}
\]

For the majority contact in $\mathcal S_B\bar B-s^2\bar B\mathcal S_B$, put $y=s^{m+1}$. The result is

\[
\left(y+1+y^{-1}-c\right)(B^2A-yAB^2)K^m,
\qquad y+1+y^{-1}-c=\frac{(y-s)(y-s^{-1})}{y}.
\tag{P10}
\]

The middle word $BAB$ cancels identically. The $A$-majority formula instead has $y=s^{1-m}$. The other minority formula is $p(s^{-m-1})$.

At generic $s$, simultaneous coefficientwise cancellation for finite Laurent polynomials therefore allows

\[
f_A=a_0+a_2K^2,\qquad f_B=b_0+b_{-2}K^{-2}.
\tag{P11}
\]

The adjective “coefficientwise” is essential: necessity assumes distinct endpoint powers and cubic words are not identified by extra relations. In a special representation further cancellations can occur. Sufficiency of (P11) for these four overlap tests does not need that assumption.

The boundary differences (P6) have precisely the permitted form. The test does not impose $a_2=-a_0$ or $b_{-2}=-b_0$; those ratios come from equal-orbit endpoint differences. It also gives no relation between $a$ and $b$. Thus the different finite normalizations of the two channels are not ruled out by the first quartic/mixed overlap. A cancellation in this calculation is not an all-order perturbative Serre proof: it checks the algebraic compatibility of the proposed leading relations if a common associative realization exists.

## 4. An explicit nonzero associative witness

One can demonstrate that the preceding consistency is nonvacuous without assigning the algebra a name. Treat $v$ as an indeterminate and set $s=v^2$. On a four-dimensional vector space with basis $e_j$, $j=0,1,2,3$, define

\[
Ke_j=v^{2j-3}e_j,\quad Ee_j=e_{j+1},\quad
Fe_j=[j]_v[4-j]_v e_{j-1},\quad
[j]_v=\frac{v^j-v^{-j}}{v-v^{-1}},
\]

with out-of-range basis vectors zero. Direct multiplication gives

\[
KE=sEK,\qquad KF=s^{-1}FK,\qquad
EF-FE=\frac{K-K^{-1}}{v-v^{-1}}.
\]

Choose arbitrary scalar symbols $a,b$ and put

\[
A=F,\quad B=E,\quad
\bar A=\frac{a(v-v^{-1})}{s}KE,\quad
\bar B=\frac{b(v-v^{-1})}{s}K^{-1}F.
\tag{P12}
\]

Then all four mixed relations (P2), the endpoint forms (P6), both chiral quartics, both antichiral quartics, and the triple conditions (P7) hold by direct symbolic matrix multiplication. Also $E^3\ne0$ and $\mathcal D\ne0$. The arbitrary independent symbols $a,b$ prove that these relations alone do not force either amplitude to vanish or tie their ratio.

This is not a Toda Hilbert-space representation. No spacetime currents, vacuum vector, positive inner product, Lorentz spin, coupling expansion, or common domain of unbounded charges has been constructed. In particular, (P12) is not evidence that the physical antichiral operators are the specified composites of the physical chiral ones. Its sole role is an explicit associative algebraic witness.

## 5. The finite-contact calculation still required

Let a prospective common product use the reported finite endpoint changes. In the endpoint normalization of this report they are

\[
\delta R_A=\frac{\nu}{2\pi i}[-\gamma\lambda_1F_A(R)],\qquad
\delta L_B=\frac{\nu}{2\pi i}[t(t-1)\lambda_2\mathscr U(L)]
\tag{P13}
\]

at their first nonzero mass orders, with $\gamma=t/(1-t)$. The mass-two dressing of $\delta L_B$ is the one displayed in mixed-regulator equation (20). Assuming the same-chirality product itself is held fixed, compatibility requires the induced change of the three-contour defects to satisfy

\[
\boxed{\delta[\mathcal D,\bar A]=(1-s^{-2})B\delta R_A,\qquad
\delta[\mathcal D,\bar B]=(s^{-1}-s)A\delta L_B.}
\tag{P14}
\]

Here $\delta$ denotes the change of the renormalized three-product defect, including its nested subtractions; it is not an ordinary commutator with an unchanged product. Equations (P13)–(P14) give the normalization target with all overall charge factors exposed. The order-$\lambda_2$ second relation is the lowest-mass nontrivial test.

A concrete field-theory calculation should therefore evaluate the four three-contour words in

\[
(BA-s^{-1}AB)\bar B-\bar B(BA-s^{-1}AB)
\]

with one $Y_2$ action insertion and all order-$\lambda_2$ current tails, using the same current-circle exclusions, ordered time separations, source cutoff, and analytic continuation for both bracketings. The binary canonical matching must be inserted on every nested mixed face. The simultaneous three-current/source face must then be extended rather than inferred from the binary computation. At the next order its $Y_0$ dressing and the endpoint counterterms must be included together. These are finite-interval calculations with endpoints outside the current support; a vacuum limit is not needed to perform them.

The reported outer-test-function restriction matters: differentiating a test function in a nested product can expose density contacts invisible in the binary integrated identity. Equations (P14) constrain the projected sum of those contacts. If it fails, a further allowed three-product contact or a change in a lower product would have to be specified and checked against the other overlaps. A failure of one unadjusted regulator would not by itself prove that no compatible associative extension exists.

For an associative reference product $\mu_0$, a first-order change $D$ must obey the direct associator condition

\[
D(\mu_0(X,Y),Z)+\mu_0(D(X,Y),Z)
-D(X,\mu_0(Y,Z))-\mu_0(X,D(Y,Z))=0.
\tag{P15}
\]

This follows by expanding associativity once. Binary endpoint matching specifies only a projected part of $D$; it does not determine all four terms. The new tests (P1) and (P14) specify useful projections that can be calculated next. Neither the existing four-fold analytic construction nor this algebraic audit evaluates that missing simultaneous face.

## Evidence and integration

The fresh script [compatibility_checks.py](../referee_scratch/sl21_product_compatibility_20261024/compatibility_checks.py) uses symbolic root products, free-word reduction, and independently constructed matrices. Its [output](../referee_scratch/sl21_product_compatibility_20261024/compatibility_checks.out) contains **60 PASS, 0 FAIL**. The script imports no prior investigation code. Symbolic indeterminates avoid accidental root-of-unity tests. The written derivation supplies the generic mode classification beyond the finite integer-mode sweep.

No book integration was performed. **Independently verified:** anchored endpoint weights, necessary triple-product identities, four quartic/mixed coefficientwise overlaps, and an explicit nonzero associative witness. **Extended:** a precise three-contour target for propagation of the canonical finite contact. **Unresolved:** equality of the actual analytic and canonically matched products with iterates of one associative product; the missing simultaneous current/source extension; physical vacuum sectors and scalar weak endpoint limits. No physical symmetry algebra is established.

Input hashes are recorded in [inputs_before.sha256](../referee_scratch/sl21_product_compatibility_20261024/inputs_before.sha256). Resumption information is in [checkpoint.md](../referee_scratch/sl21_product_compatibility_20261024/checkpoint.md).
