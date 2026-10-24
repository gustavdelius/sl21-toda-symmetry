# Sine-Liouville reflection zero and its identity resonance

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed commit: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

This is an additional audit requested by the parent. The completed global-boundary report is unchanged. The source is V. A. Fateev, *Normalization Factors, Reflection Amplitudes and Integrable Systems*, [hep-th/0103014v2](https://arxiv.org/abs/hep-th/0103014), equations (126), (130)–(132), printed pp. 32–33. The supplied PDF was read and p. 33 rendered and visually inspected as `fateev_p33.png`. No correction to the published source is proposed.

## Source reduction

Use source rank one with short root of length one. Then the positive long-root product is empty, the fundamental weight is 1/2, and the background charge is

\[
Q=\frac1{2z},\qquad \alpha^2=1+z^2,\qquad
C=\frac{\pi\mu'}{2z^2}.
\]

Equation (131) has the orientation
\(R(a,\eta)=A(2Q-a,\eta)/A(a,\eta)\), consistent with (130),
\(V_{a,\eta}=R(a,\eta)V_{2Q-a,\eta}\).
The rank-one reduction of the printed equation (132) is

\[
A(a,\eta)=C^{(a-Q)/z}
\frac{\Gamma(1-2(a-Q)z)\Gamma(1-(a-Q)/z)}
{\Gamma(1/2-(a-Q)z+\alpha\eta)
 \Gamma(1/2-(a-Q)z-\alpha\eta)}.
\tag{R1}
\]

The exponent in the source is the **barred** momentum, \((a-Q)/z\), not \(a/z\). Replacing one by the other multiplies all \(A\)'s by the same constant and therefore leaves the reflection ratio unchanged. This was a transcription issue in the proposed reduction, not an error in Fateev's formula.

For the global scope, the exact source sentence on p. 33 is: “Until now we considered spinless fields and non-compactified field Φ.” The following compactification discussion is additional. Thus using (R1) by itself does not establish the compactification or spin-sector identification of the original fermionic model.

## Simple zero at the endpoint momentum

Write \(v=z^2=-t/2\), set \(\eta=-\alpha\), and perturb the endpoint momentum by

\[
a=a_*+\varepsilon,
\qquad a_*=z+z^{-1},\qquad
a'=2Q-a=-z-\varepsilon.
\]

The four gamma arguments in the numerator and denominator of \(A(a,\eta)\), in that order, are

\[
-2v-2z\varepsilon,\quad
-\frac1{2v}-\frac\varepsilon z,\quad
-1-2v-z\varepsilon,\quad 1-z\varepsilon.
\]

The four reflected arguments are

\[
2+2v+2z\varepsilon,\quad
2+\frac1{2v}+\frac\varepsilon z,\quad
z\varepsilon,\quad 2+2v+z\varepsilon.
\]

At generic positive \(v\), away from simultaneous gamma poles,

\[
A(a_*,\eta)
=C^{1+1/(2v)}(-1-2v)\Gamma\!\left(-\frac1{2v}\right)
\]

is finite and nonzero. Since \(1/\Gamma(z\varepsilon)=z\varepsilon+O(\varepsilon^2)\),

\[
R(a_*+\varepsilon,-\alpha)
=\varepsilon\,r_1+O(\varepsilon^2),
\]
\[
r_1=
zC^{-2-1/v}
\frac{\Gamma(2+1/(2v))}
{(-1-2v)\Gamma(-1/(2v))}.
\tag{R2}
\]

Thus the generic reflection zero is simple. Exceptional couplings at intersecting gamma singularities need a separate multivariable limiting prescription. The numerical checks use irrational values of \(v\), count the order by \(R/\varepsilon\), and do not infer absence of a pole from a small printed value.

## The reflected operator has an identity resonance

The zero in (R2) cannot by itself imply a zero endpoint VEV. With canonical free covariance, the reflected vertex has momentum
\(c=(-z-\varepsilon,-i\alpha)\), while the opposite Yukawa screening has momentum \(a_0=(z,i\alpha)\). Directly,

\[
c\cdot a_0=1-z\varepsilon,
\qquad c+a_0=(-\varepsilon,0).
\]

Consequently

\[
V_{-z-\varepsilon,-\alpha}(0)Y_0(r)
\sim |r|^{-2+2z\varepsilon}V_{-\varepsilon,0}(0).
\tag{R3}
\]

The local integrated coefficient has a simple pole,

\[
\int_{|r|<R}d^2r\,|r|^{-2+2z\varepsilon}
=\frac{\pi R^{2z\varepsilon}}{z\varepsilon}.
\tag{R4}
\]

In the source normalization, the expansion of the action contributes \(-\mu'\int Y_0\). It therefore gives reflected-operator residue
\(-\pi\mu'/z\) times the identity. In the original ordinary-area convention \(S_{\rm int}=\pi^{-1}\int\lambda_0Y_0\), the residue is \(-\lambda_0/z\). The sign follows from the minus action insertion. The pole is a resonance mixing with the identity; it is not an unspecified long-distance condensate.

Within the given perturbation class, it is the unique nontrivial generic-coupling UV resonance at this reflected momentum. Let \(n_0,n_2\) be Yukawa insertion counts, \(m\) the even-root insertion count, \(n=n_0+n_2\), \(D=n_0-n_2-1\), and \(\ell\ge0\) a descendant dimension. The fused field plus coupling dimensions equal

\[
D^2+n+2m+\ell
+v\{D^2-(n-1-2m)^2+4m\}.
\]

Matching the original reflected-field dimension 1 at generic \(v\) requires the constant and \(v\) coefficients separately. The constant equation
\(D^2+n+2m+\ell=1\) leaves only the trivial no-insertion term and \((n_0,n_2,m,\ell)=(1,0,0,0)\).

For the unreflected \(U\), put \(A=1+n-2m\). Its resonance equations are

\[
D^2-A^2+4m=0,\qquad
D^2-1-n+6m+\ell=0.
\]

Eliminating \(D^2\) gives \(\ell=A(1-A)\). Since \(A\) is an integer and \(\ell\ge0\), only \(A=0,1\) and \(\ell=0\) remain. Nonzero \(m\) then forces \(D^2<0\), and \(m=0\) leaves only the no-insertion term. Thus \(U\) itself is nonresonant in this generic-coupling perturbation class. The enumeration in the script is a second check of this analytic argument.

## Conditional finite value and its normalization test

If the meromorphic reflection identity is valid for the massive operator family in the chosen completion, and the original family is its regular minimally continued Wick representative, the zero times the pole gives a finite contribution. For symmetric Yukawa coupling, it is

\[
G_U^{\rm refl}
=C^{-1-1/v}
\frac{\Gamma(1+1/(2v))}{\Gamma(-1/(2v))}.
\tag{R5}
\]

With \(C=\sqrt{\lambda_0\lambda_2}/(-t)\), unequal positive Yukawa coefficients are reduced by a constant angular shift. Its endpoint factor supplies the correct \(\lambda_0\), giving the equivalent formula

\[
G_U^{\rm refl}
=\frac{\lambda_0}{1-t}
\left(\frac{\sqrt{\lambda_0\lambda_2}}{-t}\right)^{2/t-2}
\frac{\Gamma(2-1/t)}{\Gamma(1/t)}.
\tag{R6}
\]

This is a **conditional analytic-reflection value**, not an independently constructed original-theory vacuum expectation. Its generic nonzero value shows precisely why inserting \(R=0\) before examining the reflected operator is invalid. Likewise, normalized form factors obtained by dividing by a VEV cannot decide this resonant normalization by themselves.

There is a useful direct matching check with the earlier canonical mixed endpoint. The original local OPE is

\[
U(0)Y_0(r)=|r|^{2t-4}
\bigl(F+rCF+\bar r\bar CF+r\bar r C\bar CF+\cdots\bigr).
\]

Multiplying its two angularly surviving radial primitives by \(-\lambda_0/\pi\) gives the cutoff terms
\(+\lambda_0\epsilon^{2t-2}F/(t-1)\) and
\(+\lambda_0\epsilon^{2t}C\bar CF/t\).
The minimal analytic extension therefore uses

\[
U_{\rm can}=U-rac{\lambda_0}{t-1}\epsilon^{2t-2}F
-\frac{\lambda_0}{t}\epsilon^{2t}C\bar CF,
\qquad B_{BB}=t^2\lambda_2U_{\rm can}
\tag{R7}
\]

through the displayed order. For \(-1<t<0\), these are exactly the two divergent angularly invariant Taylor terms; the next has positive radial degree. Thus the minimally continued local primary agrees with the specified canonical endpoint subtraction at this order. This makes (R6) a concrete conditional bridge, while leaving higher-order finite choices, massive reflection validity, compactification/domain identification, and infinite-volume weak limits unproved.

## Checks and scope

`reflection_zero_audit.py` independently computes every gamma argument, the source-power difference, the OPE pairing and pole, both resonance polynomials, and numerical simple-zero/finite-product scaling at three generic irrational couplings. Its output is **28 PASS, 0 FAIL**. The source page image, script, and output are all new scratch files. The completed boundary-global report remains unchanged.

Verdict: the reflection zero is **correct at generic coupling**; the inference of a zero massive endpoint VEV is **unsupported** because the reflected operator has an identity resonance. Formula (R6) is **conditional** on the analytic-reflection completion and matching described above. No no-go claim, published-source correction, or physical symmetry-algebra identification follows from this audit.
