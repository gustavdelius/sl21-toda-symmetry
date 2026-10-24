# Independent check of the 2017 normalized neutral form-factor degeneracy

Date: 2026-10-24. Model: Codex (GPT-6).
Commit of the surrounding repository: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

This is a separate follow-up note. The completed product-compatibility report and all previous files remain unchanged. No subagent was used.

## Verdict

**Correct for the normalized 2017 rational family, at generic coupling:** the three parameter values below have $J_{00}=1$ and

\[
J_{nn,p_0}=J_{nn,p_A}=J_{nn,p_U}=0\qquad(n>0).
\tag{F1}
\]

An all-$n$ induction is given below. It does not assume that a physical exponential operator has finite nonzero normalization, or that the root-system Toda endpoints equal the corresponding normalized form-factor fields.

**Unclear as a statement about the renormalized Toda endpoint limits:** the normalization factor, possible resonance limits, field identification, quantum vacuum, and operator domains are not fixed by (F1).

## Primary formulas and parameter map

Source: M. Lashkevich and Y. Pugai, [arXiv:1610.04926](https://arxiv.org/abs/1610.04926), supplied locally as `referee_scratch/sl21_vacuum_limits_20261024/lashkevich_pugai_2017.pdf`. The PDF skill was used for inspection. Printed pages 4–6 were rendered and read to check equations (2.7)–(2.17) and (3.1)–(3.9). The 2019 convention is not used.

The primary source qualifies the field identification by the exact phrase **“up to a normalization factor”**, immediately before (2.16). Its parameter map is

\[
p_\pm=1-g^{-1}+b\mp a.
\tag{F2}
\]

Here the source's $a,b$ are exponential-field labels, not the Toda coupling or endpoint amplitudes used in the other report. Set $h=g^{-1}$. The requested labels give

\[
\begin{array}{c|c}
\text{exponential label}&(p_+,p_-)\\\hline
V_{0,0,0}&p_0=(1-h,1-h)\\
V_{0,2h-2,0}&p_A=(h-1,h-1)=-p_0\\
V_{-h,-h,0}&p_U=(1-h,1-3h)=p_0+(0,-2h).
\end{array}
\tag{F3}
\]

These substitutions are direct applications of the 2017 formula. The separate proposed Toda-to-sine-Liouville field map is not audited in this note.

Let $\omega=e^{-2\pi ig}$ with powers interpreted using this displayed exponential, and put $v=\omega^{1/4}=e^{-i\pi g/2}$. Equations (2.7)–(2.9) give the elementary pair factors

\[
f_\epsilon(z)=\frac{v^{2\epsilon}z-v^{-2\epsilon}}{z-1},\qquad
g_\epsilon(z)=\frac{v^{-2\epsilon}z+v^{2\epsilon}}{z+1}.
\tag{F4}
\]

For opposite signs, two currents of the same particle type contribute $f_\epsilon(z_j/z_i)$, and a current followed by its dagger type contributes $g_\epsilon(y_j/x_i)$. Equal signs contribute one. Summing all signs with one-current coefficients $\omega^{\epsilon p_+/4}$ and $\omega^{\epsilon p_-/4}$ is a fresh direct implementation of $J_{KL,p}$.

For example, at $p_0$ both coefficients are $(iv)^\epsilon$. At arbitrary separated rapidities,

\[
J_{11,p_0}=-v^2-v^{-2}+g_+(y/x)+g_-(y/x)=0.
\tag{F5}
\]

This cancellation uses all four sign choices, including the mixed-pair factors.

## All-particle proof by a larger vanishing family

Define

\[
p_r=(1-h+r,1-h-r),\qquad r\ge0.
\]

The stronger statement is

\[
J_{K,L,p_{L-K}}(X|Y)=0
\quad\text{when }0\le K\le L,\ K+L>0.
\tag{F6}
\]

Use induction on $K+L$. At $K=0$, equation (3.6) gives

\[
J_{0,L,p_L}=\prod_{i=1}^{L}
2\cos\left[\frac{\pi g}{2}(p_{L,-}+L+1-2i)\right].
\]

The $i=1$ argument is $-\pi gh/2=-\pi/2$, so the product is zero for every $L>0$.

For $K\ge1$, apply the source's first recursion (3.3), with one of the $K$ particle variables singled out. Its constant term contains

\[
J_{K-1,L,p_r+\sigma}=J_{K-1,L,p_{r+1}},\qquad \sigma=(1,-1),
\]

which belongs to the induction family since $L-(K-1)=r+1$. Every pole term contains

\[
J_{K-1,L-1,p_r},
\]

which belongs to the same family with smaller total particle number. It therefore vanishes unless $K=L=1$, where this lower function is $J_{00}=1$. In that sole exceptional case the coefficient from (3.2) is

\[
2i\sin\pi g\left(\omega^{(p_+-p_-)/4}-\omega^{(p_--p_+)/4}\right)=0,
\]

because $p_{0,+}=p_{0,-}$. Thus both the constant and all residue terms vanish. This proves (F6) for generic separated variables. The expressions are rational functions of those variables, so it is a rational identity, and in particular proves $J_{nn,p_0}=0$ for every $n>0$.

The proof needs only recursion (3.3), its residue coefficient (3.2), and the base product (3.6). The script independently checks this recursion against generic sign sums before using the special parameters. No assumption that the $p_0$ field is the identity is used as an induction input.

## The two endpoint parameter values

First, direct substitution into (F4) gives

\[
f_\epsilon(z^{-1})=f_{-\epsilon}(z),\qquad
g_\epsilon(z^{-1})=g_{-\epsilon}(z).
\]

Changing every summation sign therefore proves the inversion identity

\[
J_{KL,-p}(X|Y)=J_{KL,p}(X^{-1}|Y^{-1}).
\tag{F7}
\]

This independently recovers source equation (3.8). Since $p_A=-p_0$, (F1) for $p_A$ follows immediately. It is unnecessary to use the more elaborate reflection relation (3.9) or any reflection amplitude. For comparison, at $K=L$ the prefactor of (3.9) is the empty-product ratio $J_{0,p_+}/J_{0,-p_-}=1$, so that relation leads to the same normalized conclusion; here $p_0^*=p_0$.

Second, the shift $p_-\mapsto p_--2h$ changes each dagger-current coefficient by

\[
\omega^{-\epsilon h/2}=e^{i\pi\epsilon}=-1.
\]

All pair factors are unchanged. Hence, for all $K,L$,

\[
J_{KL,p_U}=(-1)^L J_{KL,p_0}.
\tag{F8}
\]

This proves the last equality in (F1). It also directly verifies the appropriate instance of (2.17), with its sign and branch convention. Using principal-value powers of a separately rounded complex $\omega$ would lose this branch information and is not an appropriate check.

## The rapidity prefactor and the normalization caveat

Equation (2.15) is $F_p=J_{KL,p}R_p$. The $p$ dependence of (2.13) is

\[
R_p\ \propto\ 
\exp\left[\frac g2(p_+-p_-)
\left(\sum\theta_i-\sum\theta'_j\right)\right].
\]

Consequently

\[
R_{p_A}=R_{p_0},\qquad
R_{p_U}=\exp\left(\sum\theta_i-\sum\theta'_j\right)R_{p_0}
=\frac{\prod x_i}{\prod y_j}R_{p_0}.
\tag{F9}
\]

The periodicity of $J$ alone is therefore not equality of generic full form factors. At the special values in (F3), however, (F1) is an exact zero rational function. Multiplication by the prefactor cannot remove that zero at generic coupling and nonsingular rapidities. Kinematic-pole values are understood by meromorphic continuation from separated variables; there is no hidden nonzero residue of the identically zero function.

This statement has a precise limitation. If a physical field is represented as $G(p)F_p$ and $G(p)$ is singular at the special parameter, the physical definition may require a joint limit. For example, $G(p_*+\eta)\sim g_{-1}/\eta$ and $F_{p_*+\eta}\sim\eta F'_*$ can produce a finite nonzero limit. Derivatives of $R_p$ can matter in such parameter limits, although they cannot undo (F1) at a fixed regular parameter. Coupling resonances where the common prefactor itself requires continuation also need their own treatment. Neither finiteness of $G$ nor absence of operator mixing is established by this computation.

Thus the result is **identity-like normalized neutral form factors**, including $J_{00}=R_{00}=1$. It does not determine the normalization of either Toda endpoint, prove scalar weak boundary limits, or justify replacing a fermion bilinear endpoint by a bosonic vertex before renormalization.

## Independent checks

[ff_identity_checks.py](ff_identity_checks.py) constructs exact Laurent-polynomial sign sums in the indeterminate $v$, with rational separated rapidity variables. It checks $n=1,2,3,4$ independently for all three parameter values, and seven unequal-particle cases of the larger induction family. A distinct rational implementation compares ten generic sign sums against recursion (3.3), with independent zero-mode coefficients and a non-root-of-unity algebraic probe. It also checks the inversion kernels, the arbitrary-rapidity two-particle cancellation, the exceptional residue, and the zero-particle normalization.

The [output](ff_identity_checks.out) records **37 PASS, 0 FAIL**. These finite checks supplement the all-$n$ derivation above; they are not its replacement. No preceding check script was imported or rerun.
