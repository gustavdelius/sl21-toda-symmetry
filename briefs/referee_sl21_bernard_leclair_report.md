# Bernard–LeClair construction: correction of direction and initial OPE tests

Date: 2026-10-23. Model: Codex (GPT-6).

Book baseline: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

This note supersedes the choice of research direction in `briefs/referee_sl21_nonlocal_report.md`, without editing that report or the book. Its monodromy computations do not answer the intended Bernard–LeClair charge-construction question. Their translation and endpoint failures concern those particular transport quantities, not Bernard–LeClair currents.

## Source read

Bernard and LeClair, *Quantum group symmetries and non-local currents in 2D QFT* (1991), the editor's PDF at `/home/gustav/Zotero/storage/QJH8RCYM/Bernard and LeClair - 1991 - Quantum group symmetries and non-local currents in.pdf`.

Read sections 2b and 5a, and the ultraviolet chiral split in section 3. The printed page 123 was rendered and inspected to recover equations (5.3)–(5.5), which text extraction misses. The precise test in section 2b is that “the residues of these OPE's are total derivatives”.

In their canonical normalization, the free chiral propagator is phi^a(z) phi^b(w) ~ -delta^{ab} log(z-w), and the perturbing bosonic vertex is exp(-i b alpha.Phi), where b denotes their beta-tilde, not the book's physical beta. For a length-two simple root,

\[
J_i=:e^{i\alpha_i\cdot\phi/b}:,\qquad
H_i=-\gamma\lambda_i:
e^{i(1/b-b)\alpha_i\cdot\phi-i b\alpha_i\cdot\bar\phi}:,
\qquad \gamma=\frac{b^2}{1-b^2}.
\]

Their construction satisfies partial_bar J_i=partial H_i to first order. The charge integrates J_i dz+H_i dbar_z. The antichiral construction is analogous. First-order conservation and an all-order argument are distinct obligations; the latter cannot be assumed for a superalgebra extension.

These currents use the chiral field's disorder string. They have no explicit exponential spatial weight of the sort introduced by the vacuum-gauged transport in the preceding report. That report therefore gives no obstruction to this construction.

## Appropriate superalgebra problem

Retain the root Gram matrix

\[
(\alpha_i\cdot\alpha_j)=
\begin{pmatrix}0&-1&1\\-1&2&-1\\1&-1&0\end{pmatrix}.
\]

Use a formal canonical ultraviolet realization with perturbing chiral factors

\[
V_1=:e^{-i b\alpha_1\cdot\phi}:,\qquad
V_j=:(\alpha_j\cdot\psi)e^{-i b\alpha_j\cdot\phi}:\quad(j=0,2),
\]

and their antichiral partners. The component action also has a bosonic exponential for alpha_0+alpha_2=-alpha_1. Couplings, normal ordering, composite counterterms and the continuation to the previously used complex action must be fixed before calling this a fully specified quantum theory. In particular, b has not been identified with the source action's g or the book's beta.

Seek even and odd currents in the form

\[
J_i=:\mathcal P_i(\psi,\partial\phi,\partial\psi,\ldots)
e^{i a_i\cdot\phi}:,
\qquad |\mathcal P_i|=|\alpha_i|,
\]

enlarging the disorder dressing if necessary. For each perturbation, the residue must be the derivative of an admissible local holomorphic counterterm h_{ia}. Then H_i is the sum of lambda_a h_{ia} times the appropriate antiholomorphic perturbing factor. Mutual locality with the perturbation, higher-order conservation, braiding and charge algebra still require verification.

The odd roots have squared length zero. One must not replace them by coroots 2 alpha_i/alpha_i^2. Nor is it established that simply multiplying the bosonic current by alpha_i.psi is sufficient.

## Independent first-order checks

The new scratch program `referee_scratch/sl21_bernard_leclair/ope_checks.py` uses the canonical propagators above and independently recomputes the following.

**Even node: correct at first order in the listed perturbing couplings.** For J_1=exp(i alpha_1.phi/b), the OPE with V_1 has a double pole. Expanding the perturbation at fixed current insertion gives the residue

\[
-i b\,(\alpha_1\cdot\partial\phi)
e^{i(1/b-b)\alpha_1\cdot\phi}
=-\gamma\partial
e^{i(1/b-b)\alpha_1\cdot\phi}.
\]

This recovers Bernard–LeClair's coefficient, not just the exponential. Its OPEs with both Yukawa holomorphic factors are regular: alpha_1.alpha_j=-1 gives a positive first power of the separation, and J_1 contains no fermion to contract. The composite bosonic exponential is also regular, since alpha_1.(alpha_0+alpha_2)=-2. Consequently the even-node current has the required first-order conservation law in this formal ultraviolet setup. This is not an all-order result, particularly because two Yukawa insertions and the composite potential contribute at the same order in the original mass parameter.

**Minimal odd ansatz with plain-vertex counterterms: wrong.** Consider J_0=(p.psi) exp(i alpha_0.phi/b), with constant vector p. The node-0 Yukawa factor produces a simple-pole residue proportional to alpha_0.p, without a derivative. It must vanish in this ansatz. The node-2 factor produces a double pole and a residue containing the normal product (alpha_2.psi)(p.psi), multiplied by a bosonic exponential. A derivative of a plain bosonic vertex cannot absorb that fermion bilinear. Removing it requires det(alpha_2,p)=0. These conditions together force p=0. In particular p=alpha_0 does not pass.

This rejects only this minimal ansatz with these counterterms. Derivative prefactors, additional disorder dressing, a bosonized fermion realization and other admissible counterterms have not been excluded. The same issue occurs with 0 and 2 exchanged.

The script checks the residue coefficient, all relevant root products and the two-vector linear system exactly. Every check passes. No claim of prior-art novelty is made.

## Integration and present status

No book integration was performed. Corrected interpretation: the sought charges should be constructed by the Bernard–LeClair ultraviolet/OPE method, not inferred from auxiliary monodromy. Verified: the even-node first-order current and failure of the restricted bare odd-current ansatz. Unresolved: suitable odd currents, a complete quantum prescription, higher-order conservation, the superalgebra relations, q as a function of the action coupling, and the action on solitons. The preceding monodromy failures must not be cited as evidence against this program.
