# moment_response.md

Equivalent root bending moment of the 2-DOF wing section under the critical certification gust. Companion to `moment_response.m`.

## Physical setup

The section is a NACA 2412 airfoil with the properties of `build_ss_nl.m` (m = 2 kg, k1 = k2 = 200 N/m, knl = 5e4 N/m³, d = 0.25 m, chord 0.30 m, U = 12 m/s). The airfoil enters the model only through the thin-airfoil lift slope $C_{l\alpha}=2\pi$ and the chord — camber and thickness play no further role at this fidelity.

The wing length $L$ is given explicitly as an input argument (default $L = 0.50$ m, matching the section span). Local effects at the spring attachments are neglected: the section loads are assumed to transfer rigidly to the root, so a vertical force $F$ at the section produces a root bending moment $F \cdot L$.

## Equivalent moment equation

Each spring deflects by $\delta_i = h + r_i\theta$ and carries

$$F_i = k_i\,\delta_i + k_{nl,i}\,\delta_i^3 .$$

The structural restoring loads on the section are a root bending moment $M_b$ (vertical force over lever arm $L$) and a torsional moment $T$:

$$M_b = L\,(F_1 + F_2), \qquad T = r_1 F_1 + r_2 F_2 .$$

These are combined with the maximum normal stress (machine-design) equivalent bending moment — the pure bending moment producing the same peak normal stress as the combined loading:

$$\boxed{\;M_{eq}(t) = \tfrac{1}{2}\Big(|M_b| + \sqrt{M_b^2 + T^2}\Big)\;}$$

$|M_b|$ is used so down-bending half-cycles register with their full magnitude; the raw signed formula collapses toward zero whenever $M_b < 0$. The formula is strictly derived for circular shafts; for the airfoil section it serves as the standard scalar proxy for combined bending–torsion severity. In the symmetric linear case ($k_1 = k_2 = k$, $r_2 = -r_1$, $k_{nl}=0$),

$$M_b = 2kL\,h, \qquad T = \tfrac{kd^2}{2}\,\theta,$$

so bending is driven purely by plunge and torsion purely by pitch before they merge in the equivalent-moment formula.

## How h and theta are obtained

Both responses are driven by the §25.341 critical 1-cosine gust, reusing existing functions rather than re-deriving anything:

- **Constant k (left panel):** `impulse_function_generator('model_file', wing_nl.mat, 'plot', false)` supplies the gust family and critical column tuned against the `build_ss_nl` linearization; `lsim` on (A, B, C, D) gives $h(t),\theta(t)$. With $k_{nl}=0$ in the moment evaluation, the tangent stiffness is constant.
- **Varying k (right panel):** `nonlinear_response(4)` integrates the full cubic-spring EOM under the same gust type and returns $h(t),\theta(t)$; the full $F_i = k_i\delta_i + k_{nl,i}\delta_i^3$ is used, so the effective stiffness varies with deflection.
