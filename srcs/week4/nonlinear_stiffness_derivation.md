# Stiffness Model — From Constant Uncoupled `K` to Nonlinear Coupled `g(q)`

This documents the structural-restoring derivation used in `build_ss_nl.m` (matrix assembly) and `nonlinear_response.m` (`spring_force`). Each stage below was checked symbolically; the final restoring vector and the linear reduction match the code exactly.

## 0. Coordinates and spring geometry

Typical-section 2-DOF, generalized coordinates

$$q = \begin{bmatrix} h \\ \theta \end{bmatrix}, \qquad h = \text{plunge [m]}, \quad \theta = \text{pitch [rad]}.$$

Two discrete springs sit at moment arms $r_1, r_2$ from the reference (elastic) axis. With separation $d$ and the reference at the midpoint,

$$r_1 = +\tfrac{d}{2}, \qquad r_2 = -\tfrac{d}{2}.$$

Keeping the **geometry** linear (small angle), the vertical deflection seen by spring $i$ is

$$\delta_i = h + r_i\,\theta.$$

Only the spring **constitutive law** will become nonlinear later; the kinematics $\delta_i(h,\theta)$ stay linear throughout.

## 1. Constant, uncoupled `K` (the original `build_ss.m`)

Linear, identical springs $k_1 = k_2 = k$. Potential energy:

$$V = \tfrac{1}{2}k_1\delta_1^2 + \tfrac{1}{2}k_2\delta_2^2.$$

The stiffness matrix is the Hessian $K_{ij} = \partial^2 V / \partial q_i \partial q_j$. Differentiating,

$$\frac{\partial V}{\partial h} = \sum_i k_i\delta_i, \qquad \frac{\partial V}{\partial \theta} = \sum_i k_i r_i \delta_i,$$

so in general

$$K = \begin{bmatrix} \sum_i k_i & \sum_i k_i r_i \\[2pt] \sum_i k_i r_i & \sum_i k_i r_i^2 \end{bmatrix}.$$

Substituting $k_1=k_2=k$, $r_1=+d/2$, $r_2=-d/2$:

$$K_{hh} = 2k, \qquad K_{h\theta} = k\tfrac{d}{2} + k\!\left(-\tfrac{d}{2}\right) = 0, \qquad K_{\theta\theta} = k\tfrac{d^2}{4}+k\tfrac{d^2}{4} = k\tfrac{d^2}{2},$$

$$\boxed{\,K = \begin{bmatrix} 2k & 0 \\ 0 & \tfrac{1}{2}k d^2 \end{bmatrix}\,}$$

This is exactly the `Ks` in the original `build_ss.m`. The off-diagonal vanishes **because the springs are identical** — symmetry, not an assumption baked into the form.

## 2. Constant, coupled `K` (allow `k1 ≠ k2`)

Still linear springs, but now distinct. The same general $K$ gives

$$K = \begin{bmatrix} k_1 + k_2 & (k_1 - k_2)\tfrac{d}{2} \\[2pt] (k_1 - k_2)\tfrac{d}{2} & (k_1 + k_2)\tfrac{d^2}{4} \end{bmatrix}.$$

The coupling term $K_{h\theta} = (k_1 - k_2)\,d/2$ is now nonzero: an asymmetric pair of springs ties plunge and pitch together. $K$ stays symmetric (the structure is conservative). This is the matrix literal assembled in `build_ss_nl.m`.

## 3. Nonlinear, coupled — restoring becomes a vector `g(q)`

Give each spring a cubic term (hardening for $k_{nl,i}>0$):

$$F_i(\delta_i) = k_i\,\delta_i + k_{nl,i}\,\delta_i^3, \qquad V = \sum_i\left(\tfrac{1}{2}k_i\delta_i^2 + \tfrac{1}{4}k_{nl,i}\delta_i^4\right).$$

The generalized restoring force is $g(q) = \partial V/\partial q$. Using $\partial\delta_i/\partial h = 1$ and $\partial\delta_i/\partial\theta = r_i$,

$$g(q) = \begin{bmatrix} \sum_i F_i \\[2pt] \sum_i r_i F_i \end{bmatrix} = \begin{bmatrix} F_1 + F_2 \\ r_1 F_1 + r_2 F_2 \end{bmatrix}.$$

This is precisely `spring_force(q, P)`. The key point: the entries are assembled from each spring's **total** force $F_i$ (linear + cubic together) — the sum is what matters, not the terms separately. Once $k_{nl}\neq 0$, $g$ is no longer expressible as $K q$ for any constant $K$; the matrix picture is gone and only the vector field remains.

## 4. Tangent stiffness — recovering the matrix locally, and the `k1(t)` plot

Linearizing $g$ about an operating point recovers a matrix, the tangent stiffness $K_t(q) = \partial g/\partial q$. Each spring contributes an instantaneous stiffness

$$k_{i,t} = \frac{dF_i}{d\delta_i} = k_i + 3\,k_{nl,i}\,\delta_i^2,$$

so

$$K_t(q) = \begin{bmatrix} \sum_i k_{i,t} & \sum_i r_i\,k_{i,t} \\[2pt] \sum_i r_i\,k_{i,t} & \sum_i r_i^2\,k_{i,t} \end{bmatrix}.$$

Same structure as Stage 2, with each $k_i$ replaced by its deflection-dependent $k_{i,t}$.

**Two consequences the code relies on:**

- At $q = 0$ all $\delta_i = 0 \Rightarrow k_{i,t} = k_i$, so $K_t(0)$ equals the linear $K$ of Stage 2. That $q=0$ tangent is the $A,B,C,D$ linearization exported by `build_ss_nl.m`, and it reduces to the original model when $k_{nl}=0$.
- The 4th panel of `nonlinear_response`, $k_1(t) = k_1 + 3k_{nl,1}\delta_1(t)^2$, is exactly $k_{1,t}$ traced along the trajectory — flat at $k_1$ when linear, breathing once cubic.

A nice check on the nonlinearity: even with **symmetric** springs ($k_1=k_2$, $k_{nl,1}=k_{nl,2}=k_{nl}$), the tangent coupling no longer vanishes,

$$K_{t,\,h\theta} = 3\,k_{nl}\,d^2\,h\,\theta,$$

so the cubic springs introduce an amplitude- and motion-dependent plunge–pitch coupling that the linear model can never show.

## 5. Equation of motion (closing the loop)

The structural restoring $g(q)$ enters the full aeroelastic EOM with the aerodynamic terms kept linear:

$$M\,\ddot q + C\,\dot q + K_a\,q + g(q) = F_{\text{gust}}\,u(t),$$

with $C = C_{\text{struct}} + C_a$ and $K_a$, $C_a$, $F_{\text{gust}}$ the aero operators. This is the right-hand side integrated in `nonlinear_response.m`. Setting $k_{nl,1}=k_{nl,2}=0$ collapses $g(q)\to K q$ and returns the linear model

$$M\,\ddot q + C\,\dot q + (K + K_a)\,q = F_{\text{gust}}\,u(t),$$

which is the built-in validation case.

---

**Symbolic verification summary.** `g(q)` from $\partial V/\partial q$ matches the assembled $[F_1+F_2;\; r_1F_1+r_2F_2]$; $K_t(0)$ matches the general linear $K$; the substitution $k_1=k_2=k,\,k_{nl}=0,\,r_{1,2}=\pm d/2$ returns $\mathrm{diag}(2k,\,\tfrac12 k d^2)$; and the symmetric-spring tangent coupling evaluates to $3k_{nl}d^2 h\theta$.
