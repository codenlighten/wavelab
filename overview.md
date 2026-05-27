Below is a strong equation map for the **simulated-light / wave-based molecular geometry engine**. I would organize them into layers: **wave propagation, optical medium encoding, molecular field mapping, scoring, and learning/optimization.**

---

# 1. Core Wave Equation

The foundation is the scalar wave equation:

$$
\frac{\partial^2 u}{\partial t^2} = c^2 \nabla^2 u
$$

In 2D:

$$
\frac{\partial^2 u}{\partial t^2} = c^2 \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right)
$$

In 3D:

$$
\frac{\partial^2 u}{\partial t^2} = c^2 \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 u}{\partial z^2} \right)
$$

Where $u(x,y,z,t)$ is the wave amplitude field.

---

# 2. Discrete FDTD Update Equation

For 2D simulation:

$$
u_{i,j}^{t+1} = 2 u_{i,j}^{t} - u_{i,j}^{t-1} + \lambda^2 \left( u_{i+1,j}^{t} + u_{i-1,j}^{t} + u_{i,j+1}^{t} + u_{i,j-1}^{t} - 4 u_{i,j}^{t} \right)
$$

Where:

$$
\lambda = \frac{c \, \Delta t}{\Delta x}
$$

For 3D:

$$
u_{i,j,k}^{t+1} = 2 u_{i,j,k}^{t} - u_{i,j,k}^{t-1} + \lambda^2 \left( u_{i+1,j,k}^{t} + u_{i-1,j,k}^{t} + u_{i,j+1,k}^{t} + u_{i,j-1,k}^{t} + u_{i,j,k+1}^{t} + u_{i,j,k-1}^{t} - 6 u_{i,j,k}^{t} \right)
$$

This is the computational heartbeat of the simulated-light engine.

> **Note:** the form above assumes constant $c$. With a spatially varying medium $c(x)$ (see §9), bake $c(x)^2$ into a per-cell $\lambda^2_{i,j,k}$ inside the parenthesized stencil, or work from the more general form $\tfrac{1}{c^2(x)}\,\partial_t^2 u = \nabla^2 u$.

---

# 3. Stability / Courant Condition

To prevent the simulation from blowing up numerically:

For 2D:

$$
\frac{c \, \Delta t}{\Delta x} \leq \frac{1}{\sqrt{2}}
$$

For 3D:

$$
\frac{c \, \Delta t}{\Delta x} \leq \frac{1}{\sqrt{3}}
$$

More generally:

$$
\Delta t \leq \frac{\Delta x}{c \sqrt{d}}
$$

where $d$ is the number of spatial dimensions.

---

# 4. Harmonic Wave Source

A clean laser-like source can be modeled as:

$$
S(t) = A \sin(2\pi f t + \phi)
$$

or:

$$
S(t) = A \cos(\omega t + \phi), \qquad \omega = 2\pi f
$$

For a Gaussian pulse:

$$
S(t) = A \, e^{-\frac{(t - t_0)^2}{2 \sigma_t^2}} \sin(2\pi f t)
$$

This allows us to probe the molecular geometry with different frequencies.

---

# 5. Damped Wave Equation

Realistic simulations should include dissipation:

$$
\frac{\partial^2 u}{\partial t^2} + \gamma \frac{\partial u}{\partial t} = c^2 \nabla^2 u
$$

Discrete approximation:

$$
u^{t+1} = (2 - \gamma \, \Delta t) \, u^{t} - (1 - \gamma \, \Delta t) \, u^{t-1} + \lambda^2 \, \nabla_d^2 u^{t}
$$

Where $\gamma$ is the damping coefficient. This helps model absorption, lossy regions, and boundary decay.

---

# 6. Molecular Occupancy Field

Each atom can be projected into a grid as a density field:

$$
\rho(x) = \sum_{a=1}^{N} \exp\!\left( -\frac{\lVert x - r_a \rVert^2}{2 \sigma_a^2} \right)
$$

Where:

* $r_a$ is the atom position.
* $\sigma_a$ controls atom radius/spread.
* $\rho(x)$ becomes the molecular density field.

This turns atoms into a smooth geometry.

---

# 7. Hard Obstacle Boundary Field

A simple binary molecular mask:

$$
B(x) =
\begin{cases}
1, & \rho(x) \geq \tau \\
0, & \rho(x) < \tau
\end{cases}
$$

Where:

* $B(x) = 1$ means blocked or occupied.
* $B(x) = 0$ means open medium.
* $\tau$ is the density threshold.

This creates a protein-pocket obstacle map.

---

# 8. Soft Boundary / Permeability Field

Instead of hard obstacles, use a soft attenuation field:

$$
\alpha(x) = \alpha_0 + \alpha_1 \, \rho(x)
$$

Then the wave update includes local damping:

$$
\frac{\partial^2 u}{\partial t^2} + \alpha(x) \frac{\partial u}{\partial t} = c^2(x) \, \nabla^2 u
$$

Dense atomic regions absorb or slow the wave more strongly.

---

# 9. Refractive Index Field

To make the molecule behave like an optical medium:

$$
n(x) = n_0 + \beta_\rho \, \rho(x) + \beta_q \, \lvert Q(x) \rvert + \beta_h \, H(x)
$$

Where:

* $n(x)$ is the refractive index.
* $\rho(x)$ is atom density.
* $Q(x)$ is charge field.
* $H(x)$ is hydrophobicity field.
* $\beta_\rho, \beta_q, \beta_h$ are tunable weights.

Wave speed becomes:

$$
c(x) = \frac{c_0}{n(x)}
$$

This lets molecular chemistry alter the path of the simulated wave.

---

# 10. Charge Field

A smoothed charge field:

$$
Q(x) = \sum_{a=1}^{N} q_a \, \exp\!\left( -\frac{\lVert x - r_a \rVert^2}{2 \sigma_q^2} \right)
$$

Where $q_a$ is the partial charge of atom $a$. This allows polar regions to influence the simulated medium.

---

# 11. Hydrophobicity Field

A hydrophobicity field:

$$
H(x) = \sum_{a=1}^{N} h_a \, \exp\!\left( -\frac{\lVert x - r_a \rVert^2}{2 \sigma_h^2} \right)
$$

Where $h_a$ represents atom-level hydrophobic contribution. This lets hydrophobic pockets affect wave propagation differently from polar pockets.

---

# 12. Electrostatic Potential Field

Classical electrostatic potential:

$$
\Phi(x) = \sum_{a=1}^{N} \frac{q_a}{4 \pi \varepsilon_0 \, \lVert x - r_a \rVert}
$$

With smoothing:

$$
\Phi(x) = \sum_{a=1}^{N} \frac{q_a}{4 \pi \varepsilon_0 \, \sqrt{\lVert x - r_a \rVert^2 + \epsilon^2}}
$$

This avoids singularities near atom centers.

---

# 13. Lennard-Jones Field

The local steric/dispersion field around the molecule:

$$
V_{LJ}(x) = \sum_{a=1}^{N} 4 \varepsilon_a \left[ \left( \frac{\sigma_a}{\lVert x - r_a \rVert} \right)^{12} - \left( \frac{\sigma_a}{\lVert x - r_a \rVert} \right)^{6} \right]
$$

Smoothed version:

$$
V_{LJ}(x) = \sum_{a=1}^{N} 4 \varepsilon_a \left[ \left( \frac{\sigma_a}{\sqrt{\lVert x - r_a \rVert^2 + \epsilon^2}} \right)^{12} - \left( \frac{\sigma_a}{\sqrt{\lVert x - r_a \rVert^2 + \epsilon^2}} \right)^{6} \right]
$$

This can be used to influence the wave medium or create exclusion zones.

---

# 14. Wave Energy Density

At each point:

$$
E(x, t) = \frac{1}{2} \left( \frac{\partial u}{\partial t} \right)^{2} + \frac{1}{2} c^2 \, \lVert \nabla u \rVert^{2}
$$

Discrete approximation:

$$
E_{i,j,k}^{t} = \frac{1}{2} \left( \frac{u_{i,j,k}^{t} - u_{i,j,k}^{t-1}}{\Delta t} \right)^{2} + \frac{1}{2} c^2 \, \lVert \nabla_d u_{i,j,k}^{t} \rVert^{2}
$$

This measures where the field stores energy.

---

# 15. Total Field Energy

$$
E_{\text{total}}(t) = \sum_x E(x, t)
$$

Energy retained after ligand insertion:

$$
R_E = \frac{E_{\text{total}}^{\text{pocket}+\text{ligand}}}{E_{\text{total}}^{\text{pocket}}}
$$

This can become a wave-stability score.

---

# 16. Scattering Loss

If $E_{\text{in}}$ is injected energy and $E_{\text{out}}$ is measured exiting energy:

$$
L_{\text{scatter}} = 1 - \frac{E_{\text{out}}}{E_{\text{in}}}
$$

A ligand that causes excessive scattering may be geometrically mismatched.

---

# 17. Reflection Coefficient

$$
R = \frac{E_{\text{reflected}}}{E_{\text{incident}}}
$$

Transmission coefficient:

$$
T = \frac{E_{\text{transmitted}}}{E_{\text{incident}}}
$$

Absorption estimate:

$$
A = 1 - R - T
$$

For a stable pocket-ligand geometry, these values can become part of the fingerprint.

---

# 18. Wave Field Difference Score

Compare the empty pocket to the ligand-filled pocket:

$$
D_{\text{wave}} = \lVert W_{\text{pocket}} - W_{\text{pocket}+\text{ligand}} \rVert_2
$$

Expanded:

$$
D_{\text{wave}} = \sqrt{ \sum_x \left( W_{\text{pocket}}(x) - W_{\text{pocket}+\text{ligand}}(x) \right)^{2} }
$$

Similarity score:

$$
S_{\text{wave}} = \exp(-\eta \, D_{\text{wave}})
$$

Where $\eta$ controls sensitivity.

---

# 19. Phase Coherence Score

For complex wave representation:

$$
\psi(x, t) = A(x, t) \, e^{i \phi(x, t)}
$$

Phase coherence:

$$
C_\phi = \left\lvert \frac{1}{M} \sum_{m=1}^{M} e^{i \phi_m} \right\rvert
$$

Where:

* $C_\phi \approx 1$ means highly coherent.
* $C_\phi \approx 0$ means phase disorder.

This is useful for detecting chaotic scattering.

---

# 20. Interference Intensity

If two waves overlap:

$$
I = \lvert u_1 + u_2 \rvert^{2}
$$

Expanded:

$$
I = \lvert u_1 \rvert^{2} + \lvert u_2 \rvert^{2} + 2 \, \lvert u_1 \rvert \, \lvert u_2 \rvert \cos(\Delta \phi)
$$

Constructive interference occurs when:

$$
\Delta \phi = 2 \pi k
$$

Destructive interference occurs when:

$$
\Delta \phi = (2k + 1) \pi
$$

This lets us measure pocket resonance and destructive mismatch.

---

# 21. Standing Wave / Resonance Condition

For a cavity of approximate length $L$:

$$
L = \frac{m \lambda}{2}
$$

where $m$ is an integer mode number. Since:

$$
\lambda = \frac{c}{f}
$$

resonant frequencies are:

$$
f_m = \frac{m c}{2 L}
$$

This can be generalized to molecular cavities:

$$
f_m \sim \frac{m \, c_{\text{eff}}}{2 \, L_{\text{pocket}}}
$$

The point is not literal visible-light resonance, but geometry-driven modal behavior.

---

# 22. Spectral Fingerprint

Take the Fourier transform of the wave response:

$$
\hat{u}(x, \omega) = \int u(x, t) \, e^{-i \omega t} \, dt
$$

Frequency energy:

$$
P(\omega) = \sum_x \lvert \hat{u}(x, \omega) \rvert^{2}
$$

The vector:

$$
F_{\text{spectral}} = \big[\, P(\omega_1),\, P(\omega_2),\, \dots,\, P(\omega_m) \,\big]
$$

becomes a pocket fingerprint.

---

# 23. Entropy of Wave Energy

Normalize wave energy:

$$
p_x = \frac{E(x)}{\sum_x E(x)}
$$

Then compute entropy:

$$
H_E = -\sum_x p_x \log p_x
$$

Low entropy may indicate focused energy. High entropy may indicate diffuse scattering. A useful score:

$$
S_{\text{focus}} = 1 - \frac{H_E}{\log M}
$$

where $M$ is the number of grid cells.

---

# 24. Hotspot Concentration Score

Energy concentration in binding-region hotspots:

$$
C_{\text{hotspot}} = \frac{\sum_{x \in \Omega_{\text{hot}}} E(x)}{\sum_{x \in \Omega} E(x)}
$$

Where:

* $\Omega$ is the full simulation domain.
* $\Omega_{\text{hot}}$ is a pocket-relevant region.

This measures whether the wave focuses into meaningful pocket zones.

---

# 25. Ligand Occlusion Score

If the ligand blocks important wave pathways:

$$
O_{\text{ligand}} = \frac{\sum_{x \in \Omega_L} E_{\text{empty}}(x)}{\sum_{x \in \Omega} E_{\text{empty}}(x)}
$$

Where $\Omega_L$ is the ligand-occupied region. This measures whether the ligand occupies high-importance wave regions.

---

# 26. Shape Complementarity via Wave Fields

Let $W_P(x)$ be the pocket wave field and $W_L(x)$ be the ligand wave field. A complementarity score:

$$
S_{\text{comp}} = \frac{\sum_x W_P(x) \, W_L(x)}{\sqrt{\sum_x W_P(x)^{2}} \, \sqrt{\sum_x W_L(x)^{2}}}
$$

This is cosine similarity between wave fingerprints.

---

# 27. Hybrid Docking Score

A final scoring equation:

$$
S_{\text{total}} = w_1 S_{LJ} + w_2 S_{\text{elec}} + w_3 S_{\text{Hbond}} + w_4 S_{\text{hydro}} + w_5 S_{\text{wave}} + w_6 S_{\text{entropy}}
$$

Or in binding-energy style:

$$
\Delta G_{\text{approx}} = \Delta G_{LJ} + \Delta G_{\text{elec}} + \Delta G_{\text{solv}} + \Delta G_{\text{entropy}} + \Delta G_{\text{wave}}
$$

Where the wave term is a learned or calibrated geometric correction.

---

# 28. Learned Wave Correction Term

A machine learning model can learn the wave contribution:

$$
\Delta G_{\text{wave}} = f_\theta \!\left( F_{\text{wave}},\, F_{\text{classical}},\, F_{\text{ligand}},\, F_{\text{pocket}} \right)
$$

Where:

* $F_{\text{wave}}$ contains optical features.
* $F_{\text{classical}}$ contains LJ/electrostatic terms.
* $F_{\text{ligand}}$ contains molecular descriptors.
* $F_{\text{pocket}}$ contains pocket descriptors.

Training objective:

$$
\mathcal{L}(\theta) = \frac{1}{N} \sum_{i=1}^{N} \left( \widehat{\Delta G}_i - \Delta G_i \right)^{2}
$$

---

# 29. Multi-Frequency Wave Score

Instead of using one source frequency, scan many:

$$
S_{\text{wave}}^{\text{multi}} = \sum_{f \in F} \alpha_f \, S_{\text{wave}}(f)
$$

Where:

* $F$ is a frequency set.
* $\alpha_f$ weights each frequency.
* Different frequencies reveal different geometric scales.

This is important because one wavelength may reveal broad cavities while another reveals narrow channels.

---

# 30. Wavelength-to-Geometry Relationship

The wavelength should correspond to the structural scale being probed:

$$
\lambda \approx k \, L_{\text{feature}}, \qquad f = \frac{c}{\lambda}
$$

For a pocket feature of size $L_{\text{feature}}$, choose frequencies that satisfy:

$$
f \approx \frac{c}{k \, L_{\text{feature}}}
$$

This lets us tune simulated light to detect different molecular features.

---

# 31. Boundary Conditions

Hard reflective (Dirichlet) boundary:

$$
u(x, t) = 0 \quad \text{on} \quad \partial \Omega
$$

Neumann reflective boundary:

$$
\frac{\partial u}{\partial n} = 0 \quad \text{on} \quad \partial \Omega
$$

Absorbing boundary approximation (1st-order Mur):

$$
\frac{\partial u}{\partial t} + c \, \frac{\partial u}{\partial n} = 0
$$

The absorbing boundary prevents artificial box reflections. Note that 1st-order Mur reflects a few percent at oblique incidence — for high-fidelity docking scores, upgrade to a perfectly matched layer (PML).

---

# 32. Pocket as Operator

This is the most elegant mathematical framing.

The pocket is an operator:

$$
\mathcal{P} : S_{\text{in}} \to S_{\text{out}}
$$

The ligand-modified pocket is:

$$
\mathcal{P}_L : S_{\text{in}} \to S_{\text{out}}^{L}
$$

Then ligand fitness can be measured by how it transforms the operator:

$$
D_{\mathcal{P}} = \lVert \mathcal{P} - \mathcal{P}_L \rVert
$$

Or by whether it transforms the output toward known binder patterns:

$$
S_{\text{binder}} = \operatorname{sim}\!\left( S_{\text{out}}^{L},\, S_{\text{known binder}} \right)
$$

This is beautiful for the thesis because it turns the molecule into a computational object.

---

# 33. Green's Function / Impulse Response

The wave response to a point impulse $G(x, t; x_0)$ satisfies:

$$
\left( \frac{\partial^{2}}{\partial t^{2}} - c^{2} \nabla^{2} \right) G(x, t; x_0) = \delta(x - x_0) \, \delta(t)
$$

The field from a source is:

$$
u(x, t) = \int G(x, t; x_0) \, S(x_0, t) \, dx_0
$$

In practice, the impulse response of the pocket becomes a structural fingerprint.

---

# 34. Correlation with Known Binder Signature

If known binders generate wave signatures $F_1, F_2, \dots, F_K$, create a prototype:

$$
\bar{F}_{\text{binder}} = \frac{1}{K} \sum_{k=1}^{K} F_k
$$

Then compare a new ligand:

$$
S_{\text{known}} = \frac{F_{\text{new}} \cdot \bar{F}_{\text{binder}}}{\lVert F_{\text{new}} \rVert \, \lVert \bar{F}_{\text{binder}} \rVert}
$$

This allows empirical calibration.

---

# 35. Final Research Equation

The most complete form:

$$
\widehat{\Delta G} = f_\theta \!\left( S_{LJ},\, S_{\text{elec}},\, S_{\text{Hbond}},\, S_{\text{hydro}},\, S_{\text{shape}},\, S_{\text{wave}},\, C_\phi,\, H_E,\, R_E,\, L_{\text{scatter}},\, F_{\text{spectral}} \right)
$$

Or as a clean thesis equation:

$$
\boxed{\; \widehat{\Delta G} = \Delta G_{\text{classical}} + \lambda_w \, \Delta G_{\text{wave}} \;}
$$

Where:

$$
\Delta G_{\text{wave}} = g\!\left( \text{interference},\, \text{phase},\, \text{energy},\, \text{scattering},\, \text{resonance},\, \text{geometry} \right)
$$

---

# The Core Equation Stack

For the thesis, I would highlight these as the essential set:

| Layer             | Equation |
| ----------------- | -------- |
| Wave propagation  | $\partial_t^{2} u = c^{2} \nabla^{2} u$ |
| FDTD update       | $u^{t+1} = 2 u^{t} - u^{t-1} + \lambda^{2} \nabla_d^{2} u^{t}$ |
| Stability         | $\Delta t \leq \dfrac{\Delta x}{c \sqrt{d}}$ |
| Molecular density | $\rho(x) = \sum_a e^{-\lVert x - r_a \rVert^{2} / (2 \sigma_a^{2})}$ |
| Refractive index  | $n(x) = n_0 + \beta_\rho \rho(x) + \beta_q \lvert Q(x) \rvert + \beta_h H(x)$ |
| Local wave speed  | $c(x) = c_0 / n(x)$ |
| Energy density    | $E = \tfrac{1}{2} u_t^{2} + \tfrac{1}{2} c^{2} \lVert \nabla u \rVert^{2}$ |
| Phase coherence   | $C_\phi = \left\lvert \tfrac{1}{M} \sum_m e^{i \phi_m} \right\rvert$ |
| Energy entropy    | $H_E = -\sum_x p_x \log p_x$ |
| Wave similarity   | $S_{\text{wave}} = \exp(-\eta \lVert W_P - W_{P+L} \rVert)$ |
| Hybrid score      | $\widehat{\Delta G} = \Delta G_{\text{classical}} + \lambda_w \Delta G_{\text{wave}}$ |

This gives you a serious mathematical foundation for the simulated-light engine.
