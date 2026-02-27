# FreeTTES // Numerical Methods & Stability

This document describes the numerical methods, discretization choices, and
stability mechanisms used in the stratified TES tank model. It is intended
for engineers concerned with robustness, conservation, and long-term
simulation reliability.

---

## 1. Time Integration Strategy

### 1.1 External vs Internal Timesteps

The model is advanced using an **external timestep** $\Delta t$ provided by the
calling system (e.g., 900s).

To ensure numerical stability, $\Delta t$ is subdivided into **internal
substeps** $\Delta t_{\mathrm{sub}}$:

$$
\Delta t = N_{\mathrm{sub}} \, \Delta t_{\mathrm{sub}}
$$

with $\Delta t_{\mathrm{sub}}$ = 60 s

---

### 1.2 Rationale for Substepping

Substepping prevents:

- excessive entrainment in a single step,
- overshoot in impulse generation,
- numerical diffusion from conduction,
- instability from strong inflows.

The approach trades minimal additional cost for significant robustness gains.

---

## 2. Spatial Discretization (Vertical Grid)

### 2.1 Adaptive Layering

The water column is discretized into layers of variable height $\Delta h_i$.

Layer sizes are dynamically adjusted to:

- resolve steep temperature gradients (thermocline),
- avoid excessively thin layers,
- maintain computational efficiency.

---

### 2.2 Grid Refinement (Split)

Layers are split when:

- $\Delta h_i > \Delta h_{\max}$, or
- temperature gradients exceed a threshold.

When splitting:

- mass is divided proportionally,
- temperature is conserved by enthalpy balance,
- impulse $I$ and mixing $M$ are copied to child layers if active.

---

### 2.3 Grid Coarsening (Merge)

Adjacent layers are merged when:

- $\Delta h_i < \Delta h_{\min}$, or
- temperature difference is negligible.

During merging:

$$
T_{\mathrm{new}} =
\frac{m_1 h(T_1) + m_2 h(T_2)}{m_1 + m_2}
$$

Impulse and mixing proxies are averaged or conserved according to mass
weighting.

---

## 3. Inflow / Outflow Discretization

- Inflow volumes are distributed incrementally across substeps
- Outflow is extracted proportionally from affected layers
- Enthalpy consistency is enforced at each substep

This avoids abrupt volume or energy changes.

---

## 4. Inversion and Mixing Stability Controls

### 4.1 Entrainment Limits

Entrained volume per interaction is capped:

$$
\Delta V_e \le 0.99 \, V_{\mathrm{ambient}}
$$

This prevents full homogenization in a single interaction.

---

### 4.2 Impulse Bounding

Impulse is constrained:

$$
I \ge 0
$$

Negative impulse-squared values are clipped to zero, preventing numerical
artifacts.

---

### 4.3 Stopping Criteria

Inversion iterations terminate when:

- no density inversion remains,
- impulse decays to zero,
- domain boundaries are reached.

This ensures termination and stability.

---

## 5. Heat Conduction Solver

### 5.1 Governing Equation

The one-dimensional heat equation:

$$
\rho c_p \frac{\partial T}{\partial t} =
\frac{\partial}{\partial z}
\left(
\lambda \frac{\partial T}{\partial z}
\right)
$$

is discretized implicitly.

---

### 5.2 Implicit Scheme

The implicit scheme yields a tridiagonal linear system:

$$
\mathbf{A} \, \mathbf{T}^{n+1} = \mathbf{b}
$$

This formulation is unconditionally stable with respect to timestep size.

---

### 5.3 TDMA Solver

The tridiagonal system is solved using the Thomas algorithm (TDMA)

---

## 6. Wall and Foundation Coupling

### 6.1 Wall Capacity Nodes

Wall nodes are advanced using a lumped-capacitance ODE:

$$
C_w \frac{dT_w}{dt} =
\dot Q_{\mathrm{water}} - \dot Q_{\mathrm{amb}}
$$

Integrated consistently with substepping.

---

### 6.2 Foundation Layers

Foundation temperatures are updated via conduction coupling to the bottom water
layer, using the same implicit framework.

---

## 7. Conservation Guarantees

The numerical scheme enforces:

- mass conservation (except for inflow/outflow),
- internal enthalpy conservation (except losses),
- monotonic density stratification after inversion resolution.

---
