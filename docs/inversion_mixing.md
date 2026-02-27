# FreeTTES // Inversion & Mixing Algorithm

## 1. Purpose of the Inversion Algorithm

Discrete stratified models can develop gravitationally unstable states due to:

* inflow placement,
* horizontal mixing,
* numerical diffusion from conduction,
* grid split/merge operations.

A real fluid cannot sustain such states. The inversion algorithm:

* detects unstable stratification,
* resolves it via buoyancy-driven motion,
* mixes fluid locally until stability is restored,
* conserves mass and enthalpy (internally).

---

## 2. Conceptual Model: Moving Plug

The algorithm operates on the concept of a moving plug :

* a contiguous parcel of fluid that is buoyant relative to its surroundings,
* capable of moving vertically through adjacent layers,
* entraining ambient fluid as it moves.

Two cases are handled:

* **Rising plug** : hot, light fluid below colder fluid
* **Falling plug** : cold, dense fluid above warmer fluid

---

## 3. High-Level Algorithm Flow

1. Scan for density inversions
2. Identify a moving plug and its direction
3. While the plug remains buoyant:
   * interact with the adjacent layer
   * entrain part of that layer
   * update impulse and mixing strength
   * rewrite the two layers
   * advance the plug
4. Stop when stable stratification is restored

---

## 4. Schematic: Single Interaction Step

```
Before interaction:

z ↑
──────────────
[r]  T_r, dh_r   ← resting ambient layer
──────────────
[b]  T_b, dh_b   ← moving plug
──────────────

After interaction:

z ↑
──────────────
[b'] mixed plug  ← plug advanced upward/downward
──────────────
[r'] mixed rest  ← remaining ambient fluid
──────────────
```

The plug advances by rewriting layer pairs, not by translating geometry.

---

## 5. Detailed Step-by-Step Mechanics

### 5.1 Inversion Detection

Adjacent layers *i* and *j* are unstable if:

$$
\rho(T_i) < \rho(T_j) \quad (z_i < z_j)
$$

This condition defines the direction of motion.

---

### 5.2 Plug Initialization

The lower (rising case) or upper (falling case) layer becomes the initial moving plug with:

* temperature (T_b)
* height (dh_b)
* impulse (I_b)
* mixing proxy (M_b)

---

### 5.3 Entrainment From Ambient Layer

As the plug encounters an adjacent layer of height (dh_r), an entrained volume:

$$
\Delta V_e = f_e,dh_r
$$

is mixed into the plug.

* (f_e) is an empirical entrainment coefficient

### 5.4 Impulse Update

Impulse is updated using a work–energy balance:

[

$$
I_{new}^2 = I_{old}^2 \Bigl(1 - 2\ln\frac{m_{new}}{m_{old}}\Bigr)-2g\frac{\Delta \rho}{\rho}
dh_r
$$

Where:

* the logarithmic term represents dilution by entrainment ,
* the buoyancy term represents conversion between potential and kinetic energy .

---

### 5.5 Mixing Proxy Update

The mixing proxy is diluted as mass increases:

$$
M_{new} = M_{old}\frac{m_{old}}{m_{new}}
$$

This ensures jet-induced mixing decays naturally.

---

### 5.6 Optional Boundary Exchange (Tail Adjustment)

The plug may partially exchange fluid back into the ambient layer based on an impulse-determined boundary temperature.

This step:

* limits over-penetration,
* shapes a smooth thermal transition zone,
* prevents nonphysical sharp gradients.

---

### 5.7 Layer Rewrite and Proxy Transport

After mixing:

* the plug mixture replaces the ambient layer position,
* the remaining ambient fluid replaces the former plug position,
* impulse `I` and mixing `M` are **swapped **between the two positions.

This implements parcel-based transport of activity variables.

---

## 6. Stopping Criteria

The iteration stops when **any** of the following is met:

* no remaining density inversion
* impulse decays to zero
* tank boundary is reached

At this point, stratification is stable.

---

## 7. Conservation Properties

The algorithm enforces:

* mass conservation
* internal enthalpy conservation
* monotonic density stratification

---

## 8. Engineering Interpretation

The inversion algorithm is a  reduced-order plume solver:

* similar in spirit to integral plume models,
* adapted to a discrete 1D layered grid,
* suitable for long-term TTES simulation.

It captures:

* buoyant rise/fall,
* entrainment-driven dilution,
* decay of impulse,
* restoration of stable stratification.

---

## 9. Why This Approach Is Used

Compared to alternatives:

* Global mixing → destroys stratification unrealistically
* No inversion handling→ allows nonphysical states
* CFD → computationally infeasible for system simulation

This approach balances  physical realism and computational efficiency .
