# FreeTTES // State Definition & Data Structures

This document defines all internal state variables and data structures used by the FreeTTES model. It is written for engineers and developers who need to understand, extend, or validate the model.

The descriptions are based strictly onexecutable behavior , not inline code comments.

---

## 1. Global State Overview

At each timestep, the full physical state of the storage system is represented by three coupled state objects:

1. `Speicherzustand` – water column (primary thermal storage)
2. `Fundamentzustand` – foundation temperature profile
3. `Kapazitaeten` – wall capacity cells

These states are saved to disk between timesteps and reloaded on the next call to `main()`.

---

## 2. Speicherzustand (Water Column State)

### 2.1 Data Structure

```
Speicherzustand : dict[float, list]
Speicherzustand[h] = [T, dh, I, M]
```

* **Key `h`** : height-like coordinate [m]
* Represents the vertical position of the *bottom* of a layer (layer=cell is used interchangeably)
* Keys are sorted to reconstruct vertical order
* **Value `[T, dh, I, M]`** : state vector of a single layer

---

### 2.2 Layer State Vector

| Index   | Symbol | Name          | Units | Meaning                           |
| ------- | ------ | ------------- | ----- | --------------------------------- |
| `[0]` | `T`  | Temperature   | °C   | Mean temperature of the layer     |
| `[1]` | `dh` | Layer height  | m     | Vertical thickness of the layer   |
| `[2]` | `I`  | Impulse proxy | m/s   | Vertical penetration capability   |
| `[3]` | `M`  | Mixing proxy  | m/s   | Jet/shear-induced mixing strength |

---

### 2.3 Derived Quantities

For each layer  *i* :

* Volume: $V_i = A\,\Delta h_i$
* Mass: $m_i = \rho(T_i)\,V_i$

where *A* is the tank cross-sectional area.

---

## 3. Physical Interpretation of State Variables

### 3.1 Temperature `T`

* Represents a **fully mixed control volume**
* Updated by:
  * inflow/outflow enthalpy exchange
  * inversion mixing
  * horizontal mixing
  * thermal conduction

---

### 3.2 Layer Height `dh`

* Encodes the **geometry and volume** of the control volume
* Changes due to:
  * inflow/outflow volume addition/removal
  * entrainment during inversion
  * grid refinement (split/merge)

Mass and energy conservation are enforced whenever `dh` is modified.

---

### 3.3 Impulse Proxy `I` (`[2]`)

#### Definition

`I` is a **velocity-like scalar** representing the ability of a fluid parcel (moving/inserted cell) to move **vertically** through the stratified tank.

#### Origin

* Initialized as `I = 0` everywhere
* **Never injected directly** by inflow
* Generated only when buoyancy instability exists

Mathematically, `I` emerges from buoyancy work:

$$
I^2 \sim 2g\left|\frac{\Delta \rho}{\rho}\right|\Delta h
$$

#### Role in the Model

* Governs how far a buoyant parcel can penetrate vertically
* Enables inversion resolution and plume propagation
* Decays through entrainment and dilution

---

### 3.4 Mixing Proxy `M` (`[3]`)

#### Definition

`M` is a **velocity-like scalar** representing  **jet- or shear-induced mixing intensity** .

#### Origin

* Initialized only at inflow insertion
* Set proportional to an effective inlet velocity

#### Role in the Model

* Drives horizontal/local mixing
* Modulates entrainment strength in inversion and impulse routines
* Decays as additional mass is mixed

---

## 4. Relationship Between `I` and `M`

* There is no direct conversion `M → I`
* `I` is generated solely by buoyancy and gravity
* `M` influences `I` indirectly by modifying entrainment and dilution

Causal chain:

```
Inflow → M → Mixing / Entrainment → Stratification Change → Buoyancy → I
```

---

## 5. Transport of `I` and `M`

When a buoyant plug moves through the tank during inversion resolution:

* Temperatures and layer heights are rewritten
* `I` and `M` are swapped between layers

This implements  parcel-based transport :

* `I` and `M` are properties of the  *moving fluid* , not fixed locations
* Prevents residual impulse or mixing from remaining behind

---

## 6. Fundamental State Invariants

The following invariants are enforced by the model:

* Total water mass is conserved (except for inflow/outflow)
* Total enthalpy is conserved internally (except for losses)
* No persistent density inversion remains after inversion resolution
* `I ≥ 0`, `M ≥ 0`

---

## 7. Fundamentzustand (Foundation State)

### Structure

```
Fundamentzustand : list[float]
```

* Represents a 1D temperature profile beneath the tank
* Coupled thermally to the bottom water layer

### Role

* Prepares ground for more detailed losses to the ground model
* Influences losses

---

## 8. Kapazitaeten (Wall Capacity Nodes)

### Structure

```
Kapazitaeten : list[float]
```

Each element represents the temperature of a lumped wall control volume.

### Role

* Represents thermal inertia of tank walls
* Couples water temperature to ambient losses

---

## 9. State Persistence

At the end of each timestep:

* `Speicherzustand`
* `Fundamentzustand`
* `Kapazitaeten`

are serialized and written to disk.

On the next call to `main()`, these states are reloaded.

---

*End of state definition documentation.*
