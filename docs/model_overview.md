# FreeTTES Model

## Overview

This repository contains a **1D vertically stratified tank thermal energy storage (TTES) model** designed for system-level simulation, optimization, and academic research. The model captures the essential physics of stratified hot water tanks while remaining computationally efficient and robust for long-term simulations.

The code is  **not CFD** . Instead, it is a reduced-order, physics-informed model that resolves:

* charging and discharging via top/bottom ports,
* buoyancy-driven stratification and destratification,
* jet-induced mixing near inlets,
* thermal conduction in water and foundation,
* heat losses via tank walls and ambient coupling,
* thermocline detection and usable-energy metrics.

---

## Repository Structure

```
.
├── speicher_upgraded_deep_v2.py          # Facade / entry point (drop-in replacement)
├── speicher_upgraded_deep_v2_model.py    # Core physics and orchestration
├── speicher_upgraded_deep_v2_io.py       # File I/O and persistence
├── speicher_upgraded_deep_v2_config.py   # Configuration and path handling
├── test_speicher_regression_v2.py        # Regression tests (stillstand → charge → discharge)
└── config.json                           # Example configuration
```

---

## Conceptual Model

### 1D Layered Representation

The tank is discretized into horizontal layers (cells). Each layer represents a well-mixed water volume with a finite vertical extent.

Each layer is stored as:

```
[T, dh, I, M]
```

| Index   | Symbol | Meaning           | Physical role           |
| ------- | ------ | ----------------- | ----------------------- |
| `[0]` | `T`  | Temperature (°C) | Thermal state           |
| `[1]` | `dh` | Layer height (m)  | Geometry / volume       |
| `[2]` | `I`  | Impulse proxy     | Vertical buoyant motion |
| `[3]` | `M`  | Mixing proxy      | Jet/shear-driven mixing |

Layers are stored in `Speicherzustand`, a dictionary keyed by height-like coordinates.

---

## Governing Physical Ideas

### Separation of Mechanisms

The model  **deliberately separates vertical motion and mixing** :

* **Impulse `[2]`**
  * Represents *vertical penetration capability*
  * Is generated **only by buoyancy + gravity**
  * Is never injected directly by inflow
* **Mixing `[3]`**
  * Represents *jet/shear-induced mixing strength*
  * Is generated **only by inflow hydraulics**
  * Drives horizontal/local mixing

There is **no direct conversion** `[3] → [2]`, but **indirect coupling exists** via entrainment and dilution.

---

## Simulation Workflow (`main`)

For each external timestep:

1. Load last state from disk (or initialize at `t = 0`)
2. Determine operating mode:
   * stillstand (no flow)
   * charge
   * discharge
3. Subdivide timestep into stable substeps
4. For each substep:
   * apply inflow and outflow
   * detect density inversions
   * resolve inversions and mixing
   * apply thermal conduction
   * apply wall and ambient heat exchange
5. Compute outputs (energy, mass, outlet temperature, losses, thermocline metrics)
6. Save the new state and finish programm

---

## Core Algorithms

### Inflow (`__Modell_Zustrom`)

* Inserts incoming volume at top or bottom port
* Creates a new layer with:
  * `I = 0`
  * `M = v_effective` (from inlet hydraulics)
* Ensures inflow **does not directly cause vertical motion**

### Outflow (`__Modell_Abstrom`)

* Removes volume from top or bottom region
* Computes outlet temperature via enthalpy averaging

### Inversion Detection (`__Modell_Inversionspruefung`)

* Scans for gravitationally unstable stratification
* Triggers buoyancy correction if needed

### Inversion Resolution (`__Modell_Inversion`)

This is the  **core buoyancy solver** :

* Identifies a moving plug
* Advances it layer-by-layer
* At each interface:
  * entrains part of the encountered layer
  * updates impulse `[2]` via buoyancy and dilution
  * dilutes mixing `[3]`
  * rewrites the two layers

This enforces stable stratification while conserving enthalpy. It will keep on calculating until the inversion is resolved. This can lead to too hasty inversion resolution when compared to reality (inversions are resolved within one given timestep.

### Impulse Propagation (`__Modell_Impuls`)

* Operates only when `I > 0`
* Transports and attenuates existing impulse
* Does **not** generate impulse from zero

### Horizontal Mixing (`__Modell_Horizontalmischung`)

* Activated only when `M > 0`
* Uses a kinetic-energy vs buoyancy-resistance criterion
* Mixes neighboring layers locally
* Dilutes `M` as more mass is mixed

### Thermal Conduction (`__Modell_Waermeleitung`)

* Solves 1D heat diffusion in water + foundation
* Uses an implicit finite-difference scheme
* Solved efficiently with TDMA

### Wall & Ambient Losses (`__Modell_Kapazitaeten`)

* Exchanges heat between water and wall capacity nodes
* Applies wall-to-ambient losses
* Dominates long-term standby behavior

---

## `[2]` and `[3]` — Impulse–Mixing Subsystem

### Origin

* `[3]` (mixing proxy): born at inflow
* `[2]` (impulse proxy): born only from buoyancy instability

### Coupling

* `[3]` affects entrainment coefficients
* Entrainment affects dilution of `[2]`
* Therefore `[3]` influences the *evolution* of `[2]`, but never creates it

### Transport

When a moving plug advances, `[2]` and `[3]` are swapped between layers. This means:

* impulse and mixing are **properties of the moving cell**
* they are not tied to fixed spatial positions

---

## Energy Consistency

Internally conservative processes:

* inversion mixing
* horizontal mixing
* grid refinement
* internal conduction

Processes that intentionally change total energy:

* inflow/outflow enthalpy
* wall and ambient heat losses
* foundation coupling

---

## Outputs and Metrics

* Outlet temperature (from enthalpy)
* Usable mass and energy above a threshold temperature
* Heat loss components
* Thermocline position and thickness*

---

## Validation and Testing

A unit testtest (`test_compare_original_vs_freettes.py`) verifies:

* stillstand → charge → discharge sequence
* equality of scalar outputs
* equality of final internal state

This test should be run after  **any refactor** .

---

## Modeling Assumptions

* 1D vertical stratification
* Local mixing is instantaneous within a timestep
* Turbulence and entrainment are parameterized
* Horizontal homogeinity
* Constant water mass

---

## Intended Use

This model is well suited for:

* district heating simulations
* system optimization and control studies
* academic research on stratification dynamics

It is **not** intended to replace CFD for detailed diffuser design.

---

## Suggested Further Documentation

For larger projects, consider splitting this documentation into:

* `docs/model.md` – physical model
* `docs/state.md` – data structures
* `docs/algorithms.md` – inversion & mixing
* `docs/io.md` – file formats and persistence

---

*This documentation reflects the executable logic of the code and intentionally ignores all inline code comments.*
