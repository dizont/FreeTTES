# FreeTTES - open source stratified TTES Model

A **1D vertically stratified tank thermal energy storage (TTES) model** for **engineering system simulation**, optimization, and academic research.

This model is designed for **engineers**, not CFD specialists. It captures the dominant physics of stratified hot-water tanks while remaining fast and robust.

---

## Key Features

* Charging and discharging via top/bottom diffusers
* Buoyancy-driven stratification
* Jet-induced mixing near inlets
* Reduced-order plume/inversion solver
* 1D thermal conduction in water and foundation
* Wall heat capacity and ambient losses
* Thermocline detection and usable-energy metrics
* unittest-protected refactoring workflow

---

## What This Model Is (and Is Not)

**This model is:**

* suitable for district heating TTES simulation
* suitable for control, optimization, and MPC studies
* physically interpretable and conservative
* fast enough for long-term simulations

**This model is NOT:**

* a CFD solver
* suitable for diffuser geometry design
* intended to resolve turbulence explicitly

---

## Repository Structure

```
.

src/
    example.py           # How to use
    FreeTTES_model.py    # Core physics
    FreeTTES_io.py       # I/O and persistence
    FreeTTES_config.py   # Configuration
docs/
    model_overview.md
    governing_equations.md
    state_definition.md
    inversion_mixing.md
    numerical_methods.md
```

---

## Core Modeling Concepts

### 1D Layered Representation

The tank is discretized into horizontal layers (cells). Each layer represents a well-mixed water volume with a finite vertical extent.

Each layer is stored as:

```
[T, dh, I, M]
```

| Index   | Symbol | Meaning             | Physical role                  |
| ------- | ------ | ------------------- | ------------------------------ |
| `[0]` | `T`  | Temperature (deg C) | Thermal state                  |
| `[1]` | `dh` | Layer height (m)    | Geometry / volume              |
| `[2]` | `I`  | Impulse proxy (m/s) | Vertical buoyant motion        |
| `[3]` | `M`  | Mixing proxy (m/s)  | initial velocity-driven mixing |

Layers are stored in `Speicherzustand`, a dictionary keyed by height-like coordinates.

---

## How to Run

The model is intended to be called **once per timestep** by an external system simulation.

Typical usage pattern:

```python
# Example usage for the FreeTTES model
import FreeTTES_model as model


def run_one_timestep():
    # Predefined inputs for a single 15-minute timestep
    result = model.main(
        t=0.0,
        dt=900,
        m_VL=100.0,
        m_RL=-100.0,
        T_Zustrom=85.0,
        T_amb=10.0,
        eingabe_volumen=False,
        zustand_uebernehmen=False,
        zustand={}
    )

    print("T_Austritt:", result.get("T_Austritt"))
    print("H_WS:", result.get("H_WS"))
    print("E_nutz:", result.get("E_nutz"))
    print("Q_V_ges:", result.get("Q_V_ges"))


if __name__ == "__main__":
    run_one_timestep()
```

State is persisted automatically between calls.

---

## Documentation

Detailed engineering documentation is provided in `docs/`:

* **Model Overview** - physical assumptions and structure
* **Governing Equations** - formal mathematical formulation
* **State Definition** - data structures and state lifecycle
* **Inversion & Mixing** - core plume-based algorithm
* **Numerical Methods** - stability and discretization

Start with:

```
docs/model_overview.md
```

---

## Intended Users

* Energy system engineers
* Researchers working on TES modeling
* Control and optimization developers

---

## License / Usage

This code is intended for **engineering and research use** and is licensed under the BSD 3-Clause License.

Citation:

---

*For details, see the documentation in `docs/`.*
