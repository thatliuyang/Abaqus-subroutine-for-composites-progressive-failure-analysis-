# Abaqus-subroutine-for-composites-progressive-failure-analysis-
This repository contains an extensively optimized and robust version of the PDALAC (Progressive Damage Analysis of Laminated Composites) UMAT for Abaqus. 

Originally developed by Ammar Khallouf (TU Munich) for 3D elements, this codebase has been  improved by **Brandenburgische Technische Universität Cottbus-Senftenberg (Chair of Polymer-based Lightweight Design)** to address critical bugs, enhance computational efficiency, and prevent numerical instabilities during complex progressive damage simulations.

---

## 🚀 Key Improvements & Benefits

### 1. Performance Optimizations (Speedup)
*   **Two-Pass Fracture Angle Search (Puck IFF)**: 
    *   *Why:* The original implementation swept the fracture plane angle (θ) from -90° to +90° in 1° increments, resulting in 181 expensive trigonometric and root calculations *per integration point, per increment*.
    *   *Fix:* Introduced a two-pass search (a 5° coarse sweep followed by a ±5° fine search).
    *   *Benefit:* **Reduced mathematical evaluations from 181 to 48 (~75% reduction)** while maintaining exact 1° precision. Drastically reduces solving time for large-scale 3D models.
*   **Deferred Jacobian Update (Lazy Loading)**:
    *   *Why:* The damaged stiffness matrix (`DDSDDE`) via `ortho3D` was previously being rebuilt unconditionally, even for perfectly healthy elastic elements.
    *   *Fix:* Moved the secondary `ortho3D` call strictly inside the active damage block (`dmg < 1.0` or `e > 1.0`).
    *   *Benefit:* Skips massive amounts of redundant matrix multiplications in the pre-failure load steps, further accelerating the analysis.

### 2. Theoretical Bug Fixes (Accuracy)
*   **Puck Fiber Fracture (FF) Compressive Sign**:
    *   *Why:* The compressive fiber fracture index was evaluating to a negative number, failing to trigger the `> 1.0` failure condition.
    *   *Fix:* Corrected the sign convention to ensure e(1) remains strictly positive under compression.
*   **Puck Inter-Fiber Fracture (IFF) Self-Consistency**:
    *   *Why:* The slope parameter $P_{nt}$ was calculated as `-1/(2*tan(2*THETAF))`, which is mathematically incorrect according to Puck's theory and failed pure-transverse compression consistency checks.
    *   *Fix:* Corrected to P_nt = -1/tan(2*θ_f) and R_nt = Y_c / (2*tan(θ_f)).
    *   *Benefit:* Guarantees that under pure transverse compression (-Y_c), the failure index evaluates to exactly 1.0 at the specified fracture angle θ_f.

### 3. Numerical Stability & Robustness (No more crashes)
*   **Residual Stiffness Lower Bounds**:
    *   *Why:* Allowing stiffness degradation (`dmg`) to reach $0.0$ inevitably leads to singular stiffness matrices in Abaqus ("Too many attempts" errors).
    *   *Fix:* Implemented a strict lower bound of $0.10$ (10%) for fiber direction and $0.05$ (5%) for matrix/shear directions.
    *   *Benefit:* Allows the simulation to smoothly progress through post-failure phases and redistribute stress without crashing the implicit solver.
*   **Fiber-Driven Extreme Failure Trigger**:
    *   *Why:* The original code triggered element-level extreme failure if *any* index exceeded 2.0 (`maxval(e) > 2.0`). Because matrix strength is low, matrix stress concentrations frequently spiked > 2.0, causing premature element deletion despite healthy fibers.
    *   *Fix:* The extreme failure trigger is now strictly tied to the fiber direction (`e(1) > 2.0`). 
    *   *Benefit:* Prevents premature component failure. Matrix cracking now correctly degrades transverse stiffness while the fibers continue to carry longitudinal load.

### 4. Enhanced Outputs (State Variables - SDVs)
State variables have been restructured for clearer post-processing in Abaqus Viewer:
*   **SDV 7-12**: Failure indices (e1 to e6)
*   **SDV 14-19**: Fracture angles (Separated explicitly into Mode A, Mode B, and Mode C to track exact failure mechanisms).
*   **SDV 20-25**: Failure flags (`fflags`), indicating failure initiation time-steps.

---

## 🧪 Standalone Test Suite
To verify the integrity of the failure criteria outside of Abaqus, we have included `test_puck.f90`.
This standalone Fortran program includes **17 targeted unit tests (33 checks)** covering:
*   Biaxial fiber loading and MGF-coupling
*   IFF Mode A, Mode B, and Mode C isolation
*   Damage degradation rules (Instantaneous vs. Constant Stress)
*   Puck-specific induced shear degradation

**To compile and run the tests (Requires Intel Fortran):**
```cmd
ifx test_puck.f90 -o test_puck.exe
test_puck.exe
```

---

## Acknowledgements
* **Original Codebase**: PDALAC UMAT by Ammar Khallouf, TUM (2019).
https://github.com/ammarkh95/ABAQUS_PDALAC

* **Modifications & Optimizations**: Yang Liu, Brandenburgische Technische Universität Cottbus-Senftenberg (2026).
