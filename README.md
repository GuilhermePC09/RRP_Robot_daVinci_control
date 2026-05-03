# RRP Robot da Vinci — Modeling, Analysis & Simulation

Complete MATLAB implementation of kinematic and dynamic modeling, linearization, stability analysis, and trajectory simulation for a simplified **da Vinci surgical robot arm** modeled as an **RRP (Revolute–Revolute–Prismatic)** manipulator.

---

## Table of Contents

- [Robot Description](#robot-description)
- [Denavit-Hartenberg Parameters](#denavit-hartenberg-parameters)
- [Physical Parameters](#physical-parameters)
- [Project Structure](#project-structure)
- [Module Overview](#module-overview)
- [How to Run](#how-to-run)
- [Implemented Features](#implemented-features)
  - [Kinematics](#kinematics)
  - [Dynamics](#dynamics)
  - [Linearization and State Space](#linearization-and-state-space)
  - [System Analysis](#system-analysis)
  - [Simulations](#simulations)
- [Simulation Scenarios](#simulation-scenarios)
- [Outputs and Figures](#outputs-and-figures)
- [Dependencies](#dependencies)

---

## Robot Description

The robot represents a simplified arm of the da Vinci surgical system, consisting of three joints:

| Joint | Type      | Variable | Description                          |
|-------|-----------|----------|--------------------------------------|
| 1     | Revolute  | `q1`     | Base rotation (yaw around Z axis)    |
| 2     | Revolute  | `q2`     | Arm elevation (pitch)                |
| 3     | Prismatic | `q3`     | Linear extension along the arm axis  |

The end-effector (surgical tool tip) position is fully determined by the three joint variables `[q1, q2, q3]`. The robot operates in a workspace above a patient bed modeled as a horizontal plane at z = -0.2 m.

---

## Denavit-Hartenberg Parameters

The robot is parameterized using the standard D-H convention. The table below defines the symbolic D-H parameters stored in `Robot_dV.m`:

| Link | `d`  | `θ`  | `a`   | `α`     | Joint Type |
|------|------|------|-------|---------|------------|
| 1    | 0    | `q1` | 0     | π/2     | Revolute   |
| 2    | 0    | `q2` | `a2`  | π/2     | Revolute   |
| 3    | `q3` | 0    | 0     | 0       | Prismatic  |

Where `a2 = 0.4 m` is the length of link 2 (the arm segment).

The homogeneous transformation for each joint follows:

```
T_i = | cos(θ)  -sin(θ)cos(α)   sin(θ)sin(α)   a·cos(θ) |
      | sin(θ)   cos(θ)cos(α)  -cos(θ)sin(α)   a·sin(θ) |
      |   0        sin(α)          cos(α)            d    |
      |   0           0               0               1   |
```

---

## Physical Parameters

Numeric values are substituted in `NumSS.m` for the linearized state-space model:

### Link 1 — Cylindrical (Base Rotation)
| Parameter          | Value                |
|--------------------|----------------------|
| Mass (`m1`)        | 4.3175 kg            |
| Radius (`r1`)      | 0.05 m               |
| Height (`h1`)      | 0.2 m                |
| `Ix1 = Iz1`        | `(1/12)·m1·(3r²+h²)` |
| `Iy1`              | `(1/2)·m1·r²`        |
| CoM offset (`lc1`) | 0 (at joint origin)  |

### Link 2 — Cylindrical (Arm Segment)
| Parameter          | Value                |
|--------------------|----------------------|
| Mass (`m2`)        | 8.635 kg             |
| Radius (`r2`)      | 0.05 m               |
| Length (`h2`)      | 0.4 m                |
| `Ix2`              | `(1/2)·m2·r²`        |
| `Iy2 = Iz2`        | `(1/12)·m2·(3r²+h²)` |
| CoM offset (`lc2`) | 0.2 m (center)       |

### Link 3 — Rectangular Block (Prismatic Tool)
| Parameter          | Value                          |
|--------------------|--------------------------------|
| Mass (`m3`)        | 1.32 kg                        |
| Length (`L3`)      | 0.2 m (along x)                |
| Height (`h3`)      | 0.04 m (along y)               |
| Width (`b3`)       | 0.06 m (along z)               |
| `Ix3`              | `(1/12)·m3·(h²+b²)`           |
| `Iy3`              | `(1/12)·m3·(L²+b²)`           |
| `Iz3`              | `(1/12)·m3·(L²+h²)`           |
| CoM offset (`lc3`) | 0.1 m (center)                 |

### Joint Friction and Environment
| Parameter                        | Value       |
|----------------------------------|-------------|
| Viscous damping joint 1 (`bv1`)  | 2 N·m·s/rad |
| Viscous damping joint 2 (`bv2`)  | 2 N·m·s/rad |
| Viscous damping joint 3 (`bv3`)  | 5 N·s/m     |
| Environment stiffness (`K_env3`) | 500 N/m     |
| Gravity (`g`)                    | 9.81 m/s²   |

### Equilibrium Point
| Variable     | Value              |
|--------------|--------------------|
| `q1_bar`     | 0 rad              |
| `q2_bar`     | 0 rad              |
| `q3_bar`     | 0.2 m              |
| `dq_bar`     | [0; 0; 0]          |

---

## Project Structure

```
RRP_Robot_daVinci_control/
└── DaVinci_Robot/
    ├── MainKin.m                    # Main script: kinematics & Simulink library
    ├── MainDyn.m                    # Main script: dynamics, analysis & simulations
    │
    ├── Robot_dV.m                   # Robot definition (D-H table, masses, CoM, inertia)
    ├── DHTransf.m                   # D-H homogeneous transformation matrix
    ├── DKin.m                       # Direct kinematics (T_0_n)
    ├── IKin.m                       # Analytical inverse kinematics
    ├── GeoJac.m                     # Geometric Jacobian
    ├── NumDiff.m                    # Numerical Jacobian (finite differences)
    │
    ├── DDyn_Lagrange.m              # Direct dynamics via Lagrange formulation
    ├── DDyn_NE.m                    # Direct dynamics via Newton-Euler formulation
    ├── ExtractInMat.m               # Extract mass matrix B from symbolic tau
    ├── DecompTau.m                  # Decompose tau into B, phi (Coriolis), G (gravity)
    ├── StateSpaceFunc.m             # Builds nonlinear state-space model f(x,u)
    ├── LinModel.m                   # Linearization via Taylor expansion (A, B, C, D)
    ├── NumSS.m                      # Substitutes numeric values into state-space matrices
    ├── IDyn.m                       # Inverse dynamics (feedforward torque computation)
    │
    ├── GetBmatrix.m                 # Auto-generated: numeric mass matrix B(x)
    ├── Getnvector.m                 # Auto-generated: numeric n(x) = C·dq + G + friction
    ├── GetIKnumeric.m               # Auto-generated: numeric IK function
    │
    ├── RouthTab.m                   # Routh-Hurwitz stability table
    ├── PlotRobot_dV.m               # 3D robot skeleton visualization
    ├── PlotSim.m                    # Simulation results plotter (joint tracking + 3D path)
    ├── AnimSim.m                    # Real-time animation with MP4 video export
    │
    ├── Robot_daVinci_Lib.slx        # Simulink library (FK, IK, Jacobian blocks)
    └── Robot_daVinci_Simul_Kin.slx  # Simulink kinematic simulation model
```

---

## Module Overview

### Core Functions

| Function           | Inputs                              | Outputs                          | Description |
|--------------------|-------------------------------------|----------------------------------|-------------|
| `Robot_dV`         | —                                   | `Robot, M, CoM, I`               | Returns symbolic D-H table, mass vector, center of mass matrix, and inertia cell array |
| `DHTransf`         | `p` (D-H row)                       | `A` (4×4 matrix)                 | Builds a single homogeneous transformation from D-H parameters |
| `DKin`             | `Robot`                             | `T` (4×4 matrix)                 | Computes full forward kinematic transformation T₀ₙ |
| `IKin`             | —                                   | `IK_joints` (2×3 symbolic)       | Analytical inverse kinematics; returns two solutions `[q1, q2, q3]` for each sign of `q3` |
| `GeoJac`           | `Robot`                             | `J` (6×n symbolic)               | Geometric Jacobian using the cross-product method |
| `NumDiff`          | `Robot, Robot_T`                    | `J` (6×n symbolic)               | Numerical Jacobian via finite differences (δ = 1e-6) |
| `DDyn_Lagrange`    | `Robot, M, CoM, I, g_vec`           | `tau, B, phi, G`                 | Full Lagrangian dynamics: mass matrix B, Coriolis vector φ, gravity vector G |
| `DDyn_NE`          | `Robot, M, CoM, I, g_vec`           | `tau`                            | Newton-Euler dynamics via forward/backward recursion |
| `StateSpaceFunc`   | `B, n`                              | `f, x, u`                        | Builds nonlinear state-space: `ẋ = f(x,u)` where `x = [q; dq]` |
| `LinModel`         | `f_ss, x_ss, u_ss, G`              | `A, B, C, D`                     | Linearizes f(x,u) via Jacobian around equilibrium `[0; 0; q3_bar; 0; 0; 0]` |
| `NumSS`            | `A_lin, B_lin, C_lin, D_lin`        | `A_num, B_num, C_num, D_num, ...`| Substitutes all symbolic parameters with numeric values |
| `IDyn`             | `num_points, vec_ref, dt`           | `q_ref, dq_ref, ddq_ref, tau_ff` | Computes inverse kinematics + inverse dynamics feedforward torques along a trajectory |
| `RouthTab`         | `coeffs`                            | `rh_table`                       | Builds the Routh array from a characteristic polynomial |
| `PlotRobot_dV`     | `q, a2`                             | `h` (plot handles)               | Draws 3D robot skeleton with joints, links, and patient bed |
| `PlotSim`          | `t_traj, vec_ref, q_ref, tau, ...`  | —                                | Plots joint tracking, 3D Cartesian trajectory, and feedforward torques |
| `AnimSim`          | `t_non, x_non, vec_ref`             | —                                | Animates the simulation and exports `Robot_Animation.mp4` |

### Auto-Generated Functions

These files are generated at runtime by `matlabFunction` inside `MainDyn.m` and `MainKin.m`:

| File               | Generated By       | Purpose |
|--------------------|--------------------|---------|
| `GetBmatrix.m`     | `MainDyn.m`        | Evaluates the numeric mass matrix B(x) given the full state vector `x = [q; dq]` |
| `Getnvector.m`     | `MainDyn.m`        | Evaluates the numeric disturbance vector `n(x) = C(q,dq)·dq + G(q) + Bv·dq + K_env·q` |
| `GetIKnumeric.m`   | `MainDyn.m`        | Evaluates the numeric inverse kinematics solution for given end-effector position `(Pe_x, Pe_y, Pe_z)` |

---

## How to Run

### Prerequisites

Add the `DaVinci_Robot/` folder to the MATLAB path:

```matlab
addpath('DaVinci_Robot');
```

### Kinematics Analysis

Run `MainKin.m` to:
- Compute forward and inverse kinematics symbolically
- Generate the Geometric Jacobian
- Validate it numerically via finite differences
- Analyze kinematic singularities (`det(J) = 0`)
- Visualize the robot in 3D
- Generate and populate the Simulink library `Robot_daVinci_Lib.slx`

```matlab
run('DaVinci_Robot/MainKin.m')
```

> **Note:** `MainKin.m` requires Simulink to be open for the `matlabFunctionBlock` calls. Open `Robot_daVinci_Lib.slx` manually before running if needed.

### Dynamics Analysis and Simulation

Run `MainDyn.m` to execute the complete dynamics pipeline:

```matlab
run('DaVinci_Robot/MainDyn.m')
```

This script runs in sequential sections (use "Run Section" in the MATLAB editor for step-by-step execution):

| Section | Description |
|---------|-------------|
| Direct Dynamics (Lagrange) | Computes `B`, `φ`, `G` symbolically |
| Nonlinear State Space | Adds friction and spring; builds `f(x,u)` |
| Linearization | Computes `A, B, C, D` matrices symbolically |
| Numeric Substitution | Produces fully numeric 6×6 state-space model |
| Reduced System | Eliminates the passive `q1` state to get a 5×5 model |
| Model Analysis | Transfer functions, poles, zeros, Routh table |
| Frequency Response | Bode plots and stability margins |
| Step Response | Linear vs. nonlinear step simulations |
| Trajectory Scenarios | Surgical task simulations (circles, helix) |

> **Warning:** The symbolic computation sections (Lagrange, linearization) are computationally intensive and may take several minutes on the first run.

---

## Implemented Features

### Kinematics

- **Forward Kinematics:** Full symbolic transformation matrix `T₀ₙ(q1, q2, q3)` via chain of D-H matrices
- **Inverse Kinematics (Analytical):**
  - `q1 = atan2(Pe_y, Pe_x)` — base rotation
  - `q3 = ±sqrt(Pe_x² + Pe_y² + Pe_z² − a2²)` — reach (two solutions)
  - `q2 = atan2(Pe_z, r_xy) ± atan2(q3, a2)` — elevation angle
- **Geometric Jacobian:** 6×3 symbolic matrix distinguishing revolute (`Jv = z × (oₙ − oᵢ)`, `Jw = z`) and prismatic (`Jv = z`, `Jw = 0`) joints
- **Numerical Jacobian:** Finite-difference validation with `δ = 1e-6`
- **Singularity Analysis:** Symbolic determinant of the position Jacobian `det(Jp) = 0`

### Dynamics

**Lagrangian Method (`DDyn_Lagrange`):**
1. Computes CoM Jacobians `Jv_ci` and `Jw_i` for each link
2. Assembles total kinetic energy `K = Σ ½mᵢvᶜᵢᵀvᶜᵢ + ½ωᵢᵀIᵢωᵢ`
3. Assembles total potential energy `U = Σ −mᵢgᵀpᶜᵢ`
4. Extracts `B` (Hessian of K w.r.t. `dq`), `G` (gradient of U w.r.t. `q`), and `φ` (Christoffel symbols)

**Newton-Euler Method (`DDyn_NE`):**
1. Forward recursion: propagates angular velocities `ω`, angular accelerations `ω̇`, and linear accelerations `v̇`
2. Backward recursion: propagates forces `f` and moments `μ` from tip to base
3. Extracts joint torques/forces from `τᵢ = μᵢᵀRᵢᵀz₀` (revolute) or `τᵢ = fᵢᵀRᵢᵀz₀` (prismatic)

**Disturbance Vector `n(q, dq)` includes:**
- Coriolis and centrifugal terms `φ(q, dq)`
- Gravity terms `G(q)`
- Viscous friction `Bv · dq`
- Environmental spring force `K_env3 · (q3 − q3_bar)` on joint 3

### Linearization and State Space

The **nonlinear model** is:

```
ẋ = f(x, u)    where  x = [q1; q2; q3; dq1; dq2; dq3],  u = [τ1; τ2; F3]
```

**Linearization** around the horizontal equilibrium `x̄ = [0; 0; q3_bar; 0; 0; 0]`:

```
δẋ = A·δx + B·δu
 y  = C·δx + D·δu
```

- `A = ∂f/∂x |_{x̄,ū}` — 6×6 state matrix
- `B = ∂f/∂u |_{x̄,ū}` — 6×3 input matrix
- `C = I₃|₀` — 3×6 output matrix (positions only, `y = q`)
- `D = 0` — 3×3 feedthrough matrix

Since `q1` is passive at equilibrium (decoupled), a **5×5 reduced model** is used for analysis, keeping states `[q2, q3, dq1, dq2, dq3]` and outputs `[q2, q3]`.

### System Analysis

- **Transfer Functions:** `G22(s) = q2(s)/τ2(s)` and `G33(s) = q3(s)/F3(s)` extracted from the reduced state-space model
- **Poles and Zeros:** Open-loop poles (`pole`) and MIMO transmission zeros (`tzero`)
- **Characteristic Polynomial:** `det(sI − A) = 0` computed via `poly(A_reduced)`
- **Routh-Hurwitz Table:** Constructed column by column using the ε-method for zero pivots
- **Frequency Domain:**
  - Bode plots for each channel
  - Bandwidth (`bandwidth`)
  - Phase margin and gain margin (`margin`)
- **DC Gains:** Steady-state displacement per unit step input (`dcgain`)

### Simulations

All trajectory simulations compare the **linearized model** (`lsim`) against the **nonlinear model** (`ode45` with `RelTol = 1e-6`, `AbsTol = 1e-9`).

**ODE integration:**
```matlab
dxdt = [dq; B(x)\(u − n(x))]
```

**Inverse dynamics feedforward** computes the required torques along a desired trajectory:
```matlab
τ_ff = B(q_ref)·q̈_ref + n(q_ref, dq_ref)
```

---

## Simulation Scenarios

### Scenario 1A — Step on Joint 2 (0.1 N·m)
Step torque applied to the revolute arm joint. Compares linear and nonlinear angular displacement `q2(t)` over 5 seconds.

### Scenario 1B — Step on Joint 3 (−100 N)
Step force applied to the prismatic joint. Compares linear and nonlinear linear displacement `q3(t)` over 5 seconds.

### Scenario 2A — Sinusoidal on Joint 2: `sin(5t)` N·m
Sinusoidal torque at 5 rad/s. Tests frequency response and linear/nonlinear agreement over 20 seconds.

### Scenario 2B — Sinusoidal on Joint 3: `50·sin(2t)` N
Sinusoidal force at 2 rad/s. Larger amplitude to stress-test the linearization.

### Scenario 3A — YZ Plane Circle (Suturing Task)
- Trajectory: circle of radius **4 cm** in the YZ plane
- Center: `(x, y, z) = (0.4, 0, −0.2)` m
- Duration: 4 seconds (1 full circle)
- Inverse kinematics at each timestep → feedforward torques via inverse dynamics
- Plots: joint tracking, 3D Cartesian trajectory, computed torques, animation

### Scenario 3B — XY Plane Circle
- Trajectory: circle of radius **4 cm** in the XY plane (base rotation + reach)
- Center: same as 3A
- Duration: 4 seconds

### Scenario 3C — Helical Trajectory
- Trajectory: helix advancing **10 cm** along X while rotating in YZ plane
- Radius: **3 cm**, frequency: **0.5 Hz**
- Duration: 10 seconds (5 full turns)
- Tests full 3D motion with combined rotation and reach variation

---

## Outputs and Figures

Running `MainDyn.m` produces the following figures and console outputs:

| Output | Description |
|--------|-------------|
| Console: equilibrium torque | Torques `[τ1; τ2; F3]` required to hold `q3 = 0.2 m` |
| Console: numeric A, B matrices | 6×6 linearized state-space matrices |
| Console: eigenvalues | Open-loop poles of the full system |
| Console: reduced A, B matrices | 5×5 reduced-system matrices |
| Console: transfer functions G22, G33 | Simplified SISO transfer functions |
| Console: poles and MIMO zeros | Pole/zero locations |
| Console: Routh table | Routh-Hurwitz array for stability check |
| Console: frequency metrics | Bandwidth, phase margin, gain margin for each channel |
| Console: DC gains | Steady-state sensitivity (rad/N·m, m/N) |
| Figure: Bode plots | Frequency response for G_arm and G_prism |
| Figure: Step response (joint 2) | `q2(t)` for 1 N·m step |
| Figure: Step response (joint 3) | `q3(t)` for 1 N step |
| Figure: Initial condition response | `q3(t)` decaying from 5 cm disturbance |
| Figures: Scenarios 1A, 1B, 2A, 2B | Linear vs. nonlinear comparison plots |
| Figures: Scenarios 3A, 3B, 3C | Joint tracking + 3D path + torques |
| `Robot_Animation.mp4` | Exported video of the robot executing the last trajectory |

---

## Dependencies

| Toolbox | Required For |
|---------|-------------|
| **MATLAB Symbolic Math Toolbox** | All symbolic computations (`syms`, `simplify`, `jacobian`, `hessian`, `subs`, `matlabFunction`) |
| **MATLAB Control System Toolbox** | `ss`, `tf`, `pole`, `tzero`, `bode`, `margin`, `bandwidth`, `lsim`, `initial`, `dcgain`, `minreal` |
| **Simulink** (optional) | `MainKin.m` only — for generating function blocks in `Robot_daVinci_Lib.slx` |

Tested on **MATLAB R2023b** and later. The Symbolic Math and Control System Toolboxes are mandatory for running `MainDyn.m`.
