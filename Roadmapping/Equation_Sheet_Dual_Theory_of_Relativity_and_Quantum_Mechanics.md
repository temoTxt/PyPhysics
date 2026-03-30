# Equation Sheet: Dual Theory of Relativity and Quantum Mechanics

**Compiled from:**
1. Gill & Zachary, "Two Mathematically Equivalent Versions of Maxwell's Equations" (2011)
2. Gill, Zachary & Lindesay, "The Classical Electron Problem" (2001)
3. Gill & Ares de Parga, "FoundationsII-Classical" (2019)
4. Gill & Ares de Parga, "FOUNDATIONS FOR QED I MATHEMATICAL" (2019)
5. Gill, Ares de Parga, Morris & Wade, "Dual Relativistic Quantum Mechanics I" (2021)

---

## Table of Contents

1. [Fundamental Definitions and Identities](#1-fundamental-definitions-and-identities)
2. [Proper-Time Transformations](#2-proper-time-transformations)
3. [Maxwell's Equations - Dual Formulation](#3-maxwells-equations---dual-formulation)
4. [Wave Equations and Radiation Reaction](#4-wave-equations-and-radiation-reaction)
5. [Many-Particle Theory](#5-many-particle-theory)
6. [Canonical Proper-Time Hamiltonian](#6-canonical-proper-time-hamiltonian)
7. [Quantum Theory - Hilbert Spaces](#7-quantum-theory---hilbert-spaces)
8. [Feynman Operator Calculus](#8-feynman-operator-calculus)
9. [Dual Relativistic Wave Equations](#9-dual-relativistic-wave-equations)
10. [Dual Dirac Equation and g-Factor](#10-dual-dirac-equation-and-g-factor)
11. [Isodual Numbers and Antimatter](#11-isodual-numbers-and-antimatter)

---

## 1. Fundamental Definitions and Identities

### 1.1 Proper-Time Definitions

**First Representation (Observer Proper-Time):**
$$d\tau = \gamma^{-1}(\mathbf{w}) dt \tag{1.1}$$
$$d\tau = \gamma^{-1}(\mathbf{w}') dt' \tag{1.2}$$
where $\gamma^{-1}(\mathbf{w}) = \sqrt{1 - \mathbf{w}^2/c^2}$ and $\mathbf{w} = d\mathbf{x}/dt$

**Second Representation (Source Proper-Time):**
$$d\tau^2 = dt^2 - \frac{1}{c^2} d\mathbf{x}^2 \tag{1.3}$$
$$dt^2 = d\tau^2 + \frac{1}{c^2} d\mathbf{x}^2 \tag{1.4}$$

**Collaborative Speed of Light:**
$$b = \sqrt{c^2 + \mathbf{u}^2} \tag{1.5}$$
where $\mathbf{u} = d\mathbf{x}/d\tau$ is the proper velocity

**Key Identity:**
$$cdt = bd\tau \tag{1.6}$$

**Velocity Relationship:**
$$\frac{\mathbf{w}}{c} = \frac{\mathbf{u}}{b} \tag{1.7}$$
$$\mathbf{u} = \mathbf{w} \left[1 - \frac{\mathbf{w}^2}{c^2}\right]^{-1/2} \tag{1.8}$$
$$\mathbf{w} = \mathbf{u} \left[1 + \frac{\mathbf{u}^2}{c^2}\right]^{-1/2} \tag{1.9}$$

**Momentum Representation:**
$$\mathbf{p} = m\gamma(\mathbf{w})\mathbf{w} = m\mathbf{u} \tag{1.10}$$

### 1.2 Energy-Momentum Relations

**Total Energy:**
$$H = \gamma(\mathbf{w}) mc^2 = \sqrt{c^2\mathbf{p}^2 + m^2c^4} \tag{1.11}$$

**Proper-Time Relation:**
$$d\tau = \frac{mc^2}{H} dt \tag{1.12}$$

**Collaborative Speed Identity:**
$$H = mcb \tag{1.13}$$

---

## 2. Proper-Time Transformations

### 2.1 Transformation of Derivatives

**Time Derivative:**
$$\frac{1}{c}\frac{\partial}{\partial t} = \frac{1}{b}\frac{\partial}{\partial \tau} \tag{2.1}$$
$$\frac{1}{c}\frac{\partial}{\partial t'} = \frac{1}{b'}\frac{\partial}{\partial \tau} \tag{2.2}$$

**Spatial Gradient:**
$$\nabla = \gamma(\mathbf{v}) \left[ \nabla' - \frac{\mathbf{v}}{cb'} \frac{\partial}{\partial \tau} \right] \tag{2.3}$$
$$\nabla' = \gamma(\mathbf{v}) \left[ \nabla + \frac{\mathbf{v}}{cb} \frac{\partial}{\partial \tau} \right] \tag{2.4}$$

### 2.2 Proper-Time Group Transformations

**Scale Factor Transformation:**
$$b' = \gamma(\mathbf{v}) \left[ b - \frac{\mathbf{u} \cdot \mathbf{v}}{c} \right] \tag{2.5}$$
$$b = \gamma(\mathbf{v}) \left[ b' + \frac{\mathbf{u}' \cdot \mathbf{v}}{c} \right] \tag{2.6}$$

**Position Transformation:**
$$\mathbf{x}' = \gamma(\mathbf{v}) \left[ \mathbf{x}^{\dagger} - \frac{\mathbf{v}}{c} \bar{b}_{\tau} \tau \right] \tag{2.7}$$
$$\mathbf{x} = \gamma(\mathbf{v}) \left[ \mathbf{x}'^{\dagger} + \frac{\mathbf{v}}{c} \bar{b}'_{\tau} \tau \right] \tag{2.8}$$

where $\mathbf{x}^{\dagger} = \mathbf{x}/\gamma(\mathbf{v}) - (1 - \gamma(\mathbf{v})) \left[ \mathbf{v} \cdot \mathbf{x}/(\gamma(\mathbf{v})\mathbf{v}^2) \right] \mathbf{v}$

**Proper Velocity Transformation:**
$$\mathbf{u}' = \gamma(\mathbf{v}) \left[ \mathbf{u}^{\dagger} - \frac{\mathbf{v}}{c} b \right] \tag{2.9}$$
$$\mathbf{u} = \gamma(\mathbf{v}) \left[ \mathbf{u}'^{\dagger} + \frac{\mathbf{v}}{c} b' \right] \tag{2.10}$$

**Acceleration Transformation:**
$$\mathbf{a}' = \gamma(\mathbf{v}) \left\{ \mathbf{a}^{\dagger} - \mathbf{v} \left[ \frac{\mathbf{u} \cdot \mathbf{a}}{bc} \right] \right\} \tag{2.11}$$
$$\mathbf{a} = \gamma(\mathbf{v}) \left\{ \mathbf{a}'^{\dagger} + \mathbf{v} \left[ \frac{\mathbf{u}' \cdot \mathbf{a}'}{b'c} \right] \right\} \tag{2.12}$$

### 2.3 Charge and Current Density Transformations

**Current Density:**
$$\mathbf{J}' = \mathbf{J} + (\gamma - 1) \frac{(\mathbf{J} \cdot \mathbf{v})}{\mathbf{v}^2} \mathbf{v} - \gamma \frac{b}{c} \rho \mathbf{v} \tag{2.13}$$

**Charge Density:**
$$b'\rho' = \gamma(\mathbf{v}) \left[ b\rho - \frac{\mathbf{J} \cdot \mathbf{v}}{c} \right] \tag{2.14}$$

**Alternative Charge Density Form:**
$$\rho' = \frac{\rho - (\mathbf{J} \cdot \mathbf{v}/bc)}{1 - (\mathbf{u} \cdot \mathbf{v}/bc)} \tag{2.15}$$

**For Current $\mathbf{J} = \rho \mathbf{w}$:**
$$\rho' = \rho \gamma(\mathbf{v}) \left[ 1 - \frac{\mathbf{w} \cdot \mathbf{v}}{c^2} \right] \tag{2.16}$$

**For Current $\mathbf{J} = \rho \mathbf{u}/b$:**
$$\rho' = \rho \frac{1 - (\mathbf{u} \cdot \mathbf{v}/b^2)}{1 - (\mathbf{u} \cdot \mathbf{v}/bc)} \tag{2.17}$$

---

## 3. Maxwell's Equations - Dual Formulation

### 3.1 Conventional Maxwell's Equations (Observer Time)

$$\nabla \cdot \mathbf{B} = 0 \tag{3.1}$$
$$\nabla \cdot \mathbf{E} = 4\pi\rho \tag{3.2}$$
$$\nabla \times \mathbf{E} = -\frac{1}{c} \frac{\partial \mathbf{B}}{\partial t} \tag{3.3}$$
$$\nabla \times \mathbf{B} = \frac{1}{c} \left[ \frac{\partial \mathbf{E}}{\partial t} + 4\pi\rho \mathbf{w} \right] \tag{3.4}$$

### 3.2 Dual Maxwell's Equations (Source Proper-Time)

$$\nabla \cdot \mathbf{B} = 0 \tag{3.5}$$
$$\nabla \cdot \mathbf{E} = 4\pi\rho \tag{3.6}$$
$$\nabla \times \mathbf{E} = -\frac{1}{b} \frac{\partial \mathbf{B}}{\partial \tau} \tag{3.7}$$
$$\nabla \times \mathbf{B} = \frac{1}{b} \left[ \frac{\partial \mathbf{E}}{\partial \tau} + 4\pi\rho \mathbf{u} \right] \tag{3.8}$$

### 3.3 Four-Vector Formulation

**Electromagnetic Tensor:**
$$F = \begin{bmatrix} 0 & B_z & -B_y & -iE_x \\ -B_z & 0 & B_x & -iE_y \\ B_y & -B_x & 0 & -iE_z \\ iE_x & iE_y & iE_z & 0 \end{bmatrix} \tag{3.9}$$

**Four-Derivative:**
$$\frac{\partial}{\partial x_4} = -\frac{i}{b}\frac{\partial}{\partial \tau} \tag{3.10}$$

**Sourceless Equations:**
$$\frac{\partial F_{\alpha\beta}}{\partial x_{\gamma}} + \frac{\partial F_{\beta\gamma}}{\partial x_{\alpha}} + \frac{\partial F_{\gamma\alpha}}{\partial x_{\beta}} = 0 \tag{3.11}$$

**Equations with Sources:**
$$\frac{\partial F_{\alpha\beta}}{\partial x_{\beta}} = \frac{4\pi}{b} J_{\alpha}, \quad J_{\alpha} = (J_x, J_y, J_z, ib\rho) \tag{3.12}$$

**Lorentz Transformation Matrix:**
$$[a_{\mu\nu}] = \begin{bmatrix} 1 + (\gamma - 1)(v_x^2/v^2) & (\gamma - 1)[(v_x v_y)/v^2] & (\gamma - 1)[(v_x v_z)/v^2] & i\gamma \frac{v_x}{c} \\ (\gamma - 1)[(v_x v_y)/v^2] & 1 + (\gamma - 1)(v_y^2/v^2) & (\gamma - 1)[(v_y v_z)/v^2] & i\gamma \frac{v_y}{c} \\ (\gamma - 1)[(v_x v_z)/v^2] & (\gamma - 1)[(v_y v_z)/v^2] & 1 + (\gamma - 1)(v_z^2/v^2) & i\gamma \frac{v_z}{c} \\ -i\gamma \frac{v_x}{c} & -i\gamma \frac{v_y}{c} & -i\gamma \frac{v_z}{c} & \gamma \end{bmatrix} \tag{3.13}$$

**Field Transformations:**
$$\mathbf{E}' = \gamma \left[ \mathbf{E} + \frac{1}{c} \left( \mathbf{v} \times \mathbf{B} \right) \right] - (\gamma - 1) \frac{(\mathbf{E} \cdot \mathbf{v})}{\mathbf{v}^2} \mathbf{v} \tag{3.14}$$
$$\mathbf{B}' = \gamma \left[ \mathbf{B} - \frac{1}{c} \left( \mathbf{v} \times \mathbf{E} \right) \right] - (\gamma - 1) \frac{(\mathbf{B} \cdot \mathbf{v})}{\mathbf{v}^2} \mathbf{v} \tag{3.15}$$

---

## 4. Wave Equations and Radiation Reaction

### 4.1 Wave Equations for Potentials

**Vector Potential Wave Equation:**
$$\nabla \left[ \nabla \cdot \mathbf{A} + \frac{1}{b} \frac{\partial \Phi}{\partial \tau} \right] + \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{1}{b} \frac{\partial \mathbf{A}}{\partial \tau} \right] - \nabla^2 \mathbf{A} = \frac{1}{b} \left( 4\pi \rho \mathbf{u} \right) \tag{4.1}$$

**Scalar Potential Wave Equation:**
$$-\nabla^2 \Phi - \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \nabla \cdot \mathbf{A} \right] = 4\pi \rho \tag{4.2}$$

### 4.2 Lorentz Gauge Condition

$$\nabla \cdot \mathbf{A} + \frac{1}{b} \frac{\partial \Phi}{\partial \tau} = 0 \tag{4.3}$$

### 4.3 Wave Equations (Gauge-Independent Form)

**Electric Field Wave Equation:**
$$\frac{1}{b^2} \frac{\partial^2 \mathbf{E}}{\partial \tau^2} - \frac{\mathbf{u} \cdot \mathbf{a}}{b^4} \left[ \frac{\partial \mathbf{E}}{\partial \tau} \right] - \nabla^2 \cdot \mathbf{E} = -\nabla (4\pi\rho) - \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{4\pi(\rho\mathbf{u})}{b} \right] \tag{4.4}$$

**Magnetic Field Wave Equation:**
$$\frac{1}{b^2} \frac{\partial^2 \mathbf{B}}{\partial \tau^2} - \frac{\mathbf{u} \cdot \mathbf{a}}{b^4} \left[ \frac{\partial \mathbf{B}}{\partial \tau} \right] - \nabla^2 \cdot \mathbf{B} = \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{4\pi \nabla \times \mathbf{J}}{b} \right] \tag{4.5}$$

### 4.4 Klein-Gordon Form

**Change of Variables:**
$$\Phi = \left(\frac{b}{c}\right)^{1/2} g \tag{4.6}$$

**Klein-Gordon Equation:**
$$\frac{1}{b^2} \frac{\partial^2 g}{\partial \tau^2} - \nabla^2 g + \left[ \frac{\ddot{b}}{2b^3} - \frac{5\dot{b}^2}{4b^4} \right] g = 4\pi\rho \left( \frac{c}{b} \right)^{1/2} \tag{4.7}$$

**Effective Mass:**
$$\mu = \left\{ \frac{\hbar^2}{b^2} \left[ \frac{\ddot{b}}{2b^3} - \frac{5\dot{b}^2}{4b^4} \right] \right\}^{1/2} \tag{4.8}$$

### 4.5 Radiation Potentials

**Lienard-Wiechert Potentials (Observer Time):**
$$\mathbf{A} = \frac{q\mathbf{w}}{cs}, \quad \Phi = \frac{q}{s}, \quad s = r - \left(\frac{\mathbf{r} \cdot \mathbf{w}}{c}\right) \tag{4.9}$$

**Dual Potentials (Source Proper-Time):**
$$\mathbf{A} = \frac{q\mathbf{u}}{bs}, \quad \Phi = \frac{q}{s}, \quad s = r - \left(\frac{\mathbf{r} \cdot \mathbf{u}}{b}\right) \tag{4.10}$$

**Retarded Time Relation:**
$$c(t - t') = \int_{\tau'}^{\tau} b(s) ds \tag{4.11}$$

**Field Derivatives:**
$$\frac{1}{\bar{b}}\frac{\partial}{\partial \tau} = \frac{1}{b} \cdot \frac{r}{s}\frac{\partial}{\partial \tau'} \tag{4.12}$$
$$\nabla = \nabla_1 - \frac{\mathbf{r}}{bs} \cdot \frac{\partial}{\partial \tau'} \tag{4.13}$$

### 4.6 Electric and Magnetic Fields

**Electric Field:**
$$\mathbf{E}(\mathbf{x},\tau) = \frac{q\mathbf{r}_u(1-\mathbf{u}^2/b^2)}{s^3} + \frac{q[\mathbf{r} \times (\mathbf{r}_u \times \mathbf{a})]}{b^2s^3} + \frac{q(\mathbf{u} \cdot \mathbf{a})[\mathbf{r} \times (\mathbf{u} \times \mathbf{r})]}{b^4s^3} \tag{4.14}$$

**Magnetic Field:**
$$\mathbf{B}(\mathbf{x},\tau) = \frac{q(\mathbf{r}_u \times \mathbf{r})(1 - \mathbf{u}^2/b^2)}{rs^3} + \frac{q\mathbf{r} \times [\mathbf{r} \times (\mathbf{r}_u \times \mathbf{a})]}{rb^2s^3} + \frac{qr(\mathbf{u} \cdot \mathbf{a})(\mathbf{r} \times \mathbf{u})}{b^4s^3} \tag{4.15}$$

where $\mathbf{r}_u = \mathbf{r} - \frac{r}{b}\mathbf{u}$

---

## 5. Many-Particle Theory

### 5.1 Global Variables

**Effective Mass Energy:**
$$Mc^2 = \sqrt{H^2 - c^2 \mathbf{P}^2} \tag{5.1}$$

**Total Momentum:**
$$\mathbf{P} = \sum_{i=1}^n \mathbf{p}_i \tag{5.2}$$

**Hamiltonian Representation:**
$$H = \sqrt{c^2 \mathbf{P}^2 + M^2 c^4} \tag{5.3}$$

**Canonical Center of Mass:**
$$\mathbf{X} = \frac{1}{H} \sum_{i=1}^{n} H_i \mathbf{x}_i + \frac{c^2 (\mathbf{S} \times \mathbf{P})}{H (Mc^2 + H)} \tag{5.4}$$

### 5.2 Global Proper-Time

**Global Proper-Time Definition:**
$$d\tau = \frac{Mc^2}{H} dt \tag{5.5}$$

**Global Velocity:**
$$\mathbf{U} = \frac{d\mathbf{X}}{d\tau} = \frac{\mathbf{P}}{M} = \frac{1}{M} \sum_{i=1}^{n} m_i \mathbf{u}_i \tag{5.6}$$

**Global Collaborative Speed:**
$$b = \sqrt{\mathbf{U}^2 + c^2} \tag{5.7}$$

**Global Time Relation:**
$$cdt = bd\tau \tag{5.8}$$

### 5.3 Individual Particle Relations

**Individual Proper-Time:**
$$d\tau_i = \frac{H m_i}{H_i M} d\tau \tag{5.9}$$

**Individual Collaborative Speed:**
$$b_i = \frac{cb}{\sqrt{b^2 - \mathbf{v}_i^2}} \tag{5.10}$$

where $\mathbf{v}_i = d\mathbf{x}_i/d\tau$

**Velocity Relations:**
$$\mathbf{u}_i = \frac{c\mathbf{v}_i}{\sqrt{b^2 - \mathbf{v}_i^2}} \tag{5.11}$$
$$\frac{\mathbf{u}_i}{b_i} = \frac{\mathbf{v}_i}{b} \tag{5.12}$$

### 5.4 Center of Mass Position

**Global Position:**
$$\mathbf{X} = \frac{1}{H} \sum_{i=1}^{n} H_i \mathbf{x}_i + \mathbf{Y} \tag{5.13}$$

**Dual Form:**
$$\mathbf{X} = \frac{1}{H} \sum_{i=1}^{n} H_i \mathbf{x}_i + \frac{c^2 (\mathbf{S} \times \mathbf{P})}{H(Mc^2 + H)} \tag{5.14}$$

### 5.5 Clock Relationships

**Observable Evolution (Observer Time):**
$$\frac{dW}{dt} = \sum_{i=1}^{n} \frac{d\tau_i}{dt} \{K_i, W\} \tag{5.15}$$

**Observable Evolution (Proper-Time):**
$$\frac{dW}{d\tau} = \sum_{i=1}^{n} \frac{d\tau_i}{d\tau} \{K_i, W\} \tag{5.16}$$

**General Form:**
$$dW = \sum_{i=1}^{n} d\tau_i \{ K_i, W \} \tag{5.17}$$

### 5.6 Lagrangian Formulation

**Global Lagrangian:**
$$\mathcal{L}dt = \sum_{i=1}^{n} \mathbf{p}_{i} \cdot d\mathbf{x}_{i} - Hdt = \sum_{i=1}^{n} \mathbf{p}_{i} \cdot d\mathbf{x}_{i} - Kd\tau + dS \tag{5.18}$$

**Individual Lagrangians:**
$$\mathcal{L}dt = \sum_{i=1}^{n} \mathcal{L}_{i}d\tau_{i} + dS \tag{5.19}$$

---

## 6. Canonical Proper-Time Hamiltonian

### 6.1 Hamiltonian Definition

**Proper-Time Hamiltonian:**
$$\{K, W\} = \frac{H}{mc^2} \{H, W\}, \quad K|_{\mathbf{p}=0} = H|_{\mathbf{p}=0} = mc^2 \tag{6.1}$$

**Poisson Bracket Form:**
$$\{K,W\} = \left[\frac{H}{mc^2}\frac{\partial H}{\partial \mathbf{p}}\right]\frac{\partial W}{\partial \mathbf{x}} - \left[\frac{H}{mc^2}\frac{\partial H}{\partial \mathbf{x}}\right]\frac{\partial W}{\partial \mathbf{p}} \tag{6.2}$$

**Solution:**
$$K = \frac{H^2}{2mc^2} + \frac{mc^2}{2} \tag{6.3}$$

**General Solution:**
$$K = mc^{2} + \int_{mc^{2}}^{H} \frac{\bar{H}}{mc^{2}} d\bar{H} \tag{6.4}$$

### 6.2 Specific Forms

**Form 1 (Fixed Lorentz Frame):**
$$K = \frac{H^2}{mc^2} = \frac{\mathbf{p}^2}{mc^2} + mc^2 \tag{6.5}$$

**Form 2 (Variable Lorentz Frame):**
$$K = \frac{H^2}{2mc^2} + \frac{mc^2}{2} = \frac{\mathbf{p}^2}{2m} + mc^2 \tag{6.6}$$

**Form 3 (Fixed Momentum):**
$$K = Mc^2 = \sqrt{H^2 - c^2 \mathbf{P}^2} \tag{6.7}$$

### 6.3 Hamilton's Equations

**Velocity:**
$$\frac{d\mathbf{x}}{d\tau} = \frac{\partial K}{\partial \mathbf{p}} = \frac{H}{mc^2} \left( \frac{c^2 \boldsymbol{\pi}}{H_0} \right) = \frac{b}{c} \left( \frac{c^2 \boldsymbol{\pi}}{H_0} \right) \tag{6.8}$$

**Momentum Evolution:**
$$\frac{d\mathbf{p}}{d\tau} = \frac{b}{c} \frac{\left[ \left( c^{2}\boldsymbol{\pi} \cdot \boldsymbol{\nabla} \right) \mathbf{A} + \frac{e}{b} \left( c^{2}\boldsymbol{\pi} \times \mathbf{B} \right) \right]}{H_{0}} - \frac{b}{c} \boldsymbol{\nabla} V \tag{6.9}$$

**Simplified Form:**
$$\frac{d\mathbf{p}}{d\tau} = \frac{e}{c} \left(\mathbf{u} \cdot \nabla\right) \mathbf{A} + \frac{e}{c} \mathbf{u} \times \mathbf{B} - \nabla V \frac{b}{c} \left[1 + \frac{V}{H_0}\right] \tag{6.10}$$

### 6.4 Lorentz Force (Dual Form)

$$\frac{c}{b} \frac{d\boldsymbol{\pi}}{d\tau} = \left[ e\mathbf{E} + \frac{e}{b} \left( \mathbf{u} \times \mathbf{B} \right) \right] = \frac{d\boldsymbol{\pi}}{dt} \tag{6.11}$$

### 6.5 Particle in Potential

**Hamilton's Equations (Approximate):**
$$\mathbf{u} = \left[1 + \frac{V}{mc^2}\right] \frac{\pi}{m} \tag{6.12}$$
$$\frac{c}{b} \frac{d\mathbf{p}}{d\tau} = -\nabla V - \nabla V \frac{V}{mcb} \tag{6.13}$$

**Non-Relativistic Limit:**
$$m\mathbf{a} = -\nabla V - \nabla V \frac{V}{mc^2} \tag{6.14}$$

**Critical Point (Classical Electron Radius):**
$$-\nabla\Phi - \nabla\Phi\frac{V}{mc^2} = 0 \quad \text{at} \quad r = r_0 \tag{6.15}$$

---

## 7. Quantum Theory - Hilbert Spaces

### 7.1 KS²[ℝⁿ] Hilbert Space

**Inner Product:**
$$(f,g) = \int_{\mathbb{R}^n \times \mathbb{R}^n} f(\mathbf{x}) g(\mathbf{y})^* d\mathbf{P}_{\lambda}(\mathbf{x}, \mathbf{y}) \tag{7.1}$$

**Measure Definition:**
$$d\mathbf{P}_{\lambda}(\mathbf{x}, \mathbf{y}) = \left[ \sum_{k=1}^{\infty} t_{\lambda}^{k} \mathcal{E}_{k}(\mathbf{x}) \mathcal{E}_{k}(\mathbf{y}) \right] d\mathbf{x} d\mathbf{y} \tag{7.2}$$

where $t_{\lambda}^k = \lambda^{k-1} e^{-\lambda} / (k-1)!$

**Norm:**
$$||f||_{KS^2} = \left[ \sum_{k=1}^{\infty} t_{\lambda}^k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|^2 \right]^{1/2} \tag{7.3}$$

### 7.2 Embedding Theorems

**L^p Embedding:**
$$L^p[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n] \quad \text{for} \quad 1 \leq p \leq \infty \tag{7.4}$$

**Test Function Embedding:**
$$\mathcal{D}[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n] \tag{7.5}$$

**Distribution Embedding:**
$$\mathcal{D}^*[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n] \tag{7.6}$$

### 7.3 Fourier Transform and Convolution

**Fourier Transform Extension:**
$$\mathfrak{F}(\cdot): KS^2[\mathbb{R}^n] \to KS^2[\mathbb{R}^n] \quad \text{(bounded linear operator)} \tag{7.7}$$

**Convolution Extension:**
$$\mathfrak{C}_f(\cdot): KS^2[\mathbb{R}^n] \to KS^2[\mathbb{R}^n] \quad \text{(bounded linear operator)} \tag{7.8}$$

---

## 8. Feynman Operator Calculus

### 8.1 HK-Integral Definition

**HK-Partition:**
$\mathbf{P} = \{t_0, \tau_1, t_1, \tau_2, \cdots, \tau_n, t_n\}$ where $t_i, t_{i+1} \in (\tau_{i+1} - \delta(\tau_{i+1}), \tau_{i+1} + \delta(\tau_{i+1}))$

**HK-Integral:**
$$Q[a,b] = (HK) \int_{a}^{b} H(t)dt \tag{8.1}$$

where $\left\| \sum_{i=1}^{n} \Delta t_i H(\tau_i) - Q[a, b] \right\| < \varepsilon$

### 8.2 Feynman Kernel

**Feynman Kernel:**
$$\mathbb{K}_{\mathbf{F}}[t, \mathbf{x}; s, B] = \int_{B} (2\pi i (t - s))^{-1/2} \exp\{i |\mathbf{x} - \mathbf{y}|^{2} / 2(t - s)\} d\mathbf{y} \tag{8.2}$$

**Kernel Composition:**
$$\mathbb{K}_{\mathbf{F}}[t, \mathbf{x}; s, B] = \int_{\mathbb{R}^n} \mathbb{K}_{\mathbf{F}}[t, \mathbf{x}; \tau, d\mathbf{z}] \mathbb{K}_{\mathbf{F}}[\tau, \mathbf{z}; s, B] \tag{8.3}$$

### 8.3 Path Integral

**Path Integral Definition:**
$$\int_{\mathbb{R}^{[0,t]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}\mathbf{x}(\tau); \mathbf{x}(0)] \psi[\mathbf{x}(0)] = \lim_{\lambda \to \infty} \int_{\mathbb{R}^{n[0,t]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}_{\lambda}\mathbf{x}(\tau); \mathbf{x}(0)] \psi[\mathbf{x}(0)] \tag{8.4}$$

**Free Particle Path Integral:**
$$\int_{\mathbb{R}^{3[s,t]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}\mathbf{x}(\tau) ; \mathbf{x}(s)] = \frac{1}{\sqrt{2\pi i(t-s)}} \exp\{i|\mathbf{x} - \mathbf{y}|^2 / 2(t-s)\} \tag{8.5}$$

### 8.4 Time-Ordered Operator Calculus

**Time-Ordered Product:**
$$T\{A(t_1)B(t_2)\} = \begin{cases} A(t_1)B(t_2) & t_1 > t_2 \\ B(t_2)A(t_1) & t_2 > t_1 \end{cases} \tag{8.6}$$

**Dyson Series:**
$$U(t,s) = 1 + \sum_{n=1}^{\infty} \left(-\frac{i}{\hbar}\right)^n \int_s^t dt_1 \int_s^{t_1} dt_2 \cdots \int_s^{t_{n-1}} dt_n \, T\{H_I(t_1)H_I(t_2)\cdots H_I(t_n)\} \tag{8.7}$$

---

## 9. Dual Relativistic Wave Equations

### 9.1 General Form

**Dual Schrödinger Equation:**
$$i\hbar\frac{\partial\Phi}{\partial\tau}=K\Phi=\left[\frac{H^2}{2mc^2}+\frac{mc^2}{2}\right]\Phi \tag{9.1}$$

### 9.2 Three Dual Wave Equations

**Equation 1 (Dual Dirac):**
$$i\hbar \frac{\partial \Psi}{\partial \tau} = \left\{ \frac{\boldsymbol{\pi}^2}{2m} + \boldsymbol{\beta} V + mc^2 - \frac{e\hbar \boldsymbol{\Sigma} \cdot \mathbf{B}}{2mc} + \frac{V\boldsymbol{\alpha} \cdot \boldsymbol{\pi}}{mc} - \frac{i\hbar \boldsymbol{\alpha} \cdot \nabla V}{2mc} + \frac{V^2}{2mc^2} \right\} \Psi \tag{9.2}$$

**Equation 2 (Square-Root Form 1):**
$$i\hbar \frac{\partial \Psi}{\partial \tau} = \left\{ \frac{\pi^2}{2m} - \frac{e\hbar \mathbf{\Sigma} \cdot \mathbf{B}}{2mc} + mc^2 + \frac{V^2}{2mc^2} \right\} \Psi \tag{9.3}$$
$$+ \frac{V\beta \sqrt{c^{22} - ec\hbar \mathbf{\Sigma} \cdot \mathbf{B} + m^2c^4}}{2mc^2} \Psi \tag{9.4}$$
$$+ \frac{\beta \sqrt{c^2 \pi^2 - ec\hbar \mathbf{\Sigma} \cdot \mathbf{B} + m^2c^4}}{2mc^2} V\Psi \tag{9.5}$$

**Equation 3 (Square-Root Form 2):**
$$i\hbar \frac{\partial \Psi}{\partial \tau} = \left\{ \frac{\pi^2}{2m} + \beta V + mc^2 - \frac{e\hbar \mathbf{\Sigma} \cdot \mathbf{B}}{2mc} + \frac{V^2}{2mc^2} \right\} \Psi \tag{9.6}$$

### 9.3 Free Particle Limit

**Free Particle Equation:**
$$i\hbar\frac{\partial\Psi}{\partial\tau}=\left\{\frac{\mathbf{p^2}}{2m}+mc^2\right\}\Psi \tag{9.7}$$

---

## 10. Dual Dirac Equation and g-Factor

### 10.1 Dirac Equation

**Standard Dirac Equation:**
$$i\hbar \frac{\partial \Psi}{\partial t} = \left[ c\boldsymbol{\alpha} \cdot \mathbf{p} + \boldsymbol{\beta} mc^2 + V \right] \Psi \tag{10.1}$$

**Two-Component Form:**
$$i\hbar \frac{\partial \psi}{\partial t} - (mc^2 + V) \psi = c (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) \varphi \tag{10.2}$$
$$i\hbar \frac{\partial \varphi}{\partial t} + (mc^2 - V) \varphi = c (\sigma \cdot \pi) \psi \tag{10.3}$$

### 10.2 Dual Dirac Hamiltonian

**Dual Dirac Hamiltonian:**
$$K_D = \frac{H_D^2}{2mc^2} + \frac{mc^2}{2} = \frac{\pi^2}{2m} + V - \frac{e\hbar\Sigma \cdot \mathbf{B}}{2mc} + mc^2 + \frac{V_0^2}{2mc^2} \tag{10.4}$$

### 10.3 Eigenvalue Problem

**General Eigenvalue Equation:**
$$E\Psi = \left\{ \frac{\pi^2}{2m} + \beta V_0 + mc^2 - \frac{e\hbar\Sigma \cdot \mathbf{B}}{2mc} + \frac{V_0\alpha \cdot \pi}{mc} - \frac{i\hbar\alpha \cdot \nabla V_0}{2mc} + \frac{V_0^2}{2mc^2} \right\} \Psi \tag{10.5}$$

**Two-Component Form:**
$$E\psi_{1} = \left\{ \frac{\boldsymbol{\pi}^{2}}{2m} + V + mc^{2} - \frac{e\hbar\boldsymbol{\sigma}\cdot\mathbf{B}}{2mc} + \frac{V_{0}^{2}}{2mc^{2}} \right\} \psi_{1} \tag{10.6}$$
$$+ \left\{ \frac{V_{0}\boldsymbol{\sigma}\cdot\boldsymbol{\pi}}{mc} - \frac{i\hbar\boldsymbol{\sigma}\cdot\nabla V_{0}}{2mc} \right\} \psi_{2} \tag{10.7}$$

$$E\psi_{2} = \left\{ \frac{\boldsymbol{\pi}^{2}}{2m} - V + mc^{2} - \frac{e\hbar\boldsymbol{\sigma}\cdot\mathbf{B}}{2mc} + \frac{V_{0}^{2}}{2mc^{2}} \right\} \psi_{2} \tag{10.8}$$
$$+ \left\{ \frac{V_{0}\boldsymbol{\sigma}\cdot\boldsymbol{\pi}}{mc} - \frac{i\hbar\boldsymbol{\sigma}\cdot\nabla V_{0}}{2mc} \right\} \psi_{1} \tag{10.9}$$

### 10.4 Lower Component Relation

**Exact Relation:**
$$\psi_2 = c \left[ \lambda - V_0 + mc^2 \right]^{-1} (\sigma \cdot \boldsymbol{\pi}) \,\psi_1 \tag{10.10}$$

**Approximate Relation (with cutoff):**
$$\psi_2 = \frac{c(\boldsymbol{\sigma} \cdot \boldsymbol{\pi})}{2mc^2 \left(1 + \frac{r_0}{2r}\right)} \psi_1 \tag{10.11}$$

### 10.5 Pauli Matrix Identities

**Sigma Product Identity:**
$$(\boldsymbol{\sigma} \cdot \mathbf{X}) (\boldsymbol{\sigma} \cdot \mathbf{Y}) = \mathbf{X} \cdot \mathbf{Y} + i \, \boldsymbol{\sigma} \cdot (\mathbf{X} \times \mathbf{Y}) \tag{10.12}$$

**Momentum Product:**
$$(\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) = \boldsymbol{\pi}^2 - \frac{e\hbar}{c} \boldsymbol{\sigma} \cdot \mathbf{B} \tag{10.13}$$

**Gradient-Momentum Product:**
$$(-i\hbar\boldsymbol{\sigma}\cdot\boldsymbol{\nabla}V_0)(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}) = -i\hbar\nabla V_0\cdot\boldsymbol{\pi} + i\,\boldsymbol{\sigma}\cdot(-i\hbar\boldsymbol{\nabla}V_0\times\boldsymbol{\pi}) \tag{10.14}$$

### 10.6 g-Factor Formula

**g-Factor:**
$$g_r = 2 \left[ 1 - \frac{4r_0}{(2r + r_0)} \right] \tag{10.15}$$

**Anomalous Hamiltonian:**
$$H_a = 2\left[1 - \frac{4r_0}{(2r + r_0)}\right] \mu_B \mathbf{s} \cdot \mathbf{B} \tag{10.16}$$

**Muon g-Factor:**
$$g_{\mu}^{a} = 2 \left[ 1 - \frac{4r_{0}^{\mu}}{(2r_{\mu} + r_{0}^{\mu})} \right] \tag{10.17}$$

**Proton g-Factor:**
$$g_{p}^{a} = -2 \left[ 1 - \frac{4r_{0}^{p}}{(2r_{p} + r_{0}^{p})} \right] \tag{10.18}$$

where $r_0^{\mu} = \frac{e^2}{m_{\mu}c^2}$ and $r_0^p = \frac{e^2}{m_pc^2}$

---

## 11. Isodual Numbers and Antimatter

### 11.1 Real Number Field

**Definition:**
The real numbers is a triplet $(\mathbb{R}, +, \cdot)$, which is a field, with 0 as the additive identity (i.e., $a + 0 = a$ for all $a \in \mathbb{R}$) and 1 as the multiplicative identity (i.e., $a \cdot 1 = a$ for all $a \in \mathbb{R}$).

### 11.2 Isodual Real Numbers

**Definition:**
The isodual real numbers $(\hat{\mathbb{R}}, +, *)$ is a field, with $\hat{0}$ as the additive identity (i.e., $\hat{a} + \hat{0} = \hat{a}$ for all $-a = \hat{a} \in \hat{\mathbb{R}}$) and $\hat{1} = -1$ as the multiplicative identity (i.e., $\hat{a} * \hat{1} = (-a)(-1)(-1) = \hat{a}$ for all $\hat{a} \in \hat{\mathbb{R}}$).

### 11.3 Isodual of Physical Quantities

**Isodual Relation:**
$$A + A = 0 \tag{11.1}$$

### 11.4 Antimatter Evolution

**Particle Evolution:**
$$i\hbar\frac{\partial\psi}{\partial\tau}=K\psi \tag{11.2}$$

**Conjugate Equation:**
$$-i\hbar\frac{\partial\psi^*}{\partial\tau} = K\psi^* \tag{11.3}$$

**Isodual Evolution:**
$$\hat{i} * \hat{\hbar} * \frac{\partial \psi^*}{\partial \hat{\tau}} = \hat{K} * \psi^* \tag{11.4}$$

---

## References

1. Gill, T. L., & Zachary, W. W. (2011). Two Mathematically Equivalent Versions of Maxwell's Equations. *Foundations of Physics*, 41, 99-128.

2. Gill, T. L., Zachary, W. W., & Lindesay, J. (2001). The Classical Electron Problem. *Foundations of Physics*, 31, 1299-1354.

3. Gill, T. L., & Ares de Parga, G. (2019). Foundations for QED II: Classical Theory. *Advanced Studies in Theoretical Physics*, 13(8), 337-377.

4. Gill, T. L., & Ares de Parga, G. (2019). Foundations for QED I: Mathematical. *Advanced Studies in Theoretical Physics*, 13(8), 337-377.

5. Gill, T. L., Ares de Parga, G., Morris, T., & Wade, M. (2021). Dual Relativistic Quantum Mechanics I. *Universal Journal of Physics and Application*, 3(1), 24-40.

---

## Appendix: Key Constants and Parameters

| Symbol | Name | Value |
|--------|------|-------|
| $c$ | Speed of light | $2.99792458 \times 10^8$ m/s |
| $\hbar$ | Reduced Planck constant | $1.0545718 \times 10^{-34}$ J·s |
| $e$ | Elementary charge | $1.60217663 \times 10^{-19}$ C |
| $m_e$ | Electron mass | $9.10938356 \times 10^{-31}$ kg |