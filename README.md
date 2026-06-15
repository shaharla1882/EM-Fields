# Electromagnetic Fields - 2D Capacitor Numerical Simulation

## Overview
This repository contains a numerical simulation for solving the 2D Laplace equation ($\nabla^2 \Phi = 0$) to determine the electric potential and electric field distribution within a complex capacitor structure. The solution is implemented using the Finite Difference Method (FDM) alongside a relaxation iterative process.

This project was developed as part of the Electromagnetic Fields course at Ben-Gurion University of the Negev (BGU).

## Repository Structure
* `EM_Fields_Computer_Assignment.m`: The main MATLAB script. It handles the grid initialization, application of boundary conditions, the iterative FDM solver, electric field derivation, and capacitance computation using Gauss's Law.
* `EM_computer_final_project.pdf`: The final project report (in Hebrew). It includes the theoretical background, mathematical derivations, spatial plots of the potential and electric field, and a convergence analysis across different error tolerances.

## Geometry & Boundary Conditions
The system models a 2D capacitor consisting of two conductive electrodes ($C_1$ and $C_2$) in a vacuum/air dielectric ($\varepsilon = \varepsilon_0$):
* **Outer Electrode ($C_1$):** A $12\text{m} \times 12\text{m}$ square boundary, held at a constant potential of $\Phi = -1\text{V}$.
* **Inner Electrode ($C_2$):** A symmetric geometric structure at the center (defined by $w_{in} = 4\text{m}$ and $\Delta = 1\text{m}$), held at an equipotential of $\Phi = 1\text{V}$.

## Methodology
1. **Grid Generation:** Discretization of the continuous spatial domain into a uniform grid.
2. **Iterative Solver:** Numerical computation of the potential at internal nodes using the average of four neighboring points, derived from the Taylor series expansion of the Laplacian operator:
   $$\Phi_{i,j} = \frac{\Phi_{i+1,j} + \Phi_{i-1,j} + \Phi_{i,j+1} + \Phi_{i,j-1}}{4}$$
   The loop continues until the maximum difference between successive iterations falls below a predefined tolerance.
3. **Electric Field Derivation:** Computed numerically via central differences:
   $$\vec{E} = -\nabla\Phi = -\left(\frac{\partial\phi}{\partial x}\hat{x} + \frac{\partial\phi}{\partial y}\hat{y}\right)$$
4. **Capacitance Calculation:** Applying the integral form of Gauss's Law over a closed contour around the electrodes to find the enclosed charge per unit length ($Q$), and subsequently the capacitance:
   $$C = \frac{Q}{\Phi_{C2} - \Phi_{C1}}$$

## Convergence Analysis
The simulation evaluates the impact of the error tolerance on the accuracy of the computed capacitance ($C$ [F/m]). The results demonstrate a clear convergence trend as the grid calculations become more precise:
* **Tolerance = $10^{-1}$:** $C \approx 5.28 \times 10^{-10}\text{ F/m}$
* **Tolerance = $10^{-3}$:** $C \approx 6.27 \times 10^{-11}\text{ F/m}$
* **Tolerance = $10^{-5}$:** The system converges to a stable value of $C \approx 3.007 \times 10^{-11}\text{ F/m}$

## Usage
1. Open MATLAB and navigate to the repository directory.
2. Run the `EM_Fields_Computer_Assignment.m` script.
3. The script will execute the numerical analysis, generate 2D spatial plots for the potential distribution and the electric field vectors, and output the computed capacitance and charge values to the Command Window.

---
**Authors:** Shahar Lavi & Tomer Amran
