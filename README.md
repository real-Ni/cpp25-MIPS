# Active Matter Simulations: ABP, Run-and-Tumble, and Vicsek Models

## Overview

This project explores minimal models of active matter using numerical simulations. The aim is to understand how simple microscopic rules—self-propulsion, noise, and local interactions—give rise to collective phenomena such as persistent motion, clustering, and large-scale ordering. Three canonical models are studied: Active Brownian Particles (ABPs), one-dimensional run-and-tumble particles, and the Vicsek model.

---

## Models Studied

### 1. Active Brownian Particles (ABPs)
- Self-propelled particles moving at constant speed in two dimensions  
- Orientations evolve via rotational diffusion  
- Periodic boundary conditions  
- Repulsive interactions included at finite density  
- Mean-square displacement (MSD) used to validate dynamics  
- At high activity and density, clustering and motility-induced phase separation (MIPS) are observed  

### 2. Run-and-Tumble Particles (1D)
- Ballistic motion interrupted by stochastic tumbling events  
- Simulated using Gillespie’s algorithm for exact continuous-time dynamics  
- Two-particle simulations reveal jamming and effective interactions  
- Many-particle simulations show clustering driven by persistence  
- Comparison with symmetric random walks highlights the role of activity  

### 3. Vicsek Model
- Self-propelled particles align with local neighbors under angular noise  
- Collective ordering quantified using the velocity order parameter  
- Density and noise sweeps reveal a transition from disordered motion to collective alignment  
- Spatial snapshots show clustering and band formation  

---

## Key Results

- Persistent motion leads to a ballistic-to-diffusive crossover at the single-particle level  
- Repulsive interactions combined with activity generate clustering without attraction  
- Effective interactions emerge dynamically in run-and-tumble systems  
- The Vicsek model exhibits a nonequilibrium phase transition controlled by noise and density  
- Collective motion and spatial inhomogeneities arise from minimal interaction rules  

---

## Project Structure

- `*.cpp` — Core C++ simulation codes  
- `*.ipynb` — Data analysis and plotting notebooks  
- `figures/` — Simulation snapshots and plots  
- `report/` — Theory, derivations, and results write-up  

---

## Purpose

This project serves as a computational study of nonequilibrium statistical physics, focusing on how collective behavior emerges in active systems. It is intended both as a learning exercise and as a foundation for further investigations into active matter, phase separation, and collective dynamics.

---

## Notes

- All simulations are performed with periodic boundary conditions  
- Parameters are swept systematically to ensure fair comparisons  
- Derivations of theoretical limits are included in the report  
