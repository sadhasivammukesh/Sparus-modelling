# Sparus-modelling
# SPARUS Modeling for Underwater Robotics

## Overview
This repository contains implementations and simulations related to **SPARUS underwater vehicle modeling**, focusing on hydrodynamic modeling, added mass calculations, drag force analysis, and control system validation. The objective is to develop accurate models for underwater vehicle dynamics to ensure stable and efficient navigation.

## Features
- **Hydrodynamic Modeling:** Calculation of added mass, drag forces, and buoyancy effects.
- **Mass Matrix Analysis:** Understanding mass distribution and coupling effects in underwater vehicles.
- **MATLAB/Simulink Simulations:** Validating dynamic behavior under various flow conditions.
- **Added Mass Computation:** Applying Slender Body Theory and Lamb’s K-factors to estimate added mass.
- **Drag Force Modeling:** Deriving drag matrices for different components, including hull, antenna, and thrusters.
- **Force and Torque Analysis:** Evaluating the impact of external forces on vehicle stability.
- **Spar and Hull Interaction Studies:** Assessing the contribution of different structural elements to hydrodynamic forces.

## Installation
### Prerequisites
Ensure you have the following dependencies installed:
- **MATLAB/Simulink** with Simscape Multibody Toolbox.
- **Python (optional for data analysis):** NumPy, SciPy, Matplotlib.

## Usage
### Running the Simulink Model
1. Open `SPARUS_model.slx` in MATLAB.
2. Run the simulation and analyze the force and torque outputs.
3. Modify environmental parameters (current speed, depth) to study their effects.
4. Compute **mass and added mass matrices** using provided functions.
5. Evaluate drag matrices by adjusting flow parameters.
