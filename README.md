# CoilRailSim

A highly realistic 1D coil rail simulation wrapped in a single web application. CoilRailSim uses electro-magnetic principles, perturbation methods and realistic corrections for modelling the dynamics of ferromagnetic probes accelerated by solenoids.

The full derivation and implemention are summarized in this [paper]().

## Features

- Solenoid magnetic potential modelling with fundamental Biot-Savart approach
- Respects induction laws with first order correction, inductance reacts to changes in permeability and probe movement
- Source voltage and currents react to inductance and precise resistance modelling
- Theoretical model for dipole magnetization within the probe with proper response functions
- Demagnetization correction for the dipole magnetization and coupling to solenoid field
- Integration with superior leap-frog algorithm which respects and conserves energy over long simulation times 

## Introduction

Electro-magnets serve multiple use-cases in engineering such as magnetic adhesion, EM damping, EM rail gun, stabilizers and more. In most cases we deal with a solenoid (cylindrical coil) which accelerates a probe, target or projectile along it's axis line which allows to abstract and simplify to a 1-dimensional problem. Given the solenoid material and geometry, as well as the target probe properties and environmental settings, this simulator can create a realistic, time-resolved analysis of the resulting target dynamics, electric properties, perturbation coupling by back-EMF, self-induction and power conversion efficiency. Simulations are extremely helpful when designing a coil under specific constraints as there
Except the accelerating force, all other forces, such as friction, are ignored.