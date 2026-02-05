# Elastic Network Model (ANM) for GPCR Activation

This repository contains a Python implementation of a state-dependent Elastic Network Model (ANM) designed to analyze the conformational transition of G-Protein Coupled Receptors (GPCRs) from an inactive to an active state.

## Project Overview
Standard ANM assumes a universal spring constant for all interactions. However, biological transitions often involve specific "breaking" of contacts. This script implements a **guided ANM** where spring constants are modulated based on the structural differences between two known states (Active vs. Inactive).

**Key Features:**
- **State-Aware Physics:** Identifies contacts that break during activation and exponentially "softens" their corresponding springs.
- **Manual Hessian Construction:** Demonstrates explicit construction of the Hessian matrix ($3N \times 3N$) from pairwise potentials.
- **ProDy Integration:** Leverages the ProDy library for PDB parsing and trajectory output while injecting custom physics.

## Dependencies
* Python 3.x
* `numpy`
* `prody`

## Usage
1. Ensure your input PDB files (`inactive.pdb` and `active.pdb`) contain matching residue indices.
2. Run the script:

```bash
python gANM.py
