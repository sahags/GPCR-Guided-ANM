"""
Elastic Network Model for GPCR Activation. 
Instead of a universal spring force constant of 1, a distance-dependent force constant is used here 
to soften specific springs that differ between active and inactive structures. 
The distance-dependent softening ebables inactive to active transition.
"""

import numpy as np
import prody as pd
import sys
import os

class gpcrANM:
    """
    ANM with spring constants modulated by structural differences between active and inactive states
    """
    
    def __init__(self, inactive_pdb, active_pdb, cutoff=15.0, softening_factor=0.1, distance_threshold=2.0):
        """
        inactive_pdb : location inactive structure
        active_pdb : location of the active structure
        cutoff : contact cutoff distance (Angstroms)
        softening_factor : how much to reduce spring constants (0-1, where 0.1 = 90% softer)
        distance_threshold : Minimum Ca distance change to trigger softening (Angstroms)
        """

        self.cutoff = cutoff
        self.softening_factor = softening_factor
        self.distance_threshold = distance_threshold
        
        print(f"Loading structures: {inactive_pdb} & {active_pdb}")
        # Check if pdb files exist
        if not os.path.exists(inactive_pdb) or not os.path.exists(active_pdb):
            raise FileNotFoundError(f"Could not find input PDB files.")

        inactive = pd.parsePDB(inactive_pdb)
        active = pd.parsePDB(active_pdb)
        print("inactive and active structures are loaded.")

        self.inactive_ca = inactive.select('name CA')
        self.active_ca = active.select('name CA')
        
        self.align_structures() # if submitted structures are already aligned this step cant be commented out.
        print("Structures aligned.")

        self.distance_changes = self.calculate_distance_changes()
        print(f"Distance changes calculated for {len(self.distance_changes)} pairs.")

        print("Building resolved state-aware  force constants...")
        self.force_constants = self.build_force_constants()
        
        print("Performing ANM with revised force constants...")
        self.anm = self.build_anm()
        
    def align_structures(self):
        """Align active to inactive structure using ProDy"""
        transformation = pd.calcTransformation(self.active_ca, self.inactive_ca)
        transformation.apply(self.active_ca)
        
    def calculate_distance_changes(self):
        """
        Calculate pairwise distance changes between inactive and active states
        Returns: dict mapping (i,j) --> distance_change
        """
        coords_inactive = self.inactive_ca.getCoords()
        coords_active = self.active_ca.getCoords()
        
        n_atoms = len(coords_inactive)
        distance_changes = {}
        
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                # Distance in inactive state
                r_inactive = np.linalg.norm(coords_inactive[i] - coords_inactive[j])
                
                # Only consider contacts in inactive state which is lower than the 15 A cut-off distance
                if r_inactive < self.cutoff:
                    # Distance in active state
                    r_active = np.linalg.norm(coords_active[i] - coords_active[j])
                    
                    
                    delta_r = r_active - r_inactive
                    
                    distance_changes[(i, j)] = {
                        'delta': delta_r,
                        'inactive_dist': r_inactive,
                        'active_dist': r_active
                    }
        
        return distance_changes
    
    def build_force_constants(self):
        """
        Build force constant matrix with softening based on distance changes
        """
        n_atoms = self.inactive_ca.numAtoms()
        force_constants = np.ones((n_atoms, n_atoms))  # default force constant=1.0
        
        for (i, j), change_info in self.distance_changes.items():
            delta_r = change_info['delta']

            #if the distance difference between the contant is higher than threshold(0.2), soften the spring proportionally 
            #using an exponential decay function.

            if abs(delta_r) > self.distance_threshold:
                softening = np.exp(-abs(delta_r) / 5.0) 
                gamma_ij = self.softening_factor * softening
                
                force_constants[i, j] = gamma_ij
                force_constants[j, i] = gamma_ij
                
                # Printing major softening events
                if abs(delta_r) > 5.0:
                    resnum_i = self.inactive_ca.getResnums()[i]
                    resnum_j = self.inactive_ca.getResnums()[j]
                    print(f"Softening {resnum_i}-{resnum_j}: Δr={delta_r:.1f}Å, γ={gamma_ij:.3f}")
        
        return force_constants
    
    def build_anm(self):
        """
        Build ANM with revised force constants-manually constructing hessian
        """
        coords = self.inactive_ca.getCoords()
        n_atoms = len(coords)
        
        hessian = np.zeros((3*n_atoms, 3*n_atoms)) #initiate empty hessian 
        
        # Build Hessian with custom force constants
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                if (i, j) in self.distance_changes:
                    dist_info = self.distance_changes[(i, j)]
                    r_ij = dist_info['inactive_dist']
                    
                    gamma_ij = self.force_constants[i, j]
                    
                    diff = coords[i] - coords[j] #vector from atom j to atom i
                    diff_norm = diff / r_ij #unit vector
                    
                    # Build 3x3 block for this interaction (outer product)
                    block = gamma_ij * np.outer(diff_norm, diff_norm)
                    
                    # Fill Hessian (super-blocks)
                    for a in range(3):
                        for b in range(3):
                            # Off-diagonal super-blocks (interaction)
                            hessian[3*i+a, 3*j+b] = -block[a, b]
                            hessian[3*j+b, 3*i+a] = -block[a, b]
                            
                            # Diagonal super-blocks (self-interaction or stiffness)
                            hessian[3*i+a, 3*i+b] += block[a, b]
                            hessian[3*j+a, 3*j+b] += block[a, b]
        
        # Create ProDy ANM object
        anm = pd.ANM('transition-guided')
        anm.buildHessian(self.inactive_ca, cutoff=self.cutoff)
        
        # Replace Hessian with our custom one
        anm._hessian = hessian
        anm._n_atoms = n_atoms
        anm._dof = 3 * n_atoms
        
        # Calculate modes
        anm.calcModes(n_modes=20)
        pd.writeNMD('transition.nmd', anm[:20], self.inactive_ca)
        pd.writePDB('transition.pdb', self.inactive_ca)
        return anm
    
if __name__ == "__main__":
    
   #Example usage for GPCR inactive to active transition
   #This code requires two pdb files with exact same number of Ca-atom.
   #Replace the filenames below with your specific PDBs.
    
    inactive_pdb = "2rh1_clean.pdb"  
    active_pdb = "3sn6_clean.pdb"      
    print("--- Starting GPCR Guided ANM ---")

    try:
        g_anm = gpcrANM(
                inactive_pdb=inactive_pdb,
                active_pdb=active_pdb,
                cutoff=15.0,
                softening_factor=0.1,
                distance_threshold=2.0
        )
        print("Analysis Complete. Outputs: transition.nmd, transition.pdb")

    except FileNotFoundError:
        print("Error: PDB files not found in the current directory.")
        print(f"Please ensure '{inactive_file}' and '{active_file}' exist.")
        print("Or update the file paths in the __main__ block.")
    except Exception as e:
        print(f"An error occurred during execution: {e}")

