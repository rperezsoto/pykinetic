"""
This module contains several utilities grouped by usage:

   *   The 'Pykinetic Utilities' are functions with usage in the scripts that 
       accompany this library.
   *   Class Specializations. Subclasses of the ChemicalSystem used for the
       scripts pykinetic-model and pykinetic-scan.
"""

import warnings
import math
from itertools import chain

from .classes import ChemicalSystem, Reaction, Energy, DiffusionTS

########################## Pykinetic Utilities #################################
def write_indexfile(chemsys,file,withoutTS=True,isrelative=False):
    """
    Writes a File that summarizes both Compounds and Reactions Files as
    they were modeled. Simple division of this file should lead to two files
    that should be able reproduce at least the function and the Jacobian.

    Parameters
    ----------
    chemsys : ChemicalSystem
        A ChemicalSystem or subclass to output
    FilePath : str
        A valid filepath for the IndexFile
    """
    Out = []
    Out.append('#### compounds ####')
    for compound in chemsys.compounds:
        key, label, energy = compound.key, compound.label, compound.energy
        Out.append(f'{key})    {label}    {energy}')
    Out.append('#### reactions ####')
    if withoutTS:
        for reaction in chemsys.reactions:
            key = reaction.key
            if isrelative:
                energy = reaction.AE.as_unit(chemsys.unit)
            else:
                energy = reaction.TS.energy
            Out.append(f'{key})    {reaction}    !{energy}')
    else:
        for reaction in chemsys.reactions:
            key = reaction.key
            Out.append(f'{key})    {reaction}    !{reaction.TS.label}')
        for TS in chemsys.transitionstates: 
            Out.append(f'{TS.label}    {TS.energy.as_unit(chemsys.unit)}')
    with open(file,'w') as F :
        F.write('\n'.join(Out))
        F.write('\n')
def calc_standard_state_correction(T):
    # constants from "The NIST Reference on Constants, Units, and Uncertainty"
    R_SI = 8.314462618    # J/(mol K)
    R_atm = 0.0820573661  # atm L / (mol K)
    return Energy(R_SI*T*math.log(R_atm*T),'J/mol')

######################### Class Specializations ################################
class BiasedChemicalSystem(ChemicalSystem):
    """
    A specialization used to directly apply a constant bias to the energy of
    to all the species (as the Standard State Correction does). It is applied to
    compounds and TS except from the diffusion TSs.

    Parameters
    ----------
    bias : float, Energy
        bias to the energy of each species (the default is 0.0).
    T : float
        Temperature in K (the default is 298.15).
    unit : str
        (the default is 'kcal/mol').
    CorrUnit: str
        (the default is 'kcal/mol').
    """
    def __init__(self,bias=0.0,bias_unit='kcal/mol',T=298.15,unit='kcal/mol'):
        super().__init__(T,unit)
        if isinstance(bias,Energy):
            self._bias = bias
        else: 
            self._bias = Energy(bias,bias_unit)
    def apply_bias(self):
        for compound in self.compounds: 
            compound.energy = compound.energy + self.bias
        for TS in self.transitionstates: 
            if not hasattr(TS,'barrier'):
                TS.energy = TS.energy + self.bias
    def remove_bias(self):
        for compound in self.compounds: 
            compound.energy = compound.energy - self.bias
        for TS in self.transitionstates:
            if not hasattr(TS,'barrier'): 
                TS.energy = TS.energy - self.bias
    def change_bias(self,bias):
        self.remove_bias()
        if isinstance(bias,Energy):
            self._bias = bias
        else:
            warnings.warn("Guessing energy unit from previous bias' unit")
            self._bias = Energy(bias,self._bias.unit)
        self.apply_bias()
    
    @property
    def bias(self):
        return self._bias
    @bias.setter
    def bias(self,other):
        self.change_bias(other)
class ScannableChemicalSystem(BiasedChemicalSystem):
    """
    A specialization of BiasedChemicalSystem that provides the methods to apply
    dinamically change the energy of all the compounds and TS flagged as 
    scannable.

    Parameters
    ----------
    scan : float, Energy
        bias to the energy of each species (the default is 0.0).
    T : float
        Temperature in K (the default is 298.15).
    bias : float, Energy
        bias to the energy of each species (the default is 0.0).
    T : float
        Temperature in K (the default is 298.15).
    unit : str
        (the default is 'kcal/mol').
    CorrUnit: str
        (the default is 'kcal/mol').

    """
    def __init__(self,scan=0.0,scan_unit='kcal/mol',
                 bias=0.0,bias_unit='kcal/mol',T=298.15,unit='kcal/mol'):
        super().__init__(bias,bias_unit,T,unit)
        if isinstance(scan,Energy):
            self._scan = scan
            self._scan_unit = scan.unit
        else: 
            self._scan = Energy(scan,scan_unit)
            self._scan_unit = scan_unit

    def apply_scan(self):
        for item in chain(self.compounds,self.transitionstates): 
            if item.scannable:
                if not hasattr(item,'barrier'): 
                    item.energy = item.energy + self.scan
                else: 
                    item.barrier = item.barrier + self.scan
    def remove_scan(self):
        for item in chain(self.compounds,self.transitionstates): 
            if item.scannable:
                if not hasattr(item,'barrier'): 
                    item.energy = item.energy - self.scan
                else: 
                    item.barrier = item.barrier - self.scan
    def change_scan(self,scan):
        self.remove_scan()
        if isinstance(scan,Energy):
            self._scan = scan
        else:
            self._scan = Energy(scan,self._scan_unit)
        self.apply_scan()
    
    @property
    def scan(self):
        return self._scan
    @scan.setter
    def scan(self,other):
        self.change_scan(other)

