"""
This module contains several utilities grouped by usage:

   *   General utilities.
   *   The pykinetic Utilities are Classes and functions without chemical
       meaning but necessary to write out the classes of the 'pyKinetic' module.
   *   Class Specializations. Subclasses of the ChemicalSystem used for the
       scripts pykinetic-scan and pykinetic-scanTS. Note that New Models should
       not be place here.
"""

import os
import errno
from pkg_resources import resource_filename
from collections import UserDict

from .Classes import ChemicalSystem, Reaction, Energy, DiffusionTS

############################ General Utilities #################################

def mkdir_p(path):
    """
    Same functionality as mkdir -p in bash console
    """
    # Copied from stackoverflow:
    # stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

########################## Pykinetic Utilities #################################
class Parameters(UserDict):

    @classmethod
    def read_from(cls,file):
        """
        Reads a file line by line and transforms it into a dictionary.

        Parameters
        ----------
        File : str
            filepath

        Returns
        -------
        dict
            dictionary containing first column as key and an empty string or
            the rest of columns as a string
        """
        parameters = cls()
        with open(file,'r') as F:
            for line in F:
                if line.strip():
                    Aux = line.strip().split(maxsplit=1)
                    key = Aux[0]
                    try:
                        item = Aux[1]
                    except ValueError:
                        item = ''
                    finally:
                        parameters[Aux[0]] = Aux[1]
        return parameters
class SimulationParameters(Parameters):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self['tfin'] = self.get('tfin','')
        self['dt'] = self.get('dt','')
        self['report_t'] = self.get('report_t','')
    def read_concentrations(self):
        compounds_mark = ';'
        values_mark  = ','
        text = self.get('concentrations','')
        self['concentrations'] = dict()
        if text:
            for compound in text.strip().split(compounds_mark): 
                key,val = compound.split(values_mark)
                self['concentrations'][key.strip()] = val.strip()
    @classmethod
    def read_from(cls,file):
        parameters = super().read_from(file)
        parameters.read_concentrations()
        return parameters
class ConvergenceParameters(Parameters):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        # default values
        self['rtol'] = self.get('rtol','1E-6')
        self['atol'] = self.get('atol','1E-12')

    @property
    def variables(self):
        return [key for key in self]
    
    def as_str(self,sep=','):
        parameters = [f'{key}={val}' for key,val in self.items()] 
        return sep.join(parameters)

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
                energy = reaction.TS.energy.as_unit(chemsys.unit)
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

def DotFile(ChemSys,Rank=False):
    """
    Returns a list of lines with a graph that relates reactants with
    reactions and reactions with products for each elemental step in .dot format

    Parameters
    ----------
    ChemSys : ChemicalSystem
        An updated ChemicalSystem instance
    Rank : bool
        Wether to align or not the different initial reactants
        (the default is False).
    """
    Out = []
    Out.append('strict Digraph {')
    Out.append('\tconcentrate=true')
    for reaction in ChemSys.reactions:
        Rs = ';'.join(['"{}"'.format(i.label) for i in reaction.reactants])
        Rs = '{'+ Rs + '}'
        Ps = ';'.join(['"{}"'.format(i.label) for i in reaction.products])
        Ps = '{'+ Ps + '}'
        Out.append('\t{}-> "{}" -> {};'.format(Rs,reaction,Ps))
    rs = ';'.join(['"{}"'.format(i) for i in ChemSys.reactions])
    if Rank:
        Out.append('\t{rank = same;'+'{}'.format(rs)+'}')
    Out.append("}")
    return Out

###################### Specializations of Classes ##############################
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
        super().__init__(T)
        try:
            bias.to_unit(bias_unit)
        except AttributeError:
            bias = Energy(bias,bias_unit)
        finally:
            self._bias = bias
        self.unit = unit
    def apply_bias(self):
        for compound in self.compounds: 
            compound.energy = compound.energy + self.bias
        for TS in self.transitionstates: 
            if not hasattribute(TS,'barrier'):
                TS.energy = TS.energy + self.bias
    def remove_bias(self):
        for compound in self.compounds: 
            compound.energy = compound.energy - self.bias
        for TS in self.transitionstates: 
            if not isinstance(TS,DiffusionTS):
                TS.energy = TS.energy - self.bias
    def change_bias(self,bias):
        self.remove_bias()
        self._bias = bias
        self.apply_bias()
    
    @property
    def bias(self):
        return self._bias
    @bias.setter
    def bias(self,other):
        self.remove_bias()
        self._bias = bias
        self.apply_bias()

class MutableSystem(BiasedChemicalSystem):
    """
    A specialization of a CorrChemSys used to mutate the energy of the
    species, minima and/or TS, a fixed amount.
    """
    def MutateEnergy(self,Value,Model='ElementalStep'):
        for Comp in self.compounds:
            Comp.energy += Value
        reactions = []
        for SR in self.sreactions:
            # Recalculate the TS energy
            if SR.RType != '<d>':
                SR.energy += Value
            SR.calc_TS()
            # Calculate the Barriers for the child reactions
            NewData = SR.to_Reaction(Model)
            for child,data in zip(SR.child_reactions,NewData):
                child.set_AE_from(data)
    def MutateTSs(self,Value,Model='ElementalStep'):
        for SR in self.sreactions:
            # Only Apply mutation to indicated reactions
            if SR.RType != '<s>':
                continue
            E_Ps = sum(P.energy for P in SR.products)
            E_Rs = sum(R.energy for R in SR.reactants)
            MaxVal = max([E_Ps,E_Rs])
            if SR.TS - MaxVal < Value:
                NewVal = SR.TS - MaxVal
            else:
                NewVal = Value
            if abs(NewVal) - 1E-14 < 0 :
                NewVal = Energy(0.0,NewVal._unit)
            SR.energy -= NewVal
            SR.calc_TS()
            # Reduce the barriers of the childs the same amount
            for child in SR.child_reactions:
                child.AE -= NewVal
                child.Calc_k()
