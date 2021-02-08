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

from .Classes import ChemicalSystem, Reaction, Energy

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

class FileTemplate(object):
    """
    Object that represents the script that when executed runs the
    simulation. The parameters employed to construct the object are
    parsed and converted to strings.

    Parameters
    ----------
    Header : str
        File template where the Head of the script is located
        (the default is templates/Header.txt).
    Tail : str
        File template where the Tail of the script is located. Contains the code
        for the integration and the logic for the output generation (the default
        is templates/Tail.txt).
    Function : str
        File template where the function format of the script is located
        (the default is templates/Function.txt).
    Jacobian : str
        File template where the format of the Jacobian of the script is located
        (the default is templates/Jacobian.txt).
    Sim_Params : str
        File with the values for the parameters of the Simulation
        (the default is templates/Sim_Params.txt).
    Conv_Params : str
        File with the values for the convergence parameters of the numerical
        solver (the default is templates/Conv_Params.txt).

    Attributes
    ----------
    OFileName
    end
    T
    ks
    species
    massbalances
    rates
    partials
    Header
    Tail
    Function
    Jacobian
    Sim_Params
    Conv_Params
    txt

    """
    def __init__(self, Header=None, Function=None, Jacobian=None,
        Sim_Params=None, Conv_Params=None, Tail=None):
        templates_path = resource_filename('pykinetic','templates')
        # Assign Defaults
        if Header is None:
            Header = templates_path + '/Header.txt'
        if Tail is None:
            Tail = templates_path + '/Tail.txt'
        if Function is None:
            Function = templates_path + '/Function.txt'
        if Jacobian is None:
            Jacobian = templates_path + '/Jacobian.txt'
        if Sim_Params is None:
            Sim_Params = templates_path + '/Sim_Params.txt'
        if Conv_Params is None:
            Conv_Params = templates_path + '/Conv_Params.txt'
        self.OFileName = 'OutFile'
        self.end = '\n'
        self.fill = '{}' # Usefull for formatting
        self.T = self.ks = self.species = None
        self.massbalances = self.rates = None
        self.partials = None
        self._txt = '{0.Header}{0.Function}{0.Tail}'
        self.txt = ''
        # Read Files
        self.Header = self.ReadTemplate(Header)
        self.Tail = self.ReadTemplate(Tail)
        self.Function = self.ReadTemplate(Function)
        self.Jacobian = self.ReadTemplate(Jacobian)
        self._Conv_Params = self.ReadParams(Conv_Params)
        self.Sim_Params = self.ReadParams(Sim_Params)
        # Handle special cases
        self.Format_Cini()
        self.Format_Conv_Params()
    def __repr__(self):
        cls = type(self).__name__
        return '<{} Object>'.format(cls)

    def ReadTemplate(self,File):
        """
        Reads a file line by line and removes all ending empty lines.

        Parameters
        ----------
        File : str
            filepath

        Returns
        -------
        str
            string with the format ready to do a str.format(self)

        Raises
        -------
        RuntimeError
            If the file is empty

        """
        with open(File,'r') as F:
            Out = [line for line in F]
        while not Out[-1].strip():
            if not Out:
                raise RuntimeError('File {} is empty'.format(File))
            else:
                Out.pop(-1)
        return ''.join(Out)
    def ReadParams(self,File):
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
        Out = dict()
        with open(File,'r') as F:
            for line in F:
                if line.strip():
                    Aux = line.strip().split(maxsplit=1)
                    key = Aux[0]
                    try:
                        item = Aux[1]
                    except ValueError:
                        item = ''
                    finally:
                        Out[Aux[0]] = Aux[1]
        return Out
    def Format_Cini(self):
        """
        Parses the Initial Concentrations and separates them by lines.
        """
        Aux = []
        Sep_Comps = ';'
        Sep_Values = ','
        txt = 'xini[{0}] = {1}'
        items = self.Sim_Params['Conc_ini'].strip().split(Sep_Comps)
        for item in items:
            i,val = item.strip().split(Sep_Values)
            Aux.append(txt.format(i,val))
        self.Sim_Params['Conc_ini'] = '\n'.join(Aux)
    def Format_Conv_Params(self):
        """
        Formats self._Conv_Params dict to a single line set of kwargs and
        updates the self.Conv_Params.
        """
        txt = '{0}={1}'
        Out = []
        for key,val in self._Conv_Params.items():
            Out.append(txt.format(key,val))
        self.Conv_Params = ','.join(Out)
    def Add_Conv_Param(key,val):
        """
        Wrapper to add a certain Conv_Parameter. Forces str on parameters
        """
        self._Conv_Params[str(key)] = str(val)
        self.Format_Conv_Params()
    def Fill(self):
        """
        Fills the values of the template with the appropiated values and updates
        the 'txt' attribute.
        """
        # Construct the scheme
        self.txt = self._txt.format(self)
        # Fill in values
        self.txt = self.txt.format(self)
    def write_to(self,filepath):
        """
        Writes the contents of self.txt in a valid filepath
        """
        with open(filepath,'w') as F:
            F.write(self.txt)
    def update_with(self,ChemSys):
        """
        Uses a ChemicalSystem instance to update the 'species', 'T', 'ks',
        'rates','massbalances' and 'partials' attributes.
        """
        self.species = len(ChemSys.compounds)
        self.T = ChemSys.T
        self.ks = '\n'.join(ChemSys.k_Expr())
        self.rates = '\n'.join(ChemSys.Rates_Expr())
        self.massbalances = '\n'.join(ChemSys.MassBalances_Expr())
        self.partials = '\n'.join(ChemSys.Jacobian_Expr())

def Write_IndexFile(ChemSys,FilePath):
    """
    Writes a File that summarizes both Compounds and Reactions Files as
    they were modeled. Simple division of this file should lead to two files
    that should be able reproduce at least the function and the Jacobian.

    Parameters
    ----------
    ChemSys : ChemicalSystem
        A ChemicalSystem or subclass to output
    FilePath : str
        A valid filepath for the IndexFile
    """
    Name = FilePath
    Out = []
    Out.append('#### compounds ####')
    for i in ChemSys.compounds:
        Out.append('{0.key})\t{0.label}\t{0.energy: 0.6f}'.format(i))
    Out.append('#### reactions ####')
    oformat = '{0:})\t{1:}\t!{2:0.1f}'
    #riformat = '{0:})\t{3:}\t{2:}\t{1:}\t!{4:}'
    for i,reaction in enumerate(ChemSys.reactions):
        reaction.AE.to_unit('kcal/mol')
        AE = float(reaction.AE)
        Out.append(oformat.format(i,reaction,AE))
    with open(Name,"w") as F :
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
class CorrChemSys(ChemicalSystem):
    """
    A specialization used to directly apply an Energy Correction similar to
    Standard State Correction.

    Parameters
    ----------
    Corr : float, Energy
        Description of parameter `Corr` (the default is 0.0).
    T : float
        Temperature in K (the default is 298.15).
    unit : str
        (the default is 'kcal/mol').
    CorrUnit: str
        (the default is 'kcal/mol').
    """
    def __init__(self,Corr=0.0,T=298.15,unit='kcal/mol',CorrUnit='kcal/mol'):
        super(CorrChemSys, self).__init__(T)
        try:
            Corr.to_unit(CorrUnit)
        except AttributeError:
            Corr = Energy(Corr,CorrUnit)
        finally:
            self.Corr = Corr
        self.unit = unit
    def cadd(self,compound,update=False):
        """
        Patches the `ChemicalSystem.cadd` by adding the Corr to the
        energy of the compound added, before actually adding it
        """
        compound.energy = Energy(compound.energy,self.unit)
        compound.energy = compound.energy + self.Corr
        super(CorrChemSys, self).cadd(compound,update)
    def sradd(self,sreaction,update=True):
        """
        Patches the addition of sreactions so that it applies the correction
        to the energy of the sreaction except if it is of RType '<d>'. and then
        adds it normally.
        """
        if sreaction.RType == '<d>':
            sreaction.energy = Energy(sreaction.energy,'kcal/mol')
        else:
            sreaction.energy = Energy(sreaction.energy,self.unit)
            sreaction.energy += self.Corr
        super(CorrChemSys, self).sradd(sreaction,update)
class MutableSystem(CorrChemSys):
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
