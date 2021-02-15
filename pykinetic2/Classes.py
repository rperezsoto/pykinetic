"""
This module contains the core classes of the library: Energy, Compound,
SymbReaction, Reaction and Chemical System. It also contains the default model
for reactions: ElementalStep.
"""

from itertools import chain
from collections import Counter
from abc import abstractmethod
import warnings

import numpy

class Energy(object):
    """
    Objects that represents a number with an energy unit.

    Parameters
    ----------
    value : float
        Number or object with a __float__ (the default is 0.0).
    unit : str
        Energy unit in which the value parameter is provided see class variable
        _units for the different currently supported units
        (the default is 'hartree').

    Attributes
    ----------
    _units : list
        {}
    unit
    value

    """
    _units = 'hartree eV cm-1 kcal/mol kJ/mol J/mol K J Hz'.split()
    __doc__ = __doc__.format(_units)
    _Factors = (1.0, 27.2107, 219474.63, 627.503, 2625.5, 2625500.0,
                315777.0, 43.60E-19, 6.57966E15)
    _ConversionFactors = {key:val for key,val in zip(_units,_Factors)}

    __slots__ = ('_value','_conversion','_SI','_unit')

    def __init__(self,value=0.0,unit='hartree'):
        # Assume that value is in the specified unit
        self._conversion = type(self)._ConversionFactors[unit]
        self._SI = type(self)._ConversionFactors['J/mol']
        self.unit = unit
        self.value = float(value)
    def __repr__(self):
        cls = type(self).__name__
        return '<{0}:{1},{2}>'.format(cls,self.Val,self._unit)
    def __str__(self):
        return ' '.join([str(self.value),self.unit])
    def __format__(self,format_spec):
        return ' '.join([self.value.__format__(format_spec),self.unit])
    def __bool__(self):
        return self.value.__bool__()

    @property
    def value(self):
        return self._value * self._conversion
    @value.setter
    def value(self,other):
        self._value = float(other)/self._conversion

    @property
    def unit(self):
        return self._unit
    @unit.setter
    def unit(self,other):
        try:
            self._conversion = type(self)._ConversionFactors[other]
        except KeyError:
            raise NotImplementedError('Unit {} not implemented'.format(other))
        self._unit = other

    # Emulate Numeric Type
    def __add__(self,other):
        ''' Implements the "self + other" operation '''
        Out = Energy(unit=self.unit)
        try:
            Out._value =  self._value + other._value
        except AttributeError:
            Out.value = self.value + other
        finally:
            return Out
    def __sub__(self,other):
        ''' Implements the "self - other" operation '''
        return self.__add__(-other)
    def __mul__(self,other):
        ''' Implements the "self * other" operation '''
        if isinstance(other,type(self)):
            raise NotImplementedError('Energy*Energy not implemented')
        Out = Energy(self.value*other,unit=self.unit)
        return Out
    def __div__(self,other):
        ''' Defined for compatibility with Python2.7, Implements the
        "self/other" behaviour in python 2 '''
        return self.__truediv__(other)
    def __truediv__(self,other):
        ''' Implements the "self / other" operation '''
        Out = Energy(unit=self.unit)
        try:
            _ = other._value
        except AttributeError:
            Out.value = self.value/other
        else:
            Out = self._value/other._value
        return Out
    def __radd__(self,other):
        ''' Implements the "other + self" operation '''
        return self + other
    def __rsub__(self,other):
        ''' Implements the "self - other" operation '''
        return self.__add__(-other)
    def __rmul__(self,other):
        return self * other
    def __iadd__(self,other):
        ''' Implements the "+= other" operation '''
        try:
            self._value = self._value + other._value
        except AttributeError:
            self.value = self.value + other
        return self
    def __isub__(self,other):
        ''' Implements the "-= other" operation '''
        self.__iadd__(-other)
        return self
    def __imul__(self,other):
        ''' Implements the "*= other" operation '''
        self._value = self._value * other
        return self
    def __itruediv__(self,other):
        ''' Implements the "/= other" operation '''
        Out = self
        try:
            _ = other._value
        except AttributeError:
            self._value = self._value / other
        else:
            Out = self._value / other._value
        finally:
            return Out

    def __neg__(self):
        ''' Implements the "- self" operation '''
        return Energy(-self.value,self.unit)
    def __abs__(self):
        ''' Implements the abs() function behaviour '''
        return Energy(abs(self.value),self.unit)
    def __int__(self):
        ''' Implements the int() function behaviour '''
        return int(self.value)
    def __float__(self):
        ''' Implements the float() function behaviour '''
        return float(self.value)
    def __round__(self,ndigits):
        ''' Implements the round() function behaviour,
        Must return an integrer '''
        return round(self.value,ndigits)
    # Comparison methods
    def __eq__(self,other):
        ''' Implements the "self == other" behaviour'''
        # Implementation of equality test extracted from math.isclose
        atol = 1E-6
        rtol = 1E-9
        try:
            A,B = self._value, other._value
        except AttributeError:
            A,B = self.value, other
        test = max(rtol*max(abs(A),abs(B)),atol)
        if abs(A - B) == test:
            warnings.warn('Difference between numbers is close to tolerance')
        return abs(A - B) <= test
    def __ne__(self,other):
        ''' Implements the "self != other" behaviour'''
        return not(self == other)
    def __lt__(self,other):
        ''' Implements the "self < other" behaviour'''
        try:
            Out = self._value < other._value
        except AttributeError:
            Out = self.value < other
        return Out
    def __gt__(self,other):
        ''' Implements the "self > other" behaviour'''
        return not ((self < other) or (self == other))
    def __le__(self,other):
        ''' Implements the "self <= other" behaviour'''
        return not self > other
    def __ge__(self,other):
        ''' Implements the "self >= other" behaviour'''
        return not self < other


    def as_unit(self,unit):
        """Returns a string with the value expressed in the unit specified.

        Parameters
        ----------
        unit : str
            A valid string representing an energy unit, see `Energy._units`.

        Returns
        -------
        str

        """
        try:
            F = type(self)._ConversionFactors[unit]
        except KeyError:
            raise NotImplementedError('Unit {} not implemented'.format(unit))
        Val = str(self._value * F)
        return '{} {}'.format(Val,unit)
    def to_SI(self):
        """Returns the value in units of the SI.

        Returns
        -------
        float

        """
        return self._value*self._SI
    def to_unit(self,unit):
        """Changes the unit of the object.

        Parameters
        ----------
        unit : str
            A valid string representing an energy unit, see `Energy._units`.

        """
        self.unit = unit

class MassBalance(object):
    """
    Representation of a mass balance.

    Parameters
    ----------
    compound : Compound
        Compound of the targeted mass balance
    items : list
        list of Reactions with a net generation/consumption of 'compound' != 0.
    """
    def __init__(self,compound,items=None):
        self.compound = compound
        if items is None:
            self.items = []
        else:
            self.items = items
    def partial(self,compound):
        """
        Do the partial derivative of the mass balance with respect to 'compound'
        and return a JacobianElement object.

        Parameters
        ----------
        compound : Compound
            Compound with which the partial derivative is to be carried out.

        Returns
        -------
        JE
            JacobianElement
        """
        JE = JacobianElement(self.compound,compound)
        for coef,reaction in self.items:
            if compound in reaction:
                JE.items.append((coef,reaction))
        return JE
class JacobianElement(object):
    """
    Representation of an element of a jacobian matrix of the mass balances.

    Parameters
    ----------
    compound1 : Compound
        Compound of the targeted mass balance
    compound2 : Compound
        Compound of the partial derivative.
    items : list
        list of Reactions with a derivative with respect to compound2 != 0 for 
        that mass balance.

    """
    def __init__(self,compound1,compound2,items=None):
        self.compound1 = compound1
        self.compound2 = compound2
        if items is None:
            self.items = []
        else:
            self.items = items

class Compound(object):
    """
    Object representation of a chemical species.

    Parameters
    ----------
    label : str
        identifier/name of the species
    energy : energy
        relative or absolute free energy of the species.
    key : int
        integrer used to represent the species in the overall chemical system
        see `ChemicalSystem` (the default is None).

    """

    __slots__ = ("label","energy","key")

    def __init__(self,label,energy,key=None):
        self.label = label
        self.energy = energy
        self.key = key
    def __eq__(self,other):
        try:
            out = self.label == other.label
        except AttributeError:
            out = False
        finally:
            return out
    def __hash__(self):
        return hash(self.label)
    def __int__(self):
        return int(self.key)
    def __repr__(self):
        cls = self.__class__
        return f"{cls} <'{self.label}'>"
    def copy(self):
        """
        Creates a copy of the object with the key set as None.
        """
        return self.__class___(self.label,self.energy,key=None)
    def DotGraphRepr(self):
        """
        Returns a string for the construction of a .dot file
        """
        F0 = '"{}"'
        F1 = '"{0:}" [label="{0:}"]'.format(self.label)
        return F0.format(self.label)

class TransitionState(object):
    defaultname = 0 
    def __init__(self,energy,reactions=None,label=None,scannable=False):
        if label is None:
            label =  f'TS{self.defaultname:02d}'
            self.__class__.defaultname += 1
        self.label = label
        self.energy = energy
        if reactions is None:
            reactions = []
        self.reactions = reactions
        self.scannable = scannable
    def __eq__(self,other):
        try:
            out = self.label == other.label
        except AttributeError:
            out = False
        finally:
            return out
    def __hash__(self):
        return hash(self.label)
    def __repr__(self):
        cls = self.__class__.__name__
        return f"{cls} <'{self.label}'>"
    
    @property
    def energy(self):
        return self._energy
    @energy.setter
    def energy(self,other):
        self._energy = other

class DiffusionTS(TransitionState):
    defaultname = 0
    def __init__(self,barrier,reactions,label=None,scannable=False):
        self.barrier = barrier
        if label is None:
            label = f'DiffTS{self.defaultname:02d}'
            self.__class__.defaultname += 1
        self.label = label
        self.reactions = reactions
        self.scannable = scannable
    
    @property
    def energy(self):
        energy = max([r.reactants_energy for r in self.reactions])
        return self.barrier + energy
    @energy.setter
    def energy(self,other):
        energy = max([r.reactants_energy for r in self.reactions])
        if other < energy: 
            raise ValueError("Cannot set an energy lower than the reactants'")
        self.barrier = other - energy

class SymbReaction(object):
    """
    An object that serves as connection between the `Reaction` implementation
    and the user input reaction. Requires the existence of a dict-like in the
    namespace named 'ReactionTypes'.

    Parameters
    ----------
    Reaction : str
        Direct reading of the reaction from the user input
    EParse : str
        How the energy included in the input should be interpreted, either as
        'relative' to the reactants or as the 'absolute' energy of the TS
        (the default is "relative").

    Attributes
    ----------
    symb : str
        Symbolic user input representation
    energy : float
        see `EParse`.
    reactants
    products
    RType : string
        string representing the type of reaction
    IsUpdated : bool
        Marker that informs if the object has been read by a ChemicalSystem
        instance. (False on initialization)
    EParse
    child_reactions : Used to store a link to Reaction objects created from the
        information of each SymbReaction (Empty on initialization)
    TS : Used to store the energy of the 'TS'. (None upon initialization)

    """
    def __init__(self,Reaction,EParse='relative'):
        self.symb = Reaction.replace('\t',' '*4)
        self.EParse = EParse
        self.IsUpdated = False
        self.child_reactions = []
        self.TS = None
        Reaction,ActEnergy = Reaction.split("!")
        self.energy = float(ActEnergy.strip())
        marker = ' {} '.format # Auxiliary function for string matching
        for symbol in ReactionTypes:
            if marker(symbol) in self.symb:
                Rs,Ps = Reaction.split(symbol)
                self.RType = symbol
                break
        else:
            raise RuntimeError("No Reaction Symbol found in {}".format(Line))
        # Separe the different compounds
        key = marker('+')
        self.reactants = tuple(R.strip() for R in Rs.split(key) if R.strip())
        self.products = tuple(P.strip() for P in Ps.split(key) if P.strip())
        self._reactants = self.reactants
        self._products = self.products

    def __repr__(self):
        return self.symb
    def __str__(self):
        return self.__repr__()

    def reset_update(self):
        """
        Resets the state of the object as if it was not Updated with
        a Mapper.
        """
        if self.IsUpdated:
            self.reactants = self._reactants
            self.products = self._products
            self.TS = None
            self.IsUpdated = False
    def update(self,mapper):
        """
        Rebuilds the reactants and products lists using a dict-like object.
        to fill them with Compound Objects. After it, calculates the TS energy
        taking into account the type of Energy Parsing.

        Parameters
        ----------
        mapper : Mappable
        """
        # Map reactants and products
        self.reactants = [mapper[R] for R in self.reactants]
        self.products = [mapper[P] for P in self.products]
        # Now Calculate the TS Energy
        self.calc_TS()
        # Change State of the object
        self.IsUpdated = True
    def calc_TS(self):
        """
        Wrapper for recalculating the TS energy. Updates the value
        of the attribute 'TS'
        """
        E_Rs = sum(R.energy for R in self.reactants)
        E_Ps = sum(P.energy for P in self.products)
        self.TS = ReactionTypes[self.RType][self.EParse](E_Rs,E_Ps,self.energy)

    def child(self,child):
        """
        Adds a child to the list of child reactions of the Symbolic Reaction.

        Parameters
        ----------
        child : Reaction

        Raises
        -------
        TypeError
            If the child is not an instance of a Reaction or any Reaction
            subclass.

        """
        if not issubclass(child.__class__,Reaction):
            msg = '{} is not a instance of a Reaction class/subclasss'
            raise TypeError(msg.format(child))
        self.child_reactions.append(child)
    def reset_childs(self):
        """
        Similar to list.pop but returning the whole list of child_reactions.
        """
        out = self.child_reactions
        self.child_reactions = []
        return out

    def to_Reaction(self,Model='ElementalStep'):
        """
        Returns a list of tuples in which each tuple contains the information
        needed to initialize each Reaction Instance (See Classes.Reaction).
        Compound initialization is left to the ChemicalSystem.

        Parameters
        ----------
        Model : str
            Name of the 'modeling strategy'/subclass of 'Reaction' to use.
            (Defaults to 'ElementalStep')

        Returns
        -------
        list
            List of tuples to unpack for the init method of the Reaction
            ( or subclass), i.e. for ElementalStep is:
            (compounds, reactants, RType, AE).

        """
        if self.IsUpdated:
            Items = [item(self) for item in ReactionTypes[self.RType][Model]]
        else:
            Items = []
        return Items

class Reaction(object):
    """
    Object that represents a reaction as an elemental step.

    Parameters
    ----------
    reactants : list
        `Compound` instances that represent the reactants of the reaction
    products : list
        `Compound` instances that represent the products of the reaction
    TS : TransitionState
        Transition state that connects reactants and products. It must contain
        a 'energy' attribute

    Attributes
    ----------
    reactants
    compounds
    key : int
        number that identifies the reaction in a `ChemicalSystem` instance.
    products
    TS
    T

    """

    def __init__(self,reactants,products,TS=None,T=298.15):
        self.reactants = Counter(reactants)
        self.products = Counter(products)
        self.compounds = {i:0 for i in chain(reactants,products)}
        for i in self.reactants.elements():
            self.compounds[i] -= 1
        for i in self.products.elements():
            self.compounds[i] += 1
        self.key = None
        self.T = T # Temperature in K
        self.TS = TS # Transition state associated, may be shared
    def __str__(self):
        reactants = ' + '.join([r.label for r in self.reactants.elements()])
        products = ' + '.join([p.label for p in self.products.elements()])
        out = '{}    {}    {}'.format(reactants,'=>',products)
        return out
    def __repr__(self):
        n = len(self.reactants)
        m = len(self.products)
        i = self.key
        cls = type(self).__name__
        Rformat = '< {} {} of Type {} with {} reacts and {} products >'.format
        return Rformat(cls,i,'=>',n,m)
    def __eq__(self,other):
        return str(self) == str(other)
    def __contains__(self,item):
        return item in self.compounds

    def coefficient(self,compound):
        """
        Returns the global coefficient of the compound in the reaction if
        the compound is in the reaction, otherwise returns a 0
        This function acts as a wrapper for readability
        """
        return self.compounds.get(compound,0)

    @property
    def reactants_energy(self):
        return sum([c.energy for c in self.reactants.elements()])

    @property
    def AE(self):
        AE = self.TS.energy - self.reactants_energy
        return  AE
    @AE.setter
    def AE(self,other):
        reactants_energy = sum([c.energy for c in self.reactants.elements()])
        try:
            self.TS.energy = reactants_energy + other
        except AttributeError as e:
            raise RuntimeError('Attempting to set the activation energy to '\
                                 'a TSless reaction.')

    @property
    def k(self):
        try:
            k = self.ActE2k(self.T,self.AE)
        except ValueError as e:
            msg = str(e) + ' in reaction:\n\t{}'.format(str(self))
            e.args = (msg,)
            raise e
        return k

    @staticmethod
    def ActE2k(T,AE):
        """
        Transforms an activation energy to a kinetic constant. Uses TST and
        Erying Equation with the transmission coefficient equal to 1. The
        partition coefficients correspond to fluid phase.

        Parameters
        ----------
        T : float
            Temperature in Kelvin.
        AE : Energy
            Activation energy.

        Returns
        -------
        float

        Raises
        -------
        ValueError
            If a negative activation energy is provided

        """
        kb = 1.38064852E-023
        h =  6.62607004E-034
        R =  8.31445985
        if AE < 0:
            raise ValueError('Negative Activation energy found')
        AE_SI  = AE.to_SI()
        k = (kb*(T)/h)*numpy.exp(-AE_SI/(R*(T)))
        return k

class ChemicalSystem(object):
    """
    This class represents a chemical system, understanding
    it as a set of compounds, reactions, barriers and Free Energies
    of Reaction, it includes the utilities for human representation
    (A + B => C) and the matrix representation.

    Parameters
    ----------
    T : float
        Temperature in K (the default is 298.15).

    Attributes
    ----------
    shape : tuple
        (n_reactions, n_compounds)
    Name2Compound : dict
        Mapping used for accesing the compounds by name
    compounds : list
    reactions : list
    sreactions : list
        List of human-ish representations of the reactions. Rules and types
        are established in InputParse.py
    T

    """
    def __init__(self,T=298.15):
        self._T = T # K
        self.compounds = []
        self.reactions = []
        self.transitionstates = []
        self.Name2Compound = dict()
        self.Name2TS = dict()

    def __repr__(self):
        n,m = self._shape()
        cls = type(self)
        return '<{} with {} compounds and {} reactions>'.format(cls,n,m)
    def __contains__(self,item):
        return item in self.compounds or item in self.reactions

    @property
    def shape(self): 
        return len(self.reactions),len(self.compounds)
    
    @property
    def species(self):
        return len(self.compounds)
    
    @property
    def T(self):
        return self._T
    @T.setter
    def T(self,other):
        self._T = other
        for reaction in self.reactions:
            reaction.T = T


    # Methods related with the (c)ompounds attribute
    def cadd(self,compound,update=False):
        """
        Forcefully adds a compound object and updates the dictionaries
        overwriting it if it already existed.

        Parameters
        ----------
        compound : Compound
        update : bool
            Whether to update the compounds in the system after adding it or not
            (the default is False).
        """
        if compound.label in self.Name2Compound:
            index = self.Name2Compound[compound.label].key
            old = self.compounds[index]
            self.compounds[index] = compound
            raise Warning(f'Overwrote {old} at {self} with {compound}')
        else:
            self.compounds.append(compound)
        self.Name2Compound[compound.label] = compound
        if update:
            self.cupdate()
    def cextend(self,compounds):
        """
        Extends the self.compounds attr without overwriting.

        Parameters
        ----------
        new_compounds : list
            list of `Compound` instances

        """
        for compound in compounds:
            if compound not in self:
                self.cadd(compound)
        self.cupdate()
    def cupdate(self,start=0,UpdateDict=False):
        """
        Updates the `compounds` attr reinitializing the indices.

        Parameters
        ----------
        start : int
            Positive integrer that represents the initial index used to
            represent compounds. This parameter should only be changed when
            writing programs for non-python languages or for including the T
            and/or V in the differential equations.(the default is 0)
        UpdateDict : bool
            if True it will update the `Name2Compound` attr
            (the default is False).
        Raises
        -------
        AssertionError
            If the start value is negative the whole
        """
        assert start >= 0
        compounds = self.compounds
        for i,C in enumerate(compounds):
            C.key = i + start
            if UpdateDict:
                self.Name2Compound[C.label] = C

    def check_compounds(self):
        """
        Checks if all the compounds in the sreactions are 'defined' in the
        within the ChemicalSystem instance.

        Returns
        -------
        set
            set with the labels of the compounds not 'defined'
        """
        missing_compounds = set()
        for reaction in self.reactions:
            for compound in reaction.compounds:
                if compound not in self.Name2Compound: # Is this if overhead ?
                     missing_compounds.add(compound)
        return missing_compounds

    # Methods related with the (r)eactions attribute
    def radd(self,reaction,update=True):
        """
        Adds and updates a Reaction Instance avoiding duplicates.

        Parameters
        ----------
        reaction : Reaction
            reaction to add to the system.
        update : bool
            if True it modifies the state of the Reaction Instance
            (the default is True).
        """
        if reaction not in self:
            if update:
                reaction.key = self.reactions[-1].key + 1
                reaction.T = self.T
            self.reactions.append(reaction)
    def rextend(self,reactions):
        """
        Lazy Addition, it first adds all the reactions and then
        it updates them

        Parameters
        ----------
        reactions : list
            list of Reaction instances
        """
        for R in reactions:
            self.radd(R,update=False)
        self.rupdate()
    def rupdate(self,start=0):
        """
        Updates the `reactions` attr updating the key, T and k of each
        Reaction Instance.

        Parameters
        ----------
        start : int
            Initial index (the default is 0).
        """
        reactions = self.reactions
        for i,R in enumerate(reactions):
            R.key = i + start
            R.T = self.T

    def reactions_fromFile(self,FilePath,EParse,Model='ElementalStep'):
        """
        Reads the reactions from the file and enforces that the
        Energy Parsing mode is either 'relative' or 'absolute'.

        Parameters
        ----------
        FilePath : str
            A valid to a file that follows the specifications of the `README`.
        EParse : str
            How the energy provided with the reaction should be considered
            either 'relative' to reactants or the TS 'absolute' energy
        """
        assert EParse.lower() == 'absolute' or EParse.lower() == 'relative'
        sreactions = []
        with open(FilePath,'r') as F:
            for line in F:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                sreactions.append(SymbReaction(line,EParse))
        self.srextend(sreactions)

        reactions = []
        # Use self.sreactions to avoid the presence of duplicates in sreactions
        for SR in self.sreactions:
            new_reactions = ElementalStep.from_sreaction(SR,Model)
            for R in new_reactions:
                reactions.append(R)
                SR.child(R)
        self.rextend(reactions)

    # Utilities
    def RxCMatrix(self):
        """
        Returns a matrix that represents the system in which per each
        Matrix[i,j] corresponds to the coefficient of the j compound in the
        i reaction. Possitive if its a product and negative if its a reactant.

        Returns
        -------
        numpy.array

        Raises
        -------
        OperationError
            If the ChemicalSystem is not Updated
        """
        RxCMat = numpy.zeros(shape=self.shape)
        for i,reaction in enumerate(self.reactions):
            for R in reaction.reactants.elements(): # For each Reactant
                RxCMat[i,R.key] -= 1.0
            for P in reaction.products.elements(): # for each Product
                RxCMat[i,P.key] += 1.0
        return RxCMat

    # Main API
    def massbalance(self,compound):
        """
        returns the mass balance of the compound for the chemical system.

        Parameters
        ----------
        compound : Compound
            Compound whose mass balance is to be obtained

        Returns
        -------
        MB
            MassBalance object
        """
        MB = MassBalance(compound)
        n = len(self.reactions)
        for reaction in self.reactions:
            coef = reaction.coefficient(compound)
            # compound in reaction and does not appear as reactant and product
            if abs(coef) > 0:
                MB.items.append((coef,reaction))
        return MB
    def jacobian(self):
        """
        returns a list of Jacobian elements for the mass balances of all the 
        compounds in the system.

        Returns
        -------
        Out
            list of JacobianElement items.
        """
        Out = []
        for compound in self.compounds:
            MB = self.massbalance(compound)
            for compound2 in self.compounds:
                Out.append(MB.partial(compound2))
        return Out
