"""
This module contains the rules used for parsing the reactions.txt File.
This module is intended to be imported only by the 'Classes' module.
The module builds the 'ReactionTypes' variable from two variables:
'items' and 'keys'.

The 'keys' variable contains the string identifiers of the behaviours encoded in
items:

   *  The first position is the string corresponding to the marker used in the
      input file, such as '<=>' or '<d>'. It is not advised to use empty spaces
      or \\t in the markers.
   *  The second and third elements are functions that take a tuple of
      reactants, a tuple of products and an energy associated to the reaction
      and should return the energy of the transition state.
   *  The second element corresponds to the behaviour when the energies provided
      by the user are absolute energies.
   *  The third element corresponds to the behaviour when the energies provided
      by the user are relative energies. (By design, in this case the values
      associated to the reaction correspond to the barrier).
   *  Each of the following elements correspond to the behaviours of each model.
      The corresponding key (in keys) is a string with the name of the model.
      Every one of these items has to be a tuple of functions. These functions
      have to receive as input a SymbReaction instance and return a tuple with
      the parameters needed to create an instance of the corresponding
      'Reaction' subclass.
"""

from collections import OrderedDict

__all__ = ['ReactionTypes','EnergyInput']

####################### Rules for parsing energies #############################

def EnergyInput(item): # TODO
    if item.startswith('file'):
        pass
    elif any(item.endswith(i) for i in EnergyUnits):
        pass
    else:
        return float(item)

############################ Auxiliary Functions ###############################

# Auxiliary functions that define how to calculate the energy of the TS from   #
# the input:                                                                   #

## When the Energy is given as 'absolute' values
absolute = lambda r,p,e: e
diffusion = lambda r,p,e: e + max(r,p)

## When the Energy is given as 'relative' values
from_reacts = lambda r,p,e: e + r
from_products = lambda r,p,e: e + p # included for completeness

# Auxiliary functions that return input as needed to instantiate a certain     #
# Reaction model.                                                              #

## Model: ElementalStep (default)
def RtoP(sreaction):
    Rs,Ps = sreaction.reactants, sreaction.products
    AE = sreaction.TS - sum(R.energy for R in Rs)
    return Rs, Ps, AE
def PtoR(sreaction):
    Rs,Ps = sreaction.reactants, sreaction.products
    AE = sreaction.TS - sum(P.energy for P in Ps)
    return Ps, Rs, AE

######################## 'keys' and 'items' variables ##########################

keys  =  [ None,'absolute', 'relative','ElementalStep']
items = [['<=>', absolute, from_reacts,(RtoP,PtoR)],
         [ '=>', absolute, from_reacts,(RtoP,)],
         ['<=' , absolute, from_reacts,(PtoR,)], #By design: Erel -> from_reacts
         ['<d>', diffusion,  diffusion,(RtoP,PtoR)],
         ['<s>', absolute, from_reacts,(RtoP,PtoR)]]

######################### High Level API preparation ###########################

# Currently it should not be an OrderedDict since the current implementation
# should not have any clash between '<=>','=>' and '<='. However Better safe
# than sorry as the only cost is the usage of an OrderedDict that comes with the
# standard python library

ReactionTypes = OrderedDict()
for item in items:
    ReactionTypes[item[0]] = {key:item2 for key,item2 in zip(keys[1:],item[1:])}
