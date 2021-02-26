import re 
from itertools import chain

from .Classes import Compound,Energy,Reaction,TransitionState,DiffusionTS


# Custom errors
class MissingCompounds(ValueError):
    pass
class MissingTransitionStates(ValueError):
    pass

# GLOBALS
DIFFUSION = '<d>'
BACKWARDS = '<='
MARKERS = {'<=>': {'TS':TransitionState,},  # Reversible
           '=>' : {'TS':TransitionState,},  # Direct reaction
           BACKWARDS : {'TS':TransitionState,},  # Backwards reaction
           DIFFUSION: {'TS':DiffusionTS,}}      # Reversible diffusion

REVERSIBLE = ['<=>',DIFFUSION] # reactions modeled as 2 elemental steps

for mark in MARKERS: 
    MARKERS[mark]['matcher'] = re.compile(f'(.*)\s{mark}\s(.*)\s!(.*)')

re_is_energy = re.compile('(-?[0-9]*\.+[0-9]*)\s*?([^\s]*?)')
def is_energy(text):
    test = re_is_energy.findall(text)
    if test:
        return test[0]
    else:
        return False

def populate_chemicalsystem_fromfiles(chemicalsystem,file_c,file_r,
                                      energy_unit='J/mol',relativeE=False):

    # Process the Compound File
    raw_compounds = read_compounds(file_c)
    compounds = create_compounds(raw_compounds,energy_unit)
    chemicalsystem.cextend(compounds)

    # Process the reaction File
    raw_reactions, TS_lines = read_reactions(file_r)
    if TS_lines: 
        TS_dict = create_TS_dict(TS_lines,energy_unit)
    else:
        TS_dict = dict()

    check_missing_compounds(chemicalsystem,raw_reactions)
    check_missing_TS(raw_reactions,TS_dict)

    # Parse the reactions and create the reaction and ts instances
    for reactants,mark,products,TS_text in raw_reactions: 
        reactants = [chemicalsystem.Name2Compound[r] for r in reactants]
        products  = [chemicalsystem.Name2Compound[p] for p in products]
        reactants_energy = sum([r.energy for r in reactants])
        if mark == BACKWARDS:
            reactants, products = products, reactants
        label,energy,scannable = prepare_inline_TS(TS_text,TS_dict,energy_unit)
        # Create the reaction/s
        reactions = [Reaction(reactants,products),]
        if mark in REVERSIBLE:
            reactions.append(Reaction(products,reactants))
        # Handle relative energy specification
        if mark != DIFFUSION and relativeE: 
            energy = energy + reactants_energy
        # Create the TS
        ts_cls = MARKERS[mark]['TS']
        if label not in chemicalsystem.Name2TS: 
            TS = ts_cls(energy,reactions=reactions,label=label,
                        scannable=scannable)
        else:
            TS = chemicalsystem.Name2TS[label]
            TS.reactions.extend(reactions)
        # Add the reactions to the chemical system
        for reaction in reactions: 
            reaction.TS = TS
            chemicalsystem.radd(reaction,False)
        chemicalsystem.rupdate()

def read_compounds(file):
    raw_compounds = []
    with open(file,"r") as F:
        for line in F:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            raw_compounds.append(line.split())
    return raw_compounds
def create_compounds(raw_compounds,energy_unit='J/mol'):
    err_msg = "Unexpected number of items in \n '{}' "
    compounds = []
    for items in raw_compounds:
        scannable = 'scan' in items[-1]
        if scannable:
            _ = items.pop(-1)
        if len(items) == 3:
            label,energy,unit = items
        elif len(items) == 2: 
            label,energy = items
            unit = energy_unit
        else:
            raise RuntimeError(err_msg.format('  '.join(items)))
        compound = Compound(label,Energy(energy,unit),scannable=scannable)
        compounds.append(compound)
    
    try: # Check if there is any duplicates
        assert len(set(compounds)) == len(compounds)
    except AssertionError as e: 
        msg = 'Inconsistent number of compounds. Check for duplicates'
        raise ValueError(msg)
    return compounds

def read_reactions(file):
    reaction_lines = []
    TS_lines = [] 
    with open(file,'r') as F:
        for line in F:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            elif '!' in line:                  # [A-]  +  [B+] => C  !TS1 
                reaction_lines.append(line)
            else:   # TSNAME    Number   (Unit)
                TS_lines.append(line)
    raw_reactions = [split_reaction_line(line) for line in reaction_lines]
    return raw_reactions,TS_lines
def create_TS_dict(TS_lines,energy_unit='J/mol'):
    TS_dict = dict()
    for line in TS_lines:
        scannable = False
        items = line.strip().split()
        # Guess if it is a scannable TS
        if items[-1] == 'scan':
            scannable = True
            _ = items.pop(-1)
        # now only two cases are possible
        if len(items) == 3:  # Energy unit is specified
            label,energy,unit = items
        elif len(items) == 2: # energy unit is not specified
            label, energy = items
            unit = energy_unit
        else:
            raise RuntimeError(f"Unexpected number of items in \n '{line}' ")
        TS_dict[label] = (Energy(energy,unit),scannable)
    
    try: # Check if there is any duplicates
        assert len(TS_dict) == len(TS_lines)
    except AssertionError as e: 
        raise ValueError("Inconsistent number of TSs at the end of the file. " \
                         "Check for duplicates")
    return TS_dict
def prepare_inline_TS(TS_text,TS_dict,energy_unit='J/mol'):
    scannable = False
    # Prepare the TS data
    TS_items = TS_text.split()
    if len(TS_items) == 1 and not is_energy(TS_text):
        label = TS_items[0]
        energy,scannable = TS_dict[label]
    elif len(TS_items) == 1: 
        label = None
        energy = Energy(TS_items[0],energy_unit)
    elif len(TS_items) == 2:
        label = None
        number, unit = TS_items
        energy = Energy(number,unit)
    elif len(TS_items) == 3 and TS_items[-1] == 'scan':
        label = None
        number, unit = TS_items[:2]
        energy = Energy(number,unit)
        scannable = True
    else: 
        ts = ' '.join(TS_items)
        raise ValueError(f'Wrong TS definition at "{ts}"')
    return label,energy,scannable
def split_reaction_line(reaction_text): 
    for mark in MARKERS: 
        match = MARKERS[mark]['matcher'].findall(reaction_text)
        if match:
            break
    else:
        raise ValueError('the reaction does not contain a known marker')

    reactants_text, products_text, TS_text = match[0]
    reactants = [r.strip() for r in reactants_text.split(' + ')]
    products  = [p.strip() for p in products_text.split(' + ')]
    
    return  reactants, mark, products, TS_text.strip()

def check_missing_compounds(chemicalsystem,raw_reactions):
    """
    Checks if all the compounds in the chemicalsystem are 'defined' 
    within the ChemicalSystem instance.

    Raises
    ------
    MissingCompounds

    """
    defined_compounds = set(chemicalsystem.Name2Compound)
    reactants,mark,products,TS_text = zip(*raw_reactions)
    reaction_compounds = set(list(chain(*reactants)) + list(chain(*products)))
    missing_compounds = reaction_compounds - defined_compounds
    if missing_compounds:
        msg = 'The following compounds are not defined: \n{}\n'
        items = [f'{item}\tNUMBER\tUNIT?' for item in sorted(missing_compounds)]
        raise MissingCompounds(msg.format('\n'.join(items)))
def check_missing_TS(raw_reactions,TS_dict):
    """
    Checks if all the labeled transition states of the reactions are 'defined' 
    are defined.

    Raises
    ------
    MissingTransitionStates

    """
    defined_TS = set(TS_dict)
    reactants,mark,products,TS_text = zip(*raw_reactions)
    reaction_TS = []
    for item in TS_text :
        if not is_energy(item): 
            reaction_TS.append(item.split()[0])
    missing_TS = set(reaction_TS) - defined_TS
    if missing_TS:
        msg = 'The following labeled TS are not defined: \n{}\n'
        items = [f'{item}\tNUMBER\tUNIT?\tSCAN?' 
                 for item in sorted(missing_TS)]
        raise MissingTransitionStates(msg.format('\n'.join(items)))