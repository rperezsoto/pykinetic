import re 

from .Classes import *

def chemicalsystem_fromfiles(cls,file_c,file_r,energy_unit='J/mol',relativeE=False):
    chemicalsystem = cls()

    # Process the Compound File
    raw_compounds = read_compounds(file_c)
    compounds = create_compounds(raw_compounds)
    chemicalsystem.cextend(compounds)

    # Process the reaction File
    raw_reactions, TS_lines = read_reactions(file_r)
    if TS_lines: 
        TS_dict = create_TS_dict(TS_lines)
    else:
        TS_dict = None
    for reactants,mark,products,TS_text in raw_reactions: 
        reactants = [chemicalsystem.Name2Compound[r.label] for r in reactants]
        products  = [chemicalsystem.Name2Compound[p.label] for p in products]
        reactants_energy = sum([r.energy for r in reactants])
        if mark == '<=':
            reactants, products = products, reactants
        label,energy,scannable = prepare_inline_TS(TS_text,mark,TS_dict)
        # Create the reaction/s
        reactions = [Reaction(reactants,products),]
        if mark in ['<=>','<d>']:
            reactions.append(Reaction(products,reactants))
        # Handle relative energy specification
        if mark != '<=>' and relativeE: 
            energy = energy + reactants_energy
        # Create the TS
        ts_cls = markers[mark]['TS']
        TS = ts_cls(energy,reactions=reactions,label=label,
                    scannable=scannable)
        # Add the reactions to the chemical system
        for reaction in reactions: 
            reaction.TS = TS
            chemicalsystem.radd(reaction)
    return chemicalsystem


def read_compounds(file):
    raw_compounds = []
    with open(FilePath,"r") as F:
        for line in F:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            raw_compounds.append(line.split())
    return raw_compounds
def create_compounds(raw_compounds,energy_unit='J/mol'):
    err_msg = "Unexpected number of items in \n '{}' "
    compounds = []
    for compound in raw_compounds:
        if len(compound) == 3:
            label,energy,unit = compound
        elif len(compound) == 2: 
            label,energy = compound
            unit = energy_unit
        else:
            raise RuntimeError(err_msg.format('  '.join(compound)))
        compounds.append(Compound(label,Energy(energy,unit)))
    
    try: # Check if there is any duplicates
        assert len(set(compounds)) == len(compounds)
    except AssertionError as e: 
        msg = 'Inconsistent number of compounds. Check for duplicates'
        raise ValueError(msg)
    return compounds

def read_reactions(file):
    reaction_lines = []
    TS_lines = [] 
    with open(FilePath,'r') as F:
        for line in F:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            elif '!' in line:                  # [A-]  +  [B+] => C  !TS1 
                reaction_lines.append(line)
            elif len(line.split()) in [2,3]:   # TSNAME    Number   (Unit)
                TS_lines.append(line)
    raw_reactions = [split_reaction_line(line) for line in reaction_lines]
    return raw_reactions,TS_lines
def create_TS_dict(TS_lines,energy_unit='J/mol'):
    TS_dict = dict()
    for line in TS_lines:
        items = line.strip().split() 
        if len(items) == 3: 
            label,energy,unit = items
        elif len(items) == 2: 
            label, energy = items
            unit = energy_unit
        else:
            raise RuntimeError(f"Unexpected number of items in \n '{line}' ")
        TS_dict[label] = Energy(energy,unit)
    
    try: # Check if there is any duplicates
        assert len(TS_dict) == len(TS_lines)
    except AssertionError as e: 
        raise ValueError("Inconsistent number of TSs at the end of the file. " \
                         "Check for duplicates")
    return TS_dict
def prepare_inline_TS(TS_text,mark,TS_dict):
    scannable = False
    # Prepare the TS data
    TS_text = TS_text.strip().split()
    if len(TS_text) == 1 and not is_energy(TS_text):
        label = TS_text
        energy = TS_dict[label]
    elif len(TS_text) == 1: 
        label = None
        energy = Energy(TS_text,energy_unit)
    elif len(TS_text) == 2:
        label = None
        number, unit = TS_text
        energy = Energy(number,unit)
    elif len(TS_text) == 3 and TS_text == 'scan':
        label = None
        number, unit = TS_text
        energy = Energy(number,unit)
        scannable = True
    else: 
        ts = ' '.join(TS_text)
        raise ValueError(f'Wrong TS definition at "{ts}"')
    return label,energy,scannable

is_energy = re.compile('([0-9]*\.[0-9]*)\s*?([^\s]*?)').findall

markers = {'<=>': {'TS':TransitionState,},  # Reversible
           '=>' : {'TS':TransitionState,},  # Direct reaction
           '<=' : {'TS':TransitionState,},  # Backwards reaction
           '<d>': {'TS':DiffusionTS,}}      # Reversible diffusion

for mark in markers: 
    markers['matcher'] = re.compile(f'(.*)\s{mark}\s(.*)\s!(.*)')

def split_reaction_line(reaction_text): 
    for mark in markers: 
        match = markers[mark]['matcher'].match(reaction_text)
        if match:
            break
    else:
        raise ValueError('the reaction does not contain a known marker')

    reactants_text, products_text, TS_text = match
    reactants = [r.strip() for r in reactants_text.split(' + ')]
    products  = [p.strip() for p in products_text.split(' + ')]
    
    return  reactants, mark, products, TS_text
    
