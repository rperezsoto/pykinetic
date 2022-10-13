"""
This module contains the Writer classes that are in charge of translating the 
provided chemical system to the different languages or formats as well as the 
Base class to inherit from, 'Writer'. 
Currently only python and c++ are supported. This module depends on the 
template files that come with the libary.
"""

from collections import Counter
from ._base import TEMPLATES_PATH, Indent, Writer

class Batch(Writer):

    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        for var in ['x', 'dxdt', 'f', 'jac','jac_f']:
            self.parameters[var] = getattr(self,var)
        self.parameters['method'] = 'LSODA'

    def _load_default_header(self):
        with open(TEMPLATES_PATH / 'python_batch.head','r') as F:
            txt = F.read()
        self._header = txt
    def _load_default_tail(self):
        with open(TEMPLATES_PATH / 'python_batch.tail','r') as F:
            txt = F.read()
        self._tail = txt

    def set_parameters(self,simulation=None,convergence=None):
        if simulation is not None:
            for key in ['tfin','trep','dt']: 
                self.parameters[key] = simulation.get(key,'')
            concentrations = self._initial_concentrations(simulation)
            self.parameters['concentrations'] = concentrations
        if convergence is not None:
            self.parameters['convergence'] = convergence.as_str(sep=',')

    # methods for object -> str transformations
    def ratelaw(self,reaction): # Currently only for Elemental Steps
        var = f'r{reaction.key:02.0f}'
        Aux = []
        expr = f'k{reaction.key:02.0f}*'
        reactants = reaction.reactants.elements() # Unpack the reactants Counter
        for R in sorted(reactants,key=lambda x: x.key):
            Aux.append(f'{self.x}[{R.key}]')
        expr += '*'.join(Aux)
        return  var, expr
    def ratelaw_partial(self,reaction,compound): # Currently only for Elemental Steps
        #var = f'r{reaction.key:02.0f}' Legacy code, to remove
        if compound not in reaction.reactants:
            return ''
        Aux = []
        # Create a new Unpack the reactants Counter
        reactants = Counter(reaction.reactants.elements())
        coef = reactants[compound]
        if coef > 1:
            Aux.append(f'{coef:0.1f}*k{reaction.key:02.0f}')
        else:
            Aux.append(f'k{reaction.key:02.0f}')
        reactants[compound] -= 1
        reactants = reactants.elements()
        for R in sorted(reactants,key=lambda x: x.key):
            Aux.append(f'{self.x}[{R.key}]')
        expr = '*'.join(Aux)
        return  expr
    def constant(self,reaction,value_format='{:0.10e}'):  # Currently only for Elemental Steps
        return super().constant(reaction,value_format)
    def massbalance(self,MassBalance):
        if not MassBalance.items:
            return '', '0'
        var = f'{self.dxdt}[{MassBalance.compound.key}]'
        expr = []
        for coef,reaction in MassBalance.items:
            if coef == 1:
                Aux = '+'
            elif coef == -1:
                Aux = '-'
            else:
                Aux = f'{coef:+0.1f}*'
            r_expr = f'{Aux}r{reaction.key:02.0f}'
            expr.append(r_expr)
        return var, ''.join(expr)
    def jacobian_element(self,Jac_ij):
        C1,C2 = Jac_ij.compound1, Jac_ij.compound2
        if not Jac_ij.items:
            return '', '0'
        var = f'{self.jac}[{C1.key},{C2.key}]'
        expr = []
        for coef,reaction in Jac_ij.items:
            if coef == 1:
                Aux = '+'
            elif coef == -1:
                Aux = '-'
            else:
                Aux = f'{coef:+0.1f}*'
            r_expr = self.ratelaw_partial(reaction,C2)
            if not r_expr:
                continue
            expr.append(Aux + r_expr)
        if not expr:
            return var ,'0'
        return var, ''.join(expr)

    # auxiliary writing methods
    def _kinetic_constants(self,chemicalsys):
        constants = ['#Constants',]
        for reaction in chemicalsys.reactions:
            # vars,values = # In case we consider supporting reversible elemental steps?
            var,value = self.constant(reaction)
            constants.append(f'{var} = {value}')
        return constants
    def _ratelaws(self,chemicalsys):
        ratelaws = ['#Ratelaws',]
        for reaction in chemicalsys.reactions:
            var,law = self.ratelaw(reaction)
            ratelaws.append(f'{var} = {law}')
        return ratelaws
    def _massbalances(self,chemicalsys):
        MBs = ['#MassBalances',]
        for compound in chemicalsys.compounds:
            MB = chemicalsys.massbalance(compound)
            var,law = self.massbalance(MB)
            if var:
                MBs.append(f'{var} = {law}')
        return MBs
    def _jacobian_elements(self,chemicalsys):
        Jac = ['#Non-zero Elements',]
        jac_elements = chemicalsys.jacobian()
        for jac_ij in jac_elements: 
            var,expr = self.jacobian_element(jac_ij)
            if expr and expr != '0':
                Jac.append(f'{var} = {expr}')
        return Jac
    def _function(self,chemicalsys,level=0):
        """
        Generates the code of function 'f' in  dxdt = f(x,t) where x, dxdt and t
        are vectors. This function corresponds to the system of diferential
        equations for the mass balances of the system.
        """
        definition = f'def {self.f}(t,{self.x}):'

        lines = [ f'{self.dxdt} = np.zeros({chemicalsys.species})',]
        # Write the constants block
        lines.append('')
        constants = self._kinetic_constants(chemicalsys)
        lines.extend(constants)
        # Write the ratelaws block
        lines.append('')
        ratelaws = self._ratelaws(chemicalsys)
        lines.extend(ratelaws)
        # Write the MassBalances block
        lines.append('')
        massbalances = self._massbalances(chemicalsys)
        lines.extend(massbalances)
        # Function end
        lines.append('')
        lines.append(f'return {self.dxdt}')
        lines.append('\n')

        # Add one indentation level
        lines = Indent(lines,tab='    ',level=1)

        # Add function definition
        lines.insert(0,definition)

        # Add the indentation level of the function
        lines = Indent(lines,tab='    ',level=level)

        return '\n'.join(lines)
    def _jacobian(self,chemicalsys,level=0):
        """
        Generates the code of function 'Jac' in  dxdt = Jac(x,t) where x, dxdt
        and t are vectors. This function corresponds to the Jacobian of the
        system of differential equations for the mass balances of the system.
        """
        definition = f'def {self.jac_f}(t,{self.x}):'
        n = chemicalsys.species
        lines = [ f'{self.jac} = np.zeros(shape=({n},{n}))',]
        # Write the constants block
        lines.append('')
        constants = self._kinetic_constants(chemicalsys)
        lines.extend(constants)
        # Write the jacobian elements
        lines.append('')
        elements = self._jacobian_elements(chemicalsys)
        lines.extend(elements)
        # Function end
        lines.append('')
        lines.append(f'return {self.jac}')
        lines.append('\n')

        # Add one indentation level
        lines = Indent(lines,tab='    ',level=1)

        # Add function definition
        lines.insert(0,definition)

        # Add the indentation level of the function
        lines = Indent(lines,tab='    ',level=level)

        return '\n'.join(lines)
    def _initial_concentrations(self,simulation):
        concentrations = []
        for key,val in simulation['concentrations'].items():
            concentrations.append(f'xini[{key}] = {val}')
        return '\n'.join(Indent(concentrations,level=0))

    # main writing methods
    def fill_header(self,chemicalsys):
        out_filename = self.parameters.get('out_filename','data.txt')
        self.parameters['out_filename'] = out_filename
        self.parameters['species'] = chemicalsys.species
        self.parameters['T'] = chemicalsys.T
        super().fill_header(chemicalsys)
    def write(self,chemicalsys,filepath):
        self.fill(chemicalsys)
        # Write the constants block
        #lines.append('#Constants')
        items = [self.header,
                 self.function,
                 self.jacobian,
                 self.tail]
        with open(filepath,'w') as F:
            F.write(''.join(items))

class SemiBatch(Batch):

    def __init__(self,flow,Vini,**kwargs):
        super().__init__(**kwargs)
        self.keys = [] # Keys of the different reactants added
        self.flow = flow
        self.Vini = Vini # initial volume
        self.Cadd = None
    def _load_default_header(self):
        with open(TEMPLATES_PATH / 'python_semibatch.head','r') as F:
            txt = F.read()
        self._header = txt
    def _load_default_tail(self):
        with open(TEMPLATES_PATH / 'python_semibatch.tail','r') as F:
            txt = F.read()
        self._tail = txt
    def set_parameters(self,**kwargs):
        super().set_parameters(**kwargs)

    # methods for object -> str transformations
    def massbalance(self,MassBalance):
        if MassBalance.compound.key not in self.keys:
            return super().massbalance(MassBalance)
        
        var, expr = super().massbalance(MassBalance)
        expr = f'tflow*Cin{MassBalance.compound.key:02.0f} {expr}'
        return var, expr

    # main writing methods
    def fill_header(self,chemicalsys):
        out_filename = self.parameters.get('out_filename','data.txt')
        self.parameters['out_filename'] = out_filename
        self.parameters['species'] = chemicalsys.species
        self.parameters['T'] = chemicalsys.T
        self.parameters['Vini'] = self.Vini
        self.parameters['Cadd'] = self.Cadd
        self.parameters['flow'] = self.flow
        super().fill_header(chemicalsys)

    # auxiliary writing methods
    def _function(self,chemicalsys,level=0):
        """
        Generates the code of function 'f' in  dxdt = f(x,t) where x, dxdt and t
        are vectors. This function corresponds to the system of diferential
        equations for the mass balances of the system.
        """
        definition = f'def {self.f}(t,{self.x}):'

        lines = [ f'{self.dxdt} = np.zeros({chemicalsys.species}+1)',]
        # Write the Inlet Concentrations and Flow
        lines.append(f'tflow = {self.flow}/{self.x}[-1] # flow/V; 1/s')
        inletC = self._flask_initial_concentrations()
        lines.extend(inletC)
        # Write the constants block
        lines.append('')
        constants = self._kinetic_constants(chemicalsys)
        lines.extend(constants)
        # Write the ratelaws block
        lines.append('')
        ratelaws = self._ratelaws(chemicalsys)
        lines.extend(ratelaws)
        # Write the MassBalances block
        lines.append('')
        massbalances = self._massbalances(chemicalsys)
        lines.extend(massbalances)
        # Include the Volume 
        lines.append('')
        lines.append(f'{self.dxdt}[-1] = tflow # flow/V; 1/s')
        lines.append(f'{self.dxdt}[:] = {self.dxdt}[:]*{self.x}[-1]')
        # Function end
        lines.append('')
        lines.append(f'return {self.dxdt}')
        lines.append('\n')

        # Add one indentation level
        lines = Indent(lines,tab='    ',level=1)

        # Add function definition
        lines.insert(0,definition)

        # Add the indentation level of the function
        lines = Indent(lines,tab='    ',level=level)

        return '\n'.join(lines)
    def _jacobian(self,chemicalsys,level=0):
        """
        Generates the code of function 'Jac' in  dxdt = Jac(x,t) where x, dxdt
        and t are vectors. This function corresponds to the Jacobian of the
        system of differential equations for the mass balances of the system.
        """
        definition = f'def {self.jac_f}(t,{self.x}):'
        n = chemicalsys.species
        lines = [ f'{self.jac} = np.zeros(shape=({n}+1,{n}+1))',]
        # Write the constants block
        lines.append('')
        constants = self._kinetic_constants(chemicalsys)
        lines.extend(constants)
        # Write the jacobian elements
        lines.append('')
        elements = self._jacobian_elements(chemicalsys)
        lines.extend(elements)
        # Include the Volume 
        lines.append('')
        lines.append(f'{self.jac}[:,:] = {self.jac}[:,:]*{self.x}[-1]')

        # Function end
        lines.append('')
        lines.append(f'return {self.jac}')
        lines.append('\n')

        # Add one indentation level
        lines = Indent(lines,tab='    ',level=1)

        # Add function definition
        lines.insert(0,definition)

        # Add the indentation level of the function
        lines = Indent(lines,tab='    ',level=level)

        return '\n'.join(lines)
    def _set_concentrations(self,concentrations):
        reactor = dict()
        flask = dict()
        for key,(r,f) in concentrations.items(): 
            if r != 0: 
                reactor[key] = r
            if f != 0: 
                flask[key] = f
        self.parameters['concentrations'] = reactor
        self.Cadd = flask
        self.keys = list(flask.keys())
    def _initial_concentrations(self,simulation):
        self._set_concentrations(simulation['concentrations'])
        return super()._initial_concentrations(self.parameters)
    def _flask_initial_concentrations(self): 
        concentrations = []
        for key,val in self.Cadd.items():
            concentrations.append(f'Cin{key:02.0f} = {val}')
        return concentrations
class SemiBatchExtended(SemiBatch):
    """
    A writer to account for a SemiBatch Reactor with an influx of species 
    coupled after a 'tsection' time of the total simulation with a Batch reactor. 
    As a consequence the Simulation parameters required for this type of model 
    require the parameters 'tsection' and 'tfin' as well as the parameters for 
    the initial SemiBatch reactor. For writing it up to two Chemical Systems 
    can be provided.
    """
    def _load_default_header(self):
        with open(TEMPLATES_PATH / 'python_semibatch_extended.head','r') as F:
            txt = F.read()
        self._header = txt
    def _load_default_tail(self):
        with open(TEMPLATES_PATH / 'python_semibatch_extended.tail','r') as F:
            txt = F.read()
        self._tail = txt
    def set_parameters(self,simulation=None,convergence=None):
        if simulation is not None:
            self.parameters.update(simulation)
            concentrations = self._initial_concentrations(simulation)
            self.parameters['concentrations'] = concentrations
        if convergence is not None:
            self.parameters['convergence'] = convergence.as_str(sep=',')

    def _function(self,systems,level=0):
        """
        Generates the code of function 'f' in  dxdt = f(x,t) where x, dxdt and t
        are vectors. This function corresponds to the system of diferential
        equations for the mass balances of the system.
        """
        _f = self.f
        self.f = f'{self.f}_semi'

        fadd = super()._function(systems[0],level)
        
        self.f = f'{_f}_batch'
        # Remove the "added compounds"
        _keys = self.keys
        self.keys = []

        fnoadd = super(SemiBatch,self)._function(systems[1],level)
        
        # Return everything to normal
        self.f = _f
        self.keys = _keys


        return fadd + fnoadd
    def _jacobian(self, systems, level=0):
        _jac_f = self.jac_f 
        self.jac_f = f'{self.jac_f}_semi'

        jacadd = super()._jacobian(systems[0],level)
        
        self.jac_f = f'{_jac_f}_batch'
        # Remove the "added compounds"
        _keys = self.keys
        self.keys = []

        jacnoadd = super(SemiBatch,self)._jacobian(systems[1],level)

        self.jac_f = _jac_f
        self.keys = _keys

        return jacadd + jacnoadd
    
    # main methods for writing
    def fill(self,system1,system2=None):
        """
        Reads the information of the chemical system and updates the values
        needed for writing.

        Parameters
        ----------
        system1 : ChemicalSystem
            Chemical System to be written for the SemiBatch reactor model section
        system2 : ChemicalSystem, optional
            Chemical System to be written for the Batch reactor model section, 
            If none is provided it defaults to the system1.
        """
        if system2 is None: 
            systems = [system1,system1]
        else: 
            systems = [system1,system2]
        self.fill_header(systems[0])
        self.fill_tail(systems[0])
        self.function = self._function(systems)
        self.jacobian = self._jacobian(systems)
class PFR(Batch):
    """
    Writer class of python scripts for ideal Plug Flow Reactors (PFR) models.
    Models using this writer take inlet concentrations as floats and not tuples.
    Allows providing the parameters through a SimulationParameters instance or 
    directly upon instantiation.
    
    Parameters
    ----------
    flow : float
        Flowrate in L/s
    volume : float
        Reactor Volume in L

    """
    def __init__(self,flow=1.0,volume=1.0,**kwargs):
        super().__init__(**kwargs)
        self.flow = flow # L/s
        self.volume = volume #L
    def _load_default_header(self):
        with open(TEMPLATES_PATH / 'python_PFR.head','r') as F:
            txt = F.read()
        self._header = txt
    def _load_default_tail(self):
        with open(TEMPLATES_PATH / 'python_PFR.tail','r') as F:
            txt = F.read()
        self._tail = txt
    def set_parameters(self,**kwargs):
        super().set_parameters(**kwargs)
        self.parameters['flow'] = self.parameters.get('flow',self.flow)
        self.parameters['volume'] = self.parameters.get('volume',self.volume)

    # main writing methods
    def fill_header(self,chemicalsys):
        out_filename = self.parameters.get('out_filename','data.txt')
        self.parameters['out_filename'] = out_filename
        self.parameters['species'] = chemicalsys.species
        self.parameters['T'] = chemicalsys.T
        super().fill_header(chemicalsys)

    def _function(self,chemicalsys,level=0):
        """
        Generates the code of function 'f' in  dxdt = f(x,t) where x, dxdt and t
        are vectors. This function corresponds to the system of diferential
        equations for the mass balances of the system.
        """
        definition = f'def {self.f}(t,{self.x}):'

        lines = [ f'{self.dxdt} = np.zeros({chemicalsys.species})',]
        # Write the Inlet Concentrations and Flow
        lines.append(f'flow = {self.flow} # flow; L/s')
        # Write the constants block
        lines.append('')
        constants = self._kinetic_constants(chemicalsys)
        lines.extend(constants)
        # Write the ratelaws block
        lines.append('')
        ratelaws = self._ratelaws(chemicalsys)
        lines.extend(ratelaws)
        # Write the MassBalances block
        lines.append('')
        massbalances = self._massbalances(chemicalsys)
        lines.extend(massbalances)
        # Include the Volume 
        lines.append('')
        lines.append(f'{self.dxdt}[:] = {self.dxdt}[:]/flow')
        # Function end
        lines.append('')
        lines.append(f'return {self.dxdt}')
        lines.append('\n')

        # Add one indentation level
        lines = Indent(lines,tab='    ',level=1)

        # Add function definition
        lines.insert(0,definition)

        # Add the indentation level of the function
        lines = Indent(lines,tab='    ',level=level)

        return '\n'.join(lines)
    def _jacobian(self,chemicalsys,level=0):
        """
        Generates the code of function 'Jac' in  dxdt = Jac(x,t) where x, dxdt
        and t are vectors. This function corresponds to the Jacobian of the
        system of differential equations for the mass balances of the system.
        """
        definition = f'def {self.jac_f}(t,{self.x}):'
        n = chemicalsys.species
        lines = [ f'{self.jac} = np.zeros(shape=({n},{n}))',]
        lines.append(f'flow = {self.flow} # flow; L/s')
        # Write the constants block
        lines.append('')
        constants = self._kinetic_constants(chemicalsys)
        lines.extend(constants)
        # Write the jacobian elements
        lines.append('')
        elements = self._jacobian_elements(chemicalsys)
        lines.extend(elements)
        # Include the Volume 
        lines.append('')
        lines.append(f'{self.jac}[:,:] = {self.jac}[:,:]/flow')

        # Function end
        lines.append('')
        lines.append(f'return {self.jac}')
        lines.append('\n')

        # Add one indentation level
        lines = Indent(lines,tab='    ',level=1)

        # Add function definition
        lines.insert(0,definition)

        # Add the indentation level of the function
        lines = Indent(lines,tab='    ',level=level)

        return '\n'.join(lines)
    
