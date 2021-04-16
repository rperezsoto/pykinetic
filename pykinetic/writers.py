"""
This module contains the Writer classes that are in charge of translating the 
provided chemical system to the different languages or formats as well as the 
Base class to inherit from, 'Writer'. 
Currently only python and c++ are supported. This module depends on the 
template files that come with the libary.
"""

from abc import abstractmethod
from pathlib import Path
from collections import Counter, defaultdict
from pkg_resources import resource_filename

TEMPLATES_PATH = Path(resource_filename('pykinetic','templates'))

def Indent(lines,tab='    ',level=1):
    """
    Takes an iterable of strings, each representing a line and returns them with
    the indentation level and character specified.
    """
    Out = []
    for line in lines:
        Out.append(f'{tab*level}{line}')
    return Out

class Writer(object):
    def __init__(self, conc_var='x', mb_var='dxdt', fun_var='model',
                 jac_var='Jac',jac_fun_var='jacobian',header=None,tail=None):
        self.x = conc_var
        self.dxdt = mb_var
        self.f = fun_var
        self.jac = jac_var
        self.jac_f = jac_fun_var
        if header is None:
            self._load_default_header()
        elif header:
            self._header = header
        else:
            self._header = ''
        if tail is None:
            self._load_default_tail()
        elif tail:
            self._tail = tail
        else:
            self._tail = ''
        self.clear()
        self.parameters = defaultdict(lambda : 'MISSING')

    # Defaults Load
    @abstractmethod
    def _load_default_tail(self):
        pass
    @abstractmethod
    def _load_default_header(self):
        pass

    # methods for parameter reading
    def set_parameters(self,simulation=None,convergence=None):
        if simulation is not None:
            self.parameters.update(simulation)
        if convergence is not None:
            self.parameters.update(convergence)

    # methods for object -> str transformations
    def constant(self,reaction,value_format):
        """
        Takes in a reaction and returns 2 expressions, the string of the 
        variable and the string with the value of the constant. I.e.
        # Reaction(1 => 2)
        A = <ElementalStep 1 of Type '=>' with reacts=(1,) products(2,)>
        'k01', '1.00000000' = PythonWriter.constant(A)

        Parameters
        ----------
        reaction : Reaction
            Any reaction

        Returns
        -------
        var,expr
            variable and expresion of the reaction
        """
        return f'k{reaction.key:02.0f}', value_format.format(reaction.k)

    @abstractmethod
    def ratelaw(self,reaction):
        """
        Takes in a reaction and returns 2 expressions, the string of the 
        variable and the string of the ratelaw expression. I.e.
        # Reaction(1 => 2)
        A = <ElementalStep 1 of Type '=>' with reacts=(1,) products(2,)>
        'r01', 'k01*x[1]' = PythonWriter.ratelaw(A)

        Parameters
        ----------
        reaction : Reaction
            Any reaction

        Returns
        -------
        var,expr
            variable and expresion of the reaction
        """
        pass
    @abstractmethod
    def ratelaw_partial(self,reaction,compound):
        """
        Takes in a reaction and returns and the string of the partial derivative
         the ratelaw with respect to the concentration of compound. I.e.
        # Reaction(1 => 2), Compound 1
        A = <ElementalStep 1 of Type '=>' with reacts=(1,) products(2,)>
        # r01,'k01*x[1]'ratelaw(A)
        'k01' = PythonWriter.ratelaw_partial(A)

        Parameters
        ----------
        reaction : Reaction
            Any reaction
        compound : Compound
            Any compound

        Returns
        -------
        expr
            variable and expresion of the reaction
        """
        pass
    @abstractmethod
    def massbalance(self,Massbalance):
        """
        Takes in a mass balance and returns 2 expressions, the string of the
        variable and the string of the expression. I.e.
        # Given a Chemical system with only Reaction('A <=> B')
        # d[A]dt = -r1 + r2 = -k01[A] + k2[B]
        # being A.key = 0 and B.key = 1
        A = <ElementalStep 1 of Type '=>' with reacts=(1,) products(2,)>
        'dxdt[0]', '-k01*x[1] +k02[B]' = PythonWriter.massbalance(A)

        Parameters
        ----------
        reaction : MassBalance
            Any MassBalance

        Returns
        -------
        var,expr
            variable and expresion of the reaction
        """
        pass
    @abstractmethod
    def jacobian_element(self,Jac_ij):
        """
        Takes the partial differential of a MB and returns the variable of the
        jacobian and the string of the expression. I.e.
        # Given a Chemical system with only Reaction('A <=> B')
        # d[A]dt = -r1 + r2 = -k01[A] + k2[B]
        # being A.key = 0 and B.key = 1
        C = <JacobianElement(Reaction(1),Compound(A))>
        'Jac[0,0]', '-k01' = PythonWriter.jacobian_element(A)

        Parameters
        ----------
        Jac_ij : JacobianElement
            Any JacobianElement object

        Returns
        -------
        var,expr
            variable and expresion of the reaction
        """
        pass

    @abstractmethod
    def _function(self,chemicalsys):
        pass
    @abstractmethod
    def _jacobian(self,chemicalsys):
        pass

    # main methods for writing
    def fill_header(self,chemicalsys):
        kwargs = dict()
        kwargs.update(self.parameters)
        for attr in ['dxdt','f','jac', 'jac_f']:
            kwargs[attr] = getattr(self,attr)
        self.header = self._header.format_map(kwargs)
    def fill_tail(self,chemicalsys):
        kwargs = dict()
        kwargs.update(self.parameters)
        for attr in ['dxdt','f','jac', 'jac_f']:
            kwargs[attr] = getattr(self,attr)
        self.tail = self._tail.format_map(kwargs)

    def fill(self,chemicalsys):
        """
        Reads the information of the chemical system and updates the values
        needed for writing.
        """
        self.fill_header(chemicalsys)
        self.fill_tail(chemicalsys)
        self.function = self._function(chemicalsys)
        self.jacobian = self._jacobian(chemicalsys)
    def clear(self):
        """
        Clears the variables related to writing.
        """
        self.header = ''
        self.tail = ''
        self.jacobian = ''
        self.function = ''

    @abstractmethod
    def write(self,chemicalsys,filepath):
        """ Writes the current chemical system into the file specified """
        pass

class PythonWriter(Writer):

    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        for var in ['x', 'dxdt', 'f', 'jac','jac_f']:
            self.parameters[var] = getattr(self,var)
        self.parameters['method'] = 'LSODA'

    def _load_default_header(self):
        with open(TEMPLATES_PATH.joinpath('python_header.default'),'r') as F:
            txt = F.read()
        self._header = txt
    def _load_default_tail(self):
        with open(TEMPLATES_PATH.joinpath('python_tail.default'),'r') as F:
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

class CplusplusWriter(Writer):

    def _load_default_header(self):
        with open(TEMPLATES_PATH.joinpath('cplusplus_header.default'),'r') as F:
            txt = F.read()
        self._header = txt
    def _load_default_tail(self):
        with open(TEMPLATES_PATH.joinpath('cplusplus_tail.default'),'r') as F:
            txt = F.read()
        self._tail = txt

    def set_parameters(self,simulation=None,convergence=None):
        if simulation is not None:
            for key in ['tfin','trep','dt']: 
                self.parameters[key] = simulation.get(key,'')
            concentrations = self._initial_concentrations(simulation)
            self.parameters['concentrations'] = concentrations
        if convergence is not None:
            variables = ','.join(convergence.variables)
            values = convergence.as_str(sep=';')
            lines = [f'double {variables};',
                     f"    {values};"]
            self.parameters['convergence'] = '\n'.join(Indent(lines,level=0))

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
        var = f'{self.jac}({C1.key},{C2.key})'
        if C1.key is None or C2.key is None: 
            return '', '0'
        elif not Jac_ij.items:
            return var, '0'
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
        constants = ['// Kinetic Constants',]
        for reaction in chemicalsys.reactions:
            # vars,values = # In case we consider supporting reversible elemental steps?
            var,value = self.constant(reaction)
            constants.append(f'const double {var} = {value};')
        return constants
    def _ratelaws(self,chemicalsys):
        vars = []
        ratelaws = ['// Ratelaws',]
        for reaction in chemicalsys.reactions:
            var,law = self.ratelaw(reaction)
            if var:
                ratelaws.append(f'{var} = {law};')
            vars.append(var)
        # Add variable declaration
        ratelaws.insert(1,f'double {",".join(vars)};')
        return ratelaws
    def _massbalances(self,chemicalsys):
        MBs = ['// MassBalances',]
        for compound in chemicalsys.compounds:
            MB = chemicalsys.massbalance(compound)
            var,law = self.massbalance(MB)
            MBs.append(f'{var} = {law};')
        return MBs
    def _jacobian_elements(self,chemicalsys):
        Jac = []
        jac_elements = chemicalsys.jacobian()
        for jac_ij in jac_elements: 
            var,expr = self.jacobian_element(jac_ij)
            if expr:
                Jac.append(f'{var} = {expr};')
        return Jac
    def _function(self,chemicalsys,level=0):
        """
        Generates the code of function 'f' in  dxdt = f(x,t) where x, dxdt and t
        are vectors. This function corresponds to the system of diferential
        equations for the mass balances of the system.
        """
        definition = f'void operator()(const state_type &{self.x}, '\
                     f'state_type &{self.dxdt}, const double t)'

        lines = []

        # Write the ratelaws block
        ratelaws = self._ratelaws(chemicalsys)
        lines.extend(ratelaws)

        # Write the MassBalances block
        massbalances = self._massbalances(chemicalsys)
        lines.extend(massbalances)
        # Function end

        # Add one indentation level
        lines = Indent(lines,tab='    ',level=1)

        # Add brackets
        lines.insert(0,'{')
        lines.append('}')
        # Add function definition
        lines.insert(0,definition)

        # indent once
        lines = Indent(lines,tab='    ',level=1)

        # Add brackets
        lines.insert(0,'{')
        lines.append('};')
        # Add structure name
        lines.insert(0,f'struct {self.f}')

        # Add the indentation level of the function
        lines = Indent(lines,tab='    ',level=level)

        return '\n'.join(lines)
    def _jacobian(self,chemicalsys,level=0):
        """
        Generates the code of the jacobian of the function 'f' in dxdt = f(x,t) 
        where x, dxdt and t are vectors.
        """
        definition = f'void operator()(const state_type &{self.x}, '\
                     f'matrix_type &{self.jac}, const double t, '\
                     f'state_type &dfdt)'

        lines = []

        # Write the jacobian elements
        elements = self._jacobian_elements(chemicalsys)
        lines.extend(elements)
        # Function end
        for i in range(chemicalsys.species): 
            lines.append(f'dfdt[{i}] = 0;')

        # Add one indentation level
        lines = Indent(lines,tab='    ',level=1)

        # Add brackets
        lines.insert(0,'{')
        lines.append('}')
        # Add function definition
        lines.insert(0,definition)

        # indent once
        lines = Indent(lines,tab='    ',level=1)

        # Add brackets
        lines.insert(0,'{')
        lines.append('};')
        # Add structure name
        lines.insert(0,f'struct {self.jac_f}')

        # Ensure the final indentation level
        lines = Indent(lines,tab='    ',level=level)

        return '\n'.join(lines)
    def _initial_concentrations(self,simulation):
        concentrations = []
        for key,val in simulation['concentrations'].items():
            concentrations.append(f'{self.x}[{key}] = {val};')
        concentrations[1:] = Indent(concentrations[1:],tab='    ',level=1)
        return '\n'.join(concentrations)

    # main writing methods
    def fill_header(self,chemicalsys):
        self.header = self._header.format_map(self.parameters)
    def fill_tail(self,chemicalsys):
        self.parameters['species'] = chemicalsys.species
        self.parameters['T'] = chemicalsys.T
        super().fill_tail(chemicalsys)
    def fill(self,chemicalsys):
        super().fill(chemicalsys)
        constants = self._kinetic_constants(chemicalsys)
        self.constants = '\n'.join(constants)
    def clear(self):
        super().clear()
        self.constants = ''

    def write(self,chemicalsys,filepath):
        self.fill(chemicalsys)
        items = [self.header,
                 self.constants,
                 '',
                 self.function,
                 '',
                 self.jacobian,
                 '',
                 self.tail]
        with open(filepath,'w') as F:
            F.write('\n'.join(items))
