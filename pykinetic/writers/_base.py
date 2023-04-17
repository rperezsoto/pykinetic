"""
This module contains the Writer classes that are in charge of translating the 
provided chemical system to the different languages or formats as well as the 
Base class to inherit from, 'Writer'. 
Currently only python and c++ are supported. This module depends on the 
template files that come with the libary.
"""

from abc import abstractmethod
from pathlib import Path
from collections import defaultdict
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
    """
    Base class of a writer object. A writer object allows the customized 
    translation of a kinetic model into a script to carry out its simulation in 
    a specific language and under a specific system constraints.  
    """
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
