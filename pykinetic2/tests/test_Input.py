from pykinetic2.Classes import Compound, Energy, ChemicalSystem, Reaction, TransitionState
from pykinetic2.InputParse import *
import unittest

class InputTest(unittest.TestCase):
    def test_chemicalsystem_fromfiles(self):
        pass
    def test_read_compounds(self):
        pass
    def test_create_compounds(self):
        pass
    def test_read_reactions(self):
        pass
    def test_create_TS_dict(self):
        pass
    def test_prepare_inline_TS(self): 
        default_unit = 'kcal/mol'
        raw_data = ['TS1',
                    'TS2',
                    '25.0',
                    '-10.0',
                    '23.0 kcal/mol',
                    '0.01 hartree',
                    '2 kcal/mol scan']
        TSdict = {'TS1':(Energy(10.0,'kcal/mol'),True),
                  'TS2':(Energy(1.0,'eV'),False)}
        solutions = [('TS1',TSdict['TS1'][0],TSdict['TS1'][1]),
                     ('TS2',TSdict['TS2'][0],TSdict['TS2'][1]),
                     (None,Energy(25.0,default_unit),False),
                     (None,Energy(-10.0,default_unit),False),
                     (None,Energy(23.0,'kcal/mol'),False),
                     (None,Energy(0.01,'hartree'),False),
                     (None,Energy(2,'kcal/mol'),True)]
        for arg,sol in zip(raw_data,solutions):
            with self.subTest(arg=arg): 
                test = prepare_inline_TS(arg,'fakemark',TSdict,default_unit)
                self.assertEqual(test,sol)
    def test_split_reaction_line(self):
        reactions = """
        A  +  B    <=>   C    +   D    !25.0 kcal/mol 
        A  +  B    <=>   C    +   D    !25.0 
        A  +  B    <=>   C    +   D    !TS1
        A          <=>   C             !TS1
        A  +  B    <=>  [C+]  +   D    !TS1
        A  +  B    =>    C    +   D    !TS1
        C  +  D    <=    A    +   B    !TS1
        A  +  B    <d>   C    +   D    !25.0 
        A        <ERROR>      C        !TSError
        """.strip().split('\n')
        solutions = [(['A','B'],'<=>',['C','D'],'25.0 kcal/mol'),
                     (['A','B'],'<=>',['C','D'],'25.0'),
                     (['A','B'],'<=>',['C','D'],'TS1'),
                     (['A',],'<=>',['C',],'TS1'),
                     (['A','B'],'<=>',['[C+]','D'],'TS1'),
                     (['A','B'],'=>',['C','D'],'TS1'),
                     (['C','D'],'<=',['A','B'],'TS1'),
                     (['A','B'],'<d>',['C','D'],'25.0')]
        for reaction,solution in zip(reactions,solutions): 
            # The last reaction is not iterated here by design
            with self.subTest(reaction=reaction):
                test = split_reaction_line(reaction)
                self.assertEqual(test,solution)
        with self.subTest(test='error'):
            with self.assertRaises(ValueError):
                _ = split_reaction_line(reactions[-1])

    
    