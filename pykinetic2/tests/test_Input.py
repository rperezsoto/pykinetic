from pykinetic2.Classes import Compound, Energy, ChemicalSystem, Reaction, TransitionState
from pykinetic2.InputParse import *
import unittest
from unittest.mock import mock_open, patch

class InputTest(unittest.TestCase):
    def test_chemicalsystem_fromfiles(self):
        pass
    def test_read_compounds(self):
        module_name = 'pykinetic2.InputParse'
        raw_compounds = """
        A   20      kcal/mol
        B   10      hartree
        C   1   
        D   0.5     eV
        """
        solutions = [['A','20','kcal/mol'],
                     ['B','10','hartree'],
                     ['C','1'],
                     ['D','0.5','eV']]
        with patch(f'{module_name}.open', mock_open(read_data=raw_compounds), 
                    create=True) as m:
            tests = read_compounds('mock_file')
        for i,(test,solution) in enumerate(zip(tests,solutions)):
            with self.subTest(compound='reactions'):
                self.assertEqual(test,solution)
    def test_create_compounds(self):
        default_unit = 'kcal/mol'
        compounds = [['A','20','kcal/mol'],
                     ['B','10','hartree'],
                     ['C','1'],
                     ['D','0.5','eV']]
        solutions = [Compound('A',Energy('20','kcal/mol')),
                     Compound('B',Energy('10','hartree')),
                     Compound('C',Energy('1',default_unit)),
                     Compound('D',Energy('0.5','eV'))]
        tests = create_compounds(compounds,default_unit)
        for test,solution in zip(tests,solutions):
            with self.subTest(compound=solution.label): 
                self.assertEqual(test,solution)
                self.assertEqual(test.energy,solution.energy)
    def test_read_reactions(self):
        module_name = 'pykinetic2.InputParse'
        reactions = """
                    A  +  B    <=>   C    +   D    !25.0 kcal/mol 
                    A  +  B    <=>   C    +   D    !25.0 
                    A  +  B    <=>   C    +   D    !TS1 
                    A          <=>   C             !TS1 
                    TS1     25
                    TS2     72      J/mol
                    TS3     1.2     eV      scan
                    """
        solutions = [[(['A','B'],'<=>',['C','D'],'25.0 kcal/mol'), 
                      (['A','B'],'<=>',['C','D'],'25.0'),
                      (['A','B'],'<=>',['C','D'],'TS1'), 
                      (['A',   ],'<=>',['C',   ],'TS1')],
                      ['TS1     25',
                       'TS2     72      J/mol',
                       'TS3     1.2     eV      scan']]
        with patch(f'{module_name}.open', mock_open(read_data=reactions), 
                    create=True) as m:
            test = read_reactions('mock_file')
        with self.subTest(block='reactions'):
            self.assertEqual(test[0],solutions[0])
        with self.subTest(block='TSs'):
            self.assertEqual(test[1],solutions[1])
    def test_create_TS_dict(self):
        default_unit = 'kcal/mol'
        lines = ['TS1    -22 ',
                 'TS2    1.2    eV',
                 'TS3    25.0',
                 'TS4    -10.0  scan',
                 'TS5    23.0   kcal/mol',
                 'TS6    0.01   hartree',
                 'TS7    2  kcal/mol    scan']
        solution = {'TS1':(Energy(-22,'kcal/mol'),False),
                    'TS2':(Energy(1.2,'eV'),False),
                    'TS3':(Energy(25.0,'kcal/mol'),False),
                    'TS4':(Energy(-10.0,'kcal/mol'),True),
                    'TS5':(Energy(23.0,'kcal/mol'),False),
                    'TS6':(Energy(0.01,'hartree'),False),
                    'TS7':(Energy(2,'kcal/mol'),True)}
        TS_dict = create_TS_dict(lines,default_unit)
        for key,item in TS_dict.items():
            with self.subTest(key=key):
                self.assertEqual(item,solution.get(key,None))
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

    
    