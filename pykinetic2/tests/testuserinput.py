from pykinetic2.classes import Compound, Energy, ChemicalSystem, Reaction, TransitionState
from pykinetic2.userinput import *
import unittest
from unittest.mock import mock_open, patch, MagicMock

class InputTest(unittest.TestCase):
    def setUp(self):
        TransitionState.defaultname = 0
        DiffusionTS.defaultname = 0
    def test_populate_chemicalsystem_fromfiles(self):
        module_name = 'pykinetic2.userinput'
        default_unit = 'kcal/mol'
        # Prepare the outputs of the mocks
        compounds = [['A','0.0',default_unit],
                     ['B','0.0',default_unit],
                     ['C','1.0'],
                     ['D','-2.0',default_unit]]
        reactions1 = [[(['A','B'],'<=>',['C',],f'10.0 {default_unit}'), 
                      (['C',],'=>',['D',],'TSa')],
                     ['TSa     20',]]
        reactions2 = [[(['A','B'],'<=>',['C',],f'10.0 {default_unit}'), 
                      (['C',],'=>',['D',],'TSa')],
                     ['TSa     18',]]
        # Prepare the solution
        chemsys = ChemicalSystem()
        A = Compound('A',Energy('0.0',default_unit))
        B = Compound('B',Energy('0.0',default_unit))
        C = Compound('C',Energy('1.0',default_unit))
        D = Compound('D',Energy('-2.0',default_unit))
        TS00 = TransitionState(Energy('10.0',default_unit),label='TS00')
        TS01 = TransitionState(Energy('10.0',default_unit),label='TS01')
        TSa = TransitionState(Energy('20.0',default_unit),label='TSa')
        r1 = Reaction((A,B),(C,),TS=TS00)
        r2 = Reaction((C,),(A,B),TS=TS00)
        r3 = Reaction((C,),(D,),TS=TSa)
        for compound in [A,B,C,D]: 
            chemsys.cadd(compound)
        for reaction in [r1,r2,r3]: 
            chemsys.radd(reaction,False)
        chemsys.rupdate()
        chemicalsystem1 = ChemicalSystem()
        chemicalsystem2 = ChemicalSystem()
        # run the function with the mocks
        with patch(f'{module_name}.read_compounds',return_value=compounds):
            with patch(f'{module_name}.read_reactions',return_value=reactions1):
                populate_chemicalsystem_fromfiles(chemicalsystem1,
                                                  'file_compounds',
                                                  'file_reactions',
                                                  default_unit,
                                                  False)
            with patch(f'{module_name}.read_reactions',return_value=reactions2):
                populate_chemicalsystem_fromfiles(chemicalsystem2,
                                                  'file_compounds',
                                                  'file_reactions',
                                                  default_unit,
                                                  False)
        attributes =['compounds','reactions','Name2Compound','transitionstates']
        for attr in attributes: 
            with self.subTest(energy_parse='absolute',attribute=attr):
                test = getattr(chemicalsystem1,attr)
                sol = getattr(chemsys,attr)
                self.assertEqual(test, sol)
        _ = attributes.pop(-1)
        for attr in attributes: 
            with self.subTest(energy_parse='relative',attribute=attr):
                test = getattr(chemicalsystem2,attr)
                sol = getattr(chemsys,attr)
                self.assertEqual(test, sol)
        with self.subTest(energy_parse='relative',attribute='transitionstates'):
            test = chemicalsystem2.transitionstates
            sol = [TS01,TSa]
            self.assertEqual(test,sol)
        for compound in chemicalsystem1.compounds:
            compound2 = chemicalsystem2.Name2Compound[compound.label] 
            with self.subTest(test='energy consistency',compound=compound.label):
                self.assertEqual(compound.energy,compound2.energy)
        ts = chemicalsystem1.Name2TS['TS00']
        ts2 = chemicalsystem2.Name2TS['TS01']
        with self.subTest(test='energy consistency',TS='defaultname'):
            self.assertTrue(ts.energy,ts2.energy)
        with self.subTest(test='energy consistency',TS='TSa'):
            ts = chemicalsystem1.Name2TS['TSa']
            ts2 = chemicalsystem2.Name2TS['TSa']
            self.assertTrue(ts.energy,ts2.energy)
    def test_read_compounds(self):
        module_name = 'pykinetic2.userinput'
        raw_compounds = """
        A   20      kcal/mol    scan
        B   10      hartree
        C   1   
        D   0.5     eV
        """
        solutions = [['A','20','kcal/mol','scan'],
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
        compounds = [['A','20','kcal/mol','scan'],
                     ['B','10','hartree'],
                     ['C','1'],
                     ['D','0.5','eV']]
        solutions = [Compound('A',Energy('20','kcal/mol'),scannable=True),
                     Compound('B',Energy('10','hartree')),
                     Compound('C',Energy('1',default_unit)),
                     Compound('D',Energy('0.5','eV'))]
        tests = create_compounds(compounds,default_unit)
        for test,solution in zip(tests,solutions):
            with self.subTest(compound=solution.label): 
                self.assertEqual(test,solution)
                self.assertEqual(test.energy,solution.energy)
                self.assertEqual(test.scannable,solution.scannable)
        with self.subTest(test='raises error'):
            with self.assertRaises(ValueError):
                _ = create_compounds(compounds + compounds)
    def test_read_reactions(self):
        module_name = 'pykinetic2.userinput'
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
                test = prepare_inline_TS(arg,TSdict,default_unit)
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

    
    