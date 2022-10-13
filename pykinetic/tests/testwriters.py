import unittest
from pykinetic.classes import Compound, Energy, ChemicalSystem, Reaction, TransitionState
from pykinetic.writers import Indent, PythonWriter, CplusplusWriter

class IndentTest(unittest.TestCase):
    def test_default(self):
        items =  [(['Test',]        ,['    Test',]),
                  (['Test1','Test2'],['    Test1','    Test2']),
                  ( 'Test'          ,['    T','    e','    s','    t']),
                  (   ''            ,[]),
                  (   []            ,[]),]
        for t,s in items:
            with self.subTest(t=t):
                self.assertEqual(Indent(t),s)
    def test_level(self):
        items =  [(1, ['Test',]        ,['\tTest',]),
                  (1, ['Test1','Test2'],['\tTest1','\tTest2']),
                  (2, ['Test1','Test2'],['\t\tTest1','\t\tTest2']),
                  (2,  'Test'          ,['\t\tT','\t\te','\t\ts','\t\tt']),
                  (4,    []            ,[]),]
        for i,t,s in items:
            with self.subTest(t=t):
                self.assertEqual(Indent(t,'\t',level=i),s)
    def test_level_consistency(self):
        test = 'line0 line1 line2 line3 line4 line-1'.split()
        for i,j in zip(range(4),range(1,5)):
            s = Indent(test,level=j)
            t = Indent(Indent(test,level=i),level=1)
            with self.subTest(i=i,j=j):
                self.assertEqual(t,s)

class PythonWriterTest(unittest.TestCase):
    def setUp(self):
        self.writer = PythonWriter()
        unit = 'kcal/mol'
        A = Compound('A',Energy( 0.0,unit))
        B = Compound('B',Energy( 0.0,unit))
        C = Compound('C',Energy( 2.0,unit))
        D = Compound('D',Energy(-1.0,unit))
        E = Compound('E',Energy(  99,unit)) # Does not appear in reactions
        # TransitionState(Energy(1.0,unit))
        self.compounds = [A,B,C,D,E]
        self.chemsys = ChemicalSystem()
        for c in self.compounds:
            self.chemsys.cadd(c,update=False)
        self.chemsys.cupdate()
        self.reactions = [Reaction((A,),(C,)),
                          Reaction((A,B),(C,)),
                          Reaction((A,A),(C,)),
                          Reaction((A,),(B,C)),
                          Reaction((A,),(C,B)),
                          Reaction((A,),(C,C)),
                          Reaction((A,B),(C,D)),
                          Reaction((A,B,C),(D,)),
                          Reaction((A,),(B,C,D))]
        for r in self.reactions:
            r.TS = TransitionState(Energy(1.0,'kcal/mol') + r.reactants_energy)
            self.chemsys.radd(r,update=False)
        self.chemsys.rupdate()

    def test_ratelaw_expr(self):
        solutions = ['k00*x[0]',
                     'k01*x[0]*x[1]',
                     'k02*x[0]*x[0]',
                     'k03*x[0]',
                     'k04*x[0]',
                     'k05*x[0]',
                     'k06*x[0]*x[1]',
                     'k07*x[0]*x[1]*x[2]',
                     'k08*x[0]']
        for reaction,solution in zip(self.reactions,solutions):
            with self.subTest(reaction=reaction):
                _,test = self.writer.ratelaw(reaction)
                self.assertEqual(test,solution)
    def test_ratelaw_var(self):
        solutions = ['r00','r01','r02','r03','r04','r05','r06','r07','r08']
        for reaction,solution in zip(self.reactions,solutions):
            with self.subTest(reaction=reaction):
                test,_  = self.writer.ratelaw(reaction)
                self.assertEqual(test,solution)
    def test_ratelaw_partial(self):
        solutions = ['k00',
                     'k01*x[1]',
                     '2.0*k02*x[0]',
                     'k03',
                     'k04',
                     'k05',
                     'k06*x[1]',
                     'k07*x[1]*x[2]',
                     '',
                     '']
        compounds = [self.compounds[0] for i in solutions[:-2]]
        compounds.append(self.compounds[-2]) # Add Product
        compounds.append(self.compounds[-1]) # Add external
        for r,s,c in zip(self.reactions,solutions,compounds):
            with self.subTest(reaction=str(r),compound=c):
                test = self.writer.ratelaw_partial(r,c)
                self.assertEqual(test,s)
    def test_massbalance_expr(self):
        solutions = ['-r00-r01-2.0*r02-r03-r04-r05-r06-r07-r08',
                     '-r01+r03+r04-r06-r07+r08',
                     '+r00+r01+r02+r03+r04+2.0*r05+r06-r07+r08',
                     '+r06+r07+r08',
                     '0',
                     '0']
        compounds = [c for c in self.compounds]
        compounds.append(Compound('ClearlyNotInTheSystem',5.2))
        for s,c in zip(solutions,compounds):
            with self.subTest(compound=c):
                ## Future MB = self.chemsys.massbalance(c)
                MB = self.chemsys.massbalance(c)
                _ , test = self.writer.massbalance(MB)
                self.assertEqual(test,s)
    def test_massbalance_var(self):
        solutions = ['dxdt[0]',
                     'dxdt[1]',
                     'dxdt[2]',
                     'dxdt[3]',
                     '',
                     '']
        compounds = [c for c in self.compounds]
        compounds.append(Compound('ClearlyNotInTheSystem',5.2))
        for s,c in zip(solutions,compounds):
            with self.subTest(compound=c):
                ## Future MB = self.chemsys.massbalance(c)
                MB = self.chemsys.massbalance(c)
                test, _ = self.writer.massbalance(MB)
                self.assertEqual(test,s)
    def test_jacobian_element_expr(self):
        solutions = [['-k00-k01*x[1]-2.0*2.0*k02*x[0]-k03-k04-k05-k06*x[1]-k07*x[1]*x[2]-k08',
                      '-k01*x[0]-k06*x[0]-k07*x[0]*x[2]',
                      '-k07*x[0]*x[1]',
                      '0',
                      '0',
                      '0'],
                     ['-k01*x[1]+k03+k04-k06*x[1]-k07*x[1]*x[2]+k08',
                      '-k01*x[0]-k06*x[0]-k07*x[0]*x[2]',
                      '-k07*x[0]*x[1]',
                      '0',
                      '0',
                      '0'],
                     ['+k00+k01*x[1]+2.0*k02*x[0]+k03+k04+2.0*k05+k06*x[1]-k07*x[1]*x[2]+k08',
                      '+k01*x[0]+k06*x[0]-k07*x[0]*x[2]',
                      '-k07*x[0]*x[1]',
                      '0',
                      '0',
                      '0'],
                     ['+k06*x[1]+k07*x[1]*x[2]+k08',
                      '+k06*x[0]+k07*x[0]*x[2]',
                      '+k07*x[0]*x[1]',
                      '0',
                      '0',
                      '0'],
                     ['0']*6]
        compounds = [c for c in self.compounds]
        compounds.append(Compound('ClearlyNotInTheSystem',5.2))
        for sols,c in zip(solutions,compounds):
            ## Future MB = self.chemsys.massbalance(c)
            MB = self.chemsys.massbalance(c)
            ## Future Jacobian
            for c2,s in zip(compounds,sols):
                Jac_ij = MB.partial(c2)
                with self.subTest(compound1=c,compound2=c2):
                    _, test = self.writer.jacobian_element(Jac_ij)
                    self.assertEqual(test,s)
    def test_jacobian_element_var(self):
        solutions = [['Jac[0,0]',
                      'Jac[0,1]',
                      'Jac[0,2]',
                      'Jac[0,3]',
                      '',
                      ''],
                     ['Jac[1,0]',
                      'Jac[1,1]',
                      'Jac[1,2]',
                      'Jac[1,3]',
                      '',
                      ''],
                     ['Jac[2,0]',
                      'Jac[2,1]',
                      'Jac[2,2]',
                      'Jac[2,3]',
                      '',
                      ''],
                     ['Jac[3,0]',
                      'Jac[3,1]',
                      'Jac[3,2]',
                      'Jac[3,3]',
                      '',
                      ''],
                     ['']*6,
                     ['']*6]
        compounds = [c for c in self.compounds]
        compounds.append(Compound('ClearlyNotInTheSystem',5.2))
        for sols,c in zip(solutions,compounds):
            ## Future MB = self.chemsys.massbalance(c)
            MB = self.chemsys.massbalance(c)
            ## Future Jacobian
            for c2,s in zip(compounds,sols):
                Jac_ij = MB.partial(c2)
                with self.subTest(compound1=c,compound2=c2):
                    test, _ = self.writer.jacobian_element(Jac_ij)
                    self.assertEqual(test,s)

    # test "private" methods
    def test__kinetic_constants(self):
        solution = ['#Constants',
                    'k00 = 1.1488119097e+12',
                    'k01 = 1.1488119097e+12',
                    'k02 = 1.1488119097e+12',
                    'k03 = 1.1488119097e+12',
                    'k04 = 1.1488119097e+12',
                    'k05 = 1.1488119097e+12',
                    'k06 = 1.1488119097e+12',
                    'k07 = 1.1488119097e+12',
                    'k08 = 1.1488119097e+12']
        test = self.writer._kinetic_constants(self.chemsys)
        self.assertEqual(test,solution)
    def test__ratelaws(self):
        solution = ['#Ratelaws',
                    'r00 = k00*x[0]',
                    'r01 = k01*x[0]*x[1]',
                    'r02 = k02*x[0]*x[0]',
                    'r03 = k03*x[0]',
                    'r04 = k04*x[0]',
                    'r05 = k05*x[0]',
                    'r06 = k06*x[0]*x[1]',
                    'r07 = k07*x[0]*x[1]*x[2]',
                    'r08 = k08*x[0]']
        test = self.writer._ratelaws(self.chemsys)
        self.assertEqual(test,solution)
    def test__massbalances(self):
        solution = ['#MassBalances',
                    'dxdt[0] = -r00-r01-2.0*r02-r03-r04-r05-r06-r07-r08',
                    'dxdt[1] = -r01+r03+r04-r06-r07+r08',
                    'dxdt[2] = +r00+r01+r02+r03+r04+2.0*r05+r06-r07+r08',
                    'dxdt[3] = +r06+r07+r08']
        test = self.writer._massbalances(self.chemsys)
        self.assertEqual(test,solution)
    def test__jacobian_elements(self):
        solution = ['#Non-zero Elements',
                    'Jac[0,0] = -k00-k01*x[1]-2.0*2.0*k02*x[0]-k03-k04-k05-k06*x[1]-k07*x[1]*x[2]-k08',
                    'Jac[0,1] = -k01*x[0]-k06*x[0]-k07*x[0]*x[2]',
                    'Jac[0,2] = -k07*x[0]*x[1]',
                    'Jac[1,0] = -k01*x[1]+k03+k04-k06*x[1]-k07*x[1]*x[2]+k08',
                    'Jac[1,1] = -k01*x[0]-k06*x[0]-k07*x[0]*x[2]',
                    'Jac[1,2] = -k07*x[0]*x[1]',
                    'Jac[2,0] = +k00+k01*x[1]+2.0*k02*x[0]+k03+k04+2.0*k05+k06*x[1]-k07*x[1]*x[2]+k08',
                    'Jac[2,1] = +k01*x[0]+k06*x[0]-k07*x[0]*x[2]',
                    'Jac[2,2] = -k07*x[0]*x[1]',
                    'Jac[3,0] = +k06*x[1]+k07*x[1]*x[2]+k08',
                    'Jac[3,1] = +k06*x[0]+k07*x[0]*x[2]',
                    'Jac[3,2] = +k07*x[0]*x[1]']
        test = self.writer._jacobian_elements(self.chemsys)
        self.assertEqual(test,solution)

class CplusplusWriterTest(unittest.TestCase):
    def setUp(self):
        self.writer = CplusplusWriter()
        unit = 'kcal/mol'
        A = Compound('A',Energy( 0.0,unit))
        B = Compound('B',Energy( 0.0,unit))
        C = Compound('C',Energy( 2.0,unit))
        D = Compound('D',Energy(-1.0,unit))
        E = Compound('E',Energy(  99,unit)) # Does not appear in reactions
        self.compounds = [A,B,C,D,E]
        self.chemsys = ChemicalSystem()
        for c in self.compounds:
            self.chemsys.cadd(c,update=False)
        self.chemsys.cupdate()
        self.reactions = [Reaction((A,),(C,)),
                          Reaction((A,B),(C,)),
                          Reaction((A,A),(C,)),
                          Reaction((A,),(B,C)),
                          Reaction((A,),(C,B)),
                          Reaction((A,),(C,C)),
                          Reaction((A,B),(C,D)),
                          Reaction((A,B,C),(D,)),
                          Reaction((A,),(B,C,D))]
        for r in self.reactions:
            r.TS = TransitionState(Energy(1.0,'kcal/mol') + r.reactants_energy)
            self.chemsys.radd(r,update=False)
        self.chemsys.rupdate()

    def test_ratelaw_expr(self):
        solutions = ['k00*x[0]',
                     'k01*x[0]*x[1]',
                     'k02*x[0]*x[0]',
                     'k03*x[0]',
                     'k04*x[0]',
                     'k05*x[0]',
                     'k06*x[0]*x[1]',
                     'k07*x[0]*x[1]*x[2]',
                     'k08*x[0]']
        for reaction,solution in zip(self.reactions,solutions):
            with self.subTest(reaction=reaction):
                _,test = self.writer.ratelaw(reaction)
                self.assertEqual(test,solution)
    def test_ratelaw_var(self):
        solutions = ['r00','r01','r02','r03','r04','r05','r06','r07','r08']
        for reaction,solution in zip(self.reactions,solutions):
            with self.subTest(reaction=reaction):
                test,_  = self.writer.ratelaw(reaction)
                self.assertEqual(test,solution)
    def test_ratelaw_partial(self):
        solutions = ['k00',
                     'k01*x[1]',
                     '2.0*k02*x[0]',
                     'k03',
                     'k04',
                     'k05',
                     'k06*x[1]',
                     'k07*x[1]*x[2]',
                     '',
                     '']
        compounds = [self.compounds[0] for i in solutions[:-2]]
        compounds.append(self.compounds[-2]) # Add Product
        compounds.append(self.compounds[-1]) # Add external
        for r,s,c in zip(self.reactions,solutions,compounds):
            with self.subTest(reaction=str(r),compound=c):
                test = self.writer.ratelaw_partial(r,c)
                self.assertEqual(test,s)
    def test_massbalance_expr(self):
        solutions = ['-r00-r01-2.0*r02-r03-r04-r05-r06-r07-r08',
                     '-r01+r03+r04-r06-r07+r08',
                     '+r00+r01+r02+r03+r04+2.0*r05+r06-r07+r08',
                     '+r06+r07+r08',
                     '0',
                     '0']
        compounds = [c for c in self.compounds]
        compounds.append(Compound('ClearlyNotInTheSystem',5.2))
        for s,c in zip(solutions,compounds):
            with self.subTest(compound=c):
                ## Future MB = self.chemsys.massbalance(c)
                MB = self.chemsys.massbalance(c)
                _ , test = self.writer.massbalance(MB)
                self.assertEqual(test,s)
    def test_massbalance_var(self):
        solutions = ['dxdt[0]',
                     'dxdt[1]',
                     'dxdt[2]',
                     'dxdt[3]',
                     'dxdt[4]',
                     '']
        compounds = [c for c in self.compounds]
        compounds.append(Compound('ClearlyNotInTheSystem',5.2))
        for s,c in zip(solutions,compounds):
            with self.subTest(compound=c):
                ## Future MB = self.chemsys.massbalance(c)
                MB = self.chemsys.massbalance(c)
                test, _ = self.writer.massbalance(MB)
                self.assertEqual(test,s)
    def test_jacobian_element_expr(self):
        solutions = [['-k00-k01*x[1]-2.0*2.0*k02*x[0]-k03-k04-k05-k06*x[1]-k07*x[1]*x[2]-k08',
                      '-k01*x[0]-k06*x[0]-k07*x[0]*x[2]',
                      '-k07*x[0]*x[1]',
                      '0',
                      '0',
                      '0'],
                     ['-k01*x[1]+k03+k04-k06*x[1]-k07*x[1]*x[2]+k08',
                      '-k01*x[0]-k06*x[0]-k07*x[0]*x[2]',
                      '-k07*x[0]*x[1]',
                      '0',
                      '0',
                      '0'],
                     ['+k00+k01*x[1]+2.0*k02*x[0]+k03+k04+2.0*k05+k06*x[1]-k07*x[1]*x[2]+k08',
                      '+k01*x[0]+k06*x[0]-k07*x[0]*x[2]',
                      '-k07*x[0]*x[1]',
                      '0',
                      '0',
                      '0'],
                     ['+k06*x[1]+k07*x[1]*x[2]+k08',
                      '+k06*x[0]+k07*x[0]*x[2]',
                      '+k07*x[0]*x[1]',
                      '0',
                      '0',
                      '0'],
                     ['0']*6]
        compounds = [c for c in self.compounds]
        compounds.append(Compound('ClearlyNotInTheSystem',5.2))
        for sols,c in zip(solutions,compounds):
            ## Future MB = self.chemsys.massbalance(c)
            MB = self.chemsys.massbalance(c)
            ## Future Jacobian
            for c2,s in zip(compounds,sols):
                Jac_ij = MB.partial(c2)
                with self.subTest(compound1=c,compound2=c2):
                    _, test = self.writer.jacobian_element(Jac_ij)
                    self.assertEqual(test,s)
    def test_jacobian_element_var(self):
        solutions = [['Jac(0,0)',
                      'Jac(0,1)',
                      'Jac(0,2)',
                      'Jac(0,3)',
                      'Jac(0,4)',
                      ''],
                     ['Jac(1,0)',
                      'Jac(1,1)',
                      'Jac(1,2)',
                      'Jac(1,3)',
                      'Jac(1,4)',
                      ''],
                     ['Jac(2,0)',
                      'Jac(2,1)',
                      'Jac(2,2)',
                      'Jac(2,3)',
                      'Jac(2,4)',
                      ''],
                     ['Jac(3,0)',
                      'Jac(3,1)',
                      'Jac(3,2)',
                      'Jac(3,3)',
                      'Jac(3,4)',
                      ''],
                     ['Jac(4,0)',
                      'Jac(4,1)',
                      'Jac(4,2)',
                      'Jac(4,3)',
                      'Jac(4,4)',
                      ''],
                     ['']*6]
        compounds = [c for c in self.compounds]
        compounds.append(Compound('ClearlyNotInTheSystem',5.2))
        for sols,c in zip(solutions,compounds):
            ## Future MB = self.chemsys.massbalance(c)
            MB = self.chemsys.massbalance(c)
            ## Future Jacobian
            for c2,s in zip(compounds,sols):
                Jac_ij = MB.partial(c2)
                with self.subTest(compound1=c,compound2=c2):
                    test, _ = self.writer.jacobian_element(Jac_ij)
                    self.assertEqual(test,s)

    # test "private" methods
    def test__kinetic_constants(self):
        solution = ['// Kinetic Constants',
                    'const double k00 = 1.1488119097e+12;',
                    'const double k01 = 1.1488119097e+12;',
                    'const double k02 = 1.1488119097e+12;',
                    'const double k03 = 1.1488119097e+12;',
                    'const double k04 = 1.1488119097e+12;',
                    'const double k05 = 1.1488119097e+12;',
                    'const double k06 = 1.1488119097e+12;',
                    'const double k07 = 1.1488119097e+12;',
                    'const double k08 = 1.1488119097e+12;']
        test = self.writer._kinetic_constants(self.chemsys)
        self.assertEqual(test,solution)
    def test__ratelaws(self):
        solution = ['// Ratelaws',
                    'double r00,r01,r02,r03,r04,r05,r06,r07,r08;',
                    'r00 = k00*x[0];',
                    'r01 = k01*x[0]*x[1];',
                    'r02 = k02*x[0]*x[0];',
                    'r03 = k03*x[0];',
                    'r04 = k04*x[0];',
                    'r05 = k05*x[0];',
                    'r06 = k06*x[0]*x[1];',
                    'r07 = k07*x[0]*x[1]*x[2];',
                    'r08 = k08*x[0];']
        test = self.writer._ratelaws(self.chemsys)
        self.assertEqual(test,solution)
    def test__massbalances(self):
        solution = ['// MassBalances',
                    'dxdt[0] = -r00-r01-2.0*r02-r03-r04-r05-r06-r07-r08;',
                    'dxdt[1] = -r01+r03+r04-r06-r07+r08;',
                    'dxdt[2] = +r00+r01+r02+r03+r04+2.0*r05+r06-r07+r08;',
                    'dxdt[3] = +r06+r07+r08;',
                    'dxdt[4] = 0;']
        test = self.writer._massbalances(self.chemsys)
        self.assertEqual(test,solution)
