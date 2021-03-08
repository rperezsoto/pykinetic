from pykinetic2.classes import *
import unittest

class CompoundTest(unittest.TestCase):
    def setUp(self):
        self.label1 = 'label1'
        self.label2 = 'label2'
        self.energy1 = 0.0
        self.energy2 = 5.0
        self.energy3 = -5.0
        self.key1 = 0
        self.key2 = 'Key2'

    def test_init(self):
        C0 = Compound(self.label1,self.energy1,self.key1,True)
        C1 = Compound(self.label2,self.energy2)
        self.assertEqual(C0.label, self.label1,'Error constructing object')
        self.assertEqual(C0.energy,self.energy1,'Error constructing object')
        self.assertEqual(C0.key,self.key1,'Error constructing object')
        self.assertEqual(C1.key,None,'Error constructing object')
        self.assertTrue(C0.scannable,'Error constructing object')
        self.assertFalse(C1.scannable,'Error constructing object')

    # Test Logic
    def test_equal(self):
        C0 = Compound(self.label1,self.energy1,self.key1)
        C1 = Compound(self.label2,self.energy2)
        C2 = Compound(self.label1,self.energy3,self.key2)
        self.assertEqual(C0,C2)
        self.assertEqual(C2,C0)
        self.assertEqual(C0 == C1,False)
        self.assertEqual(C2 == C1,False)
        self.assertEqual(C1 == self.label2,False)

    # Test Hashing
    def test_hash(self):
        C0 = Compound(self.label1,self.energy1,self.key1)
        C1 = Compound(self.label2,self.energy2)
        C2 = Compound(self.label1,self.energy3,self.key2)
        C3 = Compound(self.label1,self.energy3,self.key1)
        self.assertEqual(hash(C0),hash(C0))
        self.assertNotEqual(hash(C0),hash(C1))
        self.assertEqual(hash(C0),hash(C2))
        self.assertEqual(hash(C0),hash(C3))

    # Test number conversion
    def test_int(self):
        C0 = Compound(self.label1,self.energy1,self.key1)
        C1 = Compound(self.label2,self.energy2)
        C2 = Compound(self.label1,self.energy3,self.key2)
        C3 = Compound(self.label1,self.energy3,self.key1)
        self.assertEqual(int(C0),int(self.key1))
        self.assertEqual(int(C0),int(C3))
        with self.assertRaises(ValueError):
            _ = int(C2)
        with self.assertRaises(TypeError):
            _ = int(C1)

class ReactionTest(unittest.TestCase):
    def setUp(self):
        self.Compound1 = Compound('C1',Energy(0.0,'kcal/mol'))
        self.Compound2 = Compound('C2',Energy(0.0,'kcal/mol'))
        self.Compound3 = Compound('C3',Energy(-2.0,'kcal/mol'))
        self.TS = TransitionState(Energy(10.0,'kcal/mol'),label='TS1')

    def test_init(self):
        reactants = (self.Compound1, self.Compound2)
        products  = (self.Compound3,)
        reaction = Reaction(reactants,products)
        self.assertTrue(reaction.reactants[self.Compound1],1)
        self.assertTrue(reaction.reactants[self.Compound2],1)
        self.assertTrue(reaction.products[self.Compound3],1)
        self.assertIsNone(reaction.key)
        self.assertIsNone(reaction.TS)

    def test_coefficient(self):
        reactants = (self.Compound1, self.Compound1, self.Compound2)
        products  = (self.Compound3, self.Compound1)
        reaction = Reaction(reactants,products)
        self.assertEqual(reaction.coefficient(self.Compound1),-1)
        self.assertEqual(reaction.coefficient(self.Compound2),-1)
        self.assertEqual(reaction.coefficient(self.Compound3),1)

    # Test Hashing
    def test_properties(self):
        reactants = (self.Compound1, self.Compound2)
        products  = (self.Compound3,)
        reaction = Reaction(reactants,products)
        reaction.TS = self.TS
        with self.subTest(property='reactants_energy'):
            self.assertEqual(reaction.reactants_energy,Energy(0.0,'kcal/mol'))
        with self.subTest(property='AE'):
            self.assertEqual(reaction.AE,Energy(10.0,'kcal/mol'))
            reaction.AE = Energy(0.0,'J/mol')
            self.assertEqual(reaction.AE,Energy(0.0,'J/mol'))
        with self.subTest(property='k'):
            k = 298.15*1.38064852E-023/6.62607004E-034
            self.assertAlmostEqual(k,reaction.k)

    # Test number conversion
    def test_errors(self):
        reactants = (self.Compound1, self.Compound2)
        products  = (self.Compound3,)
        reaction = Reaction(reactants,products)
        with self.assertRaises(RuntimeError):
            reaction.AE = 10
        reaction.TS = self.TS
        with self.assertRaises(ValueError):
            reaction.AE = -1.0
            _ = reaction.k

class MassBalanceTest(unittest.TestCase):
    def setUp(self):
        C1 = self.Compound1 = Compound('C1',Energy(0.0,'kcal/mol'))
        C2 = self.Compound2 = Compound('C2',Energy(0.0,'kcal/mol'))
        C3 = self.Compound3 = Compound('C3',Energy(-2.0,'kcal/mol'))
        self.reaction1 = Reaction(reactants=(C1,C2),products=(C3,))
        self.reaction2 = Reaction(reactants=(C3,),products=(C1,C2))
        self.reaction3 = Reaction(reactants=(C1,C1),products=(C1,C3))

    def test_partial(self):
        C1 = Compound('C1',Energy(0.0,'kcal/mol'))
        C2 = Compound('C2',Energy(0.0,'kcal/mol'))
        C3 = Compound('C3',Energy(-2.0,'kcal/mol'))
        r1 = Reaction(reactants=(C1,C2),products=(C3,))
        r2 = Reaction(reactants=(C3,),products=(C1,C2))
        r3 = Reaction(reactants=(C1,C1),products=(C1,C3))
        r1.key = 1
        r2.key = 2 
        r3.key = 3
        MB1 = MassBalance(C1,items=[(-1,r1),(1,r2),(-1,r3)])
        MB2 = MassBalance(C2,items=[(-1,r1),(1,r2)])
        MB3 = MassBalance(C3,items=[(1,r1),(-1,r2),(1,r3)])
        partials = [[JacobianElement(C1,C1,items=[(-1,r1),(1,r2),(-1,r3)]),
                     JacobianElement(C1,C2,items=[(-1,r1),(1,r2)]),
                     JacobianElement(C1,C3,items=[(-1,r1),(1,r2),(-1,r3)])],
                    [JacobianElement(C2,C1,items=[(-1,r1),(1,r2)]),
                     JacobianElement(C2,C2,items=[(-1,r1),(1,r2)]),
                     JacobianElement(C2,C3,items=[(-1,r1),(1,r2)])],
                    [JacobianElement(C3,C1,items=[(1,r1),(-1,r2),(1,r3)]),
                     JacobianElement(C3,C2,items=[(1,r1),(-1,r2)]),
                     JacobianElement(C3,C3,items=[(1,r1),(-1,r2),(1,r3)])]]
        for i,(MB,elements) in enumerate(zip([MB1,MB2,MB3],partials)): 
            for j,(C,element) in enumerate(zip([C1,C2,C3],elements)):
                with self.subTest(MB=i+1,partial=j+1):
                    test = MB.partial(C).items
                    self.assertEqual(test,element.items) 