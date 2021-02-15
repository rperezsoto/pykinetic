from pykinetic2.Classes import Compound, Reaction, Energy, TransitionState
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
        C0 = Compound(self.label1,self.energy1,self.key1)
        C1 = Compound(self.label2,self.energy2)
        self.assertEqual(C0.label, self.label1,'Error constructing object')
        self.assertEqual(C0.energy,self.energy1,'Error constructing object')
        self.assertEqual(C0.key,self.key1,'Error constructing object')
        self.assertEqual(C1.key,None,'Error constructing object')

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

    # Test Methods
    def test_DotGraphRepr(self):
        C0 = Compound(self.label1,self.energy1,self.key1)
        self.assertIn(self.label1,C0.DotGraphRepr())

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

