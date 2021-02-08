from pykinetic2.Classes import Compound
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
