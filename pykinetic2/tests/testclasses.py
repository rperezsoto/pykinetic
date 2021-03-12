from pykinetic.classes import *
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

class EnergyTest(unittest.TestCase):
    def setUp(self):
        u1,u2 = 'kcal/mol','kJ/mol'
        factors = Energy._ConversionFactors
        self.zero = Energy(0.0,u1)  # zero
        self.pos1 = Energy(1.0,u1) # positive
        self.pos2 = Energy(2.0,u1) # higher than pos1
        self.neg1 = Energy(-1.0,u1) # negative
        self.neg2 = Energy(-2.0,u1) # lower than neg1
        self.refpos = Energy(1.0/factors[u1])
        self.refneg = Energy(-1.0/factors[u1])
        val = factors[u2]/factors[u1]
        self.unitpos = Energy(val,u2)
        self.unitneg = Energy(-val,u2)

    def test_init(self):
        e0 = Energy(unit='kcal/mol')
        e1 = Energy( 1.2,'kcal/mol')
        e2 = Energy(-6.3)
        e3 = Energy(e0,'kJ/mol')
        self.assertEqual(e0.value, 0.0,'Error initializing default value')
        self.assertEqual(e1.value, 1.2,'Error constructing object')
        self.assertEqual(e1._unit,'kcal/mol','Error constructing object')
        self.assertEqual(e2._unit,'hartree','Error initializing default unit')
        self.assertEqual(e3.value == e0.value, True)
        self.assertEqual(e3.unit == e0.unit, False)
    # Test Properties
    def test_value(self):
        compare = lambda x,y: abs(x-y)<=1E-12
        Test = Energy(5.0,'kcal/mol')
        self.assertEqual(compare(Test.value,5.0),True)
        self.assertEqual(compare(Test._value,5.0/Test._conversion),True)
        Test.value = 0.0
        self.assertEqual(compare(Test._value,0.0),True)
        Test._value = 5.0/Test._conversion
        self.assertEqual(compare(Test.value,5.0),True)

    def test_unit(self):
        Test = Energy(5.0,'kcal/mol')
        self.assertEqual(Test.unit,'kcal/mol')
        self.assertEqual(Test._unit,'kcal/mol')
        Val = Test._value
        Test.unit = 'hartree'
        self.assertEqual(Test.value,Val)
        self.assertEqual(Test.unit,'hartree')
        self.assertEqual(Test._unit,'hartree')
        with self.assertRaises(NotImplementedError):
            Test.unit = 'FakeUnit'

    # Test Logic
    def test_bool(self):
        self.assertEqual(bool(self.zero), False)
        self.assertEqual(bool(self.pos1), True)
    def test_equal(self):
        self.assertEqual(self.pos1 == self.unitpos, True)
        self.assertEqual(self.zero == self.pos1, False)
        self.assertEqual(self.pos1 == self.pos1.value, True)
        self.assertEqual(self.pos1 == self.zero.value,False)
    def test_gt(self):
        self.assertEqual(self.pos2 > self.pos1, True)
        self.assertEqual(self.pos2 > self.pos1.value, True)
        self.assertEqual(self.pos2.value > self.pos1, True)
        self.assertEqual(self.neg1 > self.zero.value, False)
        self.assertEqual(self.neg1 > self.neg2,True)
    def test_lt(self):
        self.assertEqual(self.pos2 < self.pos1, False)
        self.assertEqual(self.pos2 < self.pos1.value, False)
        self.assertEqual(self.pos2.value < self.pos1, False)
        self.assertEqual(self.neg1 < self.zero.value, True)
        self.assertEqual(self.neg1 < self.neg2,False)
    def test_ne(self):
        self.assertEqual(self.pos1 != self.pos2, True)
        self.assertEqual(self.pos1 != abs(self.neg1),False)
        self.assertEqual(self.pos1.value != self.pos2, True)

    # Test number conversion
    def test_int(self):
        self.assertEqual(int(self.pos1),int(self.pos1.value))
    def test_float(self):
        self.assertEqual(float(self.pos1),self.pos1.value)
    def test_round(self):
        self.assertEqual(round(self.pos1,5),round(self.pos1.value,5))

    # Test Math
    def test_add_1(self):
        #Test if math addition works properly
        Out = [ (self.pos1+self.zero,self.pos1), # + 0
                (self.neg1+self.zero,self.neg1), # - 0
                (self.pos1+self.pos1,self.pos2), # + +
                (self.pos1+self.neg1,self.zero), # + -
                (self.neg1+self.neg1,self.neg2)] # - -
        for test,result in Out:
            self.assertEqual(test.value,result)
    def test_add_2(self):
        # Test effect of different units and ensure Energy Class Output
        Items = [self.pos1+self.unitpos,
                self.pos1+self.refpos,
                self.unitpos+self.pos1,
                self.pos1+self.pos1.value]
        Out = [ (Items[0],self.pos2, True,self.pos1.unit),
                (Items[1],self.pos2, True,self.pos1.unit),
                (Items[2],self.pos2, True,self.unitpos.unit),
                (Items[3],self.pos2, True,self.pos1.unit)]
        for item,test,result,unit in Out:
            self.assertEqual(item,test)
            self.assertEqual(item.unit,unit)
            self.assertIsInstance(item,Energy)

    def test_sub_1(self):
        #Test if math substraction works properly
        Out = [ (self.zero-self.pos1,self.neg1), # 0 +
                (self.zero-self.neg1,self.pos1), # 0 -
                (self.pos1-self.pos1,self.zero), # + +
                (self.pos1-self.zero,self.pos1), # + -
                (self.pos1-self.neg1,self.pos2)] # - -
        for test,result in Out:
            self.assertEqual(test.value,result)
    def test_sub_2(self):
        # Test effect of different units and ensure Energy Class Output
        Out = [ (self.pos1-self.unitpos,self.zero,self.pos1.unit),
                (self.pos1-self.refpos,self.zero,self.pos1.unit),
                (self.unitpos-self.pos1,self.zero,self.unitpos.unit),
                (self.pos1-self.pos1.value,self.zero,self.pos1.unit)]
        for test,result,unit in Out:
            self.assertEqual(test.value,result)
            self.assertEqual(test.unit,unit)
            self.assertIsInstance(test,Energy)

    def test_mul(self):
        self.assertEqual(self.pos1*2.0,self.pos2)
        self.assertEqual(self.pos1*(-1.0),self.neg1)
        self.assertEqual(self.zero*2.0,self.zero)
        with self.assertRaises(NotImplementedError):
            self.zero*self.zero

    def test_div_1(self):
        # Test correct behaviour of division
        self.assertEqual(self.pos2/2.0,self.pos1)
        self.assertEqual(self.pos1/(-1.0),self.neg1)
        self.assertEqual(self.zero/2.0,self.zero)
        self.assertEqual(self.pos1/self.pos1,1.0)
    def test_div_2(self):
        # Test instance/unit consistency
        test1 = self.pos1/self.unitpos
        test2 = self.pos1/2.0
        self.assertIsInstance(test1,float)
        self.assertIsInstance(test2,Energy)
        self.assertEqual(test2.unit,self.pos1.unit)
        self.assertEqual(test2.value,self.pos1.value/2.0)

    def test_neg(self):
        self.assertEqual(-self.pos1,self.neg1)
        self.assertEqual(-self.neg1,self.pos1)
        self.assertEqual(-self.zero,self.zero)

    def test_abs(self):
        self.assertEqual(abs(self.pos1),self.pos1)
        self.assertEqual(abs(self.neg1),self.pos1)
        self.assertEqual(abs(self.zero),self.zero)

    def test_inline(self):
        # Test +=, -=, *= and /=
        Test = Energy(self.zero.value,'kcal/mol')
        val = self.pos1.value
        Test += val
        self.assertEqual(Test.value,val,'incorrect iadd implementation')
        Test -= val
        self.assertEqual(Test == self.zero, True,'incorrect isub implementation')
        Test += self.pos2
        self.assertEqual(Test,self.pos2,'incorrect iadd implementation')
        Test -= self.pos2
        self.assertEqual(Test == self.zero, True,'incorrect isub implementation')
        Test += self.pos1
        Test *= 0.5
        self.assertEqual(Test,self.pos1*0.5,'incorrect imul implementation')
        Test /= 0.5
        self.assertEqual(Test,self.pos1,'incorrect itruediv implementation')
        Test /= self.pos1
        self.assertEqual(Test,1.0,'incorrect itruediv implementation')

    def test_commutative_eq(self):
        self.assertEqual(self.zero == self.pos1, self.pos1 == self.zero)
    def test_commutative_add(self):
        self.assertEqual(self.pos1 + self.zero, self.zero + self.pos1)
        self.assertEqual(self.pos1 + 1.0, 1.0 + self.pos1)
    def test_commutative_sub(self):
        self.assertEqual(self.pos1 - self.pos2, -self.pos2 + self.pos1)
        self.assertEqual(self.pos1 - 1.0, -1.0 - (-self.pos1))
    def test_commutative_mul(self):
        self.assertEqual(self.pos1*2.0,2.0*self.pos1)

    # Test Methods
    def test_as_unit(self):
        val,unit = self.pos1.value, self.pos1.unit
        unit2 = self.unitpos.unit
        Out1 = self.pos1.as_unit(unit2)
        self.assertIn(str(self.unitpos.value),Out1)
        self.assertIn(unit2,Out1)
        with self.assertRaises(NotImplementedError):
            self.pos1.as_unit('FakeUnit')
