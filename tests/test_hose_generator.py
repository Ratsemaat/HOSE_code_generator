import unittest
from rdkit import Chem
import sys
sys.path.append(".")
from hosegen import HoseGenerator
from hosegen.geometry import *
from hosegen.exceptions import *


class HoseGeneratorTest(unittest.TestCase):
    eth = Chem.MolFromSmiles("CC")
    mol1 = Chem.MolFromSmiles("CCCC")
    mol2 = Chem.MolFromSmiles("Cc1ccccc1")
    mol3 = Chem.MolFromSmiles("F[As@OH1-](F)(F)(F)(F)F")
    molStereo1 = Chem.MolFromMolFile("tests/stereo_tetrahedral1.mol")
    molStereo2 = Chem.MolFromMolFile("tests/stereo_squareplanar1.mol")
    molStereo3 = Chem.MolFromMolFile("tests/stereo_squareplanar2.mol")
    molStereo4 = Chem.MolFromMolFile("tests/stereo_octahedral1.mol")
    molStereo5 = Chem.MolFromMolFile("tests/stereo_octahedral2.mol")
    molStereo6 = Chem.MolFromMolFile("tests/stereo_updown.mol")
    molStereo7 = Chem.MolFromMolFile("tests/stereo_tetrahedral2.mol", removeHs=False)
    molEZ1 = Chem.MolFromMolFile("tests/ez_1.mol", removeHs=False)
    molEZ2 = Chem.MolFromMolFile("tests/ez_2.mol", removeHs=False)
    molEZ3 = Chem.MolFromMolFile("tests/ez_3.mol", removeHs=False)
    molEZ4 = Chem.MolFromMolFile("tests/ez_4.mol", removeHs=False)
    molSymmetry = Chem.MolFromMolFile("tests/symmetryandstopinhose.mol", removeHs=False)
    molZeroBond1 = Chem.MolFromMolFile("tests/zbo_1.mol")

    def setUp(self):
        self.gen = HoseGenerator()


    def test_hydrogen_as_starting_point(self):
        prop = Chem.MolFromSmiles("CCC")
        self.assertRaises(InvalidStartingAtomIndex, self.gen.get_Hose_codes, prop, 3)
        prop = Chem.AddHs(prop)
        value = self.gen.get_Hose_codes(prop, 3)
        assert value == "H-1;C(HHC/HHC/HHH)/"

    def test_simple(self):
        value = self.gen.get_Hose_codes(HoseGeneratorTest.eth, 0)
        assert value == "C-4;C(//)/"

    def test_mol1_1(self):
        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol1, 0)
        assert value == "C-4;C(C/C/)/"

    def test_mol1_2(self):
        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol1, 1)
        assert value == "C-4;CC(C,//)/"

    def test_mol1_3(self):
        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol1, 2)
        assert value == "C-4;CC(C,//)/"

    def test_mol1_4(self):
        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol1, 3)
        assert value == "C-4;C(C/C/)/"


    def test_mol2_1(self):

        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol2, 0)
        assert value == "C-4;C(*C*C/*C,*C/*C,*&)*&/"
        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol2, 0, ringsize=True)
        assert value == "C-4;C(*C*C/*C,*C/*C,*&)*&/"

        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol2, 1)
        assert value == "C-3;*C*CC(*C,*C,/*C,*&/*&)/"
        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol2, 1, ringsize=True)
        assert value == "C-3-6;*C*CC(*C,*C,/*C,*&/*&)/"

        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol2, 2)
        "C-3;*C*C(*CC,*C/*C,,*&/*&)/"
        "C-3;*C*C(*C,*C,C/*C,*&,/*&)/"
        "C-3;*C*C(*C,*C,C/*C,*&,/*&)/"
        "C-3;*C*C(*C,*C/*CC,*&/*&,)/"

        assert value == "C-3;*C*C(*CC,*C/*C,,*&/*&)/"

        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol2, 6)
        assert value == "C-3;*C*C(*CC,*C/*C,,*&/*&)/"

        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol2, 3)
        assert value == "C-3;*C*C(*C,*C/*CC,*&/*&,)/"

        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol2, 5)
        assert value == "C-3;*C*C(*C,*C/*CC,*&/*&,)/"

        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol2, 4)
        assert value == "C-3;*C*C(*C,*C/*C,*&/*&C),/"

    def test_mol3(self):
        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol3, 0)
        assert value == "F-1;As-(FFFFF/,,,,/)/"

        value = self.gen.get_Hose_codes(HoseGeneratorTest.mol3, 1)
        assert value == "As-6-;FFFFFF(,,,,,//)/"

    def test_molStereo1(self):
        value = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo1,0)
        value2 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo1,1)
        assert value == value2
        wedgemap1=create_wedgemap("tests/stereo_tetrahedral1.mol")
        value = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo1,0,usestereo=True,wedgebond=wedgemap1)
        assert value != value2
        wedgemap12=create_wedgemap("tests/stereo_tetrahedral1.mol")
        wedgemap12[20]=1
        value3= self.gen.get_Hose_codes(HoseGeneratorTest.molStereo1,0,usestereo=True,wedgebond=wedgemap12)
        value4 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo1,1,usestereo=True,wedgebond=wedgemap12)
        assert value != value3
        assert value2 != value4
        assert value4 != value3
        value3= self.gen.get_Hose_codes(HoseGeneratorTest.molStereo1,13)
        value4 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo1,13,usestereo=True,wedgebond=wedgemap1)
        assert value4 != value3
        value= self.gen.get_Hose_codes(HoseGeneratorTest.molStereo1,13)
        value2 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo1,13,usestereo=True,wedgebond=wedgemap12)
        assert value != value2
        assert value == value3
        assert value2 != value4

    def test_molStereo2(self):
        wedgemap=create_wedgemap("tests/stereo_squareplanar1.mol")
        value = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo2,0)
        value2 = self.gen.get_Hose_codes_from_file("tests/stereo_squareplanar1.mol",0,usestereo=True)
        assert value != value2
        value3 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo3,0)
        value4 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo3,0,usestereo=True,wedgebond=wedgemap)
        assert value3 != value4
        assert value == value3
        assert value2 != value4

    def test_molStereo3(self):
        wedgemap=create_wedgemap("tests/stereo_octahedral1.mol")
        value = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo4,2)
        value2 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo4,2,usestereo=True,wedgebond=wedgemap)
        assert value != value2
        value3 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo5,2)
        value4 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo5,2,usestereo=True,wedgebond=wedgemap)
        assert value3 != value4
        assert value == value3
        assert value2 != value4
        value = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo4,5)
        value2 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo4,5,usestereo=True,wedgebond=wedgemap)
        assert value != value2
        value3 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo5,5)
        value4 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo5,5,usestereo=True,wedgebond=wedgemap)
        assert value3 != value4
        assert value == value3
        assert value2 != value4
        value = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo4,6)
        value2 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo4,6,usestereo=True,wedgebond=wedgemap)
        assert value != value2
        value3 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo5,6)
        value4 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo5,6,usestereo=True,wedgebond=wedgemap)
        assert value3 != value4
        assert value == value3
        assert value2 != value4

    def test_molStereo4 (self):
        wedgemap=create_wedgemap("tests/stereo_tetrahedral2.mol")
        value = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo7,0,usestereo=True,wedgebond=wedgemap)
        value2 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo7,1,usestereo=True,wedgebond=wedgemap)
        assert value == "C-4;@FHClBr(,,//)/"
        assert value2 == "F-1;C(@ClBrH/,/)/"

    def test_molstereoupdown(self):
        wedgemap=create_wedgemap("tests/stereo_updown.mol")
        value1 = self.gen.get_Hose_codes(HoseGeneratorTest.molStereo6,2,usestereo=True)
        #no stereochemistry recognized here
        assert "@" not in value1
        value2 = self.gen.get_Hose_codes_from_file("tests/stereo_updown.mol",2,usestereo=True)
        #with the up/down mechanism it is
        assert "@" in value2
    
    def test_molez1(self):
        value1 = self.gen.get_Hose_codes(HoseGeneratorTest.molEZ1,4,usestereo=False)
        value2 = self.gen.get_Hose_codes(HoseGeneratorTest.molEZ2,3,usestereo=False)
        assert value1 == value2
        value3 = self.gen.get_Hose_codes(HoseGeneratorTest.molEZ1,4,usestereo=True)
        value4 = self.gen.get_Hose_codes(HoseGeneratorTest.molEZ2,3,usestereo=True)
        assert value3 != value4
        assert value1 != value3
        assert value2 != value4
        value3 = self.gen.get_Hose_codes(HoseGeneratorTest.molEZ1,2,usestereo=True)
        value4 = self.gen.get_Hose_codes(HoseGeneratorTest.molEZ2,2,usestereo=True)
        assert value3 != value4
        value3 = self.gen.get_Hose_codes(HoseGeneratorTest.molEZ3,2,usestereo=True)
        value4 = self.gen.get_Hose_codes(HoseGeneratorTest.molEZ4,2,usestereo=True)
        assert value3 != value4


    def test_symmetry_and_stop(self) :
        valueA = self.gen.get_Hose_codes(HoseGeneratorTest.molSymmetry, 2, 6, True)
        valueB = self.gen.get_Hose_codes(HoseGeneratorTest.molSymmetry, 3, 6, True)

        assert valueB == valueA

        valueA = self.gen.get_Hose_codes(HoseGeneratorTest.molSymmetry, 4, 6, True)
        valueB = self.gen.get_Hose_codes(HoseGeneratorTest.molSymmetry, 6, 6, True)

        assert valueB == valueA


    def test_zbo1(self):
        # Phospourus atoms should have the same HOSE codes
        value1 = self.gen.get_Hose_codes(HoseGeneratorTest.molZeroBond1,6,usestereo=False)
        value2 = self.gen.get_Hose_codes(HoseGeneratorTest.molZeroBond1,7,usestereo=False)
        assert value1 == value2

        # Carbon atoms same distance from non-planar carbon should have the same HOSE codes
        value1 = self.gen.get_Hose_codes(HoseGeneratorTest.molZeroBond1,1,2,usestereo=False)
        value2 = self.gen.get_Hose_codes(HoseGeneratorTest.molZeroBond1,2,2,usestereo=False)
        assert value1 == value2

        # Carbon atoms same distance from non-planar carbon should have the same HOSE codes
        value1 = self.gen.get_Hose_codes(HoseGeneratorTest.molZeroBond1,3,usestereo=False)
        value2 = self.gen.get_Hose_codes(HoseGeneratorTest.molZeroBond1,5,usestereo=False)
        assert value1 == value2

         # Carbon atoms with different distance from non-planar carbon should have the same HOSE codes
        value1 = self.gen.get_Hose_codes(HoseGeneratorTest.molZeroBond1,3,usestereo=False)
        value2 = self.gen.get_Hose_codes(HoseGeneratorTest.molZeroBond1,1,usestereo=False)
        assert value1 != value2



if __name__ == '__main__':
    import xmlrunner

    unittest.main(
        testRunner=xmlrunner.XMLTestRunner(output='test-reports'),
        failfast=False,
        buffer=False,
        catchbreak=False)