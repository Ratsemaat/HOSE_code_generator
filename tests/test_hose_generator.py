import unittest
from rdkit import Chem

from hosegen import HoseGenerator

eth = Chem.MolFromSmiles("CC")
mol1 = Chem.MolFromSmiles("CCCC")
mol2 = Chem.MolFromSmiles("Cc1ccccc1")


class HoseGeneratorTest(unittest.TestCase):
    def setUp(self):
        self.gen = HoseGenerator()

    def test_simple(self):
        value = self.gen.get_Hose_codes(eth, 0)
        print(value)

    def test_mol1_1(self):
        value = self.gen.get_Hose_codes(mol1, 0)
        assert value == "C-4;C(C/C/)/"

    def test_mol1_2(self):
        value = self.gen.get_Hose_codes(mol1, 1)
        assert value == "C-4;CC(C,//)/"

    def test_mol1_3(self):
        value = self.gen.get_Hose_codes(mol1, 2)
        assert value == "C-4;CC(C,//)/"

    def test_mol1_4(self):
        value = self.gen.get_Hose_codes(mol1, 3)
        assert value == "C-4;C(C/C/)/"


    def test_mol2_1(self):

        value = self.gen.get_Hose_codes(mol2, 0)
        assert value == "C-4;C(*C*C/*C,*C/*C,*&)*&/"

        value = self.gen.get_Hose_codes(mol2, 1)
        assert value == "C-3;*C*CC(*C,*C,/*C,*&/*&)/"

        value = self.gen.get_Hose_codes(mol2, 2)
        "C-3;*C*C(*CC,*C/*C,,*&/*&)/"
        "C-3;*C*C(*C,*C,C/*C,*&,/*&)/"
        "C-3;*C*C(*C,*C,C/*C,*&,/*&)/"
        "C-3;*C*C(*C,*C/*CC,*&/*&,)/"

        assert value == "C-3;*C*C(*CC,*C/*C,,*&/*&)/"

        value = self.gen.get_Hose_codes(mol2, 6)
        assert value == "C-3;*C*C(*CC,*C/*C,,*&/*&)/"

        value = self.gen.get_Hose_codes(mol2, 3)
        assert value == "C-3;*C*C(*C,*C/*CC,*&/*&,)/"

        value = self.gen.get_Hose_codes(mol2, 5)
        assert value == "C-3;*C*C(*C,*C/*CC,*&/*&,)/"

        value = self.gen.get_Hose_codes(mol2, 4)
        assert value == "C-3;*C*C(*C,*C/*C,*&/*&C),/"



if __name__ == '__main__':
    import xmlrunner

    unittest.main(
        testRunner=xmlrunner.XMLTestRunner(output='test-reports'),
        failfast=False,
        buffer=False,
        catchbreak=False)