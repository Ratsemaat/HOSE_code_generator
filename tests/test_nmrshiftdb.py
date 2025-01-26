# type: ignore
import unittest
from rdkit import Chem
import sys

sys.path.append(".")
from hosegen import HoseGenerator
from hosegen.geometry import *
import pytest


class HoseGeneratorNmrshiftdbTest(unittest.TestCase):
    def setUp(self):
        self.gen = HoseGenerator()

    @pytest.mark.skip(reason="to be run locally as nmrshiftdb2.sd is too large for git")
    def test_all(self):
        file1 = open("tests/nmrshiftdb2.sd", "r")
        count = 0
        molstring = ""
        while True:
            line = file1.readline()
            if not line:
                break
            if line.strip() == "$$$$":
                mol = Chem.MolFromMolBlock(molstring, removeHs=False)
                if mol is not None:
                    mol = Chem.rdmolops.AddHs(mol)
                    wedgemap = create_wedgemap_from_block(
                        molstring.splitlines(keepends=True)
                    )
                    i = 0
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() == "C" or atom.GetSymbol() == "H":
                            value2 = self.gen.get_Hose_codes(
                                mol, i, usestereo=True, wedgebond=wedgemap
                            )
                        i += 1
                molstring = ""
            else:
                molstring += line
        file1.close()


if __name__ == "__main__":
    import xmlrunner

    unittest.main(
        testRunner=xmlrunner.XMLTestRunner(output="test-reports"),
        failfast=False,
        buffer=False,
        catchbreak=False,
    )
