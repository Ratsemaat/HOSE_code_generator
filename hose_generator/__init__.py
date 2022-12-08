import sys

from rdkit import Chem
from helpers.Node import Node


class HoseGenerator():
    table = Chem.rdchem.GetPeriodicTable()

    def __init__(self):
        self.center_code = None
        self.symbol_rankings = [9000, 8900, 8800, 8700, 8600, 8500, 8400, 8300, 8200, 8100,
                                8000, 7900, 1200, 1100, 1000]

        self.ranked_symbols = ["C", "O", "N", "S", "P", "Si", "B", "F", "Cl", "Br", ";", "I",
                               "#", "&", ","]

        self.sphere_delimiters = ["(", "/", "/", ")", "/", "/", "/", "/", "/", "/", "/", "/"]

        self.spheres = []
        self.mol = None
        self.visited_atoms = []
        self.max_radius = 0
        self.HOSE_code = ""

    def create_code(self):
        self.visited_atoms = [False for i in range(len(self.mol.GetAtoms()))]
        for sphere in self.spheres:
            for node in sphere:
                if node.connector:
                    node.connector.ranking += node.degree

            for node_set in self.spheres:
                self.calculate_node_scores(node_set)

            for node_set in self.spheres:
                for node in node_set:
                    node.score += node.ranking
                    node.stringscore = (str(node.connector.stringscore)+ str(node.score).zfill(6))
                node_set.sort(key=lambda x: x.stringscore, reverse=True)

        self.HOSE_code += self.center_code
        for (idx, sphere) in enumerate(self.spheres):
            self.HOSE_code += self.get_sphere_code(idx, sphere)

    def get_Hose_codes(self,
                       mol: Chem.rdchem.Mol,
                       atom_idx: int,
                       min_radius: int = 2,
                       max_radius: int = 5):
        self.spheres = []
        self.mol = mol
        self.max_radius = max_radius
        self.visited_atoms = [False for i in range(len(self.mol.GetAtoms()))]
        self.visited_atoms[atom_idx] = True

        atom = mol.GetAtoms()[atom_idx]
        root = Node(self.get_element_symbol(atom.GetAtomicNum()), None, atom, 0, atom.GetDegree() + atom.GetTotalNumHs(), 0)
        self.HOSE_code = ""
        self.create_center_code(atom)
        current_layer = [root]

        for i in range(max_radius):
            next_layer = []

            for node in current_layer:
                if node.symbol in "H&;#:,":
                    continue

                bonds = node.atom.GetBonds()
                if len(bonds) == 1 and i>0:
                    next_layer.append(Node(",", node, None, 0, 0, node.score))
                else:
                    for bond in bonds:
                        to_atom = bond.GetOtherAtom(node.atom)
                        if not node.connector or to_atom.GetIdx() != node.connector.atom.GetIdx():
                            next_layer.append(Node(self.get_element_symbol(to_atom.GetAtomicNum()), node, to_atom,
                                                   bond.GetBondType(), to_atom.GetDegree() + to_atom.GetTotalNumHs(),
                                                   node.score))
            next_layer.sort(key=lambda x: x.atom.GetIdx() if x.atom else -sys.maxsize - 1, reverse=True)
            self.spheres.append(next_layer)
            current_layer = next_layer

        self.create_code()
        self.fill_up_sphere_delimiters()
        return self.HOSE_code

    def calculate_node_scores(self, sphere_nodes):
        for node in sphere_nodes:
            node.score += self.get_element_rank(node.symbol)
            node.score += self.get_bond_score(node.bond_type)

    def get_bond_score(self, bond_type):
        if bond_type == Chem.BondType.AROMATIC:
            return 100000
        if bond_type == Chem.BondType.SINGLE:
            return 0
        if bond_type == Chem.BondType.DOUBLE:
            return 200000
        if bond_type == Chem.BondType.TRIPLE:
            return 300000
        return 0

    def get_bond_symbol(self, bond_type):
        if bond_type == Chem.BondType.AROMATIC:
            return "*"
        if bond_type == Chem.BondType.SINGLE:
            return ""
        if bond_type == Chem.BondType.DOUBLE:
            return "="
        if bond_type == Chem.BondType.TRIPLE:
            return "%"
        return ""


    def get_element_rank(self, symbol):
        loc = self.ranked_symbols.index(symbol)
        if loc >= 0:
            return self.symbol_rankings[loc]
        if self.table.GetMostCommonIsotopeMass(symbol):
            return 80000 - self.table.GetMostCommonIsotopeMass(symbol)
        return 80000

    def get_sphere_code(self, radius, sphere_nodes):
        if not sphere_nodes or len(sphere_nodes) < 1:
            return self.sphere_delimiters[radius]

        branch = sphere_nodes[0].connector.atom
        code = ""
        for node in sphere_nodes:
            temp_code = ""

            if not node.connector.stopper and node.connector.atom.GetIdx() != branch.GetIdx():
                branch = node.connector.atom
                code += ","

            if not node.connector.stopper and node.connector.atom.GetIdx() == branch.GetIdx():
                temp_code += self.get_bond_symbol(node.bond_type)
                if node.atom and not self.visited_atoms[node.atom.GetIdx()]:
                    temp_code += node.symbol
                elif node.atom and self.visited_atoms[node.atom.GetIdx()]:
                    temp_code += '&'
                    node.stopper = True

                code += temp_code + self.create_charge_code(node.atom)
                # node.h

            if node.atom:
                self.visited_atoms[node.atom.GetIdx()] = True
            if node.connector.stopper:
                node.stopper = True
        code += self.sphere_delimiters[radius]
        return code

    def create_charge_code(self, atom):
        temp_code = ""
        if atom:
            formal_charge = atom.GetFormalCharge()

            if formal_charge != 0:
                if abs(formal_charge) == 1:
                    if formal_charge < 0:
                        temp_code.append('-')
                    else:
                        temp_code.append('+')
                else:
                    temp_code.append('\'')
                    if formal_charge > 0:
                        temp_code.append('+')
                    temp_code.append(formal_charge).append('\'')

        return temp_code

    def create_center_code(self, atom):
        partners_count = atom.GetTotalNumHs()+atom.GetDegree()
        self.center_code = self.table.GetElementSymbol(atom.GetAtomicNum()) + "-" + str(partners_count) \
                           + self.create_charge_code(atom) + ";"

    def get_element_symbol(self, param):
        return self.table.GetElementSymbol(param)

    def fill_up_sphere_delimiters(self):
        for i in range(len(self.spheres), 4):
            self.HOSE_code.append(self.sphere_delimiters[i]);
