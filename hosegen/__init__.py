import sys

from rdkit import Chem
from rdkit.Chem.rdchem import BondStereo
from rdkit.Chem.rdchem import BondType
from rdkit.Chem.rdchem import BondDir
from hosegen.helpers import Node
from hosegen.geometry import *

class HoseGenerator():
    table = Chem.rdchem.GetPeriodicTable()

    def __init__(self):
        self.rec=0
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
        self.rec += 1
        reorders=set()
        for sphere in self.spheres:
            for node in sphere:
                if node.connector:
                    node.connector.ranking += node.degree
                if node.reorder is not None:
                    reorders.add(tuple(node.reorder))

            for node_set in self.spheres:
                self.calculate_node_scores(node_set)

        for node_set in self.spheres:
            for reorder in reorders:
                toreorders=np.empty(len(reorder), dtype=object) 
                count=0
                for node2 in node_set:
                    if node2.reorder is not None and tuple(node2.reorder)==reorder:
                        toreorders[count]=node2
                        count=count+1
                count=0
                for node2 in node_set:
                    #that is the core of the stereochemistry handling, the scores/rankings are swapped between atoms to align with reorder
                    if node2.reorder is not None and tuple(node2.reorder)==reorder:
                        for idx, node3 in enumerate(reorder):
                            if node2.atom.GetIdx()==node3.atom.GetIdx():
                                node2.score= reorder[count].score+len(reorder)-idx   #-idx is needed to distinguish identical atoms with different wedges
                                node2.ranking= reorder[count].ranking
                        count=count+1

        for sphere in self.spheres:
            for node_set in self.spheres:
                for node in node_set:
                    node.score += node.ranking
                    node.stringscore = (str(node.connector.stringscore)+ str(node.score).zfill(6))
                node_set.sort(key=lambda x: x.stringscore, reverse=True)

        self.HOSE_code += self.center_code
        for (idx, sphere) in enumerate(self.spheres):
            self.HOSE_code += self.get_sphere_code(idx, sphere)

    def get_Hose_codes_from_file(self,
                       molfilepath: str,
                       atom_idx: int,
                       max_radius: int = 5,
                       usestereo: bool = False,
                       strict: bool = False,
                       ringsize: bool = False):
        if usestereo:
            mol = Chem.MolFromMolFile(molfilepath)
            #we add explicit hydrogens and do the up/down bonds
            mol = Chem.rdmolops.AddHs(mol)
            wedgebonds = create_wedgemap(molfilepath)
            for bond in mol.GetBonds():
                if not bond.GetIdx() in wedgebonds:
                    wedgebonds[bond.GetIdx()]=0
            makeUpDownBonds(mol, wedgebonds)
            return self.get_Hose_codes(mol, atom_idx, max_radius, usestereo, wedgebond=wedgebonds, strict=strict, ringsize=ringsize)
        else:
            return self.get_Hose_Codes(mol, atom_idx, max_radius, usestereo, strict=strict, ringsize=ringsize)


    def get_Hose_codes(self,
                       mol: Chem.rdchem.Mol,
                       atom_idx: int,
                       max_radius: int = 5,
                       usestereo: bool = False,
                       wedgebond: dict = None,
                       strict: bool = False,
                       ringsize: bool = False):
        #we use e/z info from rdkit, but we find all 4 atoms participating, to make it easier later
        self.zes=[]
        self.ezs=[]
        if usestereo:
            for bnd in mol.GetBonds():
                if bnd.GetStereo()==BondStereo.STEREOE:
                    beginNbrs = [x.GetIdx() for x in bnd.GetBeginAtom().GetNeighbors() if x.GetIdx()!=bnd.GetEndAtomIdx()]
                    if len(beginNbrs)==1:
                        beginNbrs.append(-1)
                    beginNbrs.extend([x.GetIdx() for x in bnd.GetEndAtom().GetNeighbors() if x.GetIdx()!=bnd.GetBeginAtomIdx()])
                    if len(beginNbrs)==3:
                        beginNbrs.append(-1)
                    self.ezs.append(beginNbrs)
                if bnd.GetStereo()==BondStereo.STEREOZ:
                    beginNbrs = [x.GetIdx() for x in bnd.GetBeginAtom().GetNeighbors() if x.GetIdx()!=bnd.GetEndAtomIdx()]
                    if len(beginNbrs)==1:
                        beginNbrs.append(-1)
                    beginNbrs.extend([x.GetIdx() for x in bnd.GetEndAtom().GetNeighbors() if x.GetIdx()!=bnd.GetBeginAtomIdx()])
                    if len(beginNbrs)==3:
                        beginNbrs.append(-1)
                    self.zes.append(beginNbrs)

        self.spheres = []
        self.mol = mol
        self.max_radius = max_radius
        self.visited_atoms = [False for i in range(len(self.mol.GetAtoms()))]
        self.visited_atoms[atom_idx] = True

        atom = mol.GetAtoms()[atom_idx]
        root = Node(self.get_element_symbol(atom.GetAtomicNum()), None, atom, 0, atom.GetDegree() + atom.GetTotalNumHs(), 0)
        self.HOSE_code = ""
        self.create_center_code(atom, mol, ringsize)
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
                    interim_list = []
                    for bond in bonds:
                        to_atom = bond.GetOtherAtom(node.atom)
                        if not node.connector or to_atom.GetIdx() != node.connector.atom.GetIdx():
                            interim_list.append(Node(self.get_element_symbol(to_atom.GetAtomicNum()), node, to_atom,
                                                   bond.GetBondType(), to_atom.GetDegree() + to_atom.GetTotalNumHs(),
                                                   node.score))
                    #This checks for stereochemistry and remembers it as reordered lists in nodes
                    if usestereo and is_chiral_centre(node.atom, mol, wedgebond):
                        reordered = order_stereo(self, interim_list, mol, node.atom, None if node.connector is None else node.connector, wedgebond, strict)
                        for node in interim_list:
                            node.reorder=reordered
                    next_layer.extend(interim_list)
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

        if symbol in self.ranked_symbols:
            loc = self.ranked_symbols.index(symbol)
            return self.symbol_rankings[loc]

        if self.table.GetMostCommonIsotopeMass(symbol):
            return 80000 - self.table.GetMostCommonIsotopeMass(symbol)

        return 80000

    def get_sphere_code(self, radius, sphere_nodes):
        if not sphere_nodes or len(sphere_nodes) < 1:
            return self.sphere_delimiters[radius]

        branch = sphere_nodes[0].connector.atom
        code = ""
        #do we need @?
        needat={}
        count=0
        needat[0]=False
        for node in sphere_nodes:
            if node.reorder is not None:
                needat[count]=True
            if node.connector.stopper:
                node.stopper = True
                count=count+1
                needat[count]=False


        count=0
        if needat[0]:
            code += "@"
        for node in sphere_nodes:
            temp_code = ""

            if not node.connector.stopper and node.connector.atom.GetIdx() != branch.GetIdx():
                branch = node.connector.atom
                code += ","

            if not node.connector.stopper and node.connector.atom.GetIdx() == branch.GetIdx():
                temp_code += self.get_bond_symbol(node.bond_type)

                #if an atom is part of an e/z config, we set \ or | accordingly (order of atoms in hose code is not changed)
                if node.atom is not None and node.connector is not None and node.connector.connector is not None\
                and node.connector.connector.connector is not None\
                and (node.connector.atom.GetSymbol()=="C" or node.connector.atom.GetSymbol()=="N")\
                and node.connector.bond_type==BondType.DOUBLE \
                and (node.connector.connector.atom.GetSymbol()=="C" or node.connector.connector.atom.GetSymbol()=="N"):
                    if any(x for x in self.ezs if (x[0]==node.atom.GetIdx() and x[2]==node.connector.connector.connector.atom.GetIdx())) or any(x for x in self.ezs if (x[1]==node.atom.GetIdx() and x[3]==node.connector.connector.connector.atom.GetIdx())):
                        code+="\\"
                    if any(x for x in self.zes if (x[0]==node.atom.GetIdx() and x[2]==node.connector.connector.connector.atom.GetIdx())) or any(x for x in self.zes if (x[1]==node.atom.GetIdx() and x[3]==node.connector.connector.connector.atom.GetIdx())):
                        code+="|"
                    if any(x for x in self.ezs if (x[0]==node.connector.connector.connector.atom.GetIdx() and x[2]==node.atom.GetIdx())) or any(x for x in self.ezs if (x[0]==node.connector.connector.connector.atom.GetIdx() and x[3]==node.atom.GetIdx())):
                        code+="\\"
                    if any(x for x in self.zes if (x[0]==node.connector.connector.connector.atom.GetIdx() and x[2]==node.atom.GetIdx())) or any(x for x in self.zes if (x[0]==node.connector.connector.connector.atom.GetIdx() and x[3]==node.atom.GetIdx())):
                        code+="|"

                if node.atom and not self.visited_atoms[node.atom.GetIdx()]:
                    temp_code += node.symbol
                elif node.atom and self.visited_atoms[node.atom.GetIdx()]:
                    temp_code += '&'
                    node.stopper = True

                code += temp_code + self.create_charge_code(node.atom)
            if node.atom:
                self.visited_atoms[node.atom.GetIdx()] = True
            if node.connector.stopper:
                node.stopper = True
                count=count+1
                if needat[count]:
                    code += "@"
        code += self.sphere_delimiters[radius]
        return code

    def create_charge_code(self, atom):
        temp_code = ""
        if atom:
            formal_charge = atom.GetFormalCharge()

            if formal_charge != 0:
                if abs(formal_charge) == 1:
                    if formal_charge < 0:
                        temp_code += '-'
                    else:
                        temp_code += '+'
                else:
                    temp_code += '\''
                    if formal_charge > 0:
                        temp_code += '+'
                    temp_code += str(formal_charge) + '\''

        return temp_code

    def create_center_code(self, atom, mol, ringsize):
        partners_count = atom.GetTotalNumHs()+atom.GetDegree()
        self.center_code = self.table.GetElementSymbol(atom.GetAtomicNum()) + "-" + str(partners_count) \
            + self.create_charge_code(atom) + (self.get_ringcode(atom, mol) * ringsize) + ";"

    def get_element_symbol(self, param):
        return self.table.GetElementSymbol(param)

    def fill_up_sphere_delimiters(self):
        for i in range(len(self.spheres), 4):
            self.HOSE_code+=(self.sphere_delimiters[i])

    def get_ringcode(self, root, mol):
        ri = mol.GetRingInfo()
        minsize = ri.MinAtomRingSize(root.GetIdx())
        if minsize>0:
            return "-"+str(minsize)
        else:
            return ""
