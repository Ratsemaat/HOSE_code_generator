import sys
from typing import Optional

from rdkit import Chem
from rdkit.Chem.rdchem import BondStereo, BondType, Atom, Mol, PeriodicTable, Bond
from hosegen.helpers import Node
from hosegen.exceptions import InvalidStartingAtomIndex
from hosegen.geometry import *

import numpy as np

class HoseGenerator():

    table: PeriodicTable = Chem.rdchem.GetPeriodicTable()

    def __init__(self) -> None:
        self.rec=0
        self.center_code: str =''
        self.symbol_rankings: dict[str, int] = {"C": 9000, "O": 8900, "N": 8800,"S": 8700, "P":8600, "Si":8500,  "B": 8400,  "F": 8300, "Cl": 8200, 'Br': 8100,
                                ";": 8000, "I": 7900,   "#": 1200, "&": 1100,  ",": 1000}

        self.sphere_delimiters: list[str] = ["(", "/", "/", ")", "/", "/", "/", "/", "/", "/", "/", "/"]

        self.spheres: list[list[Node]] = []
        self.visited_atoms: list[bool]
        self.HOSE_code: str = ""
        


    def create_code(self) -> None:
        self.visited_atoms = [False for i in range(len(self.mol.GetAtoms()))] #type: ignore
        self.rec += 1
        reorders: set[tuple[Node]] = set()
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
                toreorders = np.empty(len(reorder), dtype=object)
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
                            if not isinstance(node2.atom, Atom) or not isinstance(node3.atom, Atom):
                                raise Exception("Bad state")
                            if node2.atom.GetIdx()==node3.atom.GetIdx():
                                node2.score= reorder[count].score+len(reorder)-idx   #-idx is needed to distinguish identical atoms with different wedges
                                node2.ranking= reorder[count].ranking
                        count=count+1

        for sphere in self.spheres:
            for node_set in self.spheres:
                for node in node_set:
                    node.score += node.ranking
                    if not isinstance(node.connector, Node):
                            raise Exception("Bad state")

                    node.stringscore = str(node.connector.stringscore) + str(node.score).zfill(6)
                node_set.sort(key=lambda x: x.stringscore, reverse=True)

        atoms: list[Atom] = self.mol.GetAtoms() #type:ignore
        for atom in atoms:
            atom.SetProp('visited', 'false')

        self.HOSE_code += self.center_code
        for (idx, sphere) in enumerate(self.spheres):
            self.HOSE_code += self.get_sphere_code(idx, sphere)

    def get_Hose_codes_from_file(self,
                       molfilepath: str,
                       atom_idx: int,
                       max_radius: int = 5,
                       usestereo: bool = False,
                       strict: bool = False,
                       ringsize: bool = False) -> str:

        self.mol: Mol = Chem.MolFromMolFile(molfilepath)
        if usestereo:
            #we add explicit hydrogens and do the up/down bonds
            self.mol = Chem.rdmolops.AddHs(self.mol)
            wedgebonds = create_wedgemap(molfilepath) #type:ignore
            bonds: list[Bond] =  self.mol.GetBonds() #type: ignore
            for bond in bonds:
                if not bond.GetIdx() in wedgebonds:
                    wedgebonds[bond.GetIdx()]=0
            makeUpDownBonds(self.mol, wedgebonds) #type: ignore
            return self.get_Hose_codes(self.mol, atom_idx, max_radius, usestereo, wedgebond=wedgebonds, strict=strict, ringsize=ringsize)
        else:
            return self.get_Hose_codes(self.mol, atom_idx, max_radius, usestereo, strict=strict, ringsize=ringsize)


    def get_Hose_codes(self,
                       mol: Chem.rdchem.Mol,
                       atom_idx: int,
                       max_radius: int = 5,
                       usestereo: bool = False,
                       wedgebond: dict[int, int] = {},
                       strict: bool = False,
                       ringsize: bool = False) -> str:
        #we use e/z info from rdkit, but we find all 4 atoms participating, to make it easier later
        if len(mol.GetAtoms()) <= atom_idx: #type: ignore
            raise InvalidStartingAtomIndex(f"atom_idx '{atom_idx}'  is out of range. Total explicit atoms: {mol.GetAtoms()}")#type: ignore
        self.zes=[]
        self.ezs=[]
        if usestereo:
            bonds: tuple[Bond] = mol.GetBonds()#type: ignore
            for bnd in bonds:
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
        self.visited_atoms = [False for i in range(len(self.mol.GetAtoms()))] #type: ignore
        self.visited_atoms[atom_idx] = True

        atom: Atom = mol.GetAtoms()[atom_idx] #type: ignore
        root = Node(self.get_element_symbol(atom.GetAtomicNum()), None, atom, None, atom.GetDegree() + atom.GetTotalNumHs(), 0)
        self.HOSE_code = ""
        self.create_center_code(atom, mol, ringsize)
        current_layer: list[Node] = [root]


        for i in range(max_radius):
            next_layer: list[Node] = []

            for node in current_layer:
                if node.symbol in "H&;#:," and i > 0:
                    continue
                
                if not isinstance(node.atom, Atom):
                    raise Exception('Bad State')
                bonds = node.atom.GetBonds()
                if len(bonds) == 1 and i>0:
                    next_layer.append(Node(",", node, None, None, 0, node.score))
                else:
                    interim_list = []
                    for bond in bonds:
                        to_atom = bond.GetOtherAtom(node.atom)
                        if not isinstance(node.connector, Node) or isinstance(node.connector.atom, Atom) and to_atom.GetIdx() != node.connector.atom.GetIdx():
                            interim_list.append(Node(self.get_element_symbol(to_atom.GetAtomicNum()), node, to_atom,
                                                   bond.GetBondType(), to_atom.GetDegree() + to_atom.GetTotalNumHs(),
                                                   node.score))
                    #This checks for stereochemistry and remembers it as reordered lists in nodes
                    if usestereo and is_chiral_centre(node.atom, mol, wedgebond): #type: ignore
                        reordered = order_stereo(self, interim_list, mol, node.atom, None if node.connector is None else node.connector, wedgebond, strict) #type: ignore
                        for node in interim_list:
                            node.reorder=reordered
                    next_layer.extend(interim_list)
            next_layer.sort(key=lambda x: x.atom.GetIdx() if x.atom else -sys.maxsize - 1, reverse=True)
            self.spheres.append(next_layer)
            current_layer = next_layer

        self.create_code()
        self.fill_up_sphere_delimiters()
        return self.HOSE_code

    def calculate_node_scores(self, sphere_nodes: list[Node]) -> None:
        visited: set[Node] = set()

        for node in sphere_nodes:

            if isinstance(node, Node) and isinstance(node.atom, Atom) and node.atom.HasProp('visited') and node.atom.GetProp('visited') == 'true':
                node.score += self.get_element_rank("&")
            else:
                node.score += self.get_element_rank(node.symbol)
            node.score += self.get_bond_score(node.bond_type)
            visited.add(node)
        
        for node in visited:
            if isinstance(node.atom, Atom) :
                node.atom.SetProp('visited','true')

    def get_bond_score(self, bond_type: Optional[BondType]) -> int:
        if bond_type == BondType.AROMATIC:
            return 100000
        if bond_type == BondType.SINGLE:
            return 0
        if bond_type == BondType.ZERO:
            return 400000
        if bond_type == BondType.DOUBLE:
            return 200000
        if bond_type == BondType.TRIPLE:
            return 300000
        return 0

    def get_bond_symbol(self, bond_type: Optional[BondType]) -> str:
        if bond_type == BondType.AROMATIC:
            return "*"
        if bond_type == BondType.SINGLE:
            return ""
        if bond_type == BondType.ZERO:
            return "<"
        if bond_type == BondType.DOUBLE:
            return "="
        if bond_type == BondType.TRIPLE:
            return "%"
        return ""


    def get_element_rank(self, symbol:str) -> float:

        if symbol in self.symbol_rankings:
            return self.symbol_rankings[symbol]

        if HoseGenerator.table.GetMostCommonIsotopeMass(symbol):
            return 80000 - HoseGenerator.table.GetMostCommonIsotopeMass(symbol)

        return 80000


    def get_sphere_code(self, radius: int, sphere_nodes: list[Node]) -> str:
        if len(sphere_nodes) == 0:
            return self.sphere_delimiters[radius]
        
        if not isinstance(sphere_nodes[0].connector, Node):
            branch = None
        else:
            branch = sphere_nodes[0].connector.atom
        code = ""
        #do we need @?
        needat={}
        count=0
        needat[0]=False
        for node in sphere_nodes:
            if node.reorder is not None:
                needat[count]=True
            if isinstance(node.connector, Node) and node.connector.stopper:
                node.stopper = True
                count=count+1
                needat[count]=False


        count=0
        if needat[0]:
            code += "@"
        for node in sphere_nodes:
            temp_code = ""

            if (isinstance(node.connector, Node) and not node.connector.stopper):

                if (branch is None or node.connector is None\
                        or node.connector.atom is None \
                        or node.connector.atom.GetIdx() != branch.GetIdx()):
                    branch = node.connector.atom
                    code += ","

                temp_code += self.get_bond_symbol(node.bond_type)

                #if an atom is part of an e/z config, we set \ or | accordingly (order of atoms in hose code is not changed)
                if isinstance(node.atom, Atom)\
                        and isinstance(node.connector, Node) \
                        and isinstance(node.connector.atom, Atom) \
                        and isinstance(node.connector.connector, Node)\
                        and isinstance(node.connector.connector.atom, Atom)\
                        and isinstance(node.connector.connector.connector, Node)\
                        and isinstance(node.connector.connector.connector.atom, Atom)\
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

                if isinstance(node.atom, Atom):
                    code += temp_code + self.create_charge_code(node.atom)
            if node.atom:
                self.visited_atoms[node.atom.GetIdx()] = True
            if isinstance(node.connector, Node) and node.connector.stopper:
                node.stopper = True
                count=count+1
                if needat[count]:
                    code += "@"
        code += self.sphere_delimiters[radius]
        return code

    def create_charge_code(self, atom: Atom) -> str:
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

    def create_center_code(self, atom: Atom, mol: Mol, ringsize: int) -> None:
        partners_count = atom.GetTotalNumHs()+atom.GetDegree()
        self.center_code = self.table.GetElementSymbol(atom.GetAtomicNum()) + "-" + str(partners_count) \
            + self.create_charge_code(atom) + (self.get_ringcode(atom, mol) * ringsize) + ";"

    def get_element_symbol(self, param:int) -> str:
        return HoseGenerator.table.GetElementSymbol(param)

    def fill_up_sphere_delimiters(self) -> None:
        for i in range(len(self.spheres), 4):
            self.HOSE_code+=self.sphere_delimiters[i]

    def get_ringcode(self, root: Atom, mol: Mol) -> str:
        ri = mol.GetRingInfo()
        minsize = ri.MinAtomRingSize(root.GetIdx())
        if minsize>0:
            return "-"+str(minsize)
        else:
            return ""
