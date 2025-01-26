from rdkit.Chem.rdchem import Atom, BondType
from typing import Optional

class Node:

    def __init__(self, symbol:str, connector: Optional['Node'], atom: Optional[Atom], bond_type:Optional[BondType], degree:int, score:float):
        self.symbol = symbol
        self.connector = connector
        self.atom = atom
        self.degree = degree
        self.bond_type = bond_type
        self.score = score
        self.ranking = 0
        self.sortOrder = 1
        self.reorder= None
        self.stringscore: str = ''
        self.stopper: bool = False
