class Node:
    atom = None
    connector = None
    score = 0
    degree = 0
    bond_type = ""
    ranking = 0
    symbol = ""
    stringscore=""
    stopper = False
    reorder = None

    def __init__(self, symbol, connector, atom, bond_type, degree, score):
        self.symbol = symbol
        self.connector = connector
        self.atom = atom
        self.degree = degree
        self.bond_type = bond_type
        self.score = score
        ranking = 0
        sortOrder = 1
