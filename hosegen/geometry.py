import numpy as np
import collections
from hosegen.helpers import Node
import math

def is_chiral_centre(atom, mol, wedgemap):
    atoms = atom.GetNeighbors()
    if len(atoms)!=4 and len(atoms)!=6:
        return False
    for neighbour in atoms:
        if mol.GetBondBetweenAtoms(atom.GetIdx(), neighbour.GetIdx()).GetBondTypeAsDouble()>1:
            return False
    for bond in atom.GetBonds():
        if wedgemap[bond.GetIdx()]==6:
            return True
        if wedgemap[bond.GetIdx()]==1:
            return True
    return False

def is_tetrahedral(mol, atom, wedgemap, strict=False):
    atoms = atom.GetNeighbors()
    if len(atoms) != 4:
        return 0
    bonds = atom.GetBonds()
    up = 0
    down = 0
    for bond in bonds:
        if wedgemap[bond.GetIdx()] == 1 and bond.GetBeginAtom().GetIdx()==atom.GetIdx():
            up+=1
        elif wedgemap[bond.GetIdx()] == 6 and bond.GetBeginAtom().GetIdx()==atom.GetIdx():
            down+=1
    if up == 1 and down == 1:
        return 1
    if up == 2 and down == 2:
        if stereos_are_opposite(mol, atom, wedgemap):
            return 2
        else:
            return 0
    if up == 1 and down == 0 and not strict:
        return 3
    if down == 1 and up == 0 and not strict:
        return 4
    if down == 2 and up == 1 and not strict:
        return 5
    if down == 1 and up == 2 and not strict:
        return 6
    return 0


def order_stereo(gen, nodes, mol, chiralCenter, parent, wedgemap, strict=False):
    chiralNeighbours = chiralCenter.GetNeighbors()
    sorted = np.empty(3, dtype=object) 
    sortednodes = np.empty(3, dtype=object)
    #this is in case the start atom is a chiral centre
    sorted = np.empty(3, dtype=object) 
    if(parent is None):
        parent=Node(gen.get_element_symbol(chiralNeighbours[0].GetAtomicNum()), None, chiralNeighbours[0], 0, chiralNeighbours[0].GetDegree() + chiralNeighbours[0].GetTotalNumHs(), 0)
        sorted = np.empty(4, dtype=object) 
        sortednodes = np.empty(4, dtype=object)
    if is_tetrahedral(mol, chiralCenter, wedgemap, strict) == 1:
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
            first=True
            firstangle=0.0
            firstatom=0
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0 or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                        angle = give_angle(mol, chiralCenter, neighbour, parent.atom)
                        if not first:
                            if angle<firstangle:
                                sorted[2]=neighbour
                                sorted[1]=firstatom
                            else:
                                sorted[2]=neighbour
                                sorted[1]=firstatom
                        else:
                            firstangle=angle
                            first=False
                            firstatom=neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1:
                        sorted[0] = neighbour
        elif wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
            first=True
            firstangle=0.0
            firstatom=0
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0 or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                        angle = give_angle(mol, chiralCenter, neighbour, parent.atom)
                        if not first:
                            if angle<firstangle:
                                sorted[2]=neighbour
                                sorted[1]=firstatom
                            else:
                                sorted[1]=neighbour
                                sorted[2]=firstatom
                        else:
                            firstangle=angle
                            first=False
                            firstatom=neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6:
                        sorted[0] = neighbour
        else: #if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0:
            normalBindingIsLeft = False
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0 or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                        if is_left(mol, neighbour, parent.atom, chiralCenter):
                            normalBindingIsLeft = True
                            break
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if not normalBindingIsLeft:
                        if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0 or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                            sorted[0] = neighbour
                        if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                            sorted[1] = neighbour
                        if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                            sorted[2] = neighbour
                    else:
                        if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                            sorted[2] = neighbour
                        if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0 or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                            sorted[0] = neighbour
                        if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                            sorted[1] = neighbour
    if is_tetrahedral(mol, chiralCenter, wedgemap, strict) == 2:
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx())]==1 and mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6  and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                        sorted[2] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6  and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                        sorted[1] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1  and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                        sorted[0] = neighbour
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx())]==6 and mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
            angle1 = 0.0
            angle2 = 0.0
            atom1 = None
            atom2 = None
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                        if angle1 == 0:
                            angle1 = give_angle(mol, chiralCenter, parent.atom, neighbour)
                            atom1 = neighbour
                        else:
                            angle2 = give_angle(mol, chiralCenter, parent.atom, neighbour)
                            atom2 = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                            sorted[1] = neighbour
            if angle1 < angle2:
                sorted[0] = atom1
                sorted[2] = atom2
            else:
                sorted[0] = atom2
                sorted[2] = atom1
    if is_tetrahedral(mol, chiralCenter, wedgemap, strict) == 3:
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
            hm={}
            for idx,neighbour in enumerate(chiralNeighbours):
                if neighbour.GetIdx() != parent.atom.GetIdx():
                    hm[give_angle(mol, chiralCenter, parent.atom, neighbour)]= idx
            sorted_dict = collections.OrderedDict(hm)
            for idx, value in enumerate(hm.keys()):
                sorted[idx]=chiralNeighbours[sorted_dict[value]]
        elif wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0 or mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
            angle1 = 0.0
            angle2 = 0.0
            atom1 = None
            atom2 = None
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx():
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0  or wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==3  or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                        if angle1 == 0:
                            angle1 = give_angle(mol, chiralCenter, parent.atom, neighbour)
                            if angle1==0:
                                angle1=0.1
                            atom1 = neighbour
                        else:
                            angle2 = give_angle(mol, chiralCenter, parent.atom, neighbour)
                            atom2 = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and  mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                        sorted[0] = neighbour
            if angle1 < angle2:
                sorted[1] = atom2
                sorted[2] = atom1
            else:
                sorted[1] = atom1
                sorted[2] = atom2
    if is_tetrahedral(mol, chiralCenter, wedgemap, strict) == 4:
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6  and mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
            hm = {}
            idx=0
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx():
                    hm[give_angle(mol, chiralCenter, parent.atom, neighbour)]=idx
                    idx+=1
            sorted_dict = collections.OrderedDict(hm)
            for idx, value in enumerate(hm.keys()):
                sorted[idx]=chiralNeighbours[sorted_dict[value]]
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0  or mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
            angle1 = 0.0
            angle2 = 0.0
            atom1 = None
            atom2 = None
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0  or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                        if angle1 == 0:
                            angle1 = give_angle(mol, chiralCenter, parent.atom, neighbour)
                            atom1 = neighbour
                        else:
                            angle2 = give_angle(mol, chiralCenter, parent.atom, neighbour)
                            atom2 = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                        sorted[2] = neighbour
            if angle1 < angle2:
                sorted[1] = atom2
                sorted[0] = atom1
            else:
                sorted[1] = atom1
                sorted[0] = atom2
    if is_tetrahedral(mol, chiralCenter, wedgemap, strict) == 5:
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6  and mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                        sorted[0] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0  or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                        sorted[1] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                        sorted[2] = neighbour
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx() and is_left(mol, neighbour, parent.atom, chiralCenter):
                        sorted[0] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx() and not is_left(mol, neighbour, parent.atom, chiralCenter):
                        sorted[1] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0 or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                        sorted[2] = neighbour
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0 or mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx() and is_left(mol, neighbour, parent.atom, chiralCenter):
                        sorted[0] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx() and not is_left(mol, neighbour, parent.atom, chiralCenter):
                        sorted[1] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                        sorted[2] = neighbour
    if is_tetrahedral(mol, chiralCenter, wedgemap, strict) == 6:
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                        sorted[0] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0 or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                        sorted[1] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                        sorted[2] = neighbour
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1  and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx() and is_left(mol, neighbour, parent.atom, chiralCenter):
                        sorted[1] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1  and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx()and not is_left(mol, neighbour, parent.atom, chiralCenter):
                        sorted[0] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0 or mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
                        sorted[2] = neighbour
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0 or mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()!=chiralCenter.GetIdx():
            for neighbour in chiralNeighbours:
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1  and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx() and is_left(mol, neighbour, parent.atom, chiralCenter):
                        sorted[1] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx() and not is_left(mol, neighbour, parent.atom, chiralCenter):
                        sorted[0] = neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6 and mol.GetBondBetweenAtoms(neighbour.GetIdx(), chiralCenter.GetIdx()).GetBeginAtom().GetIdx()==chiralCenter.GetIdx():
                        sorted[2] = neighbour
    #this means start was chiral center, parent is part of the sphere
    if len(sorted)==4 and is_tetrahedral(mol, chiralCenter, wedgemap, strict)>0:
        sorted[3]=sorted[2]
        sorted[2]=sorted[1]
        sorted[1]=sorted[0]
        sorted[0]=parent.atom
    if is_square_planar(mol, chiralCenter, wedgemap):
        #This produces a U=SP1 order in every case
        hm = {}
        for idx,neighbour in enumerate(chiralNeighbours):
            if neighbour.GetIdx() != parent.atom.GetIdx():
                hm[give_angle(mol, chiralCenter, parent.atom, neighbour)]=idx
        sorted_dict = collections.OrderedDict(hm)
        for idx, value in enumerate(sorted_dict.keys()):
            sorted[idx]=chiralNeighbours[hm[value]]
        if len(sorted)==4:
            sorted[3]=parent.atom
    if is_trigonal_bipyramidal_or_octahedral(mol, chiralCenter, wedgemap)!=0:
        sorted = np.empty(len(chiralCenter.GetBonds())-1, dtype=object)
        hm = {}
        numbyangle=len(chiralCenter.GetBonds())-2
        sortednodes = np.empty(len(chiralCenter.GetBonds())-1, dtype=object)
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==1:
            for idx,neighbour in enumerate(chiralNeighbours):
                if wedgemap[mol.GetBondBetweenAtoms(chiralCenter.GetIdx(), neighbour.GetIdx()).GetIdx()]==0:
                    hm[give_angle(mol, chiralCenter, parent.atom, neighbour)]=idx
                if wedgemap[mol.GetBondBetweenAtoms(chiralCenter.GetIdx(), neighbour.GetIdx()).GetIdx()]==6:
                    sorted[len(sorted)-1]=neighbour
            sorted_dict = collections.OrderedDict(hm)
            for idx, value in enumerate(sorted_dict.keys()):
                sorted[idx]=chiralNeighbours[hm[value]]
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==6:
            for idx,neighbour in enumerate(chiralNeighbours):
                if wedgemap[mol.GetBondBetweenAtoms(chiralCenter.GetIdx(), neighbour.GetIdx()).GetIdx()]==0:
                    hm[give_angle(mol, chiralCenter, parent.atom, neighbour)]=idx
                if wedgemap[mol.GetBondBetweenAtoms(chiralCenter.GetIdx(), neighbour.GetIdx()).GetIdx()]==1:
                    sorted[len(sorted)-1]=neighbour
            sorted_dict = collections.OrderedDict(hm)
            for idx, value in enumerate(sorted_dict.keys()):
                sorted[idx]=chiralNeighbours[hm[value]]
        if wedgemap[mol.GetBondBetweenAtoms(parent.atom.GetIdx(), chiralCenter.GetIdx()).GetIdx()]==0:
            numbyangle-=1
            for idx,neighbour in enumerate(chiralNeighbours):
                if neighbour.GetIdx() != parent.atom.GetIdx:
                    if wedgemap[mol.GetBondBetweenAtoms(chiralCenter.GetIdx(), neighbour.GetIdx()).GetIdx()]==0:
                        hm[give_angle(mol, chiralCenter, parent.atom, neighbour)]=idx
                    if wedgemap[mol.GetBondBetweenAtoms(chiralCenter.GetIdx(), neighbour.GetIdx()).GetIdx()]==1:
                        sorted[0]=neighbour
                    if wedgemap[mol.GetBondBetweenAtoms(chiralCenter.GetIdx(), neighbour.GetIdx()).GetIdx()]==6:
                        sorted[len(sorted)-1]=neighbour
            sorted_dict = collections.OrderedDict(hm)
            if numbyangle==2:
                sorted_dict[len(sorted-1)]=chiralNeighbours[hm[0]]
                if give_angle_from_middle(mol, chiralCenter, parent,chiralNeighbours[hm[1]])<0:
                    dummy = sorted[len(sorted)-2]
                    sorted[len(sorted)-2] = sorted[0]
                    sorted[0] = dummy
            if numbyangle==3:
                for idx, value in enumerate(sorted_dict.keys()):
                    sorted[idx+1]=chiralNeighbours[hm[value]]
    if sorted[0]==None:#This is in case none of the cases applies
        return None
    for atomid, atom in enumerate(sorted):
        if parent is not None and parent.atom.GetIdx()==atom.GetIdx():
            sortednodes[atomid]=parent        
        for node in nodes:
            if node.atom.GetIdx()==atom.GetIdx():
                sortednodes[atomid]=node
                break
    return sortednodes

def give_angle(mol, fromatom, to1, to2):
    return give_angle_both(mol, fromatom, to1, to2, True)

def give_angle_from_middle(mol, fromatom, to1, to2):
    return give_angle_both(mol, fromatom, to1, to2, False)

def give_angle_both(mol, fromatom, to1, to2, bool):
    fromcoords=mol.GetConformer().GetAtomPosition(fromatom.GetIdx())
    to1coords=mol.GetConformer().GetAtomPosition(to1.GetIdx())
    to2coords=mol.GetConformer().GetAtomPosition(to2.GetIdx())
    angle1 = math.atan2((to1coords.y - fromcoords.y), (to1coords.x - fromcoords.x))
    angle2 = math.atan2((to2coords.y - fromcoords.y), (to2coords.x - fromcoords.x))
    angle = angle2 - angle1
    if angle2 < 0 and angle1 > 0 and angle2 < -(math.pi / 2):
        angle = math.pi + angle2 + math.pi - angle1
    if angle2 > 0 and angle1 < 0 and angle1 < -(math.pi / 2):
        angle = -math.pi + angle2 - math.pi - angle1
    if bool and angle < 0:
        return 2 * math.pi + angle
    else:
        return angle
    
def is_left(mol, whereIs, viewFrom, viewTo):
    angle = give_angle_both(mol, viewFrom, viewTo, whereIs, False)
    if angle < 0:
        return False
    else:
        return True

def stereos_are_opposite(container, centre, wedgemap):
    atoms = centre.GetNeighbors()
    hm = {}
    for k,atom in enumerate(atoms):
        hm[give_angle(container, atom, atoms[0], atom)]=k
    sorted_dict = collections.OrderedDict(hm)
    ohere = list(sorted_dict.values())
    stereoOne = wedgemap[container.GetBondBetweenAtoms(centre.GetIdx(), atoms[0].GetIdx()).GetIdx()]
    stereoOpposite = wedgemap[container.GetBondBetweenAtoms(centre.GetIdx(), atoms[ohere[1]].GetIdx()).GetIdx()]
    return stereoOpposite == stereoOne

def is_square_planar(container, centre, wedgemap):
    atoms = centre.GetNeighbors()
    if len(atoms) != 4:
        return False
    bonds = centre.GetBonds()
    up = 0
    down = 0
    for bond in bonds:
        if wedgemap[bond.GetIdx()]==1:
            up += 1
        elif wedgemap[bond.GetIdx()]==6:
            down += 1
    return up == 2 and down == 2 and not stereos_are_opposite(container, centre, wedgemap)

def is_trigonal_bipyramidal_or_octahedral(container, centre, wedgemap):
    atoms = centre.GetNeighbors()
    if len(atoms) < 5 or len(atoms)>6:
        return 0
    bonds = centre.GetBonds()
    up = 0
    down = 0
    for bond in bonds:
        if wedgemap[bond.GetIdx()]==1:
            up += 1
        elif wedgemap[bond.GetIdx()]==6:
            down += 1
    if up == 1 and down == 1:
        if len(atoms) == 5:
            return 1
        else:
            return 2
    return 0

def create_wedgemap_from_block(molblock):
    count = 0
    atoms = 0
    bonds = 0
    skip=False
    wedgemap={}
    for line in molblock:
        count += 1
        if skip and count>atoms+4 and count<atoms+5+bonds and not "M  END" in line:
            wedgemap[count-atoms-5]=int(line[9:12].strip())
        elif count==4:
            atoms = int(line[0:3].strip())
            bonds = int(line[3:6].strip())
            skip=True
        elif "M  END" in line:
            break
    return wedgemap


def create_wedgemap(molfile):
    file1 = open(molfile, 'r')
    lines = file1.readlines()
    return create_wedgemap_from_block(lines)

def makeUpDownBonds(container, wedgemap):
    for a in container.GetAtoms():
        if len(a.GetNeighbors()) == 4:
            up = 0
            down = 0
            hs = 0
            h = None
            for conAtom in a.GetNeighbors():
                stereo = wedgemap[container.GetBondBetweenAtoms(a.GetIdx(), conAtom.GetIdx()).GetIdx()]
                if stereo == 1:
                    up+=1
                elif stereo == 6:
                    down+=1
                elif (stereo == 0 or stereo==4) and conAtom.GetSymbol()=="H":
                    h = conAtom
                    hs+=1
            if up == 0 and down == 1 and h is not None and hs == 1:
                wedgemap[container.GetBondBetweenAtoms(a.GetIdx(), h.GetIdx()).GetIdx()]=1
            if up == 1 and down == 0 and h is not None and hs == 1:
                wedgemap[container.GetBondBetweenAtoms(a.GetIdx(), h.GetIdx()).GetIdx()]=6