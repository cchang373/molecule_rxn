#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 09:27:38 2020

@author: cchang373
"""
from rdkit import Chem
from json import load,dump
import rdkit

def add_single_metal(smi, metal_atom):
    mol = Chem.MolFromSmiles(smi)
    smi = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smi)
    mol_H = Chem.AddHs(mol)
    for num, atom_object in enumerate(mol_H.GetAtoms()):
        if atom_object.GetSymbol() != 'H' and atom_object.GetSymbol() != metal_atom:
            count_bond = 0
            saturate = 4 if atom_object.GetSymbol() == 'C' else 2
            for bonds in mol_H.GetBonds():
                if num in [bonds.GetBeginAtomIdx(), bonds.GetEndAtomIdx()]:
                    count_bond += 1 if bonds.GetBondType() == rdkit.Chem.rdchem.BondType.SINGLE else 2
            if saturate > count_bond:
                string_add = '('+metal_atom+')'
                count_heavy = 0
                for idx, char in enumerate(smi):
                    if count_heavy == num and (char == 'C' or char == 'O'):
                        element = 'C' if char == 'C' else 'O'
                        if count_heavy != 0 and smi[idx-1] == '[':
                            try:
                                a = int(smi[idx+2])
                                smi_Pt = smi[:idx+4] + string_add + smi[idx+4:]
                            except:
                                smi_Pt = smi[:idx+3] + string_add + smi[idx+3:] if smi[idx+1]=='H' else smi[:idx+2]+string_add+smi[idx+2:]
                        elif count_heavy != 0 and smi[idx-1] != '[':
                            smi_Pt = smi[:idx+1]+string_add+smi[idx+1:]
                        elif count_heavy == 0 and idx != 0:
                            try:
                                a = int(smi[idx+2])
                                smi_Pt = smi[:idx+4] + string_add + smi[idx+4:]
                            except:
                                smi_Pt = smi[:idx+3] + string_add + smi[idx+3:] if smi[idx+1] == 'H' else smi[:idx+2]+string_add+smi[idx+2:]
                        elif count_heavy == 0 and idx == 0:
                            smi_Pt = smi[:idx+1]+string_add+smi[idx+1:]
                        mol_Pt = Chem.MolFromSmiles(smi_Pt)
                        smi_Pt = Chem.MolToSmiles(mol_Pt)
                        return smi_Pt
                    if char == 'C' or char == 'O' or char == 'P' or char == 'A' or char == 'R':
                        count_heavy += 1
    return

class smi_single():
    def __init__(self):
        self.single_smi = []
    def smi_single(self,smi_list):
        new_smi = []
        for smi in smi_list:
            if 'O=' in smi or '=O' in smi:
                mol=Chem.MolFromSmiles(smi)
                smi=Chem.MolToSmiles(mol)
                mol=Chem.MolFromSmiles(smi)
                mol_H=Chem.AddHs(mol)
                for num, atom_object in enumerate(mol_H.GetAtoms()):
                    if atom_object.GetSymbol() =='C' and 'O' in [x.GetSymbol() for x in atom_object.GetNeighbors()]:
                        bonds=atom_object.GetBonds()
                        for bond in bonds:
                            if bond.GetBondType() == rdkit.Chem.rdchem.BondType.DOUBLE:
                                if mol_H.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol() == 'O' or  mol_H.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol() == 'O':
                                    O_num = bond.GetBeginAtomIdx() if  mol_H.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol() == 'O' else bond.GetEndAtomIdx()
                                    C_num = num
                                    H_count = len([x.GetSymbol() for x in atom_object.GetNeighbors()if x.GetSymbol() == 'H'])
                                    if H_count == 0:
                                        H_add = ''
                                    elif H_count == 1:
                                        H_add = 'H'
                                    else:
                                        H_add = 'H' + str(H_count)
                                    count_C, count_O = 0,0
                                    for idx,char in enumerate(smi):
                                        if count_C == C_num and char == 'C':
                                            try:
                                                h=smi[idx-1:idx+2]
                                                if h == '[C]' or smi[idx-1:idx+3] == '[CH]':
                                                    smi_C = smi
                                                    break
                                                else:
                                                    smi_C = smi[:idx]+'[C'+H_add+']'+smi[idx+1:]
                                                    break
                                            except:
                                                smi_C = smi[:idx]+'[C'+H_add+']'+smi[idx+1:]
                                                break
                                        if char == 'C' or char == 'O' or char == 'P':
                                            count_C += 1
                                    for idx,char in enumerate(smi_C):
                                        if count_O == O_num and char == 'O':
                                            try:
                                                h = smi_C[idx+1]
                                                if smi_C[idx+1] == '=':
                                                    smi_O = smi_C[:idx]+'[O]'+smi_C[idx+2:]
                                                else:
                                                    smi_O = smi_C[:idx-1]+'[O]'+smi_C[idx+1:]
                                            except:
                                                smi_O = smi_C[:idx-1]+'[O]'+smi_C[idx+1:]
                                            break
                                        if char in 'COPRA':
                                            count_O += 1
                                    mol_O = Chem.MolFromSmiles(smi_O)
                                    smi_O = Chem.MolToSmiles(mol_O)
                                    if smi_O not in self.single_smi:
                                        self.single_smi.append(smi_O)
                                        new_smi.append(smi_O)
                                    #print(new_smi)
            if len(new_smi) == 0:
                return self.single_smi
            else:
                return self.smi_single(new_smi)
        return
#print(smi_single().smi_single(['[CH]OC(C)=O']))

class add_metal():
    def __init__(self):
    	self.smis_Pt = []
    	#self.smis_Pt_list = []
    	#self.l1 = 0
    	#self.l2 = 0    
    def add_metal(self,smi,metal):
        #print('flag: ',self.smis_Pt)
        mol = Chem.MolFromSmiles(smi)
        smi = Chem.MolToSmiles(mol)
        mol = Chem.MolFromSmiles(smi)
        mol_H = Chem.AddHs(mol)
        for num, atom_object in enumerate(mol_H.GetAtoms()):
            if atom_object.GetSymbol() != 'H' and atom_object.GetSymbol() != metal:
                count_bond = 0
                saturate = 4 if atom_object.GetSymbol() == 'C' else 2
                for bonds in mol_H.GetBonds():
                    if num in [bonds.GetBeginAtomIdx(), bonds.GetEndAtomIdx()]:
                        count_bond += 1 if bonds.GetBondType() == rdkit.Chem.rdchem.BondType.SINGLE else 2
                if saturate > count_bond:
                    string_add = ('(['+metal+'])') * (saturate-count_bond)
                    count_heavy = 0
                    for idx, char in enumerate(smi):
                        if count_heavy == num and (char == 'C' or char == 'O'):
                            if count_heavy != 0 and smi[idx-1] == '[':
                                try:
                                    a = int(smi[idx+2])
                                    smi_Pt = smi[:idx+4] + string_add + smi[idx+4:]
                                    #print('flag 1: ', smi_Pt)
                                except:
                                    smi_Pt = smi[:idx+3] + string_add + smi[idx+3:] if smi[idx+1]=='H' else smi[:idx+2]+string_add+smi[idx+2:]
                                    #print('flag 2: ', smi_Pt)
                            elif count_heavy != 0 and smi[idx-1] != '[':
                                smi_Pt = smi[:idx+1]+string_add+smi[idx+1:]
                                #print('flag 3: ', smi_Pt)
                            elif count_heavy == 0 and idx != 0:
                                try:
                                    a = int(smi[idx+2])
                                    smi_Pt = smi[:idx+4] + string_add + smi[idx+4:]
                                    #print('flag 4: ', smi_Pt)
                                except:
                                    smi_Pt = smi[:idx+3] + string_add + smi[idx+3:] if smi[idx+1] == 'H' else smi[:idx+2]+string_add+smi[idx+2:]
                                    #print('flag 5: ', idx, smi_Pt, smi)
                            elif count_heavy == 0 and idx == 0:
                                smi_Pt = smi[:idx+1]+string_add+smi[idx+1:]
                                #print('flag 6: ', smi_Pt)
                            mol_Pt = Chem.MolFromSmiles(smi_Pt)
                            smi_Pt = Chem.MolToSmiles(mol_Pt)
                            if smi_Pt not in self.smis_Pt:
                                self.smis_Pt.append(smi_Pt)
                            break
                        if char == 'C' or char == 'O' or char in 'PRA':
                            count_heavy += 1
    
        return self.smis_Pt

#print(add_metal().add_metal('[O]C[C]','Rh'))

class add_metal_all():
    def __init__(self):
        self.smis_Pt_list=[]
        #self.l1=0
        #self.l2=0
        
    def add_metal_all(self,smi_list,metal):
        new_list = []
        #l1=len(self.smis_Pt_list)
        #print('l1: ',l1)
        for smi in smi_list:
            smi_n_list = add_metal().add_metal(smi, metal)
            #print(smi_n_list,smi)
            if len(smi_n_list) == 0:
                continue
            for smi_n in smi_n_list:
                if smi_n not in self.smis_Pt_list:
                    self.smis_Pt_list.append(smi_n)
                    new_list.append(smi_n)
    
        #l2=len(self.smis_Pt_list)
        if len(new_list) == 0:
            #print('l1 end: ',l)
            #print('finished: ',self.smis_Pt_list)
            return self.smis_Pt_list
        else:
            #print('still looping: ',self.smis_Pt_list)
            return self.add_metal_all(new_list, metal)
#print(add_metal_all().add_metal_all(['[O]C[C]'],'Rh'))

