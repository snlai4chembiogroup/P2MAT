######################################################
# Created: August, 2024
# By: Md Kamruzzaman (@mkamruz)
# Objective:
#          1. SMILEs string related methods
######################################################

from padelpy import from_smiles
from pubchempy import Compound
from rdkit import Chem
import rdkit, joblib
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment=None

class SMILE:
    def __init__(self, verbose=0):
        self.verbose=verbose
    
    def __checkValidity__(self, can_smile):
        retV = dict()
        for s in can_smile:
            try:
                rv = Chem.MolFromSmiles(s)
            except Exception as e:
                print(f"Error in smile string:{e}")
            finally:
                if isinstance(rv, rdkit.Chem.rdchem.Mol):
                    retV[s] = True
                else:
                    retV[s] = False

        return retV

    def isAeromatic(self, smile):
        aeromatic = False
        if self.__checkValidity__(smile):
            m = Chem.MolFromSmiles(smile)
            
            for idx in [a.GetIdx() for a in m.GetAtoms()]:
                #if m.GetAtomWithIdx(idx).GetIsAromatic():
                if m.GetRingInfo().IsAtomInRingOfSize(idx,6):
                    aeromatic = True
                    break
        return aeromatic
    
    def __fetch_smile__(self, pubChemId):
        comp = Compound.from_cid(pubChemId)
        return comp.isomeric_smiles
    
    def get_valid_smiles(self, can_smile):
        smile_dict = self.__checkValidity__(can_smile)
        vsmiles = []
        ivsmiles = []
        for k in smile_dict:
            if smile_dict[k]==True:
                vsmiles.append(k)
            else:
                ivsmiles.append(k)
        return vsmiles, ivsmiles
    
    def get_valid_smiles_from_pid(self, pubChemIdList):
        can_smile = []
        for pid in pubChemIdList:
            can_smile.append(self.__fetch_smile__(pid))
        return self.get_valid_smiles(can_smile)