import pandas as pd
import numpy as np
from utils.smile import *
from utils import helper
from rdkit import Chem
#from rdkit.Chem import AllChem
from utils.depictsmile import DepictSmile
from utils.descriptor import Descriptor
import os
import warnings
import sklearn
from sklearn import pipeline
#from sklearn.model_selection import train_test_split, KFold, GridSearchCV
#from sklearn import metrics
from xgboost import XGBRegressor
import joblib
warnings.filterwarnings('ignore')
pd.options.mode.chained_assignment=None

class MaterialStructuralPropertyPrediction:
    def __init__(self, props={"MP":True, "BP":False, "HC":False, "dH":False, "h2_uptake":False}, imgpath=None, batch=128, verbose=0):
        self.__smile__ = SMILE(verbose)
        self.__model_path__= helper.resource_path(os.path.join("include", "best_models"))
        self.__feature_path__= helper.resource_path(os.path.join("include", "feature"))

        #print(f"Model path: {self.__model_path__}")
        #print(f"Feature path: {self.__feature_path__}")

        self.__batch__ = batch
        self.__verbose__ = verbose
        self.__plot_smile__ = None
        self.props = props
        if imgpath is not None:
            self.__plot_smile__ = DepictSmile(save_file_path=self.__create_path__(imgpath))
        self.df_mol = None

    def __create_path__(self, fpath):
        if os.path.exists(fpath):
            return fpath
        
        try:
            os.makedirs(fpath, exist_ok=True)
        except Exception as e:
            print(f"Error in [__create_path__] to create file: {e}")
        
        return fpath
        
    def __read_features__(self, fileWPath):
        if not os.path.exists(fileWPath):
            print(f'Feature files not exists: {fileWPath}')
            return []

        fl = []
        with open(fileWPath, 'r') as file:
            for row in file:
                fl.append(row.rstrip('\n'))
        return fl

    def __download_features__(self, v_smile, si, ei):
        if si>=ei:
            return

        if self.__verbose__>0:
            print(f"[{si}, {ei}] processing.")
            
        try:
            desc = Descriptor()
            descriptors = pd.DataFrame(desc.from_smiles(v_smile[si:ei]))
            descriptors['smiles'] = v_smile[si:ei]
            del [desc]
            
            if self.df_mol is None:
                self.df_mol = pd.DataFrame(descriptors)
            else:
                self.df_mol = pd.concat([self.df_mol, pd.DataFrame(descriptors)], ignore_index=True)

            #if self.__verbose__>0:
            #    print(f"[{si}, {ei}] done.")
            
        except Exception as e:
            print(f"Error processing [{si}-{ei}]: {e}.")

            if si+1==ei:
                return
            
            mi = si + (ei-si)//2
            self.__download_features__(v_smile, si, mi)
            if self.__verbose__>0:
                print(f"[{si}, {mi}] done.")
            
            try:
                if self.__verbose__>0:
                    print(f"{mi} processing.")
                    
                descriptors = pd.DataFrame([from_smiles(v_smile[mi])])
                descriptors['smiles'] = v_smile[mi]

                if self.df_mol is None:
                    self.df_mol = pd.DataFrame(descriptors)
                else:
                    self.df_mol = pd.concat([self.df_mol, pd.DataFrame(descriptors)], ignore_index=True)

                if self.__verbose__>0:
                    print(f"{mi} done.")
            except Exception as e:
                print(f"Error for '{v_smile[mi]}' is: {e}")

            if self.__verbose__>0:
                print(f"[{mi+1}, {ei}] processing.")
            self.__download_features__(v_smile, mi+1, ei)
            if self.__verbose__>0:
                print(f"[{mi+1}, {ei}] done.")


    def __get_features__(self, v_smile):
        df_mol = None
        
        if len(v_smile)>0:
            si = 0
            ei = 0
            
            while ei<len(v_smile):
                si=ei
                ei += self.__batch__
                ei = min(ei, len(v_smile))

                self.__download_features__(v_smile, si, ei)

                if self.__verbose__>0:
                    print(f"[{si}, {ei}] done.")
                
                '''
                if self.__verbose__>0:
                    print(f"si={si}, ei={ei}")
                    
                try:
                    descriptors = pd.DataFrame(from_smiles(v_smile[si:ei]))
                    descriptors['smile'] = v_smile[si:ei]
                    
                    if df_mol is None:
                        df_mol = pd.DataFrame(descriptors)
                    else:
                        df_mol = pd.concat([df_mol, pd.DataFrame(descriptors)], ignore_index=True)
                    
                except Exception as e:
                    print(f"Error processing [{si}-{ei}]: {e}")
                    for i in range(si, ei):
                        try:
                            descriptors = pd.DataFrame([from_smiles(v_smile[i])])
                            descriptors['smile'] = v_smile[i]

                            if df_mol is None:
                                df_mol = pd.DataFrame(descriptors)
                            else:
                                df_mol = pd.concat([df_mol, pd.DataFrame(descriptors)], ignore_index=True)
                        except Exception as e:
                            print(f"Error for '{v_smile[i]}' is: {e}")
                '''
            
        # Replace missing values to 0
        for c in self.df_mol.columns:
            if c != 'smiles':
                self.df_mol.loc[self.df_mol[c]=='', c] = '0.0'
                self.df_mol[c] = self.df_mol[c].astype(np.float32)

    def __validate_features__(self, df):
        # Replace inf and -inf with NaN first, then handle NaN
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        # Then replace NaN values with 0, for example
        df.fillna(0, inplace=True)
    
    def __predict_MP__(self):
        model=None
        try:
            model_path = os.path.join(self.__model_path__, 'MP', 'bestmodel_LBNL.sav')
            #print(f"Loading MP model from: {model_path}")
            model = joblib.load(model_path)
        except Exception as e:
            print(f"Error to load MP model: {e}")
                
        if model is None:
            return None
        
        feature_path = os.path.join(self.__feature_path__, 'MP_features.txt')
        #print(f"Loading MP features from: {feature_path}")
        features = self.__read_features__(feature_path)
        if len(features)==0:
            return None
        
        fv = self.df_mol.loc[:, features]
        if fv is None:
            return None
            
        self.__validate_features__(fv)
        
        try:
            pv = model.predict(fv)
        except Exception as e:
            print(f'Error in prediction MP: {e}')
            return None
            
        return pv
    
    def __predict_BP__(self):
        model=None
        try:
            model = joblib.load(os.path.join(self.__model_path__, 'BP', 'bestmodel.sav'))
        except Exception as e:
            print(f"Error to load model: {e}")
                
        if model is None:
            return None
        
        features = self.__read_features__(os.path.join(self.__feature_path__, 'BP_features.txt'))
        if len(features)==0:
            return None
        
        fv = self.df_mol.loc[:, features]
        if fv is None:
            return None
            
        self.__validate_features__(fv)
        
        try:
            pv = model.predict(fv)
        except Exception as e:
            print(f'Error in prediction BP: {e}')
            return None
            
        return pv
    
    def __predict_HC__(self):
        model=None
        try:
            model = joblib.load(os.path.join(self.__model_path__, 'HC', 'bestmodel.sav'))
        except Exception as e:
            print(f"Error to load model: {e}")
                
        if model is None:
            return None
        
        features = self.__read_features__(os.path.join(self.__feature_path__, 'HC_features.txt'))
        if len(features)==0:
            return None
        
        fv = self.df_mol.loc[:, features]
        if fv is None:
            return None
            
        self.__validate_features__(fv)
        
        try:
            pv = model.predict(fv)
        except Exception as e:
            print(f'Error in prediction HC: {e}')
            return None
            
        return pv
    
    def __predict_HH__(self):
        model=None
        try:
            model = joblib.load(os.path.join(self.__model_path__, 'HH' ,'bestmodel.sav'))
        except Exception as e:
            print(f"Error to load model: {e}")
                
        if model is None:
            return None
        
        features = self.__read_features__(os.path.join(self.__feature_path__, 'HH_features.txt'))
        if len(features)==0:
            return None
        
        fv = self.df_mol.loc[:, features]
        if fv is None:
            return None
            
        self.__validate_features__(fv)
        
        try:
            pv = model.predict(fv)
        except Exception as e:
            print(f'Error in prediction HH: {e}')
            return None
            
        return pv

    def __add_hydrogens__(self, emol, atom_idx, num_h):
        for _ in range(num_h):
            new_atom_idx = emol.AddAtom(Chem.Atom(1))  # Hydrogen atom
            emol.AddBond(atom_idx, new_atom_idx, Chem.BondType.SINGLE)

    def __saturate_smiles__(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return '',0
    
        # Add explicit hydrogens to the molecule
        mol = Chem.AddHs(mol)
    
        # Create an editable molecule object
        emol = Chem.EditableMol(mol)
    
        # Saturate non-aromatic double and triple bonds
        for bond in mol.GetBonds():
            bond_type = bond.GetBondType()
            begin_atom_idx = bond.GetBeginAtomIdx()
            end_atom_idx = bond.GetEndAtomIdx()
            
            if bond_type == Chem.BondType.DOUBLE:
                emol.RemoveBond(begin_atom_idx, end_atom_idx)
                emol.AddBond(begin_atom_idx, end_atom_idx, Chem.BondType.SINGLE)
                self.__add_hydrogens__(emol, begin_atom_idx, 1)
                self.__add_hydrogens__(emol, end_atom_idx, 1)
            elif bond_type == Chem.BondType.TRIPLE:
                emol.RemoveBond(begin_atom_idx, end_atom_idx)
                emol.AddBond(begin_atom_idx, end_atom_idx, Chem.BondType.SINGLE)
                self.__add_hydrogens__(emol, begin_atom_idx, 2)
                self.__add_hydrogens__(emol, end_atom_idx, 2)
        
        # Generate the saturated molecule
        saturated_mol = emol.GetMol()
        Chem.SanitizeMol(saturated_mol)
        
        # Generate saturated SMILES
        saturated_smiles = Chem.MolToSmiles(saturated_mol)
        
        # Compute the required moles of hydrogen
        original_mol = Chem.MolFromSmiles(smiles)
        num_unsaturated_bonds = sum(1 for bond in original_mol.GetBonds() if bond.GetBondType() in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE])
        moles_of_hydrogen = num_unsaturated_bonds * 2
    
        return saturated_smiles, moles_of_hydrogen
    
    def __predict_properties__(self, v_smile):
        self.__get_features__(v_smile)
        bpos = -1

        # Predict melting point
        if self.props['MP']:
            pv = self.__predict_MP__()
            if pv is None:
                return None
            self.df_mol['MP_pred_K'] = pv
            bpos-=1
        
        # Predict boiling point
        if self.props['BP']:
            pv = self.__predict_BP__()
            if pv is None:
                return None
            self.df_mol['BP_pred_K'] = pv
            bpos-=1
        
        # Predict Heat capacity
        if self.props['HC']:
            pv = self.__predict_HC__()
            if pv is None:
                return None
            self.df_mol['HC_pred_JK'] = pv
            bpos-=1
        
        # Predict Heat of Hydrogenation
        if self.props['dH']:
            pv = self.__predict_HH__()
            if pv is None:
                return None
            self.df_mol['dH_pred_JK'] = pv
            bpos-=1

        self.df_mol = pd.DataFrame(self.df_mol.iloc[:, bpos:])

        if self.props['h2_uptake']:
            sat_smiles, h2_upkate=[],[]
            for s in self.df_mol.smiles.values.tolist():
                if self.__plot_smile__:
                    self.__plot_smile__.save_svg_to_png(s)
                ssv, h2v = self.__saturate_smiles__(s)
                sat_smiles.append(ssv)
                h2_upkate.append(h2v)
    
            self.df_mol['Saturated_smiles'] = sat_smiles
            self.df_mol['H2_uptake'] = h2_upkate
        
    def get_data(self):    
        return self.df_mol
        
    def prediction_from_smile(self, smiles):
        right_smiles, wrong_smiles = self.__smile__.get_valid_smiles(smiles)
        self.__predict_properties__(right_smiles)
        return wrong_smiles
    
    def prediction_from_pid(self, pubChemIDs):
        if isinstance(pubChemIDs, int):
            pubChemIDs = [pubChemIDs]
        can_smiles = self.__smile__.get_valid_smiles_from_pid(pubChemIDs)
        return self.prediction_from_smile(can_smiles)