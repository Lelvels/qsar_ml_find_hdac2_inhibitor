from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import Descriptors, MACCSkeys
from rdkit.ML.Descriptors import MoleculeDescriptors

#MACCS
from tqdm import tqdm
import numpy as np

class PharmacyCommon:
    def __init__(self):
        pass
    
    def gen_maccs_fpts(self, data) -> np.ndarray:
        Maccs_fpts = []
        count = 0
        with tqdm(total=len(data), desc='Progress') as pbar:
            for i in data:
                try:
                    mol = Chem.MolFromSmiles(i)
                except:
                    print("An exception occurred with " + str(count))
                    continue
                fpts = MACCSkeys.GenMACCSKeys(mol)
                mfpts = np.array(fpts)
                Maccs_fpts.append(mfpts)
                count += 1
                pbar.update(1)  # Update the progress bar
        return np.array(Maccs_fpts)
    
    def gen_ecfp4_fpts(self, data, bits) -> np.ndarray:
        if bits not in [1024, 2048]:
            raise ValueError("Invalid value for bits. Must be either 1024 or 2048.")

        result_fpts = []
        count = 0
        with tqdm(total=len(data), desc='Progress') as pbar:
            for i in data:
                try:
                    mol = Chem.MolFromSmiles(i)
                except:
                    print("An exception occurred with " + str(count))
                    continue
                fpts = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, bits))
                result_fpts.append(fpts)
                count += 1
                pbar.update(1)  # Update the progress bar
        return np.array(result_fpts) 
    
    def gen_ecfp6_fpts(self, data, bits) -> np.ndarray:
        if bits not in [1024, 2048]:
            raise ValueError("Invalid value for bits. Must be either 1024 or 2048.")
        result_fpts = []
        count = 0
        with tqdm(total=len(data), desc='Progress') as pbar:
            for i in data:
                try:
                    mol = Chem.MolFromSmiles(i)
                except:
                    print("An exception occurred with " + str(count))
                    continue
                fpts = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 3, bits))
                result_fpts.append(fpts)
                count += 1
                pbar.update(1)  # Update the progress bar
        return np.array(result_fpts) 
    
    """_summary_
    Rdkit parameters:
        minPath: the minimum path length (in bonds) to be included
        maxPath: the maximum path length (in bonds) to be included
        useHs: toggles inclusion of Hs in paths (if the molecule has explicit Hs)
        branchedPaths: toggles generation of branched subgraphs, not just linear paths
        useBondOrder: toggles inclusion of bond orders in the path hashes
        countSimulation: if set, use count simulation while generating the fingerprint
        countBounds: boundaries for count simulation, corresponding bit will be set if the count is higher than the number provided for that spot
        fpSize: size of the generated fingerprint, does not affect the sparse versions
        numBitsPerFeature: the number of bits set per path/subgraph found
        atomInvariantsGenerator: atom invariants to be used during fingerprint generation
    Default parameters:
        minPath=1, 
        maxPath=7, 
        useHs=True, 
        branchedPaths=True, 
        useBondOrder=True, 
        countSimulation=False,
        (AtomPairsParameters)countBounds=None, 
        (int)fpSize=2048, 
        (int)numBitsPerFeature=2, 
        (AtomPairsParameters)atomInvariantsGenerator=None
    """
    def gen_rdkit_fpts(self, data, bits) -> np.ndarray:
        if bits not in [1024, 2048]:
            raise ValueError("Invalid value for bits. Must be either 1024 or 2048.") 
        result_fpts = []
        count = 0
        with tqdm(total=len(data), desc='Progress') as pbar:
            for i in data:
                try:
                    mol = Chem.MolFromSmiles(i)
                except:
                    print("An exception occurred with " + str(count))
                    continue
                fpts = np.array(AllChem.GetRDKitFPGenerator(mol, 3, bits))
                result_fpts.append(fpts)
                count += 1
                pbar.update(1)  # Update the progress bar
        return np.array(result_fpts) 