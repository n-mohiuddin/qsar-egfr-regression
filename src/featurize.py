# src/featurize.py
import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs, Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator

# Morgan/ECFP4: radius=2, 2048 bits
_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

DESCRIPTOR_COLS = ["MolWt", "MolLogP", "HBD", "HBA", "TPSA", "RotatableBonds"]

def compute_lipinski_descriptors(mol) -> np.ndarray:
    """Return 6 descriptors in the exact order used in training."""
    return np.array([
        Descriptors.MolWt(mol),                       # MolWt
        Crippen.MolLogP(mol),                         # MolLogP
        Lipinski.NumHDonors(mol),                     # HBD
        Lipinski.NumHAcceptors(mol),                  # HBA
        rdMolDescriptors.CalcTPSA(mol),               # TPSA
        rdMolDescriptors.CalcNumRotatableBonds(mol),  # RotatableBonds
    ], dtype=float)

def morgan_fp(mol) -> np.ndarray:
    """Return 2048-bit Morgan fingerprint as 0/1 float array."""
    fp = _MORGAN_GEN.GetFingerprint(mol)
    arr = np.zeros((2048,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.astype(float)

def featurize_smiles(smiles: str) -> np.ndarray:
    """
    SMILES -> (2054,) feature vector = [6 descriptors | 2048 Morgan bits]
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES: could not parse with RDKit.")

    x_desc = compute_lipinski_descriptors(mol)   # (6,)
    x_fp   = morgan_fp(mol)                      # (2048,)
    x = np.hstack([x_desc, x_fp])                # (2054,)
    return x

def featurize_many(smiles_list):
    """List[str] -> (n, 2054) feature matrix."""
    feats = [featurize_smiles(s) for s in smiles_list]
    return np.vstack(feats)
