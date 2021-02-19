import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Draw import SimilarityMaps

# This code is for drawing similarity maps based on predicted fingerprints.
def get_highlight(SMILES,output_pred_fingerprint): #input is smiles and predicted fingerprints
    radi = 2
    bits = 2048
    mol = Chem.MolFromSmiles(SMILES)
    mol_H = Chem.AddHs(mol)
    mol_bi_H_QC = {0:[],1:[],2:[]}
    mol_tpls_H = {}
    weight = np.full((mol.GetNumAtoms()),-1.5,float)
    for r in range(radi+1):
        mol_bi_H = {}
        mol_fp_H = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol_H, radius=r, bitInfo=mol_bi_H, nBits = bits)

        for i in mol_fp_H.GetOnBits():
            idx = mol_bi_H[i][0][0]
            radius_list = []
            for j in range(len(mol_bi_H[i])):
                atom_radi = mol_bi_H[i][j][1]
                radius_list.append(atom_radi)
            atom = mol_H.GetAtomWithIdx(idx)
            symbol = atom.GetSymbol()
            neigbor = [x.GetAtomicNum() for x in atom.GetNeighbors()]

            if r in radius_list:# and symbol == 'C' and 1 in neigbor:#radius = 2, atom = Carbon, H possessed Carbon
                mol_bi_H_QC[r].append(i)
        for k in mol_bi_H_QC[r]:
            mol_tpls_H[2048*r+k] = {'Mol':mol_H,'Bit_info':mol_bi_H,'Radius':r}
            
    for bit in output_pred_fingerprint:
    # getting highlight of atom
        try:
            for j in range(len(mol_tpls_H[bit]['Bit_info'][bit%2048])):#for i in range(len(mol_tpls_H[[0][2][bit_n]])):
                if mol_tpls_H[bit]['Bit_info'][bit%2048][j][1] == mol_tpls_H[bit]['Radius']:
                    atom = mol_tpls_H[bit]['Bit_info'][bit%2048][j][0]
                    #if mol_tpls_H[bit]['Radius'] >=1:
                    weight[atom] += 1
                
                    

        except:
            pass
    weight = weight
    #weight = (weight ==1).a
    fig = SimilarityMaps.GetSimilarityMapFromWeights(mol,weight.tolist(),contourLines=0,step=0.01,sigma=0.03,colors=None)
    plt.close()
    return fig
