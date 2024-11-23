# Search-by-SMILES-across-the-excel-sheet
Search by SMILES across the Excel sheet (RDKit)

1. Import Panda and RDKit libraries
```
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

from io import StringIO
```

2. Input Excel sheet directory, and 
```
input_xlsx = input("Enter the directory of the Excel sheet: ")
data_df = pd.read_excel(input_xlsx)
index_flag = input("Does the Excel sheet have an index? Type 'True' or 'False': ").lower() == 'false'
```

3. Convert it to CSV
```
csv_data = data_df.to_csv(index=index_flag)
csv_file = StringIO(csv_data)
df = pd.read_csv(csv_file)
```

3. Convert SMILES to MOLS by RDKit
```
df['mols'] = [Chem.MolFromSmiles(smiles) for smiles in df['SMILES']]
```

4. Search by SMILES and save it
```
cur_scaf_smiles = input('write SMILES of interest, ex: c2nn([C@H]3CCCO3)c3ncncc23')
mol_scaffold_of_interest = Chem.MolFromSmiles(cur_scaf_smiles)
mols_to_display = []
for mol_smi in df["SMILES"]:
    mol = Chem.MolFromSmiles(mol_smi)
    if mol is None:
        continue
    hit_ats = list(mol.GetSubstructMatch(mol_scaffold_of_interest))
    if hit_ats:
        mols_to_display.append((mol, hit_ats))

mols = [item[0] for item in mols_to_display]
hit_ats_list = [item[1] for item in mols_to_display]
img1 = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(300, 300), 
                           legends=[smile for smile in df["Code"]], highlightAtomLists=hit_ats_list
                            ,returnPNG=False)

img1.save(input("/directory/../name.png"))
display(img)
```

5. Highlight the SMILES of interest
```
def get_atoms_to_highlight(mol_scaffold_of_interest, mol_compound_of_interest):
    atom_list = [i for i in range(mol_compound_of_interest.GetNumAtoms())]
    grouped_common_atoms = mol_compound_of_interest.GetSubstructMatches(mol_scaffold_of_interest)
    flattened_grouped_common_atoms = list(sum(grouped_common_atoms, ()))
    atoms_to_highlight = list(set(atom_list) - set(flattened_grouped_common_atoms))
    return atoms_to_highlight

cur_scaf_smiles = input('Write the SMILES of interest, e.g., c2nn([C@H]3CCCO3)c3ncncc23: ')

mol_scaffold_of_interest = Chem.MolFromSmiles(cur_scaf_smiles)
hl_atom_lists = []

for mol in df["mols"]:
    # Highlight atoms using the get_atoms_to_highlight function
    highlight = get_atoms_to_highlight(mol_scaffold_of_interest, mol)
    mol.__sssAtoms = highlight
    hl_atom_lists.append(highlight)

img = Draw.MolsToGridImage(df["mols"], 
                           legends=[name for name in df["Code"]], 
                           highlightAtomLists=hl_atom_lists, 
                           molsPerRow=10, maxMols= 1000,
                           subImgSize=(300, 300), useSVG=False,returnPNG=False)

# Display the image
img.save(input("/directory/.../name.png"))
display(img)
```

