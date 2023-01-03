import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def _InitialiseNeutralisationReactions():
    
    patts= (
    # Imidazoles
    ('[n+;H]','n'),
    # Amines
    ('[N+;!H0]','N'),
    # Carboxylic acids and alcohols
    ('[$([O-]);!$([O-][#7])]','O'),
    # Thiols
    ('[S-;X1]','S'),
    # Sulfonamides
    ('[$([N-;X2]S(=O)=O)]','N'),
    # Enamines
    ('[$([N-;X2][C,N]=C)]','N'),
    # Tetrazoles
    ('[n-]','[nH]'),
    # Sulfoxides
    ('[$([S-]=O)]','S'),
    # Amides
    ('[$([N-]C=O)]','N'),
    )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]


_reactions = None 
def NeutraliseCharges(smiles,reactions = None ): 
    global  _reactions 
    if  reactions  is  None : 
        if  _reactions  is  None : 
            _reactions = _InitialiseNeutralisationReactions() 
        reactions = _reactions 
    mol  =  Chem.MolFromSmiles(smiles) 
    replaced  =  False 
    for  i ,( reactant , product)  in  enumerate (reactions): 
        while  mol.HasSubstructMatch (reactant): 
            replaced  =  True 
            rms  =  AllChem.ReplaceSubstructs (mol,reactant,product ) 
            mol  =  rms [0] 
    if  replaced : 
        return (Chem.MolToSmiles(mol)) 
    else: 
        return (smiles)

def dataset_wash(task_name):
    # parameter
    # load data set
    data = pd.read_csv('{}.csv'.format(task_name))
    origin_data_num = len(data)
    data = data[data["smiles"].str.contains(r"\*") == False]

    # remove molecule can't processed by rdkit_des
    print('********dealing with compounds with rdkit_des*******')
    smiles_list = data['smiles'].values.tolist()
    cant_processed_smiles_list = []
    for index, smiles in enumerate(smiles_list):
        if index % 10000 == 0:
            print(index)
        try:
            molecule = Chem.MolFromSmiles(smiles)
            smiles_standard = Chem.MolToSmiles(molecule)
            data['smiles'][index] = smiles_standard
        except:
            cant_processed_smiles_list.append(smiles)
            data.drop(index=index, inplace=True)
    print("compounds can't be processed by rdkit_des: {} molecules, {}\n".format(len(cant_processed_smiles_list),
                                                                      cant_processed_smiles_list))

    # remove mixture and salt
    print('********dealing with inorganic compounds*******')
    data = data.reset_index(drop=True)
    smiles_list = data['smiles'].values.tolist()
    mixture_salt_list = []
    for index, smiles in enumerate(smiles_list):
        if index % 10000==0:
            print(index)
        symbol_list = list(smiles)
        if '.' in symbol_list:
            mixture_salt_list.append(smiles)
            data.drop(index=index, inplace=True)
    print('inorganic compounds: {} molecules, {}\n'.format(len(mixture_salt_list), mixture_salt_list))


    # remove inorganic compounds
    print('********dealing with inorganic compounds*******')
    data = data.reset_index(drop=True)
    smiles_list = data['smiles'].values.tolist()
    inorganics = []
    atom_list = []
    for index, smiles in enumerate(smiles_list):
        if index % 10000==0:
            print(index)
        mol = Chem.MolFromSmiles(smiles)
        count = 0
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C':
                break
            else:
                count += 1
        if count == mol.GetNumAtoms():
            inorganics.append(smiles)
            data.drop(index=index, inplace=True)
    print('inorganic compounds: {} molecules, {}\n'.format(len(inorganics), inorganics))
    return(data)
