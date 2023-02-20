import pandas as pd
from tqdm.auto import tqdm
from rdkit import Chem
import argparse
from molvs.standardize import Standardizer

def query_standardize(mol):
    try:
        mol_sta = Standardizer().charge_parent(Standardizer().fragment_parent(mol), skip_standardize=True)
        return mol_sta
    except:
        return mol

def filter (database,unwanted):
    
    #fail =[]
    matches = []
    clean = []
    error =[]
    
    ## Getting structures from database set
    
    moldb = pd.read_csv(database,index_col=0)
    substructures = pd.read_csv(unwanted, sep="\s+")
    
    try:
        
      substructures["rdkit_molecule"] = substructures.smarts.apply(Chem.MolFromSmarts)
    except :
        print("error in smarts :",substructures.smarts)
    
    for index, row in tqdm(moldb.iterrows(), total=moldb.shape[0]):
        try:
            molecule = query_standardize(Chem.MolFromSmiles(row.smiles))
            match = False
            for _, substructure in substructures.iterrows():
                if molecule.HasSubstructMatch(substructure.rdkit_molecule):
                    matches.append(
                        {
                            "rdkit_molecule": row.smiles,
                            "substructure": substructure.smarts,
                            "substructure_name": substructure["name"],
                        }
                    )
                    match = True
            if not match:
                clean.append(index)  
        except:
              error.append(row.smiles)
              print("error with smiles :",row.smiles)      
    if error:
        err = pd.DataFrame (error, columns = ['error_smiles'])
        err.to_csv("Error_smiles.csv",index=False)
        
    matches = pd.DataFrame(matches)
    nonmatches = moldb.loc[clean]
    print(f"Number of compounds with unwanted structures: {len(matches)}")
    print(f"Number of compounds without unwanted structures: {len(nonmatches)}")
    return matches , nonmatches


def main():
    
    parser = argparse.ArgumentParser(
                    prog = 'Filter Structures',
                    description = 'Filter unwanted structures from input smiles files generates two output files unmatch and matches')
    parser.add_argument("-u", "--unwanted", help="specify path of input csv smarts file for query molecules", metavar="")
    parser.add_argument("-d", "--database", help="specify path/name of the database smiles file")
    parser.add_argument("-o", "--output", help="specify name/path of output csv file",default="Filter_unwatned.csv", metavar="")
                    
    args = parser.parse_args()
    
    match,unmatch = filter(args.database,args.unwanted)                
    
    match.to_csv("match_"+args.output,index=False)
    unmatch.to_csv("unmatch_"+args.output,index=False)
    
if __name__ == "__main__":
    main()    