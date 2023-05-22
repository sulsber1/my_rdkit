from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
import rdkit
import smiles_example
import time
import pandas as pd
import redis

pool = redis.ConnectionPool(host='localhost', port=6379, db=0, decode_responses=True)
r = redis.Redis(connection_pool=pool)


SMILES = smiles_example.return_small_smiles()

data=pd.read_csv("chembl_32_chemreps.txt", sep="\t", header=0)
SMILES = data['canonical_smiles']
del data

mols = {}
def generate_mol_files():
    # Generate Mol files
    mol_time_start = time.perf_counter()
    mols = [Chem.MolFromSmiles(x) for x in SMILES]
    mols = [x for x in mols if x is not None]
    mol_time_end = time.perf_counter()
    print(f"It took {mol_time_end - mol_time_start} seconds to generate mol files")
    return mols

# Generate fingerprints
def generate_fingerprints(mols):
    fp_start_time = time.perf_counter()
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3,fpSize=512)
    fps = [(mfpgen.GetFingerprint(x)) for x in mols]
    fp_end_time = time.perf_counter()
    print(f"It took {fp_end_time - fp_start_time} seconds to generate fingerprints")
    return fps

# Convert fingerprints to binary
def fp_to_binary(fps):
    fp_start_time = time.perf_counter()
    binary = [x.ToBitString() for x in fps]
    fp_end_time = time.perf_counter()
    print(f"It took {fp_end_time - fp_start_time} seconds to convert fingerprints to binary")

    results = dict(zip(SMILES, binary))
    #print(results)
    return results

def set_to_redis(results):
    fp_start_time = time.perf_counter()
    #[r.set(x, results[x]) for x in results]
    r.hmset("smiles", results)
    fp_end_time = time.perf_counter()
    print(f"It took {fp_end_time - fp_start_time} seconds to set the values in redis")

def get_from_redis():
    fp_start_time = time.perf_counter()
    #values_in_redis = [r.get(x) for x in SMILES]
    values_in_redis = r.hgetall("smiles")
    #print(f'{values_in_redis}')
    fp_end_time = time.perf_counter()
    print(f"It took {fp_end_time - fp_start_time} seconds to get the values in redis")
    return values_in_redis

def binary_to_vector(values_in_redis):
    fp_start_time = time.perf_counter()
    converted_binary = [DataStructs.cDataStructs.CreateFromBitString(values_in_redis[x]) for x in values_in_redis.keys()]
    fp_end_time = time.perf_counter()
    print(f"It took {fp_end_time - fp_start_time} seconds to get convert values to Bit Vectors")
    return converted_binary

def generate_fingerprint(smile_arg):
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3,fpSize=512)
    fp = mfpgen.GetFingerprint(Chem.MolFromSmiles(smile_arg))
    return fp

def tanimoto_search(fp, converted_binary, threashold, values_from_redis):
    fp_start_time = time.perf_counter()
    tanimoto_results = DataStructs.BulkTanimotoSimilarity(fp, converted_binary)
    tanimoto_results = [x for x in tanimoto_results]
    fp_end_time = time.perf_counter()
    print(f"It took {fp_end_time - fp_start_time} seconds to perform stucture sim search")

    fp_start_time = time.perf_counter()
    result_df = pd.DataFrame(values_from_redis.keys(), columns=["Smiles"])
    result_df["tanimoto_results"] = tanimoto_results

    result_df = result_df.loc[result_df['tanimoto_results'] >= threashold].sort_values(by="tanimoto_results", ascending=False)
    result_df.to_csv("Tanimoto_results.csv")
    fp_end_time = time.perf_counter()
    print(f"It took {fp_end_time - fp_start_time} seconds to build the dataframes")
    #print(result_df.head(25))
    return result_df.sort_values(by="tanimoto_results", ascending=True)
    

if __name__ == "__main__":
    smile_arg = 'C1CC(C)C1'
    smile_arg = 'O=C(O)CCCCC(=O)O'
    #tautomer
    smile_arg = 'Oc1c(cccc3)c3nc2ccncc12'
    threashold = 0.2
    #mols = generate_mol_files()
    #fps = generate_fingerprints(mols)
    #binary = fp_to_binary(fps)
    #set_to_redis(binary)
    values_from_redis = get_from_redis()
    fps2 = binary_to_vector(values_from_redis)
    example_fp = generate_fingerprint(smile_arg)
    result = tanimoto_search(example_fp, fps2, threashold, values_from_redis)