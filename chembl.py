import requests
import pandas as pd

def download_chembl32():
  """Downloads the CHEMBL32 dataset and transforms the data into SMILES."""

  # Download the dataset
  url = "https://ftp.ebi.ac.uk/pub/databases/chembl/releases/CHEMBL32/chembl32.sdf"
  response = requests.get(url)
  with open("chembl32.sdf", "wb") as f:
    f.write(response.content)

  # Convert the dataset to a Pandas DataFrame
  df = pd.read_sdf("chembl32.sdf")

  # Extract the SMILES strings
  smiles = df["Molecule"].apply(lambda x: x.to_smiles())

  # Save the SMILES strings to a file
  with open("chembl32.smiles", "w") as f:
    f.write("\n".join(smiles))

if __name__ == "__main__":
  download_chembl32()