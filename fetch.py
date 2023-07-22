# provide methods to get 3D Conformer SDF files from PubChem
# https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial
import requests


def sdf_from_sid(sid) -> str:
    """Get 3D Conformer SDF file from PubChem by SID"""
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/sid/" + str(sid) + "/record/SDF/?record_type=3d"
    response = requests.get(url)
    return response.text


def sdf_from_cid(cid) -> str:
    """Get 3D Conformer SDF file from PubChem by CID"""
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + str(cid) + "/record/SDF/?record_type=3d"
    response = requests.get(url)
    return response.text


def sdf_from_name(name: str) -> str:
    """Get 3D Conformer SDF file from PubChem by name"""
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + name + "/record/SDF/?record_type=3d"
    response = requests.get(url)
    return response.text


def sdf_from_smiles(smiles: str) -> str:
    """Get 3D Conformer SDF file from PubChem by SMILES"""
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" + smiles + "/record/SDF/?record_type=3d"
    response = requests.get(url)
    return response.text


def save_sdf_to_file(sdf: str, filename: str) -> None:
    with open("molecule_sdf/" + filename, "w") as f:
        f.write(sdf)


if __name__ == "__main__":
    save_sdf_to_file(sdf_from_cid(2244), "2244.sdf")
    save_sdf_to_file(sdf_from_name("glucose"), "glucose.sdf")

