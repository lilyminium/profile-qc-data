import qcportal as ptl
import pyarrow as pa

from rdkit import Chem

QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"

def sanitize_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)

def get_client(cache_directory: str = None):
    client = ptl.PortalClient(address=QCFRACTAL_URL, cache_dir=cache_directory)
    return client


def get_canonical_smiles(entry) -> str:
    KEY = "canonical_isomeric_explicit_hydrogen_mapped_smiles"
    try:
        return entry.attributes[KEY]
    except KeyError:
        return getattr(entry.initial_molecule.identifiers, KEY)


def download_optimization(client, dataset_name: str):
    dataset = client.get_dataset("optimization", dataset_name)
    entries = {
        entry.name: entry
        for entry in dataset.iterate_entries()
    }
    data_entries = []
    for name, spec, record in dataset.iterate_records():
        if spec not in {"default", "spec_1"}:
            continue
        entry = entries[name]
        cmiles = get_canonical_smiles(entry)
        try:
            smiles = sanitize_smiles(cmiles)
        except ValueError as e:
            print(e)
            continue
        data_entry = {
            "id": record.id,
            "cmiles": cmiles,
            "smiles": smiles,
            "dataset_name": dataset_name
        }
        data_entries.append(data_entry)

    table = pa.Table.from_pylist(data_entries)
    return table


def download_torsiondrive(client, dataset_name: str):
    dataset = client.get_dataset("torsiondrive", dataset_name)
    entries = {
        entry.name: entry
        for entry in dataset.iterate_entries()
    }
    data_entries = []
    for name, spec, record in dataset.iterate_records():
        if spec not in {"default", "spec_1"}:
            continue
        entry = entries[name]
        cmiles = get_canonical_smiles(entry)
        try:
            smiles = sanitize_smiles(cmiles)
        except ValueError as e:
            print(e)
            continue
        dihedrals = record.specification.keywords.dihedrals
        # assert len(dihedrals) == 1, f"Expected 1 dihedral, got {len(dihedrals)}"
        # dihedral = dihedrals[0]
        dihedral = []
        for dih in dihedrals:
            dihedral.extend(dih)
        data_entry = {
            "id": record.id,
            "cmiles": cmiles,
            "smiles": smiles,
            "dataset_name": dataset_name,
            "dihedral": list(dihedral),
            "n_dihedrals": len(dihedrals),
        }
        data_entries.append(data_entry)

    table = pa.Table.from_pylist(data_entries)
    return table