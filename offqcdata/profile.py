import logging

import collections
import pathlib
import tqdm
import warnings

import pyarrow as pa
import numpy as np
import pyarrow.compute as pc
import pyarrow.parquet as pq
import pyarrow.dataset as ds

from openff.toolkit import Molecule, ForceField
from yammbs.checkmol import analyze_functional_groups, ChemicalEnvironment
from rdkit import Chem


logger = logging.getLogger(__name__)

ELEMENTS = ["O", "N", "S", "P", "F", "Cl", "Br", "I"]
ENVIRONMENTS = sorted([
    x.value
    for x in ChemicalEnvironment
])

def analyze_smiles(smiles: str):
    groups_ = analyze_functional_groups(smiles)
    if groups_ is None:
        warnings.warn(
            f"Could not analyze SMILES: {smiles}",
            UserWarning,
        )
        return {
            "elements": [],
            "functional-groups": [],
        }
    groups = [
        x.value for x in groups_
    ]
    rdmol = Chem.MolFromSmiles(smiles)
    elements = set([
        atom.GetSymbol() for atom in rdmol.GetAtoms()
        if atom.GetSymbol() != "H"
    ])
    return {
        "elements": sorted(elements),
        "functional-groups": sorted(groups),
    }


def get_unique_smiles(input_file: str, column: str = "smiles") -> set[str]:
    # Load the table
    table = pq.read_table(input_file)
    logger.info(f"Loaded {table.num_rows} rows from {input_file}")

    # Get the unique smiles
    unique_smiles = set(pc.unique(table.column(column)))
    logger.info(f"Loaded {len(unique_smiles)} unique {column}")
    return unique_smiles


def analyze_table(input_file: str, output_directory: str):
    unique_smiles = get_unique_smiles(input_file)
    
    file_number = 0

    # Read existing data
    output_path = pathlib.Path(output_directory)
    if not output_path.exists():
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        existing_dataset = ds.dataset(output_path)
        if existing_dataset.count_rows():
            existing_smiles = set(
                existing_dataset.to_table(columns=["smiles"]).to_pydict()["smiles"]
            )
            logger.info(f"Loaded {len(existing_smiles)} existing smiles")
            unique_smiles = unique_smiles - existing_smiles
            logger.info(f"New smiles: {len(unique_smiles)}")
            file_number = len(existing_dataset.files)
    
    entries = []

    for smiles in unique_smiles:
        smiles = str(smiles)
        result = analyze_smiles(smiles)
        mol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        mw = sum([atom.mass for atom in mol.atoms]).m
        n_atoms = len(mol.atoms)
        n_heavy_atoms = sum([1 for atom in mol.atoms if atom.atomic_number != 1])
        entry = {
            "smiles": smiles,
            "mw": mw,
            "n_atoms": n_atoms,
            "n_heavy_atoms": n_heavy_atoms,
        }
        entry.update(dict.fromkeys(ELEMENTS, False))
        entry.update(dict.fromkeys(ENVIRONMENTS, False))
        for groups in result.values():
            for group in groups:
                entry[group] = True
        entries.append(entry)
    
    # Create a new table and save
    new_table = pa.Table.from_pylist(entries)
    new_filename = output_path / f"batch-{file_number:04d}.parquet"
    assert not new_filename.exists(), f"File {new_filename} already exists"
    pq.write_table(new_table, new_filename)
    logger.info(f"Wrote {len(entries)} rows to {new_filename}")


def label_table_with_forcefield(
    forcefield_file: str,
    input_file: str,
    output_directory: str
):
    forcefield = ForceField(forcefield_file)
    forcefield_name = pathlib.Path(forcefield_file).stem
    unique_smiles = get_unique_smiles(input_file, "cmiles")

    file_number = 0

    # Read existing data
    output_path = pathlib.Path(output_directory)
    if not output_path.exists():
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        existing_dataset = ds.dataset(output_path)
        if existing_dataset.count_rows():
            subset = existing_dataset.filter(
                pc.field("forcefield") == forcefield_name
            )
            existing_smiles = set(
                subset.to_table(columns=["cmiles"]).to_pydict()["cmiles"]
            )
            logger.info(f"Loaded {len(existing_smiles)} existing cmiles")
            unique_smiles = unique_smiles - existing_smiles
            logger.info(f"New cmiles: {len(unique_smiles)}")
            file_number = len(existing_dataset.files)
    
    entries = []

    for smiles in tqdm.tqdm(unique_smiles):
        mol = Molecule.from_mapped_smiles(
            str(smiles),
            allow_undefined_stereo=True,
        )
        labels = forcefield.label_molecules(mol.to_topology())[0]

        for parameter_type in ["Bonds", "Angles", "ProperTorsions", "ImproperTorsions"]:
            parameter_indices = collections.defaultdict(list)
            for indices, parameter in labels[parameter_type].items():
                parameter_indices[parameter.id].append(list(indices))
            for k, v in parameter_indices.items():
                entry = {
                    "cmiles": smiles,
                    "forcefield": forcefield_name,
                    "parameter_type": parameter_type,
                    "parameter_id": k,
                    "parameter_indices": np.array(v).flatten().tolist(),
                }
                entries.append(entry)
    
    # Create a new table and save
    new_table = pa.Table.from_pylist(entries)
    new_filename = output_path / f"batch-{file_number:04d}.parquet"
    assert not new_filename.exists(), f"File {new_filename} already exists"
    pq.write_table(new_table, new_filename)
    logger.info(f"Wrote {len(entries)} rows to {new_filename}")



# def update_groups_from_table(table: pa.Table, output_directory):
#     unique_smiles = set(pc.unique(table.column("smiles")))
#     logger.info(f"Loaded {len(unique_smiles)} unique smiles")

#     # load existing from all
#     output_path = pathlib.Path(output_directory)
#     all_path = output_path / "functional-groups" / "All"
#     all_path.mkdir(parents=True, exist_ok=True)

#     all_dataset = ds.dataset(all_path)
#     existing_smiles = set(
#         all_dataset.to_table(columns=["smiles"]).to_pydict()["smiles"]
#     )
#     logger.info(f"Loaded {len(existing_smiles)} existing smiles")
#     new_smiles = unique_smiles - existing_smiles
#     logger.info(f"New smiles: {len(new_smiles)}")
#     if not new_smiles:
#         logger.info("No new smiles to analyze")
#         return
    
#     elements = collections.defaultdict(list)
#     functional_groups = collections.defaultdict(list)
#     GROUPS = {
#         "elements": elements,
#         "functional-groups": functional_groups,
#     }

#     for smiles in tqdm.tqdm(
#         new_smiles,
#         total=len(new_smiles),
#         desc=f"Analyzing {smiles}",
#     ):
#         result = analyze_smiles(smiles)
#         for group, values in result.items():
#             for chemgroup in values:
#                 GROUPS[group][chemgroup].append(smiles)
    

#     # update parquet tables
#     for group, values in GROUPS.items():
#         group_path = all_path / group
#         group_path.mkdir(parents=True, exist_ok=True)

#         for chemgroup, added_smiles in values.items():
#             outfile = group_path / f"{chemgroup}.parquet"
#             if outfile.exists():
#                 logger.info(f"Appending to {outfile}")
#                 table = pq.read_table(outfile)
#                 smiles = set(table.column("smiles").to_pylist())
#                 smiles.update(added_smiles)
#             else:
#                 logger.info(f"Creating {outfile}")
#                 smiles = set(added_smiles)

#             table = 


#             # smiles = sorted(smiles)
#             # table = pa.Table.from_pydict({
#             #     "smiles": smiles,
#             #     "count": [len(smiles)] * len(smiles),
#             # })
#             # pq.write_table(
#             #     table,
#             #     group_path / f"{chemgroup}.parquet",
#             #     partition_cols=["smiles"],
#             #     use_compression="snappy",
#             # )