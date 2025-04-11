import functools
import multiprocessing
import click
import tqdm
import pandas as pd

from rdkit import Chem
from openff.toolkit import Molecule, ForceField

from yammbs import MoleculeStore
from yammbs.analysis import get_internal_coordinates

MAPPING = {
    "Bond": "Bonds",
    "Angle": "Angles",
    "Dihedral": "ProperTorsions",
    "Improper": "ImproperTorsions",
}

def get_single_labels(molecule_id, store=None, ff_name: str = None):
    ff = ForceField(ff_name)

    entries = []
    smiles = store.get_smiles_by_molecule_id(molecule_id)
    mol = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
    labels = ff.label_molecules(mol.to_topology())[0]
    qcarchive_ids = store.get_qcarchive_ids_by_molecule_id(molecule_id)

    for qcarchive_id in qcarchive_ids:
        qm_conformer = store.get_qm_conformer_by_qcarchive_id(qcarchive_id)
        mm_conformer = store.get_mm_conformer_by_qcarchive_id(
            qcarchive_id,
            ff_name
        )
        internal_coordinates = get_internal_coordinates(
            mol,
            qm_conformer,
            mm_conformer,
        )
        for geometric_type, geometric_values in internal_coordinates.items():
            openff_type = MAPPING[geometric_type]
            openff_values = labels[openff_type]
            for indices, (qm_value, mm_value) in geometric_values.items():
                if indices not in openff_values and geometric_type == "Improper":
                    continue
                parameter_id = openff_values[indices].id
                entry = {
                    "parameter_type": openff_type,
                    "parameter_id": parameter_id,
                    "qcarchive_id": qcarchive_id,
                    "qm_value": qm_value,
                    "mm_value": mm_value,
                }
                entries.append(entry)
    return entries


@click.command()
@click.option(
    "--input-file",
    "-i",
    type=click.Path(exists=True, dir_okay=False),
    default="../offqcdata/data/yammbs-unsafe/optimizations.sqlite",
    help="Path to the input SQLite file.",
)
@click.option(
    "--n-processes",
    "-np",
    type=int,
    default=1,
    help="Number of processes to use for parallel processing.",
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(exists=False, dir_okay=False),
    default="../coverage/valence-parameters.csv",
    help="Path to the output CSV file.",
)
def main(
    input_file="../offqcdata/data/yammbs-unsafe/optimizations.sqlite",
    n_processes: int = 1,
    output_file="../coverage/valence-parameters.csv",
):
    store = MoleculeStore(input_file)
    molecule_ids = sorted(store.get_molecule_ids())
    print(f"Loaded {len(molecule_ids)} molecules from {input_file}")

    labeller = functools.partial(
        get_single_labels,
        store=store,
        ff_name="openff_unconstrained-2.2.1.offxml",
    )

    with multiprocessing.Pool(n_processes) as pool:
        entries = list(
            tqdm.tqdm(
                pool.imap(
                    labeller,
                    molecule_ids,
                ),
                total=len(molecule_ids),
            )
        )

    df = pd.DataFrame(entries)
    df.to_csv(output_file, index=False)
    print(f"Saved valence parameters to {output_file}")


if __name__ == "__main__":
    main()
    