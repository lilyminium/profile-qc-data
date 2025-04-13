import click
import tqdm

import pyarrow.dataset as ds
import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import SimDivFilters


def bond_distance(i, j):
    return abs(i - j)

def angle_distance(i, j):
    """Return angle distance in radians
    
    Parameters
    ----------
    i : float
        Angle in radians
    j : float
        Angle in radians
    
    Returns
    -------
    float
        Angle distance in radians, wrapped to [-pi, pi]
    """
    return abs(i - j) % (2 * np.pi)


def select_diverse_parameters(parameter_id, subdf, n_records: int = 25):
    """
    Select optimizations with the most diverse qm_values
    """
    qm_values = subdf.qm_value.values
    qcarchive_ids = subdf.qcarchive_id.values

    mmp = SimDivFilters.MaxMinPicker()
    if parameter_id.startswith("b"):
        dist_func = bond_distance
    else:
        dist_func = angle_distance
    pick_size = min([n_records, len(qm_values)])
    picked_indices = list(
        mmp.LazyPick(
            lambda i, j: dist_func(qm_values[i], qm_values[j]),
            len(qm_values), # pool size
            pick_size, # pick size
        )
    )
    print(f"{parameter_id}: {qm_values[picked_indices]}")
    return qcarchive_ids[picked_indices], qm_values[picked_indices]




@click.command()
@click.option(
    "--input-directory",
    default="../coverage/valence-parameters",
    help="Directory containing the input dataset.",
)
@click.option(
    "--n-records",
    default=25,
    help="Number of records to select for each parameter.",
)
@click.option(
    "--output-file",
    default="../coverage/valence-parameter-diverse-optimizations.csv",
    help="Path to the output CSV file.",
)
def main(
    input_directory: str = "../coverage/valence-parameters",
    n_records: int = 25,
    output_file: str = "../coverage/valence-parameter-diverse-optimizations.csv",
):
    dataset = ds.dataset(input_directory)
    df = dataset.to_table().to_pandas()

    entries = []

    for (parameter_type, parameter_id), subdf in tqdm.tqdm(
        df.groupby(by=["parameter_type", "parameter_id"]),
        desc="Selecting optimizations"
    ):
        qcarchive_ids, qm_values = select_diverse_parameters(
            parameter_id,
            subdf,
            n_records=n_records
        )
        for qcarchive_id, qm_value in zip(qcarchive_ids, qm_values):
            entry = {
                "parameter_type": parameter_type,
                "parameter_id": parameter_id,
                "qcarchive_id": qcarchive_id,
                "qm_value": qm_value,
            }
            entries.append(entry)
    
    output_df = pd.DataFrame(entries)
    output_df = output_df.sort_values("parameter_id")
    output_df.to_csv(output_file, index=False)

    print(f"Wrote {len(entries)} entries to {output_file}")
    
        
if __name__ == "__main__":
    main()
