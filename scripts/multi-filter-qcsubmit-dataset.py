import click
import pickle
import tqdm
import multiprocessing
import pathlib
import numpy as np

from qcportal import PortalClient
from openff.qcsubmit.results import OptimizationResultCollection
from openff.qcsubmit.results.filters import (
    ConnectivityFilter,
    RecordStatusFilter,
    UnperceivableStereoFilter,
    ElementFilter,
    ConformerRMSDFilter,
    
)
from openff.qcsubmit.utils import portal_client_manager


QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"

def filter_single_record(item):
    record, molecule = item
    if record.status.value.upper() != "COMPLETE":
        return None
    connectivity_filter = ConnectivityFilter(tolerance=1.2)
    connectivity = connectivity_filter._filter_function(None, None, molecule)
    if not connectivity:
        return None
    stereo_filter = UnperceivableStereoFilter()
    stereo = stereo_filter._filter_function(None, None, molecule)
    if not stereo:
        return None
    
    # check for zero coordinates
    z_coordinates = molecule.conformers[0].m[:, 2]
    if np.all(z_coordinates == 0):
        return None
    return record.id

    

@click.command()
@click.option(
    "--input-file",
    "-i",
    default="../offqcdata/data/yammbs/optimizations.json",
    help="Path to the input file containing the dataset.",
)
@click.option(
    "--output-file",
    "-o",
    default="../offqcdata/data/bad-qcarchive_ids.dat",
    help="Path to the output file for bad qcarchive ids.",
)
@click.option(
    "--n-processes",
    "-np",
    default=1,
    help="Number of processes to use for filtering.",
)
def main(
    input_file: str = "../offqcdata/data/yammbs/optimizations.json",
    output_file: str = "../offqcdata/data/bad-qcarchive_ids.dat",
    n_processes: int = 1,
):
    dataset = OptimizationResultCollection.parse_file(input_file)
    
    input_ids = set([
        entry.record_id
        for entry in dataset.entries[QCFRACTAL_URL]
    ])
    print(f"Loaded {len(input_ids)} records from {input_file}")

    output_file = pathlib.Path(output_file)
    existing_bad_ids = set()
    if output_file.exists():
        with output_file.open("r") as f:
            existing_bad_ids = set([int(x.strip()) for x in f.readlines()])

    print(f"Loaded {len(existing_bad_ids)} existing bad ids from {output_file}")

    pickle_file = pathlib.Path("records_and_molecules.pkl")
    if not pickle_file.exists():

        with portal_client_manager(lambda x: PortalClient(x, cache_dir="../_cache")):
            records_and_molecules = list(dataset.to_records())

        with open("records_and_molecules.pkl", "wb") as f:
            pickle.dump(records_and_molecules, f)

    else:
        with open("records_and_molecules.pkl", "rb") as f:
            records_and_molecules = pickle.load(f)

    print(f"Loaded {len(records_and_molecules)} records and molecules")

    with multiprocessing.Pool(n_processes) as pool:
        all_ids = [
            record_id
            for record_id in tqdm.tqdm(
                pool.imap(
                    filter_single_record,
                    records_and_molecules
                ),
                total=len(records_and_molecules)
            )
        ]

    clean_ids = [
        record_id
        for record_id in all_ids
        if record_id is not None
    ]
    removed_ids = input_ids - set(clean_ids)

    existing_bad_ids |= removed_ids
    with output_file.open("w") as f:
        f.write("\n".join(map(str, sorted(existing_bad_ids))))
    print(f"Wrote {len(existing_bad_ids)} bad qcarchive ids to {output_file}")



if __name__ == "__main__":
    main()
