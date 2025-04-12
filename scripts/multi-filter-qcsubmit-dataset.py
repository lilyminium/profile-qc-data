import click
import tqdm
import multiprocessing
import pathlib

from openff.qcsubmit.results import OptimizationResultCollection
from qcportal.models.records import RecordStatusEnum
from openff.qcsubmit.results.filters import (
    ConnectivityFilter,
    RecordStatusFilter,
    UnperceivableStereoFilter,
    ElementFilter,
    ConformerRMSDFilter,
    
)


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
    return record.record_id

    

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

    output_file = pathlib.Path(output_file)
    existing_bad_ids = set()
    if output_file.exists():
        with output_file.open("r") as f:
            existing_bad_ids = set([x.strip() for x in f.readlines()])

    records_and_molecules = list(dataset.to_records())

    with multiprocessing.Pool(n_processes) as pool:
        clean_ids = [
            record_id
            for record_id in tqdm.tqdm(
                pool.imap(
                    filter_single_record,
                    records_and_molecules
                ),
                total=len(records_and_molecules)
            )
            if record_id is not None
        ]

    removed_ids = input_ids - clean_ids

    existing_bad_ids |= removed_ids
    with output_file.open("w") as f:
        f.write(sorted(existing_bad_ids))
    print(f"Wrote {len(existing_bad_ids)} bad qcarchive ids to {output_file}")



if __name__ == "__main__":
    main()
