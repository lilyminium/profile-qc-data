import click
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


def main(
    input_file: str = "../offqcdata/data/yammbs/optimizations.json",
    output_file: str = "../offqcdata/data/bad-qcarchive_ids.dat",
):
    dataset = OptimizationResultCollection.parse_file(input_file)
    
    input_ids = set([
        entry.qcarchive_id
        for entry in dataset.entries[QCFRACTAL_URL]
    ])

    output_file = pathlib.Path(output_file)
    existing_bad_ids = set()
    if output_file.exists():
        with output_file.open("r") as f:
            existing_bad_ids = set([x.strip() for x in f.readlines()])

    clean_dataset = dataset.filter(
        RecordStatusFilter(status=RecordStatusEnum.complete),
        ConnectivityFilter(tolerance=1.2),
        UnperceivableStereoFilter(),
    )

    clean_ids = set([
        entry.qcarchive_id
        for entry in clean_dataset.entries[QCFRACTAL_URL]
    ])
    removed_ids = input_ids - clean_ids

    existing_bad_ids |= removed_ids
    with output_file.open("w") as f:
        f.write(sorted(existing_bad_ids))
    print(f"Wrote {len(existing_bad_ids)} bad qcarchive ids to {output_file}")



if __name__ == "__main__":
    main()
