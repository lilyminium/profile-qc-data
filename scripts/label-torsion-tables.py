import pathlib
import tqdm

from offqcdata.profile import label_torsion_table_with_forcefield


def main(
    root_directory="../offqcdata",
    forcefield="openff_unconstrained-2.2.1.offxml",
):
    root_directory = pathlib.Path(root_directory)
    table_files = sorted(
        (root_directory / "data/tables/torsiondrive").glob("*.parquet")
    )

    output_directory = root_directory / "data/torsion-parameters"

    for table_file in tqdm.tqdm(table_files):
        label_torsion_table_with_forcefield(
            forcefield,
            table_file,
            output_directory,
        )


if __name__ == "__main__":
    main()
