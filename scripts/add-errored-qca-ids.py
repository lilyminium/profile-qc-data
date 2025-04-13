import pathlib
import click
import pandas as pd

def main(
    input_file: str = "../coverage/valence-parameters_errors.csv",
    output_file: str = "../offqcdata/data/bad-qcarchive_ids.dat",
):
    df = pd.read_csv(input_file)
    qca_ids = df.qcarchive_id.unique()

    output_file = pathlib.Path(output_file)
    existing_bad_ids = set()
    if output_file.exists():
        with output_file.open("r") as f:
            existing_bad_ids = set([int(x.strip()) for x in f.readlines()])
    print(f"Loaded {len(existing_bad_ids)} existing bad ids from {output_file}")

    combined_bad_ids = existing_bad_ids.union(set(qca_ids))
    with output_file.open("w") as f:
        f.write("\n".join(map(str, sorted(combined_bad_ids))))
    print(f"Wrote {len(combined_bad_ids)} bad qcarchive ids to {output_file}")


if __name__ == "__main__":
    main()
