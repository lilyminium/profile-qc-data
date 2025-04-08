import collections
import tqdm
import pathlib
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
from openff.toolkit import ForceField


def main(
    dataset="../offqcdata/data/torsion-parameters",
    forcefield="openff_unconstrained-2.2.1.offxml",
    output_file: str = "../coverage/coupled-torsions.csv"
):
    dataset = ds.dataset(dataset)
    print(f"Loaded {dataset.count_rows()} rows from {dataset}")

    forcefield_name = pathlib.Path(forcefield).stem
    subset = dataset.filter(pc.field("forcefield") == forcefield_name)
    print(f"Filtered to {subset.count_rows()} rows with forcefield {forcefield_name}")

    df = subset.to_table().to_pandas()
    TORSION_PAIRS = collections.defaultdict(collections.Counter)
    for _, row in tqdm.tqdm(df.iterrows(), total=len(df)):
        parameter_ids = row["all_parameter_ids"]
        for i, parameter in enumerate(parameter_ids):
            remaining = parameter_ids[i + 1:]
            TORSION_PAIRS[parameter]["All"] += 1
            for pid in remaining:
                TORSION_PAIRS[parameter][pid] += 1
                TORSION_PAIRS[pid][parameter] += 1
    
    entries = []
    for parameter_id, values in TORSION_PAIRS.items():
        all_counts = values.pop("All")
        for counter_parameter, count in values.items():
            entries.append({
                "parameter_id": parameter_id,
                "other_parameter_id": counter_parameter,
                "all_count": all_counts,
                "count": count,
            })

    pair_df = pd.DataFrame(entries)
    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    pair_df.to_csv(output_file, index=False)

    # find same counts
    matches = pair_df[pair_df["count"] == pair_df["all_count"]]
    print("The following torsions are always trained with another:")
    print(sorted(matches.parameter_id.unique()))


if __name__ == "__main__":
    main()
