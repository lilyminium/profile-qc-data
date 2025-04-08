import tqdm
import pathlib
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
from openff.toolkit import ForceField


def main(
    dataset="../offqcdata/data/torsion-parameters",
    count_threshold: int = 10,
    forcefield="openff_unconstrained-2.2.1.offxml",
    output_file: str = "../coverage/low-coverage-torsions.csv"
):
    dataset = ds.dataset(dataset)
    print(f"Loaded {dataset.count_rows()} rows from {dataset}")

    forcefield_name = pathlib.Path(forcefield).stem
    subset = dataset.filter(pc.field("forcefield") == forcefield_name)
    print(f"Filtered to {subset.count_rows()} rows with forcefield {forcefield_name}")

    df = subset.to_table().to_pandas()
    unique_counts_per_parameter_id = {}
    for k, subdf in df.groupby("parameter_id"):
        unique_counts_per_parameter_id[k] = len(subdf.cmiles.unique())
    count_df = pd.DataFrame.from_dict(
        unique_counts_per_parameter_id, orient="index",
        columns=["Count"]
    )
    subset = count_df[count_df["Count"] < count_threshold].sort_values("Count")
    print(f"Found {len(subset)} torsions with less than {count_threshold} unique conformers")

    forcefield = ForceField(forcefield)
    handler = forcefield.get_parameter_handler("ProperTorsions")
    smiles = []
    for parameter_id in tqdm.tqdm(subset.index):
        parameter = handler.get_parameter({"id": parameter_id})[0]
        smiles.append(parameter.smirks)
    subset["smirks"] = smiles

    output_file = pathlib.Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    subset.to_csv(output_file)
    print(f"Saved low-coverage torsions to {output_file}")
    

if __name__ == "__main__":
    main()
