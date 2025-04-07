import qcportal as ptl
import tqdm
import pathlib
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.dataset as ds
from offqcdata.download import download_optimization, get_client, download_torsiondrive


IGNORE_IODINE = [
    "OpenFF Discrepancy Benchmark 1",
    "OpenFF Gen 2 Opt Set 2 Coverage",
    "OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy",
    "SMIRNOFF Coverage Set 1",
    "OpenFF Ehrman Informative Optimization v0.2",
    "FDA optimization dataset 1",
    "Kinase Inhibitors: WBO Distributions",

    # ---
    "OpenFF Gen 2 Torsion Set 2 Coverage 2",
    "OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2",
]

OPTIMIZATION_WHITELISTS = [
    "OpenFF Optimization Set 1",
    "SMIRNOFF Coverage Set 1",
    "OpenFF VEHICLe Set 1",
    "OpenFF Discrepancy Benchmark 1",
    "OpenFF Ehrman Informative Optimization v0.2",
    "Pfizer discrepancy optimization dataset 1",
    "FDA optimization dataset 1",
    "Kinase Inhibitors: WBO Distributions",
    "OpenFF Gen 2 Opt Set 1 Roche",
    "OpenFF Gen 2 Opt Set 2 Coverage",
    "OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy",
    "OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy",
    "OpenFF Gen 2 Opt Set 5 Bayer",
    "OpenFF Sandbox CHO PhAlkEthOH v1.0",
    "OpenFF Industry Benchmark Season 1 v1.1",
    "OpenFF Gen2 Optimization Dataset Protomers v1.0",
    "OpenFF Protein Capped 1-mers 3-mers Optimization Dataset v1.0",
    "OpenFF Iodine Chemistry Optimization Dataset v1.0",
    # "OpenFF multi-Br ESP Fragment Conformers v1.0",
    "XtalPi Shared Fragments OptimizationDataset v1.0",
    "XtalPi 20-percent Fragments OptimizationDataset v1.0",
    "OpenFF Torsion Benchmark Supplement v1.0",
    "OpenFF Torsion Multiplicity Optimization Training Coverage Supplement v1.0",
    "OpenFF Torsion Multiplicity Optimization Benchmarking Coverage Supplement v1.0",
    "OpenFF Iodine Fragment Opt v1.0",
    "OpenFF Sulfur Optimization Training Coverage Supplement v1.0",
    "OpenFF Sulfur Optimization Benchmarking Coverage Supplement v1.0",
    "OpenFF Lipid Optimization Training Supplement v1.0",
    "OpenFF Lipid Optimization Benchmark Supplement v1.0",
    # "SPICE DES370k Monomers Lowest E Conformer Optimization Dataset v4.0",
    "OpenFF Cresset Additional Coverage Optimizations v4.0",
    "OpenFF Protein PDB 4-mers v4.0"
]

TORSIONDRIVE_WHITELISTS = [
    "OpenFF Group1 Torsions",
    "SMIRNOFF Coverage Torsion Set 1",
    "OpenFF Substituted Phenyl Set 1",
    "Pfizer discrepancy torsion dataset 1",
    "OpenFF Primary Benchmark 1 Torsion Set",
    "OpenFF Gen 2 Torsion Set 2 Coverage 2",
    "OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2",
    "OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2",
    "OpenFF Gen 2 Torsion Set 5 Bayer 2",
    "OpenFF Gen 2 Torsion Set 6 supplemental 2",
    "OpenFF Fragmenter Validation 1.0",
    "OpenFF DANCE 1 eMolecules t142 v1.0",
    "OpenFF Rowley Biaryl v1.0",
    "OpenFF-benchmark-ligand-fragments-v1.0",
    "OpenFF Protein Fragments TorsionDrives v1.0",
    "OpenFF WBO Conjugated Series v1.0",
    "OpenFF Amide Torsion Set v1.0",
    # "OpenFF Aniline Para Opt v1.0",
    "OpenFF Gen3 Torsion Set v1.0",
    "OpenFF Aniline 2D Impropers v1.0",
    "OpenFF-benchmark-ligand-fragments-v2.0",
    "OpenFF multiplicity correction torsion drive data v1.1",
    "OpenFF Protein Capped 3-mer Omega v1.0",
    "XtalPi Shared Fragments TorsiondriveDataset v1.0",
    "OpenFF Torsion Coverage Supplement v1.0",
    "OpenFF RNA Dinucleoside Monophosphate TorsionDrives v1.0",
    "XtalPi 20-percent Fragments TorsiondriveDataset v1.0",
    "OpenFF Torsion Drive Supplement v1.0",
    "OpenFF Torsion Multiplicity Torsion Drive Coverage Supplement v1.0",
    "OpenFF Phosphate Torsion Drives v1.0",
    "OpenFF Alkane Torsion Drives v1.0",
    "OpenFF Cresset Additional Coverage TorsionDrives v4.0"
]

def main(
    root_directory=".."
):
    root_directory = pathlib.Path(root_directory)
    client = get_client((root_directory / "_cache").resolve())

    for dsname in tqdm.tqdm(OPTIMIZATION_WHITELISTS):
        table = download_optimization(client, dsname)
        if dsname in IGNORE_IODINE:
            df = table.to_pandas()
            n_df = len(df)
            mask = np.array(["I" in smi for smi in df.smiles.values])
            print(df[mask].smiles)
            df = pd.DataFrame(df[~mask])
            print(f"Found {len(mask)} in {dsname} -- filtering from {n_df} to {len(df)}")
            table = pa.Table.from_pandas(df)

        pq.write_table(
            table,
            root_directory / f"offqcdata/data/tables/optimization/{dsname}.parquet"
        )

    for dsname in tqdm.tqdm(TORSIONDRIVE_WHITELISTS):
        table = download_torsiondrive(client, dsname)
        if dsname in IGNORE_IODINE:
            df = table.to_pandas()
            n_df = len(df)
            mask = np.array(["I" in smi for smi in df.smiles.values])
            print(df[mask].smiles)
            df = pd.DataFrame(df[~mask])
            print(f"Found {len(mask)} in {dsname} -- filtering from {n_df} to {len(df)}")
            table = pa.Table.from_pandas(df)

        pq.write_table(
            table,
            root_directory / f"offqcdata/data/tables/torsiondrive/{dsname}.parquet"
        )

        
if __name__ == "__main__":
    main()
