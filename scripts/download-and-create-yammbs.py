import pathlib
import tqdm

from openff.qcsubmit.results import OptimizationResultCollection, TorsionDriveResultCollection
import qcportal as ptl
from yammbs import MoleculeStore
from yammbs.torsion import TorsionStore

QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"

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
    # "OpenFF Gen3 Torsion Set v1.0",
    # "OpenFF Aniline 2D Impropers v1.0",
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

def download_and_convert(
    client,
    dataset_names,
    qcsubmit_class,
    yammbs_class,
    output_file,
    root_directory=".."
):
    root_directory = pathlib.Path(root_directory)
    

    collections = []

    #for dsname in tqdm.tqdm(OPTIMIZATION_WHITELISTS):
    for dsname in tqdm.tqdm(dataset_names):
        try:
            qcsubmit_collection = qcsubmit_class.from_server(
                client=client,
                datasets=dsname,
                spec_name="default"
            )
        except KeyError:
            qcsubmit_collection = qcsubmit_class.from_server(
                client=client,
                datasets=dsname,
                spec_name="spec_1"
            )

        # check for iodine
        if dsname in IGNORE_IODINE:
            entries = qcsubmit_collection.entries[QCFRACTAL_URL]
            not_i_entries = [
                entry
                for entry in entries
                if "I" not in entry.cmiles
            ]
            qcsubmit_collection.entries[QCFRACTAL_URL] = not_i_entries
            print(f"Filtered {len(entries)} entries to {len(not_i_entries)} entries from {dsname}")

        collections.append(qcsubmit_collection)
    
    final_collection = collections.pop(0)
    seen_ids = set([
        entry.record_id
        for entry in final_collection.entries[QCFRACTAL_URL]
    ])
    for collection in collections:
        entries_to_add = [
            entry
            for entry in collection.entries[QCFRACTAL_URL]
            if entry.record_id not in seen_ids
        ]
        final_collection.entries[QCFRACTAL_URL].extend(entries_to_add)
        seen_ids |= set([
            entry.record_id
            for entry in entries_to_add
        ])
    
    print(f"Final collection has {len(final_collection.entries[QCFRACTAL_URL])} entries")
    output_json = output_file.with_suffix(".json")
    with output_json.open("w") as f:
        f.write(final_collection.json(indent=4))

    # convert to MoleculeStore or TorsionStore
    store = yammbs_class.from_qcsubmit_collection(
        final_collection,
        output_file
    )
    print(f"Wrote to {output_file}")


def main(
    root_directory=".."
):
    root_directory = pathlib.Path(root_directory)
    client = ptl.PortalClient(address=QCFRACTAL_URL, cache_dir=(root_directory / "_cache").resolve())

    download_and_convert(
        client=client,
        dataset_names=OPTIMIZATION_WHITELISTS,
        qcsubmit_class=OptimizationResultCollection,
        yammbs_class=MoleculeStore,
        output_file=root_directory / "offqcdata/data/yammbs/optimizations.sqlite",
        root_directory=root_directory
    )

    download_and_convert(
        client=client,
        dataset_names=TORSIONDRIVE_WHITELISTS,
        qcsubmit_class=TorsionDriveResultCollection,
        yammbs_class=TorsionStore,
        output_file=root_directory / "offqcdata/data/yammbs/torsiondrives.sqlite",
        root_directory=root_directory
    )

        
if __name__ == "__main__":
    main()
