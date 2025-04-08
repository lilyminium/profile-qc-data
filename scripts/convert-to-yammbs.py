import logging
import pathlib
import tqdm

from openff.qcsubmit.results import OptimizationResultCollection, TorsionDriveResultCollection
import qcportal as ptl
from yammbs.inputs import QCArchiveDataset
from yammbs.torsion.inputs import QCArchiveTorsionDataset
from yammbs import MoleculeStore
from yammbs.torsion import TorsionStore

QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

def create_moleculestore_verbose(dataset, dataset_name):
    from yammbs.models import MoleculeRecord, QMConformerRecord
    from yammbs._molecule import _smiles_to_inchi_key

    store = MoleculeStore(dataset_name)

    # only store unique molecule records
    input_molecule_records = []
    for qm_molecule in tqdm.tqdm(dataset.qm_molecules, desc="Creating MoleculeRecords"):
        molecule_record = MoleculeRecord(
            mapped_smiles=qm_molecule.mapped_smiles,
            inchi_key=_smiles_to_inchi_key(qm_molecule.mapped_smiles),
        )
        input_molecule_records.append(molecule_record)
    
    store.store(input_molecule_records)

    # now store QM records all at once
    input_qm_records = []
    for qm_molecule in tqdm.tqdm(dataset.qm_molecules, desc="Creating QMConformerRecords"):
        qm_conformer_record = QMConformerRecord(
            molecule_id=store.get_molecule_id_by_smiles(
                qm_molecule.mapped_smiles,
            ),
            qcarchive_id=qm_molecule.qcarchive_id,
            mapped_smiles=qm_molecule.mapped_smiles,
            coordinates=qm_molecule.coordinates,
            energy=qm_molecule.final_energy,
        )
        input_qm_records.append(qm_conformer_record)
    
    store.store_qcarchive(input_qm_records)


def main(
    root_directory=".."
):
    root_directory = pathlib.Path(root_directory)
    client = ptl.PortalClient(address=QCFRACTAL_URL, cache_dir=(root_directory / "_cache").resolve())

    optimization_file = root_directory / "offqcdata/data/yammbs/optimizations.json"
    collection = OptimizationResultCollection.parse_file(optimization_file)

    # convert to QCArchiveDataset by batching as calls to client can be really quite slow
    n_entries = len(collection.entries[QCFRACTAL_URL])
    batch_size = 1000
    n_batches = n_entries // batch_size + 1
    overall_dataset = QCArchiveDataset()
    for i in tqdm.tqdm(range(n_batches), desc="Converting optimizations to QCArchiveDataset"):
        start = i * batch_size
        end = min((i + 1) * batch_size, n_entries)
        entries = collection.entries[QCFRACTAL_URL][start:end]
        new_collection = OptimizationResultCollection(entries={QCFRACTAL_URL: entries})
        sub_qcarchive_dataset = QCArchiveDataset.from_qcsubmit_collection(new_collection)
        overall_dataset.qm_molecules.extend(sub_qcarchive_dataset.qm_molecules)
    
    store_file = root_directory / "offqcdata/data/yammbs/optimizations.sqlite"
    # optimization_store = MoleculeStore.from_qm_dataset(
    #     dataset=overall_dataset,
    #     database_name=store_file,
    # )
    create_moleculestore_verbose(
        dataset=overall_dataset,
        dataset_name=store_file,
    )
    logger.info(f"Wrote to {store_file}")

    # repeat for torsiondrives
    torsiondrive_file = root_directory / "offqcdata/data/yammbs/torsiondrives.json"
    collection = TorsionDriveResultCollection.parse_file(torsiondrive_file)
    n_entries = len(collection.entries[QCFRACTAL_URL])
    n_batches = n_entries // batch_size + 1
    overall_dataset = QCArchiveTorsionDataset()
    for i in tqdm.tqdm(range(n_batches), desc="Converting torsiondrives to QCArchiveTorsionDataset"):
        start = i * batch_size
        end = min((i + 1) * batch_size, n_entries)
        entries = collection.entries[QCFRACTAL_URL][start:end]
        new_collection = TorsionDriveResultCollection(entries={QCFRACTAL_URL:entries})
        sub_qcarchive_dataset = QCArchiveTorsionDataset.from_qcsubmit_collection(new_collection)
        overall_dataset.qm_torsions.extend(sub_qcarchive_dataset.qm_torsions)
    
    store_file = root_directory / "offqcdata/data/yammbs/torsiondrives.sqlite"
    torsiondrive_store = TorsionStore.from_qm_dataset(
        dataset=overall_dataset,
        database_name=store_file,
    )
    logger.info(f"Wrote to {store_file}")


        
if __name__ == "__main__":
    main()
