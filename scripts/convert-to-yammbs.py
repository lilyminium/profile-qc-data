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
    optimization_store = MoleculeStore.from_qm_dataset(
        dataset=overall_dataset,
        database_name=store_file,
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
