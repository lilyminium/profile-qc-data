import logging
import pickle
import pathlib
import tqdm

from openff.qcsubmit.results import OptimizationResultCollection, TorsionDriveResultCollection
import qcportal as ptl
from yammbs.inputs import QCArchiveDataset
from yammbs.torsion.inputs import QCArchiveTorsionDataset
from yammbs import MoleculeStore
from yammbs.torsion import TorsionStore
from yammbs._db import DBQMConformerRecord, DBMoleculeRecord

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
    seen_mapped_smiles = set()
    for qm_molecule in tqdm.tqdm(dataset.qm_molecules, desc="Creating MoleculeRecords"):
        if qm_molecule.mapped_smiles in seen_mapped_smiles:
            continue
        molecule_record = MoleculeRecord(
            mapped_smiles=qm_molecule.mapped_smiles,
            inchi_key=_smiles_to_inchi_key(qm_molecule.mapped_smiles),
        )
        input_molecule_records.append(molecule_record)
        seen_mapped_smiles.add(qm_molecule.mapped_smiles)

    print(f"Created {len(input_molecule_records)} MoleculeRecords")

    with store._get_session() as db:
        db_records = [
            DBMoleculeRecord(
                mapped_smiles=record.mapped_smiles,
                inchi_key=record.inchi_key,
            )
            for record in input_molecule_records
        ]
        db.db.add_all(db_records)
    
    # now store QM records all at once
    input_qm_records = []
    smiles_to_molecule_id = {
        smi: store.get_molecule_id_by_smiles(smi)
        for smi in seen_mapped_smiles
    }
    seen_qcarchive_ids = set()
    for qm_molecule in tqdm.tqdm(dataset.qm_molecules, desc="Creating QMConformerRecords"):
        if qm_molecule.qcarchive_id in seen_qcarchive_ids:
            continue
        qm_conformer_record = QMConformerRecord(
            molecule_id=smiles_to_molecule_id[qm_molecule.mapped_smiles],
            qcarchive_id=qm_molecule.qcarchive_id,
            mapped_smiles=qm_molecule.mapped_smiles,
            coordinates=qm_molecule.coordinates,
            energy=qm_molecule.final_energy,
        )
        input_qm_records.append(qm_conformer_record)
        seen_qcarchive_ids.add(qm_molecule.qcarchive_id)

    print(f"Created {len(input_qm_records)} QMConformerRecords")

    with store._get_session() as db:
        db_records = [
            DBQMConformerRecord(
                parent_id=record.molecule_id,
                qcarchive_id=record.qcarchive_id,
                mapped_smiles=record.mapped_smiles,
                coordinates=record.coordinates,
                energy=record.energy,
            )
            for record in input_qm_records
        ]
        db.db.add_all(db_records)




def main(
    root_directory=".."
):
    root_directory = pathlib.Path(root_directory)
    with open("optimizations.pkl", "rb") as f:
        overall_dataset = pickle.load(f)
    
    store_file = root_directory / "offqcdata/data/yammbs-unsafe/optimizations.sqlite"
    store_file.parent.mkdir(parents=True, exist_ok=True)

    create_moleculestore_verbose(
        dataset=overall_dataset,
        dataset_name=store_file,
    )
    logger.info(f"Wrote to {store_file}")


        
if __name__ == "__main__":
    main()
