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
from yammbs.torsion.models import TorsionRecord, QMTorsionPointRecord
from yammbs._db import DBQMConformerRecord, DBMoleculeRecord
from yammbs.torsion._db import DBMMTorsionPointRecord, DBQMTorsionPointRecord, DBTorsionRecord


QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

def create_torsionstore_verbose(dataset: QCArchiveTorsionDataset, dataset_name: str):
    from yammbs.models import MoleculeRecord, QMConformerRecord
    from yammbs._molecule import _smiles_to_inchi_key

    store = TorsionStore(dataset_name)

    # only store unique torsion records
    input_torsion_records = []
    seen_torsions = set()
    for qm_torsion in tqdm.tqdm(dataset.qm_torsions, desc="Creating TorsionRecords"):
        key = (qm_torsion.mapped_smiles, qm_torsion.dihedral_indices)
        if key in seen_torsions:
            continue
        torsion_record = TorsionRecord(
            mapped_smiles=qm_torsion.mapped_smiles,
            inchi_key=_smiles_to_inchi_key(qm_torsion.mapped_smiles),
            dihedral_indices=qm_torsion.dihedral_indices,
        )
        input_torsion_records.append(torsion_record)
        seen_torsions.add(key)

    print(f"Created {len(input_torsion_records)} TorsionRecords")
    
    with store._get_session() as db:
        db_records = [
            DBTorsionRecord(
                mapped_smiles=record.mapped_smiles,
                inchi_key=record.inchi_key,
                dihedral_indices=record.dihedral_indices,
            )
            for record in input_torsion_records
        ]
        db.db.add_all(db_records)

    smiles_and_dihedrals_to_molecule_ids = {
        k: store.get_molecule_id_by_smiles_and_dihedral_indices(
            smiles=k[0],
            dihedral_indices=k[1]
        )
        for k in seen_torsions
    }

    # now do individual angles
    input_torsion_point_records = []
    for qm_torsion in dataset.qm_torsions:
        molecule_id = smiles_and_dihedrals_to_molecule_ids[
            (torsion_record.mapped_smiles, torsion_record.dihedral_indices)
        ]

        for angle in qm_torsion.coordinates:
            qm_point_record = QMTorsionPointRecord(
                molecule_id=molecule_id,
                grid_id=angle,  # TODO: This needs to be a tuple later
                coordinates=qm_torsion.coordinates[angle],
                energy=qm_torsion.energies[angle],
            )
            input_torsion_point_records.append(qm_point_record)

    print(f"Created {len(input_torsion_point_records)} QMTorsionPointRecords")

    with store._get_session() as db:
        db_records = [
            DBQMTorsionPointRecord(
                parent_id=record.molecule_id,
                grid_id=record.grid_id,
                coordinates=record.coordinates,
                energy=record.energy,
            )
            for record in input_torsion_point_records
        ]
        db.db.add_all(db_records)




def main(
    root_directory=".."
):
    root_directory = pathlib.Path(root_directory)
    with open("torsiondrives.pkl", "rb") as f:
        overall_dataset = pickle.load(f)
    
    store_file = root_directory / "offqcdata/data/yammbs-unsafe/torsiondrives.sqlite"
    store_file.parent.mkdir(parents=True, exist_ok=True)

    create_torsionstore_verbose(
        dataset=overall_dataset,
        dataset_name=store_file,
    )
    logger.info(f"Wrote to {store_file}")


        
if __name__ == "__main__":
    main()
