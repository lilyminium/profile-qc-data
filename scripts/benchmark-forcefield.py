import pathlib
from yammbs import MoleculeStore
from yammbs.torsion import TorsionStore

def main(
    root_directory: str = "..",
    forcefield_file: str = "openff_unconstrained-2.2.1.offxml",
    n_processes: int = 16
):
    root_directory = pathlib.Path(root_directory)
    opt_store = MoleculeStore(
        root_directory / "offqcdata/data/yammbs/optimizations.sqlite",
    )
    opt_store.optimize_mm(
        forcefield_file,
        n_processes=n_processes,
    )
    
    #torsion_store = TorsionStore(
    #    root_directory / "offqcdata/data/yammbs/torsiondrives.sqlite",
    #)
    #torsion_store.optimize_mm(
    #    forcefield_file,
    #    n_processes=n_processes,
    #)


        
if __name__ == "__main__":
    main()
