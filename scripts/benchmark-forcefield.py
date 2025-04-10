import click
import pathlib
from yammbs import MoleculeStore
from yammbs.torsion import TorsionStore

@click.command()
@click.option(
    "--input",
    "-i",
    "input_file",
    type=click.Path(file_okay=True),
)
@click.option(
    "--forcefield",
    "-ff",
    "forcefield_file",
    type=str,
    default="openff_unconstrained-2.2.1.offxml",
)
@click.option(
    "--n-processes",
    "-np",
    type=int,
    default=16,
)
def main(
    input_file: str,
    forcefield_file: str = "openff_unconstrained-2.2.1.offxml",
    n_processes: int = 16
):
    if "optimization" in input_file:
        klass = MoleculeStore
    else:
        klass = TorsionStore

    store = klass(input_file)

    print(f"Benchmarking {input_file} with {forcefield_file}")

    store.optimize_mm(
        forcefield_file,
        n_processes=n_processes,
)
    print("Done!")
 


        
if __name__ == "__main__":
    main()
