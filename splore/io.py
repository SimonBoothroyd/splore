import functools
import gzip
import os.path
from typing import Generator, Literal, overload

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdmolfiles
from rdkit.Chem.Draw import rdMolDraw2D

IMAGE_UNAVAILABLE_SVG = """
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" style="isolation:isolate" viewBox="0 0 200 200" width="200pt" height="200pt"><defs><clipPath id="_clipPath_eSdCSpw1sB1xWp7flmMoZ0WjTPwPpzQh"><rect width="200" height="200"/></clipPath></defs><g clip-path="url(#_clipPath_eSdCSpw1sB1xWp7flmMoZ0WjTPwPpzQh)"><g clip-path="url(#_clipPath_LvpdWbrYj1cREqoXz8Lwbk3ZilfC6tg9)"><text transform="matrix(1,0,0,1,44.039,91.211)" style="font-family:'Open Sans';font-weight:400;font-size:30px;font-style:normal;fill:#000000;stroke:none;">Preview</text><text transform="matrix(1,0,0,1,17.342,132.065)" style="font-family:'Open Sans';font-weight:400;font-size:30px;font-style:normal;fill:#000000;stroke:none;">Unavailable</text></g><defs><clipPath id="_clipPath_LvpdWbrYj1cREqoXz8Lwbk3ZilfC6tg9"><rect x="0" y="0" width="166" height="81.709" transform="matrix(1,0,0,1,17,59.146)"/></clipPath></defs></g></svg>
"""


@overload
def molecules_from_file(
    file_path: str, as_smiles: Literal[True] = True
) -> Generator[str, None, None]:
    ...


@overload
def molecules_from_file(
    file_path: str, as_smiles: Literal[False] = False
) -> Generator[Chem.Mol, None, None]:
    ...


def molecules_from_file(file_path: str, as_smiles=False):

    from rdkit import Chem

    extension = os.path.splitext(file_path)[-1].lower()

    if extension in (".smi", ".txt"):
        supplier = Chem.SmilesMolSupplier(
            file_path, delimiter="\n", titleLine=False, nameColumn=-1
        )
    elif extension == ".sdf":
        supplier = Chem.SDMolSupplier(file_path)
    elif extension == ".gz":
        supplier = Chem.ForwardSDMolSupplier(gzip.open(file_path))
    else:
        raise NotImplementedError(f"{extension} not supported")

    for rd_molecule in supplier:

        if rd_molecule is None:
            continue

        yield Chem.MolToSmiles(rd_molecule) if as_smiles else rd_molecule


@functools.lru_cache(maxsize=4096)
def molecule_to_svg(smiles: str) -> str:

    smiles_parser = rdmolfiles.SmilesParserParams()
    smiles_parser.removeHs = True

    molecule: Chem.Mol = Chem.MolFromSmiles(smiles, smiles_parser)

    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)

    rdDepictor.Compute2DCoords(molecule)

    drawer = rdMolDraw2D.MolDraw2DSVG(200, 200, 150, 200)
    drawer.SetOffset(25, 0)
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()

    svg_content = drawer.GetDrawingText().replace(
        "<rect style='opacity:1.0", "<rect style='opacity: 0"
    )
    return svg_content
