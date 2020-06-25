# Copyright 2020 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""This module contains two molecular drawing functions to render SVGs or PNGs
given an RDKit molecule object: mol_to_svg and mol_to_png."""

import io
import base64
import numpy as np
from PIL import Image, ImageOps

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import rdDepictor
    from rdkit import __version__ as RDKIT_VERSION
    rdDepictor.SetPreferCoordGen(True)
except ImportError:
    Chem = None
    RDKIT_VERSION = None

# pylint: disable=unsubscriptable-object


def trim_image_whitespace(img, padding=0):
    """This function takes a PIL image, img, and crops it to the minimum
    rectangle based on its whiteness/transparency. 5 pixel padding used
    automatically."""

    # Convert to array
    as_array = np.array(img)  # N x N x (r,g,b,a)

    # Set previously-transparent pixels to white
    if as_array.shape[2] == 4:
        as_array[as_array[:, :, 3] == 0] = [255, 255, 255, 0]

    as_array = as_array[:, :, :3]

    # Content defined as non-white and non-transparent pixel
    has_content = np.sum(as_array, axis=2, dtype=np.uint32) != 255 * 3
    xs_nonzero, ys_nonzero = np.nonzero(has_content)

    # Crop down
    margin = 5
    x_range = (max([min(xs_nonzero) - margin,
                    0]), min([max(xs_nonzero) + margin, as_array.shape[0]]))
    y_range = (max([min(ys_nonzero) - margin,
                    0]), min([max(ys_nonzero) + margin, as_array.shape[1]]))
    as_array_cropped = as_array[x_range[0]:x_range[1], y_range[0]:y_range[1],
                                0:3]

    img = Image.fromarray(as_array_cropped, mode='RGB')
    return ImageOps.expand(img, border=padding, fill=(255, 255, 255))


def mol_to_svg(mol, max_size=300, bond_length=25):
    """Returns a string with an <svg> tag containing a rendering of the input
    molecule. Note that the SVG is not cropped."""
    rdDepictor.Compute2DCoords(mol)
    Chem.Kekulize(mol)
    drawer = Draw.MolDraw2DSVG(max_size, max_size)
    drawer.drawOptions().fixedBondLength = bond_length
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def mol_to_png(mol, max_size=1000, bond_length=25, png_quality=70):
    """Returns a string with an <img> tag containing a rendering of the input
    molecule."""
    drawer = Draw.MolDraw2DCairo(max_size, max_size)
    drawer.drawOptions().fixedBondLength = bond_length
    try:
        drawer.DrawMolecule(mol)
    except ValueError as value_error:
        raise ValueError(Chem.MolToSmiles(mol)) from value_error
    drawer.FinishDrawing()
    temp = io.BytesIO()
    temp.write(drawer.GetDrawingText())
    temp.seek(0)
    img = Image.open(temp)
    img = trim_image_whitespace(img, padding=10)
    output = io.BytesIO()
    img.save(output, format='png', quality=png_quality)
    output.seek(0)
    b64 = base64.b64encode(output.getvalue())
    return b64.decode("UTF-8")
