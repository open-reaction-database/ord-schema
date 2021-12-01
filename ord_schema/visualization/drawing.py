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
"""Drawing functions.

This module contains two molecular drawing functions to render SVGs or PNGs
given an RDKit molecule object: mol_to_svg and mol_to_png.
"""

import io
import base64
import re
from typing import Optional

import numpy as np
from PIL import Image, ImageOps
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor

rdDepictor.SetPreferCoordGen(True)

# pylint: disable=unsubscriptable-object


def trim_image_whitespace(image: Image.Image, padding: int = 0) -> Image.Image:
    """Crops and image to a minimal rectangle.

    This function takes a PIL image and crops it to the minimum rectangle based
    on its whiteness/transparency.

    Args:
        image: PIL image.
        padding: Integer number of pixels to use for padding.

    Returns:
        A new PIL image.
    """
    # Convert to array
    as_array = np.array(image)  # N x N x (r,g,b,a)

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

    image = Image.fromarray(as_array_cropped, mode='RGB')
    return ImageOps.expand(image, border=padding, fill=255)


def mol_to_svg(  # pylint: disable=inconsistent-return-statements
        mol: Chem.Mol,
        min_size: int = 50,
        max_size: int = 300,
        bond_length: int = 25,
        padding: int = 10) -> Optional[str]:
    """Creates a (cropped) SVG molecule drawing as a string.

    Args:
        mol: RDKit Mol.
        min_size: Integer minimum image size (in pixels).
        max_size: Integer maximum image size (in pixels).
        bond_length: Integer bond length (in pixels).
        padding: Integer number of padding pixels in each dimension.

    Returns:
        String SVG image.
    """
    Chem.Kekulize(mol)
    rdDepictor.Compute2DCoords(mol)
    drawer = _draw_svg(mol,
                       size_x=max_size,
                       size_y=max_size,
                       bond_length=bond_length)
    # Find the extent of the drawn image so we can crop the canvas.
    min_x, max_x, min_y, max_y = np.inf, -np.inf, np.inf, -np.inf
    for atom in mol.GetAtoms():
        canvas_x, canvas_y = drawer.GetDrawCoords(atom.GetIdx())
        min_x = min(canvas_x, min_x)
        max_x = max(canvas_x, max_x)
        min_y = min(canvas_y, min_y)
        max_y = max(canvas_y, max_y)
    drawer = _draw_svg(mol,
                       size_x=max(min_size, int(max_x - min_x + 2 * padding)),
                       size_y=max(min_size, int(max_y - min_y + 2 * padding)),
                       bond_length=bond_length)
    match = re.search(r'(<svg\s+.*</svg>)',
                      drawer.GetDrawingText(),
                      flags=re.DOTALL)
    if match:
        return match.group(1)


def _draw_svg(mol: Chem.Mol, size_x: int, size_y: int,
              bond_length: int) -> Draw.MolDraw2DSVG:
    """Creates a canvas and draws a SVG.

    Args:
        mol: RDKit Mol.
        size_x: Integer image size in the x-dimension (in pixels).
        size_y: Integer image size in the y-dimension (in pixels).
        bond_length: Integer bond length (in pixels).

    Returns:
        MolDraw2DSVG.
    """
    drawer = Draw.MolDraw2DSVG(size_x, size_y)
    drawer.drawOptions().fixedBondLength = bond_length
    drawer.drawOptions().padding = 0.0
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer


def mol_to_png(mol: Chem.Mol,
               max_size: int = 1000,
               bond_length: int = 25,
               png_quality: int = 70) -> str:
    """Creates a (cropped) PNG molecule drawing as a string.

    Args:
        mol: RDKit Mol.
        max_size: Integer maximum image size (in pixels).
        bond_length: Integer bond length (in pixels).
        png_quality: Integer PNG quality.

    Returns:
        String PNG image.
    """
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
    return b64.decode('UTF-8')
