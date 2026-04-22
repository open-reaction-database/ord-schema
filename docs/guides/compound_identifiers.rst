#####################
Compound Identifiers
#####################

**********************
Overview
**********************

The table below shows the currently available chemical identifiers in the ORD schema. A
CUSTOM option is included for capturing other chemical identifiers which are not yet
directly included in the schema. While multiple identifiers can be included for a single
chemical, users are actively encouraged to include SMILES or InChI identifiers in their 
datasets. These line notations depict the structure of the chemical and allow the ORD data
browsers to resolve images for the chemicals using RDKit.

===================  ===============================================================================================
Identifier           Description
===================  ===============================================================================================
SMILES               Simplified molecular-input line-entry system.
INCHI                IUPAC International Chemical Identifier.
MOLBLOCK             Molblock from a MDL Molfile V3000.
IUPAC_NAME           Chemical name following IUPAC nomenclature recommendations.
NAME                 Any accepted common name, trade name, etc.
CAS_NUMBER           Chemical Abstracts Service Registry Number (with hyphens).
PUBCHEM_CID          PubChem Compound ID number; `link <https://pubchem.ncbi.nlm.nih.gov/>`__.
CHEMSPIDER_ID        ChemSpider ID number; `link <https://www.chemspider.com/>`__.
CXSMILES             ChemAxon extended SMILES; `link <https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html>`__.
INCHI_KEY            IUPAC International Chemical Identifier key
XYZ                  XYZ molecule file
UNIPROT_ID           UniProt ID (for enzymes); `link <https://www.uniprot.org/>`__.
PDB_ID               Protein data bank ID (for enzymes); `link <https://www.wwpdb.org/>`__.
AMINO_ACID_SEQUENCE  Amino acid sequence (for enzymes).
HELM                 Hierarchical Editing Language for Macromolecules (HELM); `link <https://www.pistoiaalliance.org/helm-notation/>`__.
MDL                  MDL number e.g., MFCD00005972 for morpholine, often included on commercial supplier listings.
CUSTOM               Create your own custom identifiers. Include an informative name in the description field.
===================  ===============================================================================================


When preparing an ORD dataset from a paper, an often overlooked task is the collection
of SMILES (or InChI) strings for the chemicals which may have been recorded in the
paper using names, abbreviations and/or images. This guide suggests some common
workflows for collecting these identifiers.

.. IMPORTANT::

  Special care needs to be taken with the identifiers for organometallic compounds. There are some
  challenges associated with generating SMILES or InChI strings for organometallic compounds, and 
  this can cause problems with the way RDKit parses and resolves them. ORD authors should always 
  check that their SMILES or InChI strings are being resolved correctly, and please refer to :ref:`dative-bonding`
  below for a recommended workflow for including dative bonding in SMILES strings.
  
  For compounds which are hard to define using existing identifiers, it can also be helpful to include
  multiple identifiers (e.g. a NAME, CAS_NUMBER and SMILES) to remove some of the ambiguity.


************************
Draw the Molecule
************************

Most chemical drawing packages can export a SMILES or InChI string for chemical structures you have 
selected. If you have already drawn chemical structures for a published article then it can be a simple
task to export SMILE or InChI strings from your ChemDraw (or other chemical drawing tool) files. This 
`ChemSpider Blog <https://blogs.rsc.org/chemspider/2019/11/08/tips-and-tricks-generating-machine-readable-structural-data-from-a-chemdraw-structure/>`__
shows how to do this in Avogadro, ChemDoodle, ChemDraw, ChemSketch and MarvinSketch.

The following online chemical drawing tools can also be used to obtain SMILES or InChI strings
- `OpenBabel <https://www.cheminfo.org/Chemistry/Cheminformatics/FormatConverter/index.html>`__ (generate SMILES or InChI)
- `InChI web demo <https://iupac-inchi.github.io/InChI-Web-Demo/>`__ (generate InChI only)

.. IMPORTANT::

  There are some challenges associated with generating SMILES or InChI strings for organometallic 
  compounds such as transition metal catalysts. Please refer to :ref:`dative-bonding` below for 
  specific guidance on generating SMILES or InChI for organometallic compounds with dative bonding.

**************************
Look Them Up
**************************

Note about quality of online SMILES

Pubchem
Chemical Suppliers







***************************
Converting Identifiers
***************************







******************************
Image Extraction
******************************

Mathpix




.. _dative-bonding:
*******************************************
Dative Bonding in Organometallic Compounds
*******************************************

There are some challenges associated with the generation of canonical SMILES and InChI strings for
organometallic compounds such as transition metal catalysts. While the above methods can be used to 
obtain a SMILES for such compounds, extra care needs to be taken to check that the SMILES is parsed
and resolved correctly by RDKit.

Some problems to watch out for when obtaining SMILES for organometallic compounds include:

- RDKit SMILES parsing fails for ChemDraw-generated SMILES of plainly drawn metal complexes due to 
  the valence rules on N atoms.
- Using dative bonds in ChemDraw fixes this problem but generates charged species to satisfy the 
  valence rules.
- ChemDraw's syntax for ferrocene complexes leads to RDKit parsing failures.
- PubChem "canonical" SMILES for metal complexes treat each substructure as a separate component 
  (e.g., metals, ligands, ions), leading to a salt-like representation of the compounds.
- Many online sources of SMILES give the salt-like representation used by PubChem.

The current best practice for generating SMILES for transition metal catalysts and other oranometallic 
compounds, which are to be recorded in the Open Reaction Database, is as follows:

1. Use ChemDraw (or other chemical drawing tool) to draw the structure with neutral single bonds to
   ligands. Ferrocenes should be explicitly drawn as cyclopentadienyl-anion-sandwiched Fe atoms.
2. Copy the structure as SMILES (see `ChemSpider Blog <https://blogs.rsc.org/chemspider/2019/11/08/tips-and-tricks-generating-machine-readable-structural-data-from-a-chemdraw-structure/>`__)
   for how to do this in Avogadro, ChemDoodle, ChemDraw, ChemSketch and MarvinSketch.
3. Generate dative '<-' bonds programmatically using the `set_dative_bonds <https://docs.open-reaction-database.org/en/latest/ord_schema/ord_schema.html#ord_schema.message_helpers.set_dative_bonds>`__ message helper.

Here's a snippet of code to generate the dative bonds for a single SMILES string in a Jupyter/Colab 
notebook:

.. code-block:: ipython

   In [1]: # Import packages
   In [2]: from ord_schema import message_helpers
   In [3]: from rdkit import Chem

   In [4]: # Enter your SMILES string
   In [5]: smiles = 'O=S(C)O[Pd]([NH2]C(=C1C=C2)C=C2)(C(=CC=C2)C1=C2)[P](C(C=C1)=C(C(=C(C=C2)P(C(=CC=C3)C=C3)C(=CC=C3)C=C3)C(=C2C=C2)C=C2)C(=C1C=C1)C=C1)(C(C=CC1)=CC=1)C(=CC=C1)C=C1'

   In [6]: # Use RDKit to generate a mol for the SMILES string
   In [7]: mol = Chem.MolFromSmiles(smiles, sanitize=False)
   In [8]: # Use the ORD message_helper to add dative bonds to the mol
   In [9]: dative_mol = message_helpers.set_dative_bonds(mol)
   In [10]: # Use RDKit to convert the mol back to a SMILES
   In [11]: canonical_smiles = Chem.MolToSmiles(dative_mol)

   In [12]: # Print the canonical SMILES string
   In [13]: print(canonical_smiles)
   Out [13]: 'CS(=O)O[Pd]1(<-P(C2=CC=CC=C2)(C2=CC=CC=C2)C2=C(C3=C(P(C4=CC=CC=C4)C4=CC=CC=C4)C=CC4=C3C=CC=C4)C3=C(C=CC=C3)C=C2)<-NC2=C(C=CC=C2)C2=CC=CC=C21'

For a more detailed example refer to this Jupyter Notebook [add link] which shows how to use the above code
to canonicalize a single SMILES, or to batch process multiple SMILES strings input as a .csv file.
The example also shows how to use RDKit to visualize the generated SMILES so you can check that the dative
bonds have been generated correctly.


*****************************************
Checking Your SMILES or InChI Strings
*****************************************


'SMARTS Plus <https://smarts.plus/>'


Check individual SMILES by using a test reaction in ORD Reaction Editor. Does it resolve the 
compound correctly?


For batch processing rdkit can be used. Speak to Ben to get some help with this. Example code for
Pfizer dataset.