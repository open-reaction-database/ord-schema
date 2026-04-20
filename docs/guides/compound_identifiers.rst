#####################
Compound Identifiers
#####################

**********************
Overview
**********************

The table below shows the currently available chemical identifiers in the ORD schema. A
CUSTOM option is included for capturing other chemical identifiers which are not yet
directly included in ORD. Where possible submitters are actively encouraged to include 
SMILES or InChI identifiers in their datasets. These line notations depict the structure 
of the chemical and allow the ORD data browsers to resolve images for the chemicals 
using rdkit.

===================  ===============================================================================================
Identifier           Description
===================  ===============================================================================================
SMILES               Simplified molecular-input line-entry system.
INCHI                IUPAC International Chemical Identifier.
MOLBLOCK             Molblock from a MDL Molfile V3000.
IUPAC_NAME           Chemical name following IUPAC nomenclature recommendations.
NAME                 Any accepted common name, trade name, etc.
CAS_NUMBER           Chemical Abstracts Service Registry Number (with hyphens).
PUBCHEM_CID          PubChem Compound ID number.
CHEMSPIDER_ID        ChemSpider ID number.
CXSMILES             ChemAxon extended SMILES
INCHI_KEY            IUPAC International Chemical Identifier key
XYZ                  XYZ molecule file
UNIPROT_ID           UniProt ID (for enzymes)
PDB_ID               Protein data bank ID (for enzymes)
AMINO_ACID_SEQUENCE  Amino acid sequence (for enzymes).
HELM                 HELM; https://www.pistoiaalliance.org/helm-notation/.
MDL                  MDL number e.g., MFCD00005972 for morpholine, often included on commercial supplier listings.
CUSTOM               Create your own custom identifiers. Include an informative name in the description field.
===================  ===============================================================================================


When preparing an ORD dataset from a paper an often overlooked task is the collection
of SMILES (or InChI) strings for the chemicals which may have been recorded in the
paper using names, abbreviations and/or images.

.. IMPORTANT::

  Special care needs to be taken with the identifiers for organometallic compounds. Many 
  of the SMILES and InChI strings reported online for organometallics have the metal and the ligand
  disconnected. Always check that your SMILES or InChI strings are being resolved correctly.
  For best practice refer to section XX for a process to encode dative bonding in the SMILES strings.
  
  For compounds which are hard to define using existing identifiers, it can also be helpful to include
  multiple identifiers (e.g. a NAME, CAS_NUMBER and SMILES) to remove some of the ambiguity.


************************
Draw the Molecule
************************

In the  `ORD Reaction Editor <https://app.open-reaction-database.org/>`_  it is possible to define your 
chemical compounds using the Ketcher chemical drawing tool. This stores the chemical structure as a 
MOLBLOCK identifier, and the corresponding SMILES string is automatically generated and added to the 
compound record.

Many chemical drawing packages can export a SMILES or InChI string for selected structures.

Link to instructions:

ChemDraw
ACD Sketch
ChemDoodle


Links to other online chemical drawing tools.
- OpenBabel <https://www.cheminfo.org/Chemistry/Cheminformatics/FormatConverter/index.html>
- InChI web demo <https://iupac-inchi.github.io/InChI-Web-Demo/>


**************************
Look It Up
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





*******************************
Dative Bonding
*******************************

Link to example Jupyter Notebook


*****************************************
Checking Your SMILES or InChI Strings
*****************************************


'SMARTS Plus <https://smarts.plus/>'


Check individual SMILES by using a test reaction in ORD Reaction Editor. Does it resolve the 
compound correctly?


For batch processing rdkit can be used. Speak to Ben to get some help with this. Example code for
Pfizer dataset.