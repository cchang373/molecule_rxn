# molecule_rxn
Molecule_rxn is a python package used for generating possible intermediates from molecules in a recursive and brute-force way. SMILES notation can be generated from ASE trajectory file and graph representation. EC_fp enables cutomized variants for atoms in Morgan Fingerprinting.

## Installation
- Create conda environment and install `RDkit`: `conda create -n molecule_rxn -c rdkit/label/beta rdkit python=3.6`
- Activate conda environment: `conda activate molecule_rxn`
- Install dependencies:
  - imolecule: `pip install imolecule`
  - networkx: `pip install networkx`
  - ASE: `pip install --upgrade --user ase`
  - Open Babel: `conda install -c conda-forge openbabel`
- Install `molecule_rxn`: `pip install git+https://github.com/cchang373/molecule_rxn.git`

## Usage
### Generation of intermediates and SMILES notation, graph representation
See `example.py`
### Generation of surface species
```
from molecule_rxn.add_metal import add_metal, add_metal_all
species_single = add_metal().add_metal('[O]C[C]', 'Rh')
species_all = add_metal_all().add_metal_all(['[O]C[C]'], 'Rh')
```
### Generation of modified Morgan Fingerprinting
```
from molecule_rxn.EC_fp.EC_fp import ecfp

bitInfo = ecfp('CO', 2) #SMILES notation, radius
```
