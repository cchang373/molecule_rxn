# molecule_rxn
Molecule_rxn is a python package used for generating possible intermediates from molecules in a recursive and brute-force way. SMILES notation can be generated from ASE trajectory file and graph representation.

## Installation
- Create conda environment and install `RDkit`: `conda create -n molecule_rxn -c rdkit/label/beta rdkit`
- Activate conda environment: `conda activate molecule_rxn`
- Install dependencies:
  - imolecule: `pip install imolecule`
  - networkx: `pip install networkx`
  - ASE: `pip install --upgrade --user ase`
- Install `molecule_rxn`: `pip install git+https://github.com/cchang373/molecule_rxn.git`

## Usage
See `example.py`
