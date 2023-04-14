Generating HOSE codes of molecules with Python. Based on Java CDK HoseCodeGeneratorClass. Incorporate stereo-enhanced HOSE codes described in https://pubs.acs.org/doi/full/10.1021/acsomega.9b00488


### Installation
pip install git+https://github.com/Ratsemaat/HOSE_code_generator

### Usage for standard HOSE codes

The methods get_Hose_codes and get_Hose_codes_from_file can be used to geneate a standard HOSE code. get_Hose_codes accepts a Molecule object, get_Hose_codes_from_file a filepath. Both need the index atom_idx of the atom for which the HOSE code should be created. max_radius gives the number of spheres to use (default 5). usestereo, wedgebond, and strict are only used for the stereo-enhanced hose code (see below).

### Usage for stereo-enhanced HOSE codes

The code can generate stereo-enhanced HOSE codes, as described in https://pubs.acs.org/doi/full/10.1021/acsomega.9b00488 The parameter usestereo needs to be set to true for this. There are a number of additional considerations to keep in mind:

* For E/Z configurations, the detection by rdkit is used. Everything detected by rdkit is included and nothing else
* For wedge bonds, the wedge bonds as in the mol file are used. rdkit detected stereochemistry is ignored. This is because we need more than rdkit does. Take https://nmrshiftdb.nmr.uni-koeln.de/molecule/60004201 as an example. Here, atom 4 will not be considered a stereo centre by rdkit, since the stereochemistry does not matter for compound identity. For our purposes, we need to distinguish atoms 18 and 19, they can have different shifts. Since rdkit reading of molecules removes wedge bonds, the parameter wedgebond in get_Hose_codes is used to pass a map of atom ids and wedge bond values. The values are as in the mol file. The method get_Hose_codes_from_file produces this map internally, so we recommend using that.
* Explicit hydrogens are needed. Chiral centres with e.g. only three atoms around a carbon will be ignored. get_Hose_codes_from_file adds hydrogens internally.
* If hydrogens are added, we recommend using the function makeUpDownBonds once. This is because often only one wedge bond is given in cases where the hydrogen is missng, and it is assumed that the opposite wedge bond is on the hydrogen, which is made explicit by makeUpDownBonds. Again, get_Hose_codes_from_file does this internally.
* The parameter strict (default false) determines if carbon centres with one wedge bond only are considered or not. Strictly speaking this is not possible (around a carbon, three bonds can never be in a plane), but is used. With strict=true, they are ignored. Notice that we do not recommend to add mising hs automatically, not using makeUpDownBonds, and using strict=false. The order may be different in that case from what is taken from two wedge bonds.