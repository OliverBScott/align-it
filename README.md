# Align-it

Code for align-it with openbabel3 and rdkit (+Python wrappers), inspired by [pyshapeit](https://github.com/rdkit/shape-it)

## Install

* Installing from source:
```
git clone https://github.com/OliverBScott/align-it
cd align-it
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=<where you want to install> ..
make 
make install
```

* To use RDKit set BUILD_RDKIT_SUPPORT=ON
* To build Python wrappers set BUILD_PYTHON_SUPPORT=ON (only available with RDKit)
* You may also need to set RDKIT_INCLUDE_DIR and Boost_INCLUDE_DIR


* Install Python wrappers using pip (requires pip 10 or greater) (experimental)
```
git clone https://github.com/OliverBScott/align-it
pip install ./align-it
```

* It is assumed that RDKit is already installed using conda


## Original code and basic usage
- [align-it webpage](http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/align-it/1.0.4/align-it.html) (CLI usage)
- [Python library example](https://github.com/OliverBScott/align-it/tree/master/example)

## How to cite
If you use this code in your research, please include the following citation in your publication:

Taminau, J.; Thijs, G.; De Winter, H. (2008) ‘Pharao: Pharmacophore alignment and optimization’, [*J. Mol. Graph. Model.* 2008, **27**, 161-169](https://doi.org/10.1016/j.jmgm.2008.04.003)