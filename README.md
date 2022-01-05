# Align-it

Code for align-it with openbabel3 and rdkit (+Python wrappers), inspired by [pyshapeit](https://github.com/rdkit/shape-it)

Code adapted from [original](http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/align-it/1.0.4/align-it.html) and [align-it-ob3](https://github.com/iwatobipen/align-it-ob3).

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

## Original code and basic usage
- [align-it webpage](http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/align-it/1.0.4/align-it.html) (CLI usage and original code)
- [Python library example](https://github.com/OliverBScott/align-it/tree/master/example)

## Disclaimer
This repo represents my attempt at updating align-it and creating Python wrappers, in the hope that it may be useful. It is 
however likely that bugs and differences with the original code will exist. Any contributions are welcome!   

## How to cite
If you use this code in your research, please include the following citation in your publication:

Taminau, J.; Thijs, G.; De Winter, H. (2008) ‘Pharao: Pharmacophore alignment and optimization’, [*J. Mol. Graph. Model.* 2008, **27**, 161-169](https://doi.org/10.1016/j.jmgm.2008.04.003)