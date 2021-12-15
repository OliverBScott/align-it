/*******************************************************************************
options.h - Align-it

Copyright 2012-2013 by Silicos-it, a division of Imacosi BVBA

This file is part of Align-it.

	Align-it is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published
	by the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Align-it is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with Align-it.  If not, see <http://www.gnu.org/licenses/>.

Align-it can be linked against OpenBabel version 3 or the RDKit.

	OpenBabel is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation version 2 of the License.

***********************************************************************/

#ifndef __SILICOSIT_ALIGNIT_OPTIONS_H__
#define __SILICOSIT_ALIGNIT_OPTIONS_H__

// General
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

// Toolkit
#ifndef USE_RDKIT
#include <openbabel/obconversion.h>
#else
#include <GraphMol/FileParsers/MolWriters.h>
#endif

// Align-it
#include <fileType.h>
#include <rankType.h>
#include <pharmacophore.h>

class Options {
public:
    std::string refInpFile;
    std::string refInpType;
    std::ifstream *refInpStream;

    std::string dbInpFile;
    std::string dbInpType;
    std::ifstream *dbInpStream;

    std::string pharmOutFile;
    std::ofstream *pharmOutStream;
    PharmacophoreWriter *pharmOutWriter;

    std::string molOutFile;
    std::string molOutType;
    std::ofstream *molOutStream;
#ifndef USE_RDKIT
    OpenBabel::OBConversion *molOutWriter;
#else
    RDKit::SDWriter *molOutWriter;
#endif

    std::string	scoreOutFile;
    std::ofstream *scoreOutStream;

    double cutOff;
    int best;
    RankType rankby;

    std::vector<bool> funcGroupVec;
    bool noHybrid;
    double epsilon;
    bool withExclusion;
    bool merge;
    bool noNormal;
    bool scoreOnly;

    bool isQuiet;
    bool version;
    bool help;

    Options(void);
    ~Options(void);

    std::string print(void) const;
};

#endif //__SILICOSIT_ALIGNIT_OPTIONS_H__
