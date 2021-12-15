/*******************************************************************************
pharmacophore.h - Align-it

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

#ifndef __SILICOSIT_ALIGNIT_PHARMACOPHORE_H__
#define __SILICOSIT_ALIGNIT_PHARMACOPHORE_H__

// General
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <list>
#include <map>

// Align-it
#include <coordinate.h>
#include <siMath.h>

enum FuncGroup {
    AROM,   ///< Aromatic ring-system, calculated by the AromFuncCalc class
    HDON,  ///< Hydrogen donor, calculated by the HDonFuncCalc class
    HACC,  ///< Hydrogen acceptor, calculated by the HAccFuncCalc class
    LIPO,  ///< Lipophilicity, calculated by the LipoFuncCalc class
    POSC,  ///< Positive charge, calculated by the ChargeFuncCalc class
    NEGC,  ///< Negative charge, calculated by the ChargeFuncCalc class
    HYBH,  ///< Hybrid Type: HDON + HACC
    HYBL,  ///< Hybrid Type: AROM + LIPO
    EXCL,  ///< Exclusion sphere
    UNDEF, ///< Undefined value (typically used for initialisation)
};

const std::string funcName[10] = {
    "SILICOS::PHARAO::AROM",
    "SILICOS::PHARAO::HDON",
    "SILICOS::PHARAO::HACC",
    "SILICOS::PHARAO::LIPO",
    "SILICOS::PHARAO::POSC",
    "SILICOS::PHARAO::NEGC",
    "SILICOS::PHARAO::HYBH",
    "SILICOS::PHARAO::HYBL",
    "SILICOS::PHARAO::EXCL",
    "SILICOS::PHARAO::UNDEF"
};

const bool funcHasNormal[10] = {
    true,
    true,
    true,
    false,
    false,
    false,
    true,
    true,
    false,
    false,
};

const double funcSigma[10] = {
    0.7, // AROM
    1.0, // HDON
    1.0, // HACC
    0.7, // LIPO
    1.0, // POSC
    1.0, // NEGC
    1.0, // HYBH
    0.7, // HYBL
    1.6, // EXCL
    1.0, // UNDEF
};

class PharmacophorePoint {
public:
    Coordinate point;
    Coordinate normal;
    FuncGroup func;
    double alpha;
    bool hasNormal;

    PharmacophorePoint(void);
    PharmacophorePoint(const PharmacophorePoint &);
    PharmacophorePoint(const PharmacophorePoint *);

    // Dummy operators (Do not use)
    bool operator==(const PharmacophorePoint& other) {return false;}
    bool operator!=(const PharmacophorePoint& other) {return true;}
};

typedef std::vector<PharmacophorePoint> Pharmacophore;
typedef std::multimap<PharmacophorePoint*, PharmacophorePoint*> PharmacophoreMap;

class PharmacophoreReader {
private:
    void _skipPharmacophore(std::ifstream *);
public:
    PharmacophoreReader(void);
    ~PharmacophoreReader(void);

    Pharmacophore read(std::ifstream *, std::string &);
};

class PharmacophoreWriter {
public:
    PharmacophoreWriter(void);
    ~PharmacophoreWriter(void);

    void write(Pharmacophore &, std::ofstream*, const std::string &);
};

#endif //__SILICOSIT_ALIGNIT_PHARMACOPHORE_H__
