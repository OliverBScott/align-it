/*******************************************************************************
alignLib.h - Align-it

Copyright 2021 by OliverBScott and the Align-it contributors

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

#ifndef __SILICOSIT_ALIGNIT_ALIGNLIB_H__
#define __SILICOSIT_ALIGNIT_ALIGNLIB_H__

// Align-it
#include <pharmacophore.h>
#include <pharMerger.h>
#include <calcPharm.h>
#include <functionMapping.h>
#include <solutionInfo.h>
#include <alignment.h>
#include <siMath.h>
#include <result.h>

// Toolkit
#ifndef USE_RDKIT
#include <openbabel/mol.h>
using Molecule = OpenBabel::OBMol;
#else
#include <GraphMol/ROMol.h>
using Molecule = RDKit::ROMol;
#endif


namespace alignit {

// Read + write utilities for pharmacophores

Pharmacophore calcPharmacophore(
    Molecule &mol,
    bool calcArom = true,
    bool calcHDon = true,
    bool calcHAcc = true,
    bool calcLipo = true,
    bool calcCharge = true,
    bool calcHybrid = true
);

void mergePharmacophore(Pharmacophore &p);

// not sure how to align in the best way
// pharmacophores + mols

} // namespace alignit

#endif //__SILICOSIT_ALIGNIT_ALIGNLIB_H__
