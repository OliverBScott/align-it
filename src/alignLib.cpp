/*******************************************************************************
alignLib.cpp - Align-it

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

#include <alignLib.h>

namespace alignit {

Pharmacophore calcPharmacophore(
    Molecule &mol,
    bool calcArom,
    bool calcHDon,
    bool calcHAcc,
    bool calcLipo,
    bool calcCharge,
    bool calcHybrid
) {
    Pharmacophore pharm;
    if (calcArom) aromFuncCalc(&mol, &pharm);
    if (calcHDon) hDonFuncCalc(&mol, &pharm);
    if (calcHAcc) hAccFuncCalc(&mol, &pharm);
    if (calcLipo) lipoFuncCalc(&mol, &pharm);
    if (calcCharge) chargeFuncCalc(&mol, &pharm);
    if (calcHybrid) hybridCalc(&mol, &pharm);
    return pharm;
}

void mergePharmacophore(Pharmacophore &p) {
    PharMerger merger;
    merger.merge(p);
}




} // namespace alignit
