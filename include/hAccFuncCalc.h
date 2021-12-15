/*******************************************************************************
hAccFuncCalc.h - Align-it

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

#ifndef __SILICOSIT_ALIGNIT_HACCFUNCCALC_H__
#define __SILICOSIT_ALIGNIT_HACCFUNCCALC_H__

// General
#include <vector>
#include <list>

// Toolkit
#ifndef USE_RDKIT
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#else
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Conformer.h>
#endif

// Align-it
#include <pharmacophore.h>
#include <defines.h>

#ifndef USE_RDKIT
void hAccFuncCalc(OpenBabel::OBMol *, Pharmacophore *);
double _hAccCalcAccSurf(OpenBabel::OBAtom *);
std::list<OpenBabel::OBAtom *> _hAccGetNeighbors(OpenBabel::OBAtom *);
Coordinate _hAccCalcNormal(OpenBabel::OBAtom *);
bool _hAccDelocalized(OpenBabel::OBAtom *);
#else
void hAccFuncCalc(RDKit::ROMol *, Pharmacophore *);
double _hAccCalcAccSurf(RDKit::Atom *, const RDKit::Conformer &);
std::list<RDKit::Atom *> _hAccGetNeighbors(RDKit::Atom*, const RDKit::Conformer &);
Coordinate _hAccCalcNormal(RDKit::Atom *, const RDKit::Conformer &);
bool _hAccDelocalized(RDKit::Atom *);
#endif

#endif //__SILICOSIT_ALIGNIT_HACCFUNCCALC_H__
