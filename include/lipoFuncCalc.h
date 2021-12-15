/*******************************************************************************
lipoFuncCalc.h - Align-it

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

#ifndef __SILICOSIT_ALIGNIT_LIPOFUNCCALC_H__
#define __SILICOSIT_ALIGNIT_LIPOFUNCCALC_H__

// General
#include <set>
#include <list>
#include <algorithm>

// Toolkit
#ifndef USE_RDKIT
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/ring.h>
#include <openbabel/elements.h>
using Molecule = OpenBabel::OBMol;
#else
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/MolOps.h>
using Molecule = RDKit::ROMol;
#endif

// Align-it
#include <pharmacophore.h>
#include <defines.h>

void lipoFuncCalc(Molecule *, Pharmacophore *);
void _lipoLabelAtoms(Molecule *);
void _lipoGroupAtoms(Molecule *, Pharmacophore *);

#ifndef USE_RDKIT
double _lipoCalcAccSurf(OpenBabel::OBAtom *);
std::list<OpenBabel::OBAtom *> _lipoGetNeighbors(OpenBabel::OBAtom *);
void _lipoLabelNeighbors(OpenBabel::OBAtom*, double);
#else
double _lipoCalcAccSurf(RDKit::Atom *, const RDKit::Conformer &);
std::list<RDKit::Atom *> _lipoGetNeighbors(RDKit::Atom *, const RDKit::Conformer &);
void _lipoLabelNeighbors(RDKit::Atom *, double);
#endif

#endif //__SILICOSIT_ALIGNIT_LIPOFUNCCALC_H__
