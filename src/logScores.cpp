/*******************************************************************************
logScores.cpp - Align-it

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

#include <logScores.h>

void logScores(Result *res, Options &uo) {
    *uo.scoreOutStream  << res->refId << "\t"
        << std::setprecision(6) << res->refVolume << "\t"
        << res->dbId << "\t"
        << std::setprecision(6) <<  res->dbVolume  << "\t"
        << std::setprecision(6) <<  res->overlapVolume << "\t"
        << std::setprecision(6) <<  res->exclVolume << "\t"
        << std::setprecision(6) <<  (res->info).volume << "\t"
        << res->resPharSize <<  "\t"
        << std::setprecision(4) <<  res->tanimoto << "\t"
        << std::setprecision(4) <<  res->tversky_ref << "\t"
        << std::setprecision(4) <<  res->tversky_db << std::endl;
}