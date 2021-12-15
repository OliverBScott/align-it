/*******************************************************************************
functionMapping.h - Align-it

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

#ifndef __SILICOSIT_ALIGNIT_FUNCTIONMAPPING_H__
#define __SILICOSIT_ALIGNIT_FUNCTIONMAPPING_H__

// General
#include <vector>
#include <map>

// Align-it
#include <pharmacophore.h>
#include <utilities.h>

class FunctionMapping {
public:
    FunctionMapping(Pharmacophore *, Pharmacophore *, double);
    ~FunctionMapping(void);

    PharmacophoreMap getNextMap(void);

private:
    bool _hasNext;
    unsigned int _maxLevel;
    double _epsilon;

    Pharmacophore *_ref;  // Pointer to reference pharmacophore
    Pharmacophore *_db;   // Pointer to database pharmacophore

    std::vector<unsigned int> _refIndex;
    std::vector<unsigned int> _dbIndex;
    std::map<unsigned int, std::vector<std::vector<unsigned int>>*> _matchMap;
};

#endif //__SILICOSIT_ALIGNIT_FUNCTIONMAPPING_H__
