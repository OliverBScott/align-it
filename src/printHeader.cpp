/*******************************************************************************
printHeader.cpp - Align-it

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

#include <printHeader.h>

#ifndef USE_RDKIT
#include <openbabel/mol.h>
#else
#include <RDGeneral/versions.h>
#endif

void printHeader() {
    std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cerr << "  Align-it v" << ALIGNIT_VERSION << "." << ALIGNIT_RELEASE << "." << ALIGNIT_SUBRELEASE << " | ";
    std::cerr << __DATE__ " " << __TIME__ << std::endl;
    std::cerr << std::endl;
    std::cerr << "  -> GCC:         " << __VERSION__ << std::endl;
#ifndef USE_RDKIT
    std::cerr << "  -> OpenBabel: " << BABEL_VERSION << std::endl;
#else
    std::cerr << "  -> RDKit: " << RDKit::rdkitVersion << std::endl;
#endif
    std::cerr << std::endl;
    std::cerr << "  Copyright 2012 by Silicos-it, a division of Imacosi bvba" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  Align-it is free software: you can redistribute it and/or modify" << std::endl;
    std::cerr << "  it under the terms of the GNU Lesser General Public License as published" << std::endl;
    std::cerr << "  by the Free Software Foundation, either version 3 of the License, or" << std::endl;
    std::cerr << "  (at your option) any later version." << std::endl;
    std::cerr << std::endl;
    std::cerr << "  Align-it is distributed in the hope that it will be useful," << std::endl;
    std::cerr << "  but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
    std::cerr << "  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
    std::cerr << "  GNU Lesser General Public License for more details." << std::endl;
    std::cerr << std::endl;
#ifndef USE_RDKIT
    std::cerr << "  Align-it is linked against OpenBabel version 3." << std::endl;
    std::cerr << "  OpenBabel is free software; you can redistribute it and/or modify" << std::endl;
    std::cerr << "  it under the terms of the GNU General Public License as published by" << std::endl;
    std::cerr << "  the Free Software Foundation version 2 of the License." << std::endl;
#else
    std::cerr << "  Shape-it is linked against the RDKit (https://www.rdkit.org)." << std::endl;
#endif
    std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cerr << std::endl;
}
