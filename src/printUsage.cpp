/*******************************************************************************
printUsage.cpp - Align-it

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

#include "printUsage.h"

void printUsage() {
    printHeader();
    std::cerr << "TASK: " << std::endl;
    std::cerr << std::endl;
    std::cerr << "  Align-it is a tool to generate and align pharmacophores." << std::endl;
    std::cerr << std::endl;
    std::cerr << "INPUT OPTIONS: " << std::endl;
    std::cerr << std::endl;
    std::cerr << "  -r  --reference       <file>" << std::endl;
    std::cerr << "      --refType         <type>" << std::endl;
    std::cerr << "  -d  --dbase           <file>" << std::endl;
    std::cerr << "      --dbType          <type>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "OUTPUT OPTIONS: " << std::endl;
    std::cerr << std::endl;
    std::cerr << "  -p  --pharmacophore   <file>" << std::endl;
    std::cerr << "  -o  --out             <file>" << std::endl;
    std::cerr << "      --outType         <type>" << std::endl;
    std::cerr << "  -s  --scores          <file>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "      --cutOff          <double>" << std::endl;
    std::cerr << "      --best            <int>" << std::endl;
    std::cerr << "      --rankBy          <TANIMOTO|TVERSKY_REF|TVERSKY_DB>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "ALIGN-IT OPTIONS: " << std::endl;
    std::cerr << std::endl;
    std::cerr << "  -f  --funcGroup       <AROM|HDON|HACC|LIPO|CHARGE>" << std::endl;
    std::cerr << "  -e  --epsilon         <double>" << std::endl;
    std::cerr << "  -m  --merge" << std::endl;
    std::cerr << "  -n  --noNormal" << std::endl;
    std::cerr << "      --noHybrid" << std::endl;
    std::cerr << "      --withExclusion" << std::endl;
    std::cerr << "      --scoreOnly" << std::endl;
    std::cerr << std::endl;
    std::cerr << "GENERAL OPTIONS: " << std::endl;
    std::cerr << std::endl;
    std::cerr << "  -h  --help" << std::endl;
    std::cerr << "  -v" << std::endl;
    std::cerr << "      --info            <option>" << std::endl;
    std::cerr << "  -q  --quiet" << std::endl;
    std::cerr << std::endl;
}

























