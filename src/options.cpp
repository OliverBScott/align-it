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

#include <options.h>
#include <sstream>

Options::Options():
        cutOff(0.0),
        best(0),
        rankby(TANIMOTO),
        funcGroupVec(),
        noHybrid(false),
        epsilon(0.5),
        withExclusion(false),
        merge(false),
        noNormal(false),
        scoreOnly(false),
        isQuiet(false),
        version(false),
        help(false)
{
    refInpStream = NULL;
    refInpFile.clear();
    refInpType.clear();

    dbInpStream = NULL;
    dbInpFile.clear();
    dbInpType.clear();

    molOutStream = NULL;
    molOutWriter = NULL;
    molOutFile.clear();
    molOutType.clear();

    pharmOutStream = NULL;
    pharmOutWriter = NULL;
    pharmOutFile.clear();

    scoreOutStream = NULL;
    scoreOutFile.clear();
}

Options::~Options(void) {
    // Reference input
    if (!refInpFile.empty()) {
        refInpFile.clear();
    }
    if (refInpStream) {
        delete refInpStream;
        refInpStream = NULL;
    }
    // Database input
    if (!dbInpFile.empty()) {
        dbInpFile.clear();
    }
    if (dbInpStream) {
        delete dbInpStream;
        dbInpStream = NULL;
    }
    // Molecule output
    if (!molOutFile.empty()) {
        molOutFile.clear();
    }
    if (molOutStream) {
        delete molOutStream;
        molOutStream = NULL;
    }
    if (molOutWriter) {
        delete molOutWriter;
        molOutWriter = NULL;
    }
    // Pharmacophore output
    if (!pharmOutFile.empty())  {
        pharmOutFile.clear();
    }
    if (pharmOutStream) {
        delete pharmOutStream;
        pharmOutStream = NULL;
    }
    if (pharmOutWriter) {
        delete pharmOutWriter;
        pharmOutWriter = NULL;
    }
    // Score output
    if (!scoreOutFile.empty()) {
        scoreOutFile.clear();
    }
    if (scoreOutStream) {
        delete scoreOutStream;
        scoreOutStream = NULL;
    }
}

std::string Options::print(void) const {
    std::ostringstream os;
    os << std::endl;
    os << "COMMAND LINE OPTIONS:" << std::endl;
    os << std::endl;
    os << "  -> Reference file:    " << (refInpFile.empty() ? "no" : refInpFile) << std::endl;
    os << "  -> Reference type:    ";
    if (!refInpType.empty()) {
        os << refInpType << std::endl;
    } else {
        os << "unknown" << std::endl;
    }
    os << "  -> Database file:     " << (dbInpFile.empty() ? "no" : dbInpFile) << std::endl;
    os << "  -> Database type:     ";
    if (!dbInpType.empty()) {
        os << dbInpType << std::endl;
    } else {
        os << "unknown" << std::endl;
    }
    os << "  -> Mol output file:   " << (molOutFile.empty() ? "no" : molOutFile) << std::endl;
    os << "  -> Output type:       ";
    if (!molOutType.empty()) {
        os << molOutType << std::endl;
    } else {
        os << "unknown" << std::endl;
    }
    os << "  -> Pharm output file: " << (pharmOutFile.empty() ? "no" : pharmOutFile) << std::endl;
    os << "  -> Scores file:       " << (scoreOutFile.empty() ? "no" : scoreOutFile) << std::endl;
    os << "  -> Cutoff:            ";
    if (cutOff) {
        os << cutOff << std::endl;
    } else {
        os << "no" << std::endl;
    }
    os << "  -> Best hits:         ";
    if (best) {
        os << best << std::endl;
    } else {
        os << "no" << std::endl;
    }
    os << "  -> Rank by:           " << rankby << std::endl;
    os << "  -> Functional groups: ";
    if (funcGroupVec[AROM]) os << "AROM ";
    if (funcGroupVec[HDON]) os << "HDON ";
    if (funcGroupVec[HACC]) os << "HACC ";
    if (funcGroupVec[LIPO]) os << "LIPO ";
    if (funcGroupVec[NEGC]) os << "NEGC ";
    if (funcGroupVec[POSC]) os << "POSC ";
    if (funcGroupVec[HYBH]) os << "HYBH ";
    if (funcGroupVec[HYBL]) os << "HYBL ";
    os << std::endl;
    os << "  -> Hybrids:           " << (noHybrid ? "no" : "yes") << std::endl;
    os << "  -> Epsilon:           " << epsilon << std::endl;
    os << "  -> Merge pharm:       " << (merge ? "yes" : "no") << std::endl;
    os << "  -> Include normals:   " << (noNormal ? "no" : "yes") << std::endl;
    os << "  -> With exclusion:    " << (withExclusion ? "yes" : "no") << std::endl;
    os << "  -> Scores only:       " << (scoreOnly ? "yes" : "no") << std::endl;
    os << "  -> Quiet mode:        " << (isQuiet ? "yes" : "no") << std::endl;
    os << std::endl;
    std::string r = os.str();
    return r;
}
