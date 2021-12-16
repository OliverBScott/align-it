/*******************************************************************************
pharmacophore.cpp - Align-it

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

#include <pharmacophore.h>
#include <stringTokenizer.h>

PharmacophorePoint::PharmacophorePoint() {
    func = UNDEF;
    alpha = 1.0;
    normal.x = 0.0;
    normal.y = 0.0;
    normal.z = 0.0;
    point.x = 0.0;
    point.y = 0.0;
    point.z = 0.0;
    hasNormal = false;
}

PharmacophorePoint::PharmacophorePoint(const PharmacophorePoint &p) {
    point = p.point;
    func = p.func;
    alpha = p.alpha;
    normal = p.normal;
    hasNormal = p.hasNormal;
}

PharmacophorePoint::PharmacophorePoint(const PharmacophorePoint * p) {
    point = p->point;
    func = p->func;
    alpha = p->alpha;
    normal = p->normal;
    hasNormal = p->hasNormal;
}

PharmacophoreReader::PharmacophoreReader(void) {}
PharmacophoreReader::~PharmacophoreReader(void) {}

Pharmacophore PharmacophoreReader::read(std::ifstream *_input, std::string &name) {
    Pharmacophore pharmacophore;
    pharmacophore.clear();
    if (!*_input) {
        mainWar("Unable to read pharmacophore. Input stream not accessible.");
        return pharmacophore;
    }
    std::string line;
    getline(*_input, line);
    name = line;
    getline(*_input, line);

    while (line != "$$$$") {
        if (line[0]=='#') continue;
        PharmacophorePoint p;
        auto lineList = stringTokenizer(line, "\t");
        std::vector<std::string> lineVec;
        lineVec.clear();
        std::list<std::string>::iterator li;
        for (li = lineList.begin(); li != lineList.end(); ++li) {
            lineVec.push_back(*li);
        }
        if (lineVec.size() < 8) {
            _skipPharmacophore(_input);
            if (!_input->eof()) {
                mainWar("Incorrect line: " + line);
                pharmacophore.clear();
                return pharmacophore;
            } else {
                return pharmacophore;
            }
        }
        bool isOk(false);
        if (lineVec[0] == "AROM") { isOk = true; p.func = AROM; }
        if (lineVec[0] == "HDON") { isOk = true; p.func = HDON; }
        if (lineVec[0] == "HACC") { isOk = true; p.func = HACC; }
        if (lineVec[0] == "LIPO") { isOk = true; p.func = LIPO; }
        if (lineVec[0] == "POSC") { isOk = true; p.func = POSC; }
        if (lineVec[0] == "NEGC") { isOk = true; p.func = NEGC; }
        if (lineVec[0] == "HYBH") { isOk = true; p.func = HYBH; }
        if (lineVec[0] == "HYBL") { isOk = true; p.func = HYBL; }
        if (lineVec[0] == "EXCL") { isOk = true; p.func = EXCL; }
        if (!isOk) {
            _skipPharmacophore(_input);
            mainWar("Incorrect functional group:: " + line);
            pharmacophore.clear();
            return pharmacophore;
        }
        p.point.x = strtod(lineVec[1].c_str(), NULL);
        p.point.y = strtod(lineVec[2].c_str(), NULL);
        p.point.z = strtod(lineVec[3].c_str(), NULL);
        p.alpha = strtod(lineVec[4].c_str(), NULL);

        if (lineVec[5] == "1") { p.hasNormal = true; }
        if (lineVec[5] == "0") { p.hasNormal = false; }

        p.normal.x = strtod(lineVec[6].c_str(), NULL);
        p.normal.y = strtod(lineVec[7].c_str(), NULL);
        p.normal.z = strtod(lineVec[8].c_str(), NULL);

        pharmacophore.push_back(p);
        getline(*_input, line);
    }
    return pharmacophore;
}

void PharmacophoreReader::_skipPharmacophore(std::ifstream *_input) {
    std::string line;
    getline(*_input, line);
    while (line != "$$$$") {
        if (_input->eof()) return;
        getline(*_input, line);
    }
}

PharmacophoreWriter::PharmacophoreWriter(void) {}
PharmacophoreWriter::~PharmacophoreWriter(void) {}

void PharmacophoreWriter::write(Pharmacophore &p, std::ofstream *os, const std::string &name) {
    *os << name << std::endl;
    Pharmacophore::iterator itP;
    for (itP = p.begin(); itP != p.end(); ++itP) {
        switch (itP->func) {
            case AROM: *os << "AROM\t"; break;
            case HDON: *os << "HDON\t"; break;
            case HACC: *os << "HACC\t"; break;
            case LIPO: *os << "LIPO\t"; break;
            case POSC: *os << "POSC\t"; break;
            case NEGC: *os << "NEGC\t"; break;
            case HYBH: *os << "HYBH\t"; break;
            case HYBL: *os << "HYBL\t"; break;
            case EXCL: *os << "EXCL\t"; break;
            case UNDEF: break;
        }
        Coordinate c(itP->point);
        *os << c.x << '\t' << c.y << '\t' << c.z << '\t';
        *os << itP->alpha << '\t';
        itP->hasNormal ? (*os << "1\t") : (*os << "0\t");

        Coordinate n(itP->normal);
        *os << n.x << '\t' << n.y << '\t' << n.z << std::endl;
    }
    *os << "$$$$" << std::endl;
}
