/*******************************************************************************
functionMapping.cpp - Align-it

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

#include <functionMapping.h>

FunctionMapping::FunctionMapping(Pharmacophore* ref, Pharmacophore *db, double eps) {
    // Initial values
    _hasNext = true;
    _maxLevel = 0;
    _epsilon = eps;
    _ref = ref;
    _db = db;
    _refIndex.clear();
    _dbIndex.clear();

    // Create an initial set of corresponding functional groups from reference and database
    bool initCounts(false);
    for (unsigned int i(0); i < _ref->size(); ++i) {
        if ((*_ref)[i].func == EXCL) { continue; }
        for (unsigned int j(0); j < _db->size(); ++j) {
            if ((*_ref)[i].func == (*_db)[j].func) {
                _refIndex.push_back(i);
                _dbIndex.push_back(j);
            } else if (((*_ref)[i].func == HYBH) && ((*_db)[j].func == HDON)) {
                _refIndex.push_back(i);
                _dbIndex.push_back(j);
            } else if (((*_ref)[i].func == HYBH) && ((*_db)[j].func == HACC)) {
                _refIndex.push_back(i);
                _dbIndex.push_back(j);
            } else if (((*_ref)[i].func == HDON) && ((*_db)[j].func == HYBH)) {
                _refIndex.push_back(i);
                _dbIndex.push_back(j);
            }
            else if (((*_ref)[i].func == HACC) && ((*_db)[j].func == HYBH)) {
                _refIndex.push_back(i);
                _dbIndex.push_back(j);
            } else if (((*_ref)[i].func == HYBL) && ((*_db)[j].func == AROM)) {
                _refIndex.push_back(i);
                _dbIndex.push_back(j);
            } else if (((*_ref)[i].func == HYBL) && ((*_db)[j].func == LIPO)) {
                _refIndex.push_back(i);
                _dbIndex.push_back(j);
            } else if (((*_ref)[i].func == AROM) && ((*_db)[j].func == HYBL)) {
                _refIndex.push_back(i);
                _dbIndex.push_back(j);
            } else if (((*_ref)[i].func == LIPO) && ((*_db)[j].func == HYBL)) {
                _refIndex.push_back(i);
                _dbIndex.push_back(j);
            }
        }
        // If we get here in the loop, the number of properties should be counted
        initCounts = true;
    }
    // check if the map is empty
    if (_refIndex.empty() && _dbIndex.empty()) {
        return;
    }
    // number of possible combinations
    unsigned int n = _refIndex.size();
    _matchMap[0] = new std::vector<std::vector<unsigned int>>;
    _matchMap[1] = new std::vector<std::vector<unsigned int>>;
    for (unsigned int i(0); i < n; ++i) {
        std::vector<unsigned int> v(1);	// create a new vector to store the first pairs
        v[0] = i;
        _matchMap[1]->push_back(v);
    }
    // _maxLevel = (_maxLevel < ref.size()) ? _maxLevel : ref.size();
    // _maxLevel = (_maxLevel < db.size()) ? _maxLevel : db.size();
    _maxLevel = (_ref->size() < _db->size()) ? _ref->size() : _db->size();

    // create a new vector to store the first pairs
    _matchMap[2] = new std::vector<std::vector<unsigned int> >;

    // now create possible pairs
    for (unsigned int i(0); i < n-1; ++i) {
        double v1 = GCI * pow(PI/(*_ref)[_refIndex[i]].alpha, 1.5);
        double v2 = GCI * pow(PI/(*_db)[_dbIndex[i]].alpha, 1.5);
        for (unsigned int j(i+1); j < n; ++j) {
            if ((_refIndex[i] == _refIndex[j]) || (_dbIndex[i] == _dbIndex[j])) {
                // same group in one of the pairs, not possible
                continue;
            }
            // compute distances between functional groups within a molecule
            double d1 = distance((*_ref)[_refIndex[i]].point, (*_ref)[_refIndex[j]].point);
            double d2 = distance((*_db)[_dbIndex[i]].point, (*_db)[_dbIndex[j]].point);

            // check if the relative overlap between points is large enough
            d1 = (d1 - d2) * (d1 - d2);
            double o1 = GCI2 * pow(PI/((*_ref)[_refIndex[i]].alpha + (*_db)[_dbIndex[i]].alpha), 1.5) * exp(-((*_ref)[_refIndex[i]].alpha * (*_db)[_dbIndex[i]].alpha) * d1/((*_ref)[_refIndex[i]].alpha + (*_db)[_dbIndex[i]].alpha));
            double o2 = GCI2 * pow(PI/((*_ref)[_refIndex[j]].alpha + (*_db)[_dbIndex[j]].alpha), 1.5) * exp(-((*_ref)[_refIndex[j]].alpha * (*_db)[_dbIndex[j]].alpha) * d1/((*_ref)[_refIndex[j]].alpha + (*_db)[_dbIndex[j]].alpha));
            double v3 = GCI * pow(PI/(*_ref)[_refIndex[j]].alpha, 1.5);
            double v4 = GCI * pow(PI/(*_db)[_dbIndex[j]].alpha, 1.5);
            if ((o2 / (v3 + v4 - o2) > _epsilon) || (o1 / (v1 + v2 - o1) > _epsilon)) {
                std::vector<unsigned int> v(2);
                v[0] = i;
                v[1] = j;
                _matchMap[2]->push_back(v);
            }
        }
    }
    for (unsigned int level(3); level <= _maxLevel; ++level) {
        _matchMap[level] = new std::vector<std::vector<unsigned int>>;
        // loop over all stored combinations in the previous level
        std::vector<std::vector<unsigned int>>::iterator vi, vj, vk;
        for (vi = _matchMap[level-1]->begin(); vi != _matchMap[level-1]->end(); ++vi) {
            // Get index of the last pair
            unsigned int lastPairIndex = (*vi)[level-2];
            // Find in all pairs
            for (vj = _matchMap[2]->begin(); vj != _matchMap[2]->end(); ++vj) {
                if ((*vj)[0] < lastPairIndex) { continue; }
                if ((*vj)[0] > lastPairIndex) { break;}
                // try to add it to the combination
                bool possible = true;
                for (unsigned int i(0); i < level-2; ++i) {
                    bool found = false;
                    for (vk = _matchMap[2]->begin(); vk != _matchMap[2]->end(); ++vk) {
                        if ((*vk)[0] < (*vi)[i]) { continue; }
                        if ((*vk)[0] > (*vi)[i]) { break; }
                        if ((*vk)[1] == (*vj)[1]) { // found corresponding pair
                            found = true;
                            break;
                        }
                    }
                    possible &= found;
                    if (!possible) { break; }
                }
                // adding the pair if it is possible
                if (possible) {
                    std::vector<unsigned int> v(level);
                    for (unsigned int i(0); i < level-1; ++i) {
                        v[i] = (*vi)[i];
                    }
                    v[level-1] = (*vj)[1];
                    _matchMap[level]->push_back(v);
                }
            }
        } // end of loop over combinations
    }
}

FunctionMapping::~FunctionMapping() {
    std::map<unsigned int, std::vector<std::vector<unsigned int> > *>::iterator mi;
    for (mi = _matchMap.begin(); mi != _matchMap.end(); ++mi) {
        if (mi->second != NULL) {
            delete mi->second;
        }
    }
}

PharmacophoreMap FunctionMapping::getNextMap() {
    PharmacophoreMap pm;
    while ((_maxLevel != 0) && _matchMap[_maxLevel]->empty()) {
        --_maxLevel;
    }
    if (_maxLevel < 1) {
        pm.clear();
        return pm;
    }
    std::vector<std::vector<unsigned int>>::iterator vi = _matchMap[_maxLevel]->begin();
    // loop over all combinations in the current selected vector
    // and add them to the function map
    for (unsigned int i(0); i < _maxLevel; ++i) {
        unsigned int i1 = _refIndex[(*vi)[i]];
        unsigned int i2 = _dbIndex[(*vi)[i]];
        pm.insert(std::make_pair(&((*_ref)[i1]),&((*_db)[i2])));
    }
    // remove the selected combination
    _matchMap[_maxLevel]->erase(vi);
    return pm;
}
