/*******************************************************************************
pharMerger.cpp - Align-it

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

#include <pharMerger.h>

PharMerger::PharMerger() : _deltaSigma(0.7), _threshold(0.075) {}

void PharMerger::merge(Pharmacophore &phar) {
    typedef std::set<int> Group;
    Group::iterator itG;
    std::list<Group> list1;
    std::list<Group> list2;
    std::list<Group>::iterator itL1, itL2;
    Pharmacophore p;

    int n(phar.size());
    for (int i(0); i < n; ++i) {
        if (phar[i].func == EXCL) {
            p.push_back(phar[i]);
            continue;
        }
        Group g;
        g.insert(i);
        for (int j = i+1; j < n ; ++j) {
            if (phar[i].func != phar[j].func) continue;
            if (VolumeOverlap(phar[i], phar[j], false) > _threshold) {
                g.insert(j);
            }
        }
        list1.push_back(g);
    }
    for (itL1 = list1.begin(); itL1 != list1.end(); ++itL1) {
        bool partOf(false);
        for (itL2 = list1.begin(); itL2 != list1.end(); ++itL2) {
            if (itL1 == itL2) continue;
            if (includes(itL2->begin(), itL2->end(), itL1->begin(), itL1->end())) {
                partOf = true;
                break;
            }
        }
        if (!partOf) {
            if (itL1->size() == 1) {
                p.push_back(phar[*(itL1->begin())]);
            } else {
                list2.push_back(*itL1);
            }
        }
    }
    for (itL1 = list2.begin(); itL1 != list2.end(); ++itL1) {
        PharmacophorePoint pp;
        pp.func = phar[*(itL1->begin())].func;
        double sigma = 0.0;
        for (itG = itL1->begin(); itG != itL1->end(); ++itG) {
            (pp.point).x += (phar[*itG].point).x;
            (pp.point).y += (phar[*itG].point).y;
            (pp.point).z += (phar[*itG].point).z;
            sigma += _deltaSigma/phar[*itG].alpha;
        }
        (pp.point).x /= itL1->size();
        (pp.point).y /= itL1->size();
        (pp.point).z /= itL1->size();
        pp.alpha = 1.0/sigma;
        p.push_back(pp);
    }
    phar = p;
}
