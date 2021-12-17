/*******************************************************************************
alignLib.cpp - Align-it

Copyright 2021 by OliverBScott and the Align-it contributors

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

#include <alignLib.h>

namespace alignit {

Pharmacophore calcPharmacophore(
    Molecule &mol,
    bool calcArom,
    bool calcHDon,
    bool calcHAcc,
    bool calcLipo,
    bool calcCharge,
    bool calcHybrid
) {
    Pharmacophore pharm;
    if (calcArom) aromFuncCalc(&mol, &pharm);
    if (calcHDon) hDonFuncCalc(&mol, &pharm);
    if (calcHAcc) hAccFuncCalc(&mol, &pharm);
    if (calcLipo) lipoFuncCalc(&mol, &pharm);
    if (calcCharge) chargeFuncCalc(&mol, &pharm);
    if (calcHybrid) hybridCalc(&mol, &pharm);
    return pharm;
}

void mergePharmacophore(Pharmacophore &p) {
    PharMerger merger;
    merger.merge(p);
}

Result alignPharmacophores(
    Pharmacophore &ref,
    Pharmacophore &db,
    double epsilon,
    bool useNormals,
    bool useExclusion,
    Molecule *dbMol
) {
    // Prepare reference
    unsigned int exclSize = 0;
    const unsigned int refSize = ref.size();
    double refVolume = 0.0;

    for (unsigned int i(0); i < refSize; ++i) {
        if (ref[i].func == EXCL) {
            for (unsigned int j(0); j < ref.size(); ++j) {
                if (ref[j].func != EXCL) {
                    refVolume -= VolumeOverlap(ref[i], ref[j], useNormals);
                }
            }
            exclSize++;
        } else {
            refVolume += VolumeOverlap(ref[i], ref[i], useNormals);
        }
    }

    // Prepare db
    const unsigned int dbSize = db.size();
    double dbVolume = 0.0;

    for (unsigned int i(0); i < dbSize; ++i) {
        if (db[i].func == EXCL) { continue; }
        dbVolume += VolumeOverlap(db[i], db[i], useNormals);
    }

    // Initialize result
    Result res;
    res.refVolume = refVolume;
    res.dbVolume = dbVolume;
    res.overlapVolume = 0.0;
    res.exclVolume = 0.0;
    res.resPharSize = 0;
    if (dbMol) {
        res.resMol = *dbMol;
    } else {
        res.resMol.addConformer(new RDKit::Conformer(0));
    }

    // Alignment
    FunctionMapping funcMap(&ref, &db, epsilon);
    PharmacophoreMap fMap = funcMap.getNextMap();
    PharmacophoreMap bestMap;
    SiMath::Matrix rotMat(3,3,0.0);

    // Default solution
    auto best = SolutionInfo();
    best.rotor[0] = 1.0;
    double bestScore = -1000;
    int mapSize = fMap.size();
    int maxSize = mapSize - 3;

    // Alignment loop
    while (!fMap.empty()) {
        int msize = fMap.size();
        // Add exclusion spheres if requested
        if (useExclusion) {
            for (unsigned int i(0); i < refSize; ++i) {
                if (ref[i].func != EXCL) {
                    continue;
                }
                for (unsigned int j(0); j < dbSize; ++j) {
                    if (db[j].func == EXCL) {
                        continue;
                    }
                    fMap.insert(std::make_pair(&(ref[i]), &(db[j])));
                }
            }
        }
        // Only align if the expected score has any chance of being larger
        // than best score so far
        if ((msize > maxSize) && (((double) msize / (refSize - exclSize + dbSize - msize)) > bestScore)) {
            Alignment align(fMap);
            SolutionInfo r = align.align(useNormals);
            if (best.volume < r.volume) {
                best = r;
                bestScore = best.volume / (refVolume + dbVolume - best.volume);
                bestMap = fMap;
                mapSize = msize;
            }
        } else break; // Level of mapping site to low
        if (bestScore > 0.98) break;
        fMap.clear();
        fMap = funcMap.getNextMap();
    }
    // Update positions
    rotMat = quat2Rotation(best.rotor);
    positionPharmacophore(db, rotMat, best);
    positionMolecule(&res.resMol, rotMat, best);
    res.info = best;

    // Compute overlap volume between exclusion spheres and pharmacophore points
    for (int i(0); i < refSize; ++i) {
        if (ref[i].func != EXCL) continue;
        for (int j(0); j < dbSize; ++j) {
            res.exclVolume += VolumeOverlap(ref[i], db[j], useNormals);
        }
    }
    for (PharmacophoreMap::iterator itP = bestMap.begin(); itP != bestMap.end(); ++itP) {
        if (((itP->first)->func == EXCL) || ((itP->second)->func == EXCL)) continue;
        res.overlapVolume += VolumeOverlap(itP->first, itP->second, useNormals);
        PharmacophorePoint p(itP->second);
        (res.resPhar).push_back(p);
        ++res.resPharSize;
    }

    // Update scores
    res.info.volume = res.overlapVolume - res.exclVolume;
    if (res.info.volume > 0.0) {
        res.tanimoto = res.info.volume / (res.refVolume + res.dbVolume - res.info.volume);
        res.tversky_ref = res.info.volume / res.refVolume;
        res.tversky_db = res.info.volume / res.dbVolume;
    }
    return res;
}


std::tuple<Pharmacophore, Result> alignMols(
    Molecule &refMol,
    Molecule &dbMol,
    bool calcArom,
    bool calcHDon,
    bool calcHAcc,
    bool calcLipo,
    bool calcCharge,
    bool calcHybrid,
    bool merge,
    double epsilon,
    bool useNormals,
    bool useExclusion
) {
    // Prepare pharmacophores
    Pharmacophore ref = calcPharmacophore(
        refMol, calcArom, calcHDon, calcHAcc, calcLipo, calcCharge, calcHybrid);
    Pharmacophore db = calcPharmacophore(
        dbMol, calcArom, calcHDon, calcHAcc, calcLipo, calcCharge, calcHybrid);
    if (merge) {
        mergePharmacophore(ref);
        mergePharmacophore(db);
    }
    // Alignment
    Result res = alignPharmacophores(
        ref, db, epsilon, useNormals, useExclusion, &dbMol);
    return std::make_tuple(ref, res);
}

} // namespace alignit
