/*******************************************************************************
main.cpp - Align-it

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

// General
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

// Toolkit
#ifndef USE_RDKIT
// OpenBabel
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#else
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#endif


// Align-it
#include <options.h>
#include <logOut.h>
#include <logScores.h>
#include <logPharmacophores.h>
#include <addBest.h>
#include <parseCommandLine.h>
#include <pharmacophore.h>
#include <pharMerger.h>
#include <printHeader.h>
#include <printUsage.h>
#include <printProgress.h>
#include <getExt.h>
#include <calcPharm.h>
#include <result.h>
#include <alignment.h>
#include <functionMapping.h>

//*--------------------------------------------------------------------------*//
PharMerger pharMerger;

//*--------------------------------------------------------------------------*//
//* MAIN                                                                MAIN *//
//*--------------------------------------------------------------------------*//
int main(int argc, char* argv[]){

    // Initialise random number generator
    srandom(time(NULL));
    clock_t t0 = clock();

    // Read options
    Options uo = parseCommandLine(argc, argv);
    if (!uo.noHybrid) {
        if (uo.funcGroupVec[AROM] && uo.funcGroupVec[LIPO]) {
            uo.funcGroupVec[HYBL] = true;
        }
        if (uo.funcGroupVec[HDON] && uo.funcGroupVec[HACC]) {
            uo.funcGroupVec[HYBH] = true;
        }
    }
    if (uo.version){
        printHeader();
        exit(0);
    }
    if (uo.help){
        printUsage();
        exit(0);
    }

    // Print header
    printHeader();
    std::cerr << uo.print() << std::endl;

    // Db file and pharmacophore out are mandatory elements
    if (uo.dbInpFile.empty()) {
        mainErr("Missing database file. This is a required option (-d | --dbase).");
    }
    if (uo.pharmOutFile.empty() && uo.molOutFile.empty() && uo.scoreOutFile.empty()) {
        mainErr("No output file defined. So there is actually no use to compute anything at all.");
    }
    if ((uo.pharmOutFile.empty() && uo.scoreOutFile.empty()) && !uo.molOutFile.empty()) {
        mainErr("No file defined to write pharmacophore information.");
    }
    if (uo.refInpFile.empty() && uo.pharmOutFile.empty() && uo.molOutFile.empty() && !uo.scoreOutFile.empty()) {
        mainErr("Only score file requested when no reference is given. Unable to generate this output.");
    }

    // Reference variables
    Pharmacophore refPharm;
    refPharm.clear();
    std::string refId;
    double refVolume(0.0);
    int refSize(0);
    int exclSize(0);

    // Database variables
    std::vector<Result*> resList;
    Pharmacophore dbPharm;
    std::string dbId;
    double dbVolume(0.0);
    int dbSize(0);

#ifdef USE_RDKIT
    bool takeOwnership = false;
    bool sanitize = true;
    bool removeHs = false;
#endif

    //----------------------------------------------------------------------------
    //...(A).. Process the reference
    //----------------------------------------------------------------------------

    if (!uo.refInpFile.empty()) {

        //-------------------------------------------------------
        //...(1).. get reference pharmacophore
        //-------------------------------------------------------

        if (!uo.refInpFile.empty() && uo.refInpType != "phar") {
#ifndef USE_RDKIT
            OpenBabel::OBMol m;
            OpenBabel::OBConversion* reader = new OpenBabel::OBConversion();
            reader->SetInFormat(reader->FindFormat(uo.refInpType));
            if (!reader->Read(&m, uo.refInpStream)) {
                mainErr("Unable to read reference molecule");
            }
            calcPharm(&m, &refPharm, uo);
            refId = m.GetTitle();
            delete reader;
            reader = NULL;
#else
            RDKit::ForwardSDMolSupplier reader(
                    uo.refInpStream, takeOwnership, sanitize, removeHs);
            std::unique_ptr<RDKit::ROMol> refmptr(reader.next());
            if (!refmptr) {
                mainErr("Could not parse reference molecule");
            }
            RDKit::ROMol &refMol = *refmptr;
            calcPharm(&refMol, &refPharm, uo);
            refMol.getPropIfPresent("_Name", refId);
#endif
            if (refId == "") {
                refId = "Unnamed_ref";
            }
        } else if (uo.refInpType == "phar") {
            PharmacophoreReader *reader = new PharmacophoreReader();
            refPharm = reader->read(uo.refInpStream, refId);
            if (refPharm.empty()) {
                mainErr("Error reading reference pharmacophore");
            }
            delete reader;
            reader = NULL;
        } else {
            mainErr("Unknown format of reference molecule");
        }

        //-------------------------------------------------------
        //...(2).. process reference pharmacophore
        //-------------------------------------------------------

        if (uo.merge) {
            pharMerger.merge(refPharm);
        }
        refSize = refPharm.size();
        for (unsigned int i(0); i < refSize; ++i) {
            if (refPharm[i].func == EXCL) {
                // extract overlap with exclusion spheres
                for (unsigned int j(0); j < refPharm.size(); ++j) {
                    if (refPharm[j].func != EXCL) {
                        refVolume -= VolumeOverlap(refPharm[i], refPharm[j], !uo.noNormal);
                    }
                }
                exclSize++;
            } else {
                refVolume += VolumeOverlap(refPharm[i], refPharm[i], !uo.noNormal);
            }
        }
        if (!uo.isQuiet) {
            std::cerr << "Reference pharmacophore " << refId << std::endl;
            std::cerr << "   number of points:            " << refSize - exclSize << std::endl;
            std::cerr << "   number of exclusion spheres: " << exclSize << std::endl;
            std::cerr << "   total volume:                " << refVolume << std::endl;
        }
    }

    //----------------------------------------------------------------------------
    //...(B).. Process the database file
    //----------------------------------------------------------------------------

    // local storage of the rotation matrix
    SiMath::Matrix rotMat(3,3,0.0);
    unsigned int molCount(0);
    PharmacophoreReader* pharmReader = NULL;
#ifndef USE_RDKIT
    OpenBabel::OBConversion* molReader = NULL;
#else
    RDKit::ForwardSDMolSupplier* molReader = NULL;
#endif
    if (uo.dbInpType == "phar") {
        pharmReader = new PharmacophoreReader();
    } else if (!uo.dbInpType.empty() && (uo.dbInpType != "phar")) {
#ifndef USE_RDKIT
        molReader = new OpenBabel::OBConversion();
		molReader->SetInFormat(molReader->FindFormat(uo.dbInpType));
      	molReader->SetInStream(uo.dbInpStream);
#else
       molReader = new RDKit::ForwardSDMolSupplier(
                uo.dbInpStream, takeOwnership, sanitize, removeHs);
#endif
    } else {
        mainErr("Unknown format of db file.");
    }
#ifndef USE_RDKIT
    OpenBabel::OBMol m;
#else
    RDKit::RWMol m;
#endif
    bool done = false;
    while (!done) {
        dbPharm.clear();
#ifndef USE_RDKIT
        m.Clear();
#endif
        if (!uo.dbInpType.empty() && (uo.dbInpType != "phar")) {
#ifndef USE_RDKIT
            if (!molReader->Read(&m)) {
                done = true;
                break;
            } else {
                calcPharm(&m, &dbPharm, uo);
                dbId = m.GetTitle();
            }
#else
            if (molReader->atEnd()) {
                    done = true;
                    break;
                } else {
                std::unique_ptr<RDKit::ROMol> dbmptr(molReader->next());
                if (!dbmptr) {
                    continue;
                }
                m = *dbmptr;
                calcPharm(&m, &dbPharm, uo);
                m.getPropIfPresent("_Name", dbId);
            }
#endif
        } else {
            if (uo.dbInpStream->eof()) {
                done = true;
                break;
            } else {
                dbPharm = pharmReader->read(uo.dbInpStream, dbId);
#ifdef USE_RDKIT
                auto conf = RDKit::Conformer(m.getNumAtoms());
                conf.set3D(true);
                m.addConformer(&conf);
#endif
            }
        }
        if (dbPharm.empty()) {
            continue;
        }
        if (dbId == "") {
            dbId = "Unnamed_ref";
        }

        ++molCount;
        if (!uo.isQuiet ) {
            if ((molCount % 10) == 0) {
                std::cerr << "." << std::flush;
            }
            if ((molCount % 500) == 0) {
                std::cerr << molCount << std::endl << std::flush;
            }
        }
        if (uo.merge) {
            pharMerger.merge(dbPharm);
        }
        if (uo.refInpFile.empty()) {
            if (!(uo.isQuiet)) {
                printProgress(molCount);
            }
            if( !uo.pharmOutFile.empty()) {
                uo.pharmOutWriter->write(dbPharm, uo.pharmOutStream, dbId);
            }
            continue;
        }

        //-------------------------------------------------------
        //...(1).. Alignment
        //-------------------------------------------------------

        dbSize = dbPharm.size();
        dbVolume = 0.0;
        for (unsigned int i(0); i < dbSize; ++i) {
            if (dbPharm[i].func == EXCL) {
                continue;
            }
            dbVolume += VolumeOverlap(dbPharm[i], dbPharm[i], !uo.noNormal);
        }

        // Create a result structure
        Result res;
        res.refId           = refId;
        res.refVolume       = refVolume;
        res.dbId            = dbId;
        res.dbVolume        = dbVolume;
        res.overlapVolume   = 0.0;
        res.exclVolume      = 0.0;
        res.resMol          = m;
        res.resPharSize     = 0;

        if (uo.scoreOnly) {
            FunctionMapping funcMap(&refPharm, &dbPharm, uo.epsilon);
            PharmacophoreMap fMap = funcMap.getNextMap();
            double volBest(-9999.999);

            // loop over all reference points
            while (!fMap.empty()) {
                double newVol(0.0);
                double exclVol(0.0);
                for (PharmacophoreMap::iterator itP = fMap.begin(); itP != fMap.end(); ++itP) {
                    if ((itP->first)->func == EXCL) {
                        exclVol += VolumeOverlap((itP->first), (itP->second), !uo.noNormal);
                    } else if (((itP->first)->func == (itP->second)->func ) ||
                             (((itP->first)->func == HYBH ||
                               (itP->first)->func == HDON ||
                               (itP->first)->func == HACC)
                              && ((itP->second)->func == HDON ||
                                  (itP->second)->func == HACC ||
                                  (itP->second)->func == HYBH))
                             || (((itP->first)->func == HYBL ||
                                  (itP->first)->func == AROM ||
                                  (itP->first)->func == LIPO)
                                 && ((itP->second)->func == AROM ||
                                     (itP->second)->func == LIPO ||
                                     (itP->second)->func == HYBL))) {
                        newVol += VolumeOverlap((itP->first),(itP->second), !uo.noNormal);
                    }
                }
                if ((newVol - exclVol) > volBest) {
                    res.resPhar.clear();
                    res.resPharSize = 0;
                    for (PharmacophoreMap::iterator itP = fMap.begin(); itP != fMap.end(); ++itP) {
                        // add point to resulting pharmacophore
                        PharmacophorePoint p(itP->second);
                        (res.resPhar).push_back(p);
                        ++res.resPharSize;
                    }
                    res.overlapVolume = newVol;
                    res.exclVolume = exclVol;
                    volBest = newVol - exclVol;
                }
                // get the next map
                fMap.clear();
                fMap = funcMap.getNextMap();
            }
        } else {
            FunctionMapping funcMap(&refPharm, &dbPharm, uo.epsilon);
            PharmacophoreMap fMap = funcMap.getNextMap();
            PharmacophoreMap bestMap;

            // default solution
            SolutionInfo best;
            best.volume = -999.9;

            // rotor is set to no rotation
            best.rotor.resize(4);
            best.rotor = 0.0;
            best.rotor[0] = 1.0;

            double bestScore = -1000;
            int mapSize(fMap.size());
            int maxSize = mapSize - 3;


            while (!fMap.empty()) {
                int msize = fMap.size();
                // add the exclusion spheres to the alignment procedure
                if (uo.withExclusion) {
                    for (unsigned int i(0); i < refSize; ++i) {
                        if (refPharm[i].func != EXCL) {
                            continue;
                        }
                        for (unsigned int j(0); j < dbSize; ++j) {
                            if (dbPharm[j].func == EXCL) {
                                continue;
                            }
                            fMap.insert(std::make_pair(&(refPharm[i]), &(dbPharm[j])));
                        }
                    }
                }
                // Only align if the expected score has any chance of being larger
                // than best score so far
                if ((msize > maxSize)
                    && (((double) msize / (refSize - exclSize + dbSize - msize)) > bestScore)) {
                    Alignment align(fMap);
                    SolutionInfo r = align.align(!uo.noNormal);
                    if (best.volume < r.volume) {
                        best = r;
                        bestScore = best.volume / (refVolume + dbVolume - best.volume);
                        bestMap = fMap;
                        mapSize = msize;
                    }
                } else {
                    // Level of mapping site to low
                    break;
                }
                if (bestScore > 0.98) {
                    break;
                }
                // Get the next map
                fMap.clear();
                fMap = funcMap.getNextMap();
            }
            // Transform the complete pharmacophore and the molecule towards the
            // best alignment
            rotMat = quat2Rotation(best.rotor);
            positionPharmacophore(dbPharm, rotMat, best);
            positionMolecule(&res.resMol, rotMat, best);

            // Update result
            res.info = best;

            // Compute overlap volume between exclusion spheres and pharmacophore points
            for (int i(0); i < refSize; ++i) {
                if (refPharm[i].func != EXCL) {
                    continue;
                }
                for (int j(0); j < dbSize; ++j) {
                    res.exclVolume += VolumeOverlap(refPharm[i], dbPharm[j], !uo.noNormal);
                }
            }
            // make copy of the best map and compute the volume overlap
            for (PharmacophoreMap::iterator itP = bestMap.begin(); itP != bestMap.end(); ++itP) {
                if (((itP->first)->func == EXCL) || ((itP->second)->func == EXCL)) {
                    continue;
                }
                // compute overlap
                res.overlapVolume += VolumeOverlap(itP->first, itP->second, !uo.noNormal);
                // add point to resulting pharmacophore
                PharmacophorePoint p(itP->second);
                (res.resPhar).push_back(p);
                ++res.resPharSize;
            }
        }
        // update scores
        res.info.volume = res.overlapVolume - res.exclVolume;
        if (res.info.volume > 0.0) {
            res.tanimoto = res.info.volume / (res.refVolume + res.dbVolume - res.info.volume);
            res.tversky_ref = res.info.volume / res.refVolume;
            res.tversky_db = res.info.volume / res.dbVolume;
        }
        switch (uo.rankby) {
            case TANIMOTO:
                res.rankbyScore = res.tanimoto;
                break;
            case TVERSKY_REF:
                res.rankbyScore = res.tversky_ref;
                break;
            case TVERSKY_DB:
                res.rankbyScore = res.tversky_db;
                break;
        }

        //-------------------------------------------------------
        //...(5).. Generate output
        //-------------------------------------------------------
        if (uo.cutOff != 0.0) {
            if (res.rankbyScore < uo.cutOff) {
                continue;
            }
        }
        if (uo.best != 0) {
            addBest(res, uo, resList);
        } else {
            if (!uo.molOutFile.empty()) {
                logOut(&res, uo);
            }
            if (!uo.pharmOutFile.empty()) {
                logPharmacophores(&res, uo);
            }
            if (!uo.scoreOutFile.empty()) {
                logScores(&res, uo);
            }
        }
    }
    if (molReader) {
        delete molReader;
        molReader = NULL;
    }
    if (pharmReader) {
        delete pharmReader;
        pharmReader = NULL;
    }

    //----------------------------------------------------------------------------
    //...(C).. Process best list (if defined)
    //----------------------------------------------------------------------------

    if (uo.best != 0) {
        std::vector<Result*>::iterator itR;
        for (itR = resList.begin(); itR != resList.end(); ++itR) {
            Result* res(*itR);
            if (!uo.molOutFile.empty()) {
                logOut(res, uo);
            }
            if (!uo.pharmOutFile.empty()) {
                logPharmacophores(res, uo);
            }
            if (!uo.scoreOutFile.empty()) {
                logScores(res, uo);
            }
            delete res;
        }
    }

    // done processing database
    if (!uo.isQuiet) {
        if (uo.refInpFile.empty()) {
            std::cerr << std::endl;
            std::cerr << "Processed " << molCount << " molecules";
            double tt = (double)(clock() - t0 )/CLOCKS_PER_SEC;
            std::cerr << " in " << tt << " seconds (";
            std::cerr << molCount/tt << " molecules per second)." << std::endl;
        }
        else {
            std::cerr << std::endl;
            std::cerr << "Processed " << molCount << " molecules" << std::endl;
            double tt = (double)(clock() - t0 )/CLOCKS_PER_SEC;
            std::cerr << molCount << " alignments in " << tt << " seconds (";
            std::cerr << molCount/tt << " alignments per second)." << std::endl;
        }
    }
    exit(0);
}
