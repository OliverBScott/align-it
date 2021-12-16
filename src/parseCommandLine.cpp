/*******************************************************************************
parseCommandLine.cpp - Align-it

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

#include <parseCommandLine.h>

Options parseCommandLine(int argc, char *argv[]) {
    static struct option Arguments[]  = {
        {"reference", required_argument, NULL, 'r'},
        {"dbase", required_argument, NULL, 'd' },
        {"scores", required_argument, NULL, 's' },
        {"out", required_argument, NULL, 'o' },
        {"pharmacophore", required_argument, NULL, 'p' },
        {"funcGroup", required_argument, NULL, 'f' },
        {"epsilon", required_argument, NULL, 'e' },
        {"merge", no_argument, NULL, 'm' },
        {"noNormal", no_argument, NULL, 'n' },
        {"help", no_argument, NULL, 'h' },
        {"version", no_argument, NULL, 'v' },
        {"quiet", no_argument, NULL, 'q' },
        {"refType", required_argument, NULL, 1  },
        {"dbType", required_argument, NULL, 2  },
        {"cutOff", required_argument, NULL, 3  },
        {"best", required_argument, NULL, 4  },
        {"rankBy", required_argument, NULL, 5  },
        {"noHybrid", no_argument, NULL, 7  },
        {"info", required_argument, NULL, 9  },
        {"withExclusion", no_argument, NULL, 10 },
        {"scoreOnly",  no_argument, NULL, 11 },
        {"outType", required_argument, NULL,	12 },
        {NULL, 0, NULL, 0}
    };
    Options o;

    // Set defaults
    o.dbInpFile.clear();
    o.dbInpStream = NULL;
    o.dbInpType.clear();

    o.refInpFile.clear();
    o.refInpStream = NULL;
    o.refInpType.clear();

    o.molOutFile.clear();
    o.molOutStream = NULL;
    //o.molOutWriter = NULL;
    o.molOutType.clear();

    o.pharmOutFile.clear();
    o.pharmOutStream = NULL;
    o.pharmOutWriter = NULL;

    o.scoreOutFile.clear();
    o.scoreOutStream = NULL;

    o.epsilon = 0.5;
    o.cutOff = 0.0;
    o.best = 0;
    o.rankby = TANIMOTO;

    o.funcGroupVec.resize(10);
    o.funcGroupVec[AROM] = true;
    o.funcGroupVec[HDON] = true;
    o.funcGroupVec[HACC] = true;
    o.funcGroupVec[LIPO] = true;
    o.funcGroupVec[POSC] = true;
    o.funcGroupVec[NEGC] = true;
    o.funcGroupVec[HYBH] = false;
    o.funcGroupVec[HYBL] = false;

    o.isQuiet = false;
    o.noHybrid = false;
    o.merge = false;
    o.noNormal = false;
    o.withExclusion = false;
    o.scoreOnly = false;
    o.version = false;

    int choice;
    opterr = 0;
    int optionIndex = 0;
    std::string strvalue, t;
    std::list<std::string> l;
    std::list<std::string>::iterator itL;
    std::string ext;

    while ((choice = getopt_long(argc, argv, "vhqr:d:s:o:p:f:e:m", Arguments, &optionIndex)) != -1) {
        switch (choice) {
            case 'v': //......................................................version
                o.version = true;
                break;
            case 'r':  //...................................................reference
                o.refInpFile = optarg;
                o.refInpStream = new std::ifstream(optarg);
                if (!o.refInpStream->good()) {
                    mainErr("Error opening input file for reference (-r | --reference)");
                }
                break;
            case 'd':  //.......................................................dbase
                o.dbInpFile = optarg;
                o.dbInpStream = new std::ifstream(optarg);
                if (!o.dbInpStream->good()) {
                    mainErr("Error opening input file for database (-d | --dbase)");
                }
                break;
            case 's': //.......................................................scores
                o.scoreOutFile = optarg;
                o.scoreOutStream = new std::ofstream(optarg);
                if (!o.scoreOutStream->good()) {
                    mainErr("Error opening output file for scores (-s)");
                }
                break;
            case 'o': //..........................................................out
                o.molOutFile = optarg;
                o.molOutStream = new std::ofstream(optarg);
                if (!o.molOutStream->good()) {
                    mainErr("Error opening output file for molecules (-o | --out)");
                }
                break;
            case 'p': //................................................pharmacophore
                o.pharmOutFile = optarg;
                o.pharmOutStream = new std::ofstream(optarg);
                if (!o.pharmOutStream->good()) {
                    mainErr("Error opening output file for pharmacophores (-p)");
                }
                o.pharmOutWriter = new PharmacophoreWriter();
                break;
            case 'e': //......................................................epsilon
                o.epsilon = strtod(optarg,NULL);
                break;
            case 'm': //........................................................merge
                o.merge = true;
                o.noNormal = true;
                break;
            case 'n': //.....................................................noNormal
                o.noNormal = true;
                break;
            case 'f': { //.................................................funcGroup
                l = stringTokenizer(optarg, ",");;
                std::vector<bool> vec(10, false);
                for (itL = l.begin(); itL != l.end(); ++itL) {
                    if (*itL == "AROM") {
                        vec[AROM] = true;
                        continue;
                    }
                    if (*itL == "HDON") {
                        vec[HDON] = true;
                        continue;
                    }
                    if (*itL == "HACC") {
                        vec[HACC] = true;
                        continue;
                    }
                    if (*itL == "LIPO") {
                        vec[LIPO] = true;
                        continue;
                    }
                    if (*itL == "CHARGE") {
                        vec[POSC] = true;
                        vec[NEGC] = true;
                        continue;
                    }
                    mainErr("Undefined functional Group. Only AROM, HDON, HACC, LIPO and"
                            "CHARGE are allowed as argument.");
                }
                o.funcGroupVec = vec;
                break;
            }
            case 1: //........................................................refType
                o.refInpType = optarg;
                break;
            case 2: //.........................................................dbType
                o.dbInpType = optarg;
                break;
            case 3: //.........................................................cutOff
                o.cutOff = strtod(optarg, NULL);
                break;
            case 4: //...........................................................best
                o.best = strtol(optarg, NULL, 10);
                break;
            case 5: //.........................................................rankby
                t = optarg;
                if (t == "TANIMOTO") {
                    o.rankby = TANIMOTO;
                } else if (t == "TVERSKY_REF") {
                    o.rankby = TVERSKY_REF;
                } else if (t == "TVERSKY_DB") {
                    o.rankby = TVERSKY_DB;
                } else {
                    mainErr("Undefined rankby type : " + t);
                    break;
                }
            case 7: //.......................................................noHybrid
                o.noHybrid = true;
                break;
            case 'h': //.........................................................help
                o.help = true;
                break;
            case 9: //...........................................................info
                printInfo(std::string(optarg));
                break;
            case 10: //.................................................withExclusion
                o.withExclusion = true;
                break;
            case 11: //.....................................................scoreOnly
                o.scoreOnly = true;
                break;
            case 12: //.......................................................outType
                o.molOutType = optarg;
                break;
            case 'q': //........................................................quiet
                o.isQuiet = true;
                break;
            default:
                mainErr("unknown command-line option");
        }
    }
    // If no options are given print the help
    if (optind == 1) {
        o.help = true;
    } else {
        // REFERENCE
        if (!o.refInpType.empty()) {
            // Get type from specification
            std::string s = o.refInpType;
            std::string::iterator i = s.begin();
            std::string::iterator end = s.end();
            while (i != end) {
                *i = std::tolower((unsigned char)*i);
                ++i;
            }
            o.refInpType = s;
#ifndef USE_RDKIT
            if (o.refInpType != "phar") {
                OpenBabel::OBConversion* reader = new OpenBabel::OBConversion();
				if (reader->FindFormat(o.refInpType) != NULL) {
					o.refInpType = reader->FindFormat(o.refInpType)->GetID();
				} else {
					o.refInpType.clear();
					mainErr("Unknown format of reference file.");
				}
				delete reader;
            }
#else
            if (o.refInpType != "phar" && o.refInpType != "sdf") {
                mainErr("RDKit implementation currently only supports SDF or PHAR (refInpType)");
            }
#endif
        } else {
            // Get type from reference file extension
#ifndef USE_RDKIT
            OpenBabel::OBConversion* reader = new OpenBabel::OBConversion();
			if (reader->FormatFromExt(o.refInpFile.c_str()) != NULL) {
				o.refInpType = reader->FormatFromExt(o.refInpFile.c_str())->GetID();
			} else {
				o.refInpType.clear();
			}
			delete reader;
#else
            std::string extension = getExt(o.refInpFile);
            if (extension != "sdf" && extension != "phar") {
                mainErr("RDKit implementation currently only supports SDF or PHAR (refInpFile)");
            }
            o.refInpType = extension;
#endif
        }
        // DATABASE
        if (!o.dbInpType.empty()) {
            // Get type from specification
            std::string s = o.dbInpType;
            std::string::iterator i = s.begin();
            std::string::iterator end = s.end();
            while (i != end) {
                *i = std::tolower((unsigned char)*i);
                ++i;
            }
            o.dbInpType = s;
#ifndef USE_RDKIT
            if (o.dbInpType != "phar") {
                OpenBabel::OBConversion* reader = new OpenBabel::OBConversion();
				if (reader->FindFormat(o.dbInpType) != NULL) {
					o.dbInpType = reader->FindFormat(o.dbInpType)->GetID();
				} else {
					o.dbInpType.clear();
					mainErr("Unknown format of reference file.");
				}
				delete reader;
            }
#else
            if (o.dbInpType != "phar" && o.dbInpType != "sdf") {
                mainErr("RDKit implementation currently only supports SDF or PHAR (dbInpType)");
            }
#endif
        } else {
            // Get type from reference file extension
#ifndef USE_RDKIT
            OpenBabel::OBConversion* reader = new OpenBabel::OBConversion();
			if (reader->FormatFromExt(o.dbInpFile.c_str()) != NULL) {
				o.dbInpType = reader->FormatFromExt(o.dbInpFile.c_str())->GetID();
			} else {
				o.dbInpType.clear();
			}
			delete reader;
#else
            std::string extension = getExt(o.dbInpFile);
            if (extension != "sdf" && extension != "phar") {
                mainErr("RDKit implementation currently only supports SDF or PHAR (dbInpFile)");
            }
            o.dbInpType = extension;
#endif
        }
        // MOL OUTPUT
        if (!o.molOutType.empty()) {
            // Get type from specification
            std::string s = o.molOutType;
            std::string::iterator i = s.begin();
            std::string::iterator end = s.end();
            while (i != end) {
                *i = std::tolower((unsigned char) *i);
                ++i;
            }
            o.molOutType = s;
#ifndef USE_RDKIT
            OpenBabel::OBConversion* reader = new OpenBabel::OBConversion();
			if (reader->FindFormat(o.molOutType) != NULL) {
				o.molOutType = reader->FindFormat(o.molOutType)->GetID();
			} else {
				o.molOutType.clear();
				mainErr("Unknown format of output file.");
			}
			delete reader;
#else
            if (o.molOutType != "sdf") {
                mainErr("RDKit implementation currently only supports SDF (molOutType)");
            }
#endif
        } else {
            // Get type from output file extension
#ifndef  USE_RDKIT
            OpenBabel::OBConversion* reader = new OpenBabel::OBConversion();
			if (reader->FormatFromExt(o.molOutFile.c_str()) != NULL) {
				o.molOutType = reader->FormatFromExt(o.molOutFile.c_str())->GetID();
			} else {
				o.molOutType.clear();
			}
			delete reader;
#else
            if (!o.molOutFile.empty()) {
                std::string extension = getExt(o.molOutFile);
                if (extension != "sdf") {
                    mainErr("RDKit implementation currently only supports SDF (molOutFile)");
                }
                o.molOutType = "sdf";
            } else {
                o.molOutType.clear();
            }
#endif
        }
#ifndef USE_RDKIT
        if (!o.molOutFile.empty() && !o.molOutType.empty()) {
            o.molOutWriter = new OpenBabel::OBConversion();
            o.molOutWriter->SetOutFormat(o.molOutWriter->FindFormat(o.molOutType));
        } else if (!o.molOutFile.empty()) {
            o.molOutWriter = new OpenBabel::OBConversion();
            o.molOutWriter->SetOutFormat(o.molOutWriter->FormatFromExt(o.molOutFile));
        }
#else
        if (!o.molOutFile.empty() && !o.molOutType.empty()) {
            o.molOutWriter = new RDKit::SDWriter(o.molOutStream, false);
        }
#endif
    }
    argc -= optind;
    argv += optind;
    return o;
}
