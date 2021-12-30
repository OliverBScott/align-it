/*******************************************************************************
printInfo.cpp - Align-it

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

#include <printInfo.h>

void printInfo(const std::string &option) {
    bool valid(false);
    printHeader();
    if (option == "r" || option == "reference") {
        std::cerr << "  -r  --reference" << std::endl;
        std::cerr << std::endl;
        std::cerr << " Defines the reference structure that will be used to screen and/or" << std::endl;
        std::cerr << " align the database. This option is not required, so when not given," << std::endl;
        std::cerr << " only the database will be converted into pharmacophores without alignment." << std::endl;
        std::cerr << std::endl;
        std::cerr << " By default the format is deduced from the extension of the file" << std::endl;
        std::cerr << " but can be defined explicitly with the --refType option." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "refType") {
        std::cerr << "  --refType" << std::endl;
        std::cerr << std::endl;
        std::cerr << " With this option the format of the reference input file can be" << std::endl;
        std::cerr << " specified. The following formats are supported:" << std::endl;
        std::cerr << std::endl;
        std::cerr << "   <type> : The <type> argument is one of the types as specified by OpenBabel" << std::endl;
        std::cerr << "            The molecule should contain coordinate information" << std::endl;
        std::cerr << "   'PHAR' : The molecule is already transformed in pharmacophore" << std::endl;
        std::cerr << "            format." << std::endl;
        std::cerr << std::endl;
        std::cerr << " If the pharmacophore format is used, the program will generate a" << std::endl;
        std::cerr << " pharmacophore for the reference molecule. The time needed for this" << std::endl;
        std::cerr << " generation is however negligible compared with the time needed for" << std::endl;
        std::cerr << " a proper alignment." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "d" || option=="dbase") {
        std::cerr << "  -d  --dbase" << std::endl;
        std::cerr << std::endl;
        std::cerr << " Defines the database, or a collection of structures, that will be" << std::endl;
        std::cerr << " used to screen. This option is required but it is allowed to" << std::endl;
        std::cerr << " select a database containing only a single molecule." << std::endl;
        std::cerr << std::endl;
        std::cerr << " The database can consist of precalculated pharmacophores or of" << std::endl;
        std::cerr << " molecules in a format readable by the computer." << std::endl;
        std::cerr << std::endl;
        std::cerr << " By default the format is deduced from the extension of the file" << std::endl;
        std::cerr << " but can be defined explicitly with the --dbType option." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "dbType") {
        std::cerr << "  --dbType" << std::endl;
        std::cerr << std::endl;
        std::cerr << " With this option the format of the database input file can be" << std::endl;
        std::cerr << " specified. The following formats are supported: " << std::endl;
        std::cerr << std::endl;
        std::cerr << "   <type> : The <type> argument is one of the types as specified by OpenBabel" << std::endl;
        std::cerr << "            Molecules should contain coordinate information" << std::endl;
        std::cerr << "   'PHAR' : The molecules are already transformed in pharmacophore" << std::endl;
        std::cerr << "            format." << std::endl;
        std::cerr << std::endl;
        std::cerr << " If the pharmacophore format is not used, the program will generate" << std::endl;
        std::cerr << " a pharmacophore for each molecule on-the-fly. The time needed for" << std::endl;
        std::cerr << " this generation is however negligible compared with the time" << std::endl;
        std::cerr << " needed for a proper alignement." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "p" || option == "pharmacophore") {
        std::cerr << "  -p  --pharmacophore" << std::endl;
        std::cerr << std::endl;
        std::cerr << " In this file the processed pharmacophores of the structures in" << std::endl;
        std::cerr << " the input database are stored. These pharmacophores will not" << std::endl;
        std::cerr << " correspond to the original structures because they are aligned" << std::endl;
        std::cerr << " with respect to the reference input molecule and therefore can" << std::endl;
        std::cerr << " have a different orientation. Moreover, only points used in the" << std::endl;
        std::cerr << " alignment are reported." << std::endl;
        std::cerr << std::endl;
        std::cerr << " For more information about the format in which the pharmacophores" << std::endl;
        std::cerr << " are written, consult the align-it manual." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "o" || option == "out") {
        std::cerr << "  -o  --out" << std::endl;
        std::cerr << std::endl;
        std::cerr << " In this file the transformed database structures after aligning" << std::endl;
        std::cerr << " them to the reference structure can be written down. These" << std::endl;
        std::cerr << " structures correspond to the processed pharmacophores. " << std::endl;
        std::cerr << std::endl;
        std::cerr << " By default the format is deduced from the extension of the file" << std::endl;
        std::cerr << " but can be defined explicitly with the --outType option." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "outType") {
        std::cerr << "  --outType" << std::endl;
        std::cerr << std::endl;
        std::cerr << " With this option the format of the molecule output file can be" << std::endl;
        std::cerr << " specified. The formats as specified by OpenBabel are recognised." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "s" || option == "scores") {
        std::cerr << "  -s  --scores" << std::endl;
        std::cerr << std::endl;
        std::cerr << " With this option a tab-delimited text file can be generated" << std::endl;
        std::cerr << " containing the following information:" << std::endl;
        std::cerr << std::endl;
        std::cerr << "   'column 1'  : Id of the reference structure." << std::endl;
        std::cerr << "   'column 2'  : Maximum volume of the reference structure." << std::endl;
        std::cerr << "   'column 3'  : Id of the database structure." << std::endl;
        std::cerr << "   'column 4'  : Maximum volume of the database structure." << std::endl;
        std::cerr << "   'column 5'  : Maximum pharmacophore overlap of the two structures." << std::endl;
        std::cerr << "   'column 6'  : Overlap between pharmacophore and exclusion spheres" << std::endl;
        std::cerr << "                 in the reference." << std::endl;
        std::cerr << "   'column 7'  : Corrected volume overlap between database" << std::endl;
        std::cerr << "                 pharmacophore and reference." << std::endl;
        std::cerr << "   'column 8'  : Number of pharmacophore points in the processed" << std::endl;
        std::cerr << "                 pharmacophore." << std::endl;
        std::cerr << "   'column 9'  : TANIMOTO score. This value is between 0 and 1." << std::endl;
        std::cerr << "   'column 10' : TVERSKY_REF score. This value is between 0 and 1." << std::endl;
        std::cerr << "   'column 11' : TVERSKY_DB score. This value is between 0 and 1." << std::endl;
        std::cerr << std::endl;
        std::cerr << "  More information about the three different scoring schemes can be" << std::endl;
        std::cerr << "  obtained by ty typing: Align-it --info rankby " << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "cutOff") {
        std::cerr << "  --cutOff" << std::endl;
        std::cerr << std::endl;
        std::cerr << " This value should be between 0 and 1 and only structures with a" << std::endl;
        std::cerr << " score, as defined by the --rankBy option, higher than this value" << std::endl;
        std::cerr << " are reported in the three possible output files." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "best") {
        std::cerr << "  --best" << std::endl;
        std::cerr << std::endl;
        std::cerr << " With this option only a limited number of best scoring structures," << std::endl;
        std::cerr << " as defined by the --rankBy option, are reported in the three" << std::endl;
        std::cerr << " possible output files. If the --cutOff option was used, all best " << std::endl;
        std::cerr << " scoring structures must first pass this filter. The user can " << std::endl;
        std::cerr << " specify the number of best scoring structures that should be " << std::endl;
        std::cerr << " reported." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "rankBy") {
        std::cerr << "  --rankBy" << std::endl;
        std::cerr << std::endl;
        std::cerr << " This option defines the scoring used by the --cutOff and --best" << std::endl;
        std::cerr << " options. Currently three measures are implemented:" << std::endl;
        std::cerr << std::endl;
        std::cerr << "   'TANIMOTO'    : Final overlap score. This value is between 0 and" << std::endl;
        std::cerr << "                   1, and the closer this value to 1, the better the" << std::endl;
        std::cerr << "                   alignment." << std::endl;
        std::cerr << "   'TVERSKY_REF' : Overlap ref volume ratio. This value is between 0" << std::endl;
        std::cerr << "                   and 1, and the closer this value to 1, the better" << std::endl;
        std::cerr << "                   the reference pharmacophore is 'included' in the" << std::endl;
        std::cerr << "                   database pharmacophore, regardless the unmatched" << std::endl;
        std::cerr << "                   part of the database pharmacophore." << std::endl;
        std::cerr << "   'TVERSKY_DB' :  Overlap db volume ratio. This value is between 0" << std::endl;
        std::cerr << "                   and 1, and the closer this value to 1, the better" << std::endl;
        std::cerr << "                   the database pharmacophore is 'included' in the" << std::endl;
        std::cerr << "                   reference pharmacophore, regardless the unmatched" << std::endl;
        std::cerr << "                   part of the reference pharmacophore." << std::endl;
        std::cerr << std::endl;
        std::cerr << " By default the TANIMOTO score is used. For more information" << std::endl;
        std::cerr << " consult the Align-it manual." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "f" || option == "funcGroup") {
        std::cerr << "  -f  --funcGroup" << std::endl;
        std::cerr << std::endl;
        std::cerr << " By default all generated pharmacophores contain all functional" << std::endl;
        std::cerr << " groups and thus include all information that might be useful. With" << std::endl;
        std::cerr << " this option only a subset of the available functional groups can" << std::endl;
        std::cerr << " be used in the alignment. The user can define this subset by using" << std::endl;
        std::cerr << " the tags listed below with the ',' symbol as separator." << std::endl;
        std::cerr << std::endl;
        std::cerr << "   'AROM'   : Aromatic rings" << std::endl;
        std::cerr << "   'HDON'   : Hydrogen donors" << std::endl;
        std::cerr << "   'HACC'   : Hydrogen acceptor" << std::endl;
        std::cerr << "   'LIPO'   : Lipophilic spots" << std::endl;
        std::cerr << "   'CHARGE' : Charge centers (both positive and negative)" << std::endl;
        std::cerr << std::endl;
        std::cerr << " If both reference and database structures are defined in the" << std::endl;
        std::cerr << " pharmacophore format this option is discarded and has no influence" << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "e" || option == "epsilon") {
        std::cerr << "  -e  --epsilon" << std::endl;
        std::cerr << std::endl;
        std::cerr << " This option can be used to change the level of relative overlap between" << std::endl;
        std::cerr << " points to be considered as feasible combinations. As such, it is an important" << std::endl;
        std::cerr << " parameter in the 'feature mapping'." << std::endl;
        std::cerr << std::endl;
        std::cerr << " The value should be set between 0.0 and 1.0 and indicates how much overlap" << std::endl;
        std::cerr << " is required for two points to be considered as a feasible combination. The relative overlap" << std::endl;
        std::cerr << " is computed by mapping one database point on top of one reference point and measuring the " << std::endl;
        std::cerr << " overlap between the two other points. Smaller values are less stringent, but require" << std::endl;
        std::cerr << " more computing time." << std::endl;
        std::cerr << " The default value is 0.5" << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "m" || option == "merge") {
        std::cerr << "  -m  --merge" << std::endl;
        std::cerr << std::endl;
        std::cerr << " This option can be used to merge pharmacophore points that are" << std::endl;
        std::cerr << " close enough to each other together into a single pharmacophore" << std::endl;
        std::cerr << " point with the same functional group, but with an increased sigma." << std::endl;
        std::cerr << std::endl;
        std::cerr << " This flag also activates the -n or --noNormal flag because merged" << std::endl;
        std::cerr << " pharmacophore points can't have a 'direction'." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "n" || option == "noNormal") {
        std::cerr << "  -n  --noNormal" << std::endl;
        std::cerr << std::endl;
        std::cerr << " Flag to indicate that no normal information is used during alignment." << std::endl;
        std::cerr << std::endl;
        std::cerr << " By default, several pharmacophore features contain normal" << std::endl;
        std::cerr << " information to add the notion of 'direction' in addition to" << std::endl;
        std::cerr << " 'position'. AROM, HYBL, HDON, HACC and HYBH points all have such an" << std::endl;
        std::cerr << " additional constraint. Using this information makes the" << std::endl;
        std::cerr << " pharmacophore model more specific." << std::endl;
        std::cerr << std::endl;
        std::cerr << " When the -m or --merge flag is used this flag is automatically" << std::endl;
        std::cerr << " activated." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "noHybrid") {
        std::cerr << "  --noHybrid" << std::endl;
        std::cerr << std::endl;
        std::cerr << " Flag to indicate that no hybrid points should be included in the" << std::endl;
        std::cerr << " final pharmacophores. The possible hybrid pharmacophore points are:" << std::endl;
        std::cerr << std::endl;
        std::cerr << "   'HYBL' : AROM + LIPO" << std::endl;
        std::cerr << "   'HYBH' : HDON + HACC" << std::endl;
        std::cerr << std::endl;
        std::cerr << " Hybrid pharmacophore points are generated by default to reduce the" << std::endl;
        std::cerr << " number of pharmacophore points." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "scoreOnly") {
        std::cerr << "  --scoreOnly" << std::endl;
        std::cerr << std::endl;
        std::cerr << " Flag to indicate when the volume overlap should only be computed from " << std::endl;
        std::cerr << " the given poses and do not perform any translation or rotation of the " << std::endl;
        std::cerr << " pharmacophore points to optimize the overlap." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "withExclusion") {
        std::cerr << "  --withExclusion" << std::endl;
        std::cerr << std::endl;
        std::cerr << " Flag to indicate if the exclusion spheres should be part of" << std::endl;
        std::cerr << " the optimization procedure. By default the overlap between " << std::endl;
        std::cerr << " pharmacophore and exclusion spheres is only taken into account" << std::endl;
        std::cerr << " at the end of the alignment procedure. When this flag is set," << std::endl;
        std::cerr << " the exclusion spheres have an imapct on the optimization procedure." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "h" || option == "help") {
        std::cerr << "  -h  --help" << std::endl;
        std::cerr << std::endl;
        std::cerr << " With this option a general help on the usage of Align-it is provided." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "v") {
        std::cerr << "  -v" << std::endl;
        std::cerr << std::endl;
        std::cerr << " Prints the program version." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "info") {
        std::cerr << "  --info" << std::endl;
        std::cerr << std::endl;
        std::cerr << " With this option the user can get detailed information for each" << std::endl;
        std::cerr << " option." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (option == "q" || option == "quiet") {
        std::cerr << "  -q  --quiet" << std::endl;
        std::cerr << std::endl;
        std::cerr << " If no output is needed during the execution of the program this " << std::endl;
        std::cerr << " option can be used and no output, progress or warnings are shown " << std::endl;
        std::cerr << " to the user." << std::endl;
        std::cerr << std::endl;
        valid = true;
    }
    if (!valid) {
        std::cerr << " Unknown option: " << option << std::endl;
        std::cerr << std::endl;
    }
    exit(0);
}
