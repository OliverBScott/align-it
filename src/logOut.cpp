/*******************************************************************************
logOut.cpp - Align-it

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

#include <logOut.h>

#ifndef USE_RDKIT
#include <openbabel/generic.h>

void logOut(Result *res, Options &uo) {
   if((res->resMol).Empty()) {
      return;
   }

   std::ostringstream ss;

	// Add properties to the molecule
   OpenBabel::OBPairData* label1 = new OpenBabel::OBPairData();
   label1->SetAttribute("PHARAO_TANIMOTO");
   ss << res->tanimoto;
   label1->SetValue(ss.str());
   (res->resMol).SetData(label1);

   OpenBabel::OBPairData* label2 = new OpenBabel::OBPairData();
   label2->SetAttribute("PHARAO_TVERSKY_REF");
   ss.str("");
   ss << res->tversky_ref;
   label2->SetValue(ss.str());
   (res->resMol).SetData(label2);

   OpenBabel::OBPairData* label3 = new OpenBabel::OBPairData();
   label3->SetAttribute("PHARAO_TVERSKY_DB");
   ss.str("");
   ss << res->tversky_db;
   label3->SetValue(ss.str());
   (res->resMol).SetData(label3);

   // Write molecule
   uo.molOutWriter->Write(&(res->resMol), uo.molOutStream);

   // Clean up
   (res->resMol).DeleteData("PHARAO_TANIMOTO");
   (res->resMol).DeleteData("PHARAO_TVERSKY_REF");
   (res->resMol).DeleteData("PHARAO_TVERSKY_DB");
}

#else

void logOut(Result *res, Options &uo) {
    if (res->resMol.getNumAtoms() == 0) return;
    res->resMol.setProp("PHARAO_TANIMOTO", res->tanimoto);
    res->resMol.setProp("PHARAO_TVERSKY_REF", res->tversky_ref);
    res->resMol.setProp("PHARAO_TVERSKY_DB", res->tversky_db);
    uo.molOutWriter->write(res->resMol);
}

#endif
