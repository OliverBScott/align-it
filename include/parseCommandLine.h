/*******************************************************************************
parseCommandLine.h - Align-it

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

#ifndef __SILICOSIT_ALIGNIT_PARSECOMMANDLINE_H__
#define __SILICOSIT_ALIGNIT_PARSECOMMANDLINE_H__

// General
#include <list>
#include <getopt.h>
#include <stdlib.h>
#include <map>

// Align-it
#include <options.h>
#include <printInfo.h>
#include <mainErr.h>
#include <stringTokenizer.h>
#include <pharmacophore.h>
#include <getExt.h>

Options parseCommandLine(int argc, char* argv[]);

#endif //__SILICOSIT_ALIGNIT_PARSECOMMANDLINE_H__
