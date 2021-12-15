/*******************************************************************************
stringTokenizer.cpp - Align-it

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

#include <stringTokenizer.h>

std::list<std::string> stringTokenizer(const std::string &s, const std::string &spacer) {
    const std::string::size_type len = s.length();
    std::string::size_type i = 0;
    std::list<std::string> container;
    container.clear();
    while (i < len) {
        // consume leading whitespace
        i = s.find_first_not_of(spacer, i);
        if (i == std::string::npos) {
            return container;  // nothing left but whitespace
        }
        // find the end of the token
        std::string::size_type j = s.find_first_of(spacer, i);
        // push token into container
        if (j == std::string::npos) { // end of string
            container.push_back(s.substr(i));
            return container;
        } else {
            container.push_back(s.substr(i, j-i));
        }
        i = j + 1;  // move positions to next substring
    }
    return container;
}
