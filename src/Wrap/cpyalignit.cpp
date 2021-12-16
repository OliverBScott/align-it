/*********************************************************************************

Copyright 2021 by OliverBScott and the Shape-it contributors

This file is part of Align-it.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

***********************************************************************/

#include <iostream>
#include <fstream>

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "alignLib.h"
#include "config.h"

namespace python = boost::python;

namespace {

std::string getVersion() {
    std::string version = "";
    version += std::to_string(ALIGNIT_VERSION);
    version += "." + std::to_string(ALIGNIT_RELEASE);
    version += "." + std::to_string(ALIGNIT_SUBRELEASE);
    return version;
}

} // namespace

void wrap_pyalignit() {
    // Align-it version
    python::def("GetVersion", &getVersion);

    // Coordinate (the classic xyz)
    python::class_<Coordinate>("Coordinate")
        .def_readwrite("x", &Coordinate::x)
        .def_readwrite("y", &Coordinate::y)
        .def_readwrite("z", &Coordinate::z);

    // FuncGroups (enum of possible functional groups)
    python::enum_<FuncGroup>("FuncGroup")
        .value("AROM", FuncGroup::AROM)
        .value("HDON", FuncGroup::HDON)
        .value("HACC", FuncGroup::HACC)
        .value("LIPO", FuncGroup::LIPO)
        .value("POSC", FuncGroup::POSC)
        .value("NEGC", FuncGroup::NEGC)
        .value("HYBH", FuncGroup::HYBH)
        .value("HYBL", FuncGroup::HYBL)
        .value("EXCL", FuncGroup::EXCL)
        .value("UNDEF", FuncGroup::UNDEF);

    // PharmacophorePoint (represents a functional group)
    python::class_<PharmacophorePoint>("PharmacophorePoint")
        .add_property("point", &PharmacophorePoint::point)
        .add_property("normal", &PharmacophorePoint::normal)
        .def_readwrite("func", &PharmacophorePoint::func)
        .def_readwrite("alpha", &PharmacophorePoint::alpha)
        .def_readwrite("hasNormal", &PharmacophorePoint::hasNormal);

    // Pharmacophore (pharmacophore model i.e. a vector of pharmacophore points)
    python::class_<Pharmacophore>("Pharmacophore")
        .def(python::vector_indexing_suite<std::vector<PharmacophorePoint>>());

    // CalcPharm (calculate pharmacophores)
    python::def("CalcPharmacophore", &alignit::calcPharmacophore, (
            python::arg("mol"),
            python::arg("calcArom") = true,
            python::arg("calcHDon") = true,
            python::arg("calcHAcc") = true,
            python::arg("calcLipo") = true,
            python::arg("calcCharge") = true,
            python::arg("calcHybrid") = true),
            "calculate a pharmacophore model for a molecule."
    );

    // PharmMerger (merge neighboring pharmacophores)
    python::def("MergePharmacophore", &alignit::mergePharmacophore,
            (python::arg("pharmacophore")),
            "merge neighbouring pharmacophore points of the same category"
    );

    // Alignment


}

BOOST_PYTHON_MODULE(cpyalignit) { wrap_pyalignit(); }
