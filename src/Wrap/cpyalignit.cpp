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

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "alignLib.h"
#include "config.h"

namespace python = boost::python;
using namespace RDKit;

namespace {

std::string getVersion() {
  std::string version = "";
  version += std::to_string(ALIGNIT_VERSION);
  version += "." + std::to_string(ALIGNIT_RELEASE);
  version += "." + std::to_string(ALIGNIT_SUBRELEASE);
  return version;
}

python::tuple alignPharmacophore(
    Pharmacophore &ref,
    Pharmacophore &probe,
    ROMol *probeMol = nullptr,
    double epsilon = 0.5,
    bool useNormals = true,
    bool useExclusion = false
) {
  // Perform alignment
  auto res = alignit::alignPharmacophores(
      ref, probe, epsilon, useNormals, useExclusion, probeMol);
  // Update molecules conformer if given
  if (probeMol) {
      const Conformer &conf = res.resMol.getConformer();
      probeMol->clearConformers();
      probeMol->addConformer(new Conformer(conf));
  }
  // Return is (tanimoto, tversky_ref, tversky_db)
  return python::make_tuple(res.tanimoto, res.tversky_ref, res.tversky_db);
}

python::tuple alignMol(
    ROMol &refMol,
    ROMol &probeMol,
    bool calcArom = true,
    bool calcHDon = true,
    bool calcHAcc = true,
    bool calcLipo = true,
    bool calcCharge = true,
    bool calcHybrid = true,
    bool merge = false,
    double epsilon = 0.5,
    bool useNormals = true,
    bool useExclusion = false
) {
    // Perform alignment
    auto out = alignit::alignMols(refMol, probeMol, calcArom, calcHDon,
        calcHAcc, calcLipo, calcCharge, calcHybrid, merge, epsilon,
        useNormals, useExclusion);
    const Result &res = std::get<1>(out);
    // Update molecules conformer
    const Conformer &conf = res.resMol.getConformer();
    probeMol.clearConformers();
    probeMol.addConformer(new Conformer(conf));
    // Return is (tanimoto, tversky_ref, tversky_db)
    return python::make_tuple(res.tanimoto, res.tversky_ref, res.tversky_db);
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
    python::def("CalcPharmacophore", &alignit::calcPharmacophore,
                (python::arg("mol"),
                 python::arg("calcArom") = true,
                 python::arg("calcHDon") = true,
                 python::arg("calcHAcc") = true,
                 python::arg("calcLipo") = true,
                 python::arg("calcCharge") = true,
                 python::arg("calcHybrid") = true),
                "calculate a pharmacophore model for a molecule.");

    // PharmMerger (merge neighboring pharmacophores)
    python::def("MergePharmacophore", &alignit::mergePharmacophore,
                (python::arg("pharmacophore")),
                "merge neighbouring pharmacophore points of the same category");

    // Alignment (pharmacophore->pharmacophore)
    python::def("AlignPharmacophore", &alignPharmacophore,
                (python::arg("ref"), python::arg("probe"),
                 python::args("dbMol") = python::ptr((ROMol*)nullptr),
                 python::args("epsilon") = 0.5,
                 python::args("useNormals") = true,
                 python::args("useExclusion") = false),
                "aligns probe to ref, probe is modified, if the probe "
                "molecule is also given then this will also be modified");

    // Alignment (mol->mol)
    python::def("AlignMol", &alignMol,
                (python::arg("ref"), python::arg("probe"),
                 python::args("calcArom") = true,
                 python::args("calcHDon") = true,
                 python::args("calcHAcc") = true,
                 python::args("calcLipo") = true,
                 python::args("calcCharge") = true,
                 python::args("calcHybrid") = true,
                 python::args("merge") = false,
                 python::args("epsilon") = 0.5,
                 python::args("useNormals") = true,
                 python::args("useExclusion") = false),
                "aligns probe to ref, probe is modified");
}

BOOST_PYTHON_MODULE(cpyalignit) { wrap_pyalignit(); }
