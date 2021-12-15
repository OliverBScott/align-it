/*******************************************************************************
lipoFuncCalc.cpp - Align-it

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

#include <lipoFuncCalc.h>

#ifndef USE_RDKIT

void lipoFuncCalc(OpenBabel::OBMol* m, Pharmacophore* pharmacophore) {
   // Make a copy of the partial charge of each atom, since here we store the values
   std::vector<double> pcharge;
   pcharge.clear();
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai)) {
      pcharge.push_back(a->GetPartialCharge());
   }

   // Find for each atom the 'topology-dependent term': t
   _lipoLabelAtoms(m);

   // Multiply 't' with the accessible surface 's'
   for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai)) {
      if (a->GetAtomicNum() == 1) {
         continue;
      }
      double t = a->GetPartialCharge();
      if (t != 0.0) {
         a->SetPartialCharge(_lipoCalcAccSurf(a) * t);
      }
   }

   // Finally calculate the lipophilic points
   _lipoGroupAtoms(m, pharmacophore);

   // Reset partial charges
   unsigned int i(0);
   for (OpenBabel::OBAtom* atom = m->BeginAtom(ai); atom; atom = m->NextAtom(ai)) {
      atom->SetPartialCharge(pcharge[i]);
      ++i;
   }
}

void _lipoLabelAtoms(OpenBabel::OBMol* m) {
  // Give all atoms default score of 1
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* atom = m->BeginAtom(ai); atom; atom = m->NextAtom(ai)) {
      atom->SetPartialCharge(1.0);
   }

   // Decrease scores based on O,N,S and charged atoms
   std::vector<OpenBabel::OBBond*>::iterator bi;
   std::vector<OpenBabel::OBBond*>::iterator bi2;
   for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai)) {
      switch (a->GetAtomicNum()) {
         //.........................................................
         case 1: //.................................................
            a->SetPartialCharge(0.0); // category 1
            break;
         //.........................................................
         case 7: //.................................................
            a->SetPartialCharge(0.0); // category 1
            if (a->IsAromatic()) {
               break;
            }
            _lipoLabelNeighbors(a, 0.25); // category 13
            if ((a->GetImplicitHCount() - a->GetHvyDegree()) !=0 ) {
               std::vector<OpenBabel::OBBond*>::iterator bi;
               for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi)) {
                  OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
                  aa->SetPartialCharge(0.0);     // category 4
                  _lipoLabelNeighbors(aa, 0.0);  // category 4
               }
            }
            break;
         //.........................................................
         case 8: //.................................................
            a->SetPartialCharge(0.0); // category 1
            if (a->IsAromatic()) {
               break;
            }
            _lipoLabelNeighbors(a, 0.25); // category 13

            for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi)) {
               if (b->GetNbrAtom(a)->GetAtomicNum() == 1) {
                  std::vector<OpenBabel::OBBond*>::iterator bi2;
                  for (OpenBabel::OBBond* b2 = a->BeginBond(bi2); b2; b2 = a->NextBond(bi2)) {
                     OpenBabel::OBAtom* aa = b2->GetNbrAtom(a);
                     aa->SetPartialCharge(0.0);    // category 4
                     _lipoLabelNeighbors(aa, 0.0); // category 4
                  }
               }
               if (b->GetBondOrder() == 2) {
                  OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
                  aa->SetPartialCharge(0.0); // category 6
                  for (OpenBabel::OBBond* b2 = aa->BeginBond(bi2); b2; b2 = aa->NextBond(bi2)) {
                     OpenBabel::OBAtom* aaa = b2->GetNbrAtom(aa);
                     if (aaa == a) {
                        continue;
                     }
                     aaa->SetPartialCharge(0.0);    // category 6
                     _lipoLabelNeighbors(aaa, 0.6); // category 9
                  }
               }
            }
            break;
         //.........................................................
         case 16: //................................................
            for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi)) {
               if (b->GetNbrAtom(a)->GetAtomicNum() == 1) {
                  a->SetPartialCharge(0.0);    // category 2
                  _lipoLabelNeighbors(a, 0.0); // category 5
               }
               if (b->GetBondOrder() == 2) {
                  a->SetPartialCharge(0.0);    // category 8
                  _lipoLabelNeighbors(a, 0.6); // category 11
               }
            }
            // we make the assumption that every S with valence > 2 can be found
            // by counting the number of bonds.
            if ((a->GetImplicitHCount() - a->GetHvyDegree()) > 2) {
               a->SetPartialCharge(0.0); // category 7
               for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi)) {
                  OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
                  aa->SetPartialCharge(0.0);    // category 7
                  _lipoLabelNeighbors(aa, 0.6); // category 10
               }
            }
            break;
      }
      if (a->GetFormalCharge() != 0) {
         for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi)) {
            OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
            aa->SetPartialCharge(0.0);
            _lipoLabelNeighbors(aa, 0.0);
         }
      }
   }
   // Adjust combinations of scores (category 12 and 14)
   for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai)) {
      double value = a->GetPartialCharge();
      if ((value == 0.36 || value < 0.25) && value != 0.15) {
         a->SetPartialCharge(0.0); // category 12 and 14
      }
   }
}

double _lipoCalcAccSurf(OpenBabel::OBAtom* a) {
   //OpenBabel::OBElementTable et;
   double radius = OpenBabel::OBElements::GetVdwRad(a->GetAtomicNum());

   // Create sphere with uniformly distributed points
   std::vector<Coordinate> sphere;
   std::vector<Coordinate>::iterator itS;

   const double arclength(1.0 / sqrt(sqrt(3.0) * DENSITY));
   double dphi(arclength/radius);
	int nlayer(ROUND(PI/dphi) + 1);

   double phi(0.0);
   for (int i = 0; i < nlayer; ++i) {
      double rsinphi(radius*sin(phi));
      double z(radius*cos(phi));
      double dtheta((rsinphi==0) ? PI*2 : arclength/rsinphi);
      int tmpNbrPoints(ROUND(PI*2/dtheta));
      if (tmpNbrPoints <= 0) {
         tmpNbrPoints = 1;
      }
      dtheta = PI*2 / tmpNbrPoints;
      double theta((i%2) ? 0 : PI);
      for (int j = 0; j < tmpNbrPoints ; ++j) {
         Coordinate coord;
         coord.x = rsinphi*cos(theta) + a->x();
         coord.y = rsinphi*sin(theta) + a->y();
         coord.z = z + a->z();
         sphere.push_back(coord);
         theta += dtheta;
         if (theta > PI*2) {
            theta -= PI*2;
         }
      }
      phi += dphi;
   }

   // Define neighbors of atom
   std::list<OpenBabel::OBAtom*> aList(_lipoGetNeighbors(a));
   std::list<OpenBabel::OBAtom*>::iterator itA;

  // Check for every sphere-point if it is accessible
  int nbrAccSurfPoints(0);
  double delta(PROBE_RADIUS/radius);
  for (itS = sphere.begin(); itS != sphere.end(); ++itS) {
      Coordinate p;
      p.x = ((itS->x - a->x()) * delta) + itS->x;
      p.y = ((itS->y - a->y()) * delta) + itS->y;
      p.z = ((itS->z - a->z()) * delta) + itS->z;

      bool isAccessible(true);
      for (itA = aList.begin(); itA != aList.end(); ++itA) {
         OpenBabel::OBAtom* n(*itA);
         double distSq(((p.x - n->x()) * (p.x - n->x())) +
                       ((p.y - n->y()) * (p.y - n->y())) +
                       ((p.z - n->z()) * (p.z - n->z())));
         double sumSq((PROBE_RADIUS + OpenBabel::OBElements::GetVdwRad(n->GetAtomicNum())) *
                      (PROBE_RADIUS + OpenBabel::OBElements::GetVdwRad(n->GetAtomicNum())));

         if (distSq <= sumSq) {
            isAccessible = false;
            break;
         }
      }
      if (isAccessible) ++nbrAccSurfPoints;
   }
   double f(nbrAccSurfPoints/(double)sphere.size());
   return f * 4 * PI * radius*radius;
}

void _lipoGroupAtoms(OpenBabel::OBMol* m, Pharmacophore* pharmacophore) {
   std::set<OpenBabel::OBAtom*> atomSet; // keeps remaining atoms for step (3)
   std::set<OpenBabel::OBAtom*>::iterator itS;

   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai)) {
      atomSet.insert(a);
   }

   // Group rings smaller than 7
   std::vector<OpenBabel::OBRing*> allrings = m->GetSSSR();
   std::vector<OpenBabel::OBRing*>::iterator ri;
   OpenBabel::OBRing* ring;
   for (ri = allrings.begin(); ri != allrings.end(); ++ri) {
      ring = *ri;
      if (ring->Size() > 7) {
         continue;
      }

      double lipoSum(0.0);
      Coordinate center;
      for (OpenBabel::OBAtom* atom = m->BeginAtom(ai); atom; atom = m->NextAtom(ai)) {
         if (ring->IsMember(atom)) {
            atomSet.erase(atom);
            double lipo(atom->GetPartialCharge());
            lipoSum += lipo;
            center.x += lipo * atom->x();
            center.y += lipo * atom->y();
            center.z += lipo * atom->z();
         }
         else {
            continue;
         }
      }

      if (lipoSum > REF_LIPO) {
         PharmacophorePoint p;
         p.func = LIPO;
         p.hasNormal = false;
         p.normal.x = 0.0;
         p.normal.y = 0.0;
         p.normal.z = 0.0;
         p.alpha = funcSigma[LIPO];
         center.x /= lipoSum;
         center.y /= lipoSum;
         center.z /= lipoSum;
         p.point = center;
         pharmacophore->push_back(p);
      }
   }

   // Group atoms with three or more bonds
   for (itS = atomSet.begin(); itS != atomSet.end(); ++itS) {
      if (((*itS)->GetImplicitHCount() - (*itS)->GetHvyDegree()) > 2) {
         std::list<OpenBabel::OBAtom*> aList;
         aList.push_back(*itS);
         double lipoSum((*itS)->GetPartialCharge());
         Coordinate center;
         center.x += lipoSum * (*itS)->x();
         center.y += lipoSum * (*itS)->y();
         center.z += lipoSum * (*itS)->z();

         std::vector<OpenBabel::OBBond*>::iterator bi;
         for (OpenBabel::OBBond* b = (*itS)->BeginBond(bi); b; b = (*itS)->NextBond(bi)) {
            OpenBabel::OBAtom* a = b->GetNbrAtom(*itS);
            if (((a->GetImplicitHCount() - a->GetHvyDegree()) == 1) && (a->GetAtomicNum() != 1)) {
               double lipo(a->GetPartialCharge());
               lipoSum += lipo;
               aList.push_back(a);
               center.x += (lipo * a->x());
               center.y += (lipo * a->y());
               center.z += (lipo * a->z());
            }
         }

         if (lipoSum > REF_LIPO) {
            PharmacophorePoint p;
            p.func = LIPO;
            p.hasNormal = false;
            p.normal.x = 0.0;
            p.normal.y = 0.0;
            p.normal.z = 0.0;
				p.alpha = funcSigma[LIPO];
            center.x /= lipoSum;
            center.y /= lipoSum;
            center.z /= lipoSum;
            p.point = center;
            pharmacophore->push_back(p);
         }
      }
   }
}

std::list<OpenBabel::OBAtom*> _lipoGetNeighbors(OpenBabel::OBAtom* a) {
   ///OpenBabel::OBElementTable et;
   double radius = OpenBabel::OBElements::GetVdwRad(a->GetAtomicNum());
   std::list<OpenBabel::OBAtom*> aList;

   OpenBabel::OBMol* parent(a->GetParent());
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* aa = parent->BeginAtom(ai); aa; aa = parent->NextAtom(ai)) {
      if ((aa->GetAtomicNum() == 1) || (aa == a)) {
         continue;
      }

      double delta(radius + OpenBabel::OBElements::GetVdwRad(aa->GetAtomicNum()) + 2 * PROBE_RADIUS);
      double maxDistSq(delta*delta);
      double distSq((a->x() - aa->x()) * (a->x() - aa->x()) +
                    (a->y() - aa->y()) * (a->y() - aa->y()) +
                    (a->z() - aa->z()) * (a->z() - aa->z()));

      if (distSq <= maxDistSq) {
         aList.push_back(aa);
      }
   }
   return aList;
}

void _lipoLabelNeighbors(OpenBabel::OBAtom* a, double value) {
   std::vector<OpenBabel::OBBond*>::iterator bi;
   for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi)) {
      OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
      aa->SetPartialCharge(value * aa->GetPartialCharge());
   }
}

#else
#include <utilities.h>

void lipoFuncCalc(RDKit::ROMol *mol, Pharmacophore *pharmacophore) {
    // Find for each atom the 'topology-dependent term': t
    _lipoLabelAtoms(mol);
    // Multiply 't' with the accessible surface 's'
    const auto &conf = mol->getConformer();
    for (auto atom: mol->atoms()) {
        if (atom->getAtomicNum() == 1) continue;
        double t = atom->getProp<double>("_LipoContrib");
        if (t != 0.0) {
            atom->setProp("_LipoContrib", _lipoCalcAccSurf(atom, conf) * t);
        }
    }
    // Finally, calculate the lipophilic points
    _lipoGroupAtoms(mol, pharmacophore);
}

void _lipoLabelAtoms(RDKit::ROMol *mol) {
    // Give all atoms default score of 1
    for (auto atom: mol->atoms()) {
        atom->setProp("_LipoContrib", 1.0);
    }
    // Decrease scores based on O,N,S and charged atoms
    for (auto atom: mol->atoms()) {
        switch (atom->getAtomicNum()) {
            case 1:
                atom->setProp("_LipoContrib", 0.0);
                break;
            case 7:
                atom->setProp("_LipoContrib", 0.0);
                if (atom->getIsAromatic()) break;
                _lipoLabelNeighbors(atom, 0.25);
                if ((atom->getNumImplicitHs() - getHeavyDegree(atom)) != 0) {
                    for (const auto &nbri: boost::make_iterator_range(
                            atom->getOwningMol().getAtomNeighbors(atom))) {
                        const auto aa = atom->getOwningMol()[nbri];
                        aa->setProp("_LipoContrib", 0.0);
                        _lipoLabelNeighbors(aa, 0.0);
                    }
                }
                break;
            case 8:
                atom->setProp("_LipoContrib", 0.0);
                if (atom->getIsAromatic()) break;
                _lipoLabelNeighbors(atom, 0.25);
                for (const auto &nbri: boost::make_iterator_range(
                        atom->getOwningMol().getAtomBonds(atom))) {
                    const auto b = (*mol)[nbri];
                    if (b->getOtherAtom(atom)->getAtomicNum() == 1) {
                        for (const auto &nbrj: boost::make_iterator_range(
                                atom->getOwningMol().getAtomBonds(atom))) {
                            const auto b2 = (*mol)[nbrj];
                            const auto aa = b2->getOtherAtom(atom);
                            aa->setProp("_LipoContrib", 0.0);
                            _lipoLabelNeighbors(aa, 0.0);
                        }
                    }
                    if (b->getBondTypeAsDouble() == 2.0) {
                        const auto aa = b->getOtherAtom(atom);
                        aa->setProp("_LipoContrib", 0.0);
                        for (const auto &nbrj: boost::make_iterator_range(
                                atom->getOwningMol().getAtomBonds(atom))) {
                            const auto b2 = (*mol)[nbrj];
                            const auto aaa = b2->getOtherAtom(aa);
                            if (aaa == atom) continue;
                            aaa->setProp("_LipoContrib", 0.0);
                            _lipoLabelNeighbors(aaa, 0.6);
                        }
                    }
                }
                break;
            case 16:
                for (const auto &nbri: boost::make_iterator_range(
                        atom->getOwningMol().getAtomBonds(atom))) {
                    const auto b = atom->getOwningMol()[nbri];
                    if (b->getOtherAtom(atom)->getAtomicNum() == 1) {
                        atom->setProp("_LipoContrib", 0.0);
                        _lipoLabelNeighbors(atom, 0.0);
                    }
                    if (b->getBondTypeAsDouble() == 2.0) {
                        atom->setProp("_LipoContrib", 0.0);
                        _lipoLabelNeighbors(atom, 0.6);
                    }
                }
                // we make the assumption that every S with valence > 2 can be found
                // by counting the number of bonds.
                if ((atom->getNumImplicitHs() - getHeavyDegree(atom)) > 2) {
                    atom->setProp("_LipoContrib", 0.0);
                    for (const auto &nbri: boost::make_iterator_range(
                            atom->getOwningMol().getAtomNeighbors(atom))) {
                        const auto aa = atom->getOwningMol()[nbri];
                        aa->setProp("_LipoContrib", 0.0);
                        _lipoLabelNeighbors(aa, 0.6);
                    }
                }
                break;
        }
        if (atom->getFormalCharge() != 0) {
            for (const auto &nbri: boost::make_iterator_range(
                    atom->getOwningMol().getAtomNeighbors(atom))) {
                const auto aa = atom->getOwningMol()[nbri];
                aa->setProp("_LipoContrib", 0.0);
                _lipoLabelNeighbors(aa, 0.0);
            }
        }
    }
    // Adjust combinations of scores (category 12 and 14)
    for (auto atom: mol->atoms()) {
        double value = atom->getProp<double>("_LipoContrib");
        if ((value == 0.36 || value < 0.25) && value != 0.15) {
            atom->setProp("_LipoContrib", 0.0);
        }
    }
}

double _lipoCalcAccSurf(RDKit::Atom *atom, const RDKit::Conformer &conf) {
    double radius = RDKit::PeriodicTable::getTable()->getRvdw(atom->getAtomicNum());
    const auto &p = conf.getAtomPos(atom->getIdx());

    // Create sphere with uniformly distributed points
    std::vector<Coordinate> sphere;
    std::vector<Coordinate>::iterator itS;

    const double arclength(1.0 / sqrt(sqrt(3.0) * DENSITY));
    double dphi(arclength/radius);
    int nlayer(ROUND(PI/dphi) + 1);

    double phi(0.0);
    for (int i = 0; i < nlayer; ++i) {
        double rsinphi(radius*sin(phi));
        double z(radius*cos(phi));
        double dtheta((rsinphi==0) ? PI*2 : arclength/rsinphi);
        int tmpNbrPoints(ROUND(PI*2/dtheta));
        if (tmpNbrPoints <= 0) {
            tmpNbrPoints = 1;
        }
        dtheta = PI*2 / tmpNbrPoints;
        double theta((i%2) ? 0 : PI);
        for (int j = 0; j < tmpNbrPoints ; ++j) {
            Coordinate coord;
            coord.x = rsinphi*cos(theta) + p.x;
            coord.y = rsinphi*sin(theta) + p.y;
            coord.z = z + p.z;
            sphere.push_back(coord);
            theta += dtheta;
            if (theta > PI*2) {
                theta -= PI*2;
            }
        }
        phi += dphi;
    }
    // Define neighbors of atom
    std::list<RDKit::Atom *> aList(_lipoGetNeighbors(atom, conf));
    std::list<RDKit::Atom *>::iterator itA;

    // Check for every sphere-point if it is accessible
    int nbrAccSurfPoints(0);
    double delta(PROBE_RADIUS/radius);
    for (itS = sphere.begin(); itS != sphere.end(); ++itS) {
        Coordinate point;
        point.x = ((itS->x - p.x) * delta) + itS->x;
        point.y = ((itS->y - p.y) * delta) + itS->y;
        point.z = ((itS->z - p.z) * delta) + itS->z;

        bool isAccessible(true);
        for (itA = aList.begin(); itA != aList.end(); ++itA) {
            RDKit::Atom *n(*itA);
            const auto &pp = conf.getAtomPos(n->getIdx());
            double distSq(
                ((point.x - pp.x) * (point.x - pp.x)) +
                ((point.y - pp.y) * (point.y - pp.y)) +
                ((point.z - pp.z) * (point.z - pp.z))
            );
            const auto rr = RDKit::PeriodicTable::getTable()->getRvdw(n->getAtomicNum());
            double sumSq((PROBE_RADIUS + rr) * (PROBE_RADIUS + rr));
            if (distSq <= sumSq) {
                isAccessible = false;
                break;
            }
        }
        if (isAccessible) ++nbrAccSurfPoints;
    }
    double f(nbrAccSurfPoints/(double)sphere.size());
    return f * 4 * PI * radius*radius;
}

void _lipoGroupAtoms(RDKit::ROMol *m, Pharmacophore *pharmacophore) {
    const auto &conf = m->getConformer();
    if (!m->getRingInfo()->isInitialized()) {
        RDKit::MolOps::findSSSR(*m);
    }
    std::set<RDKit::Atom *> atomSet; // keeps remaining atoms for step (3)
    std::set<RDKit::Atom *>::iterator itS;
    for (auto atom: m->atoms()) {
        atomSet.insert(atom);
    }
    // Group rings smaller than 7
    const auto &ri = m->getRingInfo();
    for (auto ring: ri->atomRings()) {
        if (ring.size() > 7) continue;
        double lipoSum(0.0);
        Coordinate center;
        for (auto atom: m->atoms()) {
            if (std::find(ring.begin(), ring.end(), atom->getIdx()) != ring.end()) {
                atomSet.erase(atom);
                double lipo(atom->getProp<double>("_LipoContrib"));
                lipoSum += lipo;
                center.x += lipo * conf.getAtomPos(atom->getIdx()).x;
                center.y += lipo * conf.getAtomPos(atom->getIdx()).y;
                center.z += lipo * conf.getAtomPos(atom->getIdx()).z;
            } else {
                continue;
            }
        }
        if (lipoSum > REF_LIPO) {
            PharmacophorePoint p;
            p.func = LIPO;
            p.hasNormal = false;
            p.normal.x = 0.0;
            p.normal.y = 0.0;
            p.normal.z = 0.0;
            p.alpha = funcSigma[LIPO];
            center.x /= lipoSum;
            center.y /= lipoSum;
            center.z /= lipoSum;
            p.point = center;
            pharmacophore->push_back(p);
        }
    }
    // Group atoms with three or more bonds
    for (itS = atomSet.begin(); itS != atomSet.end(); ++itS) {
        if (((*itS)->getNumImplicitHs() - getHeavyDegree((*itS))) > 2) {
            std::list<RDKit::Atom *> aList;
            aList.push_back(*itS);
            double lipoSum((*itS)->getProp<double>("_LipoContrib"));
            Coordinate center;
            center.x += lipoSum * conf.getAtomPos((*itS)->getIdx()).x;
            center.y += lipoSum * conf.getAtomPos((*itS)->getIdx()).y;
            center.z += lipoSum * conf.getAtomPos((*itS)->getIdx()).z;
            for (const auto &nbri: boost::make_iterator_range(
                    (*itS)->getOwningMol().getAtomNeighbors((*itS)))) {
                const auto a = (*itS)->getOwningMol()[nbri];
                if (((a->getNumImplicitHs() - getHeavyDegree(a)) == 1) && (a->getAtomicNum() != 1)) {
                    double lipo(a->getProp<double>("_LipoContrib"));
                    lipoSum += lipo;
                    aList.push_back(a);
                    center.x += (lipo * conf.getAtomPos(a->getIdx()).x);
                    center.y += (lipo * conf.getAtomPos(a->getIdx()).y);
                    center.z += (lipo * conf.getAtomPos(a->getIdx()).z);
                }
            }
            if (lipoSum > REF_LIPO) {
                PharmacophorePoint p;
                p.func = LIPO;
                p.hasNormal = false;
                p.normal.x = 0.0;
                p.normal.y = 0.0;
                p.normal.z = 0.0;
                p.alpha = funcSigma[LIPO];
                center.x /= lipoSum;
                center.y /= lipoSum;
                center.z /= lipoSum;
                p.point = center;
                pharmacophore->push_back(p);
            }
        }
    }
}

std::list<RDKit::Atom *> _lipoGetNeighbors(RDKit::Atom *a, const RDKit::Conformer &conf) {
    double radius = RDKit::PeriodicTable::getTable()->getRvdw(a->getAtomicNum());
    const auto &p = conf.getAtomPos(a->getIdx());
    std::list<RDKit::Atom*> aList;
    for (auto aa: a->getOwningMol().atoms()) {
        if ((aa->getAtomicNum() == 1) || (aa == a)) continue;
        const auto &pp = conf.getAtomPos(aa->getIdx());
        double delta(radius + RDKit::PeriodicTable::getTable()->getRvdw(aa->getAtomicNum()) + 2 * PROBE_RADIUS);
        double maxDistSq(delta*delta);
        double distSq(
            (p.x - pp.x) * (p.x - pp.x) +
            (p.y - pp.y) * (p.y - pp.y) +
            (p.z - pp.z) * (p.z - pp.z)
        );
        if (distSq <= maxDistSq) {
            aList.push_back(aa);
        }
    }
    return aList;
}

void _lipoLabelNeighbors(RDKit::Atom *a, double value) {
    for (const auto &nbri: boost::make_iterator_range(a->getOwningMol().getAtomNeighbors(a))) {
        const auto aa = a->getOwningMol()[nbri];
        aa->setProp("_LipoContrib", value * aa->getProp<double>("_LipoContrib"));
    }
}

#endif
