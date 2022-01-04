/*******************************************************************************
hAccFuncCalc.cpp - Align-it

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

#include <hAccFuncCalc.h>

#ifndef USE_RDKIT
#include <openbabel/elements.h>
#include <openbabel/data.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>

void hAccFuncCalc(OpenBabel::OBMol* mol, Pharmacophore* pharmacophore) {
   // Create for every hydrogen acceptor a pharmacophore point
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* atom = mol->BeginAtom(ai); atom; atom = mol->NextAtom(ai)) {
      if (atom->GetAtomicNum() == 7 || atom->GetAtomicNum() == 8) {
         if (atom->GetFormalCharge() <= 0) {
            if(_hAccDelocalized(atom) || (_hAccCalcAccSurf(atom) < 0.02)) {
               continue;
            }
            PharmacophorePoint p;
            p.func = HACC;
            p.point.x = atom->x();
            p.point.y = atom->y();
            p.point.z = atom->z();
            p.hasNormal = true;
            p.alpha = funcSigma[HACC];
            p.normal = _hAccCalcNormal(atom);
            pharmacophore->push_back(p);
         }
      }
   }
}

double _hAccCalcAccSurf(OpenBabel::OBAtom* atom) {
   double radius(H_BOND_DIST);

   //---(1)-- create sphere with uniformly distributed points
   std::vector<Coordinate> sphere;
   std::vector<Coordinate>::iterator itS;

   const double arclength(1.0 / sqrt(sqrt(3.0) * DENSITY));
   double dphi(arclength/radius);
	int nlayer(ROUND(PI/dphi) + 1);

   double phi(0.0);
   for (int i(0); i < nlayer; ++i) {
      double rsinphi(radius*sin(phi));
      double z(radius*cos(phi));
      double dtheta((rsinphi==0) ? PI*2 : arclength/rsinphi);
      int tmpNbrPoints(ROUND(PI*2/dtheta));
      if(tmpNbrPoints <= 0) {
         tmpNbrPoints = 1;
      }
      dtheta = PI * 2.0 / tmpNbrPoints;
      double theta((i % 2) ? 0 : PI);
      for (int j(0) ; j < tmpNbrPoints ; ++j) {
         Coordinate coord;
         coord.x = rsinphi*cos(theta) + atom->x();
         coord.y = rsinphi*sin(theta) + atom->y();
         coord.z = z + atom->z();
         sphere.push_back(coord);
         theta += dtheta;
         if(theta > PI*2) {
            theta -= PI*2;
         }
      }
      phi += dphi;
   }

   //---(2)-- define neighbors of atom
   std::list<OpenBabel::OBAtom*> aList(_hAccGetNeighbors(atom));
   std::list<OpenBabel::OBAtom*>::iterator itA;

   //---(3) -- check for every sphere-point if it is accessible
   int nbrAccSurfPoints(0);
   double r;
   //OpenBabel::OBElementTable et;
   for (itS = sphere.begin(); itS != sphere.end(); ++itS) {
      bool isAccessible(true);
      for (itA = aList.begin(); itA != aList.end(); ++itA) {
         OpenBabel::OBAtom* n(*itA);
         double distSq(((itS->x - n->x()) * (itS->x - n->x())) +
                       ((itS->y - n->y()) * (itS->y - n->y())) +
                       ((itS->z - n->z()) * (itS->z - n->z())));
         r = OpenBabel::OBElements::GetVdwRad(n->GetAtomicNum());
         double sumSq((r + H_RADIUS) * (r + H_RADIUS));

         if (distSq <= sumSq) {
            isAccessible = false;
            break;
         }
      }

      if (isAccessible) {
         ++nbrAccSurfPoints;

      }
   }

   return (nbrAccSurfPoints/(double)sphere.size());
}

std::list<OpenBabel::OBAtom*> _hAccGetNeighbors(OpenBabel::OBAtom* a) {
   std::list<OpenBabel::OBAtom*> aList;
   OpenBabel::OBMol* parent(a->GetParent());

   double r;
   //OpenBabel::OBElementTable et;
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* aa = parent->BeginAtom(ai); aa; aa = parent->NextAtom(ai)) {
      if (*aa == a) {
         continue;
      }

      r = OpenBabel::OBElements::GetVdwRad(aa->GetAtomicNum());

      double delta(H_BOND_DIST + H_RADIUS + r);
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

bool _hAccDelocalized(OpenBabel::OBAtom* a) {
   if (a->GetAtomicNum() != 7) {
      return false;
   }
   //if (a->IsAromatic() && a->GetImplicitValence() == 3)
   if (a->IsAromatic() && a->GetTotalDegree() == 3) {
      return true;
   }

   std::vector<OpenBabel::OBBond*>::iterator bi1;
   for (OpenBabel::OBBond* b1 = a->BeginBond(bi1); b1; b1 = a->NextBond(bi1)) {
      OpenBabel::OBAtom* aa = b1->GetNbrAtom(a);

      //if (aa->IsAromatic() && a->GetImplicitValence() == 3)
      if (aa->IsAromatic() && a->GetTotalDegree() == 3) {
         return true;
      }

      if (aa->GetAtomicNum() == 6) {
         std::vector<OpenBabel::OBBond*>::iterator bi2;
         for (OpenBabel::OBBond* b2 = aa->BeginBond(bi2); b2; b2 = aa->NextBond(bi2)) {
            OpenBabel::OBAtom* aaa = b2->GetNbrAtom(aa);

            if (aaa == a) {
               continue;
            }
            if (b2->GetBondOrder() == 2) {
               if (aaa->GetAtomicNum() == 8)  return true;
               if (aaa->GetAtomicNum() == 7)  return true;
               if (aaa->GetAtomicNum() == 16) return true;
            }
         }
      }
      else if (aa->GetAtomicNum() == 16) {
         std::vector<OpenBabel::OBBond*>::iterator bi2;
         for (OpenBabel::OBBond* b2 = aa->BeginBond(bi2); b2; b2 = aa->NextBond(bi2)) {
            OpenBabel::OBAtom* aaa = b2->GetNbrAtom(aa);

            if (aaa == a) {
               continue;
            }
            if ((b2->GetBondOrder() == 2) && (aaa->GetAtomicNum() == 8)) {
               return true;
            }
         }
      }
   }
   return false;
}

Coordinate _hAccCalcNormal(OpenBabel::OBAtom* a) {
   Coordinate normal;
   std::vector<OpenBabel::OBBond*>::iterator bi;
   for (OpenBabel::OBBond* b = a->BeginBond(bi); b; b = a->NextBond(bi)) {
      OpenBabel::OBAtom* aa = b->GetNbrAtom(a);
      if (aa->GetAtomicNum() == 1) {
         continue;
      }
      normal.x += (aa->x() - a->x());
      normal.y += (aa->y() - a->y());
      normal.z += (aa->z() - a->z());
   }
   double length(sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z));
   normal.x /= length;
   normal.y /= length;
   normal.z /= length;

   normal.x = -normal.x;
   normal.y = -normal.y;
   normal.z = -normal.z;

   normal.x += a->x();
   normal.y += a->y();
   normal.z += a->z();

   return normal;
}

#else
#include <GraphMol/PeriodicTable.h>

void hAccFuncCalc(RDKit::ROMol *mol, Pharmacophore *pharmacophore) {
    // Create for every hydrogen acceptor a pharmacophore point
    const auto &conf = mol->getConformer();
    for (auto atom: mol->atoms()) {
        if (atom->getAtomicNum() == 7 || atom->getAtomicNum() == 8) {
            if (atom->getFormalCharge() <= 0) {
                if (_hAccDelocalized(atom) || _hAccCalcAccSurf(atom, conf) < 0.02) {
                    continue;
                }
                PharmacophorePoint p;
                p.func = HACC;
                p.point.x = conf.getAtomPos(atom->getIdx()).x;
                p.point.y = conf.getAtomPos(atom->getIdx()).y;
                p.point.z = conf.getAtomPos(atom->getIdx()).z;
                p.hasNormal = true;
                p.alpha = funcSigma[HACC];
                p.normal = _hAccCalcNormal(atom, conf);
                pharmacophore->push_back(p);
            }
        }
    }
}

double _hAccCalcAccSurf(RDKit::Atom *atom, const RDKit::Conformer &conf) {
    double radius(H_BOND_DIST);

    //---(1)-- create sphere with uniformly distributed points
    std::vector<Coordinate> sphere;
    std::vector<Coordinate>::iterator itS;

    const double arclength(1.0 / sqrt(sqrt(3.0) * DENSITY));
    double dphi(arclength/radius);
    int nlayer(ROUND(PI/dphi) + 1);

    double phi(0.0);
    for (int i(0); i < nlayer; ++i) {
        double rsinphi(radius*sin(phi));
        double z(radius*cos(phi));
        double dtheta((rsinphi==0) ? PI*2 : arclength/rsinphi);
        int tmpNbrPoints(ROUND(PI*2/dtheta));
        if(tmpNbrPoints <= 0) {
            tmpNbrPoints = 1;
        }
        dtheta = PI * 2.0 / tmpNbrPoints;
        double theta((i % 2) ? 0 : PI);
        for (int j(0) ; j < tmpNbrPoints ; ++j) {
            Coordinate coord;
            coord.x = rsinphi*cos(theta) + conf.getAtomPos(atom->getIdx()).x;
            coord.y = rsinphi*sin(theta) + conf.getAtomPos(atom->getIdx()).y;
            coord.z = z + conf.getAtomPos(atom->getIdx()).z;
            sphere.push_back(coord);
            theta += dtheta;
            if(theta > PI*2) {
                theta -= PI*2;
            }
        }
        phi += dphi;
    }
    //---(2)-- define neighbors of atom
    std::list<RDKit::Atom *> aList(_hAccGetNeighbors(atom, conf));
    std::list<RDKit::Atom *>::iterator itA;

    //---(3) -- check for every sphere-point if it is accessible
    int nbrAccSurfPoints(0);
    double r;
    for (itS = sphere.begin(); itS != sphere.end(); ++itS) {
        bool isAccessible(true);
        for (itA = aList.begin(); itA != aList.end(); ++itA) {
            RDKit::Atom* n(*itA);
            const auto &p = conf.getAtomPos(n->getIdx());
            double distSq(
                ((itS->x - p.x) * (itS->x - p.x)) +
                ((itS->y - p.y) * (itS->y - p.y)) +
                ((itS->z - p.z) * (itS->z - p.z))
            );
            r = RDKit::PeriodicTable::getTable()->getRvdw(n->getAtomicNum());
            double sumSq((r + H_RADIUS) * (r + H_RADIUS));
            if (distSq <= sumSq) {
                isAccessible = false;
                break;
            }
        }
        if (isAccessible) {
            ++nbrAccSurfPoints;
        }
    }
   return (nbrAccSurfPoints/(double)sphere.size());
}

std::list<RDKit::Atom *> _hAccGetNeighbors(RDKit::Atom *a, const RDKit::Conformer &conf) {
    std::list<RDKit::Atom *> aList;
    double r;
    for (auto atom: a->getOwningMol().atoms()) {
        if (atom == a) continue;
        const auto &p = conf.getAtomPos(a->getIdx());
        const auto &pp = conf.getAtomPos(atom->getIdx());
        r = RDKit::PeriodicTable::getTable()->getRvdw(a->getAtomicNum());
        double delta(H_BOND_DIST + H_RADIUS + r);
        double maxDistSq(delta*delta);
        double distSq(
             (p.x - pp.x) * (p.x - pp.x) +
             (p.y - pp.y) * (p.y - pp.y) +
             (p.z - pp.z) * (p.z - pp.z)
        );
        if (distSq <= maxDistSq) {
            aList.push_back(atom);
        }
    }
    return aList;
}


bool _hAccDelocalized(RDKit::Atom *a) {
    if (a->getAtomicNum() != 7) { return false; }
    if (a->getIsAromatic() && a->getTotalDegree() == 3) { return true; }
    for (const auto &nbri: boost::make_iterator_range(a->getOwningMol().getAtomNeighbors(a))) {
        const auto aa = a->getOwningMol()[nbri];
        if (aa->getIsAromatic() && a->getTotalDegree() == 3) {
            return true;
        }
        if (aa->getAtomicNum() == 6) {
            for (const auto &nbrj: boost::make_iterator_range(aa->getOwningMol().getAtomBonds(aa))) {
                const auto bnd = aa->getOwningMol()[nbrj];
                const auto aaa = bnd->getOtherAtom(aa);
                if (aaa == a) continue;
                if (bnd->getBondTypeAsDouble() == 2.0) {
                    if (aaa->getAtomicNum() == 8) return true;
                    if (aaa->getAtomicNum() == 7) return true;
                    if (aaa->getAtomicNum() == 16) return true;
                }
            }
        } else if (aa->getAtomicNum() == 16) {
            for (const auto &nbrj: boost::make_iterator_range(aa->getOwningMol().getAtomBonds(aa))) {
                const auto bnd = aa->getOwningMol()[nbrj];
                const auto aaa = bnd->getOtherAtom(aa);
                if (aaa == a) continue;
                if ((bnd->getBondTypeAsDouble() == 2.0) && (aaa->getAtomicNum() == 8)) {
                    return true;
                }
            }
        }
    }
    return false;
}

Coordinate _hAccCalcNormal(RDKit::Atom *a, const RDKit::Conformer &conf) {
    Coordinate normal;
    const auto &p = conf.getAtomPos(a->getIdx());
    for (const auto &nbri:  boost::make_iterator_range(a->getOwningMol().getAtomNeighbors(a))) {
        const auto aa = a->getOwningMol()[nbri];
        if (aa->getAtomicNum() == 1) continue;
        const auto &pp = conf.getAtomPos(aa->getIdx());
        normal.x += (pp.x - p.x);
        normal.y += (pp.y - p.y);
        normal.z += (pp.z - p.z);
    }
    double length(sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z));
    normal.x /= length;
    normal.y /= length;
    normal.z /= length;
    normal.x = -normal.x;
    normal.y = -normal.y;
    normal.z = -normal.z;
    normal.x += p.x;
    normal.y += p.y;
    normal.z += p.z;
    return normal;
}

#endif
