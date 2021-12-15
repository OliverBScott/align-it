/*******************************************************************************
utilities.cpp - Align-it

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

#include <utilities.h>

Coordinate translate(Coordinate &p, Coordinate &t) {
    Coordinate n;
    n.x = p.x + t.x;
    n.y = p.y + t.y;
    n.z = p.z + t.z;
    return n;
}

Coordinate rotate(Coordinate &p, SiMath::Matrix &U) {
    Coordinate n;
    n.x = p.x * U[0][0] + p.y * U[0][1] + p.z *  U[0][2];
    n.y = p.x * U[1][0] + p.y * U[1][1] + p.z *  U[1][2];
    n.z = p.x * U[2][0] + p.y * U[2][1] + p.z *  U[2][2];
    return n;
}

void normalise(Coordinate &p) {
    double d = p.x*p.x + p.y*p.y + p.z*p.z;
    d = sqrt(d);
    p.x /= d;
    p.y /= d;
    p.z /= d;
    return;
}

double norm(Coordinate &p) {
    return sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
}

double dotProduct(Coordinate &p1, Coordinate &p2) {
    return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
}

Coordinate crossProduct(Coordinate &p1, Coordinate &p2) {
    Coordinate p;
    p.x = (p1.y * p2.z) - (p1.z * p2.y);
    p.y = (p1.z * p2.x) - (p1.x * p2.z);
    p.z = (p1.x * p2.y) - (p1.y * p2.x);
    return p;
}

double cosine(Coordinate &p1, Coordinate &p2) {
    double c(p1.x*p2.x + p1.y*p2.y + p1.z*p2.z);
    c /= norm(p1);
    c /= norm(p2);
    return c;
}

double distance(Coordinate &p1, Coordinate &p2) {
    double d(0.0);
    d += (p1.x - p2.x)*(p1.x - p2.x);
    d += (p1.y - p2.y)*(p1.y - p2.y);
    d += (p1.z - p2.z)*(p1.z - p2.z);
    return sqrt(d);
}

void normalise(SiMath::Vector &v) {
    double d(0.0);
    for (unsigned int i = 0; i < 4; i++) d += v[i] * v[i];
    v /= sqrt(d);
    return;
}

SiMath::Matrix quat2Rotation(SiMath::Vector &Q) {
    double d1sq(Q[1]*Q[1]);
    double d2sq(Q[2]*Q[2]);
    double d3sq(Q[3]*Q[3]);

    SiMath::Matrix U(3,3);
    U[0][0] = 1.0 - 2.0 * d2sq - 2.0 * d3sq;
    U[1][0] = 2.0 * (Q[1] * Q[2] + Q[0] * Q[3]);
    U[2][0] = 2.0 * (Q[1] * Q[3] - Q[0] * Q[2]);

    U[0][1] = 2.0 * (Q[2] * Q[1] - Q[0] * Q[3]);
    U[1][1] = 1.0 - 2.0 * d1sq - 2.0 * d3sq;
    U[2][1] = 2.0 * (Q[2] * Q[3] + Q[0] * Q[1]);

    U[0][2] = 2.0 * (Q[3] * Q[1] + Q[0] * Q[2]);
    U[1][2] = 2.0 * (Q[3] * Q[2] - Q[0] * Q[1]);
    U[2][2] = 1.0 - 2.0 * d1sq - 2.0 * d2sq;

    return U;
}

void inverseHessian(SiMath::Matrix &H) {
    // define H = [H00 H01 | H10 H11]
    // R0 = inv(H00)
    SiMath::Matrix R0(2,2);
    double d = (H[0][0] * H[1][1] - H[0][1] * H[1][0]);
    if ( d > 1e-6 || d < -1e-6) d = 1.0 / d;

    R0[0][0] = d * H[1][1];
    R0[1][1] = d * H[0][0];
    R0[0][1] = -d * H[0][1];
    R0[1][0] = -d * H[1][0];

    // R1 = H10 * R0
    SiMath::Matrix R1(2,2);
    R1[0][0] = H[2][0]*R0[0][0] + H[2][1] * R0[1][0];
    R1[0][1] = H[2][0]*R0[0][1] + H[2][1] * R0[1][1];
    R1[1][0] = H[3][0]*R0[0][0] + H[3][1] * R0[1][0];
    R1[1][1] = H[3][0]*R0[0][1] + H[3][1] * R0[1][1];

    // R2 = R0 * H01
    SiMath::Matrix R2(2,2);
    R2[0][0] = R0[0][0]*H[0][2] + R0[0][1] * H[1][2];
    R2[0][1] = R0[0][0]*H[0][3] + R0[0][1] * H[1][3];
    R2[1][0] = R0[1][0]*H[0][2] + R0[1][1] * H[1][2];
    R2[1][1] = R0[1][0]*H[0][3] + R0[1][1] * H[1][3];

    // R3 = H10 * R1
    SiMath::Matrix R3(2,2);
    R3[0][0] = H[2][0]*R2[0][0] + H[2][1] * R2[1][0];
    R3[0][1] = H[2][0]*R2[0][1] + H[2][1] * R2[1][1];
    R3[1][0] = H[3][0]*R2[0][0] + H[3][1] * R2[1][0];
    R3[1][1] = H[3][0]*R2[0][1] + H[3][1] * R2[1][1];

    // R3 = R3 - A11
    R3[0][0] -= H[2][2];
    R3[0][1] -= H[2][3];
    R3[1][0] -= H[3][2];
    R3[1][1] -= H[3][3];

    // R3 = inv(R3)
    d = (R3[0][0] * R3[1][1] - R3[0][1] * R3[1][0]);
    if ( d > 1e-6 || d < -1e-6 ) R3 /= d;

    // swap [0][0] with [1][1]
    d = R3[1][1];
    R3[1][1] = R3[0][0];
    R3[0][0] = d;

    // negate [0][1] and [1][0]
    R3[1][0] = -R3[1][0];
    R3[0][1] = -R3[0][1];

    // H01 = R2 * R3
    H[0][2] = R2[0][0] * R3[0][0] + R2[0][1] * R3[1][0];
    H[0][3] = R2[0][0] * R3[0][1] + R2[0][1] * R3[1][1];
    H[1][2] = R2[1][0] * R3[0][0] + R2[1][1] * R3[1][0];
    H[1][3] = R2[1][0] * R3[0][1] + R2[1][1] * R3[1][1];

    // H10 = R3 * R1
    H[2][0] = R3[0][0] * R1[0][0] + R3[0][1] * R1[1][0];
    H[2][1] = R3[0][0] * R1[0][1] + R3[0][1] * R1[1][1];
    H[3][0] = R3[1][0] * R1[0][0] + R3[1][1] * R1[1][0];
    H[3][1] = R3[1][0] * R1[0][1] + R3[1][1] * R1[1][1];

    // R4 = R2 * H10
    SiMath::Matrix R4(2,2);
    R4[0][0] = R2[0][0] * H[2][0] + R2[0][1] * H[3][0];
    R4[0][1] = R2[0][0] * H[2][1] + R2[0][1] * H[3][1];
    R4[1][0] = R2[1][0] * H[2][0] + R2[1][1] * H[3][0];
    R4[1][1] = R2[1][0] * H[2][1] + R2[1][1] * H[3][1];

    // H00 = R0 - R4
    H[0][0] = R0[0][0] - R4[0][0];
    H[0][1] = R0[0][1] - R4[0][1];
    H[1][0] = R0[1][0] - R4[1][0];
    H[1][1] = R0[1][1] - R4[1][1];

    // H11 = -R3
    H[2][2] = -R3[0][0];
    H[2][3] = -R3[0][1];
    H[3][2] = -R3[1][0];
    H[3][3] = -R3[1][1];

    return;
}

double VolumeOverlap(PharmacophorePoint &p1, PharmacophorePoint &p2, bool n) {
    double r2 = (p1.point.x - p2.point.x) * (p1.point.x - p2.point.x);
    r2 += (p1.point.y - p2.point.y) * (p1.point.y - p2.point.y);
    r2 += (p1.point.z - p2.point.z) * (p1.point.z - p2.point.z);
    double vol(1.0);
    if (n) {
        if(((p1.func == AROM) || (p1.func == HYBL))
           && ((p2.func == AROM) || (p2.func == HYBL))
           && ( p1.hasNormal )
           && ( p2.hasNormal )) {
            vol = fabs(cosine(p1.normal, p2.normal));
        }
        else if (((p1.func == HACC) || (p1.func == HDON) || (p1.func == HYBH))
                && ((p2.func == HACC) || (p2.func == HDON) || (p2.func == HYBH))
                && ( p1.hasNormal )
                && ( p2.hasNormal )) {
            vol = cosine(p1.normal, p2.normal);
        }
    }
    vol *= GCI2 * pow(PI/(p1.alpha + p2.alpha), 1.5);
    vol *= exp(-(p1.alpha * p2.alpha) * r2/(p1.alpha + p2.alpha));
    return vol;
}

double VolumeOverlap(PharmacophorePoint *p1, PharmacophorePoint *p2, bool n) {
    double r2 = (p1->point.x - p2->point.x) * (p1->point.x - p2->point.x);
    r2 += (p1->point.y - p2->point.y) * (p1->point.y - p2->point.y);
    r2 += (p1->point.z - p2->point.z) * (p1->point.z - p2->point.z);
    double vol(1.0);
    if (n) {
        if(((p1->func == AROM) || (p1->func == HYBL))
           && ((p2->func == AROM) || (p2->func == HYBL))
           && ( p1->hasNormal )
           && ( p2->hasNormal )) {
            vol = fabs(cosine(p1->normal, p2->normal));
        }
        else if (((p1->func == HACC) || (p1->func == HDON) || (p1->func == HYBH))
                && ((p2->func == HACC) || (p2->func == HDON) || (p2->func == HYBH))
                && ( p1->hasNormal )
                && ( p2->hasNormal )) {
            vol = cosine(p1->normal, p2->normal);
        }
    }
    vol *= GCI2 * pow(PI/(p1->alpha + p2->alpha), 1.5);
    vol *= exp(-(p1->alpha * p2->alpha) * r2/(p1->alpha + p2->alpha));
    return vol;
}

void positionPharmacophore(Pharmacophore &pharm, SiMath::Matrix &U, SolutionInfo &s) {
    // transpose of rotation matrix
    SiMath::Matrix rt = s.rotation2.transpose();
    for (int i(0); i < pharm.size(); ++i) {
        // translate normal origin
        pharm[i].normal.x -= pharm[i].point.x;
        pharm[i].normal.y -= pharm[i].point.y;
        pharm[i].normal.z -= pharm[i].point.z;

        // translate pharmacophore to center of db pharm
        pharm[i].point.x -= s.center2.x;
        pharm[i].point.y -= s.center2.y;
        pharm[i].point.z -= s.center2.z;

        // align with main axes
        pharm[i].point = rotate(pharm[i].point, rt);
        pharm[i].normal = rotate(pharm[i].normal, rt);

        // rotate according to best rotor
        pharm[i].point = rotate(pharm[i].point, U);
        pharm[i].normal = rotate(pharm[i].normal, U);

        // rotate back to main axes of the reference
        pharm[i].point = rotate(pharm[i].point, s.rotation1);
        pharm[i].normal = rotate(pharm[i].normal, s.rotation1);

        // move to center of reference
        pharm[i].point.x += s.center1.x;
        pharm[i].point.y += s.center1.y;
        pharm[i].point.z += s.center1.z;

        // translate normal back from origin
        pharm[i].normal.x += pharm[i].point.x;
        pharm[i].normal.y += pharm[i].point.y;
        pharm[i].normal.z += pharm[i].point.z;
    }
    return;
}

void TransformPharmacophore(Pharmacophore &pharm, SiMath::Matrix &U, Coordinate &center1, Coordinate &center2) {
    for (int i(0); i < pharm.size(); ++i) {
        PharmacophorePoint pp(pharm[i]);

        // translate and rotate normal[0]
        pharm[i].normal.x -= pharm[i].point.x;
        pharm[i].normal.y -= pharm[i].point.y;
        pharm[i].normal.z -= pharm[i].point.z;

        // translate and rotate pharmacophore center
        pharm[i].point.x -= center2.x;
        pharm[i].point.y -= center2.y;
        pharm[i].point.z -= center2.z;
        pharm[i].point = rotate(pharm[i].point, U);
        pharm[i].point.x += center1.x;
        pharm[i].point.y += center1.y;
        pharm[i].point.z += center1.z;

        pharm[i].normal = rotate(pharm[i].normal, U);
        pharm[i].normal.x += pharm[i].point.x;
        pharm[i].normal.y += pharm[i].point.y;
        pharm[i].normal.z += pharm[i].point.z;
    }
    return;
}

#ifndef USE_RDKIT
#include <openbabel/atom.h>
#include <openbabel/mol.h>

void positionMolecule(OpenBabel::OBMol *m, SiMath::Matrix &U, SolutionInfo &s) {
    // transpose of rotation matrix
    SiMath::Matrix rt = s.rotation2.transpose();
    Coordinate point;
    std::vector<OpenBabel::OBAtom*>::iterator ai;
    for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai)) {
        point.x = a->x();
        point.y = a->y();
        point.z = a->z();
        point.x -= s.center2.x;
        point.y -= s.center2.y;
        point.z -= s.center2.z;
        point = rotate(point, rt);
        point = rotate(point, U);
        point = rotate(point, s.rotation1);
        point.x += s.center1.x;
        point.y += s.center1.y;
        point.z += s.center1.z;
        a->SetVector(point.x, point.y, point.z);
    }
    return;
}

void TransformMolecule(OpenBabel::OBMol *m, SiMath::Matrix &U, Coordinate &center1, Coordinate &center2) {
   Coordinate point;
   std::vector<OpenBabel::OBAtom*>::iterator ai;
   for (OpenBabel::OBAtom* a = m->BeginAtom(ai); a; a = m->NextAtom(ai)) {
       point.x = a->x();
       point.y = a->y();
       point.z = a->z();

       point.x -= center2.x;
       point.y -= center2.y;
       point.z -= center2.z;
       point = rotate(point, U);
       point.x += center1.x;
       point.y += center1.y;
       point.z += center1.z;
       a->SetVector(point.x, point.y, point.z);
   }
   return;
}

#else
#include <Geometry/point.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/RDKitBase.h>

void positionMolecule(RDKit::ROMol *m, SiMath::Matrix &U, SolutionInfo &s) {
    RDKit::Conformer &conf = m->getConformer();
    SiMath::Matrix rt = s.rotation2.transpose();
    Coordinate point;
    for (unsigned int i =0; i < m->getNumAtoms(); ++i) {
        RDGeom::Point3D p = conf.getAtomPos(i);
        point.x = p.x;
        point.y = p.y;
        point.z = p.z;
        point.x -= s.center2.x;
        point.y -= s.center2.y;
        point.z -= s.center2.z;
        point = rotate(point, rt);
        point = rotate(point, U);
        point = rotate(point, s.rotation1);
        point.x += s.center1.x;
        point.y += s.center1.y;
        point.z += s.center1.z;
        conf.setAtomPos(i, RDGeom::Point3D(
            point.x, point.y, point.z));
    }
}

void TransformMolecule(RDKit::ROMol *m, SiMath::Matrix &U, Coordinate &center1, Coordinate &center2) {
    RDKit::Conformer &conf = m->getConformer();
    Coordinate point;
    for (unsigned int i =0; i < m->getNumAtoms(); ++i) {
        RDGeom::Point3D p = conf.getAtomPos(i);
        point.x = p.x;
        point.y = p.y;
        point.z = p.z;
        point.x -= center2.x;
        point.y -= center2.y;
        point.z -= center2.z;
        point = rotate(point, U);
        point.x += center1.x;
        point.y += center1.y;
        point.z += center1.z;
        conf.setAtomPos(i, RDGeom::Point3D(
            point.x, point.y, point.z));
    }
}

unsigned int getHeavyDegree(RDKit::Atom *a) {
    unsigned int deg(0);
    for (const auto &nbri: boost::make_iterator_range(
            a->getOwningMol().getAtomNeighbors(a))) {
        const auto aa = a->getOwningMol()[nbri];
        if (aa->getAtomicNum() != 1) {
            ++deg;
        }
    }
    return deg;
}

#endif
