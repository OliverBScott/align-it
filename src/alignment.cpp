/*******************************************************************************
alignment.cpp - Align-it

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

#include <alignment.h>

Alignment::Alignment(PharmacophoreMap &fMap) :
    _refMap(),
    _dbMap(),
    _refCenter(),
    _dbCenter(),
    _refRotMat(3,3,0.0),
    _dbRotMat(3,3,0.0),
    _dCdq(4,0.0),
    _d2Cdq2(4,4,0.0),
    _grad(4,0.0),
    _AkA(fMap.size()),
    _nbrPoints(0),
    _nbrExcl(0)
{
    // compute centers of the two pharmacophore sets
    PharmacophoreMap::iterator mi;
    unsigned int nbrMatch(0);

    // compute the centroid of the pharmacophore sets
    double V1(0.0), v1(0.0), V2(0.0), v2(0.0);
    for (mi = fMap.begin(); mi != fMap.end(); ++mi) {
        if (mi->first->func == EXCL) {
            _nbrExcl++;
            continue; // exclusion spheres not considered
        }
        _nbrPoints++;

        v1 = GCI * pow(PI/mi->first->alpha,1.5);
        V1 += v1;
        _refCenter.x += v1 * mi->first->point.x;
        _refCenter.y += v1 * mi->first->point.y;
        _refCenter.z += v1 * mi->first->point.z;

        v2 = GCI * pow(PI/mi->second->alpha,1.5);
        V2 += v2;
        _dbCenter.x += v2 * mi->second->point.x;
        _dbCenter.y += v2 * mi->second->point.y;
        _dbCenter.z += v2 * mi->second->point.z;
        ++nbrMatch;
    }
    _refCenter.x /= V1;
    _refCenter.y /= V1;
    _refCenter.z /= V1;
    _dbCenter.x /= V2;
    _dbCenter.y /= V2;
    _dbCenter.z /= V2;

    // translate the pharmacophores to the centers
    // and compute center of mass matrix
    SiMath::Matrix mass1(3,3,0.0);
    SiMath::Matrix mass2(3,3,0.0);
    for (mi = fMap.begin(); mi != fMap.end(); ++mi) {
        PharmacophorePoint p1, p2;

        p1.point.x = mi->first->point.x - _refCenter.x;
        p1.point.y = mi->first->point.y - _refCenter.y;
        p1.point.z = mi->first->point.z - _refCenter.z;
        p1.func = mi->first->func;
        p1.alpha = mi->first->alpha;
        p1.normal.x = mi->first->normal.x - mi->first->point.x;
        p1.normal.y = mi->first->normal.y - mi->first->point.y;
        p1.normal.z = mi->first->normal.z - mi->first->point.z;

        p2.point.x = mi->second->point.x - _dbCenter.x;
        p2.point.y = mi->second->point.y - _dbCenter.y;
        p2.point.z = mi->second->point.z - _dbCenter.z;
        p2.func = mi->second->func;
        p2.alpha = mi->second->alpha;
        p2.normal.x = mi->second->normal.x - mi->second->point.x;
        p2.normal.y = mi->second->normal.y - mi->second->point.y;
        p2.normal.z = mi->second->normal.z - mi->second->point.z;

        if (mi->first->func != EXCL) {

            v1 = GCI * pow(PI/mi->first->alpha,1.5);
            mass1[0][0] += v1 * p1.point.x * p1.point.x;
            mass1[0][1] += v1 * p1.point.x * p1.point.y;
            mass1[0][2] += v1 * p1.point.x * p1.point.z;
            mass1[1][1] += v1 * p1.point.y * p1.point.y;
            mass1[1][2] += v1 * p1.point.y * p1.point.z;
            mass1[2][2] += v1 * p1.point.z * p1.point.z;

            v2 = GCI * pow(PI/mi->second->alpha,1.5);
            mass2[0][0] += v2 * p2.point.x * p2.point.x;
            mass2[0][1] += v2 * p2.point.x * p2.point.y;
            mass2[0][2] += v2 * p2.point.x * p2.point.z;
            mass2[1][1] += v2 * p2.point.y * p2.point.y;
            mass2[1][2] += v2 * p2.point.y * p2.point.z;
            mass2[2][2] += v2 * p2.point.z * p2.point.z;
        }
        // add new points to local maps
        _refMap.push_back(p1);
        _dbMap.push_back(p2);
    }
    // use SVD to compute the best set of rotations
    // set lower triangle
    mass1[1][0] = mass1[0][1];
    mass1[2][0] = mass1[0][2];
    mass1[2][1] = mass1[1][2];

    // normalize mass matrix
    mass1 /= V1;

    // compute SVD of the mass matrix
    SiMath::SVD svd(mass1, true, true);
    _refRotMat = svd.getU();

    // check if determinant is 1, otherwise it is a mirroring instead of rotation
    double det = _refRotMat[0][0]*_refRotMat[1][1]*_refRotMat[2][2]
            + _refRotMat[2][1]*_refRotMat[1][0]*_refRotMat[0][2]
            + _refRotMat[0][1]*_refRotMat[1][2]*_refRotMat[2][0]
            - _refRotMat[0][0]*_refRotMat[1][2]*_refRotMat[2][1]
            - _refRotMat[1][1]*_refRotMat[2][0]*_refRotMat[0][2]
            - _refRotMat[2][2]*_refRotMat[0][1]*_refRotMat[1][0];

    // check if it is a rotation matrix and not a mirroring
    if (det < 0) {
        // switch sign of third column
        _refRotMat[0][2] = -_refRotMat[0][2];
        _refRotMat[1][2] = -_refRotMat[1][2];
        _refRotMat[2][2] = -_refRotMat[2][2];
    }
    // set lower triangle
    mass2[1][0] = mass2[0][1];
    mass2[2][0] = mass2[0][2];
    mass2[2][1] = mass2[1][2];

    // normalize mass matrix
    mass2 /= V2;

    // compute SVD of the mass matrix
    SiMath::SVD svd2(mass2, true, true);
    _dbRotMat = svd2.getU();

    // check if determinant is 1, otherwise it is a mirroring instead of rotation
    det = _dbRotMat[0][0]*_dbRotMat[1][1]*_dbRotMat[2][2]
            + _dbRotMat[2][1]*_dbRotMat[1][0]*_dbRotMat[0][2]
            + _dbRotMat[0][1]*_dbRotMat[1][2]*_dbRotMat[2][0]
            - _dbRotMat[0][0]*_dbRotMat[1][2]*_dbRotMat[2][1]
            - _dbRotMat[1][1]*_dbRotMat[2][0]*_dbRotMat[0][2]
            - _dbRotMat[2][2]*_dbRotMat[0][1]*_dbRotMat[1][0];

    // check if it is a rotation matrix and not a mirroring
    if (det < 0) {
        // switch sign of third column
        _dbRotMat[0][2] = -_dbRotMat[0][2];
        _dbRotMat[1][2] = -_dbRotMat[1][2];
        _dbRotMat[2][2] = -_dbRotMat[2][2];
    }
    // rotate points towards main axes
    for (unsigned int i(0); i < _refMap.size(); ++i) {
        // Rotate points
        double x = _refMap[i].point.x;
        double y = _refMap[i].point.y;
        double z = _refMap[i].point.z;
        _refMap[i].point.x = _refRotMat[0][0]*x + _refRotMat[1][0]*y + _refRotMat[2][0]*z;
        _refMap[i].point.y = _refRotMat[0][1]*x + _refRotMat[1][1]*y + _refRotMat[2][1]*z;
        _refMap[i].point.z = _refRotMat[0][2]*x + _refRotMat[1][2]*y + _refRotMat[2][2]*z;

        x = _dbMap[i].point.x;
        y = _dbMap[i].point.y;
        z = _dbMap[i].point.z;
        _dbMap[i].point.x = _dbRotMat[0][0]*x + _dbRotMat[1][0]*y + _dbRotMat[2][0]*z;
        _dbMap[i].point.y = _dbRotMat[0][1]*x + _dbRotMat[1][1]*y + _dbRotMat[2][1]*z;
        _dbMap[i].point.z = _dbRotMat[0][2]*x + _dbRotMat[1][2]*y + _dbRotMat[2][2]*z;

        double dx = _refMap[i].point.x - _dbMap[i].point.x;
        double dx2 = dx * dx;
        double dy = _refMap[i].point.y - _dbMap[i].point.y;
        double dy2 = dy * dy;
        double dz = _refMap[i].point.z - _dbMap[i].point.z;
        double dz2 = dz * dz;
        double sx = _refMap[i].point.x + _dbMap[i].point.x;
        double sx2 = sx * sx;
        double sy = _refMap[i].point.y + _dbMap[i].point.y;
        double sy2 = sy * sy;
        double sz = _refMap[i].point.z + _dbMap[i].point.z;
        double sz2 = sz * sz;

        _AkA[i] = new SiMath::Matrix(4,4,0.0);
        (*_AkA[i])[0][0] = dx2 + dy2 + dz2;
        (*_AkA[i])[0][1] = dy*sz - sy*dz;
        (*_AkA[i])[0][2] = sx*dz - dx*sz;
        (*_AkA[i])[0][3] = dx*sy - sx*dy;
        (*_AkA[i])[1][0] = (*_AkA[i])[0][1];
        (*_AkA[i])[1][1] = dx2 + sy2 + sz2;
        (*_AkA[i])[1][2] = dx*dy - sx*sy;
        (*_AkA[i])[1][3] = dx*dz - sx*sz;
        (*_AkA[i])[2][0] = (*_AkA[i])[0][2];
        (*_AkA[i])[2][1] = (*_AkA[i])[1][2];
        (*_AkA[i])[2][2] = sx2 + dy2 + sz2;
        (*_AkA[i])[2][3] = dy*dz - sy*sz;
        (*_AkA[i])[3][0] = (*_AkA[i])[0][3];
        (*_AkA[i])[3][1] = (*_AkA[i])[1][3];
        (*_AkA[i])[3][2] = (*_AkA[i])[2][3];
        (*_AkA[i])[3][3] = sx2 + sy2 + dz2;

        // Rotate normals
        x = _refMap[i].normal.x;
        y = _refMap[i].normal.y;
        z = _refMap[i].normal.z;
        _refMap[i].normal.x = _refRotMat[0][0]*x + _refRotMat[1][0]*y + _refRotMat[2][0]*z;
        _refMap[i].normal.y = _refRotMat[0][1]*x + _refRotMat[1][1]*y + _refRotMat[2][1]*z;
        _refMap[i].normal.z = _refRotMat[0][2]*x + _refRotMat[1][2]*y + _refRotMat[2][2]*z;
        x = _dbMap[i].normal.x;
        y = _dbMap[i].normal.y;
        z = _dbMap[i].normal.z;
        _dbMap[i].normal.x = _dbRotMat[0][0]*x + _dbRotMat[1][0]*y + _dbRotMat[2][0]*z;
        _dbMap[i].normal.y = _dbRotMat[0][1]*x + _dbRotMat[1][1]*y + _dbRotMat[2][1]*z;
        _dbMap[i].normal.z = _dbRotMat[0][2]*x + _dbRotMat[1][2]*y + _dbRotMat[2][2]*z;
    }
    return;
}

Alignment::~Alignment(void) {
    for (unsigned int i(0); i < _AkA.size(); ++i) {
        if (_AkA[i] != NULL) {
            delete _AkA[i];
        }
    }
}

SolutionInfo Alignment::align(bool n) {
    // create initial solution
    SolutionInfo si;
    si.volume = -1000.0;
    si.iterations = 0;
    si.center1 = _refCenter;
    si.center2 = _dbCenter;
    si.rotation1 = _refRotMat;
    si.rotation2 = _dbRotMat;

    // scaling of the exclusion spheres
    double scale(1.0);
    if (_nbrExcl != 0) {
        scale /= _nbrExcl;
    }
    // try four different start orientations
    for (unsigned int _call(0); _call < 4; ++_call ) {
        // create initial rotation quaternion
        SiMath::Vector rotor(4, 0.0);
        rotor[_call] = 1.0;

        double volume(0.0), oldVolume(-999.99), v(0.0);
        SiMath::Vector dG(4, 0.0);  // gradient update
        SiMath::Matrix hessian(4, 4, 0.0), dH(4, 4, 0.0); // hessian and hessian update
        unsigned int ii(0);
        for (; ii < 100; ++ii) {
            // compute gradient of volume
            _grad = 0.0;
            volume = 0.0;
            hessian = 0.0;
            for (unsigned int i(0); i < _refMap.size(); ++i) {
                // compute the volume overlap of the two pharmacophore points
                SiMath::Vector Aq(4, 0.0);
                SiMath::Matrix *AkA = _AkA[i];
                Aq[0] = (*AkA)[0][0] * rotor[0] + (*AkA)[0][1] * rotor[1] + (*AkA)[0][2] * rotor[2] +
                        (*AkA)[0][3] * rotor[3];
                Aq[1] = (*AkA)[1][0] * rotor[0] + (*AkA)[1][1] * rotor[1] + (*AkA)[1][2] * rotor[2] +
                        (*AkA)[1][3] * rotor[3];
                Aq[2] = (*AkA)[2][0] * rotor[0] + (*AkA)[2][1] * rotor[1] + (*AkA)[2][2] * rotor[2] +
                        (*AkA)[2][3] * rotor[3];
                Aq[3] = (*AkA)[3][0] * rotor[0] + (*AkA)[3][1] * rotor[1] + (*AkA)[3][2] * rotor[2] +
                        (*AkA)[3][3] * rotor[3];
                double qAq = Aq[0] * rotor[0] + Aq[1] * rotor[1] + Aq[2] * rotor[2] + Aq[3] * rotor[3];
                v = GCI2 * pow(PI / (_refMap[i].alpha + _dbMap[i].alpha), 1.5) * exp(-qAq);
                double c(1.0);

                // add normal if AROM-AROM in this case the absolute value of the angle is needed
                if (n && (_refMap[i].func == AROM) && (_dbMap[i].func == AROM)
                    && (_refMap[i].hasNormal) && (_dbMap[i].hasNormal)) {
                    // for aromatic rings only the planar directions count
                    // therefore the absolute value of the cosine is taken
                    c = _normalContribution(_refMap[i].normal, _dbMap[i].normal, rotor);

                    // update based on the sign of the cosine
                    if (c < 0) {
                        c *= -1.0;
                        _dCdq *= -1.0;
                        _d2Cdq2 *= -1.0;
                    }
                    for (unsigned int hi(0); hi < 4; hi++) {
                        _grad[hi] += v * ( _dCdq[hi] - 2.0 * c * Aq[hi] );
                        for (unsigned int hj(0); hj < 4; hj++) {
                            hessian[hi][hj] += v * (_d2Cdq2[hi][hj] - 2.0 * _dCdq[hi]*Aq[hj] + 2.0 * c * (2.0*Aq[hi]*Aq[hj] - (*AkA)[hi][hj]));
                        }
                    }
                    v *= c;
                } else if (n
                    && ((_refMap[i].func == HACC) || (_refMap[i].func == HDON) || (_refMap[i].func == HYBH))
                    && ((_dbMap[i].func == HYBH) || (_dbMap[i].func == HACC)  || (_dbMap[i].func == HDON))
                    && (_refMap[i].hasNormal)
                    && (_dbMap[i].hasNormal)) {
                    // hydrogen donors and acceptor also have a direction
                    // in this case opposite directions have negative impact
                    c = _normalContribution(_refMap[i].normal, _dbMap[i].normal, rotor);
                    for (unsigned int hi(0); hi < 4; hi++) {
                        _grad[hi] += v * ( _dCdq[hi] - 2.0 * c * Aq[hi] );
                        for (unsigned int hj(0); hj < 4; hj++) {
                            hessian[hi][hj] += v * (_d2Cdq2[hi][hj] - 2.0 * _dCdq[hi]*Aq[hj] + 2.0 * c * (2.0*Aq[hi]*Aq[hj] - (*AkA)[hi][hj]));
                        }
                    }
                    v *= c;
                } else if (_refMap[i].func == EXCL) {
                    // scale volume overlap of exclusion sphere with a negative scaling factor
                    // => exclusion spheres have a negative impact
                    v *= -scale;
                    // update gradient and hessian directions
                    for (unsigned int hi=0; hi < 4; hi++) {
                        _grad[hi] -= 2.0 * v * Aq[hi];
                        for (unsigned int hj(0); hj < 4; hj++) {
                            hessian[hi][hj] += 2.0 * v * (2.0*Aq[hi]*Aq[hj] - (*AkA)[hi][hj]);
                        }
                    }
                } else {
                    // update gradient and hessian directions
                    for (unsigned int hi(0); hi < 4; hi++) {
                        _grad[hi] -= 2.0 * v * Aq[hi];
                        for (unsigned int hj(0); hj < 4; hj++) {
                            hessian[hi][hj] += 2.0 * v * (2.0*Aq[hi]*Aq[hj] - (*AkA)[hi][hj]);
                        }
                    }
                }
                volume += v;
            }
            // stop iterations if the increase in volume overlap is too small (gradient ascent)
            // or if the volume is not defined
            if (std::isnan(volume) || (volume - oldVolume < 1e-5)) {
                break;
            }
            // reset old volume
            oldVolume = volume;
            inverseHessian(hessian);
            // update gradient based on inverse hessian
            _grad = rowProduct(hessian,_grad);
            // small scaling of the gradient
            _grad *= 0.9;
            // update rotor based on gradient information
            rotor += _grad;
            // normalise rotor such that it has unit norm
            normalise(rotor);
        }
        // save result in info structure
        if (oldVolume > si.volume) {
            si.rotor = rotor;
            si.volume = oldVolume;
            si.iterations = ii;
        }
    }
    return si;
}

double Alignment::_normalContribution(Coordinate &n1, Coordinate &n2, SiMath::Vector &q) {
    double x = n2.x;
    double y = n2.y;
    double z = n2.z;

    double d1sq(q[1]*q[1]);
    double d2sq(q[2]*q[2]);
    double d3sq(q[3]*q[3]);


    double Ux = x * (1.0 - 2.0 * d2sq - 2.0 * d3sq )
                + y * (2.0 * (q[2] * q[1] - q[0] * q[3]))
                + z * (2.0 * (q[3] * q[1] + q[0] * q[2]));

    double Uy = x * (2.0 * (q[1] * q[2] + q[0] * q[3]))
                + y * (1.0 - 2.0 * d1sq - 2.0 * d3sq)
                + z * (2.0 * (q[3] * q[2] - q[0] * q[1]));

    double Uz = x * (2.0 * (q[1] * q[3] - q[0] * q[2]))
                + y * (2.0 * (q[2] * q[3] + q[0] * q[1]))
                + z * (1.0 - 2.0 * d1sq - 2.0 * d2sq);

    // hessian update matrix
    _d2Cdq2[0][0] = 0.0;
    _d2Cdq2[1][1] = 2.0 * (n1.y * Uy + n1.z * Uz);
    _d2Cdq2[2][2] = 2.0 * (n1.x * Ux + n1.z * Uz);
    _d2Cdq2[3][3] = 2.0 * (n1.x * Ux + n1.y * Uy);
    _d2Cdq2[0][1] = _d2Cdq2[1][0] = -2.0 * (n1.y * Uz - n1.z * Uy);
    _d2Cdq2[0][2] = _d2Cdq2[2][0] =  2.0 * (n1.x * Uz - n1.z * Ux);
    _d2Cdq2[0][3] = _d2Cdq2[3][0] = -2.0 * (n1.x * Uy - n1.y * Ux);
    _d2Cdq2[1][2] = _d2Cdq2[2][1] =  2.0 * (n1.x * Uy + n1.y * Ux);
    _d2Cdq2[1][3] = _d2Cdq2[3][1] =  2.0 * (n1.x * Uz + n1.z * Ux);
    _d2Cdq2[2][3] = _d2Cdq2[3][2] =  2.0 * (n1.y * Uz + n1.z * Uy);

    // gradient update
    _dCdq[0] = _d2Cdq2[0][0] * q[0] + _d2Cdq2[0][1] * q[1] + _d2Cdq2[0][2] * q[2] + _d2Cdq2[0][3] * q[3];
    _dCdq[1] = _d2Cdq2[1][0] * q[0] + _d2Cdq2[1][1] * q[1] + _d2Cdq2[1][2] * q[2] + _d2Cdq2[1][3] * q[3];
    _dCdq[2] = _d2Cdq2[2][0] * q[0] + _d2Cdq2[2][1] * q[1] + _d2Cdq2[2][2] * q[2] + _d2Cdq2[2][3] * q[3];
    _dCdq[3] = _d2Cdq2[3][0] * q[0] + _d2Cdq2[3][1] * q[1] + _d2Cdq2[3][2] * q[2] + _d2Cdq2[3][3] * q[3];

    // return cosine
    return n1.x * Ux + n1.y * Uy + n1.z * Uz;
}

double Alignment::_quatVolumeOverlap(double alpha1, double alpha2, const SiMath::Vector &q, const SiMath::Matrix &A) {
    // compute qTAq
    // first t = Aq
    SiMath::Vector temp = SiMath::colProduct(q,A);
    // next qT*t
    double r2 = temp.dotProd(q);
    double vol = GCI2 * pow( (PI)/(alpha1 + alpha2),1.5);
    vol *= exp(-r2);
    return vol;
}
