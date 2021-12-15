/*******************************************************************************
siMath.h - Align-it

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

#ifndef __SILICOSIT_ALIGNIT_SIMATH_H__
#define __SILICOSIT_ALIGNIT_SIMATH_H__

// General
#include <math.h>
#include <vector>

// Align-it
#include <defines.h>

#define ROUND(x) ((int) ((x) + 0.5))

#ifndef min
template <class T> inline T min(T x, T y) { return (x < y) ? x : y; }
#endif

#ifndef max
template <class T> inline T max(T x, T y) { return (x > y) ? x : y; }
#endif

#ifndef sign
template <class T> inline T sign(const T &a, const T &b) {
    return (b >= 0.0) ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
#endif

namespace SiMath {

inline double triangle(const double &a, const double &b) {
    double A(fabs(a)), B(fabs(b));
    if (A > B) {
        return A * sqrt(1.0 + (B/A)*(B/A));
    } else if (B == 0) {
        return 0;
    }
    return B * sqrt(1.0 + (A/B)*(A/B));
}

class Vector {
private:
    unsigned int _n;               ///< Number of data points in vector
    std::vector<double> _pVector;  ///< std::vector to hold values
public:
    Vector() : _n(0), _pVector(0) {};
    explicit Vector(const unsigned int n) : _n(n), _pVector(n) {};
    Vector(const unsigned int n, const double &v) : _n(n), _pVector(n, v) {};
    Vector(const unsigned int n, const double*);
    explicit Vector(const std::vector<double>&);

    Vector(const Vector&);

    ~Vector();

    void clear();
    void reset(unsigned int);

    void resize(unsigned int);

    double getValueAt(const unsigned int);
    double getValueAt(const unsigned int) const;

    double max() const;
    double max(unsigned int &) const;
    double min() const;
    double min(unsigned int &) const;
    double sum() const;
    double mean() const;
    double stDev() const;
    double stDev(double m) const;
    unsigned int size() const { return _n; };

    Vector &operator=(const Vector &);
    Vector &operator=(const double &);
    Vector &operator+=(const double &);
    Vector &operator+=(const Vector &);
    Vector &operator-=(const double &);
    Vector &operator-=(const Vector &);
    Vector &operator*=(const double &);
    Vector &operator*=(const Vector &);
    Vector &operator/=(const double &);
    Vector &operator/=(const Vector &);
    Vector &operator-();
    Vector operator+(const Vector &) const;
    Vector operator-(const Vector &) const;
    Vector operator*(const Vector &) const;
    Vector operator/(const Vector &) const;
    
    bool operator==(const Vector &) const;
    bool operator!=(const Vector &) const;

    inline double &operator[](const unsigned int i) {
        return _pVector[i];
    };
    inline double operator[](const unsigned int i) const {
        return _pVector[i];
    };

    void swap(const unsigned int, const unsigned int);

    double dotProd(const Vector &);

    const double *getArrayPointer() const {
        return &(_pVector[0]);
    };
    double *getArrayPointer() {
        return &(_pVector[0]);
    };
};

class Matrix {
private:
    unsigned int _nRows;
    unsigned int _nCols;
    double **_pMatrix;
public:
    Matrix() : _nRows(0), _nCols(0), _pMatrix(NULL){};
    Matrix(const unsigned int, const unsigned int);
    Matrix(const unsigned int, const unsigned int, const double &);
    Matrix(const Matrix &);

    Matrix(const unsigned int, const unsigned int, const Vector &);

    ~Matrix();

    void reset(const unsigned int, const unsigned int);
    void clear();

    inline unsigned int nbrRows() const { return _nRows; };
    inline unsigned int nbrColumns() const { return _nCols; };

    double getValueAt(const unsigned int, const unsigned int);
    const double getValueAt(const unsigned int, const unsigned int) const;
    Vector getRow(const unsigned int) const;
    Vector getColumn(const unsigned int) const;

    void setValueAt(const unsigned int, const unsigned int, double);
    void setRow(const unsigned int, Vector &);
    void setColumn(const unsigned int, Vector &);

    void swapRows(unsigned int, unsigned int);
    void swapColumns(unsigned int, unsigned int);
    Matrix transpose(void);

    Matrix &operator=(const Matrix&);
    Matrix &operator=(const double&);
    Matrix &operator+=(const double&);
    Matrix &operator+=(const Matrix&);
    Matrix &operator-=(const double&);
    Matrix &operator-=(const Matrix&);
    Matrix &operator*=(const double&);
    Matrix &operator*=(const Matrix&);
    Matrix &operator/=(const double&);
    Matrix &operator/=(const Matrix&);
    Matrix &operator-();

    Matrix operator+(const Matrix&) const;
    Matrix operator-(const Matrix&) const;
    Matrix operator*(const Matrix&) const;
    Matrix operator/(const Matrix&) const;

    inline double * operator[] (const unsigned int i) { return _pMatrix[i]; };
    inline const double * operator[] (const unsigned int i) const { return _pMatrix[i]; };
};

Vector rowProduct(const Matrix &A, const Vector &U);
Vector colProduct(const Vector &U, const Matrix &A);

class SVD {
public:
    SVD(const Matrix &, bool bU = true, bool bV = true);

    Vector getSingularValues() { return _S; };
    Matrix getSingularMatrix();

    Matrix getU() { return _U; };
    Matrix getV() { return _V; };

    double norm2() { return _S[0]; };
    double cond() { return _S[0] / _S[_S.size() - 1]; };

    int rank();

private:
    int _m;
    int _n;
    Matrix _U;
    Matrix _V;
    Vector _S;

    bool _computeV;
    bool _computeU;
};

double randD(double, double);

} // namespace SiMath

#endif //__SILICOSIT_ALIGNIT_SIMATH_H__
