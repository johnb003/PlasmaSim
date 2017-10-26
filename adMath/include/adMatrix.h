#ifndef MATRIX_H
#define MATRIX_H
/*
	Animadead: A Skeletal Animation Library
	Copyright (C) 2005 John C. Butterfield
	http://animadead.sourceforge.net

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
	Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

	John C. Butterfield
	johnb003@hotmail.com
*/

/** \file matrix.h
 	\brief Matrix
 	
	This file includes the class declaration for for a 3D transformation matrix.
*/

#include <iostream>

#include "adQuaternion.h"
#include "adVector.h"

namespace ad
{
/// Represents 3D Transformations.
/** Mathematical structure used to store 3D transformations.  Good for storing a
	stack of transformations and can transform points efficiently. This class uses
	row-major indexing.\n

	Since matrices are layed out in a 2D array, the indexing can be done in two ways:
	row-major, and column-major.

	Column-major indexing (like OpenGL)...
		| 0  4  8 12 |
		| 1  5  9 13 |
		| 2  6 10 14 |
		| 3  7 11 15 |

 	Row-major indexing (like DirectX)...
		| 0  1  2  3 |
		| 4  5  6  7 |
		| 8  9 10 11 |
		|12 13 14 15 |

	Regardless of how the values are layed out in memory, the math always works the same. (Row * Column)
		| - - - - |   | - - # - |   | - - - - |
		| # # # # | * | - - # - | = | - - # - |
		| - - - - |   | - - # - |   | - - - - |
		| - - - - |   | - - # - |   | - - - - |

	So, if our vectors are considered row-major, the only valid way to multiply the vector by the matrix is on the left, v * M
	and conversely if our vectors were column-major, the only valid way to multiply would be if the vector was on the right M * v

		Row-Major: ( v * M )						Column-Major: ( M * v )
					| - - # - |						| - - - - |   | x |   | - |
		[x y z w] * | - - # - | = [- - # -]			| # # # # | * | y | = | # |
					| - - # - |						| - - - - |   | z |   | - |
					| - - # - |						| - - - - |   | w |   | - |

	So when it comes to matrix multiplication order, how can you tell which way to multiply?
	Well if you want to transform a vector by a series of matrices, the first matrix it is transformed by is the most local.
	So, to transform by the chain of bones A,B,C (where A is the root):
		In column-major, it's left to right A*B*C*v.
		In row-major, it's right to left v*C*B*A


	A 2D array in C/C++ like this:
	Scalar [][4] matrix = {
		{ a, b, c, d },
		{ e, f, g, h },
		{ i, j, k, l },
		{ m, n, o, p }};

	Will be stored in memory like this: a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p.

	This library uses row major format.

	\sa Vec4,
		Quaternion
*/

	// TODO:
	// - Add an operator[][] accessor
	// - Add a GetColMajor
	// - Add more to the todo list

	class Matrix
	{
	public:
		/// Holds the data for the matrix
		Scalar m[16];

		/// Initialize the matrix to the identity Matrix.
		Matrix();
		/// Initialize the matrix from an existing Matrix.
		Matrix(const Matrix &arg);
		/// Initialize the matrix from an array, each 4 values must define an axis
		Matrix(const Scalar farray[16]);
		/// Initialize the matrix from 3 vectors and an optional position vector.
		Matrix(const Vec4 &xv, const Vec4 &yv, const Vec4 &zv, const Vec4 &pv = Vec4());
		/// Initialize the matrix from an orientation Quaternion and a position vector.
		Matrix(const Quaternion &q, const Vec4 &pv);
		~Matrix() {}

		/// Sets the matrix from an array, with the option to specify the input as either row or col-major;
		void Set(const Scalar farray[16]);
		/// Sets the matrix from 3 vectors and an optional position vector.
		void Set(const Vec4 &xv, const Vec4 &yv, const Vec4 &zv, const Vec4 &pv = Vec4());
		/// Sets the matrix from an orientation Quaternion and a position vector.
		void Set(const Quaternion &q, const Vec4 &pv);

		Vec4 &AxisX() const { return (Vec4 &)m[0]; }
		Vec4 &AxisY() const { return (Vec4 &)m[4]; }
		Vec4 &AxisZ() const { return (Vec4 &)m[8]; }
		Vec4 &Pos() const { return (Vec4 &)m[12]; }

		/// Assigns the matrix to the value of another matrix, and returns the address to this matrix.
		Matrix &operator=(const Matrix &rhs);

		/// Transforms the current matrix by an offset.
		Matrix &operator+=(const Vec4 &rhs);

		/// Applies the rhs transformation to the current matrix, stores the result in the current matrix, and returns the address.
		Matrix &operator*=(const Matrix &rhs);

		/// Transforms the current matrix by another and return the the resulting matrix.
		Matrix operator*(const Matrix &rhs) const;

		/// Transforms the current matrix by another and return the the resulting matrix.
		friend Vec4 operator*(const Vec4 &lhs, const Matrix &rhs);

		/// <x,y,z,1> * Matrix, transforms the point by the matrix, and returns the resulting point
		Vec4 Transform(const Vec4 &rhs) const;
		/// <x,y,z,0> * Matrix, transforms the vector, does not apply translation
		Vec4 Rotate(const Vec4 &v) const;

		/// Converts handedness.  For example if using right-hand and expecting X to the right, Y up, and Z would be towards you, then to convert to left hand, we'll invert Z
		/// This negates all z components of all axies, AND negates the Z axis, (so the z component of the z axis will remain positive)
		Matrix &FlipZ(Matrix &outMat) const;

		/// The transformation matrix that undoes this transformation.
		Matrix Inverse() const;

		/// Swap the rows and columns (used to get access to a column major version)
		Matrix Transpose() const;

	};
	std::ostream &operator<<(std::ostream &os, const Matrix &m);
	Vec4 operator*(const Vec4 &lhs, const Matrix &rhs);
}

#endif // MATRIX_H
