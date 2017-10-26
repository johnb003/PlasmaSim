#ifndef AD_VECTOR_H
#define AD_VECTOR_H
/*
	Animadead: A Skeletal Animation Library
	Copyright (C) 2005 John C. Butterfield

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

/** \file adVector.h
 	\brief Vector

	Details.
*/

#include <iostream>

#include "adScalar.h"

namespace ad
{
	/// 4-Component vector
	/** Mathematical structure used to hold 3D points and vectors.
	\sa ad::Matrix,
	ad::Quaternion
	*/
	class Vec4
	{
	public:
		Scalar x,y,z,w;

		/// Initialize to 0,0,0,0
		Vec4();
		/// Initialize to the values of another vector
		Vec4(const Vec4 &v);
		/// Initialize to \p x,\p y,\p z,\p w
		Vec4(Scalar x, Scalar y, Scalar z, Scalar w=0);
		/// Set to \p x,\p y,\p z,\p w
		void Set(Scalar x, Scalar y, Scalar z, Scalar w=0);

		/// Assign to the values of another vector
		Vec4 &operator=(const Vec4 &v);
		/// Add to this vector, the values of another vector.
		Vec4 &operator+=(const Vec4 &v);
		/// Subtract from this vector the values of another vector.
		Vec4 &operator-=(const Vec4 &v);
		/// Multiply and store the values of this vector by a scalar.
		Vec4 &operator*=(Scalar f);
		/// Divide and store the values of this vector by a scalar.
		Vec4 &operator/=(Scalar f);

		/// Return the values of this vector plus another
		Vec4 operator+(const Vec4 &v) const;
		/// Return the values of this vector minus another
		Vec4 operator-(const Vec4 &v) const;
		/// Return a copy of this vector, negated
		Vec4 operator-() const;
		/// Return the values of this vector times a scalar
		Vec4 operator*(const Scalar f) const;
		/// Return the values of this vector divided by a scalar
		Vec4 operator/(const Scalar f) const;

		/// Return the magnitude of this vector
		Scalar Length3() const;

		/// Return the squared magnitude of this vector
		Scalar Length3Sqr() const;

		/// Return the magnitude of this vector
		Scalar Length4() const;

		/// Return the magnitude of this vector
		Scalar Length4Sqr() const;

		/// Make this a unit-vector in the same direction
		Scalar Normalize3();

		/// Make this a unit-vector in the same direction, returns length, will check distance with threshold and return m_UnitX if less than threshold
		Scalar Normalize3Safe(Scalar theshold = 1e-7);

		/// Make this a unit-vector in the same direction
		Scalar Normalize4();

		/// Cross Product: returns the vector perpendicular to both \p v1 and \p v2 (the fourth component is ignored, and returns 0)
		static Vec4 Cross(const Vec4 &v1, const Vec4 &v2);

		/// Dot Product = | \p v1 || \p v2 | cos(a), cosine of the angle between two unit-vectors
		static Scalar Dot3(const Vec4 &v1, const Vec4 &v2);

		/// Dot Product = | \p v1 || \p v2 | cos(a), cosine of the angle between two unit-vectors
		static Scalar Dot4(const Vec4 &v1, const Vec4 &v2);

		/// Linear Interpolation between two vectors
		static Vec4 Lerp(const Vec4 &from, const Vec4 &to, Scalar t);

		/// {1, 0, 0, 0}
		static const Vec4 m_UnitX;
		/// {0, 1, 0, 0}
		static const Vec4 m_UnitY;
		/// {0, 0, 1, 0}
		static const Vec4 m_UnitZ;
		/// {0, 0, 0, 1}
		static const Vec4 m_UnitW;
		/// {0, 0, 0, 0}
		static const Vec4 m_Zero;
	};
	std::ostream &operator<<(std::ostream &os, const Vec4 &v);
	Vec4 operator*(Scalar f, const Vec4 &rhs);
}

#endif // AD_VECTOR_H
