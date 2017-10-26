#ifndef QUATERNION_H
#define QUATERNION_H
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

/** \file quaternion.h
 	\brief Quaternion
 	
	Details.
*/

#include <iostream>

#include "adScalar.h"

namespace ad
{
/// Stores a 3D rotation, free of gimbal lock
/** Mathematical structure that you shouldn't even try to visualize. These are
	good for interpolating between two rotations and applying successive
	rotations.
*/
	class Matrix;

	class Quaternion
	{
	public:
		Scalar x, y, z, w;

		/// Initialize with no rotation
		Quaternion();
		/// Initialize with the rotation of another quaternion
		Quaternion(const Quaternion &q);
		/// Initialize with \p x,\p y,\p z,\p w
		Quaternion(Scalar x, Scalar y, Scalar z, Scalar w);
		/// Initialize from Euler rotations
		Quaternion(Scalar rotX, Scalar rotY, Scalar rotZ);
		/// Set to \p x,\p y,\p z,\p w
		void Set(Scalar x, Scalar y, Scalar z, Scalar w);
		/// Set from matrix \p mat
		void FromMatrix(const Matrix &mat);
		/// Set to \p x,\p y,\p z,\p w
		void FromAxisAndAngle(Scalar x, Scalar y, Scalar z, Scalar theta);
		/// Set from Euler angles \p yaw, \p pitch, and \p roll
		void FromEuler(Scalar yaw, Scalar pitch, Scalar roll);
		/// Get the Euler angles from the quaternion
		void GetEuler(Scalar &yaw, Scalar &pitch, Scalar &roll);
		/// Assign to the value of another quaternion
		Quaternion &operator=(const Quaternion &q);

		/// Rotate this quaternion by another
		Quaternion &operator*=(const Quaternion &q);
		/// Return the rotation of this quaternion then another, \p q
		Quaternion operator*(const Quaternion &q);

		/// Smoothly interpolates between two UNIT quaternions
		static Quaternion Slerp(const Quaternion &from, const Quaternion &to, Scalar t);
		/// Linearly interpolates between two UNIT quaternions
		static Quaternion Lerp(const Quaternion &from, const Quaternion &to, Scalar t);
	};
	std::ostream &operator<<(std::ostream &os, const Quaternion &q);
}


#endif // QUATERNION_H
