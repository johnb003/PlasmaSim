#include "CurrentLineField.h"
#include "SDL.h"

#include <math.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#ifdef _WIN32
#include <windows.h>
#endif  // _WIN32
#include <GL/gl.h>
#include <GL/glu.h>
#endif  // Anything other than __APPLE__


using ad::Vec4;
using ad::Scalar;

Scalar gCurrent = 5000000; // current (current per strand * strands)
Scalar gChargePerMeterOfWire = 1.0e-30;

CurrentLineField::CurrentLineField(const ad::Vec4 &pointA, const ad::Vec4 &pointB)
: VectorFunction(), m_pointA(pointA), m_pointB(pointB)
{
	tPos = 0.0;
}

// x is the value along AtoBDir where closestPoint is the origin
Scalar LineIntegralSinDist2At(Scalar x, Scalar z)
{
	return x/sqrt(x*x + z*z);
}

void CurrentLineField::GetField(Vec4 &outBField, Vec4 &outEField, const Vec4 &position) const
{
	const Scalar mu = 1e-7;	// permeability of free space
	const Scalar Ke = 8987551787.3681764; // C^2 e-7

	// In layman's terms:
	//	the magnetic field is proportional to:
	//		the distance squared from the sample point to all points along the current carrying wire
	//		and a factor of how perpendicular each of those points are to the sample point

	// Since I'm using straight line segments, the vector component will be the same everywhere, and I only need to
	// integrate to find the magnitude of the B field. For the E field, I'll solve it component-wise with axial symmetry
	// Therefore I can break the integral into an expression in terms of the angle using a right triangle
	// formed by the closest point on the line defined by the line segment to the sample point.

	// Electric Field:

	//          * (Sample Point)
	//         /|
	//       r/ |z
	//       /  |
	//  A---*-x-*-------------------------B

	//		r = distance from a point on the line to the sample point
	//		z = distance to the line (projected point to line)
	//		x = value integrated, distance from projected point along the line
	//	theta = the angle between the line and r

	// Ez =	Integral( 1/r^2 * sin(theta) )			r^2 = (x^2 + z^2)
	//	  =	Integral( 1/r^2 * z/r )					r	= (x^2 + z^2)^(1/2)
	//	  =	z * Integral( 1/r^(3/2) )
	//	  =	z * Integral( 1/(x^2 + z^2)^(3/2) )
	//	  = 1/z * (x / sqrt(x^2 + z^2) )

	// Ex = Integral( 1/r^2 * cos(theta) dx )
	//	  = Integral( 1/r^2 * x/r dx )
	//	  = Integral( x/(x^2 + z^2)^(3/2) dx )
	//	  = -1 / sqrt(x^2 + z^2)

	// Magnetic Field:

	// using the same diagram, we can define the magnetic field with a few more variables:

	//		rHat = a unit vector where r is depicted pointing towards the sample point
	//		dI   = a vector vector pointing in the direction of the current
	//		u	 = magnetic constant

	//  B = u*I/4*pi * Integral( dI x rHat / r^2 )
	//  B =	u*I/4*pi * Integral( 1/r^2 * sin(theta) )			r^2 = (x^2 + z^2)

	// It can be observed that the cross product dI x rHat is always in the same direction
	// since it differs in magnitude by the sine of the angle, and the distance squared, we have a perfect
	// match for the integral part of Ez


	// Given 2 points to define a line segment A and B (assuming current travels from A to B)
	// find the closest point on the line

	// project A to P onto A to B
	Vec4 AtoB = m_pointB - m_pointA;
	Scalar AtoBLength = AtoB.Length3();
	Vec4 AtoBDir = AtoB / AtoBLength;

	Scalar qLambda = gChargePerMeterOfWire;// / AtoBLength; // total charge per unit of wire (guessed for now)


	Vec4 AtoP = position - m_pointA;

	Scalar dotP = Vec4::Dot3(AtoP, AtoBDir);
	Vec4 closestPoint = m_pointA + AtoBDir * dotP;

	// is the test point on the line??
	Vec4 zDir = position - closestPoint;
	Scalar z = zDir.Normalize3();
	if (z < 0.00001f)
	{
		outBField = Vec4::m_Zero;
		outEField = Vec4::m_Zero;
		return;
	}

	Scalar da = Vec4::Dot3(m_pointA - closestPoint, AtoBDir);	
	Scalar db = Vec4::Dot3(m_pointB - closestPoint, AtoBDir);	

	// now compute the field magnitude
	// 1/z * [x/sqrt(x^2 + z^2)] (a..b)
	Scalar LineIntegralSinDist2 = (db/sqrt(db*db + z*z) - da/sqrt(da*da + z*z)) / z;
	Scalar LineIntegralCosDist2 =  1/sqrt(da*da + z*z) - 1/sqrt(db*db + z*z);

	if (gCurrent > 0.0001)
	{
		Scalar Bfield = LineIntegralSinDist2 * (mu * gCurrent);
		outBField = Vec4::Cross(AtoBDir, zDir) * Bfield;
	}
	else
	{
		outBField = Vec4::m_Zero;
	}

	outEField = qLambda * Ke * (LineIntegralSinDist2 * zDir + LineIntegralCosDist2 * AtoBDir);
//	outField.w = B;
}

void CurrentLineField::DrawField() const
{
	glColor4f(1,1,0,1);
	glBegin(GL_LINES);  
	glVertex3d(m_pointA.x, m_pointA.y, m_pointA.z);
	glVertex3d(m_pointB.x, m_pointB.y, m_pointB.z);
	glEnd();

	Vec4 goopPos = m_pointA + (m_pointB-m_pointA) * tPos;

	ad::Scalar goopHalfSize = 0.001;

	glBegin(GL_QUAD_STRIP);  
	glVertex3d(goopPos.x-goopHalfSize, goopPos.y+goopHalfSize, goopPos.z-goopHalfSize);
	glVertex3d(goopPos.x-goopHalfSize, goopPos.y+goopHalfSize, goopPos.z+goopHalfSize);

	glVertex3d(goopPos.x-goopHalfSize, goopPos.y-goopHalfSize, goopPos.z-goopHalfSize);
	glVertex3d(goopPos.x-goopHalfSize, goopPos.y-goopHalfSize, goopPos.z+goopHalfSize);

	glVertex3d(goopPos.x+goopHalfSize, goopPos.y-goopHalfSize, goopPos.z-goopHalfSize);
	glVertex3d(goopPos.x+goopHalfSize, goopPos.y-goopHalfSize, goopPos.z+goopHalfSize);

	glVertex3d(goopPos.x+goopHalfSize, goopPos.y+goopHalfSize, goopPos.z-goopHalfSize);
	glVertex3d(goopPos.x+goopHalfSize, goopPos.y+goopHalfSize, goopPos.z+goopHalfSize);

	glVertex3d(goopPos.x-goopHalfSize, goopPos.y+goopHalfSize, goopPos.z-goopHalfSize);
	glVertex3d(goopPos.x-goopHalfSize, goopPos.y+goopHalfSize, goopPos.z+goopHalfSize);
	glEnd();

	glBegin(GL_QUADS);  
	glVertex3d(goopPos.x-goopHalfSize, goopPos.y+goopHalfSize, goopPos.z+goopHalfSize);
	glVertex3d(goopPos.x+goopHalfSize, goopPos.y+goopHalfSize, goopPos.z+goopHalfSize);
	glVertex3d(goopPos.x+goopHalfSize, goopPos.y-goopHalfSize, goopPos.z+goopHalfSize);
	glVertex3d(goopPos.x-goopHalfSize, goopPos.y-goopHalfSize, goopPos.z+goopHalfSize);

	glVertex3d(goopPos.x-goopHalfSize, goopPos.y-goopHalfSize, goopPos.z-goopHalfSize);
	glVertex3d(goopPos.x+goopHalfSize, goopPos.y-goopHalfSize, goopPos.z-goopHalfSize);
	glVertex3d(goopPos.x+goopHalfSize, goopPos.y+goopHalfSize, goopPos.z-goopHalfSize);
	glVertex3d(goopPos.x-goopHalfSize, goopPos.y+goopHalfSize, goopPos.z-goopHalfSize);
	glEnd();

}

void CurrentLineField::UpdateField(ad::Scalar dt)
{
	tPos += dt;
	if (tPos > 1.0)
		tPos = 0.0;
}