#include "ParticleSimulation.h"

#include "cutil_math.h"

#include <thrust/device_vector.h>
#include "adScalar.h"

#include <sstream>
#include <iomanip>




// __device__ void CurlIntegrator(const float4 &initialPosition, const float3 &initialVelocity, const float3 &Bfield, const float3 &Efield, float dt, float4 outPos, float3 outVel)
// {
// 	
// }




// void Particle::CurlIntegration( const Vec4 &position, const Vec4 &velocity, const Vec4 &bField, Scalar dt, Vec4 &outPos, Vec4 &outVel )
// {
// 	Vec4 unitB = bField;
// 	Scalar bFieldStr = unitB.Normalize3();
// 
// 	Vec4 velocityPerpToPlane = unitB * Vec4::Dot3(velocity, unitB);
// 	Vec4 velocityOnPlane = velocity - velocityPerpToPlane;
// 
// 	Scalar velPlaneMagSquared = velocityOnPlane.Length3Sqr();
// 
// 	Scalar accelMag = 0;
// 	Vec4 accel;
// 
// 	if (velPlaneMagSquared > 1e-10)
// 	{
// 		accel = Vec4::Cross(velocity, bField) * (charge / mass);
// 		accelMag = accel.Length3();
// 		Scalar distFromOrbit = velPlaneMagSquared / accelMag;
// 		if (distFromOrbit < 1000)	// this is huge for a particle, If it's more just use Euler instead
// 		{
// 			Vec4 accelDir =  accel / accelMag;
// 			Scalar velPlaneMag = sqrt(velPlaneMagSquared);
// 
// 			Quaternion quat;
// 			quat.FromAxisAndAngle(unitB.x,unitB.y,unitB.z, (accelMag/velPlaneMag)*dt);
// 			ad::Matrix mat = ad::Matrix(quat, ad::Vec4::m_Zero);
// 
// 			Vec4 offset = accelDir * -distFromOrbit;
// 
// 			outVel = mat.Rotate(velocity);	// just rotate the velocity, (no need to integrate it)
// 			outPos = position + (mat.Transform(offset) - offset) + velocityPerpToPlane * dt; // simple Euler integration for the drifting part.
// 			return;
// 		}
// 	}
// 
// 	outVel = velocity;
// 	outPos = position + velocity*dt;
// 	// Fallback to Euler integration
// 	if (accelMag > 0)
// 	{
// 		outVel += accel*dt;
// 		outPos += (0.5*dt*dt)*accel;
// 	}
// 
// }

// __device__ void integrateMotion(float3 initialPosition, float3 initialVelocity, float mass, float3 EFieldAtInitial, float3 BFieldAtInitial, float dt, float3 &outPos, float3 &outVel)
// {
// 	const float particleCharge = -1.60217646e-19;	// charge of an electron
// 
// 	// crude integration:
// 
// 	float3 force = cross(initialVelocity, BFieldAtInitial);
// 	force = (EFieldAtInitial + force) * particleCharge;
// 
// 	float3 acc = force / mass;
// 	outPos = initialPosition + initialVelocity * dt + 0.5*acc*dt*dt;
// 	outVel = initialVelocity + acc * dt;	// w = mass, F/m = a
// }

struct State
{
	float3 pos;
	float3 vel;

	__device__ State(float3 pos, float3 vel)
		: pos(pos), vel(vel)
	{
	}
};

struct Derivative
{
	float3 dPos;	// delta pos. aka velocity
	float3 dVel;	// acceleration
//	float3 turnAcc;	// turning acceleration

	__device__ Derivative(float3 dPos, float3 dVel)
		: dPos(dPos), dVel(dVel)
	{
	}
};

struct FieldParm
{
	const float3 *pairPoints;
	int numPairs;
	float currentAmperes;
	float chargePerPointOnWire;
	float mass;
//	float electronParticleDensity;
};

__device__ Derivative operator*(float f, const Derivative &rhs)
{
	return Derivative(rhs.dPos*f, rhs.dVel*f);
}

__device__ Derivative operator+(const Derivative &lhs, const Derivative &rhs)
{
	return Derivative(lhs.dPos + rhs.dPos, lhs.dVel+rhs.dVel);
}

__device__ void GetFields( const float3 &pos, const FieldParm &fieldParms, float3 &BField, float3 &EField ) 
// TODO: Explore rewriting the math for the distance to line computation using the fact that:
//
// Alternative Math:
// 
// Area = |PB x PA| / 2   -- which is useful because we need that cross product for B field direction
// Area = b * h/2    -- ie. base * height /2
// solve for h
// |PB x PA|/2 = b * h/2
// |PB x PA| = b * h
// |PB x PA| / |AB| = h       -- cause |AB| is the base.
//
// It would go something like...	
// float3 bFieldDir = cross(lineA - pos, lineB - pos);
// float area2 = sqrt(dot(bFieldDir, bFieldDir));
// float z = area2 / AtoBLen;
// 
// it'd be great to normalize the bFieldDir, but, we can't do that if we're on the line.
// if (z > minDist) {
//     bFieldDir /= area2;
// }
//
// 1/z * [x / (sqrt(x^2 + z^2) )] ( from A to B )
// 1/z * (dB / dot(PB, PB) - dA / dot(PA, PA));
// ... I'm trailing off here, dB and dA aren't available without what we had before.... Worth it???

{
	const float mu = 1e-7;	// magnetic permeability of free space
	const float Ke = 8987551787.3681764; // C^2E-7

	BField.x = 0.0f;
	BField.y = 0.0f;
	BField.z = 0.0f;
	EField.x = 0.0f;
	EField.y = 0.0f;
	EField.z = 0.0f;

	int lineIdx = 0;
	while (lineIdx < fieldParms.numPairs)
	{
		float3 lineA = fieldParms.pairPoints[lineIdx];
		float3 lineB = fieldParms.pairPoints[lineIdx+1];
		lineIdx += 2;

		// The magnetic field is proportional to:
		//    * One over distance squared from a sample point to all points along the current carrying wire
		//    * A factor of how perpendicular the test vector (sample point to a points on the line) is
		//      to the current axis

		// Since the line segments are straight, the current vector will be constant for the integration
		// and only the the magnitude of the B field needs to be integrated. So that's integrating using
		// just z in the figure below.
		//
		// For the E field, we solve each component separately, and the integral along z is the same as
		// the B field.

		//          * (Sample Point)
		//         /|
		//       r/ |z
		//       /  |
		//  A---*-x-C-------------------------B
		//           \(Closest point on line)
		//
		//		r = distance from a point on the line to the sample point
		//		z = distance to the line (projected point to line)
		//		x = value integrated, distance from projected point along the line
		//	theta = the angle between the line and r

		// Sum of distances along the axis perpendicular to the line
		// Note, sin(theta) here (which is z/r) is to select the z component.
		// Integration of 1/r^2 along z component 
		// Ez =	Integral(dx / r^2 * sin(theta))       r^2 = (x^2 + z^2)
		//    = Integral(dx / r^2 * (z/r))            r = (x^2 + z^2)^(1/2)
		//    = z * Integral(dx / r^3)                z is constant, so pull it out of the integral
		//    = z * Integral(dx / (x^2 + z^2)^(3/2))
		// Wolfram Alpha computed the integral for me.
		//    = z * [ x / (z^2 * sqrt(x^2 + z^2) )] [ from A to B in terms of distance from C ]
		//   Note: some z's cancel
		//
		//    = [x / (z * sqrt(x^2 + z^2) )] ( from A to B )
		// Ez = 1/z * [x / sqrt(x^2 + z^2)]
		//    = 1/z * [(b / sqrt(b^2 + z^2)) - (-a / sqrt(a^2 + z^2))]

		// Sum of distances along the axis parallel to the line
		// Ex = Integral( dx / r^2 * cos(theta) )
		//    = Integral( dx / r^2 * (x/r) )
		//    = Integral( x * dx / (x^2 + z^2)^(3/2) )
		// Wolfram Alpha, because math.
		//
		// Ex = [ -1 / sqrt(x^2 + z^2) ] (from A to B in terms of distance from C)

		// Magnetic Field:

		// using the same diagram, we can define the magnetic field with a few more variables:

		//  rHat = a unit vector where r is depicted pointing towards the sample point
		//  dI   = a unit vector pointing in the direction of the current
		//  u    = magnetic constant

		//  B = u*I/4*pi * Integral( dI x rHat * dx / r^2 )
		//                           magnitude of the cross product is sin(theta) from x, to AtoP
		//                           sin(theta) = z/r
		//  B =	u*I/4*pi * Integral( sin(theta) * dx / r^2 )			r^2 = (x^2 + z^2)

		// It can be observed that the cross product dI x rHat is always in the same direction
		// since it differs in magnitude by the sine of the angle, and the distance squared, we have a perfect
		// match for the integral part of Ez

		// Find the closest point on the line segment AB to P by projecting AP onto AB
		// and multiplying by the unit vector in the direction of AB.
		// technically speaking we're doing:
		// AP•AB / |AB| to get the scalar projection.
		// Yet when converted from a scalar to a vector projection we get:
		// AP•AB      AB       AP•AB * AB
		// -----  *  ----  =  -----------
		// |AB|      |AB|        AB•AB
		//
		// This requires less expensive math:
		//		AP•AB/(AB•AB) * AB
		//
		// But this gives me more useful values I will need later:
		// 		AP•(AB/|AB|) * (AB/|AB|)

		float3 AtoB = lineB - lineA;					// AB
		float AtoBLen = sqrtf(dot(AtoB, AtoB));			// |AB|
		float3 AtoBDir = AtoB / AtoBLen;				// AB/|AB|

		// project to unit vector, gives nice distance along unit vector.
		float dotP = dot(pos - lineA, AtoBDir);			// Distance from A to C
		float3 closestPoint = lineA + AtoBDir * dotP;	// C

		// In order to solve the integrals, we need the start and end values for x, which is relative to the center point
		float dA = -dotP;								// we already got the scalar for CtoA in the projection.
		float dB = dA + AtoBLen;						// And we can just compute CtoB from AtoB - AtoC

		// zDir is useful for both distance of z, and since we have AtoBdir, we now have unit vectors for the components we're working with.
		float3 zDir = pos - closestPoint;				// CP

		float distToLineSq = dot(zDir, zDir);			// z^2
		// we'll normalize zDir later, if we're sure it's not too small.

		// Ex = [ -1 / sqrt(x^2 + z^2) ] (from A to B in terms of distance from C)
		// Written backwards because it's negated.
		float LineIntegralX = -rsqrtf(dA*dA + distToLineSq) + rsqrtf(dB*dB + distToLineSq);

		// Bfield has no affect from the X 1/r^2, since it's perpendicular to that plane.
		const float minDistSq = 1e-5f;		

		// If we're effectively colliding with the line, we could say E-field attraction is REALLY high, and
		// B-field is super strong, but we can't tell in what direction, so we need to be careful here.
		// Also, this would only be true for a very small part of the particles integration, and that requries
		// an infitesimally small DT to be accurate. It'd be nice if we didn't have to deal with such a bad
		// vertical asymtote.
		// It should have a velocity to continue moving, so we'll just make this a tiny dead-zone
		float3 bFieldDir;
		float distToLineRecip;
		float LineIntegralZ;
		if (distToLineSq > minDistSq)
		{
			distToLineRecip = rsqrtf(distToLineSq+0.000001f);
			zDir *= distToLineRecip;
			// the B field only works if there is distance from the wire. We need it to form the direction.
			bFieldDir = cross(AtoBDir, zDir);
	
			// This integral only works if the distance of z > ~0
			// [x / (z * sqrt(x^2 + z^2) )] ( from A to B )
			// 1/z * [x / (sqrt(x^2 + z^2) )] ( from A to B )
			LineIntegralZ = distToLineRecip * (dB * rsqrtf(dB*dB + distToLineSq) - dA * rsqrtf(dA*dA + distToLineSq));
		}
		else
		{
			distToLineRecip = rsqrtf(distToLineSq+0.000001f);
			zDir = make_float3(0.0f);
			bFieldDir = make_float3(0.0f);

			// This integral only works if the distance of z > ~0
			// [x / (z * sqrt(x^2 + z^2) )] ( from A to B )
			// 1/z * [x / (sqrt(x^2 + z^2) )] ( from A to B )
			LineIntegralZ = distToLineRecip * (dB * rsqrtf(dB*dB + distToLineSq) - dA * rsqrtf(dA*dA + distToLineSq));
		}
		EField += (LineIntegralX * AtoBDir + LineIntegralZ * zDir) * fieldParms.chargePerPointOnWire * Ke; 

		// No point in doing this unless the B field has an effect.
		BField += bFieldDir * LineIntegralZ * (mu * fieldParms.currentAmperes);

	}
}

__device__ Derivative SampleDerivative(float dt, const State &sampleState, const FieldParm &fieldParms)
{
	float3 BField, EField;
	GetFields(sampleState.pos, fieldParms, BField, EField);

	const float particleCharge = -1.60217646e-19;	// charge of an electron (this should be read per-particle)

	// f = q(E + v x B)
	float3 force = (EField + cross(sampleState.vel, BField)) * particleCharge;

	// TODO: body-body interactions.
	// TODO: particle movement induced fields.

	return Derivative(sampleState.vel, force/fieldParms.mass);
}

__device__ State EulerStep(float dt, const State &initialState, const Derivative &initialDerivative)
{
	return State(initialState.pos + dt*initialDerivative.dPos, initialState.vel + dt*initialDerivative.dVel);
}

__device__ void RK45Integrate(const State &initialState, const FieldParm &fieldParms, float dt, float3 &outPos, float3 &outVel, float &errorOut)
{
	//const float C1=0.0f;
	static const float C2 = 0.25f;
	static const float C3 = 3.0f/8.0f;
	static const float C4 = 12.0f/13.0f;
	static const float C5 = 1.0f;
	static const float C6 = 0.5f;

	static const float A21 = 0.25f;
	static const float A31 = 3.0f/32.0f, A32 = 9.0/32.0f;
	static const float A41 = 1932.0f/2197.0f, A42 = -7200.0f/2197.0f, A43 = 7296.0f/2197.0f;
	static const float A51 = 439.0f/216.0f, A52 = -8.0f, A53 = 3680.0f/513.0f, A54 = -845.0f/4104.0f;
	static const float A61 = -8.0f/27.0f, A62 = 2.0f, A63 = -3544.0f/2565.0f, A64 = 1859.0f/4104.0f, A65 = -11.0f/40.0f;
	// -------------------------------------------------------------------------------------------------------------------------
	const float B4_1 = 25.0f/216.0f, B4_2 = 0.0f, B4_3 = 1408.0f/2565.0f, B4_4 = 2197.0f/4104.0f, B4_5 = -1.0f/5.0f, B4_6 = 0.0f;
	const float B5_1 = 16.0f/135.0f, B5_2 = 0.0f, B5_3 = 6656.0f/12825.0f,B5_4 = 28561.0f/56430.0f,B5_5 = -9.0f/50.0f, B5_6 = 2.0f/55.0f;

	Derivative k1 = SampleDerivative(0.0f, initialState, fieldParms);
	Derivative k2 = SampleDerivative(C2*dt, EulerStep(dt, initialState, A21*k1), fieldParms);
	Derivative k3 = SampleDerivative(C3*dt, EulerStep(dt, initialState, A31*k1 + A32*k2), fieldParms);
	Derivative k4 = SampleDerivative(C4*dt, EulerStep(dt, initialState, A41*k1 + A42*k2 + A43*k3), fieldParms);
	Derivative k5 = SampleDerivative(C5*dt, EulerStep(dt, initialState, A51*k1 + A52*k2 + A53*k3 + A54*k4), fieldParms);
	Derivative k6 = SampleDerivative(C6*dt, EulerStep(dt, initialState, A61*k1 + A62*k2 + A63*k3 + A64*k4 + A65*k5), fieldParms);

	// ...
	//Derivative kn = SampleDerivative(Cn*dt, EulerStep(initialState, An1*k1 + An2*k2 + ... + An_(n-1)*k_(n-1), dt));

	const Derivative deltaSum4 = B4_1*k1 + B4_2*k2+ B4_3*k3 + B4_4*k4 + B4_5*k5 + B4_6*k6;
	const Derivative deltaSum5 = B5_1*k1 + B5_2*k2 + B5_3*k3 + B5_4*k4 + B5_5*k5 + B5_6*k6;

	// For Runge-Kutta on Wikipedia, this final step is done differently because each Ki is stored as dt*derivative, instead of just the derivative as I have done.
	// What that means for the math version is that the final weighted average can just be added to the position. The only difference is that I've saved the dt for last
	// and now I must account for dt, so I can use a normal Eulerstep and do y = y0 + derivative*dt;
	State updatedOrder5 = EulerStep(dt, initialState, deltaSum5);
	State updatedOrder4 = EulerStep(dt, initialState, deltaSum4);

	float3 deltaActualApprox = updatedOrder5.pos - updatedOrder4.pos;
	errorOut = dot(deltaActualApprox, deltaActualApprox);
	outVel = updatedOrder5.vel;
	outPos = updatedOrder5.pos;
}


__global__ void	integrateBodies(
	const float4 *position,
	const float3 *velocity,
	const float3 *linePairs,
	const int numPairPoints,
	const int numBodies,
	const float currentAmperes,
	const float chargePerPointOnWire,
//	const float electronParticleDensity,
	float dt,
	float4 *outPosBuff,
	float3 *outVelBuff,
	float *outErrorBuff)
{
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	if (index >= numBodies)
		return;

	const float3 *pairPoints = linePairs;

	float4 pos4 = position[index];
	float3 vel = velocity[index];

	// each particle should iterate over the line pairs, and sum the E and B field vectors

	float3 outPos, outVel;

	FieldParm fieldParm;
	fieldParm.chargePerPointOnWire = chargePerPointOnWire;
	fieldParm.currentAmperes = currentAmperes;
	fieldParm.mass = 9.10938188e-31f;
	fieldParm.pairPoints = pairPoints;
	fieldParm.numPairs = numPairPoints;
//	fieldParm.electronParticleDensity = electronParticleDensity;

	float errorApprox;
	RK45Integrate(State(make_float3(pos4.x, pos4.y, pos4.z), vel), fieldParm, dt*pos4.w, outPos, outVel, errorApprox);

	outPosBuff[index].x = outPos.x;
	outPosBuff[index].y = outPos.y;
	outPosBuff[index].z = outPos.z;
	outPosBuff[index].w = pos4.w;
	outVelBuff[index] = outVel;
	outErrorBuff[index] = errorApprox;//dot(outPos - outPos2, outPos - outPos2);
}

static float inflicted = 1.0f;
int frame = 0;

bool IntegrateNBodySystem( DeviceData &deviceData, int numBodies, float currentAmperes, float chargePerPointOnWire,
	/*float electronParticleDensity, */int outputIndex, float dt, cudaStream_t stream)
{
	int numThreadsPerBlock = 64;
	static const int NUM_RETRIES = 10;
	int numBlocks = (numBodies-1) / numThreadsPerBlock + 1;

	dim3 dimGrid(numBlocks);
	dim3 dimBlock(numThreadsPerBlock);

	// We'll retry with successively smaller timesteps until we reach enough accuracy.
	int retryCount = 0;
	int inputIndex = deviceData.state;

	const float fltThresh = 1.0e-10f;
	for (; retryCount < NUM_RETRIES; retryCount++)
	{
		integrateBodies<<< dimGrid, dimBlock, 0, stream >>>(
			deviceData.particlePos[inputIndex],
			deviceData.particleVel[inputIndex],
			deviceData.linePairPos,
			deviceData.numPointsInPairs,
			numBodies,
			currentAmperes,
			chargePerPointOnWire,
			// electronParticleDensity,
			dt*inflicted,
			deviceData.particlePos[outputIndex],
			deviceData.particleVel[outputIndex],
			deviceData.integrationStepError);

		// wrap raw pointer with a device_ptr 
		thrust::device_ptr<float> dev_ptr(deviceData.integrationStepError);
		float max = thrust::reduce(dev_ptr, dev_ptr + numBodies, -1.0f, thrust::maximum<float>());

		if (inflicted < 1.0f)
		{
			std::stringstream ss;
			std::cout << std::setprecision(8) << "Dt reduced (" << inflicted*100 << "%) to maintain error threshold less than: " << fltThresh << ".\n";
			std::cout << "Largest integration error: " << max << std::endl;
		}

		if (max < fltThresh*1.0e-2 && inflicted < 1.0f) {
			inflicted *= 2;
		}

		if (max < fltThresh)
			break;

		// too inaccurate, try again
		inflicted /= 2;
	}
	if (retryCount == NUM_RETRIES) return false;
	return true;
}