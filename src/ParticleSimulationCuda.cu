#include "ParticleSimulation.h"

#include "cutil_math.h"

#include <thrust/device_vector.h>




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
	float3 dPos;
	float3 dVel;

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
{
	const float mu = 1e-7;	// permeability of free space
	const float Ke = 8987551787.3681764; // C^2 e-7

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

		// In layman's terms:
		//	the magnetic field is proportional to:
		//		the distance squared from the sample point to all points along the current carrying wire
		//		and a factor of how perpendicular the test vector (sample point to a points on the line) is to the current axis

		// Since I'm using straight line segments, the sample vector component will be constant for the integration
		// and I only need to integrate the magnitude of the B field.
		// The E field I can do component-wise aligned with the axis, so that my values can be shared for the B field
		// Therefore I can break the integral into an expression in terms of the angle using a right triangle
		// formed by the closest point on the line defined by the line segment to the sample point.

		// Electric Field:

		//          * (Sample Point)
		//         /|
		//       r/ |z
		//       /  |
		//  A---*-x-C-------------------------B
		//			 \(Closest point on line)

		//		r = distance from a point on the line to the sample point
		//		z = distance to the line (projected point to line)
		//		x = value integrated, distance from projected point along the line
		//	theta = the angle between the line and r

		// Sum of distances along the axis perpendicular to the line
		// Ez =	Integral( 1/r^2 * sin(theta) )			r^2 = (x^2 + z^2)
		//	  =	Integral( 1/r^2 * z/r )					r	= (x^2 + z^2)^(1/2)
		//	  =	z * Integral( 1/r^(3/2) )
		//	  =	z * Integral( 1/(x^2 + z^2)^(3/2) )
		//	  = 1/z * (x / sqrt(x^2 + z^2) ) [ from A to B in terms of distance from C ]

		// Sum of distances along the axis parallel to the line
		// Ex = Integral( 1/r^2 * cos(theta) dx )
		//	  = Integral( 1/r^2 * x/r dx )
		//	  = Integral( x/(x^2 + z^2)^(3/2) dx )
		//	  = -1 / sqrt(x^2 + z^2) [ from A to B in terms of distance from C ]

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


		// Start by computing the distance squared integrals along the axis ("x"), and perpendicular to the axis ("z")

		// Find the closest point on the line segment AB to P
		float3 AtoB = lineB - lineA;
		float AtoBLen = rsqrtf(AtoB.x*AtoB.x + AtoB.y*AtoB.y + AtoB.z*AtoB.z);
		float3 AtoBDir = AtoB / AtoBLen;
		float3 AtoP = pos - lineA;
		float dotP = dot(AtoP, AtoBDir);

		float3 closestPoint = lineA + AtoBDir * dotP;

		float3 zDir = pos - closestPoint;
		float distToLineSq = dot(zDir, zDir);

		// now to solve the integrals, we need the start and end values for x, which is relative to the center point
		float dA = -dotP;
		float dB = AtoBLen - dotP;

		// when the particle is exactly on the end points of the line segment, the strength approaches infinity.
		// As soon as it's just past the tip (on the inside), the tip starts to counter the rest of the wire so it stays within reason
		// A hack will be to fudge the terms that approach infinity

		const float minDistSq = 1e-10f;

		if (distToLineSq < minDistSq)
		{
			// Usually FALSE
			distToLineSq = minDistSq;	// fudge the calculation to prevent infinity
		}

		float distToLineRecip = rsqrtf(distToLineSq);

		// now finally compute the integrals of distance along the axis, and perpendicular to the axis
		float LineIntegralSinDistSq = (dB * rsqrtf(dB*dB + distToLineSq) - dA * rsqrtf(dA*dA + distToLineSq)) * distToLineRecip;
		float LineIntegralCosDistSq =  1.0f * rsqrtf(dA*dA + distToLineSq) - 1.0f * rsqrt(dB*dB + distToLineSq);

		float3 tmpEField = LineIntegralCosDistSq * AtoBDir;

		// we're along the line, so the z component will be 0
		if (distToLineSq > minDistSq)
		{
			// Usually TRUE
			zDir *= distToLineRecip;
			tmpEField += LineIntegralSinDistSq * zDir;

			// the B field only works if there is distance from the wire. We need it to form the direction.
			float BfieldStrength = LineIntegralSinDistSq * (mu * fieldParms.currentAmperes);
			float3 tmpBField = cross(AtoBDir, zDir);
			tmpBField *= BfieldStrength;

			BField += tmpBField;
		}
		tmpEField *= fieldParms.chargePerPointOnWire * Ke;

		EField += tmpEField;
	}
}

__device__ Derivative SampleDerivative(float dt, const State &sampleState, const FieldParm &fieldParms)
{
	float3 BField, EField;
	GetFields(sampleState.pos, fieldParms, BField, EField);

	const float particleCharge = -1.60217646e-19;	// charge of an electron

	float3 force = (EField + cross(sampleState.vel, BField)) * particleCharge;	// f = q(E + v x B)

	return Derivative(sampleState.vel, force/fieldParms.mass);
}

__device__ State EulerStep(float dt, const State &initialState, const Derivative &initialDerivative)
{
	return State(initialState.pos + dt*initialDerivative.dPos, initialState.vel + dt*initialDerivative.dVel);
}

__device__ void RK45Integrate(const State &initialState, const FieldParm &fieldParms, float dt, float3 &outPos, float3 &outVel, float &errorOut)
{
	//const float C1=0.0f;
	const float C2 = 0.25f;
	const float C3 = 3.0f/8.0f;
	const float C4 = 12.0f/13.0f;
	const float C5 = 1.0f;
	const float C6 = 0.5f;

	const float A21 = 0.25f;
	const float A31 = 3.0f/32.0f, A32 = 9.0/32.0f;
	const float A41 = 1932.0f/2197.0f, A42 = -7200.0f/2197.0f, A43 = 7296.0f/2197.0f;
	const float A51 = 439.0f/216.0f, A52 = -8.0f, A53 = 3680.0f/513.0f, A54 = -845.0f/4104.0f;
	const float A61 = -8.0f/27.0f, A62 = 2.0f, A63 = -3544.0f/2565.0f, A64 = 1859.0f/4104.0f, A65 = -11.0f/40.0f;
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
//	return updatedOrder4;
	outVel = updatedOrder4.vel;
	outPos = updatedOrder4.pos;
}



// This isn't generic, because it includes sampling and summing all of the forces
// It's here so we can call it twice out of laziness, should replace with RK4/5 to get a cheaper error metric
//__device__ void RKStep(float3 initialPos, float3 initialVel, float)

__global__ void	integrateBodies(
	const float4 *position,
	const float3 *velocity,
	const float3 *linePairs,
	const int numPairPoints,
	const int numBodies,
	const float currentAmperes,
	const float chargePerPointOnWire,
	float dt,
	float4 *outPosBuff,
	float3 *outVelBuff,
	float *outErrorBuff)
{
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	if (index >= numBodies)
		return;

	const float3 *pairPoints = linePairs;
// 	__shared__ float3 pairPoints[ParticleSimulation::NUM_THREADS_PER_BLOCK];
// 
// 	if (threadIdx.x < numPairPoints)
// 		pairPoints[threadIdx.x] = linePairs[threadIdx.x];
// 
// 	__syncthreads();

	float4 pos4 = position[index];
//	float3 myPos = make_float3(pos4.x, pos4.y, pos4.z);
	float3 vel = velocity[index];

	// each particle should iterate over the line pairs, and sum the E and B field vectors

	float3 outPos, outVel;
//	float3 outPos2, outVel2;

	FieldParm fieldParm;
	fieldParm.chargePerPointOnWire = chargePerPointOnWire;
	fieldParm.currentAmperes = currentAmperes;
	fieldParm.mass = pos4.w;
	fieldParm.pairPoints = pairPoints;
	fieldParm.numPairs = numPairPoints;

	float errorApprox;
	RK45Integrate(State(make_float3(pos4.x, pos4.y, pos4.z), vel), fieldParm, dt, outPos, outVel, errorApprox);

/*
	float3 k1E, k1B;
	float3 k2E, k2B;
	float3 k3E, k3B;
	float3 k4E, k4B;

	// regular step, outputs to outPos, outVel
	{
		// get K1 from start pos
		GetFields(myPos, pairPoints, numPairPoints, currentAmperes, chargePerPointOnWire, k1B, k1E);

		// Get K2 from projected pos using k1 at (1/2)*t
		integrateMotion(myPos, vel, pos4.w, k1E, k1B, dt/2, outPos, outVel);
		GetFields(outPos, pairPoints, numPairPoints, currentAmperes, chargePerPointOnWire, k2B, k2E);

		// Get K3 from projected pos using k2 at (1/2)*t
		integrateMotion(myPos, vel, pos4.w, k2E, k2B, dt/2, outPos, outVel);
		GetFields(outPos, pairPoints, numPairPoints, currentAmperes, chargePerPointOnWire, k3B, k3E);

		// Get K4 from the projected pos using k3 at the end point
		integrateMotion(myPos, vel, pos4.w, k3E, k3B, dt, outPos, outVel);
		GetFields(outPos, pairPoints, numPairPoints, currentAmperes, chargePerPointOnWire, k4B, k4E);

		integrateMotion(myPos, vel, pos4.w,
			(k1E + 2.0f*(k2E + k3E) + k4E)/6.0f,
			(k1B + 2.0f*(k2B + k3B) + k4B)/6.0f,
			dt, outPos, outVel);
	}*/
	// First half step, outputs to outPos2, outVel2
/*
	dt = dt/2;

	{
		// get K1 from start pos
//		GetFields(myPos, __smem, numPairPoints, currentAmperes, chargePerPointOnWire, k1B, k1E);

		// Get K2 from projected pos using k1 at (1/2)*t
		integrateMotion(myPos, vel, pos4.w, k1E, k1B, dt/2, outPos2, outVel2);
		GetFields(outPos2, (float3*)__smem, numPairPoints, currentAmperes, chargePerPointOnWire, k2B, k2E);

		// Get K3 from projected pos using k2 at (1/2)*t
		integrateMotion(myPos, vel, pos4.w, k2E, k2B, dt/2, outPos2, outVel2);
		GetFields(outPos2, (float3*)__smem, numPairPoints, currentAmperes, chargePerPointOnWire, k3B, k3E);

		// Get K4 from the projected pos using k3 at the end point
		integrateMotion(myPos, vel, pos4.w, k3E, k3B, dt, outPos2, outVel2);
		GetFields(outPos2, (float3*)__smem, numPairPoints, currentAmperes, chargePerPointOnWire, k4B, k4E);

		integrateMotion(myPos, vel, pos4.w,
			(k1E + 2.0f*(k2E + k3E) + k4E)/6.0f,
			(k1B + 2.0f*(k2B + k3B) + k4B)/6.0f,
			dt, outPos2, outVel2);
	}
	// Second half step, outputs to outPos2, outVel2
	{
		// update initial values
		myPos = outPos2;
		vel = outVel2;

		// get K1 from start pos
		GetFields(myPos, (float3*)__smem, numPairPoints, currentAmperes, chargePerPointOnWire, k1B, k1E);

		// Get K2 from projected pos using k1 at (1/2)*t
		integrateMotion(myPos, velocity[index], pos4.w, k1E, k1B, dt/2, outPos2, outVel2);
		GetFields(outPos2, (float3*)__smem, numPairPoints, currentAmperes, chargePerPointOnWire, k2B, k2E);

		// Get K3 from projected pos using k2 at (1/2)*t
		integrateMotion(myPos, velocity[index], pos4.w, k2E, k2B, dt/2, outPos2, outVel2);
		GetFields(outPos2, (float3*)__smem, numPairPoints, currentAmperes, chargePerPointOnWire, k3B, k3E);

		// Get K4 from the projected pos using k3 at the end point
		integrateMotion(myPos, velocity[index], pos4.w, k3E, k3B, dt, outPos2, outVel2);
		GetFields(outPos2, (float3*)__smem, numPairPoints, currentAmperes, chargePerPointOnWire, k4B, k4E);

		integrateMotion(myPos, velocity[index], pos4.w,
			(k1E + 2.0f*(k2E + k3E) + k4E)/6.0f,
			(k1B + 2.0f*(k2B + k3B) + k4B)/6.0f,
			dt, outPos2, outVel2);
	}
*/
	outPosBuff[index].x = outPos.x;
	outPosBuff[index].y = outPos.y;
	outPosBuff[index].z = outPos.z;
	outPosBuff[index].w = pos4.w;
	outVelBuff[index] = outVel;
	outErrorBuff[index] = 0.0f;//dot(outPos - outPos2, outPos - outPos2);

}

// 	unsigned int deviceNumBodies,
// 	float deltaTime,
// 	int totalNumBodies)


// void IntegrateNBodySystem(DeviceData *deviceData, float currentAmperes, float chargePerPointOnWire, float dt, int numBodies)
// {
// 	int numThreadsPerBlock = 64;
// 	int numBlocks = (numBodies-1) / numThreadsPerBlock + 1;
// 
// 	dim3 dimGrid(numBlocks);
// 	dim3 dimBlock(numThreadsPerBlock);
// 
// 	integrateBodies<<< dimGrid, dimBlock >>>(deviceData->particlePos, deviceData->particleVel, deviceData->linePairPos, deviceData->state, deviceData->numPointsInPairs, currentAmperes, chargePerPointOnWire, dt);
// }

float stepSize = 1.0e-10f;

void IntegrateNBodySystem( DeviceData &deviceData, int numBodies, float currentAmperes, float chargePerPointOnWire, int outputIndex, float dt, cudaStream_t stream)
{
	int numThreadsPerBlock = 64;
	int numBlocks = (numBodies-1) / numThreadsPerBlock + 1;

	dim3 dimGrid(numBlocks);
	dim3 dimBlock(numThreadsPerBlock);

	float dtRemaining = dt;

	int retryCount = 0;
	int inputIndex = deviceData.state;

	if (stepSize > dt)
		stepSize = dt;

	const float fltThresh = 1.0e-10f;
	while (dtRemaining > 0.0f/* && retryCount < 10*/)
	{
		retryCount++;
		float step = stepSize;
		if (stepSize > dtRemaining)
			step = dtRemaining;

		integrateBodies<<< dimGrid, dimBlock, 0, stream >>>(
			deviceData.particlePos[inputIndex],
			deviceData.particleVel[inputIndex],
			deviceData.linePairPos,
			deviceData.numPointsInPairs,
			numBodies,
			currentAmperes,
			chargePerPointOnWire,
			step,
			deviceData.particlePos[outputIndex],
			deviceData.particleVel[outputIndex],
			deviceData.integrationStepError);
		inputIndex = outputIndex;

		// wrap raw pointer with a device_ptr 
		thrust::device_ptr<float> dev_ptr(deviceData.integrationStepError);
		float max = thrust::reduce(dev_ptr, dev_ptr + numBodies, -1.0f, thrust::maximum<float>());

		if (max > fltThresh)
		{
			// too inaccurate, try again
			stepSize /= 2;
			continue;
		}


		if (max < fltThresh / 4)
		{
			stepSize *=2;
		}

		dtRemaining -= step;
	} 
//	printf("max: %f\n", max);

}