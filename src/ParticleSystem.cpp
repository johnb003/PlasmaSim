#include "ParticleSystem.h"
#include "VectorFunction.h"

#include "adQuaternion.h"
#include "adMatrix.h"

#include "SDL.h"

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

#include <assert.h>
#include <time.h>
#include <math.h>


using ad::Vec4;
using ad::Scalar;
using ad::Quaternion;
using ad::Matrix;

Particle::Particle(const VectorFunction &fn)
	: fn(fn)
{
	mass = 9.10938188e-31; // electron
//	mass = 3.34449428e-27; // deuterium

	charge = -1.602e-19;	// charge of electron

}

// Vec4 Particle::GetAccel( const Vec4 &position, const Vec4 &velocity )
// {
// 	Vec4 tmpBField;
// 	Scalar q = -1.602e-19;	// charge of electron
// 	return Vec4::Cross(velocity, fn.GetField(tmpBField, position)) * (q / mass);
// }

void Particle::CurlIntegration(const Vec4 &position, const Vec4 &velocity, const Vec4 &bField, Scalar dt, Vec4 &outPos, Vec4 &outVel)
{
	Vec4 unitB = bField;
	unitB.Normalize3();

	Vec4 velocityPerpToPlane = unitB * Vec4::Dot3(velocity, unitB);
	Vec4 velocityOnPlane = velocity - velocityPerpToPlane;

	Scalar velPlaneMagSquared = velocityOnPlane.Length3Sqr();

	Scalar accelMag = 0;
	Vec4 accel;

	if (velPlaneMagSquared > 1e-10)
	{
		accel = Vec4::Cross(velocity, bField) * (charge / mass);
		accelMag = accel.Length3();
		Scalar distFromOrbit = velPlaneMagSquared / accelMag;
		if (distFromOrbit < 1000)	// this is huge for a particle, If it's more just use Euler instead
		{
			Vec4 accelDir =  accel / accelMag;
			Scalar velPlaneMag = sqrt(velPlaneMagSquared);

			Quaternion quat;
			quat.FromAxisAndAngle(unitB.x,unitB.y,unitB.z, (accelMag/velPlaneMag)*dt);
			ad::Matrix mat = ad::Matrix(quat, ad::Vec4::m_Zero);

			Vec4 offset = accelDir * -distFromOrbit;
	
			outVel = mat.Rotate(velocity);	// just rotate the velocity, (no need to integrate it)
			outPos = position + (mat.Transform(offset) - offset) + velocityPerpToPlane * dt; // simple Euler integration for the drifting part.
			return;
		}
	}

	outVel = velocity;
	outPos = position + velocity*dt;
	// Fallback to Euler integration
	if (accelMag > 0)
	{
		outVel += accel*dt;
		outPos += (0.5*dt*dt)*accel;
	}
	
}

void Particle::RKStep(Scalar dt)
{
	// F = q[E + (v x B)]
	//	F is the force (in newtons)
	//	E is the electric field (in volts per meter)
	//	B is the magnetic field (in teslas)
	//	q is the electric charge of the particle (in coulombs)
	//	v is the instantaneous velocity of the particle (in meters per second)
	//	× is the vector cross product 

	// B field samples
	Vec4 k1B;
	Vec4 k2B;
	Vec4 k3B;
	Vec4 k4B;

	// E field samples
	Vec4 k1E;
	Vec4 k2E;
	Vec4 k3E;
	Vec4 k4E;

	// F = qE
	// a = q/m * E
	// v = v0 + at;
	// p = p0 + v0*t + 1/2*(aB + aE)*t^2
	// the part in the [] is handled by the curl integration, so we just need to add the second anti-derivative of acc to get position change
	// p = [p0 + v0*t + 1/2*aB*t^2] + 1/2*aE*t^2,

	Vec4 tmpVel;
	Vec4 tmpPos;
	Vec4 aE;

	// k1 is the slope at the beginning of the interval
	fn.GetField(k1B, k1E, pos);

	// k2 is the slope at the midpoint from k1
	CurlIntegration(pos, vel, k1B, dt/2, tmpPos, tmpVel);
	aE = (charge/mass)*k1E;
	tmpVel += aE*(dt/2);
	tmpPos += aE*(dt*dt/8);	// 0.5 * (t/2)^2 = t^2/(2^2 * 2) = t^2/8
	fn.GetField(k2B, k2E, tmpPos);

	// k3 is the slope at the midpoint from k2
	CurlIntegration(pos, vel, k2B, dt/2, tmpPos, tmpVel);
	aE = (charge/mass)*k2E;
	tmpVel += aE*(dt/2);
	tmpPos += aE*(dt*dt/8);	// 0.5 * (t/2)^2 = t^2/(2^2 * 2) = t^2/8
	fn.GetField(k3B, k3E, tmpPos);

	// k4 is the slope at the endpoint from k3
	CurlIntegration(pos, vel, k3B, dt, tmpPos, tmpVel);
	aE = (charge/mass)*k3E;
	tmpVel += aE*dt;
	tmpPos += 0.5*aE*dt*dt;
	fn.GetField(k4B, k4E, tmpPos);

	CurlIntegration(pos, vel, (k1B + 2*(k2B + k3B) + k4B)/6, dt, pos, vel);
	aE = (charge/mass)*(k1E + 2*(k2E + k3E) + k4E)/6;
	vel += aE*dt;
	pos += 0.5*aE*dt*dt;
}

void Particle::Draw(Scalar pSize) const
{
	ad::Scalar goopHalfSize = 0.01 * pSize;

	glColor4d(1,1,1,1);
	glBegin(GL_QUAD_STRIP);  
	glVertex3d(pos.x-goopHalfSize, pos.y+goopHalfSize, pos.z-goopHalfSize);
	glVertex3d(pos.x-goopHalfSize, pos.y+goopHalfSize, pos.z+goopHalfSize);

	glVertex3d(pos.x-goopHalfSize, pos.y-goopHalfSize, pos.z-goopHalfSize);
	glVertex3d(pos.x-goopHalfSize, pos.y-goopHalfSize, pos.z+goopHalfSize);

	glVertex3d(pos.x+goopHalfSize, pos.y-goopHalfSize, pos.z-goopHalfSize);
	glVertex3d(pos.x+goopHalfSize, pos.y-goopHalfSize, pos.z+goopHalfSize);

	glVertex3d(pos.x+goopHalfSize, pos.y+goopHalfSize, pos.z-goopHalfSize);
	glVertex3d(pos.x+goopHalfSize, pos.y+goopHalfSize, pos.z+goopHalfSize);

	glVertex3d(pos.x-goopHalfSize, pos.y+goopHalfSize, pos.z-goopHalfSize);
	glVertex3d(pos.x-goopHalfSize, pos.y+goopHalfSize, pos.z+goopHalfSize);
	glEnd();

	glBegin(GL_QUADS);  
	glVertex3d(pos.x-goopHalfSize, pos.y+goopHalfSize, pos.z+goopHalfSize);
	glVertex3d(pos.x+goopHalfSize, pos.y+goopHalfSize, pos.z+goopHalfSize);
	glVertex3d(pos.x+goopHalfSize, pos.y-goopHalfSize, pos.z+goopHalfSize);
	glVertex3d(pos.x-goopHalfSize, pos.y-goopHalfSize, pos.z+goopHalfSize);

	glVertex3d(pos.x-goopHalfSize, pos.y-goopHalfSize, pos.z-goopHalfSize);
	glVertex3d(pos.x+goopHalfSize, pos.y-goopHalfSize, pos.z-goopHalfSize);
	glVertex3d(pos.x+goopHalfSize, pos.y+goopHalfSize, pos.z-goopHalfSize);
	glVertex3d(pos.x-goopHalfSize, pos.y+goopHalfSize, pos.z-goopHalfSize);
	glEnd();
}

ParticleSystem::ParticleSystem(VectorFunction *fieldFn)
	: fieldFn(fieldFn), particles(NUM_PARTICLES, Particle(*fieldFn))
{
	srand((unsigned int)time(NULL));
	Reset();
}

float RandSpread(float center, float radius)
{
	return center - radius + (rand() / RAND_MAX) * radius * 2;
}

void ParticleSystem::Reset()
{
	int i = 0;
	for (auto it = particles.begin(); it != particles.end(); ++it)
	{
		it->pos = Vec4(RandSpread(-0.1f, 0.2), RandSpread(0, 0.2), RandSpread(0, 0.2));
		it->vel = Vec4(RandSpread(0, 0), RandSpread(0, 0), RandSpread(0, 0));

		curPoint[i] = 0;
		startPoint[i] = 0;
		i++;
	}
}

void ParticleSystem::Update(ad::Scalar dt)
{
	dt /= 100000000;
	std::vector<Particle>::iterator it;
	int i = 0;
	for (it = particles.begin(); it != particles.end(); ++it, ++i)
	{
		it->RKStep(dt);

		points[i][curPoint[i]] = it->pos;
		curPoint[i]++;

		if (curPoint[i] == MAX_POINTS)
			curPoint[i] = 0;
		if (curPoint[i] == startPoint[i])
			startPoint[i]++;
		if (startPoint[i] == MAX_POINTS)
			startPoint[i] = 0;
	}
}

#if 0  // Adaptive timestep. Useful in offline baked simulation, but not for real-time.
void ParticleSystem::Update(ad::Scalar dt)
{
	dt /= 100000000;
	std::vector<Particle>::iterator it;
	int i = 0;
	for (it = particles.begin(); it != particles.end(); ++it, i++)
	{
		static const Scalar MIN_STEP_SIZE = 0.0000000001;
		static const Scalar MAX_STEP_SIZE = 1.0;

		Scalar stepSize = 1.0f;
		Scalar ALLOWED_ERROR_THRESH = 0.0000001;

		Particle stepParticle = *it;

		Scalar dtConsumed = 0.0;
		int stepsTaken = 0;

		while (dtConsumed < 1 && stepsTaken < 1000)
		{
			Scalar step = dt * stepSize;

			Particle testParticle = stepParticle;
			testParticle.RKStep(step);

			Particle testParticleHalfStep = stepParticle;
			testParticleHalfStep.RKStep(step/2);
			testParticleHalfStep.RKStep(step/2);

			// estimate the local error
			Vec4 errorDiff = testParticleHalfStep.pos - testParticle.pos;
			Scalar errorSum = fabs(errorDiff.x) + fabs(errorDiff.y) + fabs(errorDiff.z);

			if (errorSum > ALLOWED_ERROR_THRESH && stepSize > MIN_STEP_SIZE)
			{
				// try a smaller stepsize
// 				Scalar offByPowersOfTen = max(1.0, log10(errorSum/(ALLOWED_ERROR_THRESH)));
// 				Scalar reductionFactor = pow(0.5, offByPowersOfTen);

				stepSize /= 2;
				if (stepSize < MIN_STEP_SIZE)
					stepSize = MIN_STEP_SIZE;

				// retry
				continue;
			}
			else if (errorSum < (ALLOWED_ERROR_THRESH * 0.01))
			{
				stepSize *= 2;
				if (stepSize > MAX_STEP_SIZE)
					stepSize = MAX_STEP_SIZE;
			}

			stepParticle = testParticleHalfStep;
			dtConsumed += stepSize;
			stepsTaken++;
			points[i][curPoint[i]] = stepParticle.pos;
			curPoint[i]++;
			if (curPoint[i] == MAX_POINTS)
				curPoint[i] = 0;
			if (curPoint[i] == startPoint[i])
				startPoint[i]++;
			if (startPoint[i] == MAX_POINTS)
				startPoint[i] = 0;
		}
		*it = stepParticle;
	}
}
#endif

void ParticleSystem::Draw(Scalar pSize) const
{
	std::vector<Particle>::const_iterator it;
	for (it = particles.begin(); it != particles.end(); ++it)
	{
		it->Draw(pSize);
	}

#ifdef SHOW_TRAILS
	for (int j = 0; j < NUM_PARTICLES; j++)
	{
		if (j<ELEC_ION_RATIO)
			glColor4d(1,1,1,1);
		else
			glColor4d(1,0,0,1);

		if (curPoint[j] < startPoint[j])
		{
			glBegin(GL_LINE_STRIP);
			for (int i = startPoint[j]; i < MAX_POINTS; i++)
			{
				glVertex3d(points[j][i].x, points[j][i].y, points[j][i].z);
			}
			for (int i = 0; i < curPoint[j]; i++)
			{
				glVertex3d(points[j][i].x, points[j][i].y, points[j][i].z);
			}
			glEnd();
		}
		else
		{
			glBegin(GL_LINE_STRIP);
			for (int i = startPoint[j]; i < curPoint[j]; i++)
			{
				glVertex3d(points[j][i].x, points[j][i].y, points[j][i].z);
			}
			glEnd();
		}
	}
#endif   // SHOW_TRAILS
}