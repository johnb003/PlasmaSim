#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include "adVector.h"
#include <vector>

class VectorFunction;

// GPU device data
// class DevceData
// {
// public:
// 	float4 *pos;
// // 	float4 *vel;
// };

class Particle
{
	const VectorFunction &fn;

//	ad::Vec4 Particle::GetAccel( const ad::Vec4 &position, const ad::Vec4 &velocity );
	void CurlIntegration( const ad::Vec4 &position, const ad::Vec4 &velocity, const ad::Vec4 &bField, ad::Scalar dt, ad::Vec4 &outPos, ad::Vec4 &outVel );

public:
	Particle(const VectorFunction &fn);

	Particle &operator=(const Particle &rhs)
	{
		pos = rhs.pos;
		vel = rhs.vel;
		mass = rhs.mass;
		return *this;
	}

	void RKStep(ad::Scalar stepSize);
	void Draw(ad::Scalar pSize) const;

	ad::Vec4 pos;
	ad::Vec4 vel;
	ad::Scalar mass;
	ad::Scalar charge;

};

class ParticleSystem
{
	VectorFunction *fieldFn;
	std::vector<Particle> particles;
	static const int MAX_POINTS = 1000;
	static const int NUM_PARTICLES = 100;
	static const int ELEC_ION_RATIO = 1 * NUM_PARTICLES;
	ad::Vec4 points[NUM_PARTICLES][MAX_POINTS];			// a circular history of point locations over time
	int curPoint[NUM_PARTICLES];						// the current "write" point in the circular history
	int startPoint[NUM_PARTICLES];						// the first valid point in the circular history


// 	static const int NUM_THREADS_PER_BLOCK = 256;
// 	DevceData deviceData;
public:
	ParticleSystem(VectorFunction *fieldFn);

	void Draw(ad::Scalar pSize) const;
	void Update(ad::Scalar dt);
	void Reset();
};

#endif // PARTICLE_SYSTEM_H