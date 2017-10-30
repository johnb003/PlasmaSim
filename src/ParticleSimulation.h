#ifndef PARTICLE_SIMULATION_H
#define PARTICLE_SIMULATION_H

#include <vector_types.h>  // cuda vectors
#include <driver_types.h>
#include "sim_precision.h"

class CurrentLineVecField;

struct DeviceData
{
public:
	// Device Data
	SIMREAL4 *particlePos[2];
	SIMREAL3 *particleVel[2];
	SIMREAL *integrationStepError;
	float3 *linePairPos;
	int numPointsInPairs;
	// we're using SIMREAL buffers, so the inputs are preserved.
	// this allows us to write outputs, and then if it wasn't good
	// we can try again without destruction of the device data.
	int state;  // Current Input Index
};

class ParticleSimulation
{
//	static const int MAX_POINTS = 1000;
	static const int NUM_PARTICLES = 64*256;

	CurrentLineVecField *m_vectorField;

	cudaStream_t stream1;
	cudaStream_t stream2;
	cudaEvent_t eventDoneDtoHMemCpy;
	cudaEvent_t eventDoneKernel;

	SIMREAL4 *m_particlePos;
	SIMREAL3 *m_particleVel;

	DeviceData m_deviceData;

	int m_numBlocks;
public:

	static const int NUM_THREADS_PER_BLOCK = 256;

	ParticleSimulation(CurrentLineVecField *fieldFn);
	~ParticleSimulation();

 	void Draw(float pSize) const;
 	void Update(SIMREAL dt);
 	void Reset();
	void Recycle(int index);
};

// Cuda implementation
bool IntegrateNBodySystem(DeviceData &deviceData, int numBodies, float currentAmperes, float chargePerPointOnWire, float particleDensity, int outputIndex, SIMREAL dt, cudaStream_t stream);

#endif // PARTICLE_SIMULATION_H