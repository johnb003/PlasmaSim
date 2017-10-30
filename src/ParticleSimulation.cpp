#include "ParticleSimulation.h"
#include "CurrentLineVectorField.h"
#include "adScalar.h"

#include "cuda_runtime.h"

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



ParticleSimulation::ParticleSimulation(CurrentLineVecField *vectorField)
	: m_vectorField(vectorField)
{
	// allocate host memory for the particles

	// allocate host and device memory for the particles
	m_numBlocks = (NUM_PARTICLES-1) / NUM_THREADS_PER_BLOCK + 1;

	// non-paged memory for faster copy back
	cudaHostAlloc(&m_particlePos, NUM_PARTICLES * sizeof(SIMREAL4), cudaHostAllocDefault);
	cudaHostAlloc(&m_particleVel, NUM_PARTICLES * sizeof(SIMREAL3), cudaHostAllocDefault);

	m_deviceData.numPointsInPairs = m_vectorField->lineSegs.size();

	// we need 2 buffers for alternating writes so that we can continue to compute while reading one buffer back to the host
	cudaMalloc(&m_deviceData.particlePos[0], NUM_PARTICLES * sizeof(SIMREAL4));
	cudaMalloc(&m_deviceData.particleVel[0], NUM_PARTICLES * sizeof(SIMREAL3));
	cudaMalloc(&m_deviceData.particlePos[1], NUM_PARTICLES * sizeof(SIMREAL4));
	cudaMalloc(&m_deviceData.particleVel[1], NUM_PARTICLES * sizeof(SIMREAL3));
	cudaMalloc(&m_deviceData.integrationStepError, NUM_PARTICLES * sizeof(SIMREAL));
	cudaMalloc(&m_deviceData.linePairPos, m_deviceData.numPointsInPairs * sizeof(float3));
	
	cudaStreamCreate(&stream1);
	cudaStreamCreate(&stream2);

	cudaEventCreateWithFlags(&eventDoneDtoHMemCpy, cudaEventDisableTiming);
	cudaEventCreateWithFlags(&eventDoneKernel, cudaEventDisableTiming);

	cudaError err = cudaGetLastError();
	if (err) printf("%s\n", cudaGetErrorString(err));

	// populate the host data
	Reset();
}

ParticleSimulation::~ParticleSimulation()
{
	cudaFree(m_deviceData.linePairPos);
	cudaFree(m_deviceData.integrationStepError);
	cudaFree(m_deviceData.particleVel[1]);
	cudaFree(m_deviceData.particlePos[1]);
	cudaFree(m_deviceData.particleVel[0]);
	cudaFree(m_deviceData.particlePos[0]);

	cudaEventDestroy(eventDoneKernel);
	cudaEventDestroy(eventDoneDtoHMemCpy);
	cudaStreamDestroy(stream2);
	cudaStreamDestroy(stream1);

	cudaFreeHost(m_particlePos);
	cudaFreeHost(m_particleVel);
}

float randNormal(float from, float to)
{
	return from + (to-from)*((float)rand() / RAND_MAX);
}

float randSymmetricNormal()
{
	return randNormal(-1.0f, 1.0f);
}

float randFuzzy(float center, float plusOrMinus)
{
	return center + randNormal(-plusOrMinus, plusOrMinus);
}

void ParticleSimulation::Reset()
{
	for (int i = 0; i < NUM_PARTICLES; ++i)
	{
		Recycle(i);
	}

	// copy the host data to the device memory
	cudaMemcpyAsync( m_deviceData.particlePos[0], m_particlePos, m_numBlocks * NUM_THREADS_PER_BLOCK * sizeof(SIMREAL4), cudaMemcpyHostToDevice, stream2);
	cudaMemcpyAsync( m_deviceData.particleVel[0], m_particleVel, m_numBlocks * NUM_THREADS_PER_BLOCK * sizeof(SIMREAL3), cudaMemcpyHostToDevice, stream2);
	cudaMemcpyAsync( m_deviceData.linePairPos, &m_vectorField->lineSegs[0], m_deviceData.numPointsInPairs * sizeof(float3), cudaMemcpyHostToDevice, stream2);
	m_deviceData.state = 0;


	cudaError err = cudaGetLastError();
	if (err) printf("memcpy h2d: %s\n", cudaGetErrorString(err));
}

float DT_SPREAD = 0.01f;

SIMREAL timeEvolution = 0;

void ParticleSimulation::Recycle( int index )
{
	// Get a random vector within a unit sphere
	float r1;
	float r2;
	float r3;
	do
	{
		r1 = randSymmetricNormal();
		r2 = randSymmetricNormal();
		r3 = randSymmetricNormal();
	} while (r1*r1 + r2*r2 + r3*r3 > 1.0f);

	m_particlePos[index].x = -0.195;// + r1*.002;
	m_particlePos[index].y =  0.0;// + r2*.002;
	m_particlePos[index].z =  0.0;// + r3*.002;
	m_particlePos[index].w = randFuzzy(1.0f, DT_SPREAD);//9.10938188e-31f; // mass of an electron;

	m_particleVel[index].x = randFuzzy(0.0f, 0.0f);
	m_particleVel[index].y = randFuzzy(0.0f, 0.0f);
	m_particleVel[index].z = randFuzzy(0.0f, 0.0f);
	timeEvolution = 0;
}

extern ad::Scalar gCurrent; // current (current per strand * strands)
extern ad::Scalar gChargePerMeterOfWire;

ad::Scalar gElectronDensityPerParticle = 1.0;

void ParticleSimulation::Update(SIMREAL dt)
{
	// We're counting real-time as nanoseconds here!
	dt *= 1.0e-9f;

	cudaError err;

	int outState = m_deviceData.state == 0? 1:0;

	cudaDeviceSynchronize();

	bool updated = IntegrateNBodySystem( m_deviceData, NUM_PARTICLES, (float)gCurrent, (float)gChargePerMeterOfWire, (float)gElectronDensityPerParticle, outState, dt, stream1);
	if (updated) {
		cudaMemcpyAsync( m_particlePos, m_deviceData.particlePos[m_deviceData.state], m_numBlocks * NUM_THREADS_PER_BLOCK * sizeof(SIMREAL4), cudaMemcpyDeviceToHost, stream2);
		m_deviceData.state = outState;
	}

	err = cudaGetLastError();
	if (err) printf("%s\n", cudaGetErrorString(err));
}

void ParticleSimulation::Draw(float pSize) const
{
	cudaStreamSynchronize(stream2);
	for (int i = 0; i < NUM_PARTICLES; i++)
	{
		float goopHalfSize = 0.01f * pSize;

		SIMREAL4 &pos = *(SIMREAL4 *)&m_particlePos[i];

		SIMREAL colorrange = (pos.w - 1) / (DT_SPREAD);
		colorrange /= 2.0F;
		colorrange += 0.5f;

		glColor4f(1.0f * colorrange, 1.0f * (1.0f - colorrange), 1.0f * (1.0f - colorrange), 1.0f);
		glBegin(GL_QUAD_STRIP);  
		glVertex3d(pos.x - goopHalfSize, pos.y + goopHalfSize, pos.z - goopHalfSize);
		glVertex3d(pos.x - goopHalfSize, pos.y + goopHalfSize, pos.z + goopHalfSize);

		glVertex3d(pos.x - goopHalfSize, pos.y - goopHalfSize, pos.z - goopHalfSize);
		glVertex3d(pos.x - goopHalfSize, pos.y - goopHalfSize, pos.z + goopHalfSize);

		glVertex3d(pos.x + goopHalfSize, pos.y - goopHalfSize, pos.z - goopHalfSize);
		glVertex3d(pos.x + goopHalfSize, pos.y - goopHalfSize, pos.z + goopHalfSize);

		glVertex3d(pos.x + goopHalfSize, pos.y + goopHalfSize, pos.z - goopHalfSize);
		glVertex3d(pos.x + goopHalfSize, pos.y + goopHalfSize, pos.z + goopHalfSize);

		glVertex3d(pos.x - goopHalfSize, pos.y + goopHalfSize, pos.z - goopHalfSize);
		glVertex3d(pos.x - goopHalfSize, pos.y + goopHalfSize, pos.z + goopHalfSize);
		glEnd();

		glBegin(GL_QUADS);  
		glVertex3d(pos.x - goopHalfSize, pos.y + goopHalfSize, pos.z + goopHalfSize);
		glVertex3d(pos.x + goopHalfSize, pos.y + goopHalfSize, pos.z + goopHalfSize);
		glVertex3d(pos.x + goopHalfSize, pos.y - goopHalfSize, pos.z + goopHalfSize);
		glVertex3d(pos.x - goopHalfSize, pos.y - goopHalfSize, pos.z + goopHalfSize);

		glVertex3d(pos.x - goopHalfSize, pos.y - goopHalfSize, pos.z - goopHalfSize);
		glVertex3d(pos.x + goopHalfSize, pos.y - goopHalfSize, pos.z - goopHalfSize);
		glVertex3d(pos.x + goopHalfSize, pos.y + goopHalfSize, pos.z - goopHalfSize);
		glVertex3d(pos.x - goopHalfSize, pos.y + goopHalfSize, pos.z - goopHalfSize);
		glEnd();
	}
}