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
	cudaHostAlloc(&m_particlePos, NUM_PARTICLES * sizeof(float4), cudaHostAllocDefault);
	m_particleVel = new float3[ NUM_PARTICLES ];

	m_deviceData.numPointsInPairs = m_vectorField->lineSegs.size();

	// we need 2 buffers for alternating writes so that we can continue to compute while reading one buffer back to the host
	cudaMalloc(&m_deviceData.particlePos[0], NUM_PARTICLES * sizeof(float4));
	cudaMalloc(&m_deviceData.particleVel[0], NUM_PARTICLES * sizeof(float3));
	cudaMalloc(&m_deviceData.particlePos[1], NUM_PARTICLES * sizeof(float4));
	cudaMalloc(&m_deviceData.particleVel[1], NUM_PARTICLES * sizeof(float3));
	cudaMalloc(&m_deviceData.integrationStepError, NUM_PARTICLES * sizeof(float));
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
	delete [] m_particleVel;
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
	cudaMemcpyAsync( m_deviceData.particlePos[0], m_particlePos, m_numBlocks * NUM_THREADS_PER_BLOCK * sizeof(float4), cudaMemcpyHostToDevice, stream2);
	cudaMemcpyAsync( m_deviceData.particleVel[0], m_particleVel, m_numBlocks * NUM_THREADS_PER_BLOCK * sizeof(float3), cudaMemcpyHostToDevice, stream2);
	cudaMemcpyAsync( m_deviceData.linePairPos, &m_vectorField->lineSegs[0], m_deviceData.numPointsInPairs * sizeof(float3), cudaMemcpyHostToDevice, stream2);
	m_deviceData.state = 0;


	cudaError err = cudaGetLastError();
	if (err) printf("memcpy h2d: %s\n", cudaGetErrorString(err));
}


void ParticleSimulation::Recycle( int index )
{
	float r1;
	float r2;
	float r3;

	do 
	{
		r1 = (float)rand() / RAND_MAX;
		r2 = (float)rand() / RAND_MAX;
		r3 = (float)rand() / RAND_MAX;

		r1 = (r1*2.0f - 1.0f);	// +/- 0.2
		r2 = (r2*2.0f - 1.0f);	// +/- 0.2
		r3 = (r3*2.0f - 1.0f);	// +/- 0.2
	} while (r1*r1 + r2*r2 + r3*r3 > 1.0f);

// 	float r4 = (float)rand() / RAND_MAX; 
// 	float r5 = (float)rand() / RAND_MAX;
// 	float r6 = (float)rand() / RAND_MAX;
// 
// 	r4 = (r4*2.0f - 1.0f) * 10.0f;	// +/- 0.2
// 	r5 = (r5*2.0f - 1.0f) * 0.0f;	// +/- 0.2
// 	r6 = (r6*2.0f - 1.0f) * 0.0f;	// +/- 0.2

	m_particlePos[index].x = -0.1f + r1*.002;
	m_particlePos[index].y =  0.0f + r2*.002;
	m_particlePos[index].z =  0.0f + r3*.002;
	m_particlePos[index].w = 9.10938188e-31f; // mass of an electron;

	m_particleVel[index].x = randFuzzy(1000000.f, 1000.0f);
	m_particleVel[index].y = randFuzzy(0.0f, 0.0f);
	m_particleVel[index].z = randFuzzy(0.0f, 0.0f);
}

extern ad::Scalar gCurrent; // current (current per strand * strands)
extern ad::Scalar gChargePerMeterOfWire;

void ParticleSimulation::Update(float dt)
{
// 	static int counter = 0;
// 	counter++;
// 	{
// 		static const int NUM_REFRESH_PER_FRAME = 100;
// 		if (counter * (NUM_REFRESH_PER_FRAME+1) >= NUM_PARTICLES ) counter = 0;
// 		for (int i = 0; i < NUM_REFRESH_PER_FRAME; i++)
// 		{
// 			int index = (counter * NUM_REFRESH_PER_FRAME + i) % NUM_PARTICLES;
// 			Recycle(index);
// 		}
// 		
// 		cudaMemcpyAsync( &m_deviceData.particlePos[0][index], &m_particlePos[index], sizeof(float4), cudaMemcpyHostToDevice, stream2);
// 		cudaMemcpyAsync( &m_deviceData.particleVel[0][index], &m_particleVel[index], sizeof(float3), cudaMemcpyHostToDevice, stream2);
// 	}

	cudaError err;
	dt *= 1.0e-9f; // nanoseconds

	int outState = m_deviceData.state == 0? 1:0;

	cudaDeviceSynchronize();

	IntegrateNBodySystem( m_deviceData, NUM_PARTICLES, (float)gCurrent, (float)gChargePerMeterOfWire, outState, dt, stream1);
	cudaMemcpyAsync( m_particlePos, m_deviceData.particlePos[m_deviceData.state], m_numBlocks * NUM_THREADS_PER_BLOCK * sizeof(float4), cudaMemcpyDeviceToHost, stream2);
	
	m_deviceData.state = outState;

	err = cudaGetLastError();
	if (err) printf("%s\n", cudaGetErrorString(err));
}

void ParticleSimulation::Draw(float pSize) const
{
	cudaStreamSynchronize(stream2);
	for (int i = 0; i < NUM_PARTICLES; i++)
	{
		float goopHalfSize = 0.01f * pSize;

		float3 &pos = *(float3 *)&m_particlePos[i];

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
}