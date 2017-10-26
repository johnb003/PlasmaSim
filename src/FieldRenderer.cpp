#include "FieldRenderer.h"
#include "VectorFunction.h"
#include "adVector.h"

#include "assert.h"
#include "math.h"

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


using ad::Vec4;
using ad::Scalar;
using ad::min;
using ad::max;

FieldRenderer::FieldRenderer(const VectorFunction *fieldFunction)
	: m_FieldFunction(fieldFunction)
{
	numSteps = 100;
	numStepsEuler = 10;
	resetPoint = 0;
	GenerateFieldLinesFromGrid();
}

FieldRenderer::~FieldRenderer()
{
}



void FieldRenderer::GenerateFieldLinesFromGrid()
{
	const int NUM_TRACERS_PER_AXIS = 5;
}

void FieldRenderer::DrawField() const
{
	Scalar minStr = 0.0;
	Scalar maxStr = 100.0;

	glBegin(GL_LINES);  
	glColor4f(1,0,0,1);
	glVertex3d(m_InteractiveTracerStart.x-0.01f, m_InteractiveTracerStart.y, m_InteractiveTracerStart.z);
	glVertex3d(m_InteractiveTracerStart.x+0.01f, m_InteractiveTracerStart.y, m_InteractiveTracerStart.z);
	glColor4f(0,1,0,1);
	glVertex3d(m_InteractiveTracerStart.x, m_InteractiveTracerStart.y-0.01f, m_InteractiveTracerStart.z);
	glVertex3d(m_InteractiveTracerStart.x, m_InteractiveTracerStart.y+0.01f, m_InteractiveTracerStart.z);
	glColor4f(0,0,1,1);
	glVertex3d(m_InteractiveTracerStart.x, m_InteractiveTracerStart.y, m_InteractiveTracerStart.z-0.01f);
	glVertex3d(m_InteractiveTracerStart.x, m_InteractiveTracerStart.y, m_InteractiveTracerStart.z+0.01f);
	glEnd();

	// tracer lines
	glBegin(GL_LINES);  
	for (int i = 0; i < numPoints; i++)
	{
		Scalar strength = (max(minStr, min(points[i].w, maxStr)) - minStr)/(maxStr-minStr);

		glColor4d(0.4+strength*0.6f,0.6f-strength*0.6f,0.6f-strength*0.6f,1);
		glVertex3d(points[i].x, points[i].y, points[i].z);
	}
	glEnd();

	// field vectors
// 	glBegin(GL_LINES);  
// 	for (int i = 0; i < numPoints; i+=2)
// 	{
// 		Scalar strength = (max(minStr, min(points[i].w, maxStr)) - minStr)/(maxStr-minStr);
// 
// 		glColor4d(0.4+strength*0.6f,0.6f-strength*0.6f,0.6f-strength*0.6f,1);
// 		glVertex3f(points[i].x, points[i].y, points[i].z);
// 		Vec4 outField;
// 		Vec4 field = m_FieldFunction->GetField(outField, points[i]);
// 		field /= 80;
// 		glVertex3f(points[i].x+field.x, points[i].y+field.y, points[i].z+field.z);
// 	}
// 	glEnd();
}

void FieldRenderer::IncNumSteps()
{
	numSteps *=2;
	numPoints = resetPoint;
	ExtendTracer(m_InteractiveTracerStart, 1.0);
	ExtendTracer(m_InteractiveTracerStart, -1.0);
}

void FieldRenderer::DecNumSteps()
{
	numSteps /=2;
	numPoints = resetPoint;
	ExtendTracer(m_InteractiveTracerStart, 1.0);
	ExtendTracer(m_InteractiveTracerStart, -1.0);
}

void FieldRenderer::Save()
{
	resetPoint = numPoints;
}

void FieldRenderer::NudgeInteractive(Scalar dx, Scalar dy, Scalar dz)
{
	m_InteractiveTracerStart.x += dx;
	m_InteractiveTracerStart.y += dy;
	m_InteractiveTracerStart.z += dz;

	numPoints = resetPoint;
	ExtendTracer(m_InteractiveTracerStart, 1.0);
	ExtendTracer(m_InteractiveTracerStart, -1.0);
}

Vec4 &FieldRenderer::RKStep(const Vec4 &startPos, Scalar stepSize, const VectorFunction &fn, Vec4 &rk_Sum)
{
	assert(&startPos != &rk_Sum);
	Vec4 rk_k1B;
	Vec4 rk_k2B;
	Vec4 rk_k3B;
	Vec4 rk_k4B;
	Vec4 rk_k1E;
	Vec4 rk_k2E;
	Vec4 rk_k3E;
	Vec4 rk_k4E;

	fn.GetField(rk_k1B, rk_k1E, startPos);
	rk_Sum = startPos + 0.5f * stepSize * rk_k1B;

	fn.GetField(rk_k2B, rk_k2E, rk_Sum);
	rk_Sum = startPos + 0.5f * stepSize * rk_k2B;

	fn.GetField(rk_k3B, rk_k3E, rk_Sum);
	rk_Sum = startPos +        stepSize * rk_k3B;

	fn.GetField(rk_k4B, rk_k4E, rk_Sum);
	rk_Sum = startPos + stepSize * (rk_k1B + 2*(rk_k2B + rk_k3B) + rk_k4B)/6;
	return rk_Sum;
}

void FieldRenderer::ExtendTracer(const Vec4 &startLocation, Scalar dir)
{
	static const Scalar MIN_STEP_SIZE = 0.0000000001;
	static const Scalar MAX_STEP_SIZE = 1.0;

	Scalar stepSize = 1.0f;
	Scalar ALLOWED_ERROR_THRESH = 0.000001;
	int segments = 0;

 	Vec4 pos = startLocation;

	while (numPoints+16 < MAX_POINTS && segments < numSteps)
	{
		Vec4 fieldStrB;
		Vec4 fieldStrE;
		m_FieldFunction->GetField(fieldStrB, fieldStrE, pos);
		Scalar startStr = fieldStrB.Length3();

		Scalar minStr = 0.0;
		Scalar maxStr = 100.0;
		Scalar strength = (max(minStr, min(startStr, maxStr)) - minStr)/(maxStr-minStr);

		strength = 1/(1 + strength*10);

		Vec4 rkFullStep;
		RKStep(pos, stepSize*dir, *m_FieldFunction, rkFullStep);

		Vec4 rkHalfStep1;
		Vec4 rkHalfStep2;
		RKStep(rkFullStep, stepSize*-dir/2, *m_FieldFunction, rkHalfStep1);
		RKStep(rkHalfStep1, stepSize*-dir/2, *m_FieldFunction, rkHalfStep2);

		// estimate the local error
		Vec4 errorDiff = rkHalfStep2 - pos;
		Scalar errorSum = fabs(errorDiff.x) + fabs(errorDiff.y) + fabs(errorDiff.z);

		if (errorSum > ALLOWED_ERROR_THRESH*strength && stepSize > MIN_STEP_SIZE)
		{
			// try a smaller stepsize
			Scalar offByPowersOfTen = max(1.0, log10(errorSum/(ALLOWED_ERROR_THRESH*strength)));

			Scalar reductionFactor = pow(0.6, offByPowersOfTen);

			stepSize *= reductionFactor;
			if (stepSize < MIN_STEP_SIZE)
				stepSize = MIN_STEP_SIZE;

			// retry
			continue;
		}
		else if (errorSum < (ALLOWED_ERROR_THRESH * strength * 0.01))
		{
// 			if (errorSum < DBL_MIN)
// 				errorSum = DBL_MIN;
// 			Scalar offByPowersOfTen = max(1.0, log10(ALLOWED_ERROR_THRESH*strength/errorSum));
// 			Scalar augmentationFactor = pow(1.66667, offByPowersOfTen);
			// next time we can bite off a bit more
			stepSize *= 1.66667;
			if (stepSize > MAX_STEP_SIZE)
				stepSize = MAX_STEP_SIZE;
		}

		m_FieldFunction->GetField(fieldStrB, fieldStrE, rkFullStep);
		rkFullStep.w = fieldStrB.Length3();
		pos.w = startStr;
		points[numPoints] = pos;
		points[numPoints+1] = rkFullStep;
		numPoints+=2;

// 		// mirror x
// 		points[numPoints] = pos;
// 		points[numPoints].x = -points[numPoints].x;
// 		points[numPoints+1] = rkFullStep;
// 		points[numPoints+1].x = -points[numPoints+1].x;
// 		numPoints+=2;
// 
// 		// Mirror Y
// 		points[numPoints] = pos;
// 		points[numPoints].y = -points[numPoints].y;
// 		points[numPoints+1] = rkFullStep;
// 		points[numPoints+1].y = -points[numPoints+1].y;
// 		numPoints+=2;
// 
// 		// mirror x
// 		points[numPoints] = pos;
// 		points[numPoints].y = -points[numPoints].y;
// 		points[numPoints].x = -points[numPoints].x;
// 		points[numPoints+1] = rkFullStep;
// 		points[numPoints+1].y = -points[numPoints+1].y;
// 		points[numPoints+1].x = -points[numPoints+1].x;
// 		numPoints+=2;

		// Mirror Z

// 		points[numPoints] = pos;
// 		points[numPoints].z = -points[numPoints].z;
// 		points[numPoints+1] = rkHalfStep2;
// 		points[numPoints+1].z = -points[numPoints+1].z;
// 		numPoints+=2;
// 
// 		// mirror x
// 		points[numPoints] = pos;
// 		points[numPoints].z = -points[numPoints].z;
// 		points[numPoints].x = -points[numPoints].x;
// 		points[numPoints+1] = rkFullStep;
// 		points[numPoints+1].z = -points[numPoints+1].z;
// 		points[numPoints+1].x = -points[numPoints+1].x;
// 		numPoints+=2;
// 
// 		// Mirror Y
// 		points[numPoints] = pos;
// 		points[numPoints].z = -points[numPoints].z;
// 		points[numPoints].y = -points[numPoints].y;
// 		points[numPoints+1] = rkFullStep;
// 		points[numPoints+1].z = -points[numPoints+1].z;
// 		points[numPoints+1].y = -points[numPoints+1].y;
// 		numPoints+=2;
// 
// 		// mirror x
// 		points[numPoints] = pos;
// 		points[numPoints].z = -points[numPoints].z;
// 		points[numPoints].y = -points[numPoints].y;
// 		points[numPoints].x = -points[numPoints].x;
// 		points[numPoints+1] = rkFullStep;
// 		points[numPoints+1].z = -points[numPoints+1].z;
// 		points[numPoints+1].y = -points[numPoints+1].y;
// 		points[numPoints+1].x = -points[numPoints+1].x;
//		numPoints+=2;

		pos = rkFullStep;

		segments++;
	}
}
