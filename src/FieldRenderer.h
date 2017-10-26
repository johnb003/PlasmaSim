#ifndef FIELD_RENDERER_H
#define FIELD_RENDERER_H

#include "adVector.h"

class VectorFunction;

class FieldRenderer
{
	const VectorFunction *m_FieldFunction;

	// this ought to be kind of temporary
	// a more structured approach maybe would be better for dynamic data
	static const int MAX_POINTS = 1000000;
	ad::Vec4 points[MAX_POINTS];
	int numPoints;
	int resetPoint;

	int numSteps;
	int numStepsEuler;

	ad::Vec4 m_InteractiveTracerStart;

	ad::Vec4 &RKStep(const ad::Vec4 &startPos, ad::Scalar stepSize, const VectorFunction &fn, ad::Vec4 &endPos);
	void ExtendTracer(const ad::Vec4 &startLocation, ad::Scalar dir);
	void GenerateFieldLinesFromGrid();

public:
	FieldRenderer(const VectorFunction *fieldFunction);
	~FieldRenderer();

	void IncNumSteps();
	void DecNumSteps();

	void Save();
	void DrawField() const;
	void NudgeInteractive(ad::Scalar dx, ad::Scalar dy, ad::Scalar dz);

};

#endif // FIELD_RENDERER_H


