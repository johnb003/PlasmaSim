#ifndef FIELD_ACCUMULATOR_H
#define FIELD_ACCUMULATOR_H

#include "VectorFunction.h"
#include <vector>

// Implements GetField, which returns the combined result of other `VectorFunction`s.
class FieldAccumulator : public VectorFunction
{
	std::vector<VectorFunction *> fields;

public:
	FieldAccumulator();
	virtual ~FieldAccumulator();

	virtual void GetField(ad::Vec4 &outBField, ad::Vec4 &outEField, const ad::Vec4 &position) const;
	virtual void DrawField() const;
	virtual void UpdateField(ad::Scalar dt);

	void AddField(VectorFunction *fn);
};

#endif // FIELD_ACCUMULATOR_H