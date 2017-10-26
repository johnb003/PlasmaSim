#ifndef VECTOR_FUNCTION_H
#define VECTOR_FUNCTION_H

#include "adVector.h"

// Abstract base class that allows polling a field
class VectorFunction
{
public:
	VectorFunction() { }
	virtual ~VectorFunction() {}

	virtual void GetField(ad::Vec4 &outBField, ad::Vec4 &outEField, const ad::Vec4 &position) const = 0;

	virtual void DrawField() const {}
	virtual void UpdateField(ad::Scalar dt) {}
};

#endif // VECTOR_FUNCTION_H