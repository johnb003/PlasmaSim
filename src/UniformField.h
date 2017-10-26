#ifndef UNIFORM_FIELD_H
#define UNIFORM_FIELD_H

#include "VectorFunction.h"

// Returns the same vector for any position.
class UniformField : public VectorFunction
{
	const ad::Vec4 m_bFieldDir;
	const ad::Vec4 m_eFieldDir;

public:
	UniformField(const ad::Vec4 &bField, const ad::Vec4 &eField) :
		VectorFunction(), m_bFieldDir(bField), m_eFieldDir(eField) {}
	void GetField(ad::Vec4 &outBField, ad::Vec4 &outEField, const ad::Vec4 &position) const override
	{
		outBField = m_bFieldDir;
		outEField = m_eFieldDir;
	}
};

#endif // UNIFORM_FIELD_H