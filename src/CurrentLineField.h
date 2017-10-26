#ifndef CURRENT_LINE_FIELD_H
#define CURRENT_LINE_FIELD_H

#include "VectorFunction.h"

// This draws a wire segment and a little box that indicates the charge movement on the wire.
// It also can compute the B-field at a position due to this wire segment using the VectorFunction interface.
class CurrentLineField : public VectorFunction
{
	const ad::Vec4 m_pointA;
	const ad::Vec4 m_pointB;

	ad::Scalar tPos;
public:
	CurrentLineField(const ad::Vec4 &pointA, const ad::Vec4 &pointB);

	void GetField(ad::Vec4 &outBField, ad::Vec4 &outEField, const ad::Vec4 &position) const;

	void DrawField() const;
	void UpdateField(ad::Scalar dt);
};

extern ad::Scalar gCurrent;  // current (current per strand * strands)
extern ad::Scalar gChargePerMeterOfWire;

#endif // CURRENT_LINE_FIELD_H