#include "FieldAccumulator.h"
#include "adVector.h"

using ad::Vec4;
using std::vector;

FieldAccumulator::FieldAccumulator()
	: VectorFunction()
{
}

FieldAccumulator::~FieldAccumulator()
{
	vector<VectorFunction *>::iterator it;
	for (it = fields.begin(); it != fields.end(); ++it)
	{
		delete *it;
	}
}

void FieldAccumulator::GetField(ad::Vec4 &outBField, ad::Vec4 &outEField, const ad::Vec4 &position) const
{
	outBField = outEField = Vec4::m_Zero;
	
	vector<VectorFunction *>::const_iterator it;
	for (it = fields.begin(); it != fields.end(); ++it)
	{
		Vec4 resultE;
		Vec4 resultB;
		(*it)->GetField(resultB, resultE, position);
		outBField += resultB;
		outEField += resultE;
	}
}

void FieldAccumulator::AddField(VectorFunction *vf)
{
	fields.push_back(vf);
}

void FieldAccumulator::DrawField() const
{
	vector<VectorFunction *>::const_iterator it;
	for (it = fields.begin(); it != fields.end(); ++it)
	{
		(*it)->DrawField();
	}
}

void FieldAccumulator::UpdateField(ad::Scalar dt)
{
	vector<VectorFunction *>::iterator it;
	for (it = fields.begin(); it != fields.end(); ++it)
	{
		(*it)->UpdateField(dt);
	}
}