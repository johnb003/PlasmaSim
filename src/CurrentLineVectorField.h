#ifndef CURRENT_LINE_VECTOR_FIELD_H
#define CURRENT_LINE_VECTOR_FIELD_H

#include <vector>
#include <cuda_runtime_api.h>

class CurrentLineVecField
{
	float m_tPos;
public:
	CurrentLineVecField();

	std::vector<float3> lineSegs;

	void AddField(float ax, float ay, float az, float bx, float by, float bz);
	void DrawField() const;
	void UpdateField(float dt);
};

#endif // CURRENT_LINE_VECTOR_FIELD_H