#include "CurrentLineVectorField.h"

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


CurrentLineVecField::CurrentLineVecField()
{
	m_tPos = 0.0f;
}

void CurrentLineVecField::AddField(float ax, float ay, float az, float bx, float by, float bz)
{
	lineSegs.push_back(float3());
	lineSegs.back().x = ax;
	lineSegs.back().y = ay;
	lineSegs.back().z = az;
	lineSegs.push_back(float3());
	lineSegs.back().x = bx;
	lineSegs.back().y = by;
	lineSegs.back().z = bz;
}

void CurrentLineVecField::DrawField() const
{
	glColor4f(1,1,0,1);
	glBegin(GL_LINES);  
	for (int i = 0, n = lineSegs.size(); i < n; i+=2)
	{
		glVertex3d(lineSegs[i].x, lineSegs[i].y, lineSegs[i].z);
		glVertex3d(lineSegs[i+1].x, lineSegs[i+1].y, lineSegs[i+1].z);
	}
	glEnd();

	// Draw little globs that animate around the wire, to help visualize the current direction
	for (int i = 0, n = lineSegs.size(); i < n; i+=2)
	{
		float goopX = lineSegs[i].x + (lineSegs[i+1].x - lineSegs[i].x) * m_tPos;
		float goopY = lineSegs[i].y + (lineSegs[i+1].y - lineSegs[i].y) * m_tPos;
		float goopZ = lineSegs[i].z + (lineSegs[i+1].z - lineSegs[i].z) * m_tPos;

		float goopHalfSize = 0.001f;

		glBegin(GL_QUAD_STRIP);  
		glVertex3d(goopX-goopHalfSize, goopY+goopHalfSize, goopZ-goopHalfSize);
		glVertex3d(goopX-goopHalfSize, goopY+goopHalfSize, goopZ+goopHalfSize);

		glVertex3d(goopX-goopHalfSize, goopY-goopHalfSize, goopZ-goopHalfSize);
		glVertex3d(goopX-goopHalfSize, goopY-goopHalfSize, goopZ+goopHalfSize);

		glVertex3d(goopX+goopHalfSize, goopY-goopHalfSize, goopZ-goopHalfSize);
		glVertex3d(goopX+goopHalfSize, goopY-goopHalfSize, goopZ+goopHalfSize);

		glVertex3d(goopX+goopHalfSize, goopY+goopHalfSize, goopZ-goopHalfSize);
		glVertex3d(goopX+goopHalfSize, goopY+goopHalfSize, goopZ+goopHalfSize);

		glVertex3d(goopX-goopHalfSize, goopY+goopHalfSize, goopZ-goopHalfSize);
		glVertex3d(goopX-goopHalfSize, goopY+goopHalfSize, goopZ+goopHalfSize);
		glEnd();

		glBegin(GL_QUADS);  
		glVertex3d(goopX-goopHalfSize, goopY+goopHalfSize, goopZ+goopHalfSize);
		glVertex3d(goopX+goopHalfSize, goopY+goopHalfSize, goopZ+goopHalfSize);
		glVertex3d(goopX+goopHalfSize, goopY-goopHalfSize, goopZ+goopHalfSize);
		glVertex3d(goopX-goopHalfSize, goopY-goopHalfSize, goopZ+goopHalfSize);

		glVertex3d(goopX-goopHalfSize, goopY-goopHalfSize, goopZ-goopHalfSize);
		glVertex3d(goopX+goopHalfSize, goopY-goopHalfSize, goopZ-goopHalfSize);
		glVertex3d(goopX+goopHalfSize, goopY+goopHalfSize, goopZ-goopHalfSize);
		glVertex3d(goopX-goopHalfSize, goopY+goopHalfSize, goopZ-goopHalfSize);
		glEnd();
	}

}

void CurrentLineVecField::UpdateField(float dt)
{
	m_tPos += dt;
	if (m_tPos > 1.0f)
		m_tPos = 0.0f;
}