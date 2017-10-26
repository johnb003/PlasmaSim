#include "font.h"
#include "sdl_gl_texture.h"

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

#include <string.h>

Font::Font(const char *filename, int size) : size(size)
{
	f_texture = LoadTexture(filename);
	if (!f_texture) return;

	// Creating 256 Display Lists
	unsigned int base = glGenLists(257);

	glNewList(9, GL_COMPILE);
		glTranslatef(X_SPACE_SIZE * size * 4, 0, 0);
	glEndList();

	float cx;	// Holds Our X Character Coord
	float cy;	// Holds Our Y Character Coord
	for (int loop = 0; loop < 256; loop++)
	{
		cx = (loop % 16) / 16.0f;            // X Position Of Current Character
		cy = (loop / 16) / 16.0f + 0.001f;   // Y Position Of Current Character + correction

		glNewList(32 + loop, GL_COMPILE);
			glBegin(GL_QUADS);
				glTexCoord2f(cx,           cy);           glVertex2d(0, 0);
				glTexCoord2f(cx,           cy + 0.0625f); glVertex2d(0, size);
				glTexCoord2f(cx + 0.0625f, cy + 0.0625f); glVertex2d(size, size);
				glTexCoord2f(cx + 0.0625f, cy);           glVertex2d(size, 0);
			glEnd();

			glTranslatef(X_SPACE_SIZE * size, 0, 0);
		glEndList();
	}
}

void Font::print(int x, int y, const char *text)
{
	glPushAttrib(GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glBindTexture(GL_TEXTURE_2D, f_texture);

	glColor4f(1, 1, 1, 1);
	glPushMatrix();
	glTranslated(x, y, 0);
	glCallLists(strlen(text), GL_BYTE, text);
	glPopMatrix();
	glPopAttrib();
};