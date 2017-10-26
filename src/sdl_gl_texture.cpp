#include "sdl_gl_texture.h"

#include <SDL.h>

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

#include <stdio.h>

unsigned int LoadTexture(const char *filename, bool mipmapped)
{
	SDL_Surface *surface;
	surface = SDL_LoadBMP(filename);
	if (surface == NULL) {
		fprintf(stderr, "Failed to load texture: %s", filename);
		return 0;
	}
	
	int fmt=0;
	int intFmt=0;
	switch (surface->format->BytesPerPixel)
	{
	case 1:
		fmt = GL_ALPHA;
		intFmt = GL_ALPHA;
		break;
	case 3:
		fmt = GL_BGR_EXT;
		intFmt = GL_RGB;
		break;
	case 4:
		fmt = GL_BGRA_EXT;
		intFmt = GL_RGBA;
		break;
	}

	unsigned int textureIndex = 0;
	glGenTextures(1, &textureIndex);
	glBindTexture(GL_TEXTURE_2D, textureIndex);

	if (mipmapped)
	{
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_NEAREST);
		gluBuild2DMipmaps(GL_TEXTURE_2D, intFmt, surface->w, surface->h, fmt, GL_UNSIGNED_BYTE, surface->pixels);
	}
	else
	{
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D, 0, intFmt, surface->w, surface->h, 0, fmt, GL_UNSIGNED_BYTE, surface->pixels);		
	}
	
	SDL_FreeSurface(surface);
	return textureIndex;
}
