#ifndef SDL_GL_TEXTURE_H
#define SDL_GL_TEXTURE_H

// Loads a texture from file using SDL to parse the image
// Then loads the image into texture memory and returns the GL texture index.
// Returns: GL texture index or 0 if it could not load.
// Note: Zero is a reserved texture name, so it's only returned when the load was invalid.
unsigned int LoadTexture(const char *filename, bool mipmapped = false);

#endif  // SDL_GL_TEXTURE_H
