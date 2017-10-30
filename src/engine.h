#ifndef ENGINE_H
#define ENGINE_H

#include "adScalar.h"

#include "SDL.h"
class FieldRenderer;
class VectorFunction;
class FieldAccumulator;
class ParticleSystem;
class ParticleSimulation;
class CurrentLineVecField;
class Tweaks;
class Font;

class Engine
{
	int width, height, flags;
	float xr,yr;
	float x,y,z;
	float elapsedtime;
	float objectYaw;
	float objectPitch;
	ad::Scalar dtMult;

//	FieldRenderer *fieldRenderer;
	CurrentLineVecField *field;
	ParticleSimulation *pSim;
	Tweaks *tweakSys;
	Font *font;

	SDL_Window *sdl_window;
	SDL_GLContext sdl_gl_context;

	float zoom;

	int dx,dy,dz;
	bool tweakingUI;

	void HandleKey(SDL_Keycode key, bool down);
	void Events();
	Uint32 last_time; 
public:
	Engine();
	void Update();
	void Draw();
	~Engine();
};


#endif // ENGINE_H
