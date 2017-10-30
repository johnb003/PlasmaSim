#include "engine.h"
#include "UniformField.h"
#include "CurrentLineField.h"
#include "FieldRenderer.h"
#include "FieldAccumulator.h"
#include "CurrentLineVectorField.h"
#include "ParticleSystem.h"
#include "ParticleSimulation.h"
#include "tweaks.h"
#include "font.h"

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


#include <adMatrix.h>
#define _USE_MATH_DEFINES 1
#include <math.h>
#include <iostream>

//#include "test.h"

using namespace std;
using ad::Vec4;

#define PI_OVER_180 0.01745329

//static float l0_position[] = {-1.0f, -0.3f, 1.0f, 0.0f};
static float l0_position[] = {0.0f, 0.0f, 1.0f, 0.0f};
static float l0_ambient[] =  {0.2f, 0.2f, 0.2f, 1.0f};
extern ad::Scalar gElectronDensityPerParticle;

Engine::Engine()
{
	xr = 0;
	yr = 0;
	x = 0.0f;
	y = 0.0f;
	z = 0.75f;
	dx = 0;
	dy = 0;
	dz = 0;
	objectYaw = 0.0f;
	objectPitch = 0.0f;
	tweakingUI = false;

	dtMult = 1.0;

	zoom = 1.0f;

	cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);

	if ( SDL_Init(SDL_INIT_VIDEO) < 0 )
	{
		cerr << "Unable to initialize SDL: " << SDL_GetError() << endl;
		throw(1);
	}

	width = 1280;
	height = 1024;
	flags = SDL_WINDOW_OPENGL;

	SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );

	sdl_window = SDL_CreateWindow("Field Sim",
		SDL_WINDOWPOS_CENTERED,
		SDL_WINDOWPOS_CENTERED,
		width, height, flags);

	if (sdl_window == NULL)
	{
		cerr << "Video mode set failed: " << SDL_GetError() << endl;
		throw(1);
	}

	sdl_gl_context = SDL_GL_CreateContext(sdl_window);
    
	glClearColor(0, 0, 0, 0);
	glClearDepth(1.0f);
	glShadeModel(GL_SMOOTH);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glDrawBuffer(GL_BACK);
	glEnable(GL_LIGHTING);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glEnable(GL_COLOR_MATERIAL);
//	glEnable(GL_POLYGON_SMOOTH);  // This enables shitty Anti-aliasing, and doesn't work correctly on Mac.
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_DEPTH_TEST);

	glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1.0f);
	glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05f);
	glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.01f);
	glEnable(GL_LIGHT0);
	
	glViewport(0,0,width,height);
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity( );
	gluPerspective( 60.0, (float)width/(float)height, 0.0001, 1024.0 );
	
	font = new Font("images/font.bmp", 16);

	tweakSys = &gTweaks;
	gTweaks.Init(width, height, font);
 	tweakSys->AddTweak("Dt Multiplier", &dtMult, 0.00001, 100000, Tweaks::ScaleStyle_Exponential); 
 	tweakSys->AddTweak("Current", &gCurrent, 0.000001, 1e20, Tweaks::ScaleStyle_Exponential); 
 	tweakSys->AddTweak("Charge", &gChargePerMeterOfWire, 1e-40, 1e10, Tweaks::ScaleStyle_Exponential); 
 	tweakSys->AddTweak("ParticleDensity", &gElectronDensityPerParticle, 1, 1e40, Tweaks::ScaleStyle_Exponential);


	field = new CurrentLineVecField();
	float rad = 0.2f;
	// Z planes
	field->AddField( 0,			rad,		-rad-0.01f,	-rad,		0,			-rad-0.01f );
	field->AddField(-rad,		0,			-rad-0.01f,	0,			-rad,		-rad-0.01f );
	field->AddField( 0,			-rad,		-rad-0.01f,	rad,		0,			-rad-0.01f );
	field->AddField( rad,		0,			-rad-0.01f,	0,			rad,		-rad-0.01f );
	
	field->AddField( 0,			rad,		rad+0.01f,	rad,		0,			rad+0.01f );
	field->AddField( rad,		0,			rad+0.01f,	0,			-rad,		rad+0.01f );
	field->AddField( 0,			-rad,		rad+0.01f,	-rad,		0,			rad+0.01f );
	field->AddField(-rad,		0,			rad+0.01f,	0,			rad,		rad+0.01f );

	// X Planes
	field->AddField(-rad-0.01f,	0,			rad,		-rad-0.01f,	-rad,		0 );
	field->AddField(-rad-0.01f,	-rad,		0,			-rad-0.01f,	0,			-rad );
	field->AddField(-rad-0.01f,	0,			-rad,		-rad-0.01f,	rad,		0 );
	field->AddField(-rad-0.01f,	rad,		0,			-rad-0.01f,	0,			rad );

	field->AddField( rad+0.01f,	0,			rad,		 rad+0.01f,	rad,		0 );
	field->AddField( rad+0.01f,	rad,		0,			 rad+0.01f,	0,			-rad );
	field->AddField( rad+0.01f,	0,			-rad,		 rad+0.01f,	-rad,		0 );
	field->AddField( rad+0.01f,	-rad,		0,			 rad+0.01f,	0,			rad );

	// Y Planes
	field->AddField( 0,			-rad-0.01f,	-rad,		-rad,		-rad-0.01f,	0 );
	field->AddField(-rad,		-rad-0.01f,	0,			 0,			-rad-0.01f,	rad );
	field->AddField( 0,			-rad-0.01f,	rad,		 rad,		-rad-0.01f,	0 );
	field->AddField( rad,		-rad-0.01f,	0,			 0,			-rad-0.01f,	-rad );

	field->AddField( 0,			rad+0.01f,	-rad,		 rad,		rad+0.01f,	0 );
	field->AddField( rad,		rad+0.01f,	0,			 0,			rad+0.01f,	rad );
	field->AddField( 0,			rad+0.01f,	rad,		-rad,		rad+0.01f,	0 );
	field->AddField(-rad,		rad+0.01f,	0,			 0,			rad+0.01f,	-rad );

//	pSys = new ParticleSystem(fieldAccum);
	pSim = new ParticleSimulation(field);
//	fieldRenderer = new FieldRenderer(field);
	
	last_time = SDL_GetTicks();
//    SDL_WM_GrabInput(SDL_GRAB_OFF);
    SDL_ShowCursor(SDL_ENABLE);
}

void Engine::Update()
{
	Events();
	Uint32 current_time = SDL_GetTicks();
	if (current_time == last_time) // damn we're fast
	{
		elapsedtime = 0.0;
		return;
	}
	
	elapsedtime = (double)(current_time-last_time)/1000.0;
	last_time = current_time;

	if (elapsedtime > 0.1)
		elapsedtime = 0.1;

	// TEMP HACK - use a fixed timestep for now, even though this isn't linked to framerate.
	// Too lazy to do this properly, but want consistent simulation results because I was testing the stability of the integrator.
	elapsedtime = 0.016;

	// double speed = 0.1;
	// if (dx != 0 || dy != 0 || dz != 0)
	// {
	// 	fieldRenderer->NudgeInteractive((float)dx * elapsedtime * speed, (float)dy * elapsedtime * speed, (float)dz * elapsedtime * speed);
	// }

	pSim->Update(elapsedtime * dtMult);
	field->UpdateField(elapsedtime * (float)dtMult);
}

void Engine::Draw()
{
// 	if (elapsedtime == 0.0f)
// 		return;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (float)width/(float)height, 0.0001, 1024.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glRotatef(-yr, 1.0, 0.0, 0.0);
	glRotatef(-xr, 0.0, 1.0, 0.0);
	glTranslatef(-x, -y, -z);

	glLightfv(GL_LIGHT0, GL_AMBIENT, l0_ambient);
	glLightfv(GL_LIGHT0, GL_POSITION, l0_position);

	glRotatef(objectPitch, 1, 0, 0);
	glRotatef(objectYaw, 0, 1, 0);

	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);  
	glColor4f(0.5f,0.5f,0.5f,1);
	glVertex3f(-1.0f, -1.0f, -1.0f);  glVertex3f( 1.0f, -1.0f, -1.0f); 
	glVertex3f( 1.0f, -1.0f, -1.0f);  glVertex3f( 1.0f,  1.0f, -1.0f); 
	glVertex3f( 1.0f,  1.0f, -1.0f);  glVertex3f(-1.0f,  1.0f, -1.0f); 
	glVertex3f(-1.0f,  1.0f, -1.0f);  glVertex3f(-1.0f, -1.0f, -1.0f); 

	glVertex3f(-1.0f, -1.0f, 1.0f);  glVertex3f( 1.0f, -1.0f, 1.0f); 
	glVertex3f( 1.0f, -1.0f, 1.0f);  glVertex3f( 1.0f,  1.0f, 1.0f); 
	glVertex3f( 1.0f,  1.0f, 1.0f);  glVertex3f(-1.0f,  1.0f, 1.0f); 
	glVertex3f(-1.0f,  1.0f, 1.0f);  glVertex3f(-1.0f, -1.0f, 1.0f); 

	glVertex3f(-1.0f, -1.0f, 1.0f);  glVertex3f(-1.0f, -1.0f, -1.0f);
	glVertex3f( 1.0f, -1.0f, 1.0f);  glVertex3f( 1.0f, -1.0f, -1.0f);
	glVertex3f( 1.0f,  1.0f, 1.0f);  glVertex3f( 1.0f,  1.0f, -1.0f);
	glVertex3f(-1.0f,  1.0f, 1.0f);  glVertex3f(-1.0f,  1.0f, -1.0f);
	glEnd(); 

	field->DrawField();
//	fieldRenderer->DrawField();
//	pSys->Draw(z);
	pSim->Draw(zoom);

	/*    Setup 2D View - (0,0,xres,yres)    */
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0f, width, height, 0.0f, -10, 10);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glPushAttrib(GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);
	tweakSys->Draw();
//	font->print(10, 10, "mom and dad! it works, wow... who would have guessed???");
	glPopAttrib();

	SDL_GL_SwapWindow(sdl_window);
}

Engine::~Engine()
{
//	delete fieldRenderer;
//	delete pSys;
	delete pSim;
	delete field;

 	SDL_GL_DeleteContext(sdl_gl_context);
 	SDL_DestroyWindow(sdl_window);

	SDL_Quit();
}

void Engine::HandleKey(SDL_Keycode key, bool down)
{
	switch(key)
	{
	case SDLK_w:
		if (down)
			dy += 1;
		else
			dy -= 1;
		break;
	case SDLK_a:
		if (down)
			dx -= 1;
		else
			dx += 1;
		break;
	case SDLK_s:
		if (down)
			dy -= 1;
		else
			dy += 1;
		break;
	case SDLK_d:
		if (down)
			dx += 1;
		else
			dx -= 1;
		break;
	case SDLK_q:
		if (down)
			dz -= 1;
		else
			dz += 1;
		break;
	case SDLK_e:
		if (down)
			dz += 1;
		else
			dz -= 1;
		break;
	case SDLK_z:
		gTweaks.saveTweaks();
		break;
	case SDLK_t:
		gTweaks.loadTweaks();
		break;
	case SDLK_r:
		if (down)
		{
//			pSys->Reset();
			pSim->Reset();
		}
		break;
	// case SDLK_RIGHTBRACKET:
	// 	if (down)
	// 		fieldRenderer->IncNumSteps();
	// 	break;
	// case SDLK_LEFTBRACKET:
	// 	if (down)
	// 		fieldRenderer->DecNumSteps();
	// 	break;
	case SDLK_KP_PLUS:
		if (down)
			dtMult *= 2;
		break;
	case SDLK_KP_MINUS:
		if (down)
			dtMult /= 2;
		break;
	// case SDLK_SPACE:
	// 	if (down)
	// 		fieldRenderer->Save();
	// 	break;
	case SDLK_ESCAPE:
		throw(0); 
		break;
	}
}

void Engine::Events()
{
	static bool mbuttondown1 = false;
	static bool mbuttondown2 = false;
	static bool mbuttondown3 = false;
	SDL_Event event;
	while( SDL_PollEvent(&event) )
	{
		switch( event.type )
		{
		case SDL_KEYDOWN:
			HandleKey(event.key.keysym.sym, true);
			break;
		case SDL_KEYUP:
			HandleKey(event.key.keysym.sym, false);
			break;
   		case SDL_MOUSEBUTTONDOWN:
			if (event.button.button == SDL_BUTTON_LEFT)
				if (tweakSys->mouseDown(event.button.x, event.button.y))
				{
					tweakingUI = true;
					break;
				}

			if (!mbuttondown1 && event.button.button == SDL_BUTTON_LEFT)
				mbuttondown1 = true;
			else if (!mbuttondown2 && event.button.button == SDL_BUTTON_RIGHT)
				mbuttondown2 = true;
			else if (!mbuttondown3 && event.button.button == SDL_BUTTON_MIDDLE)
				mbuttondown3 = true;
			break;
		case SDL_MOUSEBUTTONUP:
			if (tweakingUI)
			{
				tweakSys->mouseUp();
				tweakingUI = false;
			}
			if (mbuttondown1 && event.button.button == SDL_BUTTON_LEFT)
				mbuttondown1 = false;
			if (mbuttondown2 && event.button.button == SDL_BUTTON_RIGHT)
				mbuttondown2 = false;
			if (mbuttondown3 && event.button.button == SDL_BUTTON_MIDDLE)
				mbuttondown3 = false;
            break;
		case SDL_MOUSEWHEEL:
		    if (event.wheel.y < 0)
				zoom *= 0.8f;
		    else
				zoom *= 1.25f;
			break;
		case SDL_MOUSEMOTION:
			if (tweakingUI)
				tweakSys->mouseMove(event.motion.x, event.motion.y);
			else if (mbuttondown1)
			{
				objectYaw += event.motion.xrel * elapsedtime * 10.0f;
				objectPitch += event.motion.yrel * elapsedtime * 10.0f;
			}
			else if (mbuttondown2)
			{
				float delta = (float)(event.motion.xrel + event.motion.yrel)/10;
				z *= 1 - delta * elapsedtime;
				if (z < 0.0000001f)
					z = 0.0000001f;
			}
			else if (mbuttondown3)
			{
				y += (event.motion.yrel)* z * elapsedtime * 0.1f;
				x -= (event.motion.xrel)* z * elapsedtime * 0.1f;
			}
			break;
		case SDL_QUIT:
			throw(0); 
			break;
		}
    }
}
