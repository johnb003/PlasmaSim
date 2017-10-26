#ifndef TWEAKS_H
#define TWEAKS_H

#include "adScalar.h"
#include <vector>
#include <string>

class Font;

class Tweaks
{
public:
	Tweaks(){}
	void Init(int xRes, int yRes, Font *_font)
	{
		windowResX = xRes;
		windowResY = yRes;
		font = _font;
		rightAnchor = 150;
		activeID = -1;
	}

	enum ScaleStyle
	{
		ScaleStyle_Logarithmic,
		ScaleStyle_Exponential,
		ScaleStyle_Linear,
	};

	struct Tweak
	{
		std::string name;
		ad::Scalar *var;
		ad::Scalar min_val;
		ad::Scalar max_val;
		ad::Scalar defaultValue;
		ScaleStyle scaleType;

		Tweak(const char *name, ad::Scalar *var, ad::Scalar min_val, ad::Scalar max_val, ScaleStyle scaleType = ScaleStyle_Linear)
			: name(name), var(var), min_val(min_val), max_val(max_val), scaleType(scaleType)
		{
			defaultValue = *var;
		}
	};

	void AddTweak(const char *name, ad::Scalar *variable, ad::Scalar min_val, ad::Scalar max_val, ScaleStyle = ScaleStyle_Linear);

	// expects an orthographic view, with top left at 0,0 and bottom right at 
	void Draw();

	// should grab input? returns true if the click is aquired by a control.
	// subsequent mouse move events should be sent to mouseMove(x,y), until mouseUp() is called.
	bool mouseDown(int x, int y);
	void mouseMove(int x, int y);
	void mouseUp();

	// writes current tweaks to a file
	void saveTweaks();
	// loads the last saved file
	void loadTweaks();

private:
	std::vector<Tweak> tweaks;

	int windowResX, windowResY;
	int rightAnchor;
	int activeID;
	Font *font;
};

extern Tweaks gTweaks;

#endif // TWEAKS_H