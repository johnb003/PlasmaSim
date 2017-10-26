#include "tweaks.h"
#include "font.h"
#include "adScalar.h"

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


#include <math.h>
#include <sstream>
#include <iomanip>
#include <fstream>

using std::vector;
using ad::min;
using ad::max;

Tweaks gTweaks;

static const int SPACES_AFTER_LABEL = 15;
static const int TWEAKS_Y_POS = 10;

static const int BOX_WIDTH = 300;

void Tweaks::AddTweak(const char *name, ad::Scalar *var, ad::Scalar min_val, ad::Scalar max_val, ScaleStyle scaleType)
{
	tweaks.push_back(Tweak(name, var, min_val, max_val, scaleType));
	rightAnchor = max(rightAnchor, ((int)tweaks.back().name.length() + SPACES_AFTER_LABEL) * font->CharWidth());
}

// ad::Scalar logn(ad::Scalar x, ad::Scalar b)
// {
// 	return log10(x) / log10(b);
// }


void Tweaks::Draw()
{
	int y = TWEAKS_Y_POS;
	int i = 0;
	// vector<Tweak>::iterator it;
	for (auto it = tweaks.begin(); it != tweaks.end(); ++it)
	{
		std::stringstream ss;
		ss << std::setprecision(4) << (double)*(it->var);

		font->print(rightAnchor - font->CharWidth() * (ss.str().length() + 2), y, ss.str().c_str() );
		font->print(0, y, it->name.c_str() );

		ad::Scalar value = *it->var;
		ad::Scalar low = it->min_val;
		ad::Scalar high = it->max_val;

		if (it->scaleType == ScaleStyle_Exponential)
		{
			value = log(value);
			low = log(low);
			high = log(high);
		}
		else if (it->scaleType == ScaleStyle_Logarithmic)
		{
			value = exp(value);
			low = exp(value);
			high = exp(high);
		}

		// Cap the value before interpolating it.
		value = min(max(low, value), high);
		ad::Scalar percentageFull = (value - low) / (high - low);

		if (i == activeID)
			glColor4f(1, 0, 0, 0.7f);
		else
			glColor4f(1, 1, 0, 0.7f);

		glBegin(GL_QUADS);
		glVertex3i(rightAnchor, y, 0);
		glVertex3i(rightAnchor, y + font->CharHeight(), 0);
		glVertex3i(rightAnchor + int(BOX_WIDTH * percentageFull), y + font->CharHeight(), 0);
		glVertex3i(rightAnchor + int(BOX_WIDTH * percentageFull), y, 0);
		glEnd();

		glColor4f(1, 1, 1, 1);
		glBegin(GL_LINE_LOOP);
		glVertex3i(rightAnchor, y, 0);
		glVertex3i(rightAnchor, y + font->CharHeight(), 0);
		glVertex3i(rightAnchor + BOX_WIDTH, y + font->CharHeight(), 0);
		glVertex3i(rightAnchor + BOX_WIDTH, y, 0);
		glEnd();

		y += font->CharHeight();
		i++;
	}
}


bool Tweaks::mouseDown(int x, int y)
{
	int yStart = 10;
	for (unsigned int i = 0; i < tweaks.size(); i++)
	{
		if (y > yStart && y < yStart + font->CharHeight() &&
			x > rightAnchor && x < rightAnchor+BOX_WIDTH)
		{
			activeID = i;

			if (tweaks[i].scaleType == ScaleStyle_Exponential)
			{
				ad::Scalar low = log(tweaks[i].min_val);
				ad::Scalar high = log(tweaks[i].max_val);
				ad::Scalar val = (((float)x - rightAnchor) / BOX_WIDTH) * (high - low) + low;
				if (val < low)
					val = low;
				if (val > high)
					val = high;
				*tweaks[i].var = exp( val );
			}

			return true;
		}

		yStart += font->CharHeight();
	}

	return false;
}

void Tweaks::mouseMove(int x, int y)
{
	if (tweaks[activeID].scaleType == ScaleStyle_Exponential)
	{
		ad::Scalar low = log(tweaks[activeID].min_val);
		ad::Scalar high = log(tweaks[activeID].max_val);
		ad::Scalar val = (((float)x - rightAnchor) / BOX_WIDTH) * (high-low) + low;
		if (val < low)
			val = low;
		if (val > high)
			val = high;
		*tweaks[activeID].var = exp(val);
	}
}

void Tweaks::mouseUp()
{
	activeID = -1;
}

void Tweaks::saveTweaks()
{
	std::ofstream outFile("savedTweaks.txt");
	for (std::vector<Tweak>::iterator tweak = tweaks.begin(); tweak != tweaks.end(); ++tweak)
	{
		outFile << *tweak->var << std::endl;
	}
	outFile.close();
}

void Tweaks::loadTweaks()
{
	std::ifstream inFile("savedTweaks.txt", std::ios::in);
	if (inFile.good())
	{
		for (std::vector<Tweak>::iterator tweak = tweaks.begin(); tweak != tweaks.end(); ++tweak)
		{
			inFile >> *tweak->var;
		}
	}
	inFile.close();
}
