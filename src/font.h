#ifndef TEXT_H
#define TEXT_H

class Font
{
	constexpr static const float X_SPACE_SIZE = 0.6f;
	constexpr static const float Y_SPACE_SIZE = 1.25f;

	unsigned int f_texture;
	int size;

public:
	Font(const char *filename, int size = 16);
	void print(int x, int y, const char *text);

	int CharWidth() { return size * X_SPACE_SIZE; }
	int CharHeight() { return size * Y_SPACE_SIZE; }
};

#endif