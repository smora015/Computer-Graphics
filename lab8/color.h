// CS130 Fall 2013: Computer Graphics
// color.h
//
// You may modify this file and are encouraged to do so
/////////////////////////////////////////
#ifndef __COLOR3D_H__
#define __COLOR3D_H__

#include <cmath>

// This is an overly simplified container
// for an RGB Color
struct Color3d
{
	double r;
	double g;
	double b;
	
	Color3d()
		: r(0.0), g(0.0), b(0.0)
	{}
	Color3d(double r, double g, double b)
		: r(r), g(g), b(b)
	{}

	// Changes color based on the rotation by
   // the color wheel, ...
	// using steps RGB -> HSL -> transform -> RGB 
	void rotateHue(double degrees)
	{
		// hue, saturation, and luminosity
		double h = 0.0;
		double s = 0.0;
		double l = 0.0;
		double cmin = std::min(r, std::min(g, b));
		double cmax = std::max(r, std::max(g, b));

		// Convert to HSV
		l = (cmax - cmin) / 2.0;
		if (cmax - cmin != 0.0)
		{
			if (cmax == r)
			{
				h = (g - b) / (cmax - cmin);
				h = h < 0.0 ? h + 6.0 : h;
			}
			else if (cmax == g)
				h = (b - r) / (cmax - cmin) + 2.0;
			else
				h = (r - g) / (cmax - cmin) + 4.0;
			h *= 60;
			s = (cmax - cmin) / (1 - abs(2 * l - 1));
		}

		// Rotate Hue
		h += degrees;
		h = h > 360.0 ? h - 360.0 : h;
		h = h < 0.0 ? h + 360.0 : h;

		// Convert back to RGB and reassign
		double c = 1 - std::abs(2.0 * l - 1.0) * s;
		double x = c * (1 - std::abs(std::fmod(h / 60.0, 2.0) - 1.0));
		double m = l - c / 2.0;;
		
		r = 0.0;
		g = 0.0;
		b = 0.0;
		if (h < 60.0) {
			r = c;
			g = x;
		} else if (h < 120.0) {
			r = x;
			g = c;
		} else if (h < 180.0) {
			g = c;
			b = x;
		} else if (h < 240.0) {
			g = x;
			b = c;
		} else if (h < 300.0) {
			b = c;
			r = x;
		} else {
			b = x;
			r = c;
		}
		r += m;
		g += m;
		b += m;
	}
};

#endif
