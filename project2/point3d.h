// CS130 Fall 2013: Computer Graphics
// point3d.h
//
// This file does not need to be modified
/////////////////////////////////////////
#ifndef __POINT3D_H__
#define __POINT3D_H__

// A simple wrapper for to store 3D vectors
typedef struct Point3D
{
	float x;
	float y;
	float z;

	Point3D() : x(0.0), y(0.0), z(0.0)
	{}

	Point3D(float x, float y, float z)
		: x(x), y(y), z(z)
	{}

	Point3D(const Point3D & v)
		: x(v.x), y(v.y), z(v.z)
	{}

	Point3D operator+(const Point3D & rhs) const
	{ return Point3D(x + rhs.x, y + rhs.y, z + rhs.z); }
	Point3D operator-(const Point3D & rhs) const
	{ return Point3D(x - rhs.x, y + rhs.y, z - rhs.z); }
	Point3D operator*(float rhs) const
	{ return Point3D(x * rhs, y * rhs, z * rhs); }
	Point3D operator/(float rhs) const
	{ return Point3D(x / rhs, y / rhs, z / rhs); }
	Point3D operator+=(const Point3D & rhs)
	{ x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
	Point3D operator-=(const Point3D & rhs)
	{ x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
	Point3D operator*=(float rhs)
	{ x *= rhs; y *= rhs; z *= rhs; return *this; }
	Point3D operator/=(float rhs)
	{ x /= rhs; y /= rhs; z /= rhs; return *this; }

	float magnitude() const
	{ return sqrt(x * x + y * y + z * z); }
	void normalize()
	{ *this /= magnitude(); }
	float dot(const Point3D & rhs) const
	{
		return x * rhs.x + y * rhs.y + z * rhs.z;
	}
	Point3D cross(const Point3D & rhs) const
	{
		return Point3D(y * rhs.z - z * rhs.y,
					z * rhs.x - x * rhs.z,
					x * rhs.y - y * rhs.x);
	}
} Point3D;

#endif
