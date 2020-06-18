/*
Uint.hpp implements 1, 2, 3, and 4 dimensional unsigned integer vectors
The vector classes are called uint, Uint2, Uint3, and Uint4.
*/

#ifndef _UINT_HPP_INCLUDED_
#define _UINT_HPP_INCLUDED_

#include <cmath>

namespace gfd
{

typedef unsigned int uint;

class Uint2 {
public:
	Uint2() {}
	Uint2(const uint i) { x = y = i; }
	Uint2(const uint xx, const uint yy) { x = xx; y = yy; }
	Uint2 &operator+=(const Uint2 &i) { x += i.x; y += i.y; return (*this); }
	Uint2 operator+(const Uint2 &i) const { return Uint2(x + i.x, y + i.y); }
	Uint2 &operator-=(const Uint2 &i) { x -= i.x; y -= i.y; return (*this); }
	Uint2 operator-(const Uint2 &i) const { return Uint2(x - i.x, y - i.y); }
	uint dot(const Uint2 &i) const { return x * i.x + y * i.y; }
	Uint2 &operator*=(const uint j) { x *= j; y *= j; return (*this); }
	Uint2 operator*(const uint j) const { return Uint2(j * x, j * y); }
	friend Uint2 operator*(const uint j, const Uint2 &i) { return i * j; }
	bool operator==(const Uint2 &v) const { return x == v.x && y == v.y; }
	bool operator!=(const Uint2 &v) const { return x != v.x || y != v.y; }

	uint x;
	uint y;
};

class Uint3 {
public:
	Uint3() {}
	Uint3(const uint i) { x = y = z = i; }
	Uint3(const uint xx, const uint yy, const uint zz) { x = xx; y = yy; z = zz; }
	Uint3(const Uint2 &u, const uint zz) { x = u.x; y = u.y; z = zz; }
	Uint3 &operator+=(const Uint3 &i) { x += i.x; y += i.y; z += i.z; return (*this); }
	Uint3 operator+(const Uint3 &i) const { return Uint3(x + i.x, y + i.y, z + i.z); }
	Uint3 &operator-=(const Uint3 &i) { x -= i.x; y -= i.y; z -= i.z; return (*this); }
	Uint3 operator-(const Uint3 &i) const { return Uint3(x - i.x, y - i.y, z - i.z); }
	uint dot(const Uint3 &i) const { return x * i.x + y * i.y + z * i.z; }
	Uint3 &operator*=(const uint j) { x *= j; y *= j; z *= j; return (*this); }
	Uint3 operator*(const uint j) const { return Uint3(j * x, j * y, j * z); }
	friend Uint3 operator*(const uint j, const Uint3 &i) { return i * j; }
	bool operator==(const Uint3 &v) const { return x == v.x && y == v.y && z == v.z; }
	bool operator!=(const Uint3 &v) const { return x != v.x || y != v.y || z != v.z; }

	Uint2 toUint2() const {	return Uint2(x, y); }

	uint x;
	uint y;
	uint z;
};

class Uint4 {
public:
	Uint4() {}
	Uint4(const uint i) { x = y = z = i; }
	Uint4(const uint xx, const uint yy, const uint zz, const uint tt) { x = xx; y = yy; z = zz; t = tt; }
	Uint4(const Uint2 &u, const uint zz, const uint tt) { x = u.x; y = u.y; z = zz; t = tt; }
	Uint4(const Uint3 &u, const uint tt) { x = u.x; y = u.y; z = u.z; t = tt; }
	Uint4 &operator+=(const Uint4 &i) { x += i.x; y += i.y; z += i.z; t += i.t; return (*this); }
	Uint4 operator+(const Uint4 &i) const { return Uint4(x + i.x, y + i.y, z + i.z, t + i.t); }
	Uint4 &operator-=(const Uint4 &i) { x -= i.x; y -= i.y; z -= i.z; t -= i.t; return (*this); }
	Uint4 operator-(const Uint4 &i) const { return Uint4(x - i.x, y - i.y, z - i.z, t - i.t); }
	uint dot(const Uint4 &i) const { return x * i.x + y * i.y + z * i.z + t * i.t; }
	Uint4 &operator*=(const uint j) { x *= j; y *= j; z *= j; t *= j; return (*this); }
	Uint4 operator*(const uint j) const { return Uint4(j * x, j * y, j * z, j * t); }
	friend Uint4 operator*(const uint j, const Uint4 &i) { return i * j; }
	bool operator==(const Uint4 &v) const { return x == v.x && y == v.y && z == v.z && t == v.t; }
	bool operator!=(const Uint4 &v) const { return x != v.x || y != v.y || z != v.z || t != v.t; }

	Uint2 toUint2() const {	return Uint2(x, y); }
	Uint3 toUint3() const {	return Uint3(x, y, z); }

	uint x;
	uint y;
	uint z;
	uint t;
};



}

#endif //_UINT_HPP_INCLUDED_
