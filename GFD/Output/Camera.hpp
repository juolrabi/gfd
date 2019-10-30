/*
Camera is a tool for drawing meshes.
*/

#ifndef _CAMERA_HPP_INCLUDED_
#define _CAMERA_HPP_INCLUDED_

#include "../Types/Buffer.hpp"
#include "../Types/Text.hpp"
#include "Picture.hpp"

namespace gfd
{

class Camera
{
public:
	Camera() { clear(); }
	virtual ~Camera() {clear(); }
	void clear();

	void initPosition(const Vector4 &pcam, const Vector4 &ptarget, const Vector4 &x, const Vector4 &y);
	const Vector4 &getPosition() const { return m_p; }

	void initSvg(const uint width, const uint height);
	bool saveSvg(const std::string &path);

	void initPicture(Picture *const pic);
	const Vector4 getPicturePosition(const uint x, const uint y) const;

	void drawPoint(const Vector4 &p, const Vector3 &c);
	void drawLine(const Vector4 &p0, const Vector4 &p1, const Vector3 &c);
	void drawPolygon(const Buffer<Vector4> &p, const Buffer<Vector3> &c, const Vector4 &proj = Vector4(0,0,0,1));

protected:

	bool m_isvg;
	uint m_svgwid;
	uint m_svghei;
	Text m_svg;

	Picture *m_pic;

	Vector4 m_p;
	Vector4 m_x;
	Vector4 m_y;
	Vector4 m_z;

	bool insertSvgPolygon(const Buffer<Vector2> &p, const Vector3 &c = Vector3(1,1,1));
	bool insertSvgLine(const Vector2 &p0, const Vector2 &p1, const Vector3 &c = Vector3(1,1,1));
	void drawPictureTriangle(const Vector3 &t0, const Vector3 &t1, const Vector3 &t2, const Vector3 &c0, const Vector3 &c1, const Vector3 &c2);

	Vector2 getSvgCoordinates(const Vector4 &p) const;
	Vector3 getPictureCoordinates(const Vector4 &p) const;
	uint limited(const double x, const uint wid) const;

};

}

#endif //_CAMERA_HPP_INCLUDED_
