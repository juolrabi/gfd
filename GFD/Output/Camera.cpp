#include "Camera.hpp"
#include <fstream>
#include <iostream>

using namespace gfd;

void Camera::clear()
{
	m_pic = NULL;

	m_isvg = false;
	m_svg.clear();

	initPosition(Vector4(0,0,100,0), Vector4(0,0,0,0), Vector4(2,0,0,0), Vector4(0,-2,0,0));
}


void Camera::initPosition(const Vector4 &pcam, const Vector4 &ptarget, const Vector4 &x, const Vector4 &y)
{
	//pcam is the camera position
	//the line z from pcam to ptarget is at the center of the picture
	//x and y is the gradient vectors for picture coordinates
	//x and y are forced to be orthogonal with z
	m_p = pcam;
	const Vector4 d = ptarget - m_p;
	m_z = d / d.lensq();
	m_x = x - x.dot(d) * m_z;
	m_y = y - y.dot(d) * m_z;
}

void Camera::initSvg(const uint width, const uint height)
{
	if(m_isvg) m_svg.clear();
	m_isvg = true;
	m_svgwid = width;
	m_svghei = height;

	m_svg << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << std::endl;
	m_svg << "<!-- Created with Inkscape (http://www.inkscape.org/) -->" << std::endl;
	m_svg << "<svg" << std::endl;
	m_svg << "   xmlns:dc=\"http://purl.org/dc/elements/1.1/\"" << std::endl;
	m_svg << "   xmlns:cc=\"http://web.resource.org/cc/\"" << std::endl;
	m_svg << "   xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"" << std::endl;
	m_svg << "   xmlns:svg=\"http://www.w3.org/2000/svg\"" << std::endl;
	m_svg << "   xmlns=\"http://www.w3.org/2000/svg\"" << std::endl;
	m_svg << "   xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\"" << std::endl;
	m_svg << "   xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\"" << std::endl;
	m_svg << "   width=\"744.09448819\"" << std::endl;
	m_svg << "   height=\"1052.3622047\"" << std::endl;
	m_svg << "   id=\"svg2\"" << std::endl;
	m_svg << "   sodipodi:version=\"0.32\"" << std::endl;
	m_svg << "   inkscape:version=\"0.44\">" << std::endl;
	m_svg << "  <defs" << std::endl;
	m_svg << "     id=\"defs4\" />" << std::endl;
	m_svg << "  <sodipodi:namedview" << std::endl;
	m_svg << "     id=\"base\"" << std::endl;
	m_svg << "     pagecolor=\"#ffffff\"" << std::endl;
	m_svg << "     bordercolor=\"#666666\"" << std::endl;
	m_svg << "     borderopacity=\"1.0\"" << std::endl;
	m_svg << "     gridtolerance=\"10000\"" << std::endl;
	m_svg << "     guidetolerance=\"10\"" << std::endl;
	m_svg << "     objecttolerance=\"10\"" << std::endl;
	m_svg << "     inkscape:pageopacity=\"0.0\"" << std::endl;
	m_svg << "     inkscape:pageshadow=\"2\"" << std::endl;
	m_svg << "     inkscape:zoom=\"0.35\"" << std::endl;
	m_svg << "     inkscape:cx=\"375\"" << std::endl;
	m_svg << "     inkscape:cy=\"520\"" << std::endl;
	m_svg << "     inkscape:document-units=\"px\"" << std::endl;
	m_svg << "     inkscape:current-layer=\"layer1\"" << std::endl;
	m_svg << "     inkscape:window-width=\"853\"" << std::endl;
	m_svg << "     inkscape:window-height=\"573\"" << std::endl;
	m_svg << "     inkscape:window-x=\"44\"" << std::endl;
	m_svg << "     inkscape:window-y=\"44\" />" << std::endl;
	m_svg << "  <metadata" << std::endl;
	m_svg << "     id=\"metadata7\">" << std::endl;
	m_svg << "    <rdf:RDF>" << std::endl;
	m_svg << "      <cc:Work" << std::endl;
	m_svg << "         rdf:about=\"\">" << std::endl;
	m_svg << "        <dc:format>image/svg+xml</dc:format>" << std::endl;
	m_svg << "        <dc:type" << std::endl;
	m_svg << "           rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />" << std::endl;
	m_svg << "      </cc:Work>" << std::endl;
	m_svg << "    </rdf:RDF>" << std::endl;
	m_svg << "  </metadata>" << std::endl;
	m_svg << "  <g" << std::endl;
	m_svg << "     inkscape:label=\"Layer 1\"" << std::endl;
	m_svg << "     inkscape:groupmode=\"layer\"" << std::endl;
	m_svg << "     id=\"layer1\">" << std::endl;
}

bool Camera::insertSvgPolygon(const Buffer<Vector2> &p, const Vector3 &c)
{
	if(!m_isvg) return false;

	const uint col = 0x101010 + 0x010000 * uint(239.9 * c.x) + 0x000100 * uint(239.9 * c.y) + uint(239.9 * c.z);
	m_svg << "    <path" << std::endl;
	m_svg << "       style=\"fill:#" << std::hex << col << std::dec << ";fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"" << std::endl;
	m_svg << "       d=\"M " << p[0].x << "," << p[0].y;
	for(uint i=1; i<p.size(); i++)
	{
		m_svg << " L " << p[i].x << "," << p[i].y;
	}
	m_svg << " z\" />" << std::endl;

	return true;
}


bool Camera::insertSvgLine(const Vector2 &p0, const Vector2 &p1, const Vector3 &c)
{
	if(!m_isvg) return false;

	const uint col = 0x101010 + 0x010000 * uint(239.9 * c.x) + 0x000100 * uint(239.9 * c.y) + uint(239.9 * c.z);
	m_svg << "    <path" << std::endl;
	m_svg << "       style=\"fill:none;stroke:#" << std::hex << col << std::dec << ";stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"" << std::endl;
	m_svg << "       d=\"M " << p0.x << "," << p0.y << " L " << p1.x << "," << p1.y << " \" />" << std::endl;
	return true;
}

bool Camera::saveSvg(const std::string &path)
{
	if(!m_isvg) return false;

	m_svg << "  </g>" << std::endl;
	m_svg << "</svg>" << std::endl;

	return m_svg.save(path);
}


void Camera::initPicture(Picture *const pic)
{
	m_pic = pic;
}

const Vector4 Camera::getPicturePosition(const uint x, const uint y) const
{
	if(m_pic == NULL) return Vector4(0,0,0,0);
	const double px = (x + 0.5) / double(m_pic->getWidth()) - 0.5;
	const double py = (y + 0.5) / double(m_pic->getHeight()) - 0.5;
	return m_p + px * m_x / m_x.lensq() + py * m_y / m_y.lensq() + m_z / m_z.lensq();
}

void Camera::drawPoint(const Vector4 &p, const Vector3 &c)
{
	if(m_isvg) // draw svg
	{
		const Vector2 mt = getSvgCoordinates(p);
		Buffer<Vector2> t(4);
		t[0] = mt + Vector2(-1,-1);
		t[1] = mt + Vector2(1,-1);
		t[2] = mt + Vector2(1,1);
		t[3] = mt + Vector2(-1,1);
		insertSvgPolygon(t, c);
	}

	// draw picture
	if(m_pic == NULL || m_pic->empty()) return;
	const uint wid = m_pic->getWidth();
	const uint hei = m_pic->getHeight();

	Vector3 t = getPictureCoordinates(p);
	t.x += 0.5;
	t.y += 0.5;
	if(t.x < 0.0 || wid <= t.x || t.y < 0.0 || hei <= t.y) return;
	if(t.z > m_pic->getColor(uint(t.x),uint(t.y)).t) m_pic->setColor(uint(t.x), uint(t.y), Vector4(c,t.z));
}

void Camera::drawLine(const Vector4 &p0, const Vector4 &p1, const Vector3 &c)
{
	if(m_isvg) // draw svg
	{
		const Vector2 t0 = getSvgCoordinates(p0);
		const Vector2 t1 = getSvgCoordinates(p1);
		insertSvgLine(t0, t1, c);
	}

	// draw picture
	if(m_pic == NULL || m_pic->empty()) return;
	const uint wid = m_pic->getWidth();
	const uint hei = m_pic->getHeight();

	const Vector3 t0 = getPictureCoordinates(p0);
	const Vector3 t1 = getPictureCoordinates(p1);
	Vector3 d = t1 - t0;
	if(fabs(d.x) > fabs(d.y))
	{
		uint x, ex;
		if(d.x > 0.0)
		{
			x = limited(t0.x + 0.5, wid);
			ex = limited(t1.x + 0.5, wid);
		}
		else
		{
			x = limited(t1.x + 0.5, wid);
			ex = limited(t0.x + 0.5, wid);
		}
		while(x < ex)
		{
			const uint y = uint(t0.y + (x - t0.x) / d.x * d.y);
			if(y < hei)
			{
				const double dist = t0.z + (x - t0.x) / d.x * d.z;
				if(dist > m_pic->getColor(x,y).t)
				{
					m_pic->setColor(x, y, Vector4(c,dist));
				}
			}
			x++;
		}
	}
	else
	{
		uint y, ey;
		if(d.y > 0.0)
		{
			y = limited(t0.y + 0.5, hei);
			ey = limited(t1.y + 0.5, hei);
		}
		else
		{
			y = limited(t1.y + 0.5, hei);
			ey = limited(t0.y + 0.5, hei);
		}
		while(y < ey)
		{
			const uint x = uint(t0.x + (y - t0.y) / d.y * d.x);
			if(x < wid)
			{
				const double dist = t0.z + (y - t0.y) / d.y * d.z;
				if(dist > m_pic->getColor(x,y).t)
				{
					m_pic->setColor(x, y, Vector4(c,dist));
				}
			}
			y++;
		}
	}
}

void Camera::drawPolygon(const Buffer<Vector4> &p, const Buffer<Vector3> &c, const Vector4 &proj)
{
	uint i;

	// check if function parameters are ok
	if(p.size() < 3) return;

	// compute light
	Vector4 mn(0,0,0,0);
	for(i=2; i<p.size(); i++) mn += ThreeVector4(p[i-1] - p[0], p[i] - p[0], proj).dual();
	const Vector3 col = c.empty() ? Vector3(0,0,0) : c[0];
	Vector3 ms(0.0, 0.0, 0.0);
	Buffer<Vector3> s(p.size());
	const double mnsq = mn.lensq();
	bool turn = false;
	for(i=0; i<p.size(); i++)
	{
		const Vector4 d = p[i] - m_p;
		double dot = mn.dot(d) / sqrt(mnsq * d.lensq() + 1e-30);
		//if(dot < 0.0) return; // draw only front face
		if(dot < 0.0) // draw also background
		{
			turn = true;
			dot = -dot;
		}
		s[i] = (i < c.size() ? c[i] : col) * (0.2 + 0.8 * dot);
		ms += s[i];
	}
	ms /= double(p.size());

	if(m_isvg) // draw svg
	{
		Buffer<Vector2> t(p.size());
		for(i=0; i<p.size(); i++) t[i] = getSvgCoordinates(p[i]);
		insertSvgPolygon(t, ms);
	}

	// draw picture
	if(m_pic == NULL || m_pic->empty()) return;

	Vector3 mt(0,0,0);
	Buffer<Vector3> t(p.size());
	for(i=0; i<t.size(); i++)
	{
		t[i] = getPictureCoordinates(p[i]);
		mt += t[i];
	}
	mt /= double(p.size());

	// draw triangles
	for(i=0; i<p.size(); i++)
	{
		const uint j = (i > 0 ? i : p.size()) - 1;
		if(!turn) drawPictureTriangle(mt, t[i], t[j], ms, s[i], s[j]);
		else drawPictureTriangle(mt, t[j], t[i], ms, s[j], s[i]);
	}
}

void Camera::drawPictureTriangle(const Vector3 &t0, const Vector3 &t1, const Vector3 &t2, const Vector3 &c0, const Vector3 &c1, const Vector3 &c2)
{
	const Vector3 face = TwoVector3(t1 - t0, t2 - t0).dual();
	if(face.z > -1e-8) return;
	const uint wid = m_pic->getWidth();
	const uint hei = m_pic->getHeight();

	const Vector3 colx = ((t0.y - t1.y) * (c2 - c0) - (c1 - c0) * (t0.y - t2.y)) / face.z;
	const Vector3 coly = ((c0 - c1) * (t2.x - t0.x) - (t1.x - t0.x) * (c0 - c2)) / face.z;
	const Vector3 col0 = c0 - t0.x * colx - t0.y * coly;

	const double distx = -face.x / face.z;
	const double disty = -face.y / face.z;
	const double dist0 = t0.z - t0.x * distx - t0.y * disty;

	double miny = t0.y, maxy = t0.y;
	if(t1.y < miny) miny = t1.y;
	else if(t1.y > maxy) maxy = t1.y;
	if(t2.y < miny) miny = t2.y;
	else if(t2.y > maxy) maxy = t2.y;

	const uint by = limited(miny + 1.0, hei);
	const uint ey = limited(maxy + 1.0, hei);
	for(uint y=by; y<ey; y++)
	{
		double minx = 0.0;
		double maxx = 0.0;
		if(t0.y < y && y <= t1.y) minx = t0.x + (y - t0.y) / (t1.y - t0.y) * (t1.x - t0.x);
		else if(t1.y < y && y <= t0.y) maxx = t1.x + (y - t1.y) / (t0.y - t1.y) * (t0.x - t1.x);
		if(t1.y < y && y <= t2.y) minx = t1.x + (y - t1.y) / (t2.y - t1.y) * (t2.x - t1.x);
		else if(t2.y < y && y <= t1.y) maxx = t2.x + (y - t2.y) / (t1.y - t2.y) * (t1.x - t2.x);
		if(t2.y < y && y <= t0.y) minx = t2.x + (y - t2.y) / (t0.y - t2.y) * (t0.x - t2.x);
		else if(t0.y < y && y <= t2.y) maxx = t0.x + (y - t0.y) / (t2.y - t0.y) * (t2.x - t0.x);

		const uint bx = limited(minx + 1.0, wid);
		const uint ex = limited(maxx + 1.0, wid);
		for(uint x=bx; x<ex; x++)
		{
			const double dist = dist0 + x * distx + y * disty;
			if(dist > m_pic->getColor(x,y).t)
			{
				Vector3 col = col0 + x * colx + y * coly;
				if(col.x < -1.0) col.x = -1.0;
				else if(col.x > 1.0) col.x = 1.0;
				if(col.y < -1.0) col.y = -1.0;
				else if(col.y > 1.0) col.y = 1.0;
				if(col.z < -1.0) col.z = -1.0;
				else if(col.z > 1.0) col.z = 1.0;
				m_pic->setColor(x, y, Vector4(col, dist));
			}
		}
	}
}

Vector2 Camera::getSvgCoordinates(const Vector4 &p) const
{
	const Vector4 d = p - m_p;
	const double dotz = m_z.dot(p - m_p);
	if(dotz < 1e-8) return Vector2(0,0);
	const double dx = (m_x.dot(d) / dotz + 0.5) * m_svgwid + 0.5;
	const double dy = -m_y.dot(d) / dotz * m_svgwid + 0.5 * m_svghei + 0.5;
	return Vector2(dx, dy);
}

Vector3 Camera::getPictureCoordinates(const Vector4 &p) const
{
	const Vector4 d = p - m_p;
	const double dotz = m_z.dot(p - m_p);
	if(dotz < 1e-8) return Vector3(0,0,0);
	const double dx = (m_x.dot(d) / dotz + 0.5) * m_pic->getWidth() + 0.5;
	const double dy = m_y.dot(d) / dotz * m_pic->getWidth() + 0.5 * m_pic->getHeight() + 0.5;
	return Vector3(dx, dy, 1.0 / dotz);
}

uint Camera::limited(const double x, const uint wid) const
{
  if(x < 0.0) return 0;
  if(x > wid) return wid;
  return uint(x);
}


