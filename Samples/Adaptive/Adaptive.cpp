#include "Adaptive.hpp"
#include <iostream>

using namespace gfd;
using namespace std;

Adaptive::Adaptive(const Vector2 &p0, const Vector2 &p1, const double h) {
	m_dim = 2;

	const Vector2 d = p1 - p0;
	m_h = Vector3(h, h * sqrt(3.0), 1);
	m_n = Uint3(d.x / m_h.x + 0.99999, d.y / m_h.y + 0.99999, 1);
	const Vector2 scale(d.x / (m_n.x * m_h.x), d.y / (m_n.y * m_h.y));
	const double hxy = scale.x * m_h.y;
	const double hyx = scale.y * m_h.x;
	const Vector2 area((uint(d.y / hxy + 0.99999) * hxy - d.y) * d.x, 
						(uint(d.x / hyx + 0.99999) * hyx - d.x) * d.y);
	if(area.x < area.y) m_h *= scale.x;
	else m_h *= scale.y;
	m_n = Uint3(d.x / m_h.x + 0.99999, d.y / m_h.y + 0.99999, 1);
	m_p = Vector3(0.5 * (p0 + p1 - Vector2(m_n.x * m_h.x, m_n.y * m_h.y)), 0);

	m_s.resize(m_n.x * m_n.y);
	for(m_ss=0; m_ss<m_s.size(); m_ss++) m_s[m_ss] = 0;
}

Adaptive::Adaptive(const Vector3 &p0, const Vector3 &p1, const double h) {
	m_dim = 3;

	const Vector3 d = p1 - p0;
	m_h = Vector3(h, h, h);
	m_n = Uint3(d.x / m_h.x + 0.99999, d.y / m_h.y + 0.99999, d.z / m_h.z + 0.99999);
	const Vector3 scale(d.x / (m_n.x * m_h.x), d.y / (m_n.y * m_h.y), d.z / (m_n.z * m_h.z));
	const double hxy = scale.x * m_h.y;
	const double hxz = scale.x * m_h.z;
	const double hyx = scale.y * m_h.x;
	const double hyz = scale.y * m_h.z;
	const double hzx = scale.z * m_h.x;
	const double hzy = scale.z * m_h.y;
	const Vector3 area((uint(d.y / hxy + 0.99999) * hxy - d.y) * d.z * d.x + (uint(d.z / hxz + 0.99999) * hxz - d.z) * d.y * d.x, 
		(uint(d.x / hyx + 0.99999) * hyx - d.x) * d.z * d.y + (uint(d.z / hyz + 0.99999) * hyz - d.z) * d.x * d.y,
		(uint(d.x / hzx + 0.99999) * hzx - d.x) * d.y * d.z + (uint(d.y / hzy + 0.99999) * hzy - d.y) * d.x * d.z);
	if(area.x < area.y && area.x < area.z) m_h *= scale.x;
	else if(area.y < area.z) m_h *= scale.y;
	else m_h *= scale.z;
	m_n = Uint3(d.x / m_h.x + 0.99999, d.y / m_h.y + 0.99999, d.z / m_h.z + 0.99999);
	m_p = 0.5 * (p0 + p1 - Vector3(m_n.x * m_h.x, m_n.y * m_h.y, m_n.z * m_h.z));

	m_s.resize(m_n.x * m_n.y * m_n.z);
	for(m_ss=0; m_ss<m_s.size(); m_ss++) m_s[m_ss] = 0;
}

bool Adaptive::isOutOfBounds(const Vector3 &d) const {
	if(d.x < 0.0 || d.x >= double(m_n.x)) return true; 
	if(m_dim == 1) return false;
	if(d.y < 0.0 || d.y >= double(m_n.y)) return true; 
	if(m_dim == 2) return false;
	if(d.z < 0.0 || d.z >= double(m_n.z)) return true; 
	return false;
}

uint Adaptive::findLevel(const Vector3 &p) const {
	const Vector3 d = (p - m_p) / m_h;
	if(isOutOfBounds(d)) return 0;
	const Uint3 i(d.x, d.y, d.z);
	const uint ii = (i.z * m_n.y + i.y) * m_n.x + i.x; // pixel index
	if(m_s[ii] == 0) return 1; // constant pixel found
	return findLevelRecursive(Vector3(d.x - double(i.x), d.y - double(i.y), d.z - double(i.z)), m_s[ii], 2);
}
uint Adaptive::findLevelRecursive(const Vector3 &p, uint ii, const uint level) const {
	const Vector3 d = 2.0 * p;
	const Uint3 i(d.x, d.y, d.z);
	ii += (i.z * 2 + i.y) * 2 + i.x; // pixel index
	if(m_s[ii] == 0) return level; // constant pixel found
	return findLevelRecursive(Vector3(d.x - double(i.x), d.y - double(i.y), d.z - double(i.z)), m_s[ii], level + 1);
}

uint Adaptive::getBlockSize() const {
	switch(m_dim) {
	case 1: return 2;
	case 2: return 4;
	default: return 8;
	}
}
Vector3 Adaptive::getFlat(const Vector3 &v) const {
	switch(m_dim) {
	case 1: return Vector3(v.x,0,0);
	case 2: return Vector3(v.x,v.y,0);
	default: return v;
	}
}

void Adaptive::splitPixel(const uint ii) {
	if(m_s[ii] != 0) return;
	m_s[ii] = m_ss;
	const uint size = getBlockSize();
	for(uint i=0; i<size; i++) m_s.gather(0, m_ss);
}

void Adaptive::refineSphereBoundary(const Vector3 &p, const double r, const uint level) {
	if(level <= 1) return;
	Uint3 i;
	uint ii = 0;
	for(i.z=0; i.z<m_n.z; i.z++) {
		for(i.y=0; i.y<m_n.y; i.y++) {
			for(i.x=0; i.x<m_n.x; i.x++, ii++) {
				refineSphereBoundaryRecursive((p - m_p) / m_h - Vector3(i.x, i.y, i.z), Vector3(r,r,r) / m_h, ii, level);
			}
		}
	}
}
void Adaptive::refineSphereBoundaryRecursive(const Vector3 &p, const Vector3 &r, const uint ii, const uint level) {
	const Vector3 h = getFlat(Vector3(1,1,1));
	const Vector3 d = (p - 0.5 * h) / r;
	const Vector3 rh = h / r;
	if(fabs(d.len() - 1.0) > 0.7 * rh.len()) return; // do not refine
	if(m_s[ii] == 0) splitPixel(ii);
	
	if(level <= 2) return; 
	const uint s = m_s[ii];
	const Vector3 p2 = 2.0 * p;
	const Vector3 r2 = 2.0 * r;
	refineSphereBoundaryRecursive(p2, r2, s, level-1);
	refineSphereBoundaryRecursive(p2 - Vector3(1,0,0), r2, s+1, level-1);
	if(m_dim == 1) return;
	refineSphereBoundaryRecursive(p2 - Vector3(0,1,0), r2, s+2, level-1);
	refineSphereBoundaryRecursive(p2 - Vector3(1,1,0), r2, s+3, level-1);
	if(m_dim == 2) return;
	refineSphereBoundaryRecursive(p2 - Vector3(0,0,1), r2, s+4, level-1);
	refineSphereBoundaryRecursive(p2 - Vector3(1,0,1), r2, s+5, level-1);
	refineSphereBoundaryRecursive(p2 - Vector3(0,1,1), r2, s+6, level-1);
	refineSphereBoundaryRecursive(p2 - Vector3(1,1,1), r2, s+7, level-1);
}


void Adaptive::createMesh(BuilderMesh &mesh) const {
	// create base grid
	const Uint3 n((m_n.x < 3 ? m_n.x : 3), (m_n.y < 3 ? m_n.y : 3), (m_n.z < 3 ? m_n.z : 3));
	const Vector3 h = getFlat(m_h);
	mesh.createGrid(Vector4(m_p,0), Vector4(m_p+Vector3(h.x * n.x, h.y * n.y, h.z * n.z),0), Vector4(m_h,1));

	Uint3 i;
	uint node = 0;
	for(i.z=0; i.z<n.z; i.z++) {
		for(i.y=0; i.y<n.y; i.y++) {
			for(i.x=0; i.x<n.x; i.x++) {
				node = mesh.insertNode(Vector4(m_p + Vector3(h.x * (i.x + 0.5), h.y * (i.y + 0.5), h.z * (i.z + 0.5)), 0), 0.0, node, true);
			}
		}
	}
	if(m_n.x > n.x) mesh.repeatMiddle(Vector4(m_p.x + 0.99999 * m_h.x, 0,0,0), Vector4(m_h.x,0,0,0), m_n.x - n.x);
	if(m_n.y > n.y) mesh.repeatMiddle(Vector4(0,m_p.y + 0.99999 * m_h.y, 0,0), Vector4(0,m_h.y,0,0), m_n.y - n.y);
	if(m_n.z > n.z) mesh.repeatMiddle(Vector4(0,0,m_p.z + 0.99999 * m_h.z, 0), Vector4(0,0,m_h.z,0), m_n.z - n.z);
	
	// refine recursively
	bool morelevels = true;
	Buffer<bool> more(m_n.x * m_n.y * m_n.z, true);
	for(uint level=0; morelevels; level++) {
		morelevels = false;
		Uint3 i;
		uint ii = 0;
		for(i.z=0; i.z<m_n.z; i.z++) {
			for(i.y=0; i.y<m_n.y; i.y++) {
				for(i.x=0; i.x<m_n.x; i.x++, ii++) {
					if(!more[ii]) continue;
					const Vector3 p = m_p + Vector3(h.x * i.x, h.y * i.y, h.z * i.z);
					more[ii] = createMeshRecursively(mesh, level, node, p, h, m_s[ii]);
					if(more[ii]) morelevels = true;
				}
			}
		}
		cout << "iteration " << level << " " << morelevels << endl;
	}
}

bool Adaptive::createMeshRecursively(BuilderMesh &mesh, const uint level, uint &node, const Vector3 &p, const Vector3 &h, const int s) const {
	if(s == 0) return false; // constant pixel found

	// insert nodes
	const Vector3 r = 0.5 * h;
	if(level == 0) {
		const Vector3 o = p + r;
		const Vector3 rr = 0.5 * r;
		// boundaries
		if(m_dim >= 2) {
			node = mesh.insertNode(Vector4(o + Vector3(r.x,0,0), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(-r.x,0,0), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(0,r.y,0), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(0,-r.y,0), 0), 0.0, node, false);
		}
		if(m_dim >= 3) {
			node = mesh.insertNode(Vector4(o + Vector3(0,0,r.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(0,0,-r.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(-r.x,-r.y,0), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(r.x,-r.y,0), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(-r.x,r.y,0), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(r.x,r.y,0), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(-r.x,0,-r.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(r.x,0,-r.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(-r.x,0,r.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(r.x,0,r.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(0,-r.y,-r.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(0,r.y,-r.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(0,-r.y,r.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(0,r.y,r.z), 0), 0.0, node, false);
		}
		// interior
		node = mesh.insertNode(Vector4(o + Vector3(rr.x,rr.y,rr.z), 0), 0.0, node, false);
		node = mesh.insertNode(Vector4(o + Vector3(-rr.x,rr.y,rr.z), 0), 0.0, node, false);
		if(m_dim >= 2) {
			node = mesh.insertNode(Vector4(o + Vector3(rr.x,-rr.y,rr.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(-rr.x,-rr.y,rr.z), 0), 0.0, node, false);
		}
		if(m_dim >= 3) {
			node = mesh.insertNode(Vector4(o + Vector3(rr.x,rr.y,-rr.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(-rr.x,rr.y,-rr.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(-rr.x,-rr.y,-rr.z), 0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector3(rr.x,-rr.y,-rr.z), 0), 0.0, node, false);
		}
		return true;
	}

	bool morelevels = false;
	morelevels = createMeshRecursively(mesh, level-1, node, p, r, m_s[s]) || morelevels;
	morelevels = createMeshRecursively(mesh, level-1, node, p + Vector3(r.x,0,0), r, m_s[s+1]) || morelevels;
	if(m_dim == 1) return morelevels;
	morelevels = createMeshRecursively(mesh, level-1, node, p + Vector3(0,r.y,0), r, m_s[s+2]) || morelevels;
	morelevels = createMeshRecursively(mesh, level-1, node, p + Vector3(r.x,r.y,0), r, m_s[s+3]) || morelevels;
	if(m_dim == 2) return morelevels;
	morelevels = createMeshRecursively(mesh, level-1, node, p + Vector3(0,0,r.z), r, m_s[s+4]) || morelevels;
	morelevels = createMeshRecursively(mesh, level-1, node, p + Vector3(r.x,0,r.z), r, m_s[s+5]) || morelevels;
	morelevels = createMeshRecursively(mesh, level-1, node, p + Vector3(0,r.y,r.z), r, m_s[s+6]) || morelevels;
	morelevels = createMeshRecursively(mesh, level-1, node, p + Vector3(r.x,r.y,r.z), r, m_s[s+7]) || morelevels;
	return morelevels;
}
