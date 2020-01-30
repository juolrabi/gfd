/**
 * Adaptive.hpp.
 * Author: Jukka Räbinä, University of Jyväskylä, 2020.
 */

#ifndef _ADAPTIVE_HPP_INCLUDED_
#define _ADAPTIVE_HPP_INCLUDED_

#include "Surrounding.hpp"
#include <iostream>
#include <map>

using namespace std;

namespace gfd
{

class Adaptive
{
public:
	// constructors
	Adaptive(const Vector2 &p0, const Vector2 &p1, const double h);
	Adaptive(const Vector3 &p0, const Vector3 &p1, const double h);
	virtual ~Adaptive() {}

	uint findLevel(const Vector3 &p) const;
	void refineSphereBoundary(const Vector3 &p, const double r, const uint level);

/*	void refineSphere(const Vector2 &p, const double r, const double h) {
		refineSphere(p, r, h, m_n, 0, m_p, m_h);
	}
*/
/*	uint findIndex(const Vector2 &p) const {
		const Vector2 d = (p - m_p) / m_h;
		if(isOutOfBounds(d)) return m_size;
		return findIndex(d, m_n, 0);
	}
	Buffer<uint> findTrace(const Vector2 &p) const {
		const Vector2 d = (p - m_p) / m_h;
		Buffer<uint> trace;
		if(!isOutOfBounds(d)) findTrace(d, m_n, 0, trace);
		return trace;
	}
	uint getSize() const { return m_size; }
	void refineSphere(const Vector2 &p, const double r, const double h) {
		refineSphere(p, r, h, m_n, 0, m_p, m_h);
	}
	void fillLevel(const uint level) {
		fillLevel(level, m_n, 0, m_p, m_h);
	}
	Buffer<uint> getTrace(const uint index) const {
		uint i = m_ss;
		int tofind = int(index);
		while(i > 0) { // find tofind from m_s
			if(m_s[--i] == tofind) break;
		}
		const uint base = m_n.x * m_n.y;
		Buffer<uint> trace;
		while(i >= base) {
			const uint remaining = (i - base) % 4;
			trace.push_front(remaining);
			i -= remaining;
			tofind = -int(i);
			while(i > 0) { // find tofind from m_s
				if(m_s[--i] == tofind) break;
			}
		}
		trace.push_front(i);
		return trace;
	}
	void getHyperCube(const Buffer<uint> &trace, Vector2 &p, double &h) const {
		if(trace.empty()) return;
		h = m_h;
		p = m_p + Vector2(h * (trace[0] % m_n.x), h * (trace[0] / m_n.x));
		for(uint i=1; i<trace.size(); i++) {
			h *= 0.5;
			p += Vector2(h * (trace[i] % 2), h * (trace[i] / 2));
		}
	}
	void createSurroundings(Buffer<Surrounding> &stype, Buffer<uint> &slink) const {
		map<Surrounding, uint> smap;
		slink.resize(m_size);
		Uint2 i;
		Buffer<uint> trace(1, 0);
		for(i.y=0; i.y<m_n.y; i.y++) {
			for(i.x=0; i.x<m_n.x; i.x++) {
				createSurroundings(m_s[trace[0]], trace, smap, slink);
				trace[0]++;
			}
		}

		// copy from map to buffer
		stype.resize(smap.size());
		auto it = smap.begin();
		for(uint j=0; j<stype.size(); j++, it++) stype[it->second] = it->first;
	}
	void createSurroundings(const int s, Buffer<uint> &trace, map<Surrounding, uint> &smap, Buffer<uint> &slink) const {
		if(s >= 0) { // constant pixel found
			const Surrounding sur = getSurrounding(trace);
			const uint size = smap.size();
			auto it = smap.find(sur);
			if(it == smap.end()) {
				smap.insert(make_pair(sur, size));
				slink[s] = size;
			}
			else slink[s] = it->second;
			return;
		}
		const uint ii = uint(-s);
		trace.push_back(0);
		createSurroundings(m_s[ii], trace, smap, slink);
		trace.back() = 1;
		createSurroundings(m_s[ii+1], trace, smap, slink);
		trace.back() = 2;
		createSurroundings(m_s[ii+2], trace, smap, slink);
		trace.back() = 3;
		createSurroundings(m_s[ii+3], trace, smap, slink);
		trace.pop_back();
	}
	Surrounding getSurrounding(const Buffer<uint> &trace) const {
		Vector2 p(0,0);
		double h = 1.0;
		getHyperCube(trace, p, h);
		Surrounding surr;

		const double h2 = 0.5 * h;
		const Vector2 d = (p + Vector2(h2, h2) - m_p) / m_h;
		const Uint2 i(d.x, d.y);
		const Uint2 i0((i.x > 0 ? i.x - 1 : 0), (i.y > 0 ? i.y - 1 : 0));
		const Uint2 i1((i.x + 1 < m_n.x ? i.x + 1 : m_n.x - 1), (i.y + 1 < m_n.y ? i.y + 1 : m_n.y - 1));
		Uint2 j;
		for(j.y=i0.y; j.y<=i1.y; j.y++) {
			for(j.x=i0.x; j.x<=i1.x; j.x++) {
				const uint jj = j.y * m_n.x + j.x; // pixel index
				gatherSurrounding((m_p + Vector2(j.x * m_h, j.y * m_h) - p) / h, m_h / h, jj, surr);
			}
		}
		return surr;

	}
	void gatherSurrounding(const Vector2 &p, const double h, const uint ii, Surrounding &surr) const {
		if(m_s[ii] >= 0) { // constant pixel found
			surr.insertHyperCube(p, h, 1);
			return;
		}
		if(!surr.isInsertible(p, h)) return;
		const uint ii2 = uint(-m_s[ii]);
		const double h2 = 0.5 * h;
		gatherSurrounding(p, h2, ii2, surr);
		gatherSurrounding(p + Vector2(h2,0), h2, ii2+1, surr);
		gatherSurrounding(p + Vector2(0,h2), h2, ii2+2, surr);
		gatherSurrounding(p + Vector2(h2,h2), h2, ii2+3, surr);
	}

	void createMesh(Mesh &mesh, const Buffer<BuilderMesh> &smesh, const Buffer<uint> &slink) const {
		Uint2 i;
		uint ii = 0;
		for(i.y=0; i.y<m_n.y; i.y++) {
			for(i.x=0; i.x<m_n.x; i.x++) {
				createMesh(m_p + Vector2(i.x * m_h, i.y * m_h), m_h, m_s[ii++], mesh, smesh, slink);
			}
		}
	}
	void createMesh(const Vector2 &p, const double h, const int s, Mesh &mesh, const Buffer<BuilderMesh> &smesh, const Buffer<uint> &slink) const {
		if(s >= 0) { // constant pixel found
			// insert elements of smesh[slink[s]];
			const uint node0 = mesh.getNodeSize();
			const BuilderMesh &smeshs = smesh[slink[s]];
			for(uint i=0; i<smeshs.getNodeSize(); i++) {
				mesh.addNode(Vector4(p,0,0) + h * smeshs.getNodePosition(i));
			}
			const uint edge0 = mesh.getEdgeSize();
			for(uint i=0; i<smeshs.getEdgeSize(); i++) {
				const Buffer<uint> &par = smeshs.getEdgeNodes(i);
				mesh.addEdge(par[0] + node0, par[1] + node0);
			}
			for(uint i=0; i<smeshs.getFaceSize(); i++) {
				Buffer<uint> par = smeshs.getFaceEdges(i);
				for(uint j=0; j<par.size(); j++) par[j] += edge0;
				mesh.addFace(par);
			}
			return;
		}
		const uint ii = uint(-s);
		const double h2 = 0.5 * h;
		createMesh(p, h2, m_s[ii], mesh, smesh, slink);
		createMesh(p + Vector2(h2,0), h2, m_s[ii+1], mesh, smesh, slink);
		createMesh(p + Vector2(0,h2), h2, m_s[ii+2], mesh, smesh, slink);
		createMesh(p + Vector2(h2,h2), h2, m_s[ii+3], mesh, smesh, slink);
	}

	void createMesh(BuilderMesh &mesh) const {
		// create base grid
		const Uint2 n((m_n.x < 3 ? m_n.x : 3), (m_n.y < 3 ? m_n.y : 3));
		mesh.createGrid(Vector4(m_p,0,0), Vector4(m_p+m_h*Vector2(n.x,n.y),0,0), m_h);

		Uint2 i;
		uint node = 0;
		for(i.y=0; i.y<n.y; i.y++) {
			for(i.x=0; i.x<n.x; i.x++) {
				node = mesh.insertNode(Vector4(m_p + m_h * Vector2(i.x + 0.5, i.y + 0.5), 0,0), 0.0, node, true);
			}
		}
		if(m_n.x > n.x) mesh.repeatMiddle(Vector4(m_p.x + 0.99999 * m_h, 0,0,0), Vector4(m_h,0,0,0), m_n.x - n.x);
		if(m_n.y > n.y) mesh.repeatMiddle(Vector4(0,m_p.y + 0.99999 * m_h, 0,0), Vector4(0,m_h,0,0), m_n.y - n.y);
		
		bool morelevels = true;
		Buffer<bool> more(m_n.x * m_n.y, true);
		for(uint level=0; morelevels; level++) {
			morelevels = false;
			Uint2 i;
			uint ii = 0;
			for(i.y=0; i.y<m_n.y; i.y++) {
				for(i.x=0; i.x<m_n.x; i.x++, ii++) {
					if(!more[ii]) continue;
					const Vector2 p = m_p + Vector2(m_h * i.x, m_h * i.y);
					more[ii] = insertLevelNodes(mesh, level, node, p, m_h, m_s[ii]);
					if(more[ii]) morelevels = true;
				}
			}
			cout << "iteration " << level << " " << morelevels << endl;
		}
	}
*/
	void splitPixel(const uint ii);
	void createMesh(BuilderMesh &mesh) const;

protected:
	uint m_dim;
	Uint3 m_n;
	uint m_ss;
	Buffer<uint> m_s;
	Vector3 m_p;
	Vector3 m_h;



	bool isOutOfBounds(const Vector3 &d) const;
	uint findLevelRecursive(const Vector3 &p, uint ii, const uint level) const;
	uint getBlockSize() const;	
	Vector3 getFlat(const Vector3 &v) const;
	void refineSphereBoundaryRecursive(const Vector3 &p, const Vector3 &r, const uint ii, const uint level);
	bool createMeshRecursively(BuilderMesh &mesh, const uint level, uint &node, const Vector3 &p, const Vector3 &h, const int s) const;

/*	uint findIndex(const Vector2 &d, const Uint2 &n, uint ii) const {
		const Uint2 i(d.x, d.y);
		ii += i.y * n.x + i.x; // pixel index
		if(m_s[ii] >= 0) return uint(m_s[ii]); // constant pixel found
		return findIndex(2.0 * Vector2(d.x - double(i.x), d.y - double(i.y)), Uint2(2,2), uint(-m_s[ii]));
	}
	void findTrace(const Vector2 &d, const Uint2 &n, uint ii, Buffer<uint> &trace) const {
		const Uint2 i(d.x, d.y);
		const uint slot = i.y * n.x + i.x;
		trace.push_back(slot);
		ii += slot; // pixel index
		if(m_s[ii] < 0) findTrace(2.0 * Vector2(d.x - double(i.x), d.y - double(i.y)), Uint2(2,2), uint(-m_s[ii]), trace);
	}

	void refineSphere(const Vector2 &sp, const double sr, const double sh, const Uint2 &n, uint ii, const Vector2 &p, const double h) {
		if(sh > h) return;
		Uint2 i;
		for(i.y=0; i.y<n.y; i.y++) {
			for(i.x=0; i.x<n.x; i.x++, ii++) {
				const Vector2 ip = p + h * Vector2(i.x, i.y);
				const double dlen = (ip + 0.5 * Vector2(h,h) - sp).len();
				if(fabs(dlen - sr) > 0.9 * h) continue;
				if(m_s[ii] >= 0) splitPixel(ii);
				refineSphere(sp, sr, sh, Uint2(2,2), uint(-m_s[ii]), ip, 0.5 * h);
			}
		}
	}
	void fillLevel(const uint level, const Uint2 &n, uint ii, const Vector2 &p, const double h) {
		if(level <= 1) return;
		Uint2 i;
		for(i.y=0; i.y<n.y; i.y++) {
			for(i.x=0; i.x<n.x; i.x++, ii++) {
				const Vector2 ip = p + h * Vector2(i.x, i.y);
				if(m_s[ii] >= 0) splitPixel(ii);
				fillLevel(level - 1, Uint2(2,2), uint(-m_s[ii]), ip, 0.5 * h);
			}
		}
	}

	void splitPixel(const uint ii) {
		if(m_s[ii] < 0) return;
		const int oldt = m_s[ii];
		m_s[ii] = -int(m_ss);
		m_s.gather(oldt, m_ss);
		for(uint i=0; i<3; i++, m_size++) m_s.gather(m_size, m_ss);
	}

	bool insertLevelNodes(BuilderMesh &mesh, const uint level, uint &node, const Vector2 &p, const double h, const int s) const {
		const double r = 0.5 * h;
		if(s >= 0) { // constant pixel found
			return false;
		}

		// insert nodes
		if(level == 0) {
			const Vector2 o = p + Vector2(r, r);
			const double rr = 0.5 * r;
			node = mesh.insertNode(Vector4(o - Vector2(r,0), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o - Vector2(0,r), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(r,0), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(0,r), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(rr,rr), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(rr,-rr), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(-rr,-rr), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(-rr,rr), 0,0), 0.0, node, false);
			return true;
		}

		const uint ii = uint(-s);
		bool morelevels = false;
		morelevels = insertLevelNodes(mesh, level-1, node, p, r, m_s[ii]) || morelevels;
		morelevels = insertLevelNodes(mesh, level-1, node, p + Vector2(r,0), r, m_s[ii+1]) || morelevels;
		morelevels = insertLevelNodes(mesh, level-1, node, p + Vector2(0,r), r, m_s[ii+2]) || morelevels;
		morelevels = insertLevelNodes(mesh, level-1, node, p + Vector2(r,r), r, m_s[ii+3]) || morelevels;
		return morelevels;
	}
*/
};

class Adaptive2
{
public:
	// constructors
	Adaptive2(const Uint2 &n, const Vector2 &p, const double h) {
		m_n = n;
		m_p = p;
		m_h = h;
		m_s.resize(m_n.x * m_n.y);
		for(m_ss=0; m_ss<m_s.size(); m_ss++) m_s[m_ss] = int(m_ss);
		m_size = m_ss;
	}
	virtual ~Adaptive2() {}
	uint findLevel(const Vector2 &p) const {
		const Vector2 d = (p - m_p) / m_h;
		if(isOutOfBounds(d)) return 0;
		return findLevel(d, m_n, 0, 1);
	}
	uint findIndex(const Vector2 &p) const {
		const Vector2 d = (p - m_p) / m_h;
		if(isOutOfBounds(d)) return m_size;
		return findIndex(d, m_n, 0);
	}
	Buffer<uint> findTrace(const Vector2 &p) const {
		const Vector2 d = (p - m_p) / m_h;
		Buffer<uint> trace;
		if(!isOutOfBounds(d)) findTrace(d, m_n, 0, trace);
		return trace;
	}
	uint getSize() const { return m_size; }
	void refineSphere(const Vector2 &p, const double r, const double h) {
		refineSphere(p, r, h, m_n, 0, m_p, m_h);
	}
	void fillLevel(const uint level) {
		fillLevel(level, m_n, 0, m_p, m_h);
	}
	Buffer<uint> getTrace(const uint index) const {
		uint i = m_ss;
		int tofind = int(index);
		while(i > 0) { // find tofind from m_s
			if(m_s[--i] == tofind) break;
		}
		const uint base = m_n.x * m_n.y;
		Buffer<uint> trace;
		while(i >= base) {
			const uint remaining = (i - base) % 4;
			trace.push_front(remaining);
			i -= remaining;
			tofind = -int(i);
			while(i > 0) { // find tofind from m_s
				if(m_s[--i] == tofind) break;
			}
		}
		trace.push_front(i);
		return trace;
	}
	void getHyperCube(const Buffer<uint> &trace, Vector2 &p, double &h) const {
		if(trace.empty()) return;
		h = m_h;
		p = m_p + Vector2(h * (trace[0] % m_n.x), h * (trace[0] / m_n.x));
		for(uint i=1; i<trace.size(); i++) {
			h *= 0.5;
			p += Vector2(h * (trace[i] % 2), h * (trace[i] / 2));
		}
	}
	void createSurroundings(Buffer<Surrounding> &stype, Buffer<uint> &slink) const {
		map<Surrounding, uint> smap;
		slink.resize(m_size);
/*		Buffer<uint> type(m_size);

		Surrounding sur1;
//		surrounding sur2;
		mp.insert(make_pair(sur1, 0));
		map<Surrounding, uint>::iterator it = mp.find(sur1);
		if(it != mp.end()) cout << it->second << endl;
		//if(it != NULL) cout << "jee" << endl;
		//else cout << "buu" << endl;
		createSurroundings(m_p, m_n, 0, )
*/
		Uint2 i;
		Buffer<uint> trace(1, 0);
		for(i.y=0; i.y<m_n.y; i.y++) {
			for(i.x=0; i.x<m_n.x; i.x++) {
				createSurroundings(m_s[trace[0]], trace, smap, slink);
				trace[0]++;
			}
		}

		// copy from map to buffer
		stype.resize(smap.size());
		auto it = smap.begin();
		for(uint j=0; j<stype.size(); j++, it++) stype[it->second] = it->first;
	}
	void createSurroundings(const int s, Buffer<uint> &trace, map<Surrounding, uint> &smap, Buffer<uint> &slink) const {
		if(s >= 0) { // constant pixel found
			const Surrounding sur = getSurrounding(trace);
			const uint size = smap.size();
			auto it = smap.find(sur);
			if(it == smap.end()) {
				smap.insert(make_pair(sur, size));
				slink[s] = size;
			}
			else slink[s] = it->second;
			return;
		}
		const uint ii = uint(-s);
		trace.push_back(0);
		createSurroundings(m_s[ii], trace, smap, slink);
		trace.back() = 1;
		createSurroundings(m_s[ii+1], trace, smap, slink);
		trace.back() = 2;
		createSurroundings(m_s[ii+2], trace, smap, slink);
		trace.back() = 3;
		createSurroundings(m_s[ii+3], trace, smap, slink);
		trace.pop_back();
	}
	Surrounding getSurrounding(const Buffer<uint> &trace) const {
		Vector2 p(0,0);
		double h = 1.0;
		getHyperCube(trace, p, h);
		Surrounding surr;

		const double h2 = 0.5 * h;
		const Vector2 d = (p + Vector2(h2, h2) - m_p) / m_h;
		const Uint2 i(d.x, d.y);
		const Uint2 i0((i.x > 0 ? i.x - 1 : 0), (i.y > 0 ? i.y - 1 : 0));
		const Uint2 i1((i.x + 1 < m_n.x ? i.x + 1 : m_n.x - 1), (i.y + 1 < m_n.y ? i.y + 1 : m_n.y - 1));
		Uint2 j;
		for(j.y=i0.y; j.y<=i1.y; j.y++) {
			for(j.x=i0.x; j.x<=i1.x; j.x++) {
				const uint jj = j.y * m_n.x + j.x; // pixel index
				gatherSurrounding((m_p + Vector2(j.x * m_h, j.y * m_h) - p) / h, m_h / h, jj, surr);
			}
		}
		return surr;

	}
	void gatherSurrounding(const Vector2 &p, const double h, const uint ii, Surrounding &surr) const {
		if(m_s[ii] >= 0) { // constant pixel found
			surr.insertHyperCube(p, h, 1);
			return;
		}
		if(!surr.isInsertible(p, h)) return;
		const uint ii2 = uint(-m_s[ii]);
		const double h2 = 0.5 * h;
		gatherSurrounding(p, h2, ii2, surr);
		gatherSurrounding(p + Vector2(h2,0), h2, ii2+1, surr);
		gatherSurrounding(p + Vector2(0,h2), h2, ii2+2, surr);
		gatherSurrounding(p + Vector2(h2,h2), h2, ii2+3, surr);
	}

	void createMesh(Mesh &mesh, const Buffer<BuilderMesh> &smesh, const Buffer<uint> &slink) const {
		Uint2 i;
		uint ii = 0;
		for(i.y=0; i.y<m_n.y; i.y++) {
			for(i.x=0; i.x<m_n.x; i.x++) {
				createMesh(m_p + Vector2(i.x * m_h, i.y * m_h), m_h, m_s[ii++], mesh, smesh, slink);
			}
		}
	}
	void createMesh(const Vector2 &p, const double h, const int s, Mesh &mesh, const Buffer<BuilderMesh> &smesh, const Buffer<uint> &slink) const {
		if(s >= 0) { // constant pixel found
			// insert elements of smesh[slink[s]];
			const uint node0 = mesh.getNodeSize();
			const BuilderMesh &smeshs = smesh[slink[s]];
			for(uint i=0; i<smeshs.getNodeSize(); i++) {
				mesh.addNode(Vector4(p,0,0) + h * smeshs.getNodePosition(i));
			}
			const uint edge0 = mesh.getEdgeSize();
			for(uint i=0; i<smeshs.getEdgeSize(); i++) {
				const Buffer<uint> &par = smeshs.getEdgeNodes(i);
				mesh.addEdge(par[0] + node0, par[1] + node0);
			}
			for(uint i=0; i<smeshs.getFaceSize(); i++) {
				Buffer<uint> par = smeshs.getFaceEdges(i);
				for(uint j=0; j<par.size(); j++) par[j] += edge0;
				mesh.addFace(par);
			}
			return;
		}
		const uint ii = uint(-s);
		const double h2 = 0.5 * h;
		createMesh(p, h2, m_s[ii], mesh, smesh, slink);
		createMesh(p + Vector2(h2,0), h2, m_s[ii+1], mesh, smesh, slink);
		createMesh(p + Vector2(0,h2), h2, m_s[ii+2], mesh, smesh, slink);
		createMesh(p + Vector2(h2,h2), h2, m_s[ii+3], mesh, smesh, slink);
	}

	void createMesh(BuilderMesh &mesh) const {
		// create base grid
		const Uint2 n((m_n.x < 3 ? m_n.x : 3), (m_n.y < 3 ? m_n.y : 3));
		mesh.createGrid(Vector4(m_p,0,0), Vector4(m_p+m_h*Vector2(n.x,n.y),0,0), m_h);

		Uint2 i;
		uint node = 0;
		for(i.y=0; i.y<n.y; i.y++) {
			for(i.x=0; i.x<n.x; i.x++) {
				node = mesh.insertNode(Vector4(m_p + m_h * Vector2(i.x + 0.5, i.y + 0.5), 0,0), 0.0, node, true);
			}
		}
		if(m_n.x > n.x) mesh.repeatMiddle(Vector4(m_p.x + 0.99999 * m_h, 0,0,0), Vector4(m_h,0,0,0), m_n.x - n.x);
		if(m_n.y > n.y) mesh.repeatMiddle(Vector4(0,m_p.y + 0.99999 * m_h, 0,0), Vector4(0,m_h,0,0), m_n.y - n.y);
		
/*		mesh.createGrid(Vector4(m_p, 0,0), Vector4(m_p+m_h*Vector2(m_n.x,m_n.y), 0,0), m_h);
		uint node = 0;

			Uint2 i;
			uint ii = 0;
			for(i.y=0; i.y<m_n.y; i.y++) {
				for(i.x=0; i.x<m_n.x; i.x++) {
					node = mesh.insertNode(Vector4(m_p + Vector2(m_h * (i.x+0.5), m_h * (i.y+0.5)), 0,0), 0.0, node, true);
				}
			}
			cout << "vaihe 1" << endl;
*/
		bool morelevels = true;
		Buffer<bool> more(m_n.x * m_n.y, true);
		for(uint level=0; morelevels; level++) {
			morelevels = false;
			Uint2 i;
			uint ii = 0;
			for(i.y=0; i.y<m_n.y; i.y++) {
				for(i.x=0; i.x<m_n.x; i.x++, ii++) {
					if(!more[ii]) continue;
					const Vector2 p = m_p + Vector2(m_h * i.x, m_h * i.y);
					more[ii] = insertLevelNodes(mesh, level, node, p, m_h, m_s[ii]);
					if(more[ii]) morelevels = true;
				}
			}
			cout << "iteration " << level << " " << morelevels << endl;
		}
/*		Uint2 i;
		uint ii = 0;
		for(i.y=0; i.y<m_n.y; i.y++) {
			for(i.x=0; i.x<m_n.x; i.x++) {
				const Vector2 p = m_p + Vector2(m_h * i.x, m_h * i.y);
				insertConstantNodes(mesh, node, p, m_h, m_s[ii++]);
			}
		}
*/	}

protected:
	Uint2 m_n;
	uint m_ss;
	Buffer<int> m_s;
	uint m_size;
	Vector2 m_p;
	double m_h;

	bool isOutOfBounds(const Vector2 &d) const {
		if(d.x < 0.0 || d.y < 0.0) return true; // out of bounds in negative direction
		if(d.x >= double(m_n.x) || d.y >= double(m_n.y)) return true; // out of bounds in positive direction
		return false;
	}
	uint findLevel(const Vector2 &d, const Uint2 &n, uint ii, const uint level) const {
		const Uint2 i(d.x, d.y);
		ii += i.y * n.x + i.x; // pixel index
		if(m_s[ii] >= 0) return level; // constant pixel found
		return findLevel(2.0 * Vector2(d.x - double(i.x), d.y - double(i.y)), Uint2(2,2), uint(-m_s[ii]), level + 1);
	}
	uint findIndex(const Vector2 &d, const Uint2 &n, uint ii) const {
		const Uint2 i(d.x, d.y);
		ii += i.y * n.x + i.x; // pixel index
		if(m_s[ii] >= 0) return uint(m_s[ii]); // constant pixel found
		return findIndex(2.0 * Vector2(d.x - double(i.x), d.y - double(i.y)), Uint2(2,2), uint(-m_s[ii]));
	}
	void findTrace(const Vector2 &d, const Uint2 &n, uint ii, Buffer<uint> &trace) const {
		const Uint2 i(d.x, d.y);
		const uint slot = i.y * n.x + i.x;
		trace.push_back(slot);
		ii += slot; // pixel index
		if(m_s[ii] < 0) findTrace(2.0 * Vector2(d.x - double(i.x), d.y - double(i.y)), Uint2(2,2), uint(-m_s[ii]), trace);
	}

	void refineSphere(const Vector2 &sp, const double sr, const double sh, const Uint2 &n, uint ii, const Vector2 &p, const double h) {
		if(sh > h) return;
		Uint2 i;
		for(i.y=0; i.y<n.y; i.y++) {
			for(i.x=0; i.x<n.x; i.x++, ii++) {
				const Vector2 ip = p + h * Vector2(i.x, i.y);
				const double dlen = (ip + 0.5 * Vector2(h,h) - sp).len();
				if(fabs(dlen - sr) > 0.9 * h) continue;
				if(m_s[ii] >= 0) splitPixel(ii);
				refineSphere(sp, sr, sh, Uint2(2,2), uint(-m_s[ii]), ip, 0.5 * h);
			}
		}
	}
	void fillLevel(const uint level, const Uint2 &n, uint ii, const Vector2 &p, const double h) {
		if(level <= 1) return;
		Uint2 i;
		for(i.y=0; i.y<n.y; i.y++) {
			for(i.x=0; i.x<n.x; i.x++, ii++) {
				const Vector2 ip = p + h * Vector2(i.x, i.y);
				if(m_s[ii] >= 0) splitPixel(ii);
				fillLevel(level - 1, Uint2(2,2), uint(-m_s[ii]), ip, 0.5 * h);
			}
		}
	}

	void splitPixel(const uint ii) {
		if(m_s[ii] < 0) return;
		const int oldt = m_s[ii];
		m_s[ii] = -int(m_ss);
		m_s.gather(oldt, m_ss);
		for(uint i=0; i<3; i++, m_size++) m_s.gather(m_size, m_ss);
	}

	bool insertLevelNodes(BuilderMesh &mesh, const uint level, uint &node, const Vector2 &p, const double h, const int s) const {
		const double r = 0.5 * h;
		if(s >= 0) { // constant pixel found
			//if(level == 0) node = mesh.insertNode(Vector4(p + Vector2(r, r), 0,0), 0.0, node, false);
			return false;
		}

/*		// insert nodes
		if(level == 0) {
			const Vector2 o = p + Vector2(r, r);
			node = mesh.insertNode(Vector4(o, 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o - Vector2(r,0), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o - Vector2(0,r), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(r,0), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(0,r), 0,0), 0.0, node, false);
			return true;
		}
*/
		// insert nodes
		if(level == 0) {
			const Vector2 o = p + Vector2(r, r);
			const double rr = 0.5 * r;
			node = mesh.insertNode(Vector4(o - Vector2(r,0), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o - Vector2(0,r), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(r,0), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(0,r), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(rr,rr), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(rr,-rr), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(-rr,-rr), 0,0), 0.0, node, false);
			node = mesh.insertNode(Vector4(o + Vector2(-rr,rr), 0,0), 0.0, node, false);
			return true;
		}

		const uint ii = uint(-s);
		bool morelevels = false;
		morelevels = insertLevelNodes(mesh, level-1, node, p, r, m_s[ii]) || morelevels;
		morelevels = insertLevelNodes(mesh, level-1, node, p + Vector2(r,0), r, m_s[ii+1]) || morelevels;
		morelevels = insertLevelNodes(mesh, level-1, node, p + Vector2(0,r), r, m_s[ii+2]) || morelevels;
		morelevels = insertLevelNodes(mesh, level-1, node, p + Vector2(r,r), r, m_s[ii+3]) || morelevels;
		return morelevels;
	}

};

class Adaptive3
{
public:
	// constructors
	Adaptive3(const Uint3 &n, const Vector3 &p, const double h) {
		m_n = n;
		m_p = p;
		m_h = h;
		m_s.resize(m_n.x * m_n.y * m_n.z);
		for(m_ss=0; m_ss<m_s.size(); m_ss++) m_s[m_ss] = int(m_ss);
		m_size = m_ss;
	}
	virtual ~Adaptive3() {}
	uint findLevel(const Vector3 &p) const {
		const Vector3 d = (p - m_p) / m_h;
		if(isOutOfBounds(d)) return 0;
		return findLevel(d, m_n, 0, 1);
	}
	uint findIndex(const Vector3 &p) const {
		const Vector3 d = (p - m_p) / m_h;
		if(isOutOfBounds(d)) return m_size;
		return findIndex(d, m_n, 0);
	}
	uint getSize() const { return m_size; }
	void refineSphere(const Vector3 &p, const double r, const double h) {
		refineSphere(p, r, h, m_n, 0, m_p, m_h);
	}
	void fillLevel(const uint level) {
		fillLevel(level, m_n, 0, m_p, m_h);
	}


protected:
	Uint3 m_n;
	uint m_ss;
	Buffer<int> m_s;
	uint m_size;
	Vector3 m_p;
	double m_h;

	bool isOutOfBounds(const Vector3 &d) const {
		if(d.x < 0.0 || d.y < 0.0 || d.z < 0.0) return true; // out of bounds in negative direction
		if(d.x >= double(m_n.x) || d.y >= double(m_n.y) || d.z >= double(m_n.z)) return true; // out of bounds in positive direction
		return false;
	}
	uint findLevel(const Vector3 &d, const Uint3 &n, uint ii, const uint level) const {
		const Uint3 i(d.x, d.y, d.z);
		ii += (i.z * n.y + i.y) * n.x + i.x; // pixel index
		if(m_s[ii] >= 0) return level; // constant pixel found
		return findLevel(2.0 * Vector3(d.x - double(i.x), d.y - double(i.y), d.z - double(i.z)), Uint3(2,2,2), uint(-m_s[ii]), level + 1);
	}
	uint findIndex(const Vector3 &d, const Uint3 &n, uint ii) const {
		const Uint3 i(d.x, d.y, d.z);
		ii += (i.z * n.y + i.y) * n.x + i.x; // pixel index
		if(m_s[ii] >= 0) return uint(m_s[ii]); // constant pixel found
		return findIndex(2.0 * Vector3(d.x - double(i.x), d.y - double(i.y), d.z - double(i.z)), Uint3(2,2,2), uint(-m_s[ii]));
	}

	void refineSphere(const Vector3 &sp, const double sr, const double sh, const Uint3 &n, uint ii, const Vector3 &p, const double h) {
		if(sh > h) return;
		Uint3 i;
		for(i.z=0; i.z<n.z; i.z++) {
			for(i.y=0; i.y<n.y; i.y++) {
				for(i.x=0; i.x<n.x; i.x++, ii++) {
					const Vector3 ip = p + h * Vector3(i.x, i.y, i.z);
					const double dlen = (ip + 0.5 * Vector3(h,h,h) - sp).len();
					if(fabs(dlen - sr) > 0.9 * h) continue;
					if(m_s[ii] >= 0) splitPixel(ii);
					refineSphere(sp, sr, sh, Uint3(2,2,2), uint(-m_s[ii]), ip, 0.5 * h);
				}
			}
		}
	}
	void fillLevel(const uint level, const Uint3 &n, uint ii, const Vector3 &p, const double h) {
		if(level <= 1) return;
		Uint3 i;
		for(i.z=0; i.z<n.z; i.z++) {
			for(i.y=0; i.y<n.y; i.y++) {
				for(i.x=0; i.x<n.x; i.x++, ii++) {
					const Vector3 ip = p + h * Vector3(i.x, i.y, i.z);
					if(m_s[ii] >= 0) splitPixel(ii);
					fillLevel(level - 1, Uint3(2,2,2), uint(-m_s[ii]), ip, 0.5 * h);
				}
			}
		}
	}

	void splitPixel(const uint ii) {
		if(m_s[ii] < 0) return;
		const int oldt = m_s[ii];
		m_s[ii] = -int(m_ss);
		m_s.gather(oldt, m_ss);
		for(uint i=0; i<7; i++, m_size++) m_s.gather(m_size, m_ss);
	}

};


}

#endif //_ADAPTIVE_HPP_INCLUDED_
