/**
 * Surrounding.hpp implements a unique surrounding of a pixel in Adaptive.
 * The pixel is located at [0,1[X[0,1[
 * Author: Jukka Räbinä, University of Jyväskylä, 2020.
 */

#ifndef _SURROUNDING_HPP_INCLUDED_
#define _SURROUNDING_HPP_INCLUDED_

#include "../../GFD/Types/Types.hpp"
#include "../../GFD/Types/Buffer.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include <iostream>

using namespace std;

namespace gfd
{

class Surrounding
{
public:
	// constructors
	Surrounding() {
//		insertHyperCube(Vector2(-2,-2), 2.0, 0);
//		insertHyperCube(Vector2(0,0), 1.0, 1);
//		insertHyperCube(Vector2(-0.5,0), 0.5, 0);
	}
	virtual ~Surrounding() {}

	uint getSize() const { return m_data.size(); }

	bool operator<(const Surrounding &sur) const { 
		if(m_data.size() < sur.m_data.size()) return true;
		if(m_data.size() > sur.m_data.size()) return false;
		for(uint i=0; i<m_data.size(); i++) {
			if(m_data[i] < sur.m_data[i]) return true;
			if(m_data[i] > sur.m_data[i]) return false;
		}
		return false;
	}

	bool isInsertible(const double p, const double h) { return isInsertible(1, Vector4(p,0,0,0), h); }
	bool isInsertible(const Vector2 &p, const double h) { return isInsertible(2, Vector4(p,0,0), h); }
	bool isInsertible(const Vector3 &p, const double h) { return isInsertible(3, Vector4(p,0), h); }
	bool isInsertible(const Vector4 &p, const double h) { return isInsertible(4, p, h); }

	bool insertHyperCube(const double p, const double h, const uchar type) { return insertHyperCube(1, Vector4(p,0,0,0), h, type); }
	bool insertHyperCube(const Vector2 &p, const double h, const uchar type) { return insertHyperCube(2, Vector4(p,0,0), h, type); }
	bool insertHyperCube(const Vector3 &p, const double h, const uchar type) { return insertHyperCube(3, Vector4(p,0), h, type); }
	bool insertHyperCube(const Vector4 &p, const double h, const uchar type) { return insertHyperCube(4, p, h, type); }

	void createMesh(BuilderMesh &mesh) {
		uint i = 0;
		uint node = 0;
		while(i < m_data.size()) {
			uint dim;
			Vector4 p;
			double h;
			uchar type;
			i = getNextHyperCube(i, dim, p, h, type);
			mesh.insertNode(p + Vector4(0.25 * h,0.25 *h,0,0), 0.0, node, false);
			mesh.insertNode(p + Vector4(0.75 * h,0.75 *h,0,0), 0.0, node, false);
/*			mesh.insertNode(p, 0.0, node, false);
			mesh.insertNode(p + Vector4(h,0,0,0), 0.0, node, false);
			if(dim == 1) continue;
			mesh.insertNode(p + Vector4(0,h,0,0), 0.0, node, false);
			mesh.insertNode(p + Vector4(h,h,0,0), 0.0, node, false);
			if(dim == 2 && type != 0) mesh.insertNode(p + 0.5 * Vector4(h,h,0,0), 0.0, node, false);
			if(dim == 2) continue;
			mesh.insertNode(p + Vector4(0,0,h,0), 0.0, node, false);
			mesh.insertNode(p + Vector4(h,0,h,0), 0.0, node, false);
			mesh.insertNode(p + Vector4(0,h,h,0), 0.0, node, false);
			mesh.insertNode(p + Vector4(h,h,h,0), 0.0, node, false);
			if(dim == 3 && type != 0) mesh.insertNode(p + 0.5 * Vector4(h,h,h,0), 0.0, node, false);
			if(dim == 3) continue;
			mesh.insertNode(p + Vector4(0,0,0,h), 0.0, node, false);
			mesh.insertNode(p + Vector4(h,0,0,h), 0.0, node, false);
			mesh.insertNode(p + Vector4(0,h,0,h), 0.0, node, false);
			mesh.insertNode(p + Vector4(h,h,0,h), 0.0, node, false);
			mesh.insertNode(p + Vector4(0,0,h,h), 0.0, node, false);
			mesh.insertNode(p + Vector4(h,0,h,h), 0.0, node, false);
			mesh.insertNode(p + Vector4(0,h,h,h), 0.0, node, false);
			mesh.insertNode(p + Vector4(h,h,h,h), 0.0, node, false);
*/		}

		// härö: just testing
//		const double eps = 1e-8;
		for(i=mesh.getFaceSize(); i-->0; ) {
			const Vector4 p = mesh.getFacePosition(i);
			const double r = sqrt(mesh.getRadiusSq(p, mesh.getFaceNodes(i)));
			if(p.x >= r && p.y >= r) continue;
			mesh.removeFace(i);
		}
		for(i=mesh.getEdgeSize(); i-->0; ) {
//			const Vector4 p = mesh.getEdgeAverage(i);
//			if(p.x >= eps && p.x < 1.0 + eps && p.y >= eps && p.y < 1.0 + eps) continue;
			if(!mesh.getEdgeFaces(i).empty()) continue;
			mesh.removeEdge(i);
		}
		for(i=mesh.getNodeSize(); i-->0; ) {
//			const Vector4 p = mesh.getNodePosition(i);
//			if(p.x >= eps && p.x < 1.0 + eps && p.y >= eps && p.y < 1.0 + eps) continue;
			if(!mesh.getNodeEdges(i).empty()) continue;
			mesh.removeNode(i);
		}

/*		// härö: just testing
		const double eps = 1e-8;
		for(i=mesh.getFaceSize(); i-->0; ) {
			const Vector4 p = mesh.getFaceAverage(i);
			if(p.x >= eps && p.x < 1.0 + eps && p.y >= eps && p.y < 1.0 + eps) continue;
			mesh.removeFace(i);
		}
		for(i=mesh.getEdgeSize(); i-->0; ) {
			const Vector4 p = mesh.getEdgeAverage(i);
//			if(p.x >= eps && p.x < 1.0 + eps && p.y >= eps && p.y < 1.0 + eps) continue;
			if(!mesh.getEdgeFaces(i).empty()) continue;
			mesh.removeEdge(i);
		}
		for(i=mesh.getNodeSize(); i-->0; ) {
			const Vector4 p = mesh.getNodePosition(i);
//			if(p.x >= eps && p.x < 1.0 + eps && p.y >= eps && p.y < 1.0 + eps) continue;
			if(!mesh.getNodeEdges(i).empty()) continue;
			mesh.removeNode(i);
		}
*/	}

protected:
	Buffer<uchar> m_data;

	bool isInsertible(const uint dim, const Vector4 &p, const double h) {
		//cout << "try " << p.x << " " << p.y << " " << h << endl; 
		const double minp = -h - 1e-8;
		const double maxp = 1.0 + 1e-8;
		if(p.x < minp || p.x > maxp) return false; // hypercube does not touch the unit cube [0,1[^dim
		if(dim == 1) return true;
		if(p.y < minp || p.y > maxp) return false; // hypercube does not touch the unit cube [0,1[^dim
		//cout << "succeed" << endl; 
		if(dim == 2) return true;
		if(p.z < minp || p.z > maxp) return false; // hypercube does not touch the unit cube [0,1[^dim
		if(dim == 3) return true;
		if(p.t < minp || p.t > maxp) return false; // hypercube does not touch the unit cube [0,1[^dim
		return true;
/*
		const double eps = 1e-8;
		const double dr = 0.5 * h - 0.5;
		const double sr = 0.5 * h + 0.5;
		const double distx = fabs(p.x + dr) - sr;
		if(distx > eps) return false; // hypercube does not touch the unit cube [0,1[^dim
		if(dim == 1) return distx < -eps;
		const double disty = fabs(p.y + dr) - sr;
		if(disty > eps) return false; // hypercube does not touch the unit cube [0,1[^dim
		if(dim == 2) return distx + disty < -eps;
		const double distz = fabs(p.z + dr) - sr;
		if(distz > eps) return false; // hypercube does not touch the unit cube [0,1[^dim
		if(dim == 3) return distx + disty + distz < -eps;
		const double distt = fabs(p.t + dr) - sr;
		if(distt > eps) return false; // hypercube does not touch the unit cube [0,1[^dim
		return distx + disty + distz + distt < -eps;
*/	}
	bool insertHyperCube(const uint dim, const Vector4 &p, const double h, const uchar type) {
		if(!isInsertible(dim, p, h)) return false; // hypercube does not touch the unit cube [0,1[^dim
	
		// discretize hypercube data
		const double eps = 1e-8;
		int level = 0; // size level of the hypercube (-31, ..., 31)
		double hh = h;
		while(hh > 1.4 && level < 31) { level++; hh *= 0.5; }
		while(hh < 0.7 && level > -31) { level--; hh *= 2.0; }
		if(fabs(hh - 1.0) > eps) {
			cout << "Surrounding::insertHyperCube -> size of hypercube is invalid." << endl;
			return false;
		}
		const uint bytes = abs(level) / 8 + 1; // number of bytes (8bits) per position
		uint i = 0;
		uint j, k;
		Buffer<uchar> data(2 + dim * bytes);
		data[i++] = uchar(level + 32) * 4 + uchar(dim - 1);
		uint div = 1;
		for(k=1; k<bytes; k++) div *= 256;
		const double d = (h < 1.0 ? h : 1.0);
		Buffer<double> pj(4);
		pj[0] = p.x; pj[1] = p.y; pj[2] = p.z; pj[3] = p.t;
		for(j=0; j<dim; j++) {
			uint pos = uint((pj[j] + h) / d + 0.5);
			if(fabs(pj[j] + h - pos * d) > eps) {
				cout << "Surrounding::insertHyperCube -> position of hypercube is invalid." << endl;
				return false;
			}
			for(k=div; k>0; k/=256) data[i++] = uchar((pos / k) % 256);
		}
		data[i++] = type;

		// insert new
		uint size = 0;
		for(i=0; i<m_data.size(); i+=size) {
			size = getHyperCubeSize(m_data[i]);
			if(data.size() < size) {
				m_data.insert(data, i);
				return true;
			}
			if(data.size() > size) continue;
			for(j=0; j<size; j++) {
				if(data[j] < m_data[i+j]) {
					m_data.insert(data, i);
					return true;
				}
				if(data[j] > m_data[i+j]) break;
			}
			if(j == size) return false; // identical hypercube already exists
		}
		m_data.combine(data); // insert data at the back of m_data
		return true;
	}


	uint getHyperCubeSize(const uchar data) const {
		const uint dim = uint(data % 4) + 1; // dimension of the hypercube (1, ..., 4)
		const int level = int(data / 4) - 32; // size level of the hypercube (-31, ..., 31)
		const uint bytes = abs(level) / 8 + 1; // number of bytes (8bits) per position
		return 2 + dim * bytes;
	}
	uint getNextHyperCube(uint i, uint &dim, Vector4 &p, double &h, uchar &type) const {
		const uchar size = m_data[i++];
		dim = uint(size % 4) + 1; // dimension of the hypercube (1, ..., 4)
		const int level = int(size / 4) - 32; // size level of the hypercube (-31, ..., 31)
		h = 1.0; // size of the hypercube = 2^level
		for(int j=0; j<level; j++) h *= 2.0;
		for(int j=0; j>level; j--) h *= 0.5;
		const uint bytes = abs(level) / 8 + 1; // number of bytes (8bits) per position
		const double p0 = -h;
		const double p1 = level < 0 ? h : 1.0;
		Buffer<double> pj(4, 0.0);
		for(uint j=0; j<dim; j++) {
			uint pos = 0;
			for(uint k=0; k<bytes; k++) pos = 256 * pos + uint(m_data[i++]);
			pj[j] = p0 + pos * p1;
		}
		p = Vector4(pj[0], pj[1], pj[2], pj[3]); // position of lowest corner of the hypercube
		type = m_data[i++]; // type of the hypercube (0, ..., 255)
		return i;
	}

	void gatherUchar(const uchar val, uint &i) {
		m_data.gather(val, i);
	}
	uchar getUchar(uint &i) const {
		return m_data[i++];
	}
	void gatherUshort(const ushort val, uint &i) {
		if(val >= 255) {
			m_data.gather(255, i);
			m_data.gather(uchar(val / 256), i);
			m_data.gather(uchar(val % 256), i);
		}
		else m_data.gather(uchar(val), i);
	}
	uint getUshort(uint &i) const {
		const uint data0 = uint(m_data[i++]);
		if(data0 == 255) {
			const uint data2 = uint(m_data[i++]);
			const uint data1 = uint(m_data[i++]);
			return data2 * 256 + data1;
		}
		return data0;
	}
	void gatherUint(const uint val, uint &i) {
		if(val >= 255) {
			m_data.gather(255, i);
			m_data.gather(uchar(val / (256 * 256 * 256)), i);
			m_data.gather(uchar((val / (256 * 256)) % 256), i);
			m_data.gather(uchar((val / 256) % 256), i);
			m_data.gather(uchar(val % 256), i);
		}
		else m_data.gather(uchar(val), i);
	}
	uint getUint(uint &i) const {
		const uint data0 = uint(m_data[i++]);
		if(data0 == 255) {
			const uint data4 = uint(m_data[i++]);
			const uint data3 = uint(m_data[i++]);
			const uint data2 = uint(m_data[i++]);
			const uint data1 = uint(m_data[i++]);
			return ((data4 * 256 + data3) * 256 + data2) * 256 + data1;
		}
		return data0;
	}
	

};


}

#endif //_SURROUNDING_HPP_INCLUDED_
