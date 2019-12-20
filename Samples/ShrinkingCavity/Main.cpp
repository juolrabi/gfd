/**
 * Program for computing (3+1)-dimensional wave in a shrinking cavity.
 * Idea of the solution is following:
 * - create (3+1)-dimensional space-time mesh.
 * - consider discrete forms of grade 1 (values located at edges).
 * - initialize wave at t=0.
 * - solve field finding gradient zero dF=0 (using Minkowski-metric).
 * - in this occasion, dF=0 leads to explicit and causal system.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include "../../GFD/Types/Types.hpp"
#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Types/Buffer.hpp"
#include <iostream>
#include <chrono>

using namespace std;
using namespace gfd;

Vector4 getField(const Mesh &mesh, const Buffer<double> &val, const uint node) {
	const Buffer<uint> &par = mesh.getNodeEdges(node);
	SymMatrix4 A(ZEROSYMMATRIX4);
//	A += Vector4(0,0,1,0).outerProduct(); // this line is needed in 3-dimensional task
	Vector4 b(0,0,0,0);
	for(uint i=0; i<par.size(); i++) {
		const Vector4 v = mesh.getEdgeVector(par[i]);
		const double w = 1.0 / v.lensq();
		A += w * v.outerProduct();
		b += w * v * val[par[i]];
	}
	return A.inverse() * b;
}

bool saveVolumeMap(const std::string &path, const Buffer<ushort> &vol, const uint xsize, const uint ysize, const uint zsize, const Vector3 &h)
{
	Text rawpath;
	rawpath << path << ".raw";

	// save raw
	std::ofstream fs(rawpath.str().c_str(), std::ios_base::binary | std::ios::trunc);
	if(fs.fail()) return false;
	fs.write((char*)&vol[0], 2 * xsize * ysize * zsize);
	fs.close();

	// save header
	Text text;

	text <<	"ObjectType              = Image" << std::endl;
	text <<	"NDims                   = 3" << std::endl;
	text <<	"BinaryData              = True" << std::endl;
	text <<	"CompressedData          = False" << std::endl;
	text <<	"BinaryDataByteOrderMSB  = False" << std::endl;
	text <<	"TransformMatrix         = 1 0 0 0 1 0 0 0 1" << std::endl;
	text <<	"Offset                  = " << -0.5 * xsize * h.x << " " << -0.5 * ysize * h.y << " " << -0.5 * zsize * h.z << std::endl;
	text <<	"CenterOfRotation        = 0 0 0" << std::endl;
	text <<	"DimSize                 = " << xsize << " " << ysize << " " << zsize << std::endl;
	text <<	"ElementSpacing          = " << h.x << " " << h.y << " " << h.z << std::endl;
	text <<	"ElementNumberOfChannels = 1" << std::endl;
	text <<	"ElementType             = MET_USHORT" << std::endl;
	text <<	"ElementDataFile         = " << rawpath.str() << std::endl;
	text.save(path);
	return true;
}

Buffer<uint> findEdges(const uint n0, const uint n1, const uint n2, BuilderMesh &mesh) {
	Buffer<uint> e(3);
	e[0] = mesh.findEdge(n0, n1);
	e[1] = mesh.findEdge(n1, n2);
	e[2] = mesh.findEdge(n2, n0);
	return e;
}
Buffer<uint> findFaces(const uint n0, const uint n1, const uint n2, const uint n3, BuilderMesh &mesh) {
	Buffer<uint> f(4);
	f[0] = mesh.findFace(findEdges(n0, n1, n2, mesh));
	f[1] = mesh.findFace(findEdges(n0, n1, n3, mesh));
	f[2] = mesh.findFace(findEdges(n0, n2, n3, mesh));
	f[3] = mesh.findFace(findEdges(n1, n2, n3, mesh));
	return f;
}
Buffer<uint> findBodies(const uint n0, const uint n1, const uint n2, const uint n3, const uint n4, BuilderMesh &mesh) {
	Buffer<uint> b(5);
	b[0] = mesh.findBody(findFaces(n0, n1, n2, n3, mesh));
	b[1] = mesh.findBody(findFaces(n0, n1, n2, n4, mesh));
	b[2] = mesh.findBody(findFaces(n0, n1, n3, n4, mesh));
	b[3] = mesh.findBody(findFaces(n0, n2, n3, n4, mesh));
	b[4] = mesh.findBody(findFaces(n1, n2, n3, n4, mesh));
	return b;
}
uint createSpaceTimeMesh(const double h, const double dtime, const uint steps, const double slope, const bool async, BuilderMesh &mesh) {
	uint i, j, k, l, m;
	const double part = 0.5;
	const double N = 1.0 + 1e-5;
	const double n = part + 1e-5;

/*// The first (commented) block is for 3-dimensional task. The second is for 4-dimensional task.
	BuilderMesh basemesh(3);
	const double SQRT3_4 = sqrt(0.75);
	basemesh.createTriangleGrid(Vector2(-1,-SQRT3_4), Vector2(1,SQRT3_4), h, false);
	const Vector2 e0 = Vector2(0.0, SQRT3_4) / 0.75;
	const Vector2 e1 = Vector2(0.75, 0.5 * SQRT3_4) / 0.75;
	const Vector2 e2 = Vector2(0.75, -0.5 * SQRT3_4) / 0.75;
	for(i=basemesh.getNodeSize(); i-->0; ) {
		const Vector2 p = basemesh.getNodePosition2(i);
		if(e0.dot(p) > N || e0.dot(p) < -N) basemesh.removeNode(i);
		else if(e1.dot(p) > N || e1.dot(p) < -N) basemesh.removeNode(i);
		else if(e2.dot(p) > N || e2.dot(p) < -N) basemesh.removeNode(i);
	}
	const uint nodes = basemesh.getNodeSize();
	for(i=0; i<nodes; i++) {
		const Vector2 x = basemesh.getNodePosition2(i);
		basemesh.setNodeFlag(i, uint((3.0 + x.x + 1.5 * x.y / sqrt(0.75)) / h + 0.5) % 3);
	}
/*/
	BuilderMesh basemesh(3);
	basemesh.createBccGrid(Vector3(-1-h,-part-h,-part-h), Vector3(1+h,part+h,part+h), h);
	const Vector3 e0(1,1,0);
	const Vector3 e1(1,-1,0);
	const Vector3 e2(1,0,1);
	const Vector3 e3(1,0,-1);
	const Vector3 e4(0,1,1);
	const Vector3 e5(0,1,-1);
	for(i=basemesh.getNodeSize(); i-->0; ) {
		const Vector3 p = basemesh.getNodePosition3(i);
		if(e0.dot(p) > N || e0.dot(p) < -N) basemesh.removeNode(i);
		else if(e1.dot(p) > N || e1.dot(p) < -N) basemesh.removeNode(i);
		else if(e2.dot(p) > N || e2.dot(p) < -N) basemesh.removeNode(i);
		else if(e3.dot(p) > N || e3.dot(p) < -N) basemesh.removeNode(i);
		else if(e4.dot(p) > n || e4.dot(p) < -n) basemesh.removeNode(i);
		else if(e5.dot(p) > n || e5.dot(p) < -n) basemesh.removeNode(i);
	}
	const uint nodes = basemesh.getNodeSize();
	for(i=0; i<nodes; i++) {
		const Vector3 p = basemesh.getNodePosition3(i);
		basemesh.setNodeFlag(i, uint((p.x + p.y + p.z) / h + 100000.5) % 4);
	}
//*/
	// create space-time mesh
	mesh.setMetric(SymMatrix4(1,0,1,0,0,1,0,0,0,-1));
	for(i=0; i<=steps; i++) {
		for(j=0; j<nodes; j++) {
			const Vector3 x = basemesh.getNodePosition3(j);
			const double t = dtime * (async ? i - 1.0 + double(basemesh.getNodeFlag(j)) / 4.0 : i);
			const double a = exp(-slope * t);
			mesh.addNode(Vector4(a * x, (1.0 - a) / slope));
		}
	}
	for(j=0; j<basemesh.getEdgeSize(); j++) {
		const Buffer<uint> &n = basemesh.getEdgeNodes(j);
		mesh.addEdge(n[0], n[1]);
	}
	for(j=0; j<basemesh.getFaceSize(); j++) mesh.addFace(basemesh.getFaceEdges(j));
	for(j=0; j<basemesh.getBodySize(); j++) mesh.addBody(basemesh.getBodyFaces(j));
	Buffer<uint> nn(nodes);
	for(k=0; k<nodes; k++) nn[k] = k;
	for(i=1; i<=steps; i++) {
		for(j=0; j<4; j++) {
			for(k=0; k<nodes; k++) {
				if(basemesh.getNodeFlag(k) != j) continue;
				nn[k] += nodes;
				mesh.addEdge(nn[k], nn[k] - nodes);
				const Buffer<uint> &e = basemesh.getNodeEdges(k);
				for(l=0; l<e.size(); l++) {
					Buffer<uint> en = basemesh.getEdgeNodes(e[l]);
					for(m=0; m<en.size(); m++) en[m] = nn[en[m]];
					mesh.setEdgeFlag(mesh.addEdge(en[0], en[1]), 1);
					mesh.addFace(findEdges(en[0], en[1], nn[k] - nodes, mesh));
				}
				const Buffer<uint> f = basemesh.getNodeFaces(k);
				for(l=0; l<f.size(); l++) {
					Buffer<uint> fn = basemesh.getFaceNodes(f[l]);
					for(m=0; m<fn.size(); m++) fn[m] = nn[fn[m]];
					mesh.addFace(findEdges(fn[0], fn[1], fn[2], mesh));
					mesh.addBody(findFaces(fn[0], fn[1], fn[2], nn[k] - nodes, mesh));
				}
				const Buffer<uint> b = basemesh.getNodeBodies(k);
				for(l=0; l<b.size(); l++) {
					Buffer<uint> bn = basemesh.getBodyNodes(b[l]);
					for(m=0; m<bn.size(); m++) bn[m] = nn[bn[m]];
					mesh.addBody(findFaces(bn[0], bn[1], bn[2], bn[3], mesh));
					mesh.addQuad(findBodies(bn[0], bn[1], bn[2], bn[3], nn[k] - nodes, mesh));
				}
			}
		}
	}

	return basemesh.getEdgeSize() + (mesh.getEdgeSize() - basemesh.getEdgeSize()) / steps; // return number of edges in the first block
}

void drawWave(const double h, const double dtime, const uint microsteps, const double steptime, const uint steps, const double slope, const uint picwidth = 400) {
	BuilderMesh mesh(4);
	const uint fblock = createSpaceTimeMesh(h, dtime, microsteps + 1, slope, true, mesh);
	const uint nblock = microsteps * (mesh.getNodeSize() / (microsteps + 2));

	// initialize values and operator
	uint i, j, k;
	const Vector3 p0(-0.5, 0.1, 0.0);
	Buffer<double> val(mesh.getEdgeSize(), 0.0);
	Buffer< Buffer< pair<uint, double> > > buf(val.size());
	Buffer<double> hodge(val.size(), 0.0);
	Buffer<Vector4> primv(val.size());
	for(i=0; i<primv.size(); i++) {
		if((i % 100) == 0) cout << i << " / " <<primv.size() << endl;
		primv[i] = mesh.getEdgeVector(i);
	}
	Buffer<ThreeVector4> dualv(val.size(), ThreeVector4(0,0,0,0));
/*// The first (commented) block is for 3-dimensional task. The second is for 4-dimensional task.
	for(i=0; i<mesh.getBodySize(); i++) {
		if((i % 100) == 0) cout << i << " / " << mesh.getBodySize() << endl;
		Buffer<uint> ele;
		Buffer<TwoVector4> v;
		mesh.getBodyEdgeVectors(i, ele, v);
		const Vector4 sign = mesh.getBodyDualVector(i);
		for(j=0; j<ele.size(); j++) dualv[ele[j]] += ThreeVector4(v[j], sign);
	}
/*/	for(i=0; i<mesh.getQuadSize(); i++) {
		if((i % 100) == 0) cout << i << " / " << mesh.getQuadSize() << endl;
		Buffer<uint> ele;
		Buffer<ThreeVector4> v;
		mesh.getQuadEdgeVectors(i, ele, v);
		const double sign = mesh.getQuadDualVector(i);
		for(j=0; j<ele.size(); j++) dualv[ele[j]] += sign * v[j];
	}
//*/	
	const SymMatrix4 met(mesh.getMetric(0));

	for(i=0; i<val.size(); i++) {
		// set initial values
		const Vector3 p = mesh.getEdgeAverage(i).toVector3() - p0;
		const Vector3 v = primv[i].toVector3() / sqrt(12.0);
		const double plen0 = 30.0 * (p + v).len();
		if(plen0 < PI) val[i] = 3.0 * (1.0 + cos(plen0)) * primv[i].t;
		const double plen1 = 30.0 * (p - v).len();
		if(plen1 < PI) val[i] = 3.0 * (1.0 + cos(plen1)) * primv[i].t;
		
		// set equations
		if(i < fblock) { // copy value from the end
			const Buffer<uint> &n = mesh.getEdgeNodes(i);
			const uint ii = mesh.findEdge(n[0] + nblock, n[1] + nblock);
			buf[i].push_back(pair<uint,double>(ii, 1.0));
		}
		else if(mesh.getEdgeFlag(i) == 1) { // update by higher grade cell
			const Buffer<uint> &ele = mesh.getEdgeFaces(i);
			uint cell = ele[0];
			for(j=1; j<ele.size(); j++) {
				if(ele[j] < cell) cell = ele[j];
			}
			const double sign = -mesh.getFaceIncidence(cell, i);
			const Buffer<uint> &f = mesh.getFaceEdges(cell);
			buf[i].resize(f.size() - 1);
			for(j=0, k=0; j<f.size(); j++) {
				if(f[j] == i) continue;
				buf[i][k++] = pair<uint,double>(f[j], sign * mesh.getFaceIncidence(cell, f[j]));
			}
		}
		else { // update by lower grade cell
			const Buffer<uint> &ele = mesh.getEdgeNodes(i);
			uint cell = ele[0];
			for(j=1; j<ele.size(); j++) {
				if(ele[j] < cell) cell = ele[j];
			}
			const double sign = -mesh.getEdgeIncidence(i, cell) * primv[i].lensq() / FourVector4(met * primv[i], dualv[i]).xyzt;
			const Buffer<uint> &f = mesh.getNodeEdges(cell);
			buf[i].resize(f.size() - 1);
			for(j=0, k=0; j<f.size(); j++) {
				if(f[j] == i) continue;
				buf[i][k++] = pair<uint,double>(f[j], sign * mesh.getEdgeIncidence(f[j], cell) * FourVector4(met * primv[f[j]], dualv[f[j]]).xyzt / primv[f[j]].lensq());
			}
		}
	}

	// force curl of initial solution be (almost) zero
	for(k=0; k<10; k++) {
		for(i=0; i<mesh.getFaceSize(); i++) {
			const Buffer<uint> &ele = mesh.getFaceEdges(i);
			double sum = 0.0;
			for(j=0; j<ele.size(); j++) sum += mesh.getFaceIncidence(i, ele[j]) * val[ele[j]];
			sum /= double(ele.size());
			for(j=0; j<ele.size(); j++) val[ele[j]] -= mesh.getFaceIncidence(i, ele[j]) * sum;
		}
	}

	// iterate and draw fields
	double itime = 0.0;
	double idtime = dtime * microsteps;
	const double a = exp(-slope * idtime);
	Picture picture(picwidth,picwidth);
	for(uint iter=0; iter<=steps; ) {
		cout << "Time " << itime << endl;

		// update values
		for(i=0; i<buf.size(); i++) {
			val[i] = 0.0;
			for(j=0; j<buf[i].size(); j++) val[i] += buf[i][j].second * val[buf[i][j].first];
		}

		// draw pictures
		while(itime >= iter * steptime) {
			cout << "Iteration " << iter << "..." << endl;
			const double lensq = 0.6 * mesh.getEdgeVector(0).lensq();
			const double z0 = sqrt(0.125 * lensq);

			uint node = NONE;
			Vector4 color = Vector4(0,0,0, 1);
			for(j=0; j<picture.getHeight(); j++) {
				for(i=0; i<picture.getWidth(); i++) {
					const Vector4 p(2.2 * i / double(picture.getWidth()) - 1.1, 2.2 * j / double(picture.getHeight()) - 1.1, z0, 0);
					const uint newnode = mesh.findNode(p, lensq, node, false);
					if(newnode == NONE) {
						picture.setColor(i, j, Vector4(-1,-1,-1,0));
						continue;
					}
					if(newnode != node) {
						node = newnode;
						const Vector4 tv = getField(mesh, val, node);
						color = Vector4(tv.x, tv.y, tv.z, 1);
					}
					picture.setColor(i, j, color);
				}
			}
			Text pathxy;
			pathxy << "fieldxy" << iter << ".bmp";
			picture.save(pathxy.str(), true);
			
			for(j=0; j<picture.getHeight(); j++) {
				for(i=0; i<picture.getWidth(); i++) {
					const Vector4 p(2.2 * i / double(picture.getWidth()) - 1.1, z0, 2.2 * j / double(picture.getHeight()) - 1.1, 0);
					const uint newnode = mesh.findNode(p, lensq, node, false);
					if(newnode == NONE) {
						picture.setColor(i, j, Vector4(-1,-1,-1,0));
						continue;
					}
					if(newnode != node) {
						node = newnode;
						const Vector4 tv = getField(mesh, val, node);
						color = Vector4(tv.x, tv.y, tv.z, 1);
					}
					picture.setColor(i, j, color);
				}
			}
			Text pathxz;
			pathxz << "fieldxz" << iter << ".bmp";
			picture.save(pathxz.str(), true);

			for(j=0; j<picture.getHeight(); j++) {
				for(i=0; i<picture.getWidth(); i++) {
					const Vector4 p(z0, 2.2 * i / double(picture.getWidth()) - 1.1, 2.2 * j / double(picture.getHeight()) - 1.1, 0);
					const uint newnode = mesh.findNode(p, lensq, node, false);
					if(newnode == NONE) {
						picture.setColor(i, j, Vector4(-1,-1,-1,0));
						continue;
					}
					if(newnode != node) {
						node = newnode;
						const Vector4 tv = getField(mesh, val, node);
						color = Vector4(tv.x, tv.y, tv.z, 1);
					}
					picture.setColor(i, j, color);
				}
			}
			Text pathyz;
			pathyz << "fieldyz" << iter << ".bmp";
			picture.save(pathyz.str(), true);

			const double h = 0.04;
			const uint xsize = 50;
			const uint ysize = 25;
			const uint zsize = 25;
			Buffer<ushort> vol(xsize * ysize * zsize, 32767);
			node = NONE;
			uint nodevol = 32767;
			for(k=0; k<zsize; k++) {
				for(j=0; j<ysize; j++) {
					for(i=0; i<xsize; i++) {
						const Vector4 p(h * (i - 0.5 * (xsize - 1)), h * (j - 0.5 * (ysize - 1)), h * (k - 0.5 * (zsize - 1)), 0);
						const uint newnode = mesh.findNode(p, lensq, node, false);
						if(newnode == NONE) continue;
						if(newnode != node) {
							node = newnode;
							double len = getField(mesh, val, node).t;
							if(len > 1.0) len = 1.0;
							else if(len < -1.0) len = -1.0;
							nodevol = ushort(65535.0 * (0.5 * len + 0.5));
						}
						vol[(k * ysize + j) * xsize + i] = nodevol;
					}
				}
			}
			Text volpath;
			volpath << "volume_" << iter << ".mhd";
			saveVolumeMap(volpath.str(), vol, xsize, ysize, zsize, Vector3(h,h,h));

			iter++;
		}

		// update itime, idtime, and mesh node positions for the next iteration
		itime += idtime;
		idtime *= a;
		for(i=0; i<mesh.getNodeSize(); i++) {
			mesh.setNodePosition(i, a * mesh.getNodePosition(i));
		} 
	}

/*	// draw mesh
	MeshDrawer drawer;
	const Vector3 vo(0,0,5);
	const Vector3 vp(0,20,10);
	const Vector3 vx = TwoVector3(Vector3(0,0,1), vp).dual().unit() / 3.0;
	const Vector3 vy = TwoVector3(vp, vx).dual().unit() / 3.0;
	drawer.initPosition(Vector4(vp + vo,0), Vector4(vo,0), Vector4(vx,0), Vector4(vy,0));
	Picture pic(picwidth, picwidth);
	drawer.initPicture(&pic);
	drawer.initSvg(picwidth, picwidth);
	drawer.drawBoundaryFaces(mesh, Vector3(1,1,1));
	pic.save("mesh.bmp", false);
	drawer.saveSvg("mesh.svg");

	// save mesh statistics
	Text text;
	mesh.writeStatistics(text);
	text.save("stat.txt");
	cout << "Saved statistics" << endl;
*/}

void createVtkMeshes() {
	// the following is a generation routine for .vtk-mesh of cavity boundaries
	for(uint i=0; i<=300; i++) {
		const double h100 = 1.00 - 0.003 * i;
		const double h075 = 0.75 * h100;
		const double h050 = 0.5 * h100;
		const double h025 = 0.25 * h100;
		Text text;
		text << "# vtk DataFile Version 3.0" << endl;
		text << "vtk output" << endl;
		text << "ASCII" << endl;
		text << "DATASET POLYDATA" << endl;
		text << "POINTS 18 float" << endl;
		text << h100 << " 0.0 0.0" << endl;
		text << h050 << " " << h050 << " 0.0" << endl;
		text << h075 << " " << h025 << " " << h025 << endl;
		text << h050 << " 0.0 " << h050 << endl;
		text << h075 << " " << -h025 << " " << h025 << endl;
		text << h050 << " " << -h050 << " 0.0" << endl;
		text << h075 << " " << -h025 << " " << -h025 << endl;
		text << h050 << " 0.0 " << -h050 << endl;
		text << h075 << " " << h025 << " " << -h025 << endl;
		text << -h100 << " 0.0 0.0" << endl;
		text << -h050 << " " << h050 << " 0.0" << endl;
		text << -h075 << " " << h025 << " " << h025 << endl;
		text << -h050 << " 0.0 " << h050 << endl;
		text << -h075 << " " << -h025 << " " << h025 << endl;
		text << -h050 << " " << -h050 << " 0.0" << endl;
		text << -h075 << " " << -h025 << " " << -h025 << endl;
		text << -h050 << " 0.0 " << -h050 << endl;
		text << -h075 << " " << h025 << " " << -h025 << endl;
		text << "POLYGONS 12 68" << endl;
		text << "4 0 2 1 8" << endl;
		text << "4 0 4 3 2" << endl;
		text << "4 0 6 5 4" << endl;
		text << "4 0 8 7 6" << endl;
		text << "4 9 17 10 11" << endl;
		text << "4 9 11 12 13" << endl;
		text << "4 9 13 14 15" << endl;
		text << "4 9 15 16 17" << endl;
		text << "6 1 2 3 12 11 10" << endl;
		text << "6 3 4 5 14 13 12" << endl;
		text << "6 5 6 7 16 15 14" << endl;
		text << "6 7 8 1 10 17 16" << endl;
		Text path;
		path << "obj_" << i << ".vtk";
		text.save(path.str());
	}
}

int main() {
	auto starttime = chrono::system_clock::now();

	BuilderMesh mesh(4);
	drawWave(1.0 / 20.0, 1.0 / 80.0, 1, 0.01, 300, 0.3, 540);

	cout << "Elapsed time: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
	return 0;
}
