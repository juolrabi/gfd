/**
 * Program for computing (2+1)-dimensional wave in a rotating cavity.
 * Idea of the solution is following:
 * - create (2+1)-dimensional space-time mesh.
 * - consider discrete forms of grade 2 (values located at faces).
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

TwoVector3 getField(const Mesh &mesh, const Buffer<double> &val, const uint node) {
	const Buffer<uint> par = mesh.getNodeFaces(node);
	SymTwoMatrix3 A(0,0,0,0,0,0);
	TwoVector3 b(0,0,0);
	for(uint i=0; i<par.size(); i++) {
		const TwoVector3 v = mesh.getFaceVector(par[i]).toTwoVector3();
		A += v.outerProduct();
		b += v * val[par[i]];
	}
	return A.inverse() * b;
}

uint createSpaceTimeMesh(const double h, const double dtime, const uint steps, const double slope, const bool async, BuilderMesh &mesh) {
	uint i, j, k, l;
	const double N = 1.0 / h + 1e-5;
	const double n = 0.4 / h + 1e-5;
	const double SQRT3_4 = sqrt(0.75);
	BuilderMesh basemesh(2);
	basemesh.createTriangleGrid(Vector2(-1,-SQRT3_4), Vector2(1,SQRT3_4), h, false);
	const Vector2 e0 = Vector2(0.0, SQRT3_4) / (h * 0.75);
	const Vector2 e1 = Vector2(0.75, 0.5 * SQRT3_4) / (h * 0.75);
	const Vector2 e2 = Vector2(0.75, -0.5 * SQRT3_4) / (h * 0.75);
	for(i=basemesh.getNodeSize(); i-->0; ) {
		const Vector2 p = basemesh.getNodePosition2(i);
		if(e0.dot(p) > N || e0.dot(p) < -N) basemesh.removeNode(i);
		else if(e1.dot(p) > N || e1.dot(p) < -0.3 / h - 0.001) basemesh.removeNode(i);
		else if(e2.dot(p) > N || e2.dot(p) < -N) basemesh.removeNode(i);
		else if(e0.dot(p) < n && e1.dot(p) < n && e2.dot(p) < n) basemesh.removeNode(i);
	}
	const uint nodes = basemesh.getNodeSize();
	for(i=0; i<nodes; i++) {
		const Vector2 x = basemesh.getNodePosition2(i);
		basemesh.setNodeFlag(i, uint((3.0 + x.x + 1.5 * x.y / sqrt(0.75)) / h + 0.5) % 3);
	}

	// create space-time mesh
	const SymMatrix3 met(1,0,1,0,0,-1);
	mesh.setMetric(SymMatrix4(met,0,0,0,1));
	for(i=0; i<=steps; i++) {
		for(j=0; j<nodes; j++) {
			const Vector2 x = basemesh.getNodePosition2(j);
			const double z = dtime * (async ? i - 1.0 + double(basemesh.getNodeFlag(j)) / 3.0 : i);
			const Vector2 b(-0.0 * z,0);
			const double cosi = cos(slope * z);
			const double sini = sin(slope * z);
			const Matrix2 a(cosi,sini, -sini,cosi);
			mesh.addNode(Vector4(a * x + b, z, 0.0));
		}
	}
	for(j=0; j<basemesh.getEdgeSize(); j++) {
		const Buffer<uint> &n = basemesh.getEdgeNodes(j);
		mesh.addEdge(n[0], n[1]);
	}
	for(j=0; j<basemesh.getFaceSize(); j++) mesh.addFace(basemesh.getFaceEdges(j));
	for(i=1; i<=steps; i++) {
		Buffer<uint> nn(nodes);
		for(j=0; j<3; j++) {
			Buffer<uint> nn(nodes);
			for(k=0; k<nodes; k++) {
				nn[k] = (basemesh.getNodeFlag(k) > j ? i-1 : i) * nodes + k;
				if(basemesh.getNodeFlag(k) == j) mesh.addEdge(nn[k], nn[k] - nodes);
			}
			for(k=0; k<basemesh.getEdgeSize(); k++) {
				const Buffer<uint> &n = basemesh.getEdgeNodes(k);
				for(l=0; l<2; l++) {
					if(basemesh.getNodeFlag(n[l]) == j) {
						const uint n0 = nn[n[l]];
						const uint n1 = nn[n[(l+1)%2]];
						Buffer<uint> e(3);
						e[0] = mesh.findEdge(n0 - nodes, n1);
						e[1] = mesh.findEdge(n0, n0 - nodes);
						e[2] = mesh.addEdge(n1, n0);
						mesh.addFace(e);
						break;
					}
				}
			}
			for(k=0; k<basemesh.getFaceSize(); k++) {
				const Buffer<uint> n = basemesh.getFaceNodes(k); // size must be 3
				Buffer<uint> fn(4);
				for(l=0; l<3; l++) {
					fn[l] = nn[n[l]];
					if(basemesh.getNodeFlag(n[l]) == j) fn[3] = fn[l] - nodes;
				}
				Buffer<uint> fe(6);
				fe[0] = mesh.findEdge(fn[0], fn[1]);
				fe[1] = mesh.findEdge(fn[1], fn[2]);
				fe[2] = mesh.findEdge(fn[2], fn[0]);
				fe[3] = mesh.findEdge(fn[3], fn[0]);
				fe[4] = mesh.findEdge(fn[3], fn[1]);
				fe[5] = mesh.findEdge(fn[3], fn[2]);
				Buffer<uint> e(3);
				Buffer<uint> f(4);
				e[0] = fe[0]; e[1] = fe[3]; e[2] = fe[4];
				f[0] = mesh.findFace(e);
				e[0] = fe[1]; e[1] = fe[4]; e[2] = fe[5];
				f[1] = mesh.findFace(e);
				e[0] = fe[2]; e[1] = fe[3]; e[2] = fe[5];
				f[2] = mesh.findFace(e);
				e[0] = fe[0]; e[1] = fe[1]; e[2] = fe[2];
				f[3] = mesh.addFace(e);
				mesh.setFaceFlag(f[3], 1);
				mesh.addBody(f);
			}
		}
	}
	return 3 * basemesh.getFaceSize() + 2 * basemesh.getEdgeSize(); // return number of faces in a block
}

void drawMeshBar(const double h, const double dtime, const uint steps, const double slope, const uint picwidth = 400) {
	BuilderMesh mesh(3);
	createSpaceTimeMesh(h, dtime, steps, slope, false, mesh);
	for(uint i=0; i<=steps; i++)
	{
		cout << "Drawing mesh " << i << "..." << endl;
		const double height = dtime * i;
		Buffer<Vector3> col(mesh.getFaceSize());
		for(uint j=0; j<col.size(); j++) {
			if(mesh.getFaceAverage(j).z < height) col[j] = Vector3(1,0,0);
			else col[j] = Vector3(1,1,1);
		}

		MeshDrawer drawer;
		const Vector3 vo(0,0,3.7);
		const Vector3 vp(0,-100,100);
		const Vector3 vx = TwoVector3(Vector3(0,0,1), vp).dual().unit() / 6.7;
		const Vector3 vy = TwoVector3(vp, vx).dual().unit() / 6.7;
		drawer.initPosition(Vector4(vp + vo,0), Vector4(vo,0), Vector4(vx,0), Vector4(vy,0));
		Picture pic(picwidth, picwidth);
		drawer.initPicture(&pic);
		drawer.drawBoundaryFaces(mesh, col);
		Text path;
		path << "mesh" << i << ".bmp";
		pic.save(path.str(), false);
	}
}

void drawWave(const double h, const double dtime, const uint microsteps, const uint steps, const double slope, const uint picwidth = 400) {
	BuilderMesh mesh(3);
	const uint block = createSpaceTimeMesh(h, dtime, microsteps, slope, true, mesh);
	const uint block0 = mesh.getFaceSize() - (microsteps - 1) * block;

	// initialize values and operator
	uint i, j, k;
	const Vector2 p0(-0.7 * sin(PI/6.0), 0.7 * cos(PI/6.0));
	Buffer<double> val(mesh.getFaceSize(), 0.0);
	Buffer< Buffer< pair<uint, double> > > buf(val.size());
	for(i=0; i<val.size(); i++) {
		if(i < block0) { // copy value from the end
			const TwoVector3 tv = mesh.getFaceVector3(i);
			const Vector2 p = mesh.getFaceAverage3(i).toVector2() - p0;
			const double plen = 50.0 * p.len();
			if(plen < PI) val[i + (microsteps-1) * block] = 3.0 * (1.0 + cos(plen)) * tv.xy;
			buf[i].push_back(pair<uint,double>(i + (microsteps-1) * block, 1.0));
		}
		else if(mesh.getFaceFlag(i) == 1) { // update by a body
			const uint body = mesh.getFaceBodies(i)[0];
			const double sign = -mesh.getBodyIncidence(body, i);
			const Buffer<uint> &f = mesh.getBodyFaces(body);
			buf[i].resize(f.size() - 1);
			for(j=0, k=0; j<f.size(); j++) {
				if(f[j] == i) continue;
				buf[i][k++] = pair<uint,double>(f[j], sign * mesh.getBodyIncidence(body, f[j]));
			}
		}
		else { // update by an edge
			const uint edge = mesh.getFaceEdges(i)[0];
			const double sign = -mesh.getFaceIncidence(i, edge) / mesh.getFaceHodge(i);
			const Buffer<uint> &f = mesh.getEdgeFaces(edge);
			buf[i].resize(f.size() - 1);
			for(j=0, k=0; j<f.size(); j++) {
				if(f[j] == i) continue;
				buf[i][k++] = pair<uint,double>(f[j], sign * mesh.getFaceIncidence(f[j], edge) * mesh.getFaceHodge(f[j]));
			}
		}
	}

	// iterate and draw fields
	const uint node0 = mesh.findNode(Vector4(p0,0,0), 0.5 * h * h, 0, false);
	const uint node1 = mesh.findNode(Vector4(0.7 * cos(PI/6.0), -0.7 * sin(PI/6.0),0,0), 0.5 * h * h, 0, false);
	const double cosi = cos(slope * dtime * microsteps);
	const double sini = sin(slope * dtime * microsteps);
	const Matrix4 A(cosi,sini,0,0, -sini,cosi,0,0, 0,0,1,0, 0,0,0,0);
	Picture picture(picwidth,picwidth);
	picture.fillColor(Vector4(-1,-1,-1,0));
	for(uint iter=0; iter<=steps; iter++) {
		cout << "Iteration " << iter << "..." << endl;

		// update values
		for(i=0; i<buf.size(); i++) {
			val[i] = 0.0;
			for(j=0; j<buf[i].size(); j++) val[i] += buf[i][j].second * val[buf[i][j].first];
		}

		// draw picture
		uint node = node0;
		const TwoVector3 tv = getField(mesh, val, node);
		Vector4 color(tv.xy, tv.xz, tv.yz, 1);
		for(j=0; j<picture.getHeight(); j++) {
			for(i=0; i<picture.getWidth(); i++) {
				if(iter > 0 && picture.getColor(i,j).t == 0.0 && i >= 5 && picture.getColor(i-5,j).t == 0.0 &&
					i + 5 < picture.getWidth() && picture.getColor(i+5,j).t == 0.0 &&
					j >= 5 && picture.getColor(i,j-5).t == 0.0 &&
					j + 5 < picture.getHeight() && picture.getColor(i,j+5).t == 0.0) continue;

				const Vector4 p(2.1 * i / double(picture.getWidth()) - 1.05, 2.1 * j / double(picture.getHeight()) - 1.05, 0,0);
				uint newnode = mesh.findNode(p, 0.5 * h * h, node, false);
				if(newnode == NONE) {
					const double sq0 = (mesh.getNodePosition(node0) - p).lensq();
					const double sq1 = (mesh.getNodePosition(node1) - p).lensq();
					if(sq0 < sq1) {
						if(sq0 < (mesh.getNodePosition(node) - p).lensq()) newnode = mesh.findNode(p, 0.5 * h * h, node0, false);
					}
					else if(sq1 < (mesh.getNodePosition(node) - p).lensq()) newnode = mesh.findNode(p, 0.5 * h * h, node1, false);
					if(newnode == NONE) {
						picture.setColor(i, j, Vector4(-1,-1,-1,0));
						continue;
					}
				}
				if(newnode != node) {
					node = newnode;
					const TwoVector3 tv = getField(mesh, val, node);
					color = Vector4(tv.xy, tv.xz, tv.yz, 1);
				}
				picture.setColor(i, j, color);
			}
		}
		Text path;
		path << "field" << iter << ".bmp";
		picture.save(path.str(), true);

		// update mesh node positions for the next iteration
		for(i=0; i<mesh.getNodeSize(); i++) {
			mesh.setNodePosition(i, A * mesh.getNodePosition(i));
		} 
	}

	// draw and save statistics
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
}

int main() {
	auto starttime = chrono::system_clock::now();

	drawMeshBar(1.0 / 4.0, 1.0 / 10.0, 80, PI / 12.0, 400);
	drawWave(1.0 / 40.0, 1.0 / 100.0, 10, 80, PI / 12.0, 400);

	cout << "Elapsed time: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
	return 0;
}
