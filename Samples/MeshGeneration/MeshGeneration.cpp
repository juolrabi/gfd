/**
 * Program for generating a Delaunay-mesh with variable metric and weights.
 * This sample is under construction.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include "../../GFD/Types/Random.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace gfd;

int main() {
	uint i, j;
	Random rnd(134); // initialize random variables with seed
	BuilderMesh mesh(4);

	//mesh.createTriangleGrid(Vector2(-1,-1), Vector2(1,1), 0.5, true);
	//mesh.setMetric(SymMatrix4(-1,0,1,0,0,1,0,0,0,1), 1);
	for(i=0; i<20; i++) {
		mesh.insertNode(Vector4(rnd.getUniformSphere2(),0,0), 0.0, 0, false);
	}
	mesh.fillBoundaryFlags(1);
	//mesh.optimizeNodes(0, 100, true, false);

	//mesh.createGrid(Vector4(-1,-1,-1,-1), Vector4(1,1,1,1), 1);
	for(i=0; i<mesh.getEdgeSize(); i++) {
		Buffer<Vector4> p;
		mesh.getEdgeSimplices(i, p);
		Vector4 sum(0,0,0,0);
		for(j=0; j<p.size(); j+=2) {
			sum += p[j+1] - p[j];
		}
		const double diff = (mesh.getEdgeVector(i) - sum).len();
		if(fabs(diff) > 1e-12) cout << "Edge " << diff << " " << sum.len() << endl;
	}
	for(i=0; i<mesh.getFaceSize(); i++) {
		Buffer<Vector4> p;
		mesh.getFaceSimplices(i, p);
		TwoVector4 sum(0,0,0,0,0,0);
		for(j=0; j<p.size(); j+=3) {
			sum += 0.5 * TwoVector4(p[j+1] - p[j], p[j+2] - p[j]);
		}
		const double diff = (mesh.getFaceVector(i) - sum).len();
		if(fabs(diff) > 1e-12) cout << "Face " << diff << " " << sum.len() << endl;
	}
	for(i=0; i<mesh.getBodySize(); i++) {
		Buffer<Vector4> p;
		mesh.getBodySimplices(i, p);
		ThreeVector4 sum(0,0,0,0);
		for(j=0; j<p.size(); j+=4) {
			sum += ThreeVector4(p[j+1] - p[j], p[j+2] - p[j], p[j+3] - p[j]) / 6.0;
		}
		const double diff = (mesh.getBodyVector(i) - sum).len();
		if(fabs(diff) > 1e-12) cout << "Body " << diff << " " << sum.len() << endl;
	}
	for(i=0; i<mesh.getQuadSize(); i++) {
		Buffer<Vector4> p;
		mesh.getQuadSimplices(i, p);
		FourVector4 sum(0);
		for(j=0; j<p.size(); j+=5) {
			sum += FourVector4(p[j+1] - p[j], p[j+2] - p[j], p[j+3] - p[j], p[j+4] - p[j]) / 24.0;
		}
		const double diff = (mesh.getQuadVector(i) - sum).len();
		if(fabs(diff) > 1e-12) cout << "Quad " << diff << " " << sum.len() << endl;
	}

	Buffer<FourVector4> nd(mesh.getNodeSize(), FourVector4(0));
	if(mesh.getFaceSize() == 0) { 
		for(i=0; i<nd.size(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getNodeEdgeSimplices(i, s, p);
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[2 * j];
				const ThreeVector4 dv = mesh.getEdgeVector(s[j]).unit().dual();
				nd[i] += FourVector4(pj[1] - pj[0], dv);
			}
		}
/*		for(i=0; i<mesh.getEdgeSize(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getEdgeNodeSimplices(i, s, p);
			const ThreeVector4 dv = mesh.getEdgeVector(i).unit().dual();
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[2 * j];
				nd[s[j]] += FourVector4(pj[1] - pj[0], dv);
			}
		}
*/	}
	else if(mesh.getBodySize() == 0) {
		for(i=0; i<nd.size(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getNodeFaceSimplices(i, s, p);
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[3 * j];
				const TwoVector4 dv = mesh.getFaceVector(s[j]).unit().dual() / 2.0;
				nd[i] += FourVector4(pj[1] - pj[0], pj[2] - pj[0], dv);
			}
		}
/*		for(i=0; i<mesh.getFaceSize(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getFaceNodeSimplices(i, s, p);
			const TwoVector4 dv = mesh.getFaceVector(i).unit().dual() / 2.0;
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[3 * j];
				nd[s[j]] += FourVector4(pj[1] - pj[0], pj[2] - pj[0], dv);
			}
		}
*/	}
	else if(mesh.getQuadSize() == 0) {
		for(i=0; i<nd.size(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getNodeBodySimplices(i, s, p);
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[4 * j];
				const Vector4 dv = mesh.getBodyVector(s[j]).unit().dual() / 6.0;
				nd[i] += FourVector4(pj[1] - pj[0], pj[2] - pj[0], pj[3] - pj[0], dv);
			}
		}
/*		for(i=0; i<mesh.getBodySize(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getBodyNodeSimplices(i, s, p);
			const Vector4 dv = mesh.getBodyVector(i).unit().dual() / 6.0;
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[4 * j];
				nd[s[j]] += FourVector4(pj[1] - pj[0], pj[2] - pj[0], pj[3] - pj[0], dv);
			}
		}
*/	}
	else {
		for(i=0; i<nd.size(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getNodeQuadSimplices(i, s, p);
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[5 * j];
				const double dv = mesh.getQuadVector(s[j]).unit().dual() / 24.0;
				nd[i] += FourVector4(pj[1] - pj[0], pj[2] - pj[0], pj[3] - pj[0], pj[4] - pj[0]) * dv;
			}
		}
/*		for(i=0; i<mesh.getQuadSize(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getQuadNodeSimplices(i, s, p);
			const double dv = mesh.getQuadVector(i).unit().dual() / 24.0;
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[5 * j];
				nd[s[j]] += FourVector4(pj[1] - pj[0], pj[2] - pj[0], pj[3] - pj[0], pj[4] - pj[0]) * dv;
			}
		}
*/	}
	for(i=0; i<nd.size(); i++) {
		const double diff = (mesh.getNodeDualVector(i) - nd[i]).len();
		if(fabs(diff) > 1e-12) cout << "Node dual " << diff << " " << nd[i].len() << endl;
		else cout << "Node jee" << endl;
	}

	Buffer<ThreeVector4> ed(mesh.getEdgeSize(), ThreeVector4(0,0,0,0));
	if(mesh.getFaceSize() == 0) { 
		for(i=0; i<mesh.getEdgeSize(); i++) {
			ed[i] = mesh.getEdgeVector(i).unit().dual();
		}
	}
	else if(mesh.getBodySize() == 0) {
		for(i=0; i<ed.size(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getEdgeFaceSimplices(i, s, p);
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[2 * j];
				const TwoVector4 dv = mesh.getFaceVector(s[j]).unit().dual();
				ed[i] += ThreeVector4(pj[1] - pj[0], dv);
			}
		}
/*		for(i=0; i<mesh.getFaceSize(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getFaceEdgeSimplices(i, s, p);
			const TwoVector4 dv = mesh.getFaceVector(i).unit().dual();
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[2 * j];
				ed[s[j]] += ThreeVector4(pj[1] - pj[0], dv);
			}
		}
*/	}
	else if(mesh.getQuadSize() == 0) {
		for(i=0; i<ed.size(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getEdgeBodySimplices(i, s, p);
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[3 * j];
				const Vector4 dv = mesh.getBodyVector(s[j]).unit().dual() / 2.0;
				ed[i] += ThreeVector4(pj[1] - pj[0], pj[2] - pj[0], dv);
			}
		}
/*		for(i=0; i<mesh.getBodySize(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getBodyEdgeSimplices(i, s, p);
			const Vector4 dv = mesh.getBodyVector(i).unit().dual() / 2.0;
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[3 * j];
				ed[s[j]] += ThreeVector4(pj[1] - pj[0], pj[2] - pj[0], dv);
			}
		}
*/	}
	else {
		for(i=0; i<ed.size(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getEdgeQuadSimplices(i, s, p);
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[4 * j];
				const double dv = mesh.getQuadVector(s[j]).unit().dual() / 6.0;
				ed[i] += ThreeVector4(pj[1] - pj[0], pj[2] - pj[0], pj[3] - pj[0]) * dv;
			}
		}
/*		for(i=0; i<mesh.getQuadSize(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getQuadEdgeSimplices(i, s, p);
			const double dv = mesh.getQuadVector(i).unit().dual() / 6.0;
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[4 * j];
				ed[s[j]] += ThreeVector4(pj[1] - pj[0], pj[2] - pj[0], pj[3] - pj[0]) * dv;
			}
		}
*/	}
	for(i=0; i<ed.size(); i++) {
		const double diff = (mesh.getEdgeDualVector(i) - ed[i]).len();
		if(fabs(diff) > 1e-12) cout << "Edge dual " << diff << " " << ed[i].len() << endl;
		else cout << "Edge jee" << endl;
	}

	Buffer<TwoVector4> fd(mesh.getFaceSize(), TwoVector4(0,0,0,0,0,0));
	if(mesh.getBodySize() == 0) {
		for(i=0; i<mesh.getFaceSize(); i++) {
			fd[i] = mesh.getFaceVector(i).unit().dual();
		}
	}
	else if(mesh.getQuadSize() == 0) {
		for(i=0; i<fd.size(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getFaceBodySimplices(i, s, p);
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[2 * j];
				const Vector4 dv = mesh.getBodyVector(s[j]).unit().dual();
				fd[i] += TwoVector4(pj[1] - pj[0], dv);
			}
		}

/*		for(i=0; i<mesh.getBodySize(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getBodyFaceSimplices(i, s, p);
			const Vector4 dv = mesh.getBodyVector(i).unit().dual();
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[2 * j];
				fd[s[j]] += TwoVector4(pj[1] - pj[0], dv);
			}
		}
*/	}
	else {
		for(i=0; i<fd.size(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getFaceQuadSimplices(i, s, p);
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[3 * j];
				const double dv = mesh.getQuadVector(s[j]).unit().dual() / 2.0;
				fd[i] += TwoVector4(pj[1] - pj[0], pj[2] - pj[0]) * dv;
			}
		}
/*		for(i=0; i<mesh.getQuadSize(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getQuadFaceSimplices(i, s, p);
			const double dv = mesh.getQuadVector(i).unit().dual() / 2.0;
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[3 * j];
				fd[s[j]] += TwoVector4(pj[1] - pj[0], pj[2] - pj[0]) * dv;
			}
		}
*/	}
	for(i=0; i<fd.size(); i++) {
		const double diff = (mesh.getFaceDualVector(i) - fd[i]).len();
		if(fabs(diff) > 1e-12) cout << "Face dual " << diff << " " << fd[i].len() << endl;
		else cout << "Face jee" << endl;
	}

	Buffer<Vector4> bd(mesh.getBodySize(), Vector4(0,0,0,0));
	if(mesh.getQuadSize() == 0) {
		for(i=0; i<mesh.getBodySize(); i++) {
			bd[i] = mesh.getBodyVector(i).unit().dual();
		}
	}
	else {
		for(i=0; i<bd.size(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getBodyQuadSimplices(i, s, p);
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[2 * j];
				const double dv = mesh.getQuadVector(s[j]).unit().dual();
				bd[i] += (pj[1] - pj[0]) * dv;
			}
		}
/*		for(i=0; i<mesh.getQuadSize(); i++) {
			Buffer<uint> s;
			Buffer<Vector4> p;
			mesh.getQuadBodySimplices(i, s, p);
			const double dv = mesh.getQuadVector(i).unit().dual();
			for(j=0; j<s.size(); j++) {
				const Vector4 *pj = &p[2 * j];
				bd[s[j]] += (pj[1] - pj[0]) * dv;
			}
		}
*/	}
	for(i=0; i<bd.size(); i++) {
		const double diff = (mesh.getBodyDualVector(i) - bd[i]).len();
		if(fabs(diff) > 1e-12) cout << "Body dual " << diff << " " << bd[i].len() << endl;
		else cout << "Body jee" << endl;
	}


/*	mesh.insertNode(Vector4(0,0,0,0), 0.0, 1, 0, true);
	mesh.insertNode(Vector4(1,0,0,0), 0.0, 1, 0, true);
	mesh.insertNode(Vector4(0,1,0,0), 0.0, 1, 0, true);
	//mesh.insertNode(Vector4(1,1,0,0), 0.0, 0, 0, true);
	//mesh.setNodeMetric(1, 1);
*/
/*	for(i=0; i<mesh.getEdgeSize(); i++) {
		const Vector4 p = mesh.getEdgePosition(i);
		const Buffer<uint> &n = mesh.getEdgeNodes(i);
		for(j=0; j<n.size(); j++) cout << mesh.getRadiusSq(p, n[j]) << " ";
		cout << endl;
	}
	for(i=0; i<mesh.getFaceSize(); i++) {
		const Vector4 p = mesh.getFacePosition(i);
		const Buffer<uint> n = mesh.getFaceNodes(i);
		for(j=0; j<n.size(); j++) cout << mesh.getRadiusSq(p, n[j]) << " ";
		cout << endl;
	}
*/
	// draw the mesh
	MeshDrawer drawer;
	const Vector3 vo(0,0,0);
	const Vector3 vp(0,0,10);
	const Vector3 vx = TwoVector3(Vector3(0,1,0), vp-vo).dual().unit() / 3.0;
	const Vector3 vy = TwoVector3(vp-vo, vx).dual().unit() / 3.0;
	drawer.initPosition(Vector4(vp,0), Vector4(vo,0), Vector4(vx,0), Vector4(vy,0));
	drawer.initSvg(500, 500);
	drawer.drawPrimalEdges(mesh, Vector3(0,0,0));
	drawer.drawDualEdges(mesh, Vector3(1,0,0));
	drawer.drawPrimalNodes(mesh, Vector3(0,0,1));
	//drawer.drawPrimalEdges(mesh, Vector3(1,0,0));
	//drawer.drawBoundaryFaces(mesh, Vector3(1,0.5,0.5));
	drawer.saveSvg("mesh.svg");

	// save statistics
	Text stat;
	mesh.writeStatistics(stat, 0);
	stat.save("stat.txt");

	cout << "Ready" << endl;
	return 0;
}

