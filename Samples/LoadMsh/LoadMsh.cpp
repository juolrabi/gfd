/**
 * Program for loading a mesh file of format .mhs and save to .jrm file format.
 * Call with parameters '.msh_input_file_name' '.jrm_output_file_name'
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace gfd;

double getDouble(const std::string &str, const double def) {
	char *rest;
	const double val = strtod(str.c_str(), &rest);
	if(rest == NULL) return def;
	return val;
}
Buffer<double> getDoubles(const std::string &str) {
	Buffer<double> res;
	const string space = " \t";
	size_t pos = 0;
	pos = str.find_first_not_of(space, pos);
	while(pos != string::npos) {
		const size_t until = str.find_first_of(space, pos);
		res.push_back(getDouble(str.substr(pos, until), 0.0));
		pos = str.find_first_not_of(space, until);
	}
	return res;
}
uint getUint(const std::string &str, const uint def) {
	const double val = getDouble(str, double(def));
	if(val < 0.0) return 0;
	return uint(val + 0.5);
}
Buffer<uint> getUints(const std::string &str) {
	Buffer<uint> res;
	const string space = " \t";
	size_t pos = 0;
	pos = str.find_first_not_of(space, pos);
	while(pos != string::npos) {
		const size_t until = str.find_first_of(space, pos);
		res.push_back(getUint(str.substr(pos, until), 0));
		pos = str.find_first_not_of(space, until);
	}
	return res;
}

int main(int argc, const char* argv[]) {
	if(argc <= 2) {
		std::cout << "Please insert input mesh file path (*.msh) and an output mesh file path (*.jrm)." << std::endl;
		return 0;
	}

	// load text file
	Text input;
	if(input.load(argv[1])) std::cout << ".msh loaded succesfully." << std::endl;
	else {
		std::cout << "Unable to load input mesh file " << argv[1] << std::endl;
		return 0;
	}

	uint i, j;
	Mesh mesh(3);

	// create nodes
	while(input.hasRow()) {
		if(input.getRow().substr(0, 6).compare("$Nodes") == 0) break;
	}
	const Buffer<uint> nodes = getUints(input.getRow());
	if(nodes.size() != 4) { cout << "FAILED reading nodes" << endl; return 0; }
	uint ids = 0;
	Buffer<uint> id(nodes[3]+1);
	mesh.resizeNodeBuffer(nodes[1]);
	for(i=0; i<nodes[0]; i++) {
		const Buffer<uint> buf = getUints(input.getRow());
		if(buf.size() != 4) { cout << "FAILED reading nodes block" << endl; return 0; }
		for(j=0; j<buf[3]; j++) id[getUint(input.getRow(), 0)] = ids++;
		for(j=0; j<buf[3]; j++) {
			const Buffer<double> pos = getDoubles(input.getRow());
			const uint node = mesh.addNode(Vector4(pos[0], (pos.size() > 1 ? pos[1] : 0.0), 
													(pos.size() > 2 ? pos[2] : 0.0), 0.0));
			mesh.setNodeFlag(node, i);
		}
	}
	if(input.getRow().substr(0, 9).compare("$EndNodes") == 0) cout << "Nodes read succesfully." << endl;
	else cout << "Unable to load nodes." << endl;

	// create elements
	while(input.hasRow()) {
		if(input.getRow().substr(0, 9).compare("$Elements") == 0) break;
	}
	const Buffer<uint> elems = getUints(input.getRow());
	if(elems.size() != 4) { cout << "FAILED reading elements" << endl; return 0; }
	for(i=0; i<elems[0]; i++) {
		const Buffer<uint> buf = getUints(input.getRow());
		if(buf.size() != 4) { cout << "FAILED reading elements block" << endl; return 0; }
		for(j=0; j<buf[3]; j++) {
			const Buffer<uint> link = getUints(input.getRow());
			if(link.size() == 3) { // line segment
				const uint edge = mesh.addEdge(id[link[1]], id[link[2]]);
				mesh.setEdgeFlag(edge, i);
			}
			else if(link.size() == 4) { // triangle
				Buffer<uint> edge(3);
				edge[0] = mesh.addEdge(id[link[1]], id[link[2]]);
				edge[1] = mesh.addEdge(id[link[2]], id[link[3]]);
				edge[2] = mesh.addEdge(id[link[3]], id[link[1]]);
				const uint face = mesh.addFace(edge);
				mesh.setFaceFlag(face, i);
			}
			else if(link.size() == 5) { // tetrahedra
				Buffer<uint> edge(6);
				edge[0] = mesh.addEdge(id[link[1]], id[link[2]]);
				edge[1] = mesh.addEdge(id[link[2]], id[link[3]]);
				edge[2] = mesh.addEdge(id[link[3]], id[link[1]]);
				edge[3] = mesh.addEdge(id[link[1]], id[link[4]]);
				edge[4] = mesh.addEdge(id[link[2]], id[link[4]]);
				edge[5] = mesh.addEdge(id[link[3]], id[link[4]]);
				Buffer<uint> face(4);
				Buffer<uint> fe(3);
				fe[0] = edge[0]; fe[1] = edge[1]; fe[2] = edge[2];
				face[0] = mesh.addFace(fe);
				fe[0] = edge[0]; fe[1] = edge[3]; fe[2] = edge[4];
				face[1] = mesh.addFace(fe);
				fe[0] = edge[1]; fe[1] = edge[4]; fe[2] = edge[5];
				face[2] = mesh.addFace(fe);
				fe[0] = edge[2]; fe[1] = edge[5]; fe[2] = edge[3];
				face[3] = mesh.addFace(fe);
				const uint body = mesh.addBody(face);
				mesh.setBodyFlag(body, i);
			}
			else cout << "Unknown element dimension." << endl;
		}
	}
	if(input.getRow().substr(0, 12).compare("$EndElements") == 0) cout << "Elements read succesfully." << endl;
	else cout << "Unable to load elements." << endl;

	mesh.saveJRMesh(argv[2]);

	if(argc > 3) {
		// draw the mesh
		MeshDrawer drawer;
		const Vector3 vo(0,0,0);
		const Vector3 vp(-3,5,10);
		const Vector3 vx = TwoVector3(Vector3(0,1,0), vp-vo).dual().unit() / 1.0;
		const Vector3 vy = TwoVector3(vp-vo, vx).dual().unit() / 1.0;
		drawer.initPosition(Vector4(vp,0), Vector4(vo,0), Vector4(vx,0), Vector4(vy,0));
		drawer.initSvg(500, 500);
		//drawer.drawPrimalEdges(mesh, Vector3(1,0,0));
		drawer.drawBoundaryFaces(mesh, Vector3(1,0.5,0.5));
		drawer.saveSvg(argv[3]);
	}

	if(argc > 4) {
		// save statistics
		Text stat;
		mesh.writeStatistics(stat);
		stat.save(argv[4]);
	}
	return 0;
}

