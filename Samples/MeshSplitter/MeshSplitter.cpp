/**
 * Program for splitting a mesh by domain decomposition method.
 * Call with parameters 'mesh_file_name' 'split_command_1' 'split_command_2' ...
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#include "../../GFD/Mesh/PartMesh.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Types/Types.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include <fstream>
#include <iostream>
#include <map>

using namespace std;
using namespace gfd;

uint minimum(const Buffer<uint> &part, const Buffer<uint> &ele) {
	uint res = NONE;
	for(uint i=0; i<ele.size(); i++) {
		if(res > part[ele[i]]) res = part[ele[i]];
	}
	return res;
}

void splitElements(const uint newparts, const Mesh &mesh, const Buffer<uint> &npart, Buffer<uint> &epart, Buffer<uint> &fpart, Buffer<uint> &bpart, Buffer<uint> &qpart, uint &parts) {
	uint i;
	for(i=0; i<epart.size(); i++) {
		if(epart[i] < parts)  epart[i] = minimum(npart, mesh.getEdgeNodes(i));
		else if(epart[i] != NONE) epart[i] += newparts - parts;
	}
	for(i=0; i<fpart.size(); i++) {
		if(fpart[i] < parts)  fpart[i] = minimum(epart, mesh.getFaceEdges(i));
		else if(fpart[i] != NONE) fpart[i] += newparts - parts;
	}
	for(i=0; i<bpart.size(); i++) {
		if(bpart[i] < parts)  bpart[i] = minimum(fpart, mesh.getBodyFaces(i));
		else if(bpart[i] != NONE) bpart[i] += newparts - parts;
	}
	for(i=0; i<qpart.size(); i++) {
		if(qpart[i] < parts)  qpart[i] = minimum(bpart, mesh.getQuadBodies(i));
		else if(qpart[i] != NONE) qpart[i] += newparts - parts;
	}
	parts = newparts;
}

bool split(const uint divide, double (*func)(const Vector4 &), const Mesh &mesh, Buffer<uint> &npart, Buffer<uint> &epart, Buffer<uint> &fpart, Buffer<uint> &bpart, Buffer<uint> &qpart, uint &parts) {
	if(divide == 0) return false;
	if(divide == 1) return true;

	// order nodes by func
	uint i, j;
	const uint newparts = parts * divide;
	Buffer< multimap<double,uint> > tree(parts);
	for(i=0; i<npart.size(); i++) {
		if(npart[i] < parts) tree[npart[i]].insert(pair<double,uint>(func(mesh.getNodePosition(i)), i));
		else if(npart[i] != NONE) npart[i] += newparts - parts;
	}
	for(i=0; i<parts; i++) {
		const uint treesize = uint(tree[i].size());
		auto it = tree[i].begin();
		for(j=0; j<treesize; j++,it++) {
			npart[it->second] = divide * npart[it->second] + uint(j * divide / treesize);
		}
	}
	splitElements(newparts, mesh, npart, epart, fpart, bpart, qpart, parts);
	return true;
}

bool flagSplit(const Mesh &mesh, Buffer<uint> &npart, Buffer<uint> &epart, Buffer<uint> &fpart, Buffer<uint> &bpart, Buffer<uint> &qpart, uint &parts) {
	// order nodes by func
	uint i;
	Buffer<uint> flags(parts, 0);
	Buffer< Buffer<uint> > flag(parts);
	for(i=0; i<npart.size(); i++) {
		if(npart[i] < parts) flag[npart[i]].gatherOnce(mesh.getNodeFlag(i), flags[npart[i]]);
		if(npart[i] < parts) cout << "flag " << mesh.getNodeFlag(i) << endl;
	}
	Buffer<uint> part0(parts, 0);
	for(i=1; i<parts; i++) part0[i] = part0[i-1] + flags[i-1];
	const uint newparts = part0.back() + flags.back();

	for(i=0; i<npart.size(); i++) {
		if(npart[i] < parts) npart[i] = part0[npart[i]] + flag[npart[i]].findFirst(mesh.getNodeFlag(i));
		else if(npart[i] != NONE) npart[i] += newparts - parts;
	}
	splitElements(newparts, mesh, npart, epart, fpart, bpart, qpart, parts);
	return true;
}

void chain(Buffer<uint> &part, const uint parts) {
	for(uint i=0; i<part.size(); i++) {
		if(part[i] == NONE || part[i] < parts) continue;
		const uint parti = part[part[i] - parts];
		if(parti == NONE || parti < parts) continue;
		part[i] = parti;
	}
}

bool repeat(const double v0, const double v1, double (*func)(const Vector4 &), Vector4 (*_func)(const Vector4 &, const double), const Mesh &mesh, Buffer<uint> &npart, Buffer<uint> &epart, Buffer<uint> &fpart, Buffer<uint> &bpart, Buffer<uint> &qpart, const uint parts) {
	if(v0 >= v1) {
		cout << "Repeat failed: Second parameter should be larger than the first." << endl;
		return false;
	}
	if(parts != 1) cout << "Please apply repeat before split-commands." << endl;
	
	uint i, j;
	const double vM = 2.0 * v0 - v1 - 1e-8;
	Buffer<uchar> nslot(npart.size(), 0);
	for(i=0; i<npart.size(); i++) {
		if(npart[i] == NONE) continue;
		const Vector4 p = mesh.getNodePosition(i);
		const double v = func(p);
		if(v > v1 || v < vM) npart[i] = NONE;
		else if(v <= v0) npart[i] = parts + i;
	}
	for(i=0; i<epart.size(); i++) {
		if(epart[i] == NONE) continue;
		uint sum = 0;
		const Buffer<uint> &ele = mesh.getEdgeNodes(i);
		for(j=0; j<ele.size(); j++) {
			const uint jpart = npart[ele[j]];
			if(jpart == NONE) {
				epart[i] = NONE; 
				break;
			}
			else if(jpart == parts + ele[j]) sum++;
		}
		if(sum == j) epart[i] = parts + i;
	}
	for(i=0; i<fpart.size(); i++) {
		if(fpart[i] == NONE) continue;
		uint sum = 0;
		const Buffer<uint> &ele = mesh.getFaceEdges(i);
		for(j=0; j<ele.size(); j++) {
			const uint jpart = epart[ele[j]];
			if(jpart == NONE) {
				fpart[i] = NONE; 
				break;
			}
			else if(jpart == parts + ele[j]) sum++;
		}
		if(sum == j) fpart[i] = parts + i;
	}
	for(i=0; i<bpart.size(); i++) {
		if(bpart[i] == NONE) continue;
		uint sum = 0;
		const Buffer<uint> &ele = mesh.getBodyFaces(i);
		for(j=0; j<ele.size(); j++) {
			const uint jpart = fpart[ele[j]];
			if(jpart == NONE) {
				bpart[i] = NONE; 
				break;
			}
			else if(jpart == parts + ele[j]) sum++;
		}
		if(sum == j) bpart[i] = parts + i;
	}
	for(i=0; i<qpart.size(); i++) {
		if(qpart[i] == NONE) continue;
		uint sum = 0;
		const Buffer<uint> &ele = mesh.getQuadBodies(i);
		for(j=0; j<ele.size(); j++) {
			const uint jpart = bpart[ele[j]];
			if(jpart == NONE) {
				qpart[i] = NONE; 
				break;
			}
			else if(jpart == parts + ele[j]) sum++;
		}
		if(sum == j) qpart[i] = parts + i;
	}

	// ignore unnecessary elements (set part = NONE, where possible)
	for(i=0; i<qpart.size(); i++) {
		if(qpart[i] != parts + i) continue;
		qpart[i] = NONE;
	}
	for(i=0; i<bpart.size(); i++) {
		if(bpart[i] != parts + i) continue;
		const Buffer<uint> &ele = mesh.getBodyQuads(i);
		for(j=0; j<ele.size(); j++) {
			if(qpart[ele[j]] != NONE) break;
		}
		if(j == ele.size()) bpart[i] = NONE;
	}
	for(i=0; i<fpart.size(); i++) {
		if(fpart[i] != parts + i) continue;
		const Buffer<uint> &ele = mesh.getFaceBodies(i);
		for(j=0; j<ele.size(); j++) {
			if(bpart[ele[j]] != NONE) break;
		}
		if(j == ele.size()) fpart[i] = NONE;
	}
	for(i=0; i<epart.size(); i++) {
		if(epart[i] != parts + i) continue;
		const Buffer<uint> &ele = mesh.getEdgeFaces(i);
		for(j=0; j<ele.size(); j++) {
			if(fpart[ele[j]] != NONE) break;
		}
		if(j == ele.size()) epart[i] = NONE;
	}
	for(i=0; i<npart.size(); i++) {
		if(npart[i] != parts + i) continue;
		const Buffer<uint> &ele = mesh.getNodeEdges(i);
		for(j=0; j<ele.size(); j++) {
			if(epart[ele[j]] != NONE) break;
		}
		if(j == ele.size()) npart[i] = NONE;
	}

	// find repeating matched elements
	for(i=0; i<npart.size(); i++) {
		if(npart[i] != parts + i) continue;
		const Vector4 p = _func(mesh.getNodePosition(i), v1 - v0);
		const uint match = mesh.findNode(p, 1e-8, i, true);
		if(match == NONE) cout << "Repeat failed: Can not find matched node." << endl;
		else npart[i] = parts + match;
	}
	for(i=0; i<epart.size(); i++) {
		if(epart[i] != parts + i) continue;
		Buffer<uint> ele = mesh.getEdgeNodes(i);
		for(j=0; j<ele.size(); j++) ele[j] = npart[ele[j]]-parts;
		const uint match = mesh.findEdge(ele[0], ele[1]);
		if(match == NONE) cout << "Repeat failed: Can not find matched edge." << endl;
		else epart[i] = parts + match;
	}
	for(i=0; i<fpart.size(); i++) {
		if(fpart[i] != parts + i) continue;
		Buffer<uint> ele = mesh.getFaceEdges(i);
		for(j=0; j<ele.size(); j++) ele[j] = epart[ele[j]]-parts;
		const uint match = mesh.findFace(ele);
		if(match == NONE) cout << "Repeat failed: Can not find matched face." << endl;
		else fpart[i] = parts + match;
	}
	for(i=0; i<bpart.size(); i++) {
		if(bpart[i] != parts + i) continue;
		Buffer<uint> ele = mesh.getBodyFaces(i);
		for(j=0; j<ele.size(); j++) ele[j] = fpart[ele[j]]-parts;
		const uint match = mesh.findBody(ele);
		if(match == NONE) cout << "Repeat failed: Can not find matched body." << endl;
		else bpart[i] = parts + match;
	}

	// chain links
	chain(npart, parts);
	chain(epart, parts);
	chain(fpart, parts);
	chain(bpart, parts);
	return true;
}

double valuex(const Vector4 &p) { return p.x; }
Vector4 _valuex(const Vector4 &p, const double d) { return p + Vector4(d,0,0,0); }
double valuey(const Vector4 &p) { return p.y; }
Vector4 _valuey(const Vector4 &p, const double d) { return p + Vector4(0,d,0,0); }
double valuez(const Vector4 &p) { return p.z; }
Vector4 _valuez(const Vector4 &p, const double d) { return p + Vector4(0,0,d,0); }
double valuet(const Vector4 &p) { return p.t; }
Vector4 _valuet(const Vector4 &p, const double d) { return p + Vector4(0,0,0,d); }
double radius(const Vector4 &p) { return p.len(); }
Vector4 _radius(const Vector4 &p, const double d) { return p + d * p.unit(); }

inline double polar(const Vector2 &p) {
	if(p.y >= 0.0) {
		if(p.x > p.y) return atan(p.y / p.x); // 0..1
		if(p.x > -p.y) return 0.5 * PI - atan(p.x / p.y); // 1..3
		return PI + atan(p.y / p.x); // 3..4
	}
	if(p.x < p.y) return -PI + atan(p.y / p.x); // -4..-3
	if(p.x < -p.y) return -0.5 * PI - atan(p.x / p.y); // -3..-1
	return atan(p.y / p.x); // -1..0
}
inline Vector2 _polar(const Vector2 &p, const double d) {
	const double cosi = cos(d);
	const double sini = sin(d);
	return Matrix2(cosi, -sini, sini, cosi) * p;
}
double polarxy(const Vector4 &p) { return polar(Vector2(p.x, p.y)); }
Vector4 _polarxy(const Vector4 &p, const double d) { const Vector2 v = _polar(Vector2(p.x, p.y), d); return Vector4(v.x, v.y, p.z, p.t); }
double polarxz(const Vector4 &p) { return polar(Vector2(p.x, p.z)); }
Vector4 _polarxz(const Vector4 &p, const double d) { const Vector2 v = _polar(Vector2(p.x, p.z), d); return Vector4(v.x, p.y, v.y, p.t); }
double polarxt(const Vector4 &p) { return polar(Vector2(p.x, p.t)); }
Vector4 _polarxt(const Vector4 &p, const double d) { const Vector2 v = _polar(Vector2(p.x, p.t), d); return Vector4(v.x, p.y, p.z, v.y); }
double polaryz(const Vector4 &p) { return polar(Vector2(p.y, p.z)); }
Vector4 _polaryz(const Vector4 &p, const double d) { const Vector2 v = _polar(Vector2(p.y, p.z), d); return Vector4(p.x, v.x, v.y, p.t); }
double polaryt(const Vector4 &p) { return polar(Vector2(p.y, p.t)); }
Vector4 _polaryt(const Vector4 &p, const double d) { const Vector2 v = _polar(Vector2(p.y, p.t), d); return Vector4(p.x, v.x, p.z, v.y); }
double polarzt(const Vector4 &p) { return polar(Vector2(p.z, p.t)); }
Vector4 _polarzt(const Vector4 &p, const double d) { const Vector2 v = _polar(Vector2(p.z, p.t), d); return Vector4(p.x, p.y, v.x, v.y); }

inline double azimuth(const double x, const Vector4 &p) { return asin(x / p.len()); }
inline Vector4 _azimuth(const Vector4 &x, const Vector4 &p, const double d) { 
	const double psq = p.lensq();
	const Vector4 q = x.dot(p) * p + psq * x; 
	return cos(d) * p + sin(d) * sqrt(psq / q.lensq()) * q;
}
double azimuthx(const Vector4 &p) { return azimuth(p.x, p); }
Vector4 _azimuthx(const Vector4 &p, const double d) { return _azimuth(Vector4(1,0,0,0), p, d); }
double azimuthy(const Vector4 &p) { return azimuth(p.y, p); }
Vector4 _azimuthy(const Vector4 &p, const double d) { return _azimuth(Vector4(0,1,0,0), p, d); }
double azimuthz(const Vector4 &p) { return azimuth(p.z, p); }
Vector4 _azimuthz(const Vector4 &p, const double d) { return _azimuth(Vector4(0,0,1,0), p, d); }
double azimutht(const Vector4 &p) { return azimuth(p.t, p); }
Vector4 _azimutht(const Vector4 &p, const double d) { return _azimuth(Vector4(0,0,0,1), p, d); }


double getDouble(const string &str, const double def) {
	char *rest;
	const double val = strtod(str.c_str(), &rest);
	if(rest == NULL) return def;
	return val;
}
uint getUint(const string &str, const uint def) {
	const double val = getDouble(str, double(def));
	if(val < 0.0) return 0;
	return uint(val + 0.5);
}

bool checkSplit(const string &line, const string &command, double (*func)(const Vector4 &), const Mesh &mesh, Buffer<uint> &npart, Buffer<uint> &epart, Buffer<uint> &fpart, Buffer<uint> &bpart, Buffer<uint> &qpart, uint &parts){
	const uint commsize = command.size();
	if(line.substr(0, commsize).compare(command) != 0) return false;
	const uint v = getUint(line.substr(commsize), 2);
	if(split(v, func, mesh, npart, epart, fpart, bpart, qpart, parts)) 
		cout << "Command " << line << " succeeded." << endl;
	else cout << "Command " << line << " failed!" << endl;
	return true;
}
bool checkRepeat(const string &line, const string &command, double (*func)(const Vector4 &), Vector4 (*_func)(const Vector4 &, const double), const Mesh &mesh, Buffer<uint> &npart, Buffer<uint> &epart, Buffer<uint> &fpart, Buffer<uint> &bpart, Buffer<uint> &qpart, const uint parts) {
	const uint commsize = command.size();
	if(line.substr(0, commsize).compare(command) != 0) return false;
	const size_t separ = line.find_first_of('_', commsize);
	if(separ == string::npos) cout << "Please type " << command << "{min}_{max}." << endl;
	const double v0 = getDouble(line.substr(commsize, separ), 0.0);
	const double v1 = getDouble(line.substr(separ+1), 0.0);
	if(repeat(v0, v1, func, _func, mesh, npart, epart, fpart, bpart, qpart, parts))
		cout << "Command " << line << " succeeded." << endl;
	else cout << "Command " << line << " failed!" << endl;
	return true;
}

int main(int argc, const char* argv[])
{
	if(argc <= 2)
	{
		cout << "Please insert mesh file path as the first argument and at least one split argument." << endl;
		return 0;
	}

	// run by command arguments
	// first load mesh
	Mesh mesh;
	if(mesh.loadJRMesh(argv[1])) cout << "Mesh loaded succesfully." << endl;
	else
	{
		cout << "Unable to load mesh " << argv[1] << endl;
		cout << "We create a cubic mesh for testing purposes. Please try again." << endl;
		BuilderMesh bmesh(3);
		bmesh.createGrid(Vector4(-1,-1,0,0), Vector4(1,1,0,0), 0.1);
		bmesh.fillRectangleFlags(Vector4(-0.25,-0.25,-0.25,-0.25), Vector4(0.25,0.25,0.25,0.25), 1);
		bmesh.fillRectangleFlags(Vector4(0.45,-0.05,-0.05,-0.05), Vector4(0.95,0.45,0.45,0.45), 2);
		bmesh.saveJRMesh(argv[1]);
		return 0;
	}

	uint i;
	uint parts = 1;
	Buffer<uint> npart(mesh.getNodeSize(), 0);
	Buffer<uint> epart(mesh.getEdgeSize(), 0);
	Buffer<uint> fpart(mesh.getFaceSize(), 0);
	Buffer<uint> bpart(mesh.getBodySize(), 0);
	Buffer<uint> qpart(mesh.getQuadSize(), 0);

	for(int ai=2; ai<argc; ai++) {
		string line = argv[ai];
		
		if(line.substr(0, 9).compare("flagsplit") == 0) {
			if(flagSplit(mesh, npart, epart, fpart, bpart, qpart, parts)) 
				cout << "Command " << line << " succeeded." << endl;
			else cout << "Command " << line << " failed!" << endl;
			continue;
		}

		if(checkSplit(line, "splitx_", valuex, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splity_", valuey, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitz_", valuez, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitt_", valuet, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitr_", radius, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitpolarxy_", polarxy, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitpolarxz_", polarxz, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitpolarxt_", polarxt, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitpolaryz_", polaryz, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitpolaryt_", polaryt, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitpolarzt_", polarzt, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitazimuthx_", azimuthx, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitazimuthy_", azimuthy, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitazimuthz_", azimuthz, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkSplit(line, "splitazimutht_", azimutht, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;

		if(checkRepeat(line, "repeatx_", valuex, _valuex, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeaty_", valuey, _valuey, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatz_", valuez, _valuez, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatt_", valuet, _valuet, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatr_", radius, _radius, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatpolarxy_", polarxy, _polarxy, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatpolarxz_", polarxz, _polarxz, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatpolarxt_", polarxt, _polarxt, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatpolaryz_", polaryz, _polaryz, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatpolaryt_", polaryt, _polaryt, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatpolarzt_", polarzt, _polarzt, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatazimuthx_", azimuthx, _azimuthx, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatazimuthy_", azimuthy, _azimuthy, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatazimuthz_", azimuthz, _azimuthz, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		if(checkRepeat(line, "repeatazimutht_", azimutht, _azimutht, mesh, npart, epart, fpart, bpart, qpart, parts)) continue;
		cout << "Unknown argument: " << line << "." << endl;
	}

	// save part meshes
	for(i=0; i<parts; i++) {
		PartMesh pmesh(i, parts);
		pmesh.createPart(mesh, npart, epart, fpart, bpart, qpart);
		Text path;
		path << argv[1] << "." << i;
		if(pmesh.saveJRMesh(path.str())) cout << "Saved part mesh " << path.str() << "." << endl;
		else cout << "Unable to save part mesh " << path.str() << "." << endl;

		// draw the mesh
		MeshDrawer drawer;
		const Vector3 vo(0,0,0);
		const Vector3 vp(-3,5,10);
		const Vector3 vx = TwoVector3(Vector3(0,1,0), vp-vo).dual().unit() / 2.0;
		const Vector3 vy = TwoVector3(vp-vo, vx).dual().unit() / 2.0;
		drawer.initPosition(Vector4(vp,0), Vector4(vo,0), Vector4(vx,0), Vector4(vy,0));
		drawer.initSvg(500, 500);
		//drawer.drawPrimalEdges(mesh, Vector3(1,0,0));
		drawer.drawBoundaryFaces(pmesh, Vector3(1,0.5,0.5));
		Text picpath;
		picpath << "mesh_" << i << ".svg";
		drawer.saveSvg(picpath.str());
	}
	return 0;
}

