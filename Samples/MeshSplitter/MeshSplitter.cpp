/**
 * Program for splitting a mesh by domain decomposition method.
 * Call with parameters 'mesh_file_name' 'split_command_1' 'split_command_2' ...
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#include "../../GFD/Mesh/PartMesh.hpp"
#include "../../GFD/Mesh/BuilderMesh.hpp"
#include "../../GFD/Types/Types.hpp"
#include <fstream>
#include <iostream>
#include <map>

using namespace gfd;

bool divide(const uint divide, double (*func)(const Vector4 &), const Mesh &mesh, Buffer<uint> &npart, uint &parts) {
	if(divide == 0) return false;
	if(divide == 1) return true;

	multimap<double,uint> tree;
	const uint nodes = mesh.getNodeSize();
	for(uint i=0; i<nodes; i++) tree.insert(pair<double,uint>(func(mesh.getNodePosition(i)), i));

	auto it = tree.begin();
	for(uint i=0; i<nodes; i++, it++) npart[it->second] += parts * uint(i * divide / nodes);
	parts *= divide;
	return true;
}

double lensq(const Vector4 &p) {
	return p.lensq();
}

double polarxy(const Vector4 &p) {
	if(p.x * p.x + p.y * p.y < 1e-13) return 0;
	if(p.y >= 0.0) {
		if(p.x > p.y) return p.y / p.x; // 0..1
		if(p.x > -p.y) return 2.0 - p.x / p.y; // 1..3
		return 4.0 + p.y / p.x; // 3..4
	}
	if(p.x < p.y) return 4.0 + p.y / p.x; // 4..5
	if(p.x < -p.y) return 6.0 - p.x / p.y; // 5..7
	return 8.0 + p.y / p.x; // 7..8
}
double polarxz(const Vector4 &p) {
	if(p.x * p.x + p.z * p.z < 1e-13) return 0;
	if(p.z >= 0.0) {
		if(p.x > p.z) return p.z / p.x; // 0..1
		if(p.x > -p.z) return 2.0 - p.x / p.z; // 1..3
		return 4.0 + p.z / p.x; // 3..4
	}
	if(p.x < p.z) return 4.0 + p.z / p.x; // 4..5
	if(p.x < -p.z) return 6.0 - p.x / p.z; // 5..7
	return 8.0 + p.z / p.x; // 7..8
}
double polarxt(const Vector4 &p) {
	if(p.x * p.x + p.t * p.t < 1e-13) return 0;
	if(p.t >= 0.0) {
		if(p.x > p.t) return p.t / p.x; // 0..1
		if(p.x > -p.t) return 2.0 - p.x / p.t; // 1..3
		return 4.0 + p.t / p.x; // 3..4
	}
	if(p.x < p.t) return 4.0 + p.t / p.x; // 4..5
	if(p.x < -p.t) return 6.0 - p.x / p.t; // 5..7
	return 8.0 + p.t / p.x; // 7..8
}
double polaryz(const Vector4 &p) {
	if(p.y * p.y + p.z * p.z < 1e-13) return 0;
	if(p.z >= 0.0) {
		if(p.y > p.z) return p.z / p.y; // 0..1
		if(p.y > -p.z) return 2.0 - p.y / p.z; // 1..3
		return 4.0 + p.z / p.y; // 3..4
	}
	if(p.y < p.z) return 4.0 + p.z / p.y; // 4..5
	if(p.y < -p.z) return 6.0 - p.y / p.z; // 5..7
	return 8.0 + p.z / p.y; // 7..8
}
double polaryt(const Vector4 &p) {
	if(p.y * p.y + p.t * p.t < 1e-13) return 0;
	if(p.t >= 0.0) {
		if(p.y > p.t) return p.t / p.y; // 0..1
		if(p.y > -p.t) return 2.0 - p.y / p.t; // 1..3
		return 4.0 + p.t / p.y; // 3..4
	}
	if(p.y < p.t) return 4.0 + p.t / p.y; // 4..5
	if(p.y < -p.t) return 6.0 - p.y / p.t; // 5..7
	return 8.0 + p.t / p.y; // 7..8
}
double polarzt(const Vector4 &p) {
	if(p.z * p.z + p.t * p.t < 1e-13) return 0;
	if(p.t >= 0.0) {
		if(p.z > p.t) return p.t / p.z; // 0..1
		if(p.z > -p.t) return 2.0 - p.z / p.t; // 1..3
		return 4.0 + p.t / p.z; // 3..4
	}
	if(p.z < p.t) return 4.0 + p.t / p.z; // 4..5
	if(p.z < -p.t) return 6.0 - p.z / p.t; // 5..7
	return 8.0 + p.t / p.z; // 7..8
}

double azimuthx(const Vector4 &p) {
	return p.x / sqrt(p.y * p.y + p.z * p.z + p.t * p.t);
}
double azimuthy(const Vector4 &p) {
	return p.y / sqrt(p.x * p.x + p.z * p.z + p.t * p.t);
}
double azimuthz(const Vector4 &p) {
	return p.z / sqrt(p.x * p.x + p.y * p.y + p.t * p.t);
}
double azimutht(const Vector4 &p) {
	return p.t / sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

double valuex(const Vector4 &p) {
	return p.x;
}
double valuey(const Vector4 &p) {
	return p.y;
}
double valuez(const Vector4 &p) {
	return p.z;
}
double valuet(const Vector4 &p) {
	return p.t;
}

uint getUint(const std::string &str, const uint def)
{
	char *rest;
	const double val = strtod(str.c_str(), &rest);
	if(rest == NULL) return def;
	if(val < 0.0) return 0;
	return uint(val + 0.5);
}

double getDouble(const std::string &str, const double def)
{
	char *rest;
	const double val = strtod(str.c_str(), &rest);
	if(rest == NULL) return def;
	return val;
}


int main(int argc, const char* argv[])
{
	if(argc <= 2)
	{
		std::cout << "Please insert mesh file path as the first argument and at least one split argument." << std::endl;
		return 0;
	}

	// run by command arguments
	// first load mesh
	Mesh mesh;
	if(mesh.loadJRMesh(argv[1])) std::cout << "Mesh loaded succesfully." << std::endl;
	else
	{
		std::cout << "Unable to load mesh " << argv[1] << std::endl;
		std::cout << "We create a cubic mesh for testing purposes. Please try again." << std::endl;
		BuilderMesh bmesh(3);
		bmesh.createGrid(Vector4(-1,-1,-1,0), Vector4(1,1,1,0), 0.1);
		bmesh.saveJRMesh(argv[1]);
		return 0;
	}

	uint parts = 1;
	Buffer<uint> npart(mesh.getNodeSize(), 0);

	uint i, j;
	//Vector3 pmin(-1e30, -1e30, -1e30);
	//Vector3 pmax(1e30, 1e30, 1e30);
	for(int ai=2; ai<argc; ai++) {
		std::string line = argv[ai];
		if(line.substr(0, 4).compare("flag") == 0) {
			uint iparts = 1;
			for(i=0; i<npart.size(); i++) {
				const uint flag = mesh.getNodeFlag(i);
				if(iparts <= flag) iparts = flag + 1; 
				npart[i] += parts * flag;
				mesh.setNodeFlag(i, 0);
			}
			parts = iparts;
			std::cout << "Division by node flags." << std::endl;
		}
		else if(line.substr(0, 2).compare("x_") == 0) {
			const uint val = getUint(line.substr(2), 2);
			if(divide(val, valuex, mesh, npart, parts)) std::cout << "X-split succeeded into " << val << " components." << std::endl;
			else std::cout << "X-split failed!" << std::endl;
		}
		else if(line.substr(0, 2).compare("y_") == 0) {
			const uint val = getUint(line.substr(2), 2);
			if(divide(val, valuey, mesh, npart, parts)) std::cout << "Y-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Y-split failed!" << std::endl;
		}
		else if(line.substr(0, 2).compare("z_") == 0) {
			const uint val = getUint(line.substr(2), 2);
			if(divide(val, valuez, mesh, npart, parts)) std::cout << "Z-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Z-split failed!" << std::endl;
		}
		else if(line.substr(0, 2).compare("t_") == 0) {
			const uint val = getUint(line.substr(2), 2);
			if(divide(val, valuet, mesh, npart, parts)) std::cout << "T-split succeeded into " << val << " components." << std::endl;
			else std::cout << "T-split failed!" << std::endl;
		}
		else if(line.substr(0, 9).compare("azimuthx_") == 0) {
			const uint val = getUint(line.substr(9), 2);
			if(divide(val, azimuthx, mesh, npart, parts)) std::cout << "Azimuth X-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Azimuth X-split failed!" << std::endl;
		}
		else if(line.substr(0, 9).compare("azimuthy_") == 0) {
			const uint val = getUint(line.substr(9), 2);
			if(divide(val, azimuthy, mesh, npart, parts)) std::cout << "Azimuth Y-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Azimuth Y-split failed!" << std::endl;
		}
		else if(line.substr(0, 9).compare("azimuthz_") == 0) {
			const uint val = getUint(line.substr(9), 2);
			if(divide(val, azimuthz, mesh, npart, parts)) std::cout << "Azimuth Z-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Azimuth Z-split failed!" << std::endl;
		}
		else if(line.substr(0, 9).compare("azimutht_") == 0) {
			const uint val = getUint(line.substr(9), 2);
			if(divide(val, azimutht, mesh, npart, parts)) std::cout << "Azimuth T-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Azimuth T-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("polarxy_") == 0) {
			const uint val = getUint(line.substr(8), 2);
			if(divide(val, polarxy, mesh, npart, parts)) std::cout << "Polar XY-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar XY-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("polarxz_") == 0) {
			const uint val = getUint(line.substr(8), 2);
			if(divide(val, polarxz, mesh, npart, parts)) std::cout << "Polar XZ-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar XZ-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("polarxt_") == 0) {
			const uint val = getUint(line.substr(8), 2);
			if(divide(val, polarxt, mesh, npart, parts)) std::cout << "Polar XT-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar XT-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("polaryz_") == 0) {
			const uint val = getUint(line.substr(8), 2);
			if(divide(val, polaryz, mesh, npart, parts)) std::cout << "Polar YZ-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar YZ-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("polaryt_") == 0) {
			const uint val = getUint(line.substr(8), 2);
			if(divide(val, polaryt, mesh, npart, parts)) std::cout << "Polar YT-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar YT-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("polarzt_") == 0) {
			const uint val = getUint(line.substr(8), 2);
			if(divide(val, polarzt, mesh, npart, parts)) std::cout << "Polar ZT-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar ZT-split failed!" << std::endl;
		}
		else if(line.substr(0, 7).compare("radius_") == 0) {
			const uint val = getUint(line.substr(7), 2);
			if(divide(val, lensq, mesh, npart, parts)) std::cout << "Radius-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Radius-split failed!" << std::endl;
		}

/*		else if(line.substr(0, 6).compare("polarx") == 0)
		{
			const uint val = getUint(line.substr(6), 2);
			if(mesh.divide(val, polarx)) std::cout << "Polar X-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar X-split failed!" << std::endl;
		}
		else if(line.substr(0, 6).compare("polary") == 0)
		{
			const uint val = getUint(line.substr(6), 2);
			if(mesh.divide(val, polary)) std::cout << "Polar Y-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar Y-split failed!" << std::endl;
		}
		else if(line.substr(0, 6).compare("polarz") == 0)
		{
			const uint val = getUint(line.substr(6), 2);
			if(mesh.divide(val, polarz)) std::cout << "Polar Z-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar Z-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("azimuthx") == 0)
		{
			const uint val = getUint(line.substr(8), 2);
			if(mesh.divide(val, polarx)) std::cout << "Azimuth X-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Azimuth X-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("azimuthy") == 0)
		{
			const uint val = getUint(line.substr(8), 2);
			if(mesh.divide(val, polary)) std::cout << "Azimuth Y-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Azimuth Y-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("azimuthz") == 0)
		{
			const uint val = getUint(line.substr(8), 2);
			if(mesh.divide(val, polarz)) std::cout << "Azimuth Z-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Azimuth Z-split failed!" << std::endl;
		}
		else if(line.substr(0, 4).compare("minx") == 0)
		{
			pmin.x = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 4).compare("maxx") == 0)
		{
			pmax.x = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 4).compare("miny") == 0)
		{
			pmin.y = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 4).compare("maxy") == 0)
		{
			pmax.y = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 4).compare("minz") == 0)
		{
			pmin.z = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 4).compare("maxz") == 0)
		{
			pmax.z = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 6).compare("repeat") == 0)
		{
			if(mesh.restrictArea(pmin, pmax)) std::cout << "Repeat area set between (" << pmin.x << ", " << pmin.y << ", " << pmin.z << ") and (" << pmax.x << ", " << pmax.y << ", " << pmax.z << ")." << std::endl;
			else std::cout << "Failed to set repeat area!" << std::endl;
		}
*/	}

	// parts must be greater that 1 for a split to occur
	if(parts == 1) {
		cout << "No split occurred." << endl;
		return 0;
	}


	// compute higher element parts from node parts
	Buffer<uint> epart(mesh.getEdgeSize(), 0);
	for(i=0; i<epart.size(); i++) { 
		const Buffer<uint> &ele = mesh.getEdgeNodes(i);
		for(j=0; j<ele.size(); j++) {
			if(epart[i] < npart[ele[j]]) epart[i] = npart[ele[j]];
		}
	}
	Buffer<uint> fpart(mesh.getFaceSize(), 0);
	for(i=0; i<fpart.size(); i++) { 
		const Buffer<uint> &ele = mesh.getFaceEdges(i);
		for(j=0; j<ele.size(); j++) {
			if(fpart[i] < epart[ele[j]]) fpart[i] = epart[ele[j]];
		}
	}
	Buffer<uint> bpart(mesh.getBodySize(), 0);
	for(i=0; i<bpart.size(); i++) { 
		const Buffer<uint> &ele = mesh.getBodyFaces(i);
		for(j=0; j<ele.size(); j++) {
			if(bpart[i] < fpart[ele[j]]) bpart[i] = fpart[ele[j]];
		}
	}
	Buffer<uint> qpart(mesh.getQuadSize(), 0);
	for(i=0; i<qpart.size(); i++) { 
		const Buffer<uint> &ele = mesh.getQuadBodies(i);
		for(j=0; j<ele.size(); j++) {
			if(qpart[i] < bpart[ele[j]]) qpart[i] = bpart[ele[j]];
		}
	}

	// save part meshes
	for(i=0; i<parts; i++) {
		PartMesh pmesh(i, parts);
		pmesh.createPart(mesh, npart, epart, fpart, bpart, qpart);
		Text path;
		path << argv[1] << "." << i;
		if(pmesh.saveJRMesh(path.str())) cout << "Saved part mesh " << path.str() << "." << endl;
		else cout << "Unable to save part mesh " << path.str() << "." << endl;
	}



/*	// run by command arguments
	// first load mesh
	MeshSplitter mesh;
	if(mesh.loadJRMesh(argv[1])) std::cout << "Mesh loaded succesfully." << std::endl;
	else
	{
		std::cout << "Unable to load mesh " << argv[1] << std::endl;
		return 0;
	}

	Vector3 pmin(-1e30, -1e30, -1e30);
	Vector3 pmax(1e30, 1e30, 1e30);
	for(int i=2; i<argc; i++)
	{
		std::string line = argv[i];
		if(line.substr(0, 4).compare("flag") == 0)
		{
			const uint val = getUint(line.substr(4), 1);
			uint pow = 1;
			while(pow < val) pow *= 2;
			if(pow > 1024 || pow != val)
			{
				std::cout << "Flag-split value must be a power of 2 and at most 1024." << std::endl;
				return 0;
			}
			if(mesh.divideByFlag(val)) std::cout << "Flag-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Flag-split failed!" << std::endl;
		}
		else if(line.substr(0, 1).compare("x") == 0)
		{
			const uint val = getUint(line.substr(1), 2);
			if(mesh.divide(val, valuex)) std::cout << "X-split succeeded into " << val << " components." << std::endl;
			else std::cout << "X-split failed!" << std::endl;
		}
		else if(line.substr(0, 1).compare("y") == 0)
		{
			const uint val = getUint(line.substr(1), 2);
			if(mesh.divide(val, valuey)) std::cout << "Y-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Y-split failed!" << std::endl;
		}
		else if(line.substr(0, 1).compare("z") == 0)
		{
			const uint val = getUint(line.substr(1), 2);
			if(mesh.divide(val, valuez)) std::cout << "Z-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Z-split failed!" << std::endl;
		}
		else if(line.substr(0, 6).compare("polarx") == 0)
		{
			const uint val = getUint(line.substr(6), 2);
			if(mesh.divide(val, polarx)) std::cout << "Polar X-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar X-split failed!" << std::endl;
		}
		else if(line.substr(0, 6).compare("polary") == 0)
		{
			const uint val = getUint(line.substr(6), 2);
			if(mesh.divide(val, polary)) std::cout << "Polar Y-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar Y-split failed!" << std::endl;
		}
		else if(line.substr(0, 6).compare("polarz") == 0)
		{
			const uint val = getUint(line.substr(6), 2);
			if(mesh.divide(val, polarz)) std::cout << "Polar Z-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Polar Z-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("azimuthx") == 0)
		{
			const uint val = getUint(line.substr(8), 2);
			if(mesh.divide(val, polarx)) std::cout << "Azimuth X-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Azimuth X-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("azimuthy") == 0)
		{
			const uint val = getUint(line.substr(8), 2);
			if(mesh.divide(val, polary)) std::cout << "Azimuth Y-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Azimuth Y-split failed!" << std::endl;
		}
		else if(line.substr(0, 8).compare("azimuthz") == 0)
		{
			const uint val = getUint(line.substr(8), 2);
			if(mesh.divide(val, polarz)) std::cout << "Azimuth Z-split succeeded into " << val << " components." << std::endl;
			else std::cout << "Azimuth Z-split failed!" << std::endl;
		}
		else if(line.substr(0, 4).compare("minx") == 0)
		{
			pmin.x = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 4).compare("maxx") == 0)
		{
			pmax.x = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 4).compare("miny") == 0)
		{
			pmin.y = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 4).compare("maxy") == 0)
		{
			pmax.y = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 4).compare("minz") == 0)
		{
			pmin.z = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 4).compare("maxz") == 0)
		{
			pmax.z = getDouble(line.substr(4), 0.0);
		}
		else if(line.substr(0, 6).compare("repeat") == 0)
		{
			if(mesh.restrictArea(pmin, pmax)) std::cout << "Repeat area set between (" << pmin.x << ", " << pmin.y << ", " << pmin.z << ") and (" << pmax.x << ", " << pmax.y << ", " << pmax.z << ")." << std::endl;
			else std::cout << "Failed to set repeat area!" << std::endl;
		}
	}

	if(mesh.saveSplit(argv[1])) std::cout << "Splitted mesh saved succesfully." << std::endl;
	else std::cout << "Unable to save mesh " << argv[1] << std::endl;
*/
	return 0;
}

