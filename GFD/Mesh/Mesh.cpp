#include "Mesh.hpp"
#include "../Types/Text.hpp"
#include <string>
#include <fstream>
#include <iostream>

using namespace gfd;

Mesh::Mesh(const uint dim)
{
	if(dim < 1) m_dim = 1;
	else if(dim > 4) m_dim = 4;
	else m_dim = dim;

	m_nsize = 0;
	m_esize = 0;
	m_fsize = 0;
	m_bsize = 0;
	m_qsize = 0;
}

void Mesh::clear()
{
	resizeNodeBuffer(0);
	m_nsize = 0;
	resizeEdgeBuffer(0);
	m_esize = 0;
	resizeFaceBuffer(0);
	m_fsize = 0;
	resizeBodyBuffer(0);
	m_bsize = 0;
	resizeQuadBuffer(0);
	m_qsize = 0;

	// flags
	m_nflag.clear(); // node flags (optional)
	m_eflag.clear(); // edge flags (optional)
	m_fflag.clear(); // face flags (optional)
	m_bflag.clear(); // body flags (optional)
	m_qflag.clear(); // body flags (optional)

	// circumcenter computation
	m_m.clear(); // symmetric matrix to determine dot product (optional)
	m_w.clear(); // node weights (optional)
}

void Mesh::swap(Mesh &mesh)
{
	const uint mdim = m_dim;
	m_dim = mesh.m_dim;
	mesh.m_dim = mdim;

	m_p.swap(mesh.m_p);

	const uint mnsize = m_nsize;
	m_nsize = mesh.m_nsize;
	mesh.m_nsize = mnsize;
	m_n.swap(mesh.m_n);

	const uint mesize = m_esize;
	m_esize = mesh.m_esize;
	mesh.m_esize = mesize;
	m_e.swap(mesh.m_e);

	const uint mfsize = m_fsize;
	m_fsize = mesh.m_fsize;
	mesh.m_fsize = mfsize;
	m_f.swap(mesh.m_f);

	const uint mbsize = m_bsize;
	m_bsize = mesh.m_bsize;
	mesh.m_bsize = mbsize;
	m_b.swap(mesh.m_b);

	const uint mqsize = m_qsize;
	m_qsize = mesh.m_qsize;
	mesh.m_qsize = mqsize;
	m_q.swap(mesh.m_q);

	// flags
	m_nflag.swap(mesh.m_nflag);
	m_eflag.swap(mesh.m_eflag);
	m_fflag.swap(mesh.m_fflag);
	m_bflag.swap(mesh.m_bflag);
	m_qflag.swap(mesh.m_qflag);

	// circumcenter computation
	m_m.swap(mesh.m_m);
	m_w.swap(mesh.m_w);
}

/*bool Mesh::loadMesh(const std::string &path)
{
	std::ifstream fs(path.c_str(), std::ios::in);
	if(fs.fail() != 0) return false;

	clear();
	uint i, j;

	const uint strlen = 1000;
	char str[strlen];
	std::vector<double> data;

  // read nodes
	while(fs.getline(str, strlen))
	{
		const std::string line = str;
		if(line.find("#") == line.npos && line.find("nodes") != line.npos) break;
	}
	fs.getline(str, strlen);
	data = readLine(str);
	if(!data.empty())
	{
		m_n.resize(uint(data[0] + 0.5));
		for(i=0; i<m_n.size(); i++)
		{
			fs.getline(str, strlen);
			data = readLine(str);
			if(data.size() < 5)
			{
				fs.close();
				return false;
			}
			m_n[i] = new Node(Vector3(data[0], data[1], data[2]));
			m_n[i]->setIndex(i);
			m_n[i]->setWeight(data[3]);
			m_n[i]->setFlag(uint(data[4] + 0.5));
		}
	}

	// read edges
	while(fs.getline(str, strlen))
	{
		const std::string line = str;
		if(line.find("#") == line.npos && line.find("edges") != line.npos) break;
	}
	fs.getline(str, strlen);
	data = readLine(str);
	if(!data.empty())
	{
		m_e.resize(uint(data[0] + 0.5));
		for(i=0; i<m_e.size(); i++)
		{
			fs.getline(str, strlen);
			data = readLine(str);
			if(data.size() < 3)
			{
				fs.close();
				return false;
			}
			Node *const n0 = m_n[uint(data[0] - 0.5)];
			Node *const n1 = m_n[uint(data[1] - 0.5)];
			m_e[i] = new Edge(n0, n1);
			m_e[i]->setIndex(i);
			m_e[i]->setFlag(uint(data[2] + 0.5));
			n0->addEdge(m_e[i]);
			n1->addEdge(m_e[i]);
		}
	}

  // read faces
	while(fs.getline(str, strlen))
	{
		const std::string line = str;
		if(line.find("#") == line.npos && line.find("faces") != line.npos) break;
	}
	fs.getline(str, strlen);
	data = readLine(str);
	if(!data.empty())
	{
		m_f.resize(uint(data[0] + 0.5));
		for(i=0; i<m_f.size(); i++)
		{
			fs.getline(str, strlen);
			data = readLine(str);
			if(data.size() < 4)
			{
				fs.close();
				return false;
			}
			Buffer<Edge *> e(data.size() - 1);
			for(j=0; j<e.size(); j++)
			{
				e[j] = m_e[uint(data[j] - 0.5)];
			}
			m_f[i] = new Face(e);
			m_f[i]->setIndex(i);
			m_f[i]->setFlag(uint(data[j] + 0.5));
			for(j=0; j<e.size(); j++)
			{
				e[j]->addFace(m_f[i]);
			}
		}
	}

  // read bodies
	while(fs.getline(str, strlen))
	{
		const std::string line = str;
		if(line.find("#") == line.npos && line.find("bodies") != line.npos) break;
	}
	fs.getline(str, strlen);
	data = readLine(str);
	if(!data.empty())
	{
		m_b.resize(uint(data[0] + 0.5));
		for(i=0; i<m_b.size(); i++)
		{
			fs.getline(str, strlen);
			data = readLine(str);
			if(data.size() < 5)
			{
				fs.close();
				return false;
			}
			Buffer<Face *> f(data.size() - 1);
			for(j=0; j<f.size(); j++)
			{
				f[j] = m_f[uint(data[j] - 0.5)];
			}
			m_b[i] = new Body(f);
			m_b[i]->setIndex(i);
			m_b[i]->setFlag(uint(data[j] + 0.5));
			for(j=0; j<f.size(); j++)
			{
				f[j]->addBody(m_b[i]);
			}
		}
	}

	fs.close();
	return true;
}
bool Mesh::saveMesh(const std::string &path) const
{
	std::ofstream fs(path.c_str(), std::ios_base::trunc);
	if(fs.fail() != 0) return false;

	uint i, j;

	fs.precision(10);

	fs << "nodes" << std::endl;
	fs << getNodeSize() << std::endl;
	for(i=0; i<getNodeSize(); i++)
	{
		const Vector4 &p = getNodePosition(i);
		fs << p.x << "\t" << p.y << "\t" << p.z << "\t" << p.t << "\t";
		fs << getNodeWeight(i) << "\t" << getNodeFlag(i) << std::endl;
	}
	fs << std::endl;

	fs << "edges" << std::endl;
	fs << getEdgeSize() << std::endl;
	for(i=0; i<getEdgeSize(); i++)
	{
		const Buffer<uint> &par = getEdgeNodes(i);
		for(j=0; j<par.size(); j++) fs << par[j] + 1 << "\t";
		fs << getEdgeFlag(i) << std::endl;
	}
	fs << std::endl;

	fs << "faces" << std::endl;
	fs << getFaceSize() << std::endl;
	for(i=0; i<getFaceSize(); i++)
	{
		const Buffer<uint> &par = getFaceEdges(i);
		for(j=0; j<par.size(); j++) fs << par[j] + 1 << "\t";
		fs << getFaceFlag(i) << std::endl;
	}
	fs << std::endl;

	fs << "bodies" << std::endl;
	fs << getBodySize() << std::endl;
	for(i=0; i<getBodySize(); i++)
	{
		const Buffer<uint> &par = getBodyFaces(i);
		for(j=0; j<par.size(); j++) fs << par[j] + 1 << "\t";
		fs << getBodyFlag(i) << std::endl;
	}
	fs << std::endl;

	fs << "quads" << std::endl;
	fs << getQuadSize() << std::endl;
	for(i=0; i<getQuadSize(); i++)
	{
		const Buffer<uint> &par = getQuadBodies(i);
		for(j=0; j<par.size(); j++) fs << par[j] + 1 << "\t";
		fs << getQuadFlag(i) << std::endl;
	}
	fs << std::endl;

	fs.close();

	return true;
}
*/
bool Mesh::loadJRMesh(const std::string &path)
{
	std::ifstream fs(path.c_str(), std::ios::binary | std::ios::in);
	if(fs.fail()) return false;

	// header
	char id[4];
	fs.read(id, 4);
	if(std::string(id).compare(0, 4, "JRM1") != 0)
	{
		fs.close();
		return false;
	}

	clear();
	uint i, size, index;

	// dimension
	fs.read((char*)&m_dim, sizeof(uint));

	// node positions
	fs.read((char*)&m_nsize, sizeof(uint));
	resizeNodeBuffer(m_nsize);
	fs.read((char*)&m_p[0], m_dim * m_nsize * sizeof(double));

	// edges
	fs.read((char*)&m_esize, sizeof(uint));
	resizeEdgeBuffer(m_esize);
	for(i=0; i<m_esize; i++)
	{
		fs.read((char*)&index, sizeof(uint));
		const uint n0 = index;
		fs.read((char*)&index, sizeof(uint));
		const uint n1 = index;
		m_e[i].n.resize(2);
		m_e[i].n[0] = n0;
		m_n[n0].e.push_back(i);
		m_e[i].n[1] = n1;
		m_n[n1].e.push_back(i);
	}

	// faces
	if(m_dim > 1)
	{
		fs.read((char*)&m_fsize, sizeof(uint));
		resizeFaceBuffer(m_fsize);
		for(i=0; i<m_fsize; i++)
		{
			fs.read((char*)&size, sizeof(uint));
			m_f[i].e.resize(uint(size));
			for(uint j=0; j<m_f[i].e.size(); j++)
			{
				fs.read((char*)&index, sizeof(uint));
				m_f[i].e[j] = index;
				m_e[index].f.push_back(i);
			}
			orderFaceEdges(i);
		}
	}

	// bodies
	if(m_dim > 2)
	{
		fs.read((char*)&m_bsize, sizeof(uint));
		resizeBodyBuffer(m_bsize);
		for(i=0; i<m_bsize; i++)
		{
			fs.read((char*)&size, sizeof(uint));
			m_b[i].f.resize(uint(size));
			for(uint j=0; j<m_b[i].f.size(); j++)
			{
				fs.read((char*)&index, sizeof(uint));
				m_b[i].f[j] = index;
				m_f[index].b.push_back(i);
			}
			orderBodyFaces(i);
		}
	}

	// quads
	if(m_dim > 3)
	{
		fs.read((char*)&m_qsize, sizeof(uint));
		resizeQuadBuffer(m_qsize);
		for(i=0; i<m_qsize; i++)
		{
			fs.read((char*)&size, sizeof(uint));
			m_q[i].b.resize(uint(size));
			for(uint j=0; j<m_q[i].b.size(); j++)
			{
				fs.read((char*)&index, sizeof(uint));
				m_q[i].b[j] = index;
				m_b[index].q.push_back(i);
			}
			orderQuadBodies(i);
		}
	}

	// flags
	fs.read((char*)&size, sizeof(uint));
	m_nflag.resize(size);
	if(size > 0) fs.read((char*)&m_nflag[0], size * sizeof(uint));
	fs.read((char*)&size, sizeof(uint));
	m_eflag.resize(size);
	if(size > 0) fs.read((char*)&m_eflag[0], size * sizeof(uint));
	fs.read((char*)&size, sizeof(uint));
	m_fflag.resize(size);
	if(size > 0) fs.read((char*)&m_fflag[0], size * sizeof(uint));
	fs.read((char*)&size, sizeof(uint));
	m_bflag.resize(size);
	if(size > 0) fs.read((char*)&m_bflag[0], size * sizeof(uint));
	fs.read((char*)&size, sizeof(uint));
	m_qflag.resize(size);
	if(size > 0) fs.read((char*)&m_qflag[0], size * sizeof(uint));

	// circumcenter computation
	fs.read((char*)&size, sizeof(uint));
	m_m.resize(size);
	if(size > 0) fs.read((char*)&m_m[0], size * sizeof(double));
	fs.read((char*)&size, sizeof(uint));
	m_w.resize(size);
	if(size > 0) fs.read((char*)&m_w[0], size * sizeof(double));

	fs.close();
	return true;
}

bool Mesh::saveJRMesh(const std::string &path) const
{
	std::ofstream fs(path.c_str(), std::ios_base::binary | std::ios::trunc);
	if(fs.fail()) return false;

	// header
	const char *id = "JRM1";
	fs.write(id, 4);

	uint i, j, size;

	// dimension
	fs.write((char*)&m_dim, sizeof(uint));

	// node positions
	fs.write((char*)&m_nsize, sizeof(uint));
	fs.write((char*)&m_p[0], m_dim * m_nsize * sizeof(double));

	// edges
	fs.write((char*)&m_esize, sizeof(uint));
	for(i=0; i<m_esize; i++)
	{
		const Buffer<uint> &n = getEdgeNodes(i);
		fs.write((char*)&n[0], sizeof(uint));
		fs.write((char*)&n[1], sizeof(uint));
	}

	// faces
	if(m_dim > 1)
	{
		fs.write((char*)&m_fsize, sizeof(uint));
		for(i=0; i<m_fsize; i++)
		{
			const Buffer<uint> &e = getFaceEdges(i);
			size = e.size();
			fs.write((char*)&size, sizeof(uint));
			for(j=0; j<size; j++) fs.write((char*)&e[j], sizeof(uint));
		}
	}

	// bodies
	if(m_dim > 2)
	{
		fs.write((char*)&m_bsize, sizeof(uint));
		for(i=0; i<m_bsize; i++)
		{
			const Buffer<uint> &f = getBodyFaces(i);
			size = f.size();
			fs.write((char*)&size, sizeof(uint));
			for(j=0; j<size; j++) fs.write((char*)&f[j], sizeof(uint));
		}
	}

	// quads
	if(m_dim > 3)
	{
		fs.write((char*)&m_qsize, sizeof(uint));
		for(i=0; i<m_qsize; i++)
		{
			const Buffer<uint> &b = getQuadBodies(i);
			size = b.size();
			fs.write((char*)&size, sizeof(uint));
			for(j=0; j<size; j++) fs.write((char*)&b[j], sizeof(uint));
		}
	}

	// flags
	size = m_nflag.size();
	while(size > 0 && m_nflag[size - 1] == 0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_nflag[0], size * sizeof(uint));
	size = m_eflag.size();
	while(size > 0 && m_eflag[size - 1] == 0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_eflag[0], size * sizeof(uint));
	size = m_fflag.size();
	while(size > 0 && m_fflag[size - 1] == 0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_fflag[0], size * sizeof(uint));
	size = m_bflag.size();
	while(size > 0 && m_bflag[size - 1] == 0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_bflag[0], size * sizeof(uint));
	size = m_qflag.size();
	while(size > 0 && m_qflag[size - 1] == 0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_qflag[0], size * sizeof(uint));

	// circumcenter computation
	size = m_m.size();
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_m[0], size * sizeof(double));
	size = m_w.size();
	while(size > 0 && m_w[size - 1] == 0.0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_w[0], size * sizeof(double));

	fs.close();
	return true;
}

void Mesh::writeStatistics(Text &text, const UintSet &flag) const
{
	uint i;

	// neighbour elements
	if(m_nsize > 0) text << "nodes:  " << m_nsize << std::endl;
	if(m_esize > 0) text << "edges:  " << m_esize << std::endl;
	if(m_fsize > 0) text << "faces:  " << m_fsize << std::endl;
	if(m_bsize > 0) text << "bodies: " << m_bsize << std::endl;
	if(m_qsize > 0) text << "quads:  " << m_qsize << std::endl;
	text << "euler:  " << int(m_nsize - m_esize + m_fsize - m_bsize + m_qsize) << std::endl << std::endl;

	Buffer<uint> links(8, 0);
	text << "links:  \t0\t1\t2\t3\t4\t5\t6\tmore\taverage";
	if(m_esize > 0)
	{
		uint sum = 0;
		uint num = 0;
		for(i=0; i<m_esize; i++)
		{
			if(!flag.includes(getEdgeFlag(i))) continue;
			const uint li = getEdgeNodes(i).size();
			links[li < links.size() ? li : links.size() - 1]++;
			num += li;
			sum++;
		}
		text << std::endl << "edge->node:";
		for(i=0; i<links.size(); i++)
		{
			text << "\t" << links[i];
			links[i] = 0;
		}
		text << "\t" << double(num) / double(sum);
	}
	if(m_fsize > 0)
	{
		uint sum = 0;
		uint num = 0;
		for(i=0; i<m_fsize; i++)
		{
			if(!flag.includes(getFaceFlag(i))) continue;
			const uint li = getFaceEdges(i).size();
			links[li < links.size() ? li : links.size() - 1]++;
			num += li;
			sum++;
		}
		text << std::endl << "face->edge:";
		for(i=0; i<links.size(); i++)
		{
			text << "\t" << links[i];
			links[i] = 0;
		}
		text << "\t" << double(num) / double(sum);
	}
	if(m_bsize > 0)
	{
		uint sum = 0;
		uint num = 0;
		for(i=0; i<m_bsize; i++)
		{
			if(!flag.includes(getBodyFlag(i))) continue;
			const uint li = getBodyFaces(i).size();
			links[li < links.size() ? li : links.size() - 1]++;
			num += li;
			sum++;
		}
		text << std::endl << "body->face:";
		for(i=0; i<links.size(); i++)
		{
			text << "\t" << links[i];
			links[i] = 0;
		}
		text << "\t" << double(num) / double(sum);
	}
	if(m_qsize > 0)
	{
		uint sum = 0;
		uint num = 0;
		for(i=0; i<m_qsize; i++)
		{
			if(!flag.includes(getQuadFlag(i))) continue;
			const uint li = getQuadBodies(i).size();
			links[li < links.size() ? li : links.size() - 1]++;
			num += li;
			sum++;
		}
		text << std::endl << "quad->body:";
		for(i=0; i<links.size(); i++)
		{
			text << "\t" << links[i];
			links[i] = 0;
		}
		text << "\t" << double(num) / double(sum);
	}
	if(m_qsize > 0)
	{
		uint sum = 0;
		uint num = 0;
		for(i=0; i<m_bsize; i++)
		{
			if(!flag.includes(getBodyFlag(i))) continue;
			const uint li = getBodyQuads(i).size();
			links[li < links.size() ? li : links.size() - 1]++;
			num += li;
			sum++;
		}
		text << std::endl << "body->quad:";
		for(i=0; i<links.size(); i++)
		{
			text << "\t" << links[i];
			links[i] = 0;
		}
		text << "\t" << double(num) / double(sum);
	}
	if(m_bsize > 0)
	{
		uint sum = 0;
		uint num = 0;
		for(i=0; i<m_fsize; i++)
		{
			if(!flag.includes(getFaceFlag(i))) continue;
			const uint li = getFaceBodies(i).size();
			links[li < links.size() ? li : links.size() - 1]++;
			num += li;
			sum++;
		}
		text << std::endl << "face->body:";
		for(i=0; i<links.size(); i++)
		{
			text << "\t" << links[i];
			links[i] = 0;
		}
		text << "\t" << double(num) / double(sum);
	}
	if(m_fsize > 0)
	{
		uint sum = 0;
		uint num = 0;
		for(i=0; i<m_esize; i++)
		{
			if(!flag.includes(getEdgeFlag(i))) continue;
			const uint li = getEdgeFaces(i).size();
			links[li < links.size() ? li : links.size() - 1]++;
			num += li;
			sum++;
		}
		text << std::endl << "edge->face:";
		for(i=0; i<links.size(); i++)
		{
			text << "\t" << links[i];
			links[i] = 0;
		}
		text << "\t" << double(num) / double(sum);
	}
	if(m_esize > 0)
	{
		uint sum = 0;
		uint num = 0;
		for(i=0; i<m_nsize; i++)
		{
			if(!flag.includes(getNodeFlag(i))) continue;
			const uint li = getNodeEdges(i).size();
			links[li < links.size() ? li : links.size() - 1]++;
			num += li;
			sum++;
		}
		text << std::endl << "node->edge:";
		for(i=0; i<links.size(); i++)
		{
			text << "\t" << links[i];
			links[i] = 0;
		}
		text << "\t" << double(num) / double(sum);
	}
	text << std::endl << std::endl;

	// Cell statistics
	text << "Cell statistics:  min\tmean\tmax" << std::endl;
	double vole = 0.0;
	double volf = 0.0;
	double volb = 0.0;
	double volq = 0.0;
	Vector3 sumx = Vector3(1e30,0,-1e30);
	Vector3 sumy = Vector3(1e30,0,-1e30);
	Vector3 sumz = Vector3(1e30,0,-1e30);
	Vector3 sumt = Vector3(1e30,0,-1e30);
	Vector3 summ = Vector3(1e30,0,-1e30);
	Vector3 sumd = Vector3(1e30,0,-1e30);
	Vector3 sumh = Vector3(1e30,0,-1e30);
	Vector3 sumw = Vector3(1e30,0,-1e30);
	uint count = 0;
	for(i=0; i<m_nsize; i++)
	{
		if(!flag.includes(getNodeFlag(i))) continue;
		count++;

		const Vector4 p = getNodePosition(i);
		const double w = getNodeWeight(i);
		if(p.x < sumx.x) sumx.x = p.x;
		sumx.y += p.x;
		if(p.x > sumx.z) sumx.z = p.x;
		if(p.y < sumy.x) sumy.x = p.y;
		sumy.y += p.y;
		if(p.y > sumy.z) sumy.z = p.y;
		if(p.z < sumz.x) sumz.x = p.z;
		sumz.y += p.z;
		if(p.z > sumz.z) sumz.z = p.z;
		if(p.t < sumt.x) sumt.x = p.t;
		sumt.y += p.t;
		if(p.t > sumt.z) sumt.z = p.t;
		if(w < summ.x) summ.x = w;
		summ.y += w;
		if(w > summ.z) summ.z = w;

		const double dv = getNodeDualVector(i).xyzt;
		double len = fabs(dv);
		if(len < sumd.x) sumd.x = len;
		sumd.y += len;
		if(len > sumd.z) sumd.z = len;

		len = getNodeHodge(i);
		if(len < sumh.x) sumh.x = len;
		sumh.y += len;
		if(len > sumh.z) sumh.z = len;

		if(len < sumw.x) sumw.x = len;
		sumw.y += len;
		if(len > sumw.z) sumw.z = len;
	}
	if(count > 0)
	{
		text << "node pos.x:       " << sumx.x << " \t" << sumx.y / double(count) << " \t" << sumx.z << std::endl;
		if(m_dim > 1) text << "node pos.y:       " << sumy.x << " \t" << sumy.y / double(count) << " \t" << sumy.z << std::endl;
		if(m_dim > 2) text << "node pos.z:       " << sumz.x << " \t" << sumz.y / double(count) << " \t" << sumz.z << std::endl;
		if(m_dim > 3) text << "node pos.t:       " << sumt.x << " \t" << sumt.y / double(count) << " \t" << sumt.z << std::endl;
		if(!m_w.empty()) text << "node weight:      " << summ.x << " \t" << summ.y / double(count) << " \t" << summ.z << std::endl;
		text << std::endl;
		text << "node unit:        1 \t1 \t1" << std::endl;
		if(m_esize == 0) text << "node dual unit:   ";
		else if(m_fsize == 0) text << "node dual length: ";
		else if(m_bsize == 0) text << "node dual area:   ";
		else text << "node dual volume: ";
		text << sumd.x << " \t" << sumd.y / double(count) << " \t" << sumd.z << std::endl;
		text << "node Hodge:       " << sumh.x << " \t" << sumh.y / double(count) << " \t" << sumh.z << std::endl;
		text << "node wedge:       " << sumw.x << " \t" << sumw.y / double(count) << " \t" << sumw.z << std::endl << std::endl;
	}
	Vector3 sum = Vector3(1e30,0,-1e30);
	sumd = Vector3(1e30,0,-1e30);
	sumh = Vector3(1e30,0,-1e30);
	sumw = Vector3(1e30,0,-1e30);
	count = 0;
	for(i=0; i<m_esize; i++)
	{
		if(!flag.includes(getEdgeFlag(i))) continue;
		count++;

		const Vector4 v = getEdgeVector(i);
		double len = v.len();
		if(len < sum.x) sum.x = len;
		sum.y += len;
		if(len > sum.z) sum.z = len;
		if(getEdgeFaces(i).size() < 2) vole += len;

		const ThreeVector4 dv = getEdgeDualVector(i);
		len = dv.len();
		if(len < sumd.x) sumd.x = len;
		sumd.y += len;
		if(len > sumd.z) sumd.z = len;

		len = getEdgeHodge(i);
		if(len < sumh.x) sumh.x = len;
		sumh.y += len;
		if(len > sumh.z) sumh.z = len;

		len = FourVector4(v, dv).xyzt;
		if(len < sumw.x) sumw.x = len;
		sumw.y += len;
		if(len > sumw.z) sumw.z = len;
	}
	if(count > 0)
	{
		text << "edge length:      " << sum.x << " \t" << sum.y / double(count) << " \t" << sum.z << std::endl;
		if(m_fsize == 0) text << "edge dual unit:   ";
		else if(m_bsize == 0) text << "edge dual length: ";
		else if(m_qsize == 0) text << "edge dual area:   ";
		else text << "edge dual volume: ";
		text << sumd.x << " \t" << sumd.y / double(count) << " \t" << sumd.z << std::endl;
		text << "edge Hodge:       " << sumh.x << " \t" << sumh.y / double(count) << " \t" << sumh.z << std::endl;
		text << "edge wedge:       " << sumw.x << " \t" << sumw.y / double(count) << " \t" << sumw.z << std::endl << std::endl;
	}
	sum = Vector3(1e30,0,-1e30);
	sumd = Vector3(1e30,0,-1e30);
	sumh = Vector3(1e30,0,-1e30);
	sumw = Vector3(1e30,0,-1e30);
	count = 0;
	for(i=0; i<m_fsize; i++)
	{
		if(!flag.includes(getFaceFlag(i))) continue;
		count++;

		const TwoVector4 v = getFaceVector(i);
		double len = v.len();
		if(len < sum.x) sum.x = len;
		sum.y += len;
		if(len > sum.z) sum.z = len;
		if(getFaceBodies(i).size() < 2) volf += len;

		const TwoVector4 dv = getFaceDualVector(i);
		len = dv.len();
		if(len < sumd.x) sumd.x = len;
		sumd.y += len;
		if(len > sumd.z) sumd.z = len;

		len = getFaceHodge(i);
		if(len < sumh.x) sumh.x = len;
		sumh.y += len;
		if(len > sumh.z) sumh.z = len;

		len = FourVector4(v, dv).xyzt;
		if(len < sumw.x) sumw.x = len;
		sumw.y += len;
		if(len > sumw.z) sumw.z = len;
	}
	if(count > 0)
	{
		text << "face area:        " << sum.x << " \t" << sum.y / double(count) << " \t" << sum.z << std::endl;
		if(m_bsize == 0) text << "face dual unit:   ";
		else if(m_qsize == 0) text << "face dual length: ";
		else text << "face dual area:   ";
		text << sumd.x << " \t" << sumd.y / double(count) << " \t" << sumd.z << std::endl;
		text << "face Hodge:       " << sumh.x << " \t" << sumh.y / double(count) << " \t" << sumh.z << std::endl;
		text << "face wedge:       " << sumw.x << " \t" << sumw.y / double(count) << " \t" << sumw.z << std::endl << std::endl;
	}
	sum = Vector3(1e30,0,-1e30);
	sumd = Vector3(1e30,0,-1e30);
	sumh = Vector3(1e30,0,-1e30);
	sumw = Vector3(1e30,0,-1e30);
	count = 0;
	for(i=0; i<m_bsize; i++)
	{
		if(!flag.includes(getBodyFlag(i))) continue;
		count++;

		const ThreeVector4 v = getBodyVector(i);
		double len = v.len();
		if(len < sum.x) sum.x = len;
		sum.y += len;
		if(len > sum.z) sum.z = len;
		if(getBodyQuads(i).size() < 2) volb += len;

		const Vector4 dv = getBodyDualVector(i);
		len = dv.len();
		if(len < sumd.x) sumd.x = len;
		sumd.y += len;
		if(len > sumd.z) sumd.z = len;

		len = getBodyHodge(i);
		if(len < sumh.x) sumh.x = len;
		sumh.y += len;
		if(len > sumh.z) sumh.z = len;

		len = FourVector4(v, dv).xyzt;
		if(len < sumw.x) sumw.x = len;
		sumw.y += len;
		if(len > sumw.z) sumw.z = len;
	}
	if(count > 0)
	{
		text << "body volume:      " << sum.x << " \t" << sum.y / double(count) << " \t" << sum.z << std::endl;
		if(m_qsize == 0) text << "body dual unit:   ";
		else text << "body dual length: ";
		text << sumd.x << " \t" << sumd.y / double(count) << " \t" << sumd.z << std::endl;
		text << "body Hodge:       " << sumh.x << " \t" << sumh.y / double(count) << " \t" << sumh.z << std::endl;
		text << "body wedge:       " << sumw.x << " \t" << sumw.y / double(count) << " \t" << sumw.z << std::endl << std::endl;
	}
	sum = Vector3(1e30,0,-1e30);
	sumh = Vector3(1e30,0,-1e30);
	count = 0;
	for(i=0; i<m_qsize; i++)
	{
		if(!flag.includes(getQuadFlag(i))) continue;
		count++;

		double len = fabs(getQuadVector(i).xyzt);
		if(len < sum.x) sum.x = len;
		sum.y += len;
		if(len > sum.z) sum.z = len;
		volq += len;

		len = getQuadHodge(i);
		if(len < sumh.x) sumh.x = len;
		sumh.y += len;
		if(len > sumh.z) sumh.z = len;
	}
	sumw = sum;
	if(count > 0)
	{
		text << "quad volume:      " << sum.x << " \t" << sum.y / double(count) << " \t" << sum.z << std::endl;
		text << "quad dual unit:   1 \t1 \t1" << std::endl;
		text << "quad Hodge:       " << sumh.x << " \t" << sumh.y / double(count) << " \t" << sumh.z << std::endl;
		text << "quad wedge:       " << sumw.x << " \t" << sumw.y / double(count) << " \t" << sumw.z << std::endl << std::endl;
	}

	if(m_esize > 0) text << "total edge:       " << vole << std::endl;
	if(m_fsize > 0) text << "total face:       " << volf << std::endl;
	if(m_bsize > 0) text << "total body:       " << volb << std::endl;
	if(m_qsize > 0) text << "total quad:       " << volq << std::endl;
}

Buffer<uint> Mesh::getNodeFaces(const uint n) const
{
	const Buffer<uint> &ne = getNodeEdges(n);
	Buffer<uint> f(ne.size());
	uint l = 0;
	for(uint i=0; i<ne.size(); i++)
	{
		const Buffer<uint> &ef = getEdgeFaces(ne[i]);
		for(uint j=0; j<ef.size(); j++) f.gatherOnce(ef[j], l);
	}
	f.resize(l);
	return f;
}

Buffer<uint> Mesh::getNodeBodies(const uint n) const
{
	const Buffer<uint> f = getNodeFaces(n);
	Buffer<uint> b(f.size());
	uint l = 0;
	for(uint i=0; i<f.size(); i++)
	{
		const Buffer<uint> &fb = getFaceBodies(f[i]);
		for(uint j=0; j<fb.size(); j++) b.gatherOnce(fb[j], l);
	}
	b.resize(l);
	return b;
}

Buffer<uint> Mesh::getNodeQuads(const uint n) const
{
	const Buffer<uint> b = getNodeBodies(n);
	Buffer<uint> q(b.size());
	uint l = 0;
	for(uint i=0; i<b.size(); i++)
	{
		const Buffer<uint> &bq = getBodyQuads(b[i]);
		for(uint j=0; j<bq.size(); j++) q.gatherOnce(bq[j], l);
	}
	q.resize(l);
	return q;
}

Buffer<uint> Mesh::getEdgeBodies(const uint e) const
{
	const Buffer<uint> &f = getEdgeFaces(e);
	Buffer<uint> b(f.size());
	uint l = 0;
	for(uint i=0; i<f.size(); i++)
	{
		const Buffer<uint> &fb = getFaceBodies(f[i]);
		for(uint j=0; j<fb.size(); j++) b.gatherOnce(fb[j], l);
	}
	b.resize(l);
	return b;
}

Buffer<uint> Mesh::getEdgeQuads(const uint e) const
{
	const Buffer<uint> b = getEdgeBodies(e);
	Buffer<uint> q(b.size());
	uint l = 0;
	for(uint i=0; i<b.size(); i++)
	{
		const Buffer<uint> &bq = getBodyQuads(b[i]);
		for(uint j=0; j<bq.size(); j++) q.gatherOnce(bq[j], l);
	}
	q.resize(l);
	return q;
}

Buffer<uint> Mesh::getFaceNodes(const uint f) const
{
	const Buffer<uint> &e = m_f[f].e;
	Buffer<uint> n(e.size());
	for(uint i=0; i<e.size(); i++)
	{
		const uint j = (i > 0 ? i : e.size()) - 1;
		n[i] = getEdgeIntersection(e[i], e[j]);
	}
	return n;
}

Buffer<uint> Mesh::getFaceQuads(const uint f) const
{
	const Buffer<uint> &b = getFaceBodies(f);
	Buffer<uint> q(b.size());
	uint l = 0;
	for(uint i=0; i<b.size(); i++)
	{
		const Buffer<uint> &bq = getBodyQuads(b[i]);
		for(uint j=0; j<bq.size(); j++) q.gatherOnce(bq[j], l);
	}
	q.resize(l);
	return q;
}

Buffer<uint> Mesh::getBodyNodes(const uint b) const
{
	const Buffer<uint> e = getBodyEdges(b);
	Buffer<uint> n(e.size());
	uint l = 0;
	for(uint i=0; i<e.size(); i++)
	{
		const Buffer<uint> &en = getEdgeNodes(e[i]);
		for(uint j=0; j<en.size(); j++) n.gatherOnce(en[j], l);
	}
	n.resize(l);
	return n;
}

Buffer<uint> Mesh::getBodyEdges(const uint b) const
{
	const Buffer<uint> &f = getBodyFaces(b);
	Buffer<uint> e(f.size());
	uint l = 0;
	for(uint i=0; i<f.size(); i++)
	{
		const Buffer<uint> &fe = getFaceEdges(f[i]);
		for(uint j=0; j<fe.size(); j++) e.gatherOnce(fe[j], l);
	}
	e.resize(l);
	return e;
}

Buffer<uint> Mesh::getQuadNodes(const uint q) const
{
	const Buffer<uint> e = getQuadEdges(q);
	Buffer<uint> n(e.size());
	uint l = 0;
	for(uint i=0; i<e.size(); i++)
	{
		const Buffer<uint> &en = getEdgeNodes(e[i]);
		for(uint j=0; j<en.size(); j++) n.gatherOnce(en[j], l);
	}
	n.resize(l);
	return n;
}

Buffer<uint> Mesh::getQuadEdges(const uint q) const
{
	const Buffer<uint> f = getQuadFaces(q);
	Buffer<uint> e(f.size());
	uint l = 0;
	for(uint i=0; i<f.size(); i++)
	{
		const Buffer<uint> &fe = getFaceEdges(f[i]);
		for(uint j=0; j<fe.size(); j++) e.gatherOnce(fe[j], l);
	}
	e.resize(l);
	return e;
}

Buffer<uint> Mesh::getQuadFaces(const uint q) const
{
	const Buffer<uint> &b = getQuadBodies(q);
	Buffer<uint> f(b.size());
	uint l = 0;
	for(uint i=0; i<b.size(); i++)
	{
		const Buffer<uint> &bf = getBodyFaces(b[i]);
		for(uint j=0; j<bf.size(); j++) f.gatherOnce(bf[j], l);
	}
	f.resize(l);
	return f;
}

uint Mesh::getEdgeOtherNode(const uint e, const uint n) const
{
	const Buffer<uint> &en = getEdgeNodes(e);
	if(n == en[0]) return en[1];
	if(n == en[1]) return en[0];
	return NONE;
}

sign Mesh::getEdgeIncidence(const uint e, const uint n) const
{
	const Buffer<uint> &en = getEdgeNodes(e);
	if(n == en[0]) return -1;
	if(n == en[1]) return 1;
	return 0;
}

void Mesh::orderFaceEdges(const uint f) // ensure the circular order of edges
{
	uint i, j;
	Buffer<uint> &e = m_f[f].e;
	if(e.size() <= 3) return;
	for(i=1; i+1<e.size(); i++)
	{
		const uint e0 = e[i-1];
		if(getEdgeIntersection(e0, e[i]) != NONE) continue;

		for(j=i+1; getEdgeIntersection(e0, e[j]) == NONE; j++);
		const uint ei = e[i];
		e[i] = e[j];
		e[j] = ei;
	}
}

sign Mesh::getFaceIncidence(const uint f, const uint e) const
{
	const Buffer<uint> &fe = getFaceEdges(f);
	if(e == fe[0])
	{
		const uint n = getEdgeNodes(e)[0];
		const Buffer<uint> &en = getEdgeNodes(fe[1]);
		if(n == en[0] || n == en[1]) return -1;
		return 1;
	}
	for(uint i=1; i<fe.size(); i++)
	{
		if(e == fe[i])
		{
			const uint n = getEdgeNodes(e)[0];
			const Buffer<uint> &en = getEdgeNodes(fe[i - 1]);
			if(n == en[0] || n == en[1]) return 1;
			return -1;
		}
	}
	return 0;
}

void Mesh::orderBodyFaces(const uint b) // order faces to optimize getBodyIncidence
{
	Buffer<uint> &f = m_b[b].f;
	if(f.size() <= 4) return;

	uint i, j;
	uint level0 = 0;
	uint level1 = 1;
	uint level2 = f.size();
	for(i=level1; i+1<f.size(); )
	{
		const uint fi = f[i];
		for(j=level0; j<level1 && getFaceIntersection(fi,f[j])==NONE; j++);
		if(j < level1) ++i; // connection found
		else // connection not found -> put f[i] to the back
		{
			f[i] = f[--level2];
			f[level2] = fi;
		}
		if(i == level2) // next level
		{
			level0 = level1;
			level1 = level2;
			level2 = f.size();
		}
	}
}

sign Mesh::getBodyIncidence(const uint b, const uint f) const
{
	const Buffer<uint> &bf = getBodyFaces(b);

	// is f the first face in the list?
	if(f == bf[0]) return 1;

	// find f from the list
	uint i = 1;
	while(f != bf[i]) { if(++i >= bf.size()) return 0; }

	// find the incidence recursively
	sign res = 1;
	for(uint j=0; j<i; )
	{
		const uint e = getFaceEdges(bf[i]).getFirstIntersection(getFaceEdges(bf[j]), NONE);
		if(e == NONE) { ++j; continue; }

		res *= -getFaceIncidence(bf[i], e) * getFaceIncidence(bf[j], e);
		if(j == 0) return res;
		i = j;
		j = 0;
	}
	return 0;
}

void Mesh::orderQuadBodies(const uint q) // order bodies to optimize getQuadIncidence
{
	Buffer<uint> &b = m_q[q].b;
	if(b.size() <= 5) return;

	uint i, j;
	uint level0 = 0;
	uint level1 = 1;
	uint level2 = b.size();
	for(i=level1; i+1<b.size(); )
	{
		const uint bi = b[i];
		for(j=level0; j<level1 && getBodyIntersection(bi,b[j])==NONE; j++);
		if(j < level1) ++i; // connection found
		else // connection not found -> put f[i] to the back
		{
			b[i] = b[--level2];
			b[level2] = bi;
		}
		if(i == level2) // next level
		{
			level0 = level1;
			level1 = level2;
			level2 = b.size();
		}
	}
}

sign Mesh::getQuadIncidence(const uint q, const uint b) const
{
	const Buffer<uint> &qb = getQuadBodies(q);

	// is b the first body in the list?
	if(b == qb[0]) return 1;

	// find b from the list
	uint i = 1;
	while(b != qb[i]) { if(++i >= qb.size()) return 0; }

	// find the incidence recursively
	sign res = 1.0;
	for(uint j=0; j<i; )
	{
		const uint f = getBodyFaces(qb[i]).getFirstIntersection(getBodyFaces(qb[j]), NONE);
		if(f == NONE) { ++j; continue; }

		res *= -getBodyIncidence(qb[i], f) * getBodyIncidence(qb[j], f);
		if(j == 0) return res;
		i = j;
		j = 0;
	}
	return 0;
}

Vector4 Mesh::getNodePosition(const uint n) const
{
	if(m_dim == 1) return Vector4(getNodePosition1(n),0,0,0);
	if(m_dim == 2) return Vector4(getNodePosition2(n),0,0);
	if(m_dim == 3) return Vector4(getNodePosition3(n),0);
	return getNodePosition4(n);
}

double Mesh::getEdgePosition1(const uint e) const
{
	const Buffer<uint> &n = m_e[e].n;
	const double w = getNodeWeight(n[1]) - getNodeWeight(n[0]);
	if(w == 0.0) return 0.5 * (getNodePosition1(n[0]) + getNodePosition1(n[1]));
	const double p = getNodePosition1(n[0]);
	const double d = getNodePosition1(n[1]) - p;
	return p + 0.5 * (d + w / getTransformed1(d));
}
Vector2 Mesh::getEdgePosition2(const uint e) const
{
	const Buffer<uint> &n = m_e[e].n;
	const double w = getNodeWeight(n[1]) - getNodeWeight(n[0]);
	if(w == 0.0) return 0.5 * (getNodePosition2(n[0]) + getNodePosition2(n[1]));
	const Vector2 p = getNodePosition2(n[0]);
	const Vector2 d = getNodePosition2(n[1]) - p;
	return p + 0.5 * (1.0 + w / d.dot(getTransformed2(d))) * d;
}
Vector3 Mesh::getEdgePosition3(const uint e) const
{
	const Buffer<uint> &n = m_e[e].n;
	const double w = getNodeWeight(n[1]) - getNodeWeight(n[0]);
	if(w == 0.0) return 0.5 * (getNodePosition3(n[0]) + getNodePosition3(n[1]));
	const Vector3 p = getNodePosition3(n[0]);
	const Vector3 d = getNodePosition3(n[1]) - p;
	return p + 0.5 * (1.0 + w / d.dot(getTransformed3(d))) * d;
}
Vector4 Mesh::getEdgePosition4(const uint e) const
{
	const Buffer<uint> &n = m_e[e].n;
	const double w = getNodeWeight(n[1]) - getNodeWeight(n[0]);
	if(w == 0.0) return 0.5 * (getNodePosition4(n[0]) + getNodePosition4(n[1]));
	const Vector4 p = getNodePosition4(n[0]);
	const Vector4 d = getNodePosition4(n[1]) - p;
	return p + 0.5 * (1.0 + w / d.dot(getTransformed4(d))) * d;
}
Vector4 Mesh::getEdgePosition(const uint e) const
{
	if(m_dim == 1) return Vector4(getEdgePosition1(e),0,0,0);
	if(m_dim == 2) return Vector4(getEdgePosition2(e),0,0);
	if(m_dim == 3) return Vector4(getEdgePosition3(e),0);
	return getEdgePosition4(e);
}

Vector2 Mesh::getFacePosition2(const uint f) const
{
	const Buffer<uint> &e = getFaceEdges(f);
	const Vector2 p = getEdgePosition2(e[0]);
	SymMatrix2 mA(ZEROSYMMATRIX2);
	Vector2 vB(0,0);
	for(uint i=0; i<e.size(); i++)
	{
		const Vector2 vi = getTransformed2(getEdgeVector2(e[i]));
		mA += vi.outerProduct();
		if(i > 0) vB += vi * vi.dot(getEdgePosition2(e[i]) - p);
	}
	return p + mA.inverse() * vB;
}
Vector3 Mesh::getFacePosition3(const uint f) const
{
	const Buffer<uint> &e = getFaceEdges(f);
	const Vector3 p = getEdgePosition3(e[0]);
	SymMatrix3 mA(ZEROSYMMATRIX3);
	Vector3 vB(0,0,0);
	for(uint i=0; i<e.size(); i++)
	{
		const Vector3 vi = getTransformed3(getEdgeVector3(e[i]));
		mA += vi.outerProduct();
		if(i > 0) vB += vi * vi.dot(getEdgePosition3(e[i]) - p);
	}
	mA += getFaceVector3(f).dual().outerProduct();
	return p + mA.inverse() * vB;
}
Vector4 Mesh::getFacePosition4(const uint f) const
{
	const Buffer<uint> &e = getFaceEdges(f);
	const Vector4 p = getEdgePosition4(e[0]);
	SymMatrix4 mA(ZEROSYMMATRIX4);
	Vector4 vB(0,0,0,0);
	for(uint i=0; i<e.size(); i++)
	{
		const Vector4 vi = getTransformed4(getEdgeVector4(e[i]));
		mA += vi.outerProduct();
		if(i > 0) vB += vi * vi.dot(getEdgePosition4(e[i]) - p);
	}
	const TwoVector4 v = getFaceVector(f);
	mA += SymMatrix4(v.yz*v.yz+v.yt*v.yt+v.zt*v.zt,-v.xz*v.yz-v.xt*v.yt,v.xz*v.xz+v.xt*v.xt+v.zt*v.zt,
		v.xy*v.yz-v.xt*v.zt,-v.yt*v.zt-v.xy*v.xz,v.xy*v.xy+v.xt*v.xt+v.yt*v.yt,
		v.xz*v.zt+v.xy*v.yt,v.yz*v.zt-v.xy*v.xt,-v.yz*v.yt-v.xz*v.xt,v.xy*v.xy+v.xz*v.xz+v.yz*v.yz);
	return p + mA.inverse() * vB;
}
Vector4 Mesh::getFacePosition(const uint f) const
{
	if(m_dim == 2) return Vector4(getFacePosition2(f),0,0);
	if(m_dim == 3) return Vector4(getFacePosition3(f),0);
	return getFacePosition4(f);
}

Vector3 Mesh::getBodyPosition3(const uint b) const
{
	const Buffer<uint> e = getBodyEdges(b);
	const Vector3 p = getEdgePosition3(e[0]);
	SymMatrix3 mA(ZEROSYMMATRIX3);
	Vector3 vB(0,0,0);
	for(uint i=0; i<e.size(); i++)
	{
		const Vector3 vi = getTransformed3(getEdgeVector3(e[i]));
		mA += vi.outerProduct();
		if(i > 0) vB += vi * vi.dot(getEdgePosition3(e[i]) - p);
	}
	return p + mA.inverse() * vB;
}
Vector4 Mesh::getBodyPosition4(const uint b) const
{
	const Buffer<uint> e = getBodyEdges(b);
	const Vector4 p = getEdgePosition4(e[0]);
	SymMatrix4 mA(ZEROSYMMATRIX4);
	Vector4 vB(0,0,0,0);
	for(uint i=0; i<e.size(); i++)
	{
		const Vector4 vi = getTransformed4(getEdgeVector4(e[i]));
		mA += vi.outerProduct();
		if(i > 0) vB += vi * vi.dot(getEdgePosition4(e[i]) - p);
	}
	mA += getBodyVector4(b).dual().outerProduct();
	return p + mA.inverse() * vB;
}
Vector4 Mesh::getBodyPosition(const uint b) const
{
	if(m_dim == 3) return Vector4(getBodyPosition3(b),0);
	return getBodyPosition4(b);
}

Vector4 Mesh::getQuadPosition(const uint q) const
{
	const Buffer<uint> e = getQuadEdges(q);
	const Vector4 p = getEdgePosition4(e[0]);
	SymMatrix4 mA(ZEROSYMMATRIX4);
	Vector4 vB(0,0,0,0);
	for(uint i=0; i<e.size(); i++)
	{
		const Vector4 vi = getTransformed4(getEdgeVector4(e[i]));
		mA += vi.outerProduct();
		if(i > 0) vB += vi * vi.dot(getEdgePosition4(e[i]) - p);
	}
	return p + mA.inverse() * vB;
}

double Mesh::getEdgeAverage1(const uint e) const
{
	const Buffer<uint> &n = getEdgeNodes(e);
	return 0.5 * (getNodePosition1(n[0]) + getNodePosition1(n[1]));
}
Vector2 Mesh::getEdgeAverage2(const uint e) const
{
	const Buffer<uint> &n = getEdgeNodes(e);
	return 0.5 * (getNodePosition2(n[0]) + getNodePosition2(n[1]));
}
Vector3 Mesh::getEdgeAverage3(const uint e) const
{
	const Buffer<uint> &n = getEdgeNodes(e);
	return 0.5 * (getNodePosition3(n[0]) + getNodePosition3(n[1]));
}
Vector4 Mesh::getEdgeAverage4(const uint e) const
{
	const Buffer<uint> &n = getEdgeNodes(e);
	return 0.5 * (getNodePosition4(n[0]) + getNodePosition4(n[1]));
}
Vector4 Mesh::getEdgeAverage(const uint e) const
{
	if(m_dim == 1) return Vector4(getEdgeAverage1(e),0,0,0);
	if(m_dim == 2) return Vector4(getEdgeAverage2(e),0,0);
	if(m_dim == 3) return Vector4(getEdgeAverage3(e),0);
	return getEdgeAverage4(e);
}

Vector2 Mesh::getFaceAverage2(const uint f) const
{
	const Buffer<uint> n = getFaceNodes(f);
	const Vector2 p = getNodePosition2(n[0]);
	if(n.size() == 3) return (p + getNodePosition2(n[1]) + getNodePosition2(n[2])) / 3.0;

	uint i;
	Buffer<TwoVector2> vv(n.size()-2);
	Buffer<Vector2> dd(n.size()-1);
	TwoVector2 v(0);
	dd[0] = getNodePosition2(n[1]) - p;
	for(i=0; i<vv.size(); i++)
	{
		dd[i+1] = getNodePosition2(n[i+2]) - p;
		vv[i] = TwoVector2(dd[i], dd[i+1]);
		v += vv[i];
	}
	v.xy = 1.0 / v.xy;
	Vector2 d(0,0);
	for(i=0; i<vv.size(); i++) d += vv[i].dot(v) * (dd[i] + dd[i+1]);
	return p + d / 3.0;
}
Vector3 Mesh::getFaceAverage3(const uint f) const
{
	const Buffer<uint> n = getFaceNodes(f);
	const Vector3 p = getNodePosition3(n[0]);
	if(n.size() == 3) return (p + getNodePosition3(n[1]) + getNodePosition3(n[2])) / 3.0;

	uint i;
	Buffer<TwoVector3> vv(n.size()-2);
	Buffer<Vector3> dd(n.size()-1);
	TwoVector3 v(0,0,0);
	dd[0] = getNodePosition3(n[1]) - p;
	for(i=0; i<vv.size(); i++)
	{
		dd[i+1] = getNodePosition3(n[i+2]) - p;
		vv[i] = TwoVector3(dd[i], dd[i+1]);
		v += vv[i];
	}
	v /= v.lensq();
	Vector3 d(0,0,0);
	for(i=0; i<vv.size(); i++) d += vv[i].dot(v) * (dd[i] + dd[i+1]);
	return p + d / 3.0;
}
Vector4 Mesh::getFaceAverage4(const uint f) const
{
	const Buffer<uint> n = getFaceNodes(f);
	const Vector4 p = getNodePosition4(n[0]);
	if(n.size() == 3) return (p + getNodePosition4(n[1]) + getNodePosition4(n[2])) / 3.0;

	uint i;
	Buffer<TwoVector4> vv(n.size()-2);
	Buffer<Vector4> dd(n.size()-1);
	TwoVector4 v(0,0,0,0,0,0);
	dd[0] = getNodePosition4(n[1]) - p;
	for(i=0; i<vv.size(); i++)
	{
		dd[i+1] = getNodePosition4(n[i+2]) - p;
		vv[i] = TwoVector4(dd[i], dd[i+1]);
		v += vv[i];
	}
	v /= v.lensq();
	Vector4 d(0,0,0,0);
	for(i=0; i<vv.size(); i++) d += vv[i].dot(v) * (dd[i] + dd[i+1]);
	return p + d / 3.0;
}
Vector4 Mesh::getFaceAverage(const uint f) const
{
	if(m_dim == 2) return Vector4(getFaceAverage2(f),0,0);
	if(m_dim == 3) return Vector4(getFaceAverage3(f),0);
	return getFaceAverage4(f);
}

Vector3 Mesh::getBodyAverage3(const uint b) const
{
	const Buffer<uint> &f = getBodyFaces(b);
	if(f.size() == 4)
	{
		const Buffer<uint> n = getBodyNodes(b);
		return 0.25 * (getNodePosition3(n[0]) + getNodePosition3(n[1]) + getNodePosition3(n[2]) + getNodePosition3(n[3]));
	}

	uint i;
	const uint node = getFaceAnyNode(f.back());
	const Buffer<uint> cf = f.getComplement(getNodeFaces(node));
	const Vector3 p = getNodePosition3(node);
	Buffer<ThreeVector3> vv(cf.size());
	Buffer<Vector3> dd(cf.size());
	ThreeVector3 v(0);
	for(i=0; i<vv.size(); i++)
	{
		dd[i] = getFaceAverage3(cf[i]) - p;
		vv[i] = getBodyIncidence(b, cf[i]) * ThreeVector3(getFaceVector3(cf[i]), dd[i]);
		v += vv[i];
	}
	v.xyz = 1.0 / v.xyz;
	Vector3 d(0,0,0);
	for(i=0; i<vv.size(); i++) d += vv[i].dot(v) * dd[i];
	return p + 0.75 * d;
}
Vector4 Mesh::getBodyAverage4(const uint b) const
{
	const Buffer<uint> &f = getBodyFaces(b);
	if(f.size() == 4)
	{
		const Buffer<uint> n = getBodyNodes(b);
		return 0.25 * (getNodePosition4(n[0]) + getNodePosition4(n[1]) + getNodePosition4(n[2]) + getNodePosition4(n[3]));
	}

	uint i;
	const uint node = getFaceAnyNode(f.back());
	const Buffer<uint> cf = f.getComplement(getNodeFaces(node));
	const Vector4 p = getNodePosition4(node);
	Buffer<ThreeVector4> vv(cf.size());
	Buffer<Vector4> dd(cf.size());
	ThreeVector4 v(0,0,0,0);
	for(i=0; i<vv.size(); i++)
	{
		dd[i] = getFaceAverage4(cf[i]) - p;
		vv[i] = getBodyIncidence(b, cf[i]) * ThreeVector4(getFaceVector4(cf[i]), dd[i]);
		v += vv[i];
	}
	v /= v.lensq();
	Vector4 d(0,0,0,0);
	for(i=0; i<vv.size(); i++) d += vv[i].dot(v) * dd[i];
	return p + 0.75 * d;
}
Vector4 Mesh::getBodyAverage(const uint b) const
{
	if(m_dim == 3) return Vector4(getBodyAverage3(b),0);
	return getBodyAverage4(b);
}

Vector4 Mesh::getQuadAverage(const uint q) const
{
	const Buffer<uint> &b = getQuadBodies(q);
	if(b.size() == 5)
	{
		const Buffer<uint> n = getQuadNodes(q);
		return 0.2 * (getNodePosition4(n[0]) + getNodePosition4(n[1]) + getNodePosition4(n[2]) + getNodePosition4(n[3]) + getNodePosition4(n[4]));
	}

	uint i;
	const uint node = getBodyAnyNode(b.back());
	const Buffer<uint> cb = b.getComplement(getNodeBodies(node));
	const Vector4 p = getNodePosition4(node);
	Buffer<FourVector4> vv(cb.size());
	Buffer<Vector4> dd(cb.size());
	FourVector4 v(0);
	for(i=0; i<vv.size(); i++)
	{
		dd[i] = getBodyAverage4(cb[i]) - p;
		vv[i] = getQuadIncidence(q, cb[i]) * FourVector4(getBodyVector4(cb[i]), dd[i]);
		v += vv[i];
	}
	v.xyzt = 1.0 / v.xyzt;
	Vector4 d(0,0,0,0);
	for(i=0; i<vv.size(); i++) d += vv[i].dot(v) * dd[i];
	return p + 0.8 * d;
}

double Mesh::getNodeDualAverage1(const uint n) const
{
	const Buffer<uint> &e = getNodeEdges(n);
	if(e.empty()) return getNodePosition1(n);
	if(e.size() == 1) return 0.5 * (getNodePosition1(n) + getEdgePosition1(e[0]));
	return 0.5 * (getEdgePosition1(e[1]) + getEdgePosition1(e[0]));
}
Vector2 Mesh::getNodeDualAverage2(const uint n) const
{
	const Buffer<uint> &e = getNodeEdges(n);
	if(e.empty()) return getNodePosition2(n);

	uint i;
	Buffer<TwoVector2> vv(e.size());
	Buffer<Vector2> dd(e.size());
	TwoVector2 v(0);
	const Vector2 p = getNodePosition2(n);
	for(i=0; i<e.size(); i++)
	{
		dd[i] = getEdgeDualAverage2(e[i]) - p;
		vv[i] = getEdgeIncidence(e[i], n) * TwoVector2(getEdgeDualVector2(e[i]), dd[i]);
		v += vv[i];
	}
	v.xy = 1.0 / v.xy;
	Vector2 d(0,0);
	for(i=0; i<e.size(); i++) d += vv[i].dot(v) * dd[i];
	if(m_fsize == 0) return p + 0.5 * d; // assume 1-dimensional mesh
	return p + 2.0 / 3.0 * d; // assume 2-dimensional mesh
}
Vector3 Mesh::getNodeDualAverage3(const uint n) const
{
	const Buffer<uint> &e = getNodeEdges(n);
	if(e.empty()) return getNodePosition3(n);

	uint i;
	Buffer<ThreeVector3> vv(e.size());
	Buffer<Vector3> dd(e.size());
	ThreeVector3 v(0);
	const Vector3 p = getNodePosition3(n);
	for(i=0; i<e.size(); i++)
	{
		dd[i] = getEdgeDualAverage3(e[i]) - p;
		vv[i] = getEdgeIncidence(e[i], n) * ThreeVector3(getEdgeDualVector3(e[i]), dd[i]);
		v += vv[i];
	}
	v.xyz = 1.0 / v.xyz;
	Vector3 d(0,0,0);
	for(i=0; i<e.size(); i++) d += vv[i].dot(v) * dd[i];
	if(m_fsize == 0) return p + 0.5 * d; // assume 1-dimensional mesh
	if(m_bsize == 0) return p + 2.0 / 3.0 * d; // assume 2-dimensional mesh
	return p + 0.75 * d; // assume 3-dimensional mesh
}
Vector4 Mesh::getNodeDualAverage4(const uint n) const
{
	const Buffer<uint> &e = getNodeEdges(n);
	if(e.empty()) return getNodePosition4(n);

	uint i;
	Buffer<FourVector4> vv(e.size());
	Buffer<Vector4> dd(e.size());
	FourVector4 v(0);
	const Vector4 p = getNodePosition4(n);
	for(i=0; i<e.size(); i++)
	{
		dd[i] = getEdgeDualAverage4(e[i]) - p;
		vv[i] = getEdgeIncidence(e[i], n) * FourVector4(getEdgeDualVector4(e[i]), dd[i]);
		v += vv[i];
	}
	v.xyzt = 1.0 / v.xyzt;
	Vector4 d(0,0,0,0);
	for(i=0; i<e.size(); i++) d += vv[i].dot(v) * dd[i];
	if(m_fsize == 0) return p + 0.5 * d; // assume 1-dimensional mesh
	if(m_bsize == 0) return p + 2.0 / 3.0 * d; // assume 2-dimensional mesh
	if(m_qsize == 0) return p + 0.75 * d; // assume 3-dimensional mesh
	return p + 0.8 * d; // assume 4-dimensional mesh
}
Vector4 Mesh::getNodeDualAverage(const uint n) const
{
	if(m_dim == 1) return Vector4(getNodeDualAverage1(n),0,0,0);
	if(m_dim == 2) return Vector4(getNodeDualAverage2(n),0,0);
	if(m_dim == 3) return Vector4(getNodeDualAverage3(n),0);
	return getNodeDualAverage4(n);
}

Vector2 Mesh::getEdgeDualAverage2(const uint e) const
{
	const Buffer<uint> &f = getEdgeFaces(e);
	if(f.empty()) return getEdgePosition2(e);
	if(f.size() == 1) return 0.5 * (getEdgePosition2(e) + getFacePosition2(f[0]));
	return 0.5 * (getFacePosition2(f[1]) + getFacePosition2(f[0]));
}
Vector3 Mesh::getEdgeDualAverage3(const uint e) const
{
	const Buffer<uint> &f = getEdgeFaces(e);
	if(f.empty()) return getEdgePosition3(e);

	uint i;
	Buffer<TwoVector3> vv(f.size());
	Buffer<Vector3> dd(f.size());
	TwoVector3 v(0,0,0);
	const Vector3 p = getEdgePosition3(e);
	for(i=0; i<f.size(); i++)
	{
		dd[i] = getFaceDualAverage3(f[i]) - p;
		vv[i] = getFaceIncidence(f[i], e) * TwoVector3(getFaceDualVector3(f[i]), dd[i]);
		v += vv[i];
	}
	v /= v.lensq();
	Vector3 d(0,0,0);
	for(i=0; i<f.size(); i++) d += vv[i].dot(v) * dd[i];
	if(m_bsize == 0) return p + 0.5 * d; // assume 2-dimensional mesh
	return p + 2.0 / 3.0 * d; // assume 3-dimensional mesh
}
Vector4 Mesh::getEdgeDualAverage4(const uint e) const
{
	const Buffer<uint> &f = getEdgeFaces(e);
	if(f.empty()) return getEdgePosition4(e);

	uint i;
	Buffer<ThreeVector4> vv(f.size());
	Buffer<Vector4> dd(f.size());
	ThreeVector4 v(0,0,0,0);
	const Vector4 p = getEdgePosition4(e);
	for(i=0; i<f.size(); i++)
	{
		dd[i] = getFaceDualAverage4(f[i]) - p;
		vv[i] = getFaceIncidence(f[i], e) * ThreeVector4(getFaceDualVector4(f[i]), dd[i]);
		v += vv[i];
	}
	v /= v.lensq();
	Vector4 d(0,0,0,0);
	for(i=0; i<f.size(); i++) d += vv[i].dot(v) * dd[i];
	if(m_bsize == 0) return p + 0.5 * d; // assume 2-dimensional mesh
	if(m_qsize == 0) return p + 2.0 / 3.0 * d; // assume 3-dimensional mesh
	return p + 0.75 * d; // assume 4-dimensional mesh
}
Vector4 Mesh::getEdgeDualAverage(const uint e) const
{
	if(m_dim == 1) return Vector4(getEdgePosition1(e),0,0,0);
	if(m_dim == 2) return Vector4(getEdgeDualAverage2(e),0,0);
	if(m_dim == 3) return Vector4(getEdgeDualAverage3(e),0);
	return getEdgeDualAverage4(e);
}

Vector3 Mesh::getFaceDualAverage3(const uint f) const
{
	const Buffer<uint> &b = getFaceBodies(f);
	if(b.empty()) return getFacePosition3(f);
	if(b.size() == 1) return 0.5 * (getFacePosition3(f) + getBodyPosition3(b[0]));
	return 0.5 * (getBodyPosition3(b[1]) + getBodyPosition3(b[0]));
}
Vector4 Mesh::getFaceDualAverage4(const uint f) const
{
	const Buffer<uint> &b = getFaceBodies(f);
	if(b.empty()) return getFacePosition4(f);

	uint i;
	Buffer<TwoVector4> vv(b.size());
	Buffer<Vector4> dd(b.size());
	TwoVector4 v(0,0,0,0,0,0);
	const Vector4 p = getFacePosition4(f);
	for(i=0; i<b.size(); i++)
	{
		dd[i] = getBodyDualAverage4(b[i]) - p;
		vv[i] = getBodyIncidence(b[i], f) * TwoVector4(getBodyDualVector4(b[i]), dd[i]);
		v += vv[i];
	}
	v /= v.lensq();
	Vector4 d(0,0,0,0);
	for(i=0; i<b.size(); i++) d += vv[i].dot(v) * dd[i];
	if(m_qsize == 0) return p + 0.5 * d; // assume 3-dimensional mesh
	return p + 2.0 / 3.0 * d; // assume 4-dimensional mesh
}
Vector4 Mesh::getFaceDualAverage(const uint f) const
{
	if(m_dim == 2) return Vector4(getFacePosition2(f),0,0);
	if(m_dim == 3) return Vector4(getFaceDualAverage3(f),0);
	return getFaceDualAverage4(f);
}

Vector4 Mesh::getBodyDualAverage4(const uint b) const
{
	const Buffer<uint> &q = getBodyQuads(b);
	if(q.empty()) return getBodyPosition4(b);
	if(q.size() == 1) return 0.5 * (getBodyPosition4(b) + getQuadPosition(q[0]));
	return 0.5 * (getQuadPosition(q[1]) + getQuadPosition(q[0]));
}
Vector4 Mesh::getBodyDualAverage(const uint b) const
{
	if(m_dim == 3) return Vector4(getBodyPosition3(b),0);
	return getBodyDualAverage4(b);
}

double Mesh::getEdgeVector1(const uint e) const
{
	const Buffer<uint> &n = getEdgeNodes(e);
	return getNodePosition1(n[1]) - getNodePosition1(n[0]);
}
Vector2 Mesh::getEdgeVector2(const uint e) const
{
	const Buffer<uint> &n = getEdgeNodes(e);
	return getNodePosition2(n[1]) - getNodePosition2(n[0]);
}
Vector3 Mesh::getEdgeVector3(const uint e) const
{
	const Buffer<uint> &n = getEdgeNodes(e);
	return getNodePosition3(n[1]) - getNodePosition3(n[0]);
}
Vector4 Mesh::getEdgeVector4(const uint e) const
{
	const Buffer<uint> &n = getEdgeNodes(e);
	return getNodePosition4(n[1]) - getNodePosition4(n[0]);
}
Vector4 Mesh::getEdgeVector(const uint e) const
{
	if(m_dim == 1) return Vector4(getEdgeVector1(e),0,0,0);
	if(m_dim == 2) return Vector4(getEdgeVector2(e),0,0);
	if(m_dim == 3) return Vector4(getEdgeVector3(e),0);
	return getEdgeVector4(e);
}

TwoVector2 Mesh::getFaceVector2(const uint f) const
{
	const Buffer<uint> n = getFaceNodes(f);
	const Vector2 p = getNodePosition2(n[0]);
	if(n.size() == 3) return 0.5 * TwoVector2(getNodePosition2(n[1]) - p, getNodePosition2(n[2]) - p);

	TwoVector2 v(0);
	Vector2 b = getNodePosition2(n[1]) - p;
	for(uint i=2; i<n.size(); i++)
	{
		const Vector2 a = b;
		b = getNodePosition2(n[i]) - p;
		v += TwoVector2(a,b);
	}
	return 0.5 * v;
}
TwoVector3 Mesh::getFaceVector3(const uint f) const
{
	const Buffer<uint> n = getFaceNodes(f);
	const Vector3 p = getNodePosition3(n[0]);
	if(n.size() == 3) return 0.5 * TwoVector3(getNodePosition3(n[1]) - p, getNodePosition3(n[2]) - p);

	TwoVector3 v(0,0,0);
	Vector3 b = getNodePosition3(n[1]) - p;
	for(uint i=2; i<n.size(); i++)
	{
		const Vector3 a = b;
		b = getNodePosition3(n[i]) - p;
		v += TwoVector3(a,b);
	}
	return 0.5 * v;
}
TwoVector4 Mesh::getFaceVector4(const uint f) const
{
	const Buffer<uint> n = getFaceNodes(f);
	const Vector4 p = getNodePosition4(n[0]);
	if(n.size() == 3) return 0.5 * TwoVector4(getNodePosition4(n[1]) - p, getNodePosition4(n[2]) - p);

	TwoVector4 v(0,0,0,0,0,0);
	Vector4 b = getNodePosition4(n[1]) - p;
	for(uint i=2; i<n.size(); i++)
	{
		const Vector4 a = b;
		b = getNodePosition4(n[i]) - p;
		v += TwoVector4(a,b);
	}
	return 0.5 * v;
}
TwoVector4 Mesh::getFaceVector(const uint f) const
{
	if(m_dim == 2) return TwoVector4(getFaceVector2(f),0,0,0,0,0);
	if(m_dim == 3) return TwoVector4(getFaceVector3(f),0,0,0);
	return getFaceVector4(f);
}

ThreeVector3 Mesh::getBodyVector3(const uint b) const
{
	const Buffer<uint> &f = getBodyFaces(b);
	if(f.size() == 4)
	{
		const Buffer<uint> n = getFaceNodes(f[0]);
		return ThreeVector3(getFaceVector3(f[0]), getNodePosition3(n[0]) - getNodePosition3(getBodyNodes(b).getFirstComplement(n, NONE))) / 3.0;
	}

	const uint node = getFaceAnyNode(f.back());
	const Buffer<uint> cf = f.getComplement(getNodeFaces(node));
	const Vector3 p = getNodePosition3(node);
	ThreeVector3 v(0);
	for(uint i=0; i<cf.size(); i++)
	{
		v += getBodyIncidence(b, cf[i]) * ThreeVector3(getFaceVector3(cf[i]), getNodePosition3(getFaceAnyNode(cf[i])) - p);
	}
	return v / 3.0;
}
ThreeVector4 Mesh::getBodyVector4(const uint b) const
{
	const Buffer<uint> &f = getBodyFaces(b);
	if(f.size() == 4)
	{
		const Buffer<uint> n = getFaceNodes(f[0]);
		return ThreeVector4(getFaceVector4(f[0]), getNodePosition4(n[0]) - getNodePosition4(getBodyNodes(b).getFirstComplement(n, NONE))) / 3.0;
	}

	const uint node = getFaceAnyNode(f.back());
	const Buffer<uint> cf = f.getComplement(getNodeFaces(node));
	const Vector4 p = getNodePosition4(node);
	ThreeVector4 v(0,0,0,0);
	for(uint i=0; i<cf.size(); i++)
	{
		v += getBodyIncidence(b, cf[i]) * ThreeVector4(getFaceVector4(cf[i]), getNodePosition4(getFaceAnyNode(cf[i])) - p);
	}
	return v / 3.0;
}
ThreeVector4 Mesh::getBodyVector(const uint b) const
{
	if(m_dim == 3) return ThreeVector4(getBodyVector3(b),0,0,0);
	return getBodyVector4(b);
}

FourVector4 Mesh::getQuadVector(const uint q) const
{
	const Buffer<uint> &b = getQuadBodies(q);
	if(b.size() == 4)
	{
		const Buffer<uint> n = getBodyNodes(b[0]);
		return 0.25 * FourVector4(getBodyVector4(b[0]), getNodePosition4(n[0]) - getNodePosition4(getQuadNodes(q).getFirstComplement(n, NONE)));
	}

	const uint node = getBodyAnyNode(b.back());
	const Buffer<uint> cb = b.getComplement(getNodeBodies(node));
	const Vector4 p = getNodePosition4(node);
	FourVector4 v(0);
	for(uint i=0; i<cb.size(); i++)
	{
		v += getQuadIncidence(q, cb[i]) * FourVector4(getBodyVector4(cb[i]), getNodePosition4(getBodyAnyNode(cb[i])) - p);
	}
	return 0.25 * v;
}

double Mesh::getNodeDualVector1(const uint n) const
{
	if(m_esize == 0) return 1.0;

	double v = 0.0;
	const double p = getNodePosition1(n);
	const Buffer<uint> &e = getNodeEdges(n);
	for(uint i=0; i<e.size(); i++)
	{
		v += getEdgeIncidence(e[i], n) * getEdgeDualVector1(e[i]) * (p - getEdgePosition1(e[i]));
	}
	return v; // assume 1-dimensional mesh
}
TwoVector2 Mesh::getNodeDualVector2(const uint n) const
{
	if(m_esize == 0) return TwoVector2(1.0);

	TwoVector2 v(0.0);
	const Vector2 p = getNodePosition2(n);
	const Buffer<uint> &e = getNodeEdges(n);
	for(uint i=0; i<e.size(); i++)
	{
		v += getEdgeIncidence(e[i], n) * TwoVector2(p - getEdgePosition2(e[i]), getEdgeDualVector2(e[i]));
	}
	if(m_fsize == 0) return v; // assume 1-dimensional mesh
	return 0.5 * v; // assume 2-dimensional mesh
}
ThreeVector3 Mesh::getNodeDualVector3(const uint n) const
{
	if(m_esize == 0) return ThreeVector3(1.0);

	ThreeVector3 v(0.0);
	const Vector3 p = getNodePosition3(n);
	const Buffer<uint> &e = getNodeEdges(n);
	for(uint i=0; i<e.size(); i++)
	{
		v += getEdgeIncidence(e[i], n) * ThreeVector3(p - getEdgePosition3(e[i]), getEdgeDualVector3(e[i]));
	}
	if(m_fsize == 0) return v; // assume 1-dimensional mesh
	if(m_bsize == 0) return 0.5 * v; // assume 2-dimensional mesh
	return v / 3.0; // assume 3-dimensional mesh
}
FourVector4 Mesh::getNodeDualVector4(const uint n) const
{
	if(m_esize == 0) return FourVector4(1.0);

	FourVector4 v(0.0);
	const Vector4 p = getNodePosition4(n);
	const Buffer<uint> &e = getNodeEdges(n);
	for(uint i=0; i<e.size(); i++)
	{
		v += getEdgeIncidence(e[i], n) * FourVector4(p - getEdgePosition4(e[i]), getEdgeDualVector4(e[i]));
	}
	if(m_fsize == 0) return v; // assume 1-dimensional mesh
	if(m_bsize == 0) return 0.5 * v; // assume 2-dimensional mesh
	if(m_qsize == 0) return v / 3.0; // assume 3-dimensional mesh
	return 0.25 * v; // assume 4-dimensional mesh
}
FourVector4 Mesh::getNodeDualVector(const uint n) const
{
	if(m_dim == 1) return FourVector4(getNodeDualVector1(n));
	if(m_dim == 2) return FourVector4(getNodeDualVector2(n).xy);
	if(m_dim == 3) return FourVector4(getNodeDualVector3(n).xyz);
	return getNodeDualVector4(n);
}

double Mesh::getEdgeDualVector1(const uint e) const
{
	return (getEdgeVector1(e) > 0.0 ? 1.0 : -1.0);
}
Vector2 Mesh::getEdgeDualVector2(const uint e) const
{
	if(m_fsize == 0) return getEdgeVector2(e).dual().unit();

	Vector2 v(0,0);
	const Vector2 p = getEdgePosition2(e);
	const Buffer<uint> &f = getEdgeFaces(e);
	for(uint i=0; i<f.size(); i++)
	{
		v += getFaceIncidence(f[i], e) * getFaceDualVector2(f[i]) * (getFacePosition2(f[i]) - p);
	}
	return v; // assume 2-dimensional mesh
}
TwoVector3 Mesh::getEdgeDualVector3(const uint e) const
{
	if(m_fsize == 0) return getEdgeVector3(e).dual().unit();

	TwoVector3 v(0,0,0);
	const Vector3 p = getEdgePosition3(e);
	const Buffer<uint> &f = getEdgeFaces(e);
	for(uint i=0; i<f.size(); i++)
	{
		v += getFaceIncidence(f[i], e) * TwoVector3(getFacePosition3(f[i]) - p, getFaceDualVector3(f[i]));
	}
	if(m_bsize == 0) return v; // assume 2-dimensional mesh
	return 0.5 * v; // assume 3-dimensional mesh
}
ThreeVector4 Mesh::getEdgeDualVector4(const uint e) const
{
	if(m_fsize == 0) return getEdgeVector4(e).dual().unit();

	ThreeVector4 v(0,0,0,0);
	const Vector4 p = getEdgePosition4(e);
	const Buffer<uint> &f = getEdgeFaces(e);
	for(uint i=0; i<f.size(); i++)
	{
		v += getFaceIncidence(f[i], e) * ThreeVector4(getFacePosition4(f[i]) - p, getFaceDualVector4(f[i]));
	}
	if(m_bsize == 0) return v; // assume 2-dimensional mesh
	if(m_qsize == 0) return 0.5 * v; // assume 3-dimensional mesh
	return v / 3.0; // assume 4-dimensional mesh
}
ThreeVector4 Mesh::getEdgeDualVector(const uint e) const
{
	if(m_dim == 1) return ThreeVector4(0,0,0,getEdgeDualVector1(e));
	if(m_dim == 2) return ThreeVector4(0,0,getEdgeDualVector2(e));
	if(m_dim == 3) return ThreeVector4(0,getEdgeDualVector3(e));
	return getEdgeDualVector4(e);
}

double Mesh::getFaceDualVector2(const uint f) const
{
	return (getFaceVector2(f).dual() > 0.0 ? 1.0 : -1.0);
}
Vector3 Mesh::getFaceDualVector3(const uint f) const
{
	if(m_bsize == 0) return getFaceVector3(f).dual().unit();

	Vector3 v(0,0,0);
	const Vector3 p = getFacePosition3(f);
	const Buffer<uint> &b = getFaceBodies(f);
	for(uint i=0; i<b.size(); i++)
	{
		v += getBodyIncidence(b[i], f) * getBodyDualVector3(b[i]) * (p - getBodyPosition3(b[i]));
	}
	return v; // assume 3-dimensional mesh
}
TwoVector4 Mesh::getFaceDualVector4(const uint f) const
{
	if(m_bsize == 0) return getFaceVector4(f).dual().unit();

	TwoVector4 v(0,0,0,0,0,0);
	const Vector4 p = getFacePosition4(f);
	const Buffer<uint> &b = getFaceBodies(f);
	for(uint i=0; i<b.size(); i++)
	{
		v += getBodyIncidence(b[i], f) * TwoVector4(p - getBodyPosition4(b[i]), getBodyDualVector4(b[i]));
	}
	if(m_qsize == 0) return v; // assume 3-dimensional mesh
	return 0.5 * v; // assume 4-dimensional mesh
}
TwoVector4 Mesh::getFaceDualVector(const uint f) const
{
	if(m_dim == 2) return TwoVector4(0,0,0,0,0,getFaceDualVector2(f));
	if(m_dim == 3) return TwoVector4(0,0,0,getFaceDualVector3(f));
	return getFaceDualVector4(f);
}

double Mesh::getBodyDualVector3(const uint b) const
{
	return (getBodyVector3(b).dual() > 0.0 ? 1.0 : -1.0);
}
Vector4 Mesh::getBodyDualVector4(const uint b) const
{
	if(m_qsize == 0) return getBodyVector4(b).dual().unit();

	Vector4 v(0,0,0,0);
	const Vector4 p = getBodyPosition4(b);
	const Buffer<uint> &q = getBodyQuads(b);
	for(uint i=0; i<q.size(); i++)
	{
		v += getQuadIncidence(q[i], b) * getQuadDualVector(q[i]) * (p - getQuadPosition(q[i]));
	}
	return v; // assume 3-dimensional mesh
}
Vector4 Mesh::getBodyDualVector(const uint b) const
{
	if(m_dim == 3) return Vector4(0,0,0,getBodyDualVector3(b));
	return getBodyDualVector4(b);
}

double Mesh::getQuadDualVector(const uint q) const
{
	return (getQuadVector(q).dual() > 0.0 ? 1.0 : -1.0);
}

double Mesh::getNodeHodge(const uint n) const
{
	if(m_dim == 1) return getNodeDualVector1(n);
	if(m_dim == 2) return getNodeDualVector2(n).xy;
	if(m_dim == 3) return getNodeDualVector3(n).xyz;
	return getNodeDualVector4(n).xyzt;
}
double Mesh::getNodeHodge(const uint n, const double &metric) const
{
	if(m_dim == 1) return metric * getNodeDualVector1(n);
	if(m_dim == 2) return metric * getNodeDualVector2(n).xy;
	if(m_dim == 3) return metric * getNodeDualVector3(n).xyz;
	return metric * getNodeDualVector4(n).xyzt;
}

double Mesh::getEdgeHodge(const uint e) const
{
	if(m_dim == 1)
	{
		const double ev = getEdgeVector1(e);
		const double edv = getEdgeDualVector1(e);
		if(m_m.empty()) return edv / ev;
		return getMetric1() * edv / ev;
	}
	if(m_dim == 2)
	{
		const Vector2 ev = getEdgeVector2(e);
		const Vector2 edv = getEdgeDualVector2(e);
		if(m_m.empty()) return TwoVector2(ev, edv).xy / ev.lensq();
		return TwoVector2(getMetric2() * ev, edv).xy / ev.lensq();
	}
	if(m_dim == 3)
	{
		const Vector3 ev = getEdgeVector3(e);
		const TwoVector3 edv = getEdgeDualVector3(e);
		if(m_m.empty()) return ThreeVector3(ev, edv).xyz / ev.lensq();
		return ThreeVector3(getMetric3() * ev, edv).xyz / ev.lensq();
	}
	const Vector4 ev = getEdgeVector4(e);
	const ThreeVector4 edv = getEdgeDualVector4(e);
	if(m_m.empty()) return FourVector4(ev, edv).xyzt / ev.lensq();
	return FourVector4(getMetric4() * ev, edv).xyzt / ev.lensq();
}
double Mesh::getEdgeHodge(const uint e, const SymMatrix4 &metric) const
{
	if(m_dim == 1)
	{
		const double ev = getEdgeVector1(e);
		const double edv = getEdgeDualVector1(e);
		return metric.xx * edv / ev;
	}
	if(m_dim == 2)
	{
		const Vector2 ev = getEdgeVector2(e);
		const Vector2 edv = getEdgeDualVector2(e);
		return TwoVector2(metric.toSymMatrix2() * ev, edv).xy / ev.lensq();
	}
	if(m_dim == 3)
	{
		const Vector3 ev = getEdgeVector3(e);
		const TwoVector3 edv = getEdgeDualVector3(e);
		return ThreeVector3(metric.toSymMatrix3() * ev, edv).xyz / ev.lensq();
	}
	const Vector4 ev = getEdgeVector4(e);
	const ThreeVector4 edv = getEdgeDualVector4(e);
	return FourVector4(metric * ev, edv).xyzt / ev.lensq();
}

double Mesh::getFaceHodge(const uint f) const
{
	if(m_dim == 2)
	{
		const TwoVector2 fv = getFaceVector2(f);
		const double fdv = getFaceDualVector2(f);
		if(m_m.empty()) return fdv / fv.xy;
		return SymTwoMatrix2(getMetric2()).xyxy * fdv / fv.xy;
	}
	if(m_dim == 3)
	{
		const TwoVector3 fv = getFaceVector3(f);
		const Vector3 fdv = getFaceDualVector3(f);
		if(m_m.empty()) return ThreeVector3(fv, fdv).xyz / fv.lensq();
		return ThreeVector3(SymTwoMatrix3(getMetric3()) * fv, fdv).xyz / fv.lensq();
	}
	const TwoVector4 fv = getFaceVector4(f);
	const TwoVector4 fdv = getFaceDualVector4(f);
	if(m_m.empty()) return FourVector4(fv, fdv).xyzt / fv.lensq();
	return FourVector4(SymTwoMatrix4(getMetric4()) * fv, fdv).xyzt / fv.lensq();
}
double Mesh::getFaceHodge(const uint f, const SymTwoMatrix4 &metric) const
{
	if(m_dim == 2)
	{
		const TwoVector2 fv = getFaceVector2(f);
		const double fdv = getFaceDualVector2(f);
		return metric.toSymTwoMatrix2().xyxy * fdv / fv.xy;
	}
	if(m_dim == 3)
	{
		const TwoVector3 fv = getFaceVector3(f);
		const Vector3 fdv = getFaceDualVector3(f);
		return ThreeVector3(metric.toSymTwoMatrix3() * fv, fdv).xyz / fv.lensq();
	}
	const TwoVector4 fv = getFaceVector4(f);
	const TwoVector4 fdv = getFaceDualVector4(f);
	return FourVector4(metric * fv, fdv).xyzt / fv.lensq();
}

double Mesh::getBodyHodge(const uint b) const
{
	if(m_dim == 3)
	{
		const ThreeVector3 bv = getBodyVector3(b);
		const double bdv = getBodyDualVector3(b);
		if(m_m.empty()) return bdv / bv.xyz;
		return SymThreeMatrix3(getMetric3()).xyzxyz * bdv / bv.xyz;
	}
	const ThreeVector4 bv = getBodyVector4(b);
	const Vector4 bdv = getBodyDualVector4(b);
	if(m_m.empty()) return FourVector4(bv, bdv).xyzt / bv.lensq();
	return FourVector4(SymThreeMatrix4(getMetric4()) * bv, bdv).xyzt / bv.lensq();
}
double Mesh::getBodyHodge(const uint b, const SymThreeMatrix4 &metric) const
{
	if(m_dim == 3)
	{
		const ThreeVector3 bv = getBodyVector3(b);
		const double bdv = getBodyDualVector3(b);
		return metric.xyzxyz * bdv / bv.xyz;
	}
	const ThreeVector4 bv = getBodyVector4(b);
	const Vector4 bdv = getBodyDualVector4(b);
	return FourVector4(metric * bv, bdv).xyzt / bv.lensq();
}

double Mesh::getQuadHodge(const uint q) const
{
	const FourVector4 qv = getQuadVector(q);
	const double qdv = getQuadDualVector(q);
	if(m_m.empty()) return qdv / qv.xyzt;
	return SymFourMatrix4(getMetric4()).xyztxyzt * qdv / qv.xyzt;
}
double Mesh::getQuadHodge(const uint q, const SymFourMatrix4 &metric) const
{
	const FourVector4 qv = getQuadVector(q);
	const double qdv = getQuadDualVector(q);
	return metric.xyztxyzt * qdv / qv.xyzt;
}

Vector4 Mesh::getEdgeDeviation(const uint e, const Vector4 &p) const
{
	if(m_dim == 1) return Vector4(0,0,0,0);
	const Buffer<uint> &n = getEdgeNodes(e);
	const Vector4 p0 = getNodePosition(n[0]);
	const Vector4 v0 = getNodePosition(n[1]) - p0;
	const Vector4 d = p - p0;
	return d - d.dot(v0) / v0.lensq() * v0;
}

Vector4 Mesh::getFaceDeviation(const uint f, const Vector4 &p) const
{
	if(m_dim <= 2) return Vector4(0,0,0,0);

	uint i;
	const Buffer<uint> n = getFaceNodes(f);
	const Vector4 p0 = getNodePosition(n[0]);
	Buffer<Vector4> v(n.size() - 1);
	for(i=0; i<v.size(); i++) v[i] = getNodePosition(n[i + 1]) - p0;
	const double sq0 = v[0].lensq();

	v[1] -= v[1].dot(v[0]) / sq0 * v[0];
	double sq1 = v[1].lensq();
	for(i=2; i<v.size(); i++)
	{
		v[i] -= v[i].dot(v[0]) / sq0 * v[0];
		const double sqi = v[i].lensq();
		if(sqi > sq1)
		{
			sq1 = sqi;
			v[1] = v[i];
		}
	}

	const Vector4 d = p - p0;
	return d - d.dot(v[0]) / sq0 * v[0] - d.dot(v[1]) / sq1 * v[1];
}

Vector4 Mesh::getBodyDeviation(const uint b, const Vector4 &p) const
{
	if(m_dim <= 3) return Vector4(0,0,0,0);

	uint i;
	const Buffer<uint> n = getBodyNodes(b);
	const Vector4 p0 = getNodePosition(n[0]);
	Buffer<Vector4> v(n.size() - 1);
	for(i=0; i<v.size(); i++) v[i] = getNodePosition(n[i + 1]) - p0;
	const double sq0 = v[0].lensq();

	v[1] -= v[1].dot(v[0]) / sq0 * v[0];
	double sq1 = v[1].lensq();
	for(i=2; i<v.size(); i++)
	{
		v[i] -= v[i].dot(v[0]) / sq0 * v[0];
		const double sqi = v[i].lensq();
		if(sqi > sq1)
		{
			sq1 = sqi;
			const Vector4 v1 = v[1];
			v[1] = v[i];
			v[i] = v1;
		}
	}

	v[2] -= v[2].dot(v[1]) / sq1 * v[1];
	double sq2 = v[2].lensq();
	for(i=3; i<v.size(); i++)
	{
		v[i] -= v[i].dot(v[1]) / sq1 * v[1];
		const double sqi = v[i].lensq();
		if(sqi > sq2)
		{
			sq2 = sqi;
			v[2] = v[i];
		}
	}

	const Vector4 d = p - p0;
	return d - d.dot(v[0]) / sq0 * v[0] - d.dot(v[1]) / sq1 * v[1] - d.dot(v[2]) / sq2 * v[2];
}


Vector4 Mesh::getEdgeOrthogonal(const uint e, const Vector4 &d) const
{
	if(m_dim == 1) return Vector4(0,0,0,0);
	const Buffer<uint> &n = getEdgeNodes(e);
	const Vector4 p0 = getNodePosition(n[0]);
	const Vector4 v0 = getNodePosition(n[1]) - p0;
	const Vector4 t0 = getTransformed(v0);
	return d - d.dot(t0) / v0.dot(t0) * v0;
}

Vector4 Mesh::getFaceOrthogonal(const uint f, const Vector4 &d) const
{
	if(m_dim <= 2) return Vector4(0,0,0,0);

	uint i;
	const Buffer<uint> n = getFaceNodes(f);
	const Vector4 p0 = getNodePosition(n[0]);
	Buffer<Vector4> v(n.size() - 1);
	for(i=0; i<v.size(); i++) v[i] = getNodePosition(n[i + 1]) - p0;
	const Vector4 t0 = getTransformed(v[0]);
	const double sq0 = v[0].dot(t0);

	v[1] -= v[1].dot(t0) / sq0 * v[0];
	Vector4 t1 = getTransformed(v[1]);
	double sq1 = v[1].dot(t1);
	for(i=2; i<v.size(); i++)
	{
		v[i] -= v[i].dot(t0) / sq0 * v[0];
		const Vector4 ti = getTransformed(v[i]);
		const double sqi = v[i].dot(ti);
		if(sqi > sq1)
		{
			t1 = ti;
			sq1 = sqi;
			v[1] = v[i];
		}
	}

	return d - d.dot(t0) / sq0 * v[0] - d.dot(t1) / sq1 * v[1];
}

Vector4 Mesh::getBodyOrthogonal(const uint b, const Vector4 &d) const
{
	if(m_dim <= 3) return Vector4(0,0,0,0);

	uint i;
	const Buffer<uint> n = getBodyNodes(b);
	const Vector4 p0 = getNodePosition(n[0]);
	Buffer<Vector4> v(n.size() - 1);
	for(i=0; i<v.size(); i++) v[i] = getTransformed(getNodePosition(n[i + 1]) - p0);
	const Vector4 t0 = getTransformed(v[0]);
	const double sq0 = v[0].dot(t0);

	v[1] -= v[1].dot(t0) / sq0 * v[0];
	Vector4 t1 = getTransformed(v[1]);
	double sq1 = v[1].dot(t1);
	for(i=2; i<v.size(); i++)
	{
		v[i] -= v[i].dot(t0) / sq0 * v[0];
		const Vector4 ti = getTransformed(v[i]);
		const double sqi = v[i].dot(ti);
		if(sqi > sq1)
		{
			t1 = ti;
			sq1 = sqi;
			const Vector4 v1 = v[1];
			v[1] = v[i];
			v[i] = v1;
		}
	}

	v[2] -= v[2].dot(t1) / sq1 * v[1];
	Vector4 t2 = getTransformed(v[2]);
	double sq2 = v[2].dot(t2);
	for(i=3; i<v.size(); i++)
	{
		v[i] -= v[i].dot(t1) / sq1 * v[1];
		const Vector4 ti = getTransformed(v[i]);
		const double sqi = v[i].dot(ti);
		if(sqi > sq2)
		{
			t2 = ti;
			sq2 = sqi;
			v[2] = v[i];
		}
	}

	return d - d.dot(t0) / sq0 * v[0] - d.dot(t1) / sq1 * v[1] - d.dot(t2) / sq2 * v[2];
}

uint Mesh::findNode(const Vector4 &p, const double zerolensq, uint curr, const bool assured) const
{
	if(m_nsize == 0) return NONE;

	// try searching the node
	if(curr >= m_nsize) curr = 0;

	// travel and find the nearest node
	double currsq = (getNodePosition(curr) - p).lensq();
	while(true)
	{
		const uint prev = curr;
		const Buffer<uint> &e = getNodeEdges(prev);
		for(uint i=0; i<e.size(); i++)
		{
			const uint next = getEdgeOtherNode(e[i], prev);
			const double nextsq = (getNodePosition(next) - p).lensq();
			if(nextsq < currsq)
			{
				curr = next;
				currsq = nextsq;
			}
		}
		if(prev == curr) break;
	}
	if(currsq < zerolensq) return curr;

	// not found by search -> check all nodes (this can be computationally heavy)
	if(!assured) return NONE;
	for(uint i=0; i<m_nsize; i++)
	{
		if((getNodePosition(i) - p).lensq() < zerolensq) return i;
	}
	return NONE;
}

uint Mesh::findEdge(const uint n0, const uint n1) const
{
	const Buffer<uint> &e = getNodeEdges(n0);
	for(uint i=0; i<e.size(); i++)
	{
		const Buffer<uint> &n = getEdgeNodes(e[i]);
		if(n[0] == n0 && n[1] == n1) return e[i];
		if(n[1] == n0 && n[0] == n1) return e[i];
	}
	return NONE;
}

uint Mesh::findFace(const Buffer<uint> &e) const
{
	const Buffer<uint> &f = getEdgeFaces(e[0]);
	for(uint i=0; i<f.size(); i++)
	{
		if(e.isAnagram(getFaceEdges(f[i]))) return f[i];
	}
	return NONE;
}

uint Mesh::findBody(const Buffer<uint> &f) const
{
	const Buffer<uint> &b = getFaceBodies(f[0]);
	for(uint i=0; i<b.size(); i++)
	{
		if(f.isAnagram(getBodyFaces(b[i]))) return b[i];
	}
	return NONE;
}

uint Mesh::findQuad(const Buffer<uint> &b) const
{
	const Buffer<uint> &q = getBodyQuads(b[0]);
	for(uint i=0; i<q.size(); i++)
	{
		if(b.isAnagram(getQuadBodies(q[i]))) return q[i];
	}
	return NONE;
}

void Mesh::setNodePosition(const uint n, const Vector4 &p)
{
	const uint nn = m_dim * n;
	m_p[nn] = p.x;
	if(m_dim == 1) return;
	m_p[nn+1] = p.y;
	if(m_dim == 2) return;
	m_p[nn+2] = p.z;
	if(m_dim == 3) return;
	m_p[nn+3] = p.t;
}

void Mesh::setNodeWeight(const uint n, const double w)
{
	if(n >= m_w.size())
	{
		if(w == 0.0) return;
		uint i = m_w.size();
		m_w.resize(m_nsize);
		while(i<m_nsize) m_w[i++] = 0.0;
	}
	m_w[n] = w;
}

void Mesh::setMetric(const SymMatrix4 &m)
{
	m_m.clear();
	if(m_dim == 1)
	{
		if(m.xx == 1.0) return;
		m_m.resize(1);
		m_m[0] = m.xx;
		return;
	}
	if(m_dim == 2)
	{
		if(m.toSymMatrix2() == IDENTITYSYMMATRIX2) return;
		m_m.resize(3);
		m_m[0] = m.xx;
		m_m[1] = m.xy; m_m[2] = m.yy;
		return;
	}
	if(m_dim == 3)
	{
		if(m.toSymMatrix3() == IDENTITYSYMMATRIX3) return;
		m_m.resize(6);
		m_m[0] = m.xx;
		m_m[1] = m.xy; m_m[2] = m.yy;
		m_m[3] = m.xz; m_m[4] = m.yz; m_m[5] = m.zz;
		return;
	}
	if(m == IDENTITYSYMMATRIX4) return;
	m_m.resize(10);
	m_m[0] = m.xx;
	m_m[1] = m.xy; m_m[2] = m.yy;
	m_m[3] = m.xz; m_m[4] = m.yz; m_m[5] = m.zz;
	m_m[6] = m.xt; m_m[7] = m.yt; m_m[8] = m.zt; m_m[9] = m.tt;
}

SymMatrix4 Mesh::getMetric() const
{
	if(m_m.empty()) return IDENTITYSYMMATRIX4;
	if(m_dim == 1) return SymMatrix4(getMetric1(),0,1,0,0,1,0,0,0,1);
	if(m_dim == 2) return SymMatrix4(getMetric2(),0,0,1,0,0,0,1);
	if(m_dim == 3) return SymMatrix4(getMetric3(),0,0,0,1);
	return getMetric4();
}

Vector4 Mesh::getTransformed(const Vector4 &r) const
{
	if(m_m.empty()) return r;
	if(m_dim == 1) return Vector4(getMetric1() * r.x,0,0,0);
	if(m_dim == 2) return Vector4(getMetric2() * r.toVector2(),0,0);
	if(m_dim == 3) return Vector4(getMetric3() * r.toVector3(),0);
	return getMetric4() * r;
}

void Mesh::transform(const Matrix4 &mat)
{
	for(uint i=0; i<m_nsize; i++)
	{
		setNodePosition(i, mat * getNodePosition(i));
	}
	setMetric(getMetric() * (mat * mat.transpose()).toSymMatrix4().inverse());
}

void Mesh::move(const Vector4 &vec)
{
	for(uint i=0; i<m_nsize; i++)
	{
		setNodePosition(i, vec + getNodePosition(i));
	}
}

uint Mesh::addNode(const Vector4 &p)
{
	// check if m_n is full -> resize the table
	const uint res = m_nsize++;
	if(m_nsize > m_n.size()) resizeNodeBuffer(2 * m_nsize);

	// create new node
	setNodePosition(res, p);
	return res;
}
uint Mesh::addEdge(const uint n0, const uint n1)
{
	// check if edge already exists
	uint res = findEdge(n0, n1);
	if(res != NONE) return res;

	// check if m_e is full -> resize the table
	res = m_esize++;
	if(m_esize > m_e.size()) resizeEdgeBuffer(2 * m_esize);

	// create new edge
	m_e[res].n.resize(2);
	m_e[res].n[0] = n0;
	m_e[res].n[1] = n1;
	m_n[n0].e.push_back(res);
	m_n[n1].e.push_back(res);
	return res;
}

uint Mesh::addFace(const Buffer<uint> &e)
{
	// check if face already exists
	uint res = findFace(e);
	if(res != NONE)	return res;

	// check if m_f is full -> resize the table
	res = m_fsize++;
	if(m_fsize > m_f.size()) resizeFaceBuffer(2 * m_fsize);

	// create new face
	m_f[res].e = e;
	for(uint i=0; i<e.size(); i++) m_e[e[i]].f.push_back(res);
	orderFaceEdges(res);
	return res;
}

uint Mesh::addBody(const Buffer<uint> &f)
{
	// check if body already exists
	uint res = findBody(f);
	if(res != NONE)	return res;

	// check if m_b is full -> resize the table
	res = m_bsize++;
	if(m_bsize > m_b.size()) resizeBodyBuffer(2 * m_bsize);

	// create new body
	m_b[res].f = f;
	for(uint i=0; i<f.size(); i++) m_f[f[i]].b.push_back(res);
	orderBodyFaces(res);
	return res;
}

uint Mesh::addQuad(const Buffer<uint> &b)
{
	// check if quad already exists
	uint res = findQuad(b);
	if(res != NONE)	return res;

	// check if m_q is full -> resize the table
	res = m_qsize++;
	if(m_qsize > m_q.size()) resizeQuadBuffer(2 * m_qsize);

	// create new quad
	m_q[res].b = b;
	for(uint i=0; i<b.size(); i++) m_b[b[i]].q.push_back(res);
	orderQuadBodies(res);
	return res;
}

void Mesh::removeNode(const uint n)
{
	uint i;

	// remove linked edges
	Buffer<uint> &e = m_n[n].e;
	for(i=e.size(); i-->0; ) removeEdge(e[i]);

	// replace this node by the last one
	--m_nsize;
	if(n != m_nsize)
	{
		setNodePosition(n, getNodePosition(m_nsize));
		setNodeFlag(n, getNodeFlag(m_nsize));
		setNodeWeight(n, getNodeWeight(m_nsize));

		Buffer<uint> &e = m_n[n].e;
		e.swap(m_n[m_nsize].e);
		for(i=e.size(); i-->0; ) m_e[e[i]].n.replaceFirst(m_nsize, n);
	}
	setNodeFlag(m_nsize, 0);
	setNodeWeight(m_nsize, 0.0);
}

void Mesh::removeEdge(const uint e)
{
	uint i;

	// remove linked faces
	Buffer<uint> &f = m_e[e].f;
	for(i=f.size(); i-->0; ) removeFace(f[i]);

	// remove links from nodes
	Buffer<uint> &n = m_e[e].n;
	for(i=n.size(); i-->0; ) m_n[n[i]].e.eraseFirst(e);
	n.clear();

	// replace this edge by the last one
	--m_esize;
	if(e != m_esize)
	{
		setEdgeFlag(e, getEdgeFlag(m_esize));

		Buffer<uint> &n = m_e[e].n;
		n.swap(m_e[m_esize].n);
		for(i=n.size(); i-->0; ) m_n[n[i]].e.replaceFirst(m_esize, e);

		Buffer<uint> &f = m_e[e].f;
		f.swap(m_e[m_esize].f);
		for(i=f.size(); i-->0; ) m_f[f[i]].e.replaceFirst(m_esize, e);
	}
	setEdgeFlag(m_esize, 0);
}

void Mesh::removeFace(const uint f)
{
	uint i;

	// remove linked bodies
	Buffer<uint> &b = m_f[f].b;
	for(i=b.size(); i-->0; ) removeBody(b[i]);

	// remove links from edges
	Buffer<uint> &e = m_f[f].e;
	for(i=e.size(); i-->0; ) m_e[e[i]].f.eraseFirst(f);
	e.clear();

	// replace this face by the last one
	--m_fsize;
	if(f != m_fsize)
	{
		setFaceFlag(f, getFaceFlag(m_fsize));

		Buffer<uint> &e = m_f[f].e;
		e.swap(m_f[m_fsize].e);
		for(i=e.size(); i-->0; ) m_e[e[i]].f.replaceFirst(m_fsize, f);

		Buffer<uint> &b = m_f[f].b;
		b.swap(m_f[m_fsize].b);
		for(i=b.size(); i-->0; ) m_b[b[i]].f.replaceFirst(m_fsize, f);
	}
	setFaceFlag(m_fsize, 0);
}

void Mesh::removeBody(const uint b)
{
	uint i;

	// remove linked quads
	Buffer<uint> &q = m_b[b].q;
	for(i=q.size(); i-->0; ) removeQuad(q[i]);

	// remove links from faces
	Buffer<uint> &f = m_b[b].f;
	for(i=f.size(); i-->0; ) m_f[f[i]].b.eraseFirst(b);
	f.clear();

	// replace this body by the last one
	--m_bsize;
	if(b != m_bsize)
	{
		setBodyFlag(b, getBodyFlag(m_bsize));

		Buffer<uint> &f = m_b[b].f;
		f.swap(m_b[m_bsize].f);
		for(i=f.size(); i-->0; ) m_f[f[i]].b.replaceFirst(m_bsize, b);

		Buffer<uint> &q = m_b[b].q;
		q.swap(m_b[m_bsize].q);
		for(i=q.size(); i-->0; ) m_q[q[i]].b.replaceFirst(m_bsize, b);
	}
	setBodyFlag(m_bsize, 0);
}

void Mesh::removeQuad(const uint q)
{
	uint i;

	// remove links from bodies
	Buffer<uint> &b = m_q[q].b;
	for(i=b.size(); i-->0; ) m_b[b[i]].q.eraseFirst(q);
	b.clear();

	// replace this quad by the last one
	--m_qsize;
	if(q != m_qsize)
	{
		setQuadFlag(q, getQuadFlag(m_qsize));

		Buffer<uint> &b = m_q[q].b;
		b.swap(m_q[m_qsize].b);
		for(i=b.size(); i-->0; ) m_b[b[i]].q.replaceFirst(m_qsize, q);
	}
	setQuadFlag(m_qsize, 0);
}

void Mesh::setNodeFlag(const uint n, const uint flag)
{
	if(n >= m_nflag.size())
	{
		if(flag == 0) return;
		uint i = m_nflag.size();
		m_nflag.resize(m_n.size());
		while(i<m_n.size()) m_nflag[i++] = 0;
	}
	m_nflag[n] = flag;
}

void Mesh::setEdgeFlag(const uint e, const uint flag)
{
	if(e >= m_eflag.size())
	{
		if(flag == 0) return;
		uint i = m_eflag.size();
		m_eflag.resize(m_e.size());
		while(i<m_e.size()) m_eflag[i++] = 0;
	}
	m_eflag[e] = flag;
}

void Mesh::setFaceFlag(const uint f, const uint flag)
{
	if(f >= m_fflag.size())
	{
		if(flag == 0) return;
		uint i = m_fflag.size();
		m_fflag.resize(m_f.size());
		while(i<m_f.size()) m_fflag[i++] = 0;
	}
	m_fflag[f] = flag;
}

void Mesh::setBodyFlag(const uint b, const uint flag)
{
	if(b >= m_bflag.size())
	{
		if(flag == 0) return;
		uint i = m_bflag.size();
		m_bflag.resize(m_b.size());
		while(i<m_b.size()) m_bflag[i++] = 0;
	}
	m_bflag[b] = flag;
}

void Mesh::setQuadFlag(const uint q, const uint flag)
{
	if(q >= m_qflag.size())
	{
		if(flag == 0) return;
		uint i = m_qflag.size();
		m_qflag.resize(m_q.size());
		while(i<m_q.size()) m_qflag[i++] = 0;
	}
	m_qflag[q] = flag;
}

void Mesh::resizeNodeBuffer(const uint size)
{
	if(size == m_n.size()) return;
	m_p.resize(m_dim * size);
	const uint minsize = size < m_n.size() ? size : m_n.size();
	if(minsize == 0)
	{
		m_n.resize(size);
		return;
	}

	// resize buffer keeping the old sub-buffers
	uint i;
	Buffer< Buffer<uint> > e(minsize);
	for(i=0; i<minsize; i++) m_n[i].e.swap(e[i]);
	m_n.resize(size);
	for(i=0; i<minsize; i++) m_n[i].e.swap(e[i]);
}

void Mesh::resizeEdgeBuffer(const uint size)
{
	if(size == m_e.size()) return;
	const uint minsize = size < m_e.size() ? size : m_e.size();
	if(minsize == 0)
	{
		m_e.resize(size);
		return;
	}

	// resize buffer keeping the old sub-buffers
	uint i;
	Buffer< Buffer<uint> > n(minsize);
	Buffer< Buffer<uint> > f(minsize);
	for(i=0; i<minsize; i++) { m_e[i].n.swap(n[i]); m_e[i].f.swap(f[i]); }
	m_e.resize(size);
	for(i=0; i<minsize; i++) { m_e[i].n.swap(n[i]); m_e[i].f.swap(f[i]); }
}

void Mesh::resizeFaceBuffer(const uint size)
{
	if(size == m_f.size()) return;
	const uint minsize = size < m_f.size() ? size : m_f.size();
	if(minsize == 0)
	{
		m_f.resize(size);
		return;
	}

	// resize buffer keeping the old sub-buffers
	uint i;
	Buffer< Buffer<uint> > e(minsize);
	Buffer< Buffer<uint> > b(minsize);
	for(i=0; i<minsize; i++) { m_f[i].e.swap(e[i]); m_f[i].b.swap(b[i]); }
	m_f.resize(size);
	for(i=0; i<minsize; i++) { m_f[i].e.swap(e[i]); m_f[i].b.swap(b[i]); }
}

void Mesh::resizeBodyBuffer(const uint size)
{
	if(size == m_b.size()) return;
	const uint minsize = size < m_b.size() ? size : m_b.size();
	if(minsize == 0)
	{
		m_b.resize(size);
		return;
	}

	// resize buffer keeping the old sub-buffers
	uint i;
	Buffer< Buffer<uint> > f(minsize);
	Buffer< Buffer<uint> > q(minsize);
	for(i=0; i<minsize; i++) { m_b[i].f.swap(f[i]); m_b[i].q.swap(q[i]); }
	m_b.resize(size);
	for(i=0; i<minsize; i++) { m_b[i].f.swap(f[i]); m_b[i].q.swap(q[i]); }
}

void Mesh::resizeQuadBuffer(const uint size)
{
	if(size == m_q.size()) return;
	const uint minsize = size < m_q.size() ? size : m_q.size();
	if(minsize == 0)
	{
		m_q.resize(size);
		return;
	}

	// resize buffer keeping the old sub-buffers
	uint i;
	Buffer< Buffer<uint> > b(minsize);
	for(i=0; i<minsize; i++) m_q[i].b.swap(b[i]);
	m_q.resize(size);
	for(i=0; i<minsize; i++) m_q[i].b.swap(b[i]);
}
/*
Buffer<double> Mesh::getNodeQuadrature(const uint n) const
{
	uint qdrs = 0;
	Buffer<double> qdr;
	insertQuadrature(getNodePosition(n), 1.0, qdr, qdrs);
	qdr.resize(qdrs);
	return qdr;
}
Buffer<double> Mesh::getEdgeQuadrature(const uint e, const double h) const
{
	const Buffer<uint> &n = getEdgeNodes(e);
	const Vector4 v = getEdgeVector(e);

	uint qdrs = 0;
	Buffer<double> qdr;
	insertQuadrature(getNodePosition(n[0]), getNodePosition(n[1]), v / v.lensq(), h, qdr, qdrs);
	qdr.resize(qdrs);
	return qdr;
}
Buffer<double> Mesh::getFaceQuadrature(const uint f, const double h) const
{
	uint ps = 0;
	Buffer<Vector4> p;
	gatherFaceSimplices(f, p, ps);

	uint i;
	TwoVector4 v(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	for(i=0; i<ps; i+=3) v += TwoVector4(p[i+1] - p[i], p[i+2] - p[i]);
	v *= 2.0 / v.lensq();

	uint qdrs = 0;
	Buffer<double> qdr;
	for(i=0; i<ps; i+=3) insertQuadrature(p[i], p[i+1], p[i+2], v, h, qdr, qdrs);
	qdr.resize(qdrs);
	return qdr;
}
Buffer<double> Mesh::getBodyQuadrature(const uint b, const double h) const
{
	uint ps = 0;
	Buffer<Vector4> p;
	gatherBodySimplices(b, p, ps);

	uint i;
	ThreeVector4 v(0.0, 0.0, 0.0, 0.0);
	for(i=0; i<ps; i+=4) v += ThreeVector4(p[i+1] - p[i], p[i+2] - p[i], p[i+3] - p[i]);
	v *= 6.0 / v.lensq();

	uint qdrs = 0;
	Buffer<double> qdr;
	for(i=0; i<ps; i+=4) insertQuadrature(p[i], p[i+1], p[i+2], p[i+3], v, h, qdr, qdrs);
	qdr.resize(qdrs);
	return qdr;
}
Buffer<double> Mesh::getQuadQuadrature(const uint q, const double h) const
{
	uint ps = 0;
	Buffer<Vector4> p;
	gatherQuadSimplices(q, p, ps);

	uint i;
	FourVector4 v(0.0);
	for(i=0; i<ps; i+=5) v += FourVector4(p[i+1] - p[i], p[i+2] - p[i], p[i+3] - p[i], p[i+4] - p[i]);
	v *= 24.0 / v.lensq();

	uint qdrs = 0;
	Buffer<double> qdr;
	for(i=0; i<ps; i+=5) insertQuadrature(p[i], p[i+1], p[i+2], p[i+3], p[i+4], v, h, qdr, qdrs);
	qdr.resize(qdrs);
	return qdr;
}
Buffer<double> Mesh::getNodeDualQuadrature(const uint n, const double h) const
{
	uint qdrs = 0;
	Buffer<double> qdr;
	if(m_esize == 0) insertQuadrature(getNodePosition(n), 1.0, qdr, qdrs);
	else
	{
		uint i;
		uint ps = 0;
		Buffer<Vector4> p;
		uint vs = 0;
		Buffer<FourVector4> v;
		gatherNodeDualSimplices(n, p, ps, v, vs);

		Buffer<Vector4>



		if(m_fsize == 0)
		{
			Vector4 v(0.0, 0.0, 0.0, 0.0);
			for(i=0; i<ps; i+=2) v += p[i+1] - p[i];
			v *= 1.0 / v.lensq();
			for(i=0; i<ps; i+=2) insertQuadrature(p[i], p[i+1], v, h, qdr, qdrs);
		}
		else if(m_bsize == 0)
		{
			TwoVector4 v(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
			for(i=0; i<ps; i+=3) v += TwoVector4(p[i+1] - p[i], p[i+2] - p[i]);
			v *= 2.0 / v.lensq();
			for(i=0; i<ps; i+=3) insertQuadrature(p[i], p[i+1], p[i+2], v, h, qdr, qdrs);
		}
		else if(m_qsize == 0)
		{
			ThreeVector4 v(0.0, 0.0, 0.0, 0.0);
			for(i=0; i<ps; i+=4) v += ThreeVector4(p[i+1] - p[i], p[i+2] - p[i], p[i+3] - p[i]);
			v *= 6.0 / v.lensq();
			for(i=0; i<ps; i+=4) insertQuadrature(p[i], p[i+1], p[i+2], p[i+3], v, h, qdr, qdrs);
		}
		else
		{
			FourVector4 v(0.0);
			for(i=0; i<ps; i+=5) v += FourVector4(p[i+1] - p[i], p[i+2] - p[i], p[i+3] - p[i], p[i+4] - p[i]);
			v *= 24.0 / v.lensq();
			for(i=0; i<ps; i+=5) insertQuadrature(p[i], p[i+1], p[i+2], p[i+3], p[i+4], v, h, qdr, qdrs);
		}
	}
	qdr.resize(qdrs);
	return qdr;
}
Buffer<double> Mesh::getEdgeDualQuadrature(const uint e, const double h) const
{
	uint qdrs = 0;
	Buffer<double> qdr;
	if(m_fsize == 0) insertQuadrature(getEdgePosition(e), 1.0, qdr, qdrs);
	else
	{
		uint i;
		uint ps = 0;
		Buffer<Vector4> p;
		gatherEdgeDualSimplices(e, p, ps);

		if(m_bsize == 0)
		{
			Vector4 v(0.0, 0.0, 0.0, 0.0);
			for(i=0; i<ps; i+=2) v += p[i+1] - p[i];
			v *= 1.0 / v.lensq();
			for(i=0; i<ps; i+=2) insertQuadrature(p[i], p[i+1], v, h, qdr, qdrs);
		}
		else if(m_qsize == 0)
		{
			TwoVector4 v(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
			for(i=0; i<ps; i+=3) v += TwoVector4(p[i+1] - p[i], p[i+2] - p[i]);
			v *= 2.0 / v.lensq();
			for(i=0; i<ps; i+=3) insertQuadrature(p[i], p[i+1], p[i+2], v, h, qdr, qdrs);
		}
		else
		{
			ThreeVector4 v(0.0, 0.0, 0.0, 0.0);
			for(i=0; i<ps; i+=4) v += ThreeVector4(p[i+1] - p[i], p[i+2] - p[i], p[i+3] - p[i]);
			v *= 6.0 / v.lensq();
			for(i=0; i<ps; i+=4) insertQuadrature(p[i], p[i+1], p[i+2], p[i+3], v, h, qdr, qdrs);
		}
	}
	qdr.resize(qdrs);
	return qdr;
}
Buffer<double> Mesh::getFaceDualQuadrature(const uint f, const double h) const
{
	uint qdrs = 0;
	Buffer<double> qdr;
	if(m_bsize == 0) insertQuadrature(getFacePosition(f), 1.0, qdr, qdrs);
	else
	{
		uint i;
		uint ps = 0;
		Buffer<Vector4> p;
		gatherFaceDualSimplices(f, p, ps);

		if(m_qsize == 0)
		{
			Vector4 v(0.0, 0.0, 0.0, 0.0);
			for(i=0; i<ps; i+=2) v += p[i+1] - p[i];
			v *= 1.0 / v.lensq();
			for(i=0; i<ps; i+=2) insertQuadrature(p[i], p[i+1], v, h, qdr, qdrs);
		}
		else
		{
			TwoVector4 v(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
			for(i=0; i<ps; i+=3) v += TwoVector4(p[i+1] - p[i], p[i+2] - p[i]);
			v *= 2.0 / v.lensq();
			for(i=0; i<ps; i+=3) insertQuadrature(p[i], p[i+1], p[i+2], v, h, qdr, qdrs);
		}
	}
	qdr.resize(qdrs);
	return qdr;
}
Buffer<double> Mesh::getBodyDualQuadrature(const uint b, const double h) const
{
	uint qdrs = 0;
	Buffer<double> qdr;
	if(m_qsize == 0) insertQuadrature(getBodyPosition(b), 1.0, qdr, qdrs);
	else
	{
		uint i;
		uint ps = 0;
		Buffer<Vector4> p;
		gatherBodyDualSimplices(b, p, ps);

		Vector4 v(0.0, 0.0, 0.0, 0.0);
		for(i=0; i<ps; i+=2) v += p[i+1] - p[i];
		v *= 1.0 / v.lensq();
		for(i=0; i<ps; i+=2) insertQuadrature(p[i], p[i+1], v, h, qdr, qdrs);
	}
	qdr.resize(qdrs);
	return qdr;
}
Buffer<double> Mesh::getQuadDualQuadrature(const uint q) const
{
	uint qdrs = 0;
	Buffer<double> qdr;
	insertQuadrature(getQuadPosition(q), 1.0, qdr, qdrs);
	qdr.resize(qdrs);
	return qdr;
}
*/
void Mesh::gatherNodeSimplices(const uint n, Buffer<Vector4> &p, uint &ps, Buffer<double> &v, uint &vs) const
{
	p.gather(getNodePosition(n), ps);
	v.gather(1.0, vs);
}
void Mesh::gatherEdgeSimplices(const uint e, Buffer<Vector4> &p, uint &ps, Buffer<Vector4> &v, uint &vs) const
{
	const Buffer<uint> &n = getEdgeNodes(e);
	p.gather(getNodePosition(n[0]), ps);
	p.gather(getNodePosition(n[1]), ps);
	v.gather(getEdgeVector(e), vs);
}
void Mesh::gatherFaceSimplices(const uint f, Buffer<Vector4> &p, uint &ps, Buffer<TwoVector4> &v, uint &vs) const
{
	const Buffer<uint> n = getFaceNodes(f);
	const Vector4 p0 = getNodePosition(n[0]);
	Vector4 p1 = getNodePosition(n[1]);
	for(uint i=2; i<n.size(); i++)
	{
		const Vector4 p2 = getNodePosition(n[i]);
		p.gather(p0, ps);
		p.gather(p1, ps);
		p.gather(p2, ps);
		v.gather(TwoVector4(0.5 * (p1 - p0), p2 - p0), vs);
		p1 = p2;
	}
}
void Mesh::gatherBodySimplices(const uint b, Buffer<Vector4> &p, uint &ps, Buffer<ThreeVector4> &v, uint &vs) const
{
	const Buffer<uint> &f = getBodyFaces(b);
	const uint n0 = getFaceAnyNode(f.back());
	const Buffer<uint> cf = f.getComplement(getNodeFaces(n0));
	const Vector4 p0 = getNodePosition(n0);
	Buffer<Vector4> fp(3);
	Buffer<TwoVector4> fv(1);
	for(uint i=0; i<cf.size(); i++)
	{
		uint fps = 0;
		uint fvs = 0;
		gatherFaceSimplices(cf[i], fp, fps, fv, fvs);
		const double inc = getBodyIncidence(b, cf[i]) / 3.0;
		for(uint j=0; j<fvs; j++)
		{
			const uint jj = 3 * j;
			p.gather(p0, ps);
			p.gather(fp[jj], ps);
			p.gather(fp[jj+1], ps);
			p.gather(fp[jj+2], ps);
			v.gather(ThreeVector4(fv[j], inc * (fp[jj] - p0)), vs);
		}
	}
}
void Mesh::gatherQuadSimplices(const uint q, Buffer<Vector4> &p, uint &ps, Buffer<FourVector4> &v, uint &vs) const
{
	const Buffer<uint> &b = getQuadBodies(q);
	const uint n0 = getBodyAnyNode(b.back());
	const Buffer<uint> cb = b.getComplement(getNodeBodies(n0));
	const Vector4 p0 = getNodePosition(n0);
	Buffer<Vector4> bp(4);
	Buffer<ThreeVector4> bv(1);
	for(uint i=0; i<cb.size(); i++)
	{
		uint bps = 0;
		uint bvs = 0;
		gatherBodySimplices(cb[i], bp, bps, bv, bvs);
		const double inc = 0.25 * getQuadIncidence(q, cb[i]);
		for(uint j=0; j<bvs; j++)
		{
			const uint jj = 4 * j;
			p.gather(p0, ps);
			p.gather(bp[jj], ps);
			p.gather(bp[jj+1], ps);
			p.gather(bp[jj+2], ps);
			p.gather(bp[jj+3], ps);
			v.gather(inc * FourVector4(bv[j], bp[jj] - p0), vs);
		}
	}
}

void Mesh::gatherNodeDualSimplices(const uint n, Buffer<Vector4> &p, uint &ps, Buffer<FourVector4> &v, uint &vs) const
{
	const Vector4 p0 = getNodePosition(n);
	if(m_esize == 0)
	{
		p.gather(p0, ps);
		v.gather(FourVector4(1.0), vs);
		return;
	}
	const uint num = (m_bsize == 0 ? (m_fsize == 0 ? 1 : 2) : (m_qsize == 0 ? 3 : 4));
	Buffer<Vector4> ep;
	Buffer<ThreeVector4> ev;
	const Buffer<uint> &e = getNodeEdges(n);
	for(uint i=0; i<e.size(); i++)
	{
		uint eps = 0;
		uint evs = 0;
		gatherEdgeDualSimplices(e[i], ep, eps, ev, evs);
		if(eps == 0) continue;

		const Vector4 vi = getEdgeIncidence(e[i], n) / double(num) * (ep[0] - p0);
		for(uint j=0; j<evs; j++)
		{
			const uint jj = num * j;
			p.gather(p0, ps);
			for(uint k=0; k<num; k++) p.gather(ep[jj+k], ps);
			v.gather(FourVector4(ev[j], vi), vs);
		}
	}
}
void Mesh::gatherEdgeDualSimplices(const uint e, Buffer<Vector4> &p, uint &ps, Buffer<ThreeVector4> &v, uint &vs) const
{
	const Vector4 p0 = getEdgePosition(e);
	if(m_fsize == 0)
	{
		p.gather(p0, ps);
		v.gather(getEdgeVector(e).dual().unit(), vs);
		return;
	}

	const uint num = (m_bsize == 0 ? 1 : (m_qsize == 0 ? 2 : 3));
	Buffer<Vector4> fp;
	Buffer<TwoVector4> fv;
	const Buffer<uint> &f = getEdgeFaces(e);
	for(uint i=0; i<f.size(); i++)
	{
		uint fps = 0;
		uint fvs = 0;
		gatherFaceDualSimplices(f[i], fp, fps, fv, fvs);
		if(fps == 0) continue;

		const Vector4 vi = getFaceIncidence(f[i], e) / double(num) * (fp[0] - p0);
		for(uint j=0; j<fvs; j++)
		{
			const uint jj = num * j;
			p.gather(p0, ps);
			for(uint k=0; k<num; k++) p.gather(fp[jj+k], ps);
			v.gather(ThreeVector4(fv[j], vi), vs);
		}
	}
}
void Mesh::gatherFaceDualSimplices(const uint f, Buffer<Vector4> &p, uint &ps, Buffer<TwoVector4> &v, uint &vs) const
{
	const Vector4 p0 = getFacePosition(f);
	if(m_bsize == 0)
	{
		p.gather(p0, ps);
		v.gather(getFaceVector(f).dual().unit(), vs);
		return;
	}

	const uint num = (m_qsize == 0 ? 1 : 2);
	Buffer<Vector4> bp;
	Buffer<Vector4> bv;
	const Buffer<uint> &b = getFaceBodies(f);
	for(uint i=0; i<b.size(); i++)
	{
		uint bps = 0;
		uint bvs = 0;
		gatherBodyDualSimplices(b[i], bp, bps, bv, bvs);
		if(bps == 0) continue;

		const Vector4 vi = getBodyIncidence(b[i], f) / double(num) * (p0 - bp[0]);
		for(uint j=0; j<bvs; j++)
		{
			const uint jj = num * j;
			p.gather(p0, ps);
			for(uint k=0; k<num; k++) p.gather(bp[jj+k], ps);
			v.gather(TwoVector4(vi, bv[j]), vs);
		}
	}
}
void Mesh::gatherBodyDualSimplices(const uint b, Buffer<Vector4> &p, uint &ps, Buffer<Vector4> &v, uint &vs) const
{
	const Vector4 p0 = getBodyPosition(b);
	if(m_qsize == 0)
	{
		p.gather(p0, ps);
		v.gather(getBodyVector(b).dual().unit(), vs);
		return;
	}

	const Buffer<uint> &q = getBodyQuads(b);
	for(uint i=0; i<q.size(); i++)
	{
		p.gather(p0, ps);
		p.gather(getQuadPosition(q[i]), ps);
		v.gather(getQuadIncidence(q[i], b) * getQuadDualVector(q[i]) * (p0 - getQuadPosition(q[i])), vs);
	}
}
void Mesh::gatherQuadDualSimplices(const uint q, Buffer<Vector4> &p, uint &ps, Buffer<double> &v, uint &vs) const
{
	p.gather(getQuadPosition(q), ps);
	v.gather(getQuadDualVector(q), vs);
}

/*
void Mesh::insertQuadrature(const Vector4 &p0, const double w, Buffer<double> &q, uint &qs) const
{
	q.gather(w, qs);
	q.gather(p0.x, qs);
	if(m_dim > 1) q.gather(p0.y, qs);
	if(m_dim > 2) q.gather(p0.z, qs);
	if(m_dim > 3) q.gather(p0.t, qs);
}
void Mesh::insertQuadrature(const Vector4 &p0, const Vector4 &p1, const Vector4 &w, const double h, Buffer<double> &q, uint &qs) const
{
	const Vector4 v1 = p1 - p0;
	const double vol = v1.len() / h + 1.0;
	const uint n = (vol > 2.0 ? uint(vol) : 2);
	const double div = 1.0 / sqrt((n-1) * (n+1));
	const Vector4 d0 = p0 + (1.0 - (n - 1) * div) / 2.0 * v1;
	const double dw = w.dot(v1) / double(n);
	for(uint i=0; i<n; i++) insertQuadrature(d0 + i * div * v1, dw, q, qs);
}
void Mesh::insertQuadrature(const Vector4 &p0, const Vector4 &p1, const Vector4 &p2, const TwoVector4 &w, const double h, Buffer<double> &q, uint &qs) const
{
	const Vector4 v1 = p1 - p0;
	const Vector4 v2 = p2 - p0;
	const TwoVector4 v12(v1, v2);
	const double vol = pow(v12.lensq(), 0.25) / h;
	const uint n = (vol < 2.0 ? 2 : uint(vol) + 1);
	const double div = 1.0 / sqrt((n-1) * (n+2));
	const Vector4 d0 = p0 + (1.0 - (n - 1) * div) / 3.0 * (v1 + v2);
	const double dw = w.dot(v12) / (double(n) * double(n + 1));
	for(uint i=0; i<n; i++) {
		const Vector4 d1 = d0 + i * div * v1;
		for(uint j=0; j<n-i; j++) insertQuadrature(d1 + j * div * v2, dw, q, qs);
	}
}
void Mesh::insertQuadrature(const Vector4 &p0, const Vector4 &p1, const Vector4 &p2, const Vector4 &p3, const ThreeVector4 &w, const double h, Buffer<double> &q, uint &qs) const
{
	const Vector4 v1 = p1 - p0;
	const Vector4 v2 = p2 - p0;
	const Vector4 v3 = p3 - p0;
	const ThreeVector4 v123(v1, v2, v3);
	const double vol = pow(v123.lensq(), 1.0 / 6.0) / h;
	const uint n = (vol < 2.0 ? 2 : uint(vol) + 1);
	const double div = 1.0 / sqrt((n-1) * (n+3));
	const Vector4 d0 = p0 + (1.0 - (n - 1) * div) / 4.0 * (v1 + v2 + v3);
	const double dw = w.dot(v123) / (double(n) * double(n + 1) * double(n + 2));
	for(uint i=0; i<n; i++) {
		const Vector4 d1 = d0 + i * div * v1;
		for(uint j=0; j<n-i; j++) {
			const Vector4 d2 = d1 + j * div * v2;
			for(uint k=0; k<n-i-j; k++) insertQuadrature(d2 + k * div * v3, dw, q, qs);
		}
	}
}
void Mesh::insertQuadrature(const Vector4 &p0, const Vector4 &p1, const Vector4 &p2, const Vector4 &p3, const Vector4 &p4, const FourVector4 &w, const double h, Buffer<double> &q, uint &qs) const
{
	const Vector4 v1 = p1 - p0;
	const Vector4 v2 = p2 - p0;
	const Vector4 v3 = p3 - p0;
	const Vector4 v4 = p4 - p0;
	const FourVector4 v1234(v1, v2, v3, v4);
	const double vol = pow(v1234.lensq(), 0.125) / h;
	const uint n = (vol < 2.0 ? 2 : uint(vol) + 1);
	const double div = 1.0 / sqrt((n-1) * (n+4));
	const Vector4 d0 = p0 + (1.0 - (n - 1) * div) / 5.0 * (v1 + v2 + v3 + v4);
	const double dw = w.dot(v1234) / (double(n) * double(n + 1) * double(n + 2) * double(n + 3));
	for(uint i=0; i<n; i++) {
		const Vector4 d1 = d0 + i * div * v1;
		for(uint j=0; j<n-i; j++) {
			const Vector4 d2 = d1 + j * div * v2;
			for(uint k=0; k<n-i-j; k++) {
				const Vector4 d3 = d2 + k * div * v3;
				for(uint l=0; l<n-i-j-k; l++) insertQuadrature(d3 + l * div * v4, dw, q, qs);
			}
		}
	}
}
*/
/*
void Mesh::insertQuadrature(const Buffer<Vector4> &p, const uint n, const double w, Buffer<double> &q, uint &qs) const
{
	if(p.size() == 1)
	{
		gatherWeightAndPosition(w, p[0], q, qs);
		return;
	}
	const double div = 1.0 / sqrt((n-1) * (n+p.size()-1));
	const Vector4 v1 = p[1] - p[0];
	if(p.size() == 2)
	{
		const Vector4 p0 = p[0] + (1.0 - (n - 1) * div) * 0.5 * v1;
		const double dw = w / double(n);
		for(uint i=0; i<n; i++) gatherWeightAndPosition(dw, p0 + i * div * v1, q, qs);
		return;
	}
	const Vector4 v2 = p[2] - p[0];
	if(p.size() == 3)
	{
		const Vector4 p0 = p[0] + (1.0 - (n - 1) * div) / 3.0 * (v1 + v2);
		const double dw = w / (double(n) * double(n + 1));
		for(uint i=0; i<n; i++) {
			const Vector4 p1 = p0 + i * div * v1;
			for(uint j=0; j<n-i; j++) gatherWeightAndPosition(dw, p1 + j * div * v2, q, qs);
		}
		return TwoVector4(v1, v2).len() / (double(n) * double(n + 1)) * sum;
	}



	double sum = 0.0;
	const double div = (n > 1 ? 1.0 / sqrt((n-1) * (n+p.size()-1)) : 0.0);
	const Vector4 v1 = p[1] - p[0];
	if(p.size() == 2)
	{
		const Vector4 p0 = p[0] + (1.0 - (n - 1) * div) * 0.5 * v1;
		for(uint i=0; i<n; i++) sum += func(p0 + i * div * v1);
		return v1.len() / double(n) * sum;
	}
	const Vector4 v2 = p[2] - p[0];
	if(p.size() == 3)
	{
		const Vector4 p0 = p[0] + (1.0 - (n - 1) * div) / 3.0 * (v1 + v2);
		for(uint i=0; i<n; i++) {
			const Vector4 vi = p0 + i * div * v1;
			for(uint j=0; j<n-i; j++) sum += func(vi + j * div * v2);
		}
		return TwoVector4(v1, v2).len() / (double(n) * double(n + 1)) * sum;
	}
	const Vector4 v3 = p[3] - p[0];
	if(p.size() == 4)
	{
		const Vector4 p0 = p[0] + (1.0 - (n - 1) * div) * 0.25 * (v1 + v2 + v3);
		for(uint i=0; i<n; i++) {
			const Vector4 vi = p0 + i * div * v1;
			for(uint j=0; j<n-i; j++) {
				const Vector4 vj = vi + j * div * v2;
				for(uint k=0; k<n-i-j; k++) sum += func(vj + k * div * v3);
			}
		}
		return ThreeVector4(v1, v2, v3).len() / (double(n) * double(n + 1) * double(n + 2)) * sum;
	}
	const Vector4 v4 = p[4] - p[0];
	const Vector4 p0 = p[0] + (1.0 - (n - 1) * div) * 0.2 * (v1 + v2 + v3 + v4);
	for(uint i=0; i<n; i++) {
		const Vector4 vi = p0 + i * div * v1;
		for(uint j=0; j<n-i; j++) {
			const Vector4 vj = vi + j * div * v2;
			for(uint k=0; k<n-i-j; k++) {
				const Vector4 vk = vj + k * div * v3;
				for(uint l=0; l<n-i-j-k; l++) sum += func(vk + l * div * v4);
			}
		}
	}
	return FourVector4(v1, v2, v3, v4).len() / (double(n) * double(n + 1) * double(n + 2) * double(n + 3)) * sum;
}
*/
