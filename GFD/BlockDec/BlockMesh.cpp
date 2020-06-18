#include "BlockMesh.hpp"
#include <iostream>
using namespace std;
using namespace gfd;

const uint PARTS = 81;

BlockMesh::BlockMesh()
{
	m_dim = 0;
	m_p = Vector4(0,0,0,0);
	m_d = Vector4(1,1,1,1);
	m_s = Uint4(0,0,0,0);
	m_0 = Uint4(0,0,0,0);

	m_nsize = 0;
	m_nlocs = 0;
	m_esize = 0;
	m_elocs = 0;
	m_fsize = 0;
	m_flocs = 0;
	m_bsize = 0;
	m_blocs = 0;
	m_qsize = 0;
	m_qlocs = 0;
}

void BlockMesh::clear()
{
	m_prev.clear();
	m_next.clear();

	m_p = Vector4(0,0,0,0);
	m_d = Vector4(1,1,1,1);
	m_0 = Uint4(0,0,0,0);
	m_type.clear();

	m_nsize = 0;
	m_nlocs = 0;
	m_nbeg.clear();
	m_esize = 0;
	m_elocs = 0;
	m_ebeg.clear();
	m_fsize = 0;
	m_flocs = 0;
	m_fbeg.clear();
	m_bsize = 0;
	m_blocs = 0;
	m_bbeg.clear();
	m_qsize = 0;
	m_qlocs = 0;
	m_qbeg.clear();

	uint i;
	for(i=0; i<m_part.size(); i++) clearMap(m_part[i]);
	m_part.clear();
	for(i=0; i<fg_num; i++) clearMap(m_poly[i]);
}

void BlockMesh::init(const double p, const double d, const uint sx, const uchar *type)
{
	clear();
	m_dim = 1;
	m_part.resize(2);
	m_p = Vector4(p,0,0,0);
	m_d = Vector4(d,1e10,1e10,1e10);
	m_0 = Uint4(1,0,0,0);
	finalizeInit(type, Uint4(sx, 1, 1, 1));
}
void BlockMesh::init(const Vector2 &p, const Vector2 &d, const uint sx, const uint sy, const uchar *type)
{
	clear();
	m_dim = 2;
	m_part.resize(4);
	m_p = Vector4(p,0,0);
	m_d = Vector4(d,1e10,1e10);
	m_0 = Uint4(1,1,0,0);
	finalizeInit(type, Uint4(sx, sy, 1, 1));
}
void BlockMesh::init(const Vector3 &p, const Vector3 &d, const uint sx, const uint sy, const uint sz, const uchar *type)
{
	clear();
	m_dim = 3;
	m_part.resize(8);
	m_p = Vector4(p,0);
	m_d = Vector4(d,1e10);
	m_0 = Uint4(1,1,1,0);
	finalizeInit(type, Uint4(sx, sy, sz, 1));
}
void BlockMesh::init(const Vector4 &p, const Vector4 &d, const uint sx, const uint sy, const uint sz, const uint st, const uchar *type)
{
	clear();
	m_dim = 4;
	m_part.resize(16);
	m_p = p;
	m_d = d;
	m_0 = Uint4(1,1,1,1);
	finalizeInit(type, Uint4(sx, sy, sz, st));
}

void BlockMesh::finalizeInit(const uchar *type, const Uint4 &s)
{
	// divide the area for ranks
	const uint rank = getMPIrank();
	uint rank0 = 0;
	uint rank1 = getMPIranks();
	Uint4 i0(0,0,0,0);
	Uint4 i1 = s;
	ullong ns = getWeight(i0, i1, type, s, rank0, rank1);
	while(rank0 + 1 < rank1)
	{
		uint dori = i1.x - i0.x; uint iori = 0;
		if(m_dim > 1 && i1.y - i0.y > dori) { dori = i1.y - i0.y; iori = 1; }
		if(m_dim > 2 && i1.z - i0.z > dori) { dori = i1.z - i0.z; iori = 2; }
		if(m_dim > 3 && i1.t - i0.t > dori) { dori = i1.t - i0.t; iori = 3; }

		uint mid;
		ullong nsMid = 0;
		const uint rankMid = (rank0 + rank1) / 2;
		const ullong nsMax = ns * (rankMid - rank0) / (rank1 - rank0);
		if(iori == 0) { // divide in x-direction
			for(mid=i0.x; nsMid<nsMax; mid++) nsMid += getWeight(Uint4(mid,i0.y,i0.z,i0.t), Uint4(mid+1,i1.y,i1.z,i1.t), type, s, rank0, rank1);
			if(rank < rankMid) i1.x = mid;
			else i0.x = mid;
		}
		else if(iori == 1) { // divide in y-direction
			for(mid=i0.y; nsMid<nsMax; mid++) nsMid += getWeight(Uint4(i0.x,mid,i0.z,i0.t), Uint4(i1.x,mid+1,i1.z,i1.t), type, s, rank0, rank1);
			if(rank < rankMid) i1.y = mid;
			else i0.y = mid;
		}
		else if(iori == 2) { // divide in z-direction
			for(mid=i0.z; nsMid<nsMax; mid++) nsMid += getWeight(Uint4(i0.x,i0.y,mid,i0.t), Uint4(i1.x,i1.y,mid+1,i1.t), type, s, rank0, rank1);
			if(rank < rankMid) i1.z = mid;
			else i0.z = mid;
		}
		else { // divide in t-direction
			for(mid=i0.t; nsMid<nsMax; mid++) nsMid += getWeight(Uint4(i0.x,i0.y,i0.z,mid), Uint4(i1.x,i1.y,i1.z,mid+1), type, s, rank0, rank1);
			if(rank < rankMid) i1.t = mid;
			else i0.t = mid;
		}

		if(rank < rankMid) {
			rank1 = rankMid;
			ns = nsMid;
		}
		else {
			rank0 = rankMid;
			ns -= nsMid;
		}
	}
	m_s = m_0 + i1 - i0;
	if(i1.x <= i0.x || i1.y <= i0.y || i1.z <= i0.z || i1.t <= i0.t) {
		cout << "BlockMesh::finalizeInit -> the domain is too narrow for rank " << rank << "." << endl;
	}

	// initialize m_type
	Uint4 i;
	m_type.resize(m_s.x * m_s.y * m_s.z * m_s.t);
	for(i.t=0; i.t<m_s.t; i.t++) {
		const bool zerot = (i.t + i0.t < m_0.t);
		for(i.z=0; i.z<m_s.z; i.z++) {
			const bool zeroz = (zerot || i.z + i0.z < m_0.z);
			for(i.y=0; i.y<m_s.y; i.y++) {
				const bool zeroy = (zeroz || i.y + i0.y < m_0.y);
				for(i.x=0; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					if(zeroy || i.x + i0.x < m_0.x) m_type[ii] = 0;
					else m_type[ii] = type[(((i.t + i0.t - m_0.t) * s.z + i.z + i0.z - m_0.z) *
											s.y + i.y + i0.y - m_0.y) * s.x + i.x + i0.x - m_0.x];
				}
			}
		}
	}
	m_p += Vector4((i0.x-double(m_0.x))*m_d.x, (i0.y-double(m_0.y))*m_d.y, (i0.z-double(m_0.z))*m_d.z, (i0.t-double(m_0.t))*m_d.t);
	synchronizeNeighbors(i0, i1);

	// create parts
	Uint4 j;
	const Uint4 k(1,2,4,8);
	for(i.t=0; i.t<m_s.t; i.t++) {
		for(i.z=0; i.z<m_s.z; i.z++) {
			for(i.y=0; i.y<m_s.y; i.y++) {
				for(i.x=0; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);

					const Uint4 ss(uint(i.x > 0), uint(i.y > 0), uint(i.z > 0), uint(i.t > 0));
					for(j.t=0; j.t<=ss.t; j.t++) {
						for(j.z=0; j.z<=ss.z; j.z++) {
							for(j.y=0; j.y<=ss.y; j.y++) {
								for(j.x=0; j.x<=ss.x; j.x++) {
									createPart(ii, k.dot(j));
								}
							}
						}
					}
				}
			}
		}
	}

	computeNodes();
	computeEdges();
	computeFaces();
	computeBodies();
	computeQuads();
}

void BlockMesh::synchronizeNeighbors(const Uint4 &i0, const Uint4 &i1)
{
	// usage of m_prev and m_send is following:
	// send necessary ids to m_prev
	// receive ids from m_next
	// send knowledge to m_next in the same order
	// receive the knowledge from m_prev in given order

	const uint ranks = getMPIranks();
	if(ranks <= 1) return;
	const uint rank = getMPIrank();
	uint i;

	// send own rectangle to others
	Buffer<uint> bnds(2 * m_dim);
	bnds[0] = i0.x, bnds[1] = i1.x;
	if(m_dim >= 2) bnds[2] = i0.y, bnds[3] = i1.y;
	if(m_dim >= 3) bnds[4] = i0.z, bnds[5] = i1.z;
	if(m_dim >= 4) bnds[6] = i0.t, bnds[7] = i1.t;
	for(i=0; i<ranks; i++)
	{
		if(i == rank) continue;
		sendMPI(&bnds[0], bnds.size() * sizeof(uint), i, 0);
	}

	// receive rectangles from all
	const Uint4 ones(1,1,1,1);
	uint prevs = 0;
	m_prev.resize(rank);
	uint nexts = 0;
	m_next.resize(ranks - rank - 1);
	for(i=0; i<ranks; i++)
	{
		if(i == rank) continue;
		recvMPI(&bnds[0], bnds.size() * sizeof(uint), i, 0);

		if(i < rank)
		{
			Uint4 a0((i0.x > bnds[0] ? i0.x : bnds[0] + 1), 0,0,0);
			Uint4 a1((i1.x < bnds[1] ? i1.x : bnds[1]), 0,0,0);
			if(a0.x > a1.x) continue;
			if(m_dim >= 2) {
				a0.y = (i0.y > bnds[2] ? i0.y : bnds[2] + 1);
				a1.y = (i1.y < bnds[3] ? i1.y : bnds[3]);
				if(a0.y > a1.y) continue;
			}
			if(m_dim >= 3) {
				a0.z = (i0.z > bnds[4] ? i0.z : bnds[4] + 1);
				a1.z = (i1.z < bnds[5] ? i1.z : bnds[5]);
				if(a0.z > a1.z) continue;
			}
			if(m_dim >= 4) {
				a0.t = (i0.t > bnds[6] ? i0.t : bnds[6] + 1);
				a1.t = (i1.t < bnds[7] ? i1.t : bnds[7]);
				if(a0.t > a1.t) continue;
			}
			m_prev[prevs++] = RankRectangle(i, a0-i0, ones+a1-i0);
		}
		else
		{
			Uint4 a0((bnds[0] > i0.x ? bnds[0] : i0.x + 1), 0,0,0);
			Uint4 a1((bnds[1] < i1.x ? bnds[1] : i1.x), 0,0,0);
			if(a0.x > a1.x) continue;
			if(m_dim >= 2) {
				a0.y = (bnds[2] > i0.y ? bnds[2] : i0.y + 1);
				a1.y = (bnds[3] < i1.y ? bnds[3] : i1.y);
				if(a0.y > a1.y) continue;
			}
			if(m_dim >= 3) {
				a0.z = (bnds[4] > i0.z ? bnds[4] : i0.z + 1);
				a1.z = (bnds[5] < i1.z ? bnds[5] : i1.z);
				if(a0.z > a1.z) continue;
			}
			if(m_dim >= 4) {
				a0.t = (bnds[6] > i0.t ? bnds[6] : i0.t + 1);
				a1.t = (bnds[7] < i1.t ? bnds[7] : i1.t);
				if(a0.t > a1.t) continue;
			}
			m_next[nexts++] = RankRectangle(i, a0-i0, ones+a1-i0);
		}
	}
	m_prev.resize(prevs);
	m_next.resize(nexts);
}

ullong BlockMesh::getWeight(Uint4 i0, Uint4 i1, const uchar *type, const Uint4 &s, const uint rank0, const uint rank1)
{
	// divide are for ranks
	const uint rank = getMPIrank();
	const uint drank = rank1 - rank0;
	if(drank > 1)
	{
		uint dori = i1.x - i0.x; uint iori = 0;
		if(m_dim > 1 && i1.y - i0.y > dori) { dori = i1.y - i0.y; iori = 1; }
		if(m_dim > 2 && i1.z - i0.z > dori) { dori = i1.z - i0.z; iori = 2; }
		if(m_dim > 3 && i1.t - i0.t > dori) { dori = i1.t - i0.t; iori = 3; }
		if(iori == 0) {
			i1.x = i0.x + dori * (rank + 1 - rank0) / drank;
			i0.x = i0.x + dori * (rank - rank0) / drank;
		}
		else if(iori == 1) {
			i1.y = i0.y + dori * (rank + 1 - rank0) / drank;
			i0.y = i0.y + dori * (rank - rank0) / drank;
		}
		else if(iori == 2) {
			i1.z = i0.z + dori * (rank + 1 - rank0) / drank;
			i0.z = i0.z + dori * (rank - rank0) / drank;
		}
		else {
			i1.t = i0.t + dori * (rank + 1 - rank0) / drank;
			i0.t = i0.t + dori * (rank - rank0) / drank;
		}
	}

	// compute area weight
	ullong size = 0;
	const ullong notinitialized = ullong(-1);
	Buffer<ullong> isize(256, notinitialized);

	Uint4 i;
	for(i.t=i0.t; i.t<i1.t; i.t++) {
		for(i.z=i0.z; i.z<i1.z; i.z++) {
			for(i.y=i0.y; i.y<i1.y; i.y++) {
				for(i.x=i0.x; i.x<i1.x; i.x++) {
					const uchar itype = type[((i.t * s.z + i.z) * s.y + i.y) * s.x + i.x];
					if(isize[itype] == notinitialized) {
						BuilderMesh mesh(m_dim);
						insertNodes(mesh, itype, Vector4(0,0,0,0));
						isize[itype] = mesh.getNodeSize();
					}
					size += isize[itype];
				}
			}
		}
	}

	// sum size with other ranks
	if(drank > 1)
	{
		const uint n = rank1 - rank;
		ullong nsize;
		if(2 * n <= drank) { recvMPI(&nsize, sizeof(ullong), rank1 - 2 * n, 0); size += nsize; }
		if(2 * n + 1 <= drank) { recvMPI(&nsize, sizeof(ullong), rank1 - 2 * n - 1, 0); size += nsize; }
		if(n >= 2) {
			sendMPI(&size, sizeof(ullong), rank1 - n / 2, 0);
			recvMPI(&size, sizeof(ullong), rank1 - n / 2, 1);
		}
		if(2 * n <= drank) sendMPI(&size, sizeof(ullong), rank1 - 2 * n, 1);
		if(2 * n + 1 <= drank) sendMPI(&size, sizeof(ullong), rank1 - 2 * n - 1, 1);
	}
	return size;
}

void BlockMesh::insertNodes(BuilderMesh &mesh, const uint type, const Vector4 &p) const
{
	uint node = 0;
	if(m_dim == 1)
	{
		for(uint i=0; i<type; i++) node = mesh.insertNode(p + Vector4(((i + 0.5) / double(type) - 0.5) * m_d.x, 0,0,0), 0.0, node, true);
		return;
	}
	if(m_dim == 2)
	{
/*		for(uint i=0; i<type; i++) {
			const double px0 = ((i + 0.25) / double(type) - 0.5) * m_d.x;
			const double px1 = ((i + 0.75) / double(type) - 0.5) * m_d.x;
			for(uint j=0; j<type; j++) {
				node = mesh.insertNode(p + Vector4(px0, ((j + 0.25) / double(type) - 0.5) * m_d.y, 0,0), 0.0, node, true);
				node = mesh.insertNode(p + Vector4(px1, ((j + 0.75) / double(type) - 0.5) * m_d.y, 0,0), 0.0, node, true);
			}
		}
*/
		if(type < 10)
		{
			for(uint i=0; i<type; i++) {
				const double px = ((i + 0.5) / double(type) - 0.5) * m_d.x;
				for(uint j=0; j<type; j++) node = mesh.insertNode(p + Vector4(px, ((j + 0.5) / double(type) - 0.5) * m_d.y, 0,0), 0.0, node, true);
			}
		}
		else if(type == 10)
		{
			node = mesh.insertNode(p + Vector4(0.125*m_d.x, 0.25*m_d.y, 0,0), 0.0, node, true);
			node = mesh.insertNode(p + Vector4(-0.125*m_d.x, -0.25*m_d.y, 0,0), 0.0, node, true);
			node = mesh.insertNode(p + Vector4(0.375*m_d.x, -0.25*m_d.y, 0,0), 0.0, node, true);
			node = mesh.insertNode(p + Vector4(-0.375*m_d.x, 0.25*m_d.y, 0,0), 0.0, node, true);
		}
		else if(type == 11)
		{
			node = mesh.insertNode(p + Vector4(m_d.x / 4.0, m_d.y / 12.0, 0,0), 0.0, node, true);
			node = mesh.insertNode(p + Vector4(-m_d.x / 12.0, m_d.y / 4.0, 0,0), 0.0, node, true);
			node = mesh.insertNode(p + Vector4(-m_d.x / 4.0, -m_d.y / 12.0, 0,0), 0.0, node, true);
			node = mesh.insertNode(p + Vector4(m_d.x / 12.0, -m_d.y / 4.0, 0,0), 0.0, node, true);
			node = mesh.insertNode(p + Vector4(m_d.x / 4.0, m_d.y * 5.0 / 12.0, 0,0), 0.0, node, true);
			node = mesh.insertNode(p + Vector4(-m_d.x * 5.0 / 12.0, m_d.y / 4.0, 0,0), 0.0, node, true);
			node = mesh.insertNode(p + Vector4(-m_d.x / 4.0, -m_d.y * 5.0 / 12.0, 0,0), 0.0, node, true);
			node = mesh.insertNode(p + Vector4(m_d.x * 5.0 / 12.0, -m_d.y / 4.0, 0,0), 0.0, node, true);
		}
		return;
	}

	if(type < 10)
	{
		for(uint i=0; i<type; i++) {
			const double px = ((i + 0.5) / double(type) - 0.5) * m_d.x;
			for(uint j=0; j<type; j++) {
				const double py = ((j + 0.5) / double(type) - 0.5) * m_d.y;
				for(uint k=0; k<type; k++) {
					node = mesh.insertNode(p + Vector4(px, py, ((k + 0.5) / double(type) - 0.5) * m_d.z, 0), 0.0, node, true);
				}
			}
		}
	}
/*	else if(type <= 20) // bcc grid
	{
		const uint typei = type - 10;
		const Vector4 d(m_d.toVector3() / double(4 * typei), 0);
		for(uint i=0; i<typei; i++) {
			const double px = ((i + 0.5) / double(typei) - 0.5) * m_d.x;
			for(uint j=0; j<typei; j++) {
				const double py = ((j + 0.5) / double(typei) - 0.5) * m_d.y;
				for(uint k=0; k<typei; k++) {
					const Vector4 p0 = p + Vector4(px, py, ((k + 0.5) / double(typei) - 0.5) * m_d.z, 0);
					node = mesh.insertNode(p0 - d, 0.0, node, true);
					node = mesh.insertNode(p0 + d, 0.0, node, true);
				}
			}
		}
	}
	else if(type <= 30) // fcc grid
	{
		const uint typei = type - 20;
		const Vector4 d(m_d.toVector3() / double(4 * typei), 0);
		for(uint i=0; i<typei; i++) {
			const double px = ((i + 0.5) / double(typei) - 0.5) * m_d.x;
			for(uint j=0; j<typei; j++) {
				const double py = ((j + 0.5) / double(typei) - 0.5) * m_d.y;
				for(uint k=0; k<typei; k++) {
					const Vector4 p0 = p + Vector4(px, py, ((k + 0.5) / double(typei) - 0.5) * m_d.z, 0);
					node = mesh.insertNode(p0 - d, 0.0, node, true);
					node = mesh.insertNode(p0 + Vector4(-d.x, d.y, d.z, 0), 0.0, node, true);
					node = mesh.insertNode(p0 + Vector4(d.x, -d.y, d.z, 0), 0.0, node, true);
					node = mesh.insertNode(p0 + Vector4(d.x, d.y, -d.z, 0), 0.0, node, true);
				}
			}
		}
	}
*//*	else if(type <= 127) // x-y-plane grid
	{
		const uint typei = 2;
		const uint typex = (type - 100) % 3;
		const uint typey = ((type - 100) / 3) % 3;
		const uint typez = (type - 100) / 9;
		const Vector4 d(m_d.toVector3() / double(4 * typei), 0);
		for(uint i=0; i<typei; i++) {
			const double px = ((i + 0.5) / double(typei) - 0.5) * m_d.x;
			for(uint j=0; j<typei; j++) {
				const double py = ((j + 0.5) / double(typei) - 0.5) * m_d.y;
				for(uint k=0; k<typei; k++) {
					Vector4 v0 = Vector4(px, py, ((k + 0.5) / double(typei) - 0.5) * m_d.z, 0) + d;
					if(typex == 0 && v0.x < -d.x) v0.x = -d.x;
					else if(typex == 2 && v0.x > d.x) v0.x = d.x;
					if(typey == 0 && v0.y < -d.y) v0.y = -d.y;
					else if(typey == 2 && v0.y > d.y) v0.y = d.y;
					if(typez == 0 && v0.z < -d.z) v0.z = -d.z;
					else if(typez == 2 && v0.z > d.z) v0.z = d.z;
					node = mesh.insertNode(p + v0, 0.0, node, false);

					Vector4 v1 = Vector4(px, py, ((k + 0.5) / double(typei) - 0.5) * m_d.z, 0) - d;
					if(typex == 0 && v1.x < -d.x) v1.x = -d.x;
					else if(typex == 2 && v1.x > d.x) v1.x = d.x;
					if(typey == 0 && v1.y < -d.y) v1.y = -d.y;
					else if(typey == 2 && v1.y > d.y) v1.y = d.y;
					if(typez == 0 && v1.z < -d.z) v1.z = -d.z;
					else if(typez == 2 && v1.z > d.z) v1.z = d.z;
					node = mesh.insertNode(p + v1, 0.0, node, false);

				}
			}
		}
	}
	else if(type <= 144) // x-y-plane grid
	{
		const uint typei = 2;
		const uint typex = (type - 130) % 3;
		const uint typey = (type - 130) / 3;
		const Vector4 d(m_d.toVector3() / double(4 * typei), 0);
		for(uint i=0; i<typei; i++) {
			const double px = ((i + 0.5) / double(typei) - 0.5) * m_d.x;
			for(uint j=0; j<typei; j++) {
				const double py = ((j + 0.5) / double(typei) - 0.5) * m_d.y;
				for(uint k=0; k<typei; k++) {
					Vector4 v0 = Vector4(px, py, ((k + 0.5) / double(typei) - 0.5) * m_d.z, 0) + d;
					double z0 = -1e30;
					if(typey <= 2 && v0.x > z0) z0 = v0.x;
					if(typey >= 2 && v0.y > z0) z0 = v0.y;
					if(typex == 0 && d.z > z0) z0 = d.z;
					else if(typex == 2 && d.z < z0) z0 = d.z;

					if(v0.z > z0) v0.z = z0;
					if(typey == 0 && v0.y < -d.y) v0.y = -d.y;
					else if(typey == 4 && v0.x < -d.x) v0.x = -d.x;
					node = mesh.insertNode(p + v0, 0.0, node, false);

					Vector4 v1 = Vector4(px, py, ((k + 0.5) / double(typei) - 0.5) * m_d.z, 0) - d;
					double z1 = -1e30;
					if(typey <= 2 && v1.x > z1) z1 = v1.x;
					if(typey >= 2 && v1.y > z1) z1 = v1.y;
					if(typex == 0 && d.z > z1) z1 = d.z;
					else if(typex == 2 && d.z < z1) z1 = d.z;

					if(v1.z > z1) v1.z = z1;
					if(typey == 0 && v1.y < -d.y) v1.y = -d.y;
					else if(typey == 4 && v1.x < -d.x) v1.x = -d.x;
					node = mesh.insertNode(p + v1, 0.0, node, false);
				}
			}
		}
	}
*/	else if(type >= 10 && type < 118)
	{
		const uint num = 5;
		const uint typei = (type - 10) % 4;

		const uint typex = ((type - 10) / 4) % 3;
		const uint typey = ((type - 10) / 12) % 3;
		const uint typez = (type - 10) / 36;
//		const Vector4 d(m_d.toVector3() / double(4 * typei), 0);
		for(uint i=0; i<num; i++) {
			const double px = ((i + 0.5) / double(num) - 0.5) * m_d.x;
			for(uint j=0; j<num; j++) {
				const double py = ((j + 0.5) / double(num) - 0.5) * m_d.y;
				for(uint k=0; k<num; k++) {
					const uint ijk = (i % 2) + 2 * (j % 2) + 4 * (k % 2);
					if(ijk != typei && ijk != 7 - typei) continue;
					Vector4 v = Vector4(px, py, ((k + 0.5) / double(num) - 0.5) * m_d.z, 0);
					if(typex == 0 && v.x < 0.0) v.x = 0.0;
					else if(typex == 2 && v.x > 0.0) v.x = 0.0;
					if(typey == 0 && v.y < 0.0) v.y = 0.0;
					else if(typey == 2 && v.y > 0.0) v.y = 0.0;
					if(typez == 0 && v.z < 0.0) v.z = 0.0;
					else if(typez == 2 && v.z > 0.0) v.z = 0.0;
					node = mesh.insertNode(p + v, 0.0, node, false);
				}
			}
		}
	}
	else if(type >= 118 && type < 226)
	{
		const uint num = 5;
		const uint typei = (type - 118) % 4;

		const uint typex = ((type - 118) / 4) % 3;
		const uint typey = ((type - 118) / 12) % 3;
		const uint typez = (type - 118) / 36;
		for(uint i=0; i<num; i++) {
			const double px = ((i + 0.5) / double(num) - 0.5) * m_d.x;
			for(uint j=0; j<num; j++) {
				const double py = ((j + 0.5) / double(num) - 0.5) * m_d.y;
				for(uint k=0; k<num; k++) {
					const uint ijk = (i % 2) + 2 * (j % 2) + 4 * (k % 2);
					if(ijk != typei && ijk != 7 - typei) continue;
					Vector4 v = Vector4(px, py, ((k + 0.5) / double(num) - 0.5) * m_d.z, 0);

					double z = 1e30;
					if(typex == 0 && m_d.x * z > m_d.z * v.x) z = m_d.z / m_d.x * v.x;
					else if(typex == 2 && m_d.x * z > -m_d.z * v.x) z = -m_d.z / m_d.x * v.x;
					if(typey == 0 && m_d.y * z > m_d.z * v.y) z = m_d.z / m_d.y * v.y;
					else if(typey == 2 && m_d.y * z > -m_d.z * v.y) z = -m_d.z / m_d.y * v.y;
					if(typez == 0 && z < 0.0) z = 0.0;
					else if(typez == 2 && z > 0.0) z = 0.0;

					if(v.z < z) v.z = z;
					node = mesh.insertNode(p + v, 0.0, node, false);
				}
			}
		}
	}
}

void BlockMesh::removeNonconvex(BuilderMesh &mesh, const uint type, const Vector4 &p) const
{
	if(m_dim == 1) return;
	if(m_dim == 2) return;

	if(type >= 118 && type < 226)
	{
		uint i;
		const uint typex = ((type - 118) / 4) % 3;
		const uint typey = ((type - 118) / 12) % 3;
		const uint typez = (type - 118) / 36;

		if(typez < 2 && (typex == 1 || typey == 1)) return;

		for(i=mesh.getEdgeSize(); i-->0; ) {
			const Vector4 v = mesh.getEdgeAverage(i) - p;
			if(v.x < -0.5 * m_d.x || 0.5 * m_d.x < v.x || v.y < -0.5 * m_d.y || 0.5 * m_d.y < v.y) continue;

			double z = 1e30;
			if(typex == 0 && m_d.x * z > m_d.z * v.x) z = m_d.z / m_d.x * v.x;
			else if(typex == 2 && m_d.x * z > -m_d.z * v.x) z = -m_d.z / m_d.x * v.x;
			if(typey == 0 && m_d.y * z > m_d.z * v.y) z = m_d.z / m_d.y * v.y;
			else if(typey == 2 && m_d.y * z > -m_d.z * v.y) z = -m_d.z / m_d.y * v.y;
			if(typez == 0 && z < 0.0) z = 0.0;
			else if(typez == 2 && z > 0.0) z = 0.0;

			if(v.z < z - 1e-8) mesh.removeEdge(i);
		}
	}
}

const PartMesh *BlockMesh::createPart(const uint ii, const uint part)
{
	const Type type = getType(ii, part);
	map<Type, PartMesh*>::iterator it = m_part[part].find(type);
	if(it != m_part[part].end()) return it->second;

	// create new mesh
	const Uint4 s(part % 2, (part / 2) % 2, (part / 4) % 2, part / 8);
	PartMesh *pm = new PartMesh(Uint4(1,3,9,27).dot(s), PARTS, m_dim);
	m_part[part].insert(std::make_pair(type, pm));

	// create empty mesh, if any of the corners is of type zero
	Uint4 i;
/*	for(i.t=0; i.t<=s.t; i.t++) {
		for(i.z=0; i.z<=s.z; i.z++) {
			for(i.y=0; i.y<=s.y; i.y++) {
				for(i.x=0; i.x<=s.x; i.x++) {
					if(m_type[ii - getBlockId(i)] == 0) return pm;
				}
			}
		}
	}
*/
	// create mesh and insert nodes
	const Uint4 k(2,6,18,54);
	BuilderMesh mesh(m_dim);
	Buffer<uint> nbeg(PARTS, 0);
	for(i.t=0; i.t<=s.t; i.t++) {
		for(i.z=0; i.z<=s.z; i.z++) {
			for(i.y=0; i.y<=s.y; i.y++) {
				for(i.x=0; i.x<=s.x; i.x++) {
					nbeg[k.dot(i)] = mesh.getNodeSize();
					insertNodes(mesh, m_type[ii - getBlockId(i)], -Vector4(i.x*m_d.x,i.y*m_d.y,i.z*m_d.z,i.t*m_d.t));
				}
			}
		}
	}
	for(i.t=0; i.t<=s.t; i.t++) {
		for(i.z=0; i.z<=s.z; i.z++) {
			for(i.y=0; i.y<=s.y; i.y++) {
				for(i.x=0; i.x<=s.x; i.x++) {
					removeNonconvex(mesh, m_type[ii - getBlockId(i)], -Vector4(i.x*m_d.x,i.y*m_d.y,i.z*m_d.z,i.t*m_d.t));
				}
			}
		}
	}
	mesh.fillFlags(pm->getPart());

	// mark existing elements
	const uint dim = s.x + s.y + s.z + s.t;
	if(dim > 0)	{
		for(i.t=0; i.t<=s.t; i.t++) {
			for(i.z=0; i.z<=s.z; i.z++) {
				for(i.y=0; i.y<=s.y; i.y++) {
					for(i.x=0; i.x<=s.x; i.x++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 0), nbeg, k.dot(i), PARTS);
					}
				}
			}
		}
	}
	if(dim > 1)	{
		if(s.x > 0)	{ i.x = 0;
			for(i.t=0; i.t<=s.t; i.t++) {
				for(i.z=0; i.z<=s.z; i.z++) {
					for(i.y=0; i.y<=s.y; i.y++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 1), nbeg, k.dot(i), PARTS);
					}
				}
			}
		}
		if(s.y > 0)	{ i.y = 0;
			for(i.t=0; i.t<=s.t; i.t++) {
				for(i.z=0; i.z<=s.z; i.z++) {
					for(i.x=0; i.x<=s.x; i.x++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 2), nbeg, k.dot(i), PARTS);
					}
				}
			}
		}
		if(s.z > 0)	{ i.z = 0;
			for(i.t=0; i.t<=s.t; i.t++) {
				for(i.y=0; i.y<=s.y; i.y++) {
					for(i.x=0; i.x<=s.x; i.x++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 4), nbeg, k.dot(i), PARTS);
					}
				}
			}
		}
		if(s.t > 0)	{ i.t = 0;
			for(i.z=0; i.z<=s.z; i.z++) {
				for(i.y=0; i.y<=s.y; i.y++) {
					for(i.x=0; i.x<=s.x; i.x++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 8), nbeg, k.dot(i), PARTS);
					}
				}
			}
		}
	}
	if(dim > 2)	{
		if(s.x > 0)	{ i.x = 0;
			if(s.y > 0)	{ i.y = 0;
				for(i.t=0; i.t<=s.t; i.t++) {
					for(i.z=0; i.z<=s.z; i.z++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 3), nbeg, k.dot(i), PARTS);
					}
				}
			}
			if(s.z > 0)	{ i.z = 0;
				for(i.t=0; i.t<=s.t; i.t++) {
					for(i.y=0; i.y<=s.y; i.y++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 5), nbeg, k.dot(i), PARTS);
					}
				}
			}
			if(s.t > 0)	{ i.t = 0;
				for(i.z=0; i.z<=s.z; i.z++) {
					for(i.y=0; i.y<=s.y; i.y++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 9), nbeg, k.dot(i), PARTS);
					}
				}
			}
		}
		if(s.y > 0) { i.y = 0;
			if(s.z > 0)	{ i.z = 0;
				for(i.t=0; i.t<=s.t; i.t++) {
					for(i.x=0; i.x<=s.x; i.x++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 6), nbeg, k.dot(i), PARTS);
					}
				}
			}
			if(s.t > 0)	{ i.t = 0;
				for(i.z=0; i.z<=s.z; i.z++) {
					for(i.x=0; i.x<=s.x; i.x++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 10), nbeg, k.dot(i), PARTS);
					}
				}
			}
		}
		if(s.z > 0) { i.z = 0;
			if(s.t > 0)	{ i.t = 0;
				for(i.y=0; i.y<=s.y; i.y++) {
					for(i.x=0; i.x<=s.x; i.x++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 12), nbeg, k.dot(i), PARTS);
					}
				}
			}
		}
	}
	if(dim > 3)	{
		if(s.x > 0)	{ i.x = 0;
			if(s.y > 0)	{ i.y = 0;
				if(s.z > 0)	{ i.z = 0;
					for(i.t=0; i.t<=s.t; i.t++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 7), nbeg, k.dot(i), PARTS);
					}
				}
				if(s.t > 0)	{ i.t = 0;
					for(i.z=0; i.z<=s.z; i.z++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 11), nbeg, k.dot(i), PARTS);
					}
				}
			}
			if(s.z > 0)	{ i.z = 0;
				if(s.t > 0)	{ i.t = 0;
					for(i.y=0; i.y<=s.y; i.y++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 13), nbeg, k.dot(i), PARTS);
					}
				}
			}
		}
		if(s.y > 0)	{ i.y = 0;
			if(s.z > 0)	{ i.z = 0;
				if(s.t > 0)	{ i.t = 0;
					for(i.x=0; i.x<=s.x; i.x++) {
						setPartFlags(mesh, getPart(ii - getBlockId(i), 14), nbeg, k.dot(i), PARTS);
					}
				}
			}
		}
	}

	cutSecure(mesh, -0.5*m_d-Vector4(s.x*m_d.x,s.y*m_d.y,s.z*m_d.z,s.t*m_d.t), 0.5 * m_d, pm->getPart());
/*	for(i.t=0; i.t<=s.t; i.t++) {
		for(i.z=0; i.z<=s.z; i.z++) {
			for(i.y=0; i.y<=s.y; i.y++) {
				for(i.x=0; i.x<=s.x; i.x++) {
					if(m_type[ii - getBlockId(i)] == 0) {
						const Vector4 d(i.x*m_d.x,i.y*m_d.y,i.z*m_d.z,i.t*m_d.t);
						cutZero(mesh, -0.5*m_d-d, 0.5 * m_d-d, pm->getPart());
					}
				}
			}
		}
	}
*/	pm->createPartFromFlags(mesh);
	return pm;
}

Type BlockMesh::getType(const uint ii) const
{
	const uint x = 1;
	if(m_dim == 1) return Type(m_type[ii],m_type[ii-x]); // part == 1
	const uint y = m_s.x;
	if(m_dim == 2) return Type(m_type[ii],m_type[ii-x],m_type[ii-y],m_type[ii-x-y]); // part == 3
	const uint z = y * m_s.y;
	if(m_dim == 3) return Type(m_type[ii],m_type[ii-x],m_type[ii-y],m_type[ii-x-y],m_type[ii-z],m_type[ii-x-z],m_type[ii-y-z],m_type[ii-x-y-z]); // part == 7
	const uint t = z * m_s.z;
	return Type(m_type[ii], m_type[ii-x],m_type[ii-y],m_type[ii-x-y],m_type[ii-z],m_type[ii-x-z],m_type[ii-y-z],m_type[ii-x-y-z],
		m_type[ii-t], m_type[ii-x-t],m_type[ii-y-t],m_type[ii-x-y-t],m_type[ii-z-t],m_type[ii-x-z-t],m_type[ii-y-z-t],m_type[ii-x-y-z-t]); // part == 15
}
Type BlockMesh::getType(const uint ii, const uint part) const
{
	const uint x = 1;
	const uint y = m_s.x;
	if(part < 8) {
		if(part < 4) {
			if(part < 2) {
				if(part == 0) return Type(m_type[ii]); // part == 0
				return Type(m_type[ii],m_type[ii-x]); // part == 1
			}
			if(part == 2) return Type(m_type[ii],m_type[ii-y]); // part == 2
			return Type(m_type[ii],m_type[ii-x],m_type[ii-y],m_type[ii-x-y]); // part == 3
		}
		const uint z = m_s.x * m_s.y;
		if(part < 6) {
			if(part == 4) return Type(m_type[ii],m_type[ii-z]); // part == 4
			return Type(m_type[ii],m_type[ii-x],m_type[ii-z],m_type[ii-x-z]); // part == 5
		}
		if(part == 6) return Type(m_type[ii],m_type[ii-y],m_type[ii-z],m_type[ii-y-z]); // part == 6
		return Type(m_type[ii],m_type[ii-x],m_type[ii-y],m_type[ii-x-y],m_type[ii-z],m_type[ii-x-z],m_type[ii-y-z],m_type[ii-x-y-z]); // part == 7
	}
	const uint t = m_s.x * m_s.y * m_s.z;
	if(part < 12) {
		if(part < 10) {
			if(part == 8) return Type(m_type[ii], m_type[ii-t]); // part == 8
			return Type(m_type[ii],m_type[ii-x],m_type[ii-t],m_type[ii-x-t]); // part == 9
		}
		if(part == 10) return Type(m_type[ii],m_type[ii-y],m_type[ii-t],m_type[ii-y-t]); // part == 10
		return Type(m_type[ii],m_type[ii-x],m_type[ii-y],m_type[ii-x-y],m_type[ii-t],m_type[ii-x-t],m_type[ii-y-t],m_type[ii-x-y-t]); // part == 11
	}
	const uint z = m_s.x * m_s.y;
	if(part < 14) {
		if(part == 12) return Type(m_type[ii],m_type[ii-z],m_type[ii-t],m_type[ii-z-t]); // part == 12
		return Type(m_type[ii],m_type[ii-x],m_type[ii-z],m_type[ii-x-z],m_type[ii-t],m_type[ii-x-t],m_type[ii-z-t],m_type[ii-x-z-t]); // part == 13
	}
	if(part == 14) return Type(m_type[ii],m_type[ii-y],m_type[ii-z],m_type[ii-y-z],m_type[ii-t],m_type[ii-y-t],m_type[ii-z-t],m_type[ii-y-z-t]); // part == 14
	return Type(m_type[ii], m_type[ii-x],m_type[ii-y],m_type[ii-x-y],m_type[ii-z],m_type[ii-x-z],m_type[ii-y-z],m_type[ii-x-y-z],
		m_type[ii-t], m_type[ii-x-t],m_type[ii-y-t],m_type[ii-x-y-t],m_type[ii-z-t],m_type[ii-x-z-t],m_type[ii-y-z-t],m_type[ii-x-y-z-t]); // part == 15
}

void BlockMesh::setPartFlags(BuilderMesh &mesh, const PartMesh *part, const Buffer<uint> &nbeg, const uint flag, const uint flags) const
{
	uint i, j;

	// nodes
	Buffer<uint> ns(part->getNodeSize());
	for(i=0; i<ns.size(); i++)
	{
		ns[i] = nbeg[part->getNodePart(i) + flag] + part->getNodeLink(i);
		if(i < part->getNodeLocals()) mesh.setNodeFlag(ns[i], part->getPart() + flag + flags * i);
	}

	// edges
	Buffer<uint> es(part->getEdgeSize());
	for(i=0; i<es.size(); i++)
	{
		Buffer<uint> n = part->getEdgeNodes(i);
		for(j=0; j<n.size(); j++) n[j] = ns[n[j]];
		es[i] = mesh.findEdge(n[0], n[1]);
		if(es[i] == NONE) cout << "BlockMesh::setPartFlags -> Generation failed: Edge not found." << endl;
		else if(i < part->getEdgeLocals()) {
			mesh.setEdgeNodes(es[i], n);
			mesh.setEdgeFlag(es[i], part->getPart() + flag + flags * i);
		}
	}

	// faces
	Buffer<uint> fs(part->getFaceSize());
	for(i=0; i<fs.size(); i++)
	{
		Buffer<uint> e = part->getFaceEdges(i);
		for(j=0; j<e.size(); j++) e[j] = es[e[j]];
		fs[i] = mesh.findFace(e);
		if(fs[i] == NONE) cout << "BlockMesh::setPartFlags -> Generation failed: Face not found." << endl;
		else if(i < part->getFaceLocals()) {
			mesh.setFaceEdges(fs[i], e);
			mesh.setFaceFlag(fs[i], part->getPart() + flag + flags * i);
		}
	}

	// bodies
	Buffer<uint> bs(part->getBodySize());
	for(i=0; i<bs.size(); i++)
	{
		Buffer<uint> f = part->getBodyFaces(i);
		for(j=0; j<f.size(); j++) f[j] = fs[f[j]];
		bs[i] = mesh.findBody(f);
		if(bs[i] == NONE) cout << "BlockMesh::setPartFlags -> Generation failed: Body not found." << endl;
		else if(i < part->getBodyLocals()) {
			mesh.setBodyFaces(bs[i], f);
			mesh.setBodyFlag(bs[i], part->getPart() + flag + flags * i);
		}
	}

	// quads
	for(i=0; i<part->getQuadLocals(); i++)
	{
		Buffer<uint> b = part->getQuadBodies(i);
		for(j=0; j<b.size(); j++) b[j] = bs[b[j]];
		const uint qs = mesh.findQuad(b);
		if(qs == NONE) cout << "BlockMesh::setPartFlags -> Generation failed: Quad not found." << endl;
		else {
			mesh.setQuadBodies(qs, b);
			mesh.setQuadFlag(qs, part->getPart() + flag + flags * i);
		}
	}
}

void BlockMesh::cutSecure(BuilderMesh &mesh, const Vector4 &p0, const Vector4 &p1, const uint flag) const
{
	uint i;
	for(i=mesh.getEdgeSize(); i-->0; ) {
		const Vector4 p = mesh.getEdgePosition(i);
		const double r = sqrt(mesh.getRadiusSq(p, mesh.getEdgeNodes(i))) - 1e-13;
		if(p.x - r < p0.x || p.x + r > p1.x || p.y - r < p0.y || p.y + r > p1.y ||
			p.z - r < p0.z || p.z + r > p1.z || p.t - r < p0.t || p.t + r > p1.t) mesh.removeEdge(i); // the edge is out of bounds
	}
	for(i=mesh.getFaceSize(); i-->0; ) {
		const Vector4 p = mesh.getFacePosition(i);
		const double r = sqrt(mesh.getRadiusSq(p, mesh.getFaceNodes(i))) - 1e-13;
		if(p.x - r < p0.x || p.x + r > p1.x || p.y - r < p0.y || p.y + r > p1.y ||
			p.z - r < p0.z || p.z + r > p1.z || p.t - r < p0.t || p.t + r > p1.t) mesh.removeFace(i); // the face is out of bounds
	}
	for(i=mesh.getBodySize(); i-->0; ) {
		const Vector4 p = mesh.getBodyPosition(i);
		const double r = sqrt(mesh.getRadiusSq(p, mesh.getBodyNodes(i))) - 1e-13;
		if(p.x - r < p0.x || p.x + r > p1.x || p.y - r < p0.y || p.y + r > p1.y ||
			p.z - r < p0.z || p.z + r > p1.z || p.t - r < p0.t || p.t + r > p1.t) mesh.removeBody(i); // the body is out of bounds
	}
	for(i=mesh.getQuadSize(); i-->0; ) {
		const Vector4 p = mesh.getQuadPosition(i);
		const double r = sqrt(mesh.getRadiusSq(p, mesh.getQuadNodes(i))) - 1e-13;
		if(p.x - r < p0.x || p.x + r > p1.x || p.y - r < p0.y || p.y + r > p1.y ||
			p.z - r < p0.z || p.z + r > p1.z || p.t - r < p0.t || p.t + r > p1.t) mesh.removeQuad(i); // the quad is out of bounds
	}
	for(i=mesh.getQuadSize(); i-->0; ) {
		if(mesh.getQuadFlag(i) != flag) mesh.removeQuad(i); // the quad belongs to another mesh
	}
	for(i=mesh.getBodySize(); i-->0; ) {
		if(mesh.getBodyFlag(i) != flag && mesh.getBodyQuads(i).empty()) mesh.removeBody(i); // the body belongs to another mesh
	}
	for(i=mesh.getFaceSize(); i-->0; ) {
		if(mesh.getFaceFlag(i) != flag && mesh.getFaceBodies(i).empty()) mesh.removeFace(i); // the face belongs to another mesh
	}
	for(i=mesh.getEdgeSize(); i-->0; ) {
		if(mesh.getEdgeFlag(i) != flag && mesh.getEdgeFaces(i).empty()) mesh.removeEdge(i); // the edge belongs to another mesh
	}
	for(i=mesh.getNodeSize(); i-->0; ) {
		if(mesh.getNodeFlag(i) != flag && mesh.getNodeEdges(i).empty()) mesh.removeNode(i); // the node belongs to another mesh
	}
/*
	uint i;
	for(i=mesh.getQuadSize(); i-->0; )
	{
		if(mesh.getQuadFlag(i) != flag) mesh.removeQuad(i); // the quad belongs to another mesh
		else
		{
			const Vector4 p = mesh.getQuadPosition(i);
			const double r = sqrt(mesh.getRadiusSq(p, mesh.getQuadNodes(i)))+1e-13;
			if(p.x - r < p0.x || p.x + r > p1.x || p.y - r < p0.y || p.y + r > p1.y ||
				p.z - r < p0.z || p.z + r > p1.z || p.t - r < p0.t || p.t + r > p1.t) mesh.removeQuad(i); // the quad is out of bounds
		}
	}
	for(i=mesh.getBodySize(); i-->0; )
	{
		if(!mesh.getBodyQuads(i).empty()) continue; // the body is necessary for a quad
		if(mesh.getBodyFlag(i) != flag) mesh.removeBody(i); // the body belongs to another mesh
		else
		{
			const Vector4 p = mesh.getBodyPosition(i);
			const double r = sqrt(mesh.getRadiusSq(p, mesh.getBodyNodes(i)))+1e-13;
			if(p.x - r < p0.x || p.x + r > p1.x || p.y - r < p0.y || p.y + r > p1.y ||
				p.z - r < p0.z || p.z + r > p1.z || p.t - r < p0.t || p.t + r > p1.t) mesh.removeBody(i); // the body is out of bounds
		}
	}
	for(i=mesh.getFaceSize(); i-->0; )
	{
		if(!mesh.getFaceBodies(i).empty()) continue; // the face is necessary for a body
		if(mesh.getFaceFlag(i) != flag) mesh.removeFace(i); // the face belongs to another mesh
		else
		{
			const Vector4 p = mesh.getFacePosition(i);
			const double r = sqrt(mesh.getRadiusSq(p, mesh.getFaceNodes(i)))+1e-13;
			if(p.x - r < p0.x || p.x + r > p1.x || p.y - r < p0.y || p.y + r > p1.y ||
				p.z - r < p0.z || p.z + r > p1.z || p.t - r < p0.t || p.t + r > p1.t) mesh.removeFace(i); // the face is out of bounds
		}
	}
	for(i=mesh.getEdgeSize(); i-->0; )
	{
		if(!mesh.getEdgeFaces(i).empty()) continue; // the edge is necessary for a face
		if(mesh.getEdgeFlag(i) != flag) mesh.removeEdge(i); // the edge belongs to another mesh
		else
		{
			const Vector4 p = mesh.getEdgePosition(i);
			const double r = sqrt(mesh.getRadiusSq(p, mesh.getEdgeNodes(i)))+1e-13;
			if(p.x - r < p0.x || p.x + r > p1.x || p.y - r < p0.y || p.y + r > p1.y ||
				p.z - r < p0.z || p.z + r > p1.z || p.t - r < p0.t || p.t + r > p1.t) mesh.removeEdge(i); // the edge is out of bounds
		}
	}
	for(i=mesh.getNodeSize(); i-->0; )
	{
		if(!mesh.getNodeEdges(i).empty()) continue; // the node is necessary for a edge
		if(mesh.getNodeFlag(i) != flag) mesh.removeNode(i); // the node belongs to another mesh
	}
*/
}
/*
void BlockMesh::cutZero(BuilderMesh &mesh, const Vector4 &p0, const Vector4 &p1, const uint flag) const
{
	uint i;
	for(i=mesh.getEdgeSize(); i-->0; ) {
		const Vector4 p = mesh.getEdgePosition(i);
		const double r = 0.0;//sqrt(mesh.getRadiusSq(p, mesh.getEdgeNodes(i))) - 1e-13;
		if(p0.x - r < p.x && p.x < p1.x + r && p0.y - r < p.y && p.y < p1.y + r &&
			p0.z - r < p.z && p.z < p1.z + r && p0.t - r < p.t && p.t < p1.t + r) mesh.removeEdge(i); // the edge is out of bounds
	}
	for(i=mesh.getFaceSize(); i-->0; ) {
		const Vector4 p = mesh.getFacePosition(i);
		const double r = 0.0;//sqrt(mesh.getRadiusSq(p, mesh.getFaceNodes(i))) - 1e-13;
		if(p0.x - r < p.x && p.x < p1.x + r && p0.y - r < p.y && p.y < p1.y + r &&
			p0.z - r < p.z && p.z < p1.z + r && p0.t - r < p.t && p.t < p1.t + r) mesh.removeFace(i); // the face is out of bounds
	}
	for(i=mesh.getBodySize(); i-->0; ) {
		const Vector4 p = mesh.getBodyPosition(i);
		const double r = 0.0;//sqrt(mesh.getRadiusSq(p, mesh.getBodyNodes(i))) - 1e-13;
		if(p0.x - r < p.x && p.x < p1.x + r && p0.y - r < p.y && p.y < p1.y + r &&
			p0.z - r < p.z && p.z < p1.z + r && p0.t - r < p.t && p.t < p1.t + r) mesh.removeBody(i); // the body is out of bounds
	}
	for(i=mesh.getQuadSize(); i-->0; ) {
		const Vector4 p = mesh.getQuadPosition(i);
		const double r = 0.0;//sqrt(mesh.getRadiusSq(p, mesh.getQuadNodes(i))) - 1e-13;
		if(p0.x - r < p.x && p.x < p1.x + r && p0.y - r < p.y && p.y < p1.y + r &&
			p0.z - r < p.z && p.z < p1.z + r && p0.t - r < p.t && p.t < p1.t + r) mesh.removeQuad(i); // the quad is out of bounds
	}
	for(i=mesh.getQuadSize(); i-->0; ) {
		if(mesh.getQuadFlag(i) != flag) mesh.removeQuad(i); // the quad belongs to another mesh
	}
	for(i=mesh.getBodySize(); i-->0; ) {
		if(mesh.getBodyFlag(i) != flag && mesh.getBodyQuads(i).empty()) mesh.removeBody(i); // the body belongs to another mesh
	}
	for(i=mesh.getFaceSize(); i-->0; ) {
		if(mesh.getFaceFlag(i) != flag && mesh.getFaceBodies(i).empty()) mesh.removeFace(i); // the face belongs to another mesh
	}
	for(i=mesh.getEdgeSize(); i-->0; ) {
		if(mesh.getEdgeFlag(i) != flag && mesh.getEdgeFaces(i).empty()) mesh.removeEdge(i); // the edge belongs to another mesh
	}
	for(i=mesh.getNodeSize(); i-->0; ) {
		if(mesh.getNodeFlag(i) != flag && mesh.getNodeEdges(i).empty()) mesh.removeNode(i); // the node belongs to another mesh
	}
}
*/

void BlockMesh::computeNodes()
{
	if(!m_nbeg.empty()) return; // already computed

	Uint4 i;
	uint j;
	m_nsize = 0;
	m_nbeg.resize(m_s.x * m_s.y * m_s.z * m_s.t);
	m_nbeg.fill(0);

	// local nodes
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					m_nbeg[ii] = m_nsize;
					m_nsize += getPart(ii, 0)->getNodeLocals();
				}
			}
		}
	}
	m_nlocs = m_nsize;

	// external nodes
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						const uint ii = getBlockId(i);
						m_nbeg[ii] = m_nsize;
						m_nsize += getNodeLinks(i).size();
					}
				}
			}
		}
	}
}

void BlockMesh::computeEdges()
{
	if(!m_ebeg.empty()) return; // already computed

	Uint4 i;
	uint j, k;
	m_esize = 0;
	m_ebeg.resize(m_s.x * m_s.y * m_s.z * m_s.t);
	m_ebeg.fill(0);

	// local edges
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					m_ebeg[ii] = m_esize;
					for(k=0; k<m_part.size(); k++) m_esize += getPart(ii, k)->getEdgeLocals();
				}
			}
		}
	}
	m_elocs = m_esize;

	// external edges
	if(m_dim < 2) return;
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						const uint ii = getBlockId(i);
						m_ebeg[ii] = m_esize;
						for(k=0; k+1<m_part.size(); k++) m_esize += getEdgeLinks(i, k).size();
					}
				}
			}
		}
	}
}

void BlockMesh::computeFaces()
{
	if(m_dim < 2) return; // no faces exist
	if(!m_fbeg.empty()) return; // already computed

	Uint4 i;
	uint j, k;
	m_fsize = 0;
	m_fbeg.resize(m_s.x * m_s.y * m_s.z * m_s.t);
	m_fbeg.fill(0);

	// local faces
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					m_fbeg[ii] = m_fsize;
					for(j=0; j<m_part.size(); j++) m_fsize += getPart(ii, j)->getFaceLocals();
				}
			}
		}
	}
	m_flocs = m_fsize;

	// external faces
	if(m_dim < 3) return;
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						const uint ii = getBlockId(i);
						m_fbeg[ii] = m_fsize;
						for(k=0; k+1<m_part.size(); k++) m_fsize += getFaceLinks(i, k).size();
					}
				}
			}
		}
	}
}

void BlockMesh::computeBodies()
{
	if(m_dim < 3) return; // no bodies exist
	if(!m_bbeg.empty()) return; // already computed

	Uint4 i;
	uint j, k;
	m_bsize = 0;
	m_bbeg.resize(m_s.x * m_s.y * m_s.z * m_s.t);
	m_bbeg.fill(0);

	// local bodies
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					m_bbeg[ii] = m_bsize;
					for(j=0; j<m_part.size(); j++) m_bsize += getPart(ii, j)->getBodyLocals();
				}
			}
		}
	}
	m_blocs = m_bsize;

	// external faces
	if(m_dim < 4) return;
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						const uint ii = getBlockId(i);
						m_fbeg[ii] = m_bsize;
						for(k=0; k+1<m_part.size(); k++) m_bsize += getBodyLinks(i, k).size();
					}
				}
			}
		}
	}
}

void BlockMesh::computeQuads()
{
	if(m_dim < 4) return; // no quads exist
	if(!m_qbeg.empty()) return; // already computed

	Uint4 i;
	uint j;
	m_qsize = 0;
	m_qbeg.resize(m_s.x * m_s.y * m_s.z * m_s.t);
	m_qbeg.fill(0);

	// local bodies
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					m_qbeg[ii] = m_qsize;
					for(j=0; j<m_part.size(); j++) m_qsize += getPart(ii, j)->getQuadLocals();
				}
			}
		}
	}
	m_qlocs = m_qsize;
}


void BlockMesh::toMesh(PartMesh &mesh) const
{
	mesh.clear();
	Uint4 i;
	uint j, k, l;
	Buffer<uint> link;
	const PartMesh *part = NULL;

	// add local nodes
	mesh.resizeNodeBuffer(m_nsize);
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					const PartMesh *part = getPart(ii, 0);
					for(k=0; k<part->getNodeLocals(); k++) addNodeToMesh(mesh, part, k, i);
				}
			}
		}
	}

	// add external nodes
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						// add external nodes
						const uint ii = getBlockId(i);
						Buffer<uint> link = getNodeLinks(i);
						const PartMesh *part = getPart(ii, 0);
						for(k=0; k<link.size(); k++) addNodeToMesh(mesh, part, link[k], i);
					}
				}
			}
		}
	}
	mesh.setExternalNodes(getExternalNodes());

	// add local edges
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					for(j=0; j<m_part.size(); j++) {
						const PartMesh *part = getPart(ii, j);
						for(k=0; k<part->getEdgeLocals(); k++) addEdgeToMesh(mesh, part, k, i);
					}
				}
			}
		}
	}

	// add external edges
	if(m_dim < 2) return;
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						const uint ii = getBlockId(i);
						for(k=0; k+1<m_part.size(); k++) {
							link = getEdgeLinks(i, k);
							if(!link.empty()) part = getPart(ii, k);
							for(l=0; l<link.size(); l++) addEdgeToMesh(mesh, part, link[l], i);
						}
					}
				}
			}
		}
	}
	mesh.setExternalEdges(getExternalEdges());

	// add local faces
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					for(j=0; j<m_part.size(); j++) {
						const PartMesh *part = getPart(ii, j);
						for(k=0; k<part->getFaceLocals(); k++) addFaceToMesh(mesh, part, k, i);
					}
				}
			}
		}
	}

	// add external faces
	if(m_dim < 3) return;
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						const uint ii = getBlockId(i);
						for(k=0; k+1<m_part.size(); k++) {
							link = getFaceLinks(i, k);
							if(!link.empty()) part = getPart(ii, k);
							for(l=0; l<link.size(); l++) addFaceToMesh(mesh, part, link[l], i);
						}
					}
				}
			}
		}
	}
	mesh.setExternalFaces(getExternalFaces());

	// add local bodies
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					for(j=0; j<m_part.size(); j++) {
						const PartMesh *part = getPart(ii, j);
						for(k=0; k<part->getBodyLocals(); k++) addBodyToMesh(mesh, part, k, i);
					}
				}
			}
		}
	}

	// add external bodies
	if(m_dim < 4) return;
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						const uint ii = getBlockId(i);
						for(k=0; k+1<m_part.size(); k++) {
							link = getBodyLinks(i, k);
							if(!link.empty()) part = getPart(ii, k);
							for(l=0; l<link.size(); l++) addBodyToMesh(mesh, part, link[l], i);
						}
					}
				}
			}
		}
	}
	mesh.setExternalBodies(getExternalBodies());

	// add local quads
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					for(j=0; j<m_part.size(); j++) {
						const PartMesh *part = getPart(ii, j);
						for(k=0; k<part->getQuadLocals(); k++) addQuadToMesh(mesh, part, k, i);
					}
				}
			}
		}
	}
}

Buffer< pair<uint,uint> > BlockMesh::getExternalNodes() const
{
	// send m_nbeg to neighbors
	Uint4 i;
	uint j, k;
	for(j=0; j<m_next.size(); j++) {
		for(i.t=m_next[j].i0.t; i.t<m_next[j].i1.t; i.t++) {
			for(i.z=m_next[j].i0.z; i.z<m_next[j].i1.z; i.z++) {
				for(i.y=m_next[j].i0.y; i.y<m_next[j].i1.y; i.y++) {
					for(i.x=m_next[j].i0.x; i.x<m_next[j].i1.x; i.x++) {
						const uint ii = getBlockId(i);
						sendMPI(&m_nbeg[ii], sizeof(uint), m_next[j].rank, 0);
					}
				}
			}
		}
	}

	// receive and compute external nodes
	uint extns = 0;
	Buffer< pair<uint,uint> > extn(m_nsize - m_nlocs);
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						uint beg;
						recvMPI(&beg, sizeof(uint), m_prev[j].rank, 0);
						Buffer<uint> link = getNodeLinks(i);
						for(k=0; k<link.size(); k++) extn[extns++] = pair<uint,uint>(m_prev[j].rank, beg + link[k]);
					}
				}
			}
		}
	}
	return extn;
}
Buffer< pair<uint,uint> > BlockMesh::getExternalEdges() const
{
	// send ebeg to neighbors
	Uint4 i;
	uint j, k, l;
	Buffer<uint> beg(m_part.size() - 1);
	for(j=0; j<m_next.size(); j++) {
		for(i.t=m_next[j].i0.t; i.t<m_next[j].i1.t; i.t++) {
			for(i.z=m_next[j].i0.z; i.z<m_next[j].i1.z; i.z++) {
				for(i.y=m_next[j].i0.y; i.y<m_next[j].i1.y; i.y++) {
					for(i.x=m_next[j].i0.x; i.x<m_next[j].i1.x; i.x++) {
						const uint ii = getBlockId(i);
						beg[0] = m_ebeg[ii];
						for(k=0; k+1<beg.size(); k++) {
							beg[k+1] = beg[k];
							if(i.x >= (k % 2) && i.y >= ((k / 2) % 2) && i.z >= ((k / 4) % 2) && i.t >= k / 8) beg[k+1] += getPart(ii, k)->getEdgeLocals();
						}
						sendMPI(&beg[0], beg.size() * sizeof(uint), m_next[j].rank, 0);
					}
				}
			}
		}
	}

	// receive and compute external edges
	uint extns = 0;
	Buffer< pair<uint,uint> > extn(m_esize - m_elocs);
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						recvMPI(&beg[0], beg.size() * sizeof(uint), m_prev[j].rank, 0);
						for(k=0; k<beg.size(); k++) {
							Buffer<uint> link = getEdgeLinks(i, k);
							for(l=0; l<link.size(); l++) extn[extns++] = pair<uint,uint>(m_prev[j].rank, beg[k] + link[l]);
						}
					}
				}
			}
		}
	}
	return extn;
}
Buffer< pair<uint,uint> > BlockMesh::getExternalFaces() const
{
	// send ebeg to neighbors
	Uint4 i;
	uint j, k, l;
	Buffer<uint> beg(m_part.size() - 1);
	for(j=0; j<m_next.size(); j++) {
		for(i.t=m_next[j].i0.t; i.t<m_next[j].i1.t; i.t++) {
			for(i.z=m_next[j].i0.z; i.z<m_next[j].i1.z; i.z++) {
				for(i.y=m_next[j].i0.y; i.y<m_next[j].i1.y; i.y++) {
					for(i.x=m_next[j].i0.x; i.x<m_next[j].i1.x; i.x++) {
						const uint ii = getBlockId(i);
						beg[0] = m_fbeg[ii];
						for(k=0; k+1<beg.size(); k++) {
							beg[k+1] = beg[k];
							if(i.x >= (k % 2) && i.y >= ((k / 2) % 2) && i.z >= ((k / 4) % 2) && i.t >= k / 8) beg[k+1] += getPart(ii, k)->getFaceLocals();
						}
						sendMPI(&beg[0], beg.size() * sizeof(uint), m_next[j].rank, 0);
					}
				}
			}
		}
	}

	// receive and compute external edges
	uint extns = 0;
	Buffer< pair<uint,uint> > extn(m_fsize - m_flocs);
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						recvMPI(&beg[0], beg.size() * sizeof(uint), m_prev[j].rank, 0);
						for(k=0; k<beg.size(); k++) {
							Buffer<uint> link = getFaceLinks(i, k);
							for(l=0; l<link.size(); l++) extn[extns++] = pair<uint,uint>(m_prev[j].rank, beg[k] + link[l]);
						}
					}
				}
			}
		}
	}
	return extn;
}
Buffer< pair<uint,uint> > BlockMesh::getExternalBodies() const
{
	// send ebeg to neighbors
	Uint4 i;
	uint j, k, l;
	Buffer<uint> beg(m_part.size() - 1);
	for(j=0; j<m_next.size(); j++) {
		for(i.t=m_next[j].i0.t; i.t<m_next[j].i1.t; i.t++) {
			for(i.z=m_next[j].i0.z; i.z<m_next[j].i1.z; i.z++) {
				for(i.y=m_next[j].i0.y; i.y<m_next[j].i1.y; i.y++) {
					for(i.x=m_next[j].i0.x; i.x<m_next[j].i1.x; i.x++) {
						const uint ii = getBlockId(i);
						beg[0] = m_bbeg[ii];
						for(k=0; k+1<beg.size(); k++) {
							beg[k+1] = beg[k];
							if(i.x >= (k % 2) && i.y >= ((k / 2) % 2) && i.z >= ((k / 4) % 2) && i.t >= k / 8) beg[k+1] += getPart(ii, k)->getBodyLocals();
						}
						sendMPI(&beg[0], beg.size() * sizeof(uint), m_next[j].rank, 0);
					}
				}
			}
		}
	}

	// receive and compute external edges
	uint extns = 0;
	Buffer< pair<uint,uint> > extn(m_bsize - m_blocs);
	for(j=0; j<m_prev.size(); j++) {
		for(i.t=m_prev[j].i0.t; i.t<m_prev[j].i1.t; i.t++) {
			for(i.z=m_prev[j].i0.z; i.z<m_prev[j].i1.z; i.z++) {
				for(i.y=m_prev[j].i0.y; i.y<m_prev[j].i1.y; i.y++) {
					for(i.x=m_prev[j].i0.x; i.x<m_prev[j].i1.x; i.x++) {
						recvMPI(&beg[0], beg.size() * sizeof(uint), m_prev[j].rank, 0);
						for(k=0; k<beg.size(); k++) {
							Buffer<uint> link = getBodyLinks(i, k);
							for(l=0; l<link.size(); l++) extn[extns++] = pair<uint,uint>(m_prev[j].rank, beg[k] + link[l]);
						}
					}
				}
			}
		}
	}
	return extn;
}
Buffer< pair<uint,uint> > BlockMesh::getMyExternals(const Buffer< pair<uint,uint> > &ext) const
{
	uint i, j, num;
	for(j=0; j<m_prev.size(); j++)
	{
		num = 0;
		for(i=0; i<ext.size(); i++)
		{
			if(ext[i].first == m_prev[j].rank) num++;
		}
		sendMPI(&num, sizeof(uint), m_prev[j].rank, 0);
		if(num == 0) continue;

		Buffer<uint> link(num);
		num = 0;
		for(i=0; i<ext.size(); i++)
		{
			if(ext[i].first == m_prev[j].rank) link[num++] = ext[i].second;
		}
		sendMPI(&link[0], num * sizeof(uint), m_prev[j].rank, 1);
	}

	num = 0;
	Buffer<uint> nnum(m_next.size());
	for(j=0; j<m_next.size(); j++)
	{
		recvMPI(&nnum[j], sizeof(uint), m_next[j].rank, 0);
		num += nnum[j];
	}
	Buffer< pair<uint,uint> > mext(num);
	num = 0;
	for(j=0; j<m_next.size(); j++)
	{
		if(nnum[j] == 0) continue;
		Buffer<uint> link(nnum[j]);
		recvMPI(&link[0], nnum[j] * sizeof(uint), m_next[j].rank, 1);
		for(i=0; i<link.size(); i++) mext[num++] = pair<uint,uint>(m_next[j].rank, link[i]);
	}
	return mext;
}


Buffer<uint> BlockMesh::getNodeLinks(const Uint4 &i) const
{
	uint links = 0;
	Buffer<uint> link;
	const uint ii = getBlockId(i);
	Uint4 j, k;
	const Uint4 k2(1,2,4,8);
	const Uint4 k23(2,6,18,54);
	const Uint4 s0(uint(i.x<m_0.x), uint(i.y<m_0.y), uint(i.z<m_0.z), uint(i.t<m_0.t));
	const Uint4 s1(uint(i.x+1<m_s.x), uint(i.y+1<m_s.y), uint(i.z+1<m_s.z), uint(i.t+1<m_s.t));
	for(j.t=s0.t; j.t<=s1.t; j.t++) {
		for(j.z=s0.z; j.z<=s1.z; j.z++) {
			for(j.y=s0.y; j.y<=s1.y; j.y++) {
				for(j.x=s0.x; j.x<=s1.x; j.x++) {
					const uint jj = ii + getBlockId(j);
					const uint part = k23.dot(j);
					for(k.t=j.t; k.t<=m_0.t; k.t++) {
						for(k.z=j.z; k.z<=m_0.z; k.z++) {
							for(k.y=j.y; k.y<=m_0.y; k.y++) {
								for(k.x=j.x; k.x<=m_0.x; k.x++) {
									const PartMesh *mesh = getPart(jj, k2.dot(k));
									for(uint l=mesh->getNodeLocals(); l<mesh->getNodeSize(); l++) {
										if(mesh->getNodePart(l) == part) link.gatherOnce(mesh->getNodeLink(l), links);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	link.resize(links);
	return link;
}

Buffer<uint> BlockMesh::getEdgeLinks(const Uint4 &i, const uint part) const
{
	uint links = 0;
	Buffer<uint> link;
	const uint ii = getBlockId(i);
	Uint4 j, k;
	const Uint4 k2(1,2,4,8);
	const Uint4 k3(1,3,9,27);
	const Uint4 k23(2,6,18,54);
	const Uint4 p(part % 2, (part / 2) % 2, (part / 4) % 2, part / 8);
	const Uint4 s0(uint(i.x<m_0.x), uint(i.y<m_0.y), uint(i.z<m_0.z), uint(i.t<m_0.t));
	const Uint4 s1(uint(p.x==0&&i.x+1<m_s.x), uint(p.y==0&&i.y+1<m_s.y), uint(p.z==0&&i.z+1<m_s.z), uint(p.t==0&&i.t+1<m_s.t));
	for(j.t=s0.t; j.t<=s1.t; j.t++) {
		for(j.z=s0.z; j.z<=s1.z; j.z++) {
			for(j.y=s0.y; j.y<=s1.y; j.y++) {
				for(j.x=s0.x; j.x<=s1.x; j.x++) {
					const uint jj = ii + getBlockId(j);
					const uint jpart = k23.dot(j) + k3.dot(p);
					for(k.t=j.t+p.t; k.t<=m_0.t; k.t++) {
						for(k.z=j.z+p.z; k.z<=m_0.z; k.z++) {
							for(k.y=j.y+p.y; k.y<=m_0.y; k.y++) {
								for(k.x=j.x+p.x; k.x<=m_0.x; k.x++) {
									const PartMesh *mesh = getPart(jj, k2.dot(k));
									for(uint l=mesh->getEdgeLocals(); l<mesh->getEdgeSize(); l++) {
										if(mesh->getEdgePart(l) == jpart) link.gatherOnce(mesh->getEdgeLink(l), links);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	link.resize(links);
	return link;
}

Buffer<uint> BlockMesh::getFaceLinks(const Uint4 &i, const uint part) const
{
	uint links = 0;
	Buffer<uint> link;
	const uint ii = getBlockId(i);
	Uint4 j, k;
	const Uint4 k2(1,2,4,8);
	const Uint4 k3(1,3,9,27);
	const Uint4 k23(2,6,18,54);
	const Uint4 p(part % 2, (part / 2) % 2, (part / 4) % 2, part / 8);
	const Uint4 s0(uint(i.x<m_0.x), uint(i.y<m_0.y), uint(i.z<m_0.z), uint(i.t<m_0.t));
	const Uint4 s1(uint(p.x==0&&i.x+1<m_s.x), uint(p.y==0&&i.y+1<m_s.y), uint(p.z==0&&i.z+1<m_s.z), uint(p.t==0&&i.t+1<m_s.t));
	for(j.t=s0.t; j.t<=s1.t; j.t++) {
		for(j.z=s0.z; j.z<=s1.z; j.z++) {
			for(j.y=s0.y; j.y<=s1.y; j.y++) {
				for(j.x=s0.x; j.x<=s1.x; j.x++) {
					const uint jj = ii + getBlockId(j);
					const uint jpart = k23.dot(j) + k3.dot(p);
					for(k.t=j.t+p.t; k.t<=m_0.t; k.t++) {
						for(k.z=j.z+p.z; k.z<=m_0.z; k.z++) {
							for(k.y=j.y+p.y; k.y<=m_0.y; k.y++) {
								for(k.x=j.x+p.x; k.x<=m_0.x; k.x++) {
									const PartMesh *mesh = getPart(jj, k2.dot(k));
									for(uint l=mesh->getFaceLocals(); l<mesh->getFaceSize(); l++) {
										if(mesh->getFacePart(l) == jpart) link.gatherOnce(mesh->getFaceLink(l), links);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	link.resize(links);
	return link;
}

Buffer<uint> BlockMesh::getBodyLinks(const Uint4 &i, const uint part) const
{
	uint links = 0;
	Buffer<uint> link;
	const uint ii = getBlockId(i);
	Uint4 j, k;
	const Uint4 k2(1,2,4,8);
	const Uint4 k3(1,3,9,27);
	const Uint4 k23(2,6,18,54);
	const Uint4 p(part % 2, (part / 2) % 2, (part / 4) % 2, part / 8);
	const Uint4 s0(uint(i.x<m_0.x), uint(i.y<m_0.y), uint(i.z<m_0.z), uint(i.t<m_0.t));
	const Uint4 s1(uint(p.x==0&&i.x+1<m_s.x), uint(p.y==0&&i.y+1<m_s.y), uint(p.z==0&&i.z+1<m_s.z), uint(p.t==0&&i.t+1<m_s.t));
	for(j.t=s0.t; j.t<=s1.t; j.t++) {
		for(j.z=s0.z; j.z<=s1.z; j.z++) {
			for(j.y=s0.y; j.y<=s1.y; j.y++) {
				for(j.x=s0.x; j.x<=s1.x; j.x++) {
					const uint jj = ii + getBlockId(j);
					const uint jpart = k23.dot(j) + k3.dot(p);
					for(k.t=j.t+p.t; k.t<=m_0.t; k.t++) {
						for(k.z=j.z+p.z; k.z<=m_0.z; k.z++) {
							for(k.y=j.y+p.y; k.y<=m_0.y; k.y++) {
								for(k.x=j.x+p.x; k.x<=m_0.x; k.x++) {
									const PartMesh *mesh = getPart(jj, k2.dot(k));
									for(uint l=mesh->getBodyLocals(); l<mesh->getBodySize(); l++) {
										if(mesh->getBodyPart(l) == jpart) link.gatherOnce(mesh->getBodyLink(l), links);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	link.resize(links);
	return link;
}

void BlockMesh::addNodeToMesh(Mesh &mesh, const PartMesh *part, const uint n, const Uint4 &i) const
{
	const Vector4 p = Vector4(m_p.x + i.x * m_d.x, m_p.y + i.y * m_d.y, m_p.z + i.z * m_d.z, m_p.t + i.t * m_d.t) + part->getNodePosition(n);
	const uint node = mesh.addNode(p);
	mesh.setNodeWeight(node, part->getNodeWeight(n));
}
void BlockMesh::addEdgeToMesh(Mesh &mesh, const PartMesh *part, const uint e, const Uint4 &i) const
{
	Buffer<uint> n = part->getEdgeNodes(e);
	n[0] = getNode(part, n[0], i);
	n[1] = getNode(part, n[1], i);
	mesh.addEdge(n[0], n[1]);
}
void BlockMesh::addFaceToMesh(Mesh &mesh, const PartMesh *part, const uint f, const Uint4 &i) const
{
	Buffer<uint> e = part->getFaceEdges(f);
	for(uint j=0; j<e.size(); j++) e[j] = getEdge(part, e[j], i);
	mesh.addFace(e);
}
void BlockMesh::addBodyToMesh(Mesh &mesh, const PartMesh *part, const uint b, const Uint4 &i) const
{
	Buffer<uint> f = part->getBodyFaces(b);
	for(uint j=0; j<f.size(); j++) f[j] = getFace(part, f[j], i);
	mesh.addBody(f);
}
void BlockMesh::addQuadToMesh(Mesh &mesh, const PartMesh *part, const uint q, const Uint4 &i) const
{
	Buffer<uint> b = part->getQuadBodies(q);
	for(uint j=0; j<b.size(); j++) b[j] = getBody(part, b[j], i);
	mesh.addQuad(b);
}

void BlockMesh::createFormIntegrator(const FormGrade grade, map<Type, BlockIntegrator*> &blocks, const int num) const
{
	PartMesh mesh(0, PARTS, m_dim);
	Buffer<const PartMesh *> part(m_part.size());
	uint j;
	Uint4 i;
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					const Type type = getType(ii);
					map<Type, BlockIntegrator*>::iterator it = blocks.find(type);
					if(it != blocks.end()) continue;

					for(j=0; j<m_part.size(); j++) part[j] = getPart(ii, j);
					mesh.createCombined(part);

					BlockIntegrator *block = new BlockIntegrator();
					block->init(grade, mesh, num);
					blocks.insert(make_pair(type, block));
				}
			}
		}
	}
}
void BlockMesh::createWedgeIntegrator(const FormGrade grade, map<Type, BlockIntegrator*> &blocks, const int num) const
{
	PartMesh mesh(0, PARTS, m_dim);
	Buffer<const PartMesh *> part(m_part.size());
	uint j;
	Uint4 i;
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					const Type type = getType(ii);
					map<Type, BlockIntegrator*>::iterator it = blocks.find(type);
					if(it != blocks.end()) continue;

					for(j=0; j<m_part.size(); j++) part[j] = getPart(ii, j);
					mesh.createCombined(part);

					BlockIntegrator *block = new BlockIntegrator();
					block->initWedge(grade, mesh, num);
					blocks.insert(make_pair(type, block));
				}
			}
		}
	}
}
void BlockMesh::createInterpolator(const FormGrade grade)
{
	if(!m_poly[grade].empty()) return;

	PartMesh mesh(0, PARTS, m_dim);
	Buffer<const PartMesh *> part(m_part.size());
	uint j;
	Uint4 i;
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					const Type type = getType(ii);
					map<Type, BlockInterpolator*>::iterator it = m_poly[grade].find(type);
					if(it != m_poly[grade].end()) continue;

					for(j=0; j<m_part.size(); j++) part[j] = getPart(ii, j);
					mesh.createCombined(part);

					BlockInterpolator *poly = new BlockInterpolator();
					poly->init(grade, mesh, m_d, 2);
					m_poly[grade].insert(std::make_pair(type, poly));
				}
			}
		}
	}
}
/*
uint BlockMesh::integrateVectors(const FormGrade grade, Buffer<double> &result) const
{
	// initialize result
	const uint fields = FormGradeVectorDimension(grade, m_dim);
	const uint gdim = FormGradeDimension(grade);
	const uint locs = getLocals(gdim);
	if(result.size() != fields * locs) {
		result.resize(fields * locs);
		result.fill(0.0);
	}

	// integrate over mesh elements
	map<Type, BlockIntegrator*> blocks;
	createFormIntegrator(grade, blocks, 0);

	// reserve space for external terms
	const bool dual = FormGradeIsDual(grade);
	Buffer< pair<uint,uint> > ext;
	if(dual) ext = getExternals(gdim);
	Buffer<double> add(fields * ext.size(), 0.0);

	// use blocks to integrate over mesh elements
	Uint4 i;
	uint j, k;
	Buffer<double *> term;
	for(i.t=m_0.t; i.t<m_s.t; i.t++) {
		for(i.z=m_0.z; i.z<m_s.z; i.z++) {
			for(i.y=m_0.y; i.y<m_s.y; i.y++) {
				for(i.x=m_0.x; i.x<m_s.x; i.x++) {
					const uint ii = getBlockId(i);
					const BlockIntegrator *block = blocks.find(getType(ii))->second;
					if(dual)
					{
						const Buffer< pair<uint,uint> > &pext = block->getExternals();
						if(term.size() < pext.size()) term.resize(pext.size());
						for(j=0; j<pext.size(); j++)
						{
							const uint cell = getCell(gdim, pext[j], i);
							term[j] = (cell < locs ? &result[fields * cell] : &add[fields * (cell - locs)]);
						}
					}
					block->integrateVectors(&result[fields * getBegin(gdim, ii)], fields, term);
//						block->integrate(func, getPosition(i), &result[getBegin(gdim, ii)], term);
				}
			}
		}
	}
	clearMap(blocks);

	// communicate with external terms
	if(dual) {
		const Buffer< pair<uint,uint> > mext = getMyExternals(ext);
		for(j=0; j<ext.size(); j++) {
			sendMPI(&add[fields * j], fields * sizeof(double), ext[j].first, 0);
		}
		Buffer<double> addition(fields);
		for(j=0; j<mext.size(); j++) {
			recvMPI(&addition[0], fields * sizeof(double), mext[j].first, 0);
			const uint jj = fields * mext[j].second;
			for(k=0; k<fields; k++) result[jj + k] += addition[k];
		}
	}
	return fields;
}
*/
Sparse<sign> &BlockMesh::integrateDerivative(const FormGrade grade, Sparse<sign> &d) const
{
	uint j, k, l;
	Uint4 i;
	if(grade == fg_prim0 || grade == fg_dual1)
	{
		Buffer< Buffer< pair<uint,sign> > > inc(m_elocs);
		l = 0;
		for(i.t=m_0.t; i.t<m_s.t; i.t++) {
			for(i.z=m_0.z; i.z<m_s.z; i.z++) {
				for(i.y=m_0.y; i.y<m_s.y; i.y++) {
					for(i.x=m_0.x; i.x<m_s.x; i.x++) {
						const uint ii = getBlockId(i);
						for(j=0; j<m_part.size(); j++) {
							const PartMesh *part = getPart(ii, j);
							for(k=0; k<part->getEdgeLocals(); k++, l++) addEdgeIncidences(part, k, i, inc[l]);
						}
					}
				}
			}
		}
		d.setFull(m_nlocs, inc, getExternalNodes(), getNextRanks(), getPrevRanks());
		if(grade == fg_dual1) d.setTranspose(d);
		return d;
	}
	if(grade == fg_prim1 || grade == fg_dual2)
	{
		Buffer< Buffer< pair<uint,sign> > > inc(m_flocs);
		l = 0;
		for(i.t=m_0.t; i.t<m_s.t; i.t++) {
			for(i.z=m_0.z; i.z<m_s.z; i.z++) {
				for(i.y=m_0.y; i.y<m_s.y; i.y++) {
					for(i.x=m_0.x; i.x<m_s.x; i.x++) {
						const uint ii = getBlockId(i);
						for(j=0; j<m_part.size(); j++) {
							const PartMesh *part = getPart(ii, j);
							for(k=0; k<part->getFaceLocals(); k++, l++) addFaceIncidences(part, k, i, inc[l]);
						}
					}
				}
			}
		}
		d.setFull(m_elocs, inc, getExternalEdges(), getNextRanks(), getPrevRanks());
		if(grade == fg_dual2) d.setTranspose(d);
		return d;
	}
	if(grade == fg_prim2 || grade == fg_dual3)
	{
		Buffer< Buffer< pair<uint,sign> > > inc(m_blocs);
		l = 0;
		for(i.t=m_0.t; i.t<m_s.t; i.t++) {
			for(i.z=m_0.z; i.z<m_s.z; i.z++) {
				for(i.y=m_0.y; i.y<m_s.y; i.y++) {
					for(i.x=m_0.x; i.x<m_s.x; i.x++) {
						const uint ii = getBlockId(i);
						for(j=0; j<m_part.size(); j++) {
							const PartMesh *part = getPart(ii, j);
							for(k=0; k<part->getBodyLocals(); k++, l++) addBodyIncidences(part, k, i, inc[l]);
						}
					}
				}
			}
		}
		d.setFull(m_flocs, inc, getExternalFaces(), getNextRanks(), getPrevRanks());
		if(grade == fg_dual3) d.setTranspose(d);
		return d;
	}
	if(grade == fg_prim3 || grade == fg_dual4)
	{
		Buffer< Buffer< pair<uint,sign> > > inc(m_qlocs);
		l = 0;
		for(i.t=m_0.t; i.t<m_s.t; i.t++) {
			for(i.z=m_0.z; i.z<m_s.z; i.z++) {
				for(i.y=m_0.y; i.y<m_s.y; i.y++) {
					for(i.x=m_0.x; i.x<m_s.x; i.x++) {
						const uint ii = getBlockId(i);
						for(j=0; j<m_part.size(); j++) {
							const PartMesh *part = getPart(ii, j);
							for(k=0; k<part->getQuadLocals(); k++, l++) addQuadIncidences(part, k, i, inc[l]);
						}
					}
				}
			}
		}
		d.setFull(m_blocs, inc, getExternalBodies(), getNextRanks(), getPrevRanks());
		if(grade == fg_dual4) d.setTranspose(d);
		return d;
	}
	cout << "BlockMesh::integrate(Derivative &) -> Cannot integrate derivative due to unsuitable grade." << endl;
	return d;
}

void BlockMesh::addEdgeIncidences(const PartMesh *part, const uint e, const Uint4 &i, Buffer< pair<uint,sign> > &inc) const
{
	const Buffer<uint> &n = part->getEdgeNodes(e);
	inc.resize(n.size());
	for(uint j=0; j<n.size(); j++)
	{
		const uint nj = getNode(part, n[j], i);
		inc[j] = pair<uint,sign>(nj, sign(part->getEdgeIncidence(e, n[j])));
	}
}
void BlockMesh::addFaceIncidences(const PartMesh *part, const uint f, const Uint4 &i, Buffer< pair<uint,sign> > &inc) const
{
	const Buffer<uint> &e = part->getFaceEdges(f);
	inc.resize(e.size());
	for(uint j=0; j<e.size(); j++)
	{
		const uint ej = getEdge(part, e[j], i);
		inc[j] = pair<uint,sign>(ej, sign(part->getFaceIncidence(f, e[j])));
	}
}
void BlockMesh::addBodyIncidences(const PartMesh *part, const uint b, const Uint4 &i, Buffer< pair<uint,sign> > &inc) const
{
	const Buffer<uint> &f = part->getBodyFaces(b);
	inc.resize(f.size());
	for(uint j=0; j<f.size(); j++)
	{
		const uint fj = getFace(part, f[j], i);
		inc[j] = pair<uint,sign>(fj, sign(part->getBodyIncidence(b, f[j])));
	}
}
void BlockMesh::addQuadIncidences(const PartMesh *part, const uint q, const Uint4 &i, Buffer< pair<uint,sign> > &inc) const
{
	const Buffer<uint> &b = part->getQuadBodies(q);
	inc.resize(b.size());
	for(uint j=0; j<b.size(); j++)
	{
		const uint bj = getBody(part, b[j], i);
		inc[j] = pair<uint,sign>(bj, sign(part->getQuadIncidence(q, b[j])));
	}
}

uint BlockMesh::getNode(const PartMesh *part, const uint n, const Uint4 &i) const
{
	return getNode(pair<uint,uint>(part->getNodePart(n), part->getNodeLink(n)), i);
}
uint BlockMesh::getNode(const pair<uint,uint> &ext, const Uint4 &i) const
{
	const Uint4 j(i.x - (ext.first % 3) / 2, i.y - (ext.first % 9) / 6, i.z - (ext.first % 27) / 18, i.t - (ext.first % 81) / 54);
	if(j.x < m_0.x || j.y < m_0.y || j.z < m_0.z || j.t < m_0.t) return m_nbeg[getBlockId(j)] + getNodeLinks(j).findFirst(ext.second);
	return m_nbeg[getBlockId(j)] + ext.second;
}
uint BlockMesh::getEdge(const PartMesh *part, const uint e, const Uint4 &i) const
{
	return getEdge(pair<uint,uint>(part->getEdgePart(e), part->getEdgeLink(e)), i);
}
uint BlockMesh::getEdge(const pair<uint,uint> &ext, const Uint4 &i) const
{
	const Uint4 j(i.x - (ext.first % 3) / 2, i.y - (ext.first % 9) / 6, i.z - (ext.first % 27) / 18, i.t - (ext.first % 81) / 54);
	const uint jp = ((ext.first % 3) % 2) + 2 * (((ext.first / 3) % 3) % 2) + 4 * (((ext.first / 9) % 3) % 2) + 8 * (((ext.first / 27) % 3) % 2);
	const uint jj = getBlockId(j);
	if(j.x < m_0.x || j.y < m_0.y || j.z < m_0.z || j.t < m_0.t) {
		uint res = m_ebeg[jj] + getEdgeLinks(j, jp).findFirst(ext.second);
		for(uint k=0; k<jp; k++) res += getEdgeLinks(j, k).size();
		return res;
	}
	uint res = m_ebeg[jj] + ext.second;
	for(uint k=0; k<jp; k++) res += getPart(jj, k)->getEdgeLocals();
	return res;
}
uint BlockMesh::getFace(const PartMesh *part, const uint f, const Uint4 &i) const
{
	return getFace(pair<uint,uint>(part->getFacePart(f), part->getFaceLink(f)), i);
}
uint BlockMesh::getFace(const pair<uint,uint> &ext, const Uint4 &i) const
{
	const Uint4 j(i.x - (ext.first % 3) / 2, i.y - (ext.first % 9) / 6, i.z - (ext.first % 27) / 18, i.t - (ext.first % 81) / 54);
	const uint jp = ((ext.first % 3) % 2) + 2 * (((ext.first / 3) % 3) % 2) + 4 * (((ext.first / 9) % 3) % 2) + 8 * (((ext.first / 27) % 3) % 2);
	const uint jj = getBlockId(j);
	if(j.x < m_0.x || j.y < m_0.y || j.z < m_0.z || j.t < m_0.t) {
		uint res = m_fbeg[jj] + getFaceLinks(j, jp).findFirst(ext.second);
		for(uint k=0; k<jp; k++) res += getFaceLinks(j, k).size();
		return res;
	}
	uint res = m_fbeg[jj] + ext.second;
	for(uint k=0; k<jp; k++) res += getPart(jj, k)->getFaceLocals();
	return res;
}
uint BlockMesh::getBody(const PartMesh *part, const uint b, const Uint4 &i) const
{
	return getBody(pair<uint,uint>(part->getBodyPart(b), part->getBodyLink(b)), i);
}
uint BlockMesh::getBody(const pair<uint,uint> &ext, const Uint4 &i) const
{
	const Uint4 j(i.x - (ext.first % 3) / 2, i.y - (ext.first % 9) / 6, i.z - (ext.first % 27) / 18, i.t - (ext.first % 81) / 54);
	const uint jp = ((ext.first % 3) % 2) + 2 * (((ext.first / 3) % 3) % 2) + 4 * (((ext.first / 9) % 3) % 2) + 8 * (((ext.first / 27) % 3) % 2);
	const uint jj = getBlockId(j);
	if(j.x < m_0.x || j.y < m_0.y || j.z < m_0.z || j.t < m_0.t) {
		uint res = m_bbeg[jj] + getBodyLinks(j, jp).findFirst(ext.second);
		for(uint k=0; k<jp; k++) res += getBodyLinks(j, k).size();
		return res;
	}
	uint res = m_bbeg[jj] + ext.second;
	for(uint k=0; k<jp; k++) res += getPart(jj, k)->getBodyLocals();
	return res;
}

Buffer<uint> BlockMesh::getPrevRanks() const
{
	Buffer<uint> buf(m_prev.size());
	for(uint i=0; i<buf.size(); i++) buf[i] = m_prev[i].rank;
	return buf;
}
Buffer<uint> BlockMesh::getNextRanks() const
{
	Buffer<uint> buf(m_next.size());
	for(uint i=0; i<buf.size(); i++) buf[i] = m_next[i].rank;
	return buf;
}

