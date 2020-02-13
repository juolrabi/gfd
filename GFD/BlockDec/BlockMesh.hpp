/*
BlockMesh is a mesh consisting of blocks of cells.
*/

#ifndef _BLOCKMESH_HPP_INCLUDED_
#define _BLOCKMESH_HPP_INCLUDED_

#include "BlockInterpolator.hpp"
#include "../Mesh/BuilderMesh.hpp"
#include "../Mesh/PartMesh.hpp"
#include "../Types/Types.hpp"
#include "../Types/Buffer.hpp"
#include <map>
//#include <cstdint>
#include "../Types/MpiEasy.hpp"

using namespace std;

namespace gfd
{

class Type
{
public:
	Type(const uchar a){
		c0 = ullong(a);
		c1 = 0;
	}
	Type(const uchar a, const uchar a1) {
		c0 = ullong(a) | (ullong(a1) << 8);
		c1 = 0;
	}
	Type(const uchar a, const uchar a1, const uchar a10, const uchar a11) {
		c0 = ullong(a) | (ullong(a1) << 8) | (ullong(a10) << 16) | (ullong(a11) << 24);
		c1 = 0;
	}
	Type(const uchar a, const uchar a1, const uchar a10, const uchar a11,
		const uchar a100, const uchar a101, const uchar a110, const uchar a111) {
		c0 = ullong(a) | (ullong(a1) << 8) | (ullong(a10) << 16) | (ullong(a11) << 24) |
			(ullong(a100) << 32) | (ullong(a101) << 40) | (ullong(a110) << 48) | (ullong(a111) << 56);
		c1 = 0;
	}
	Type(const uchar a, const uchar a1, const uchar a10, const uchar a11,
		const uchar a100, const uchar a101, const uchar a110, const uchar a111,
		const uchar a1000, const uchar a1001, const uchar a1010, const uchar a1011,
		const uchar a1100, const uchar a1101, const uchar a1110, const uchar a1111) {
		c0 = ullong(a) | (ullong(a1) << 8) | (ullong(a10) << 16) | (ullong(a11) << 24) |
			(ullong(a100) << 32) | (ullong(a101) << 40) | (ullong(a110) << 48) | (ullong(a111) << 56);
		c1 = ullong(a1000) | (ullong(a1001) << 8) | (ullong(a1010) << 16) | (ullong(a1011) << 24) |
			(ullong(a1100) << 32) | (ullong(a1101) << 40) | (ullong(a1110) << 48) | (ullong(a1111) << 56);
	}
	bool operator<(const Type &t) const { return c0<t.c0 || (c0==t.c0 && c1<t.c1); }
	uchar get() const { return uchar(c0); }
	uchar get1() const { return uchar(c0 >> 8); }
	uchar get10() const { return uchar(c0 >> 16); }
	uchar get11() const { return uchar(c0 >> 24); }
	uchar get100() const { return uchar(c0 >> 32); }
	uchar get101() const { return uchar(c0 >> 40); }
	uchar get110() const { return uchar(c0 >> 48); }
	uchar get111() const { return uchar(c0 >> 56); }
	uchar get1000() const { return uchar(c1); }
	uchar get1001() const { return uchar(c1 >> 8); }
	uchar get1010() const { return uchar(c1 >> 16); }
	uchar get1011() const { return uchar(c1 >> 24); }
	uchar get1100() const { return uchar(c1 >> 32); }
	uchar get1101() const { return uchar(c1 >> 40); }
	uchar get1110() const { return uchar(c1 >> 48); }
	uchar get1111() const { return uchar(c1 >> 56); }
private:
	ullong c0;
	ullong c1;
};

class BlockMesh
{
public:
	BlockMesh();
	virtual ~BlockMesh() { clear(); }
	void clear();

	void init(const double p, const double d, const uint sx, const uchar *type);
	void init(const Vector2 &p, const Vector2 &d, const uint sx, const uint sy, const uchar *type);
	void init(const Vector3 &p, const Vector3 &d, const uint sx, const uint sy, const uint sz, const uchar *type);
	void init(const Vector4 &p, const Vector4 &d, const uint sx, const uint sy, const uint sz, const uint st, const uchar *type);

	void toMesh(PartMesh &mesh) const;
	Sparse<sign> &integrateDerivative(const FormGrade grade, Sparse<sign> &d) const;
	template<typename T> Diagonal<T> &integrateForm(T func(const Buffer<double> &), const int num, const FormGrade grade, Diagonal<T> &result) const {
		// initialize form
		const uint gdim = FormGradeDimension(grade);
		const uint locs = getLocals(gdim);
		if(!result.m_full || result.m_height != locs) result.setFullOfZeros(locs);

		// integrate over mesh elements
		map<Type, BlockIntegrator*> blocks;
		createFormIntegrator(grade, blocks, num);
		integrate(func, blocks, gdim, FormGradeIsDual(grade), result);
		clearMap(blocks);
		return result;
	}
	template<typename T> Diagonal<T> &integrateWedge(T func(const Buffer<double> &), const int num, const FormGrade grade, Diagonal<T> &result) const {
		// initialize hodge
		const uint gdim = FormGradeDimension(grade);
		const uint locs = getLocals(gdim);
		if(!result.m_full || result.m_height != locs) result.setFullOfZeros(locs);

		// integrate over mesh elements
		map<Type, BlockIntegrator*> blocks;
		createWedgeIntegrator(grade, blocks, num);
		integrate(func, blocks, gdim, true, result);
		clearMap(blocks);
		return result;
	}
	template<typename T> Diagonal<T> &integrateHodge(T func(const Buffer<double> &), const int num, const FormGrade grade, Diagonal<T> &result) const {
		integrateWedge(func, num, grade, result);
		switch(FormGradeVectorDimension(grade, m_dim)) {
		case 1: return divideByVectorSquare(FormTwoVector2, TwoVector2(0), grade, result);
		case 2: return divideByVectorSquare(FormVector2, Vector2(0,0), grade, result);
		case 3: return divideByVectorSquare(FormVector3, Vector3(0,0,0), grade, result);
		case 4: return divideByVectorSquare(FormVector4, Vector4(0,0,0,0), grade, result);
		default: return divideByVectorSquare(FormTwoVector4, TwoVector4(0,0,0,0,0,0), grade, result);
		}
	}
	template<typename T> Buffer<T> interpolateForm(const FormGrade grade, const Column<T> &form, const Vector4 &p) { return interpolateForm(form, p, Buffer<T>()); }
	template<typename T> Buffer<T> &interpolateForm(const FormGrade grade, const Column<T> &form, const Vector4 &p, Buffer<T> &result) {
		// check that m_poly is initialized
		createInterpolator(grade);

		// initialize result
		result.resize(FormGradeVectorDimension(grade, m_dim));
		result.fill(form.m_zero);

		// find block
		const Vector4 r = (p - m_p + m_d) / m_d;
		if(r.x < 0.0 || r.y < 0.0 || r.z < 0.0 || r.t < 0.0) return result;
		const Uint4 i(uint(r.x), uint(r.y), uint(r.z), uint(r.t));
		if(i.x < m_0.x || i.y < m_0.y || i.z < m_0.z || i.t < m_0.t) return result;
		if(i.x > m_s.x || i.y > m_s.y || i.z > m_s.z || i.t > m_s.t) return result;

		// interpolate from block corners
		if(i.t < m_s.t) {
			if(i.z < m_s.z) {
				if(i.y < m_s.y) {
					if(i.x < m_s.x) sumInterpolation(grade, form, p, Uint4(i.x, i.y, i.z, i.t), &result[0]);
					if(i.x > m_0.x) sumInterpolation(grade, form, p, Uint4(i.x-1, i.y, i.z, i.t), &result[0]);
				}
				if(i.y > m_0.y) {
					if(i.x < m_s.x) sumInterpolation(grade, form, p, Uint4(i.x, i.y-1, i.z, i.t), &result[0]);
					if(i.x > m_0.x) sumInterpolation(grade, form, p, Uint4(i.x-1, i.y-1, i.z, i.t), &result[0]);
				}
			}
			if(i.z > m_0.z) {
				if(i.y < m_s.y) {
					if(i.x < m_s.x) sumInterpolation(grade, form, p, Uint4(i.x, i.y, i.z-1, i.t), &result[0]);
					if(i.x > m_0.x) sumInterpolation(grade, form, p, Uint4(i.x-1, i.y, i.z-1, i.t), &result[0]);
				}
				if(i.y > m_0.y) {
					if(i.x < m_s.x) sumInterpolation(grade, form, p, Uint4(i.x, i.y-1, i.z-1, i.t), &result[0]);
					if(i.x > m_0.x) sumInterpolation(grade, form, p, Uint4(i.x-1, i.y-1, i.z-1, i.t), &result[0]);
				}
			}
		}
		if(i.t > m_0.t) {
			if(i.z < m_s.z) {
				if(i.y < m_s.y) {
					if(i.x < m_s.x) sumInterpolation(grade, form, p, Uint4(i.x, i.y, i.z, i.t-1), &result[0]);
					if(i.x > m_0.x) sumInterpolation(grade, form, p, Uint4(i.x-1, i.y, i.z, i.t-1), &result[0]);
				}
				if(i.y > m_0.y) {
					if(i.x < m_s.x) sumInterpolation(grade, form, p, Uint4(i.x, i.y-1, i.z, i.t-1), &result[0]);
					if(i.x > m_0.x) sumInterpolation(grade, form, p, Uint4(i.x-1, i.y-1, i.z, i.t-1), &result[0]);
				}
			}
			if(i.z > m_0.z) {
				if(i.y < m_s.y) {
					if(i.x < m_s.x) sumInterpolation(grade, form, p, Uint4(i.x, i.y, i.z-1, i.t-1), &result[0]);
					if(i.x > m_0.x) sumInterpolation(grade, form, p, Uint4(i.x-1, i.y, i.z-1, i.t-1), &result[0]);
				}
				if(i.y > m_0.y) {
					if(i.x < m_s.x) sumInterpolation(grade, form, p, Uint4(i.x, i.y-1, i.z-1, i.t-1), &result[0]);
					if(i.x > m_0.x) sumInterpolation(grade, form, p, Uint4(i.x-1, i.y-1, i.z-1, i.t-1), &result[0]);
				}
			}
		}
		return result;
	}


protected:

	uint m_dim;
	Buffer< map<Type, PartMesh*> > m_part;

	struct RankRectangle
	{
		uint rank;
		Uint4 i0;
		Uint4 i1;
		RankRectangle() {}
		RankRectangle(const uint arank, const Uint4 &ai0, const Uint4 &ai1) {
			rank = arank;
			i0 = ai0;
			i1 = ai1;
		}
	};
	Buffer<RankRectangle> m_prev; // previous neighbor ranks and their areas
	Buffer<RankRectangle> m_next; // next neighbor ranks and their areas

	Uint4 m_s; // number of blocks in each dimension
	Uint4 m_0; // begin block in each dimension (0 or 1)
	Buffer<uchar> m_type;

	Vector4 m_p;
	Vector4 m_d;

	uint m_nsize; // number of nodes
	uint m_nlocs; // number of local nodes
	Buffer<uint> m_nbeg; // begin index for nodes

	uint m_esize; // number of edges
	uint m_elocs; // number of local edges
	Buffer<uint> m_ebeg; // begin index for edges

	uint m_fsize; // number of faces
	uint m_flocs; // number of local faces
	Buffer<uint> m_fbeg; // begin index for faces

	uint m_bsize; // number of bodies
	uint m_blocs; // number of local bodies
	Buffer<uint> m_bbeg; // begin index for bodies

	uint m_qsize; // number of quads
	uint m_qlocs; // number of local quads
	Buffer<uint> m_qbeg; // begin index for quads

	// lists of polynomials
	map<Type, BlockInterpolator*> m_poly[fg_num]; // list of interpolator polynomials

	void finalizeInit(const uchar *type, const Uint4 &s);
	uint getLocals(const uint dim) const {
		if(dim == 0) return m_nlocs;
		if(dim == 1) return m_elocs;
		if(dim == 2) return m_flocs;
		if(dim == 3) return m_blocs;
		return m_qlocs;
	}
	uint getSize(const uint dim) const {
		if(dim == 0) return m_nsize;
		if(dim == 1) return m_esize;
		if(dim == 2) return m_fsize;
		if(dim == 3) return m_bsize;
		return m_qsize;
	}
	uint getBegin(const uint dim, const uint ii) const {
		if(dim == 0) return m_nbeg[ii];
		if(dim == 1) return m_ebeg[ii];
		if(dim == 2) return m_fbeg[ii];
		if(dim == 3) return m_bbeg[ii];
		return m_qbeg[ii];
	}

	const PartMesh *createPart(const uint ii, const uint part);
	const PartMesh *getPart(const uint ii, const uint part) const {
		map<Type, PartMesh*>::iterator it = m_part[part].find(getType(ii, part));
		if(it != m_part[part].end()) return it->second;
		return NULL;
	}
	Type getType(const uint ii) const;
	Type getType(const uint ii, const uint part) const;

	const Vector4 getPosition(const Uint4 &i) const { return Vector4(m_p.x + i.x * m_d.x, m_p.y + i.y * m_d.y, m_p.z + i.z * m_d.z, m_p.t + i.t * m_d.t); }
	uint getBlockId(const Uint4 &i) const { return ((i.t * m_s.z + i.z) * m_s.y + i.y) * m_s.x + i.x; }

//	uint integrateVectors(const FormGrade grade, Buffer<double> &result) const;

	void synchronizeNeighbors(const Uint4 &i0, const Uint4 &i1);
	ullong getWeight(Uint4 i0, Uint4 i1, const uchar *type, const Uint4 &s, const uint rank0, const uint rank1);

	void insertNodes(BuilderMesh &mesh, const uint type, const Vector4 &p) const;
	void removeNonconvex(BuilderMesh &mesh, const uint type, const Vector4 &p) const;
	void setPartFlags(BuilderMesh &mesh, const PartMesh *part, const Buffer<uint> &nbeg, const uint flag, const uint flags) const;
	void cutSecure(BuilderMesh &mesh, const Vector4 &p0, const Vector4 &p1, const uint flag) const;
//	void cutZero(BuilderMesh &mesh, const Vector4 &p0, const Vector4 &p1, const uint flag) const;

	void computeNodes();
	void computeEdges();
	void computeFaces();
	void computeBodies();
	void computeQuads();

	Buffer< pair<uint,uint> > getExternalNodes() const;
	Buffer< pair<uint,uint> > getExternalEdges() const;
	Buffer< pair<uint,uint> > getExternalFaces() const;
	Buffer< pair<uint,uint> > getExternalBodies() const;
	Buffer< pair<uint,uint> > getExternals(const uint dim) const {
		if(dim == 0) return getExternalNodes();
		if(dim == 1) return getExternalEdges();
		if(dim == 2) return getExternalFaces();
		if(dim == 3) return getExternalBodies();
		return Buffer< pair<uint,uint> >();
	}
	Buffer< pair<uint,uint> > getMyExternals(const Buffer< pair<uint,uint> > &ext) const;

	void addNodeToMesh(Mesh &mesh, const PartMesh *part, const uint n, const Uint4 &i) const;
	void addEdgeToMesh(Mesh &mesh, const PartMesh *part, const uint e, const Uint4 &i) const;
	void addFaceToMesh(Mesh &mesh, const PartMesh *part, const uint f, const Uint4 &i) const;
	void addBodyToMesh(Mesh &mesh, const PartMesh *part, const uint b, const Uint4 &i) const;
	void addQuadToMesh(Mesh &mesh, const PartMesh *part, const uint q, const Uint4 &i) const;

	Buffer<uint> getNodeLinks(const Uint4 &i) const;
	Buffer<uint> getEdgeLinks(const Uint4 &i, const uint part) const;
	Buffer<uint> getFaceLinks(const Uint4 &i, const uint part) const;
	Buffer<uint> getBodyLinks(const Uint4 &i, const uint part) const;

	void createFormIntegrator(const FormGrade grade, map<Type, BlockIntegrator*> &blocks, const int num) const;
	void createWedgeIntegrator(const FormGrade grade, map<Type, BlockIntegrator*> &blocks, const int num) const;
	void createInterpolator(const FormGrade grade);

	void addEdgeIncidences(const PartMesh *part, const uint e, const Uint4 &i, Buffer< pair<uint,sign> > &inc) const;
	void addFaceIncidences(const PartMesh *part, const uint f, const Uint4 &i, Buffer< pair<uint,sign> > &inc) const;
	void addBodyIncidences(const PartMesh *part, const uint b, const Uint4 &i, Buffer< pair<uint,sign> > &inc) const;
	void addQuadIncidences(const PartMesh *part, const uint q, const Uint4 &i, Buffer< pair<uint,sign> > &inc) const;

	uint getNode(const PartMesh *part, const uint n, const Uint4 &i) const;
	uint getNode(const pair<uint,uint> &ext, const Uint4 &i) const;
	uint getEdge(const PartMesh *part, const uint e, const Uint4 &i) const;
	uint getEdge(const pair<uint,uint> &ext, const Uint4 &i) const;
	uint getFace(const PartMesh *part, const uint f, const Uint4 &i) const;
	uint getFace(const pair<uint,uint> &ext, const Uint4 &i) const;
	uint getBody(const PartMesh *part, const uint b, const Uint4 &i) const;
	uint getBody(const pair<uint,uint> &ext, const Uint4 &i) const;
	uint getCell(const uint dim, const pair<uint,uint> &ext, const Uint4 &i) const {
		if(dim == 0) return getNode(ext, i);
		if(dim == 1) return getEdge(ext, i);
		if(dim == 2) return getFace(ext, i);
		return getBody(ext, i);
	}

	Buffer<uint> getPrevRanks() const;
	Buffer<uint> getNextRanks() const;

	template<typename T> void integrate(T func(const Buffer<double> &), const map<Type, BlockIntegrator*> &blocks, const uint gdim, const bool dual, Discrete<T> &result) const
	{
		// reserve space for external terms
		Buffer< pair<uint,uint> > ext;
		if(dual) ext = getExternals(gdim);
		Buffer<T> add(ext.size(), result.m_zero);

		// use blocks to integrate over mesh elements
		Uint4 i;
		uint j;
		const uint locs = result.m_height;
		Buffer<T *> term;
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
								term[j] = (cell < locs ? &result.m_val[cell] : &add[cell - locs]);
							}
						}
						block->integrate(func, getPosition(i), &result.m_val[getBegin(gdim, ii)], term);
					}
				}
			}
		}

		// communicate with external terms
		if(dual)
		{
			const Buffer< pair<uint,uint> > mext = getMyExternals(ext);
			for(j=0; j<ext.size(); j++)
			{
				sendMPI(&add[j], sizeof(T), ext[j].first, 0);
			}
			T addition;
			for(j=0; j<mext.size(); j++)
			{
				recvMPI(&addition, sizeof(T), mext[j].first, 0);
				result.m_val[mext[j].second] += addition;
			}
		}
	}
	template<typename T, typename V> Diagonal<T> &divideByVectorSquare(V func(const Buffer<double> &), const V &zero, const FormGrade grade, Diagonal<T> &result) const {
		Column<V> vec(result.m_height, zero);
		integrateForm(func, 0, grade, vec);
		for(uint j=0; j<result.m_height; j++) result.m_val[j] /= vec.m_val[j].lensq();
		return result;
	}
	template<typename T> void sumInterpolation(const FormGrade grade, const Column<T> &form, const Vector4 &p, const Uint4 &i, T *result) const {
		const uint ii = getBlockId(i);
		const Vector4 d(p - getPosition(i));
		const uint gdim = FormGradeDimension(grade);
		const BlockInterpolator *poly = m_poly[grade].find(getType(ii))->second;
		if(gdim == 0) poly->interpolate(d, &form.m_val[m_nbeg[ii]], result);
		else if(gdim == 1) poly->interpolate(d, &form.m_val[m_ebeg[ii]], result);
		else if(gdim == 2) poly->interpolate(d, &form.m_val[m_fbeg[ii]], result);
		else if(gdim == 3) poly->interpolate(d, &form.m_val[m_bbeg[ii]], result);
		else poly->interpolate(d, &form.m_val[m_qbeg[ii]], result);
	}
	template<typename Key, typename Link> void clearMap(map<Key,Link*> &mapi) const {
		for(auto it=mapi.begin(); it!= mapi.end(); it++) {
			if(it->second != NULL) delete it->second;
		}
		mapi.clear();
	}

};

}

#endif //_BLOCKMESH_HPP_INCLUDED_
