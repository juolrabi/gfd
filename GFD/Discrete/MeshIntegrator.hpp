/**
 * MeshIntegrator implements tools for integration over mesh elements.
 * variable num defines the number of integration points per element:
 *   num > 0 leads to speed-optimized integration
 *   num = 0 leads to position-independent integration (constant field)
 *   num < 0 leads to primal-dual-partitioned integration
 * Author: Jukka Räbinä, University of Jyväskylä, 2020.
 */

#ifndef _MESHINTEGRATOR_HPP_INCLUDED_
#define _MESHINTEGRATOR_HPP_INCLUDED_

#include "Form.hpp"
#include "Quadrature.hpp"
#include "../Mesh/PartMesh.hpp"

namespace gfd
{

class MeshIntegrator
{
public:
	MeshIntegrator(const PartMesh &mesh, const FormGrade grade, const int num, const uint lowdim, const uint highdim);
	virtual ~MeshIntegrator() { }

	Quadrature getEmptyQuadrature() const;
	Quadrature &getQuadrature(const uint i, Quadrature &q) const;
	Buffer< pair<uint, Quadrature> > &getBaseQuadrature(const uint i, Buffer< pair<uint,Quadrature> > &q) const;

	uint getLocals() const { return getLocals(FormGradeDimension(m_grade)); }
	uint getSize() const { return getSize(FormGradeDimension(m_grade)); }
	const Buffer< pair<uint,uint> > &getExternals() const { return getExternals(FormGradeDimension(m_grade)); }
	uint getFlag(const uint i) const { return getFlag(i, FormGradeDimension(m_grade)); }
	Vector4 getPosition(const uint i) const { return getPosition(i, FormGradeDimension(m_grade)); }

	uint getBaseLocals() const { return getLocals(getBaseDimension()); }
	uint getBaseSize() const { return getSize(getBaseDimension()); }
	const Buffer< pair<uint,uint> > &getBaseExternals() const { return getExternals(getBaseDimension()); }
	uint getBaseFlag(const uint i) const { return getFlag(i, getBaseDimension()); }
	Vector4 getBasePosition(const uint i) const { return getPosition(i, getBaseDimension()); }


protected:
	const PartMesh &m_mesh;
	FormGrade m_grade;
	int m_num;
	uint m_lowdim;
	uint m_highdim;

	Quadrature &getVector(const uint i, Quadrature &q) const;
	Buffer< pair<uint,Quadrature> > &getBaseVector(const uint i, Buffer< pair<uint,Quadrature> > &q) const;

	template<typename V, typename B, typename R> Quadrature &createVector(const Buffer<uint> &j, const Buffer<V> &v, R vector(const V &, const B &), B base(const PartMesh &, const uint), Quadrature &q) const;
	template<typename V, typename B, typename R> Buffer< pair<uint,Quadrature> > &createBaseVector(const Buffer<uint> &j, const Buffer<V> &v, R vector(const V &, const B &), const B &b, Buffer< pair<uint,Quadrature> > &q) const;
	template<typename V> Quadrature &createSingleVector(const V &v, Quadrature &q) const;
	template<typename V> Buffer< pair<uint,Quadrature> > &createSingleBaseVector(const uint i, const V &v, Buffer< pair<uint,Quadrature> > &q) const;
	void removeExternalSimplices(const uint locs, Buffer<uint> &j, Buffer<Vector4> &p) const;
	template<typename V> void removeExternalVectors(const uint locs, Buffer<uint> &j, Buffer<V> &v) const;
	template<typename V> Quadrature &createEntryQuadrature(const uint pdim, const Buffer<Vector4> &p, V func(const Vector4 *), Quadrature &q) const;
	template<typename V, typename B> Quadrature &createQuadrature(const Buffer<uint> &j, const Buffer<Vector4> &p, V vector(const Vector4 *, const B &), B base(const PartMesh &, const uint), Quadrature &q) const;
	template<typename V, typename B> Buffer< pair<uint,Quadrature> > &createBaseQuadrature(const Buffer<uint> &j, const Buffer<Vector4> &p, V vector(const Vector4 *, const B &), const B &b, Buffer< pair<uint,Quadrature> > &q) const;
	template<typename V> Quadrature &createSingleQuadrature(const Vector4 &p, const V &v, Quadrature &q) const;
	template<typename V> Buffer< pair<uint,Quadrature> > &createSingleBaseQuadrature(const uint i, const Vector4 &p, const V &v, Buffer< pair<uint,Quadrature> > &q) const;
	Quadrature &createMultiQuadrature(const Buffer<Vector4> &p, const Quadrature &v, Quadrature &q) const;

	uint getBaseDimension() const { return (FormGradeIsPrim(m_grade) ? m_lowdim : m_highdim); }
	uint getLocals(const uint gdim) const;
	uint getSize(const uint gdim) const;
	const Buffer< pair<uint,uint> > &getExternals(const uint gdim) const;
	uint getFlag(const uint i, const uint gdim) const;
	Vector4 getPosition(const uint i, const uint gdim) const;

};

}

#endif //_MESHINTEGRATOR_HPP_INCLUDED_
