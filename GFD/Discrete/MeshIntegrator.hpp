/**
 * MeshIntegrator implements tools for integration over mesh elements.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _MESHINTEGRATOR_HPP_INCLUDED_
#define _MESHINTEGRATOR_HPP_INCLUDED_

#include "Form.hpp"
#include "../Mesh/PartMesh.hpp"

namespace gfd
{

class MeshIntegrator
{
public:
	MeshIntegrator(const PartMesh &mesh, const FormGrade grade, const uint num, const uint lowdim, const uint highdim);
	virtual ~MeshIntegrator() { }

	Buffer<double> &getVector(const uint i, Buffer<double> &q) const;
	Buffer< pair<uint, Buffer<double> > > &getBaseVector(const uint i, Buffer< pair<uint, Buffer<double> > > &q) const;
	Buffer<double> &getQuadrature(const uint i, Buffer<double> &q) const;
	Buffer< pair<uint, Buffer<double> > > &getBaseQuadrature(const uint i, Buffer< pair<uint, Buffer<double> > > &q) const;

	Buffer<double> &invertVector(Buffer<double> &q) const;
	Buffer<double> &invertQuadrature(const Vector4 &p, Buffer<double> &q) const;

	template<typename T> void integrateVector(const Buffer<double> &q, const T *f, T &val) const {
		for(uint i=0; i<q.size(); i++) val += q[i] * f[i];
	}
	template<typename T> void integrateQuadrature(const Buffer<double> &q, void func(const Vector4 &, T *), T &val) const {
		const uint dim = m_mesh.getDimension();
		const uint fields = getFields();
		Buffer<T> f(fields);
		for(uint i=0; i<q.size(); ) {
			Vector4 p(q[i++],0,0,0);
			if(dim >= 2) p.y = q[i++];
			if(dim >= 3) p.z = q[i++];
			if(dim >= 4) p.t = q[i++];
			func(p, &f[0]);
			for(uint j=0; j<f.size(); j++) val += q[i++] * f[j];
		}
	}
	template<typename T> void integrateUnitProductVector(const Buffer<double> &q0, const Buffer<double> &q1, T &val) const {
		if(q0.size() != q1.size()) return;
		for(uint i=0; i<q0.size(); i++) val += q0[i] * q1[i];
	}
	template<typename T> void integrateProductVector(const Buffer<double> &q0, const Buffer<double> &q1, const T *f, T &val) const {
		if(q0.size() != q1.size()) return;
		uint n = 0;
		for(uint i=0; i<q0.size(); i++) {
			for(uint j=0; j<i; j++) val += f[n++] * (q0[i] * q1[j] + q0[i] * q1[j]);
			val += f[n++] * (q0[i] * q1[i]);
		}
	}
	template<typename T> void integrateProductQuadrature(const Buffer<double> &q0, const Buffer<double> &q1, void func(const Vector4 &, T *), T &val) const {
		const uint dim = m_mesh.getDimension();
		const uint fields = getFields();
		const uint size = dim + fields;
		Buffer<T> f(fields * (fields + 1) / 2);
		for(uint i=0; i<q0.size(); i+=size) {
			const double *v0 = &q0[i + dim];
			for(uint j=0; j<q1.size(); j+=size) {
				const double *v1 = &q1[j + dim];
				Vector4 p(q0[i] + q1[j],0,0,0);
				if(dim >= 2) p.y = q0[i+1] + q1[j+1];
				if(dim >= 3) p.z = q0[i+2] + q1[j+2];
				if(dim >= 4) p.t = q0[i+3] + q1[j+3];
				func(p, &f[0]);
				uint n = 0;
				for(uint k=0; k<fields; k++) {
					for(uint l=0; l<k; l++) val += f[n++] * (v0[k] * v1[l] + v0[l] * v1[k]);
					val += f[n++] * (v0[k] * v1[k]);
				}
			}
		}
	}

	uint getFields() const;
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
	uint m_num;
	uint m_lowdim;
	uint m_highdim;

	template<typename V, typename B, typename R> Buffer<double> &createVector(const Buffer<uint> &j, const Buffer<V> &v, R vector(const V &, const B &), B base(const PartMesh &, const uint), Buffer<double> &q) const;
	template<typename V, typename B, typename R> Buffer< pair<uint, Buffer<double> > > &createBaseVector(const Buffer<uint> &j, const Buffer<V> &v, R vector(const V &, const B &), const B &b, Buffer< pair<uint, Buffer<double> > > &q) const;
	template<typename V> Buffer<double> &createSingleVector(const V &v, Buffer<double> &q) const;
	template<typename V> Buffer< pair<uint, Buffer<double> > > &createSingleBaseVector(const uint i, const V &v, Buffer< pair<uint, Buffer<double> > > &q) const;
	void removeExternalSimplices(const uint locs, Buffer<uint> &j, Buffer<Vector4> &p) const;
	template<typename V> void removeExternalVectors(const uint locs, Buffer<uint> &j, Buffer<V> &v) const;
	template<typename V> Buffer<double> &createEntryQuadrature(const uint pdim, const Buffer<Vector4> &p, V func(const Vector4 *), Buffer<double> &q) const;
	template<typename V, typename B> Buffer<double> &createQuadrature(const Buffer<uint> &j, const Buffer<Vector4> &p, V vector(const Vector4 *, const B &), B base(const PartMesh &, const uint), Buffer<double> &q) const;
	template<typename V, typename B> Buffer< pair<uint, Buffer<double> > > &createBaseQuadrature(const Buffer<uint> &j, const Buffer<Vector4> &p, V vector(const Vector4 *, const B &), const B &b, Buffer< pair<uint, Buffer<double> > > &q) const;
	template<typename V> Buffer<double> &createSingleQuadrature(const Vector4 &p, const V &v, Buffer<double> &q) const;
	template<typename V> Buffer< pair<uint, Buffer<double> > > &createSingleBaseQuadrature(const uint i, const Vector4 &p, const V &v, Buffer< pair<uint, Buffer<double> > > &q) const;
	Buffer<double> &createMultiQuadrature(const uint elems, const Buffer<Vector4> &p, const Buffer<double> &v, Buffer<double> &q) const;

	uint getBaseDimension() const { return (FormGradeIsPrim(m_grade) ? m_lowdim : m_highdim); }
	uint getLocals(const uint gdim) const;
	uint getSize(const uint gdim) const;
	const Buffer< pair<uint,uint> > &getExternals(const uint gdim) const;
	uint getFlag(const uint i, const uint gdim) const;
	Vector4 getPosition(const uint i, const uint gdim) const;

	void setVector(const double &v, Buffer<double> &q, uint &qs) const;
	void setVector(const Vector4 &v, Buffer<double> &q, uint &qs) const;
	void setVector(const TwoVector4 &v, Buffer<double> &q, uint &qs) const;
	void setVector(const ThreeVector4 &v, Buffer<double> &q, uint &qs) const;
	void setVector(const FourVector4 &v, Buffer<double> &q, uint &qs) const;

};

}

#endif //_MESHINTEGRATOR_HPP_INCLUDED_
