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
	MeshIntegrator(const PartMesh &mesh, const FormGrade grade, const uint num);
	virtual ~MeshIntegrator() { }

	uint getFields() const;
	uint getWedgeFields() const;

	void gatherSetter0(const uint i, Buffer<double> &q, uint &qs) const;
	void gatherWedgeSetter0(const uint i, Buffer<double> &q, uint &qs) const;

	void gatherSetter(const uint i, Buffer<double> &q, uint &qs) const;
	void gatherWedgeSetter(const uint i, Buffer<double> &q, uint &qs) const;

	uint getLocals() const;
	const Buffer< pair<uint,uint> > &getExternals() const;

protected:
	const PartMesh &m_mesh;
	FormGrade m_grade;
	uint m_num;
	
	void gatherQuadrature(const Vector4 &p, const double w, Buffer<double> &q, uint &qs) const;
	void gatherQuadrature(const Buffer<Vector4> &p, const double w, Buffer<double> &q, uint &qs) const;
	void gatherVector(const Vector4 &v, Buffer<double> &q, uint &qs) const;
	void gatherVector(const TwoVector4 &v, Buffer<double> &q, uint &qs) const;
	void gatherVector(const ThreeVector4 &v, Buffer<double> &q, uint &qs) const;
	void gatherVector(const FourVector4 &v, Buffer<double> &q, uint &qs) const;
	template<typename V> void gatherQuadrature(const Buffer<Vector4> &p, const uint ps, const Buffer<V> &v, const uint vs, const bool dual, Buffer<double> &q, uint &qs) const;
	void gatherWedgeQuadrature(const Buffer<double> &prim, const uint prims, const Buffer<double> &dual, const uint duals, const Vector4 &p0, Buffer<double> &q, uint &qs) const;

};

}

#endif //_MESHINTEGRATOR_HPP_INCLUDED_
