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

	void getSetter(const uint i, Buffer<double> &q) const;
	uint getFields() const;
	
	void getWedgeSetter(const uint i, Buffer<double> &q) const;
	uint getWedgeFields() const;

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
	template<typename V> void createQuadrature(Buffer<double> &q, const Buffer<Vector4> &p, const uint ps, const Buffer<V> &v, const uint vs, const bool dual) const;
	void createWedgeQuadrature(Buffer<double> &q, const Buffer<double> &prim, const Buffer<double> &dual, const Vector4 &p0) const;

};

}

#endif //_MESHINTEGRATOR_HPP_INCLUDED_
