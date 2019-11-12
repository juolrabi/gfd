/**
 * class Dec implements tools of discrete exterior calculus.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _DEC_HPP_INCLUDED_
#define _DEC_HPP_INCLUDED_

#include "MeshIntegrator.hpp"
#include "../Mesh/PartMesh.hpp"
#include "../Types/MpiEasy.hpp"

namespace gfd
{

class Dec
{
public:
	Dec(const PartMesh &mesh) : m_mesh(mesh) { }
	virtual ~Dec() { }

	Sparse<sign> &integrateDerivative(const FormGrade grade, Sparse<sign> &d) const;
	template<typename T> Column<T> &integrateForm(const FormGrade grade, void func(const Vector4 &, T *), const uint num, Column<T> &result) const {
		// initialize form
		MeshIntegrator intg(m_mesh, grade, num);
		const uint locs = intg.getLocals();
		if(!result.m_full || result.m_height != locs) result.setFullOfZeros(locs);
		if(num == 0) return result;

		// integrate over mesh elements
		uint i;
		Buffer<T> f(intg.getFields());
		Buffer<double> q;
		for(i=0; i<locs; i++) {
			intg.getSetter(i, q);
			sumIntegration(result.m_val[i], q, func, f);
		}
		if(FormGradeIsPrim(grade)) return result;
		// communicate with external terms
		T val;
		const uint rank = getMPIrank();
		const Buffer< pair<uint,uint> > &ext = intg.getExternals();
		const Buffer< pair<uint,uint> > mext = getMyExternals(ext);
		for(i=0; i<ext.size(); i++) {
			val = result.m_zero;
			intg.getSetter(locs + i, q);
			sumIntegration(val, q, func, f);
			if(ext[i].first == rank) result.m_val[ext[i].second] += val;
			else sendMPI(&val, sizeof(T), ext[i].first, 0);
		}
		for(i=0; i<mext.size(); i++) {
			recvMPI(&val, sizeof(T), mext[i].first, 0);
			result.m_val[mext[i].second] += val;
		}
		return result;
	}
	template<typename T> Sparse<T> &integrateHodge(const FormGrade grade, void func(const Vector4 &, T *), const uint num, Sparse<T> &result) const {
		// initialize form
		MeshIntegrator intg(m_mesh, grade, num);
		const uint locs = intg.getLocals();
		if(!result.m_full || result.m_height != locs) result.setFullOfZeros(locs);
		if(num == 0) return result;

		// integrate over mesh elements
		uint i, j, k;
		const uint dim = m_mesh.getDimension();
		Buffer<T> f(intg.getWedgeFields());
		Buffer<double> q;
		for(i=0; i<locs; i++) {
			intg.getWedgeSetter(i, q);
			sumIntegration(result.m_val[i], q, func, f);
		}
		return result;
	}

protected:
	const PartMesh &m_mesh;

	template<typename T> void sumIntegration(T &val, const Buffer<double> &q, void func(const Vector4 &, T *), Buffer<T> &f) const {
		uint i, j;
		const uint dim = m_mesh.getDimension();
		for(i=f.size(); i<q.size(); ) {
			const double fac = q[i++];
			Vector4 p(q[i++],0,0,0);
			if(dim >= 2) p.y = q[i++];
			if(dim >= 3) p.z = q[i++];
			if(dim >= 4) p.t = q[i++];
			func(p, &f[0]);
			for(j=0; j<f.size(); j++) val += fac * q[j] * f[j];
		}
	}

	Buffer< pair<uint,uint> > getMyExternals(const Buffer< pair<uint,uint> > &ext) const;

};

}

#endif //_DEC_HPP_INCLUDED_
