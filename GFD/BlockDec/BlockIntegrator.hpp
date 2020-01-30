/**
 * BlockIntegrator.hpp implements DeRham map (integrates continuous function into discrete forms)
 * Author: Jukka Räbinä, University of Jyväskylä, 2020.
 */

#ifndef _BLOCKINTEGRATOR_HPP_INCLUDED_
#define _BLOCKINTEGRATOR_HPP_INCLUDED_

#include "../Discrete/Quadrature.hpp"
#include "../Discrete/Form.hpp"
#include "../Mesh/PartMesh.hpp"

namespace gfd
{

class BlockIntegrator
{
public:
	BlockIntegrator();
	virtual ~BlockIntegrator() { clear(); }
	void clear();

	void init(const FormGrade grade, const PartMesh &mesh, const int num);
	void initWedge(const FormGrade grade, const PartMesh &mesh, const int num);
	template<typename T> void integrate(T func(const Buffer<double> &), const Vector4 &p, T *result, const Buffer<T *> &exterm) const {
		for(uint i=0; i<m_q.size(); i++) {
			T &ival = (i < m_values ? result[i] : *exterm[i - m_values]);
			m_q[i].integrateRelocated(func, p, ival);
		}
	}

	uint getLocals() const { return m_values; }
	const Buffer< pair<uint,uint> > &getExternals() const { return m_ext; }
	const Quadrature &getQuadrature(const uint i) const { return m_q[i]; }

protected:

	Buffer<Quadrature> m_q; // defines numerical integration
	uint m_values; // number of local terms
	Buffer< pair<uint,uint> > m_ext; // external term data
};

}

#endif //_BLOCKINTEGRATOR_HPP_INCLUDED_
