/*
BlockInterpolator.hpp interpolates discrete Form into continuous function. Works in association with BlockMesh
*/

#ifndef _BLOCKINTERPOLATOR_HPP_INCLUDED_
#define _BLOCKINTERPOLATOR_HPP_INCLUDED_

#include "BlockIntegrator.hpp"

namespace gfd
{

class BlockInterpolator
{
public:
	BlockInterpolator();
	virtual ~BlockInterpolator() { clear(); }
	void clear();

	void init(const FormGrade grade, const PartMesh &mesh, const Vector4 &d, const uint num);
	template<typename T> void interpolate(const Vector4 &p, const T *val, T *result) const {
		uint l = 0;
		const VectorN fac = getWeight(p) * getOrderFactors(p);
		for(uint i=0; i<m_values; i++) {
			for(uint j=0; j<m_order; j++) {
				const T ival = fac[j] * val[i];
				for(uint k=0; k<m_fields; k++) result[k] += m_getter[l++] * ival;
			}
		}
	}

protected:

	uint m_dim; // dimension of the space
	uint m_fields; // size of the field variable
	uint m_values; // number of discrete values to be used
	uint m_order; // number of terms in the polynomial
	Buffer<double> m_getter; // coded factors for interpolation
	Vector4 m_d; // block size

	Buffer<bool> getBoundaries(const Mesh &mesh, const FormGrade grade) const;
	uint getOrder(const uint order) const;
	VectorN getOrderFactors(const Vector4 &p) const;
	double getWeight(const Vector4 &p) const;

};

}

#endif //_BLOCKINTERPOLATOR_HPP_INCLUDED_
