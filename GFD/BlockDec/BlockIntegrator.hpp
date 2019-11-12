/*
BlockIntegrator.hpp implements DeRham map (integrates continuous function into discrete form)
*/

#ifndef _BLOCKINTEGRATOR_HPP_INCLUDED_
#define _BLOCKINTEGRATOR_HPP_INCLUDED_

#include "../Discrete/Form.hpp"
#include "../Mesh/PartMesh.hpp"

//#include "../Types/MpiEasy.hpp"

namespace gfd
{

class BlockIntegrator
{
public:
	static uint getFieldDimension(const FormGrade grade, const uint dim) {
		const uint gdim = FormGradeDimension(grade);
		if(gdim == 0 || gdim == 4) return 1;
		if(gdim == 1) return dim;
		if(gdim == 2) return dim * (dim - 1) / 2;
		return dim * (dim - 1) * (dim - 2) / 6;
	}

	BlockIntegrator();
	virtual ~BlockIntegrator() { clear(); }
	void clear();

	void init(const FormGrade grade, const PartMesh &mesh, const uint num = 2);
	void initWedge(const FormGrade grade, const PartMesh &mesh, const uint num = 2);
	const Buffer< pair<uint,uint> > &getExternal() const { return m_ext; }
	template<typename T> void integrate(void func(const Vector4 &, T *result), const Vector4 &p0, T *result, const Buffer<T *> &exterm) const {
		Buffer<T> f(m_fields);
		for(uint i=0; i<m_setter.size(); i++)
		{
			T &ival = (i < m_values ? result[i] : *exterm[i - m_values]);
			const Buffer<double> &setteri = m_setter[i];
			for(uint j=m_fields; j<setteri.size(); )
			{
				const double fac = setteri[j++];
				Vector4 p(setteri[j++] + p0.x,0,0,0);
				if(m_dim >= 2) p.y = setteri[j++] + p0.y;
				if(m_dim >= 3) p.z = setteri[j++] + p0.z;
				if(m_dim >= 4) p.t = setteri[j++] + p0.t;
				func(p, &f[0]);
				for(uint k=0; k<m_fields; k++) ival += fac * setteri[k] * f[k];
			}
		}
	}
	void getVectors(double *result, const Buffer<double *> &exterm) const;

	uint getNumberOfLocalValues() const { return m_values; }
	const Buffer<double> &getIntegrationData(const uint i) const { return m_setter[i]; }

protected:

	uint m_dim; // dimension of the space
	uint m_fields; // size of the field variable
	uint m_values; // number of discrete values to be used
	Buffer< Buffer<double> > m_setter;
	Buffer< pair<uint,uint> > m_ext;
};

}

#endif //_BLOCKINTEGRATOR_HPP_INCLUDED_
