/**
 * Quadrature is a tool for numerical integration over a cell.
 * Author: Jukka Räbinä, University of Jyväskylä, 2020.
 */

#ifndef _QUADRATURE_HPP_INCLUDED_
#define _QUADRATURE_HPP_INCLUDED_

#include "../Types/Types.hpp"

namespace gfd
{

// default functions for integration
inline double FormVector1(const Buffer<double> &q) { return q[0]; }
inline Vector2 FormVector2(const Buffer<double> &q) { return Vector2(q[0], q[1]); }
inline Vector3 FormVector3(const Buffer<double> &q) { return Vector3(q[0], q[1], q[2]); }
inline Vector4 FormVector4(const Buffer<double> &q) { return Vector4(q[0], q[1], q[2], q[3]); }
inline TwoVector2 FormTwoVector2(const Buffer<double> &q) { return TwoVector2(q[0]); }
inline TwoVector3 FormTwoVector3(const Buffer<double> &q) { return TwoVector3(q[0], q[1], q[2]); }
inline TwoVector4 FormTwoVector4(const Buffer<double> &q) { return TwoVector4(q[0], q[1], q[2], q[3], q[4], q[5]); }
inline ThreeVector3 FormThreeVector3(const Buffer<double> &q) { return ThreeVector3(q[0]); }
inline ThreeVector4 FormThreeVector4(const Buffer<double> &q) { return ThreeVector4(q[0], q[1], q[2], q[3]); }
inline FourVector4 FormFourVector4(const Buffer<double> &q) { return FourVector4(q[0]); }
inline VectorN FormVectorN(const Buffer<double> &q) { return VectorN(q); }
inline double HodgeUnit1(const Buffer<double> &q) { return q[0]; }
inline double HodgeUnit2(const Buffer<double> &q) { return q[0] + q[2]; }
inline double HodgeUnit3(const Buffer<double> &q) { return q[0] + q[2] + q[5]; }
inline double HodgeUnit4(const Buffer<double> &q) { return q[0] + q[2] + q[5] + q[9]; }
inline double HodgeUnit6(const Buffer<double> &q) { return q[0] + q[2] + q[5] + q[9] + q[14]; }


class Quadrature : public Buffer<double>
{
public:
	Quadrature() { }
	Quadrature(const Quadrature &q);
	Quadrature(const uint vdim, const uint pdim, const uint pcount) { init(vdim, pdim, pcount); }
	Quadrature(const uint vdim, const uint pdim, const uint pcount, const uint vcount);
	Quadrature(const Quadrature &l, const Quadrature &r) { init(l, r); }
	virtual ~Quadrature() { }

	void init(const uint vdim, const uint pdim, const uint pcount);
	void init(const Quadrature &l, const Quadrature &r); // product of two quadratures
	void reserve(const uint vcount);

	Quadrature &relocate(const Vector4 &p);
	Quadrature &invert(const Vector4 &p);
	void push(const double &v, uint &size) const;
	void push(const Vector4 &v, uint &size) const;
	void push(const TwoVector4 &v, uint &size) const;
	void push(const ThreeVector4 &v, uint &size) const;
	void push(const FourVector4 &v, uint &size) const;
	void push(const Buffer<double> &v, uint &size) const;
	void push(const double *v, const uint vs, uint &size) const;

	template<typename T> T &integrate(T func(const Buffer<double> &), T &val) const {
		Buffer<double> data(m_vdim + m_pdim);
		for(uint i=0; i<m_size; ) {
			memcpy(&data[0], &m_data[i], m_vdim * sizeof(double));
			i += m_vdim;
			for(uint j=0; j<m_pcount; j++) { // integrate over positions
				memcpy(&data[m_vdim], &m_data[i], m_pdim * sizeof(double));
				i += m_pdim;
				val += func(data);
			}
			if(m_pcount == 0) val += func(data); // integrate constant
		}
		return val;
	}
	template<typename T> T &integrateRelocated(T func(const Buffer<double> &), const Vector4 &p, T &val) const {
		Buffer<double> data(m_vdim + m_pdim);
		for(uint i=0; i<m_size; ) {
			memcpy(&data[0], &m_data[i], m_vdim * sizeof(double));
			i += m_vdim;
			for(uint j=0; j<m_pcount; j++) { // integrate over positions
				data[m_vdim] = m_data[i++] + p.x;
				if(m_pdim >= 2) data[m_vdim + 1] = m_data[i++] + p.y;
				if(m_pdim >= 3) data[m_vdim + 2] = m_data[i++] + p.z;
				if(m_pdim >= 4) data[m_vdim + 3] = m_data[i++] + p.t;
				val += func(data);
			}
			if(m_pcount == 0) val += func(data); // integrate constant
		}
		return val;
	}
	template<typename T> T &integrateProduct(const Quadrature &r, T func(const Buffer<double> &), T &val) const {
		if(m_vdim != r.m_vdim) return val;
		if(m_pdim != r.m_pdim) return val;
		const uint pcount = m_pcount * r.m_pcount;
		Buffer<double> data(m_vdim * (m_vdim + 1) / 2 + m_pdim);
		for(uint li=0; li<m_size; li+=vsize()) {
			const double *lv = &m_data[li];
			for(uint ri=0; ri<r.m_size; ri+=r.vsize()) {
				const double *rv = &r.m_data[ri];
				uint size = 0;
				for(uint lj=0; lj<m_vdim; lj++) {
					for(uint rj=0; rj<lj; rj++) data[size++] = lv[lj] * rv[rj] + lv[rj] * rv[lj];
					data[size++] = lv[lj] * rv[lj];
				}
				for(uint lj=0; lj<m_pcount; lj++) {
					const double *lp = &lv[m_vdim + m_pdim * lj];
					for(uint rj=0; rj<r.m_pcount; rj++) {
						const double *rp = &rv[m_vdim + m_pdim * rj];
						for(uint k=0; k<m_pdim; k++) data[size + k] = lp[k] + rp[k];
						val += func(data);
					}
				}
				if(pcount == 0) val += func(data);
			}
		}
		return val;
	}

	uint vdim() const { return m_vdim; } // number of doubles per vector
	uint pdim() const { return m_pdim; } // number of doubles per position
	uint pcount() const { return m_pcount; } // number of positions per vector block
	uint vsize() const { return m_vdim + m_pdim * m_pcount; } // number of doubles per vector block
	uint psize() const { return m_pdim * m_pcount; } // number of position doubles per vector block
	uint vcount() const { return m_size / vsize(); } // number of vector blocks

protected:
	uint m_vdim; // number of doubles per vector
	uint m_pdim; // number of doubles per position
	uint m_pcount; // number of positions per vector block

};

}

#endif //_QUADRATURE_HPP_INCLUDED_
