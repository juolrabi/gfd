#include "Quadrature.hpp"
#include <iostream>

using namespace gfd;
using namespace std;

Quadrature::Quadrature(const Quadrature &q) 
: Buffer<double>(q) {
	m_vdim = q.m_vdim;
	m_pdim = q.m_pdim;
	m_pcount = q.m_pcount;
}
Quadrature::Quadrature(const uint vdim, const uint pdim, const uint pcount, const uint vcount) { 
	init(vdim, pdim, pcount); 
	reserve(vcount);
}

void Quadrature::init(const uint vdim, const uint pdim, const uint pcount) {
	m_vdim = vdim;
	m_pdim = pdim;
	m_pcount = pcount;
	Buffer<double>::clear();
}
void Quadrature::init(const Quadrature &l, const Quadrature &r) { // product of two quadratures
	if(l.m_vdim != r.m_vdim) return;
	if(l.m_pdim != r.m_pdim) return;
	init(l.m_vdim * (l.m_vdim + 1) / 2, l.m_pdim, l.m_pcount * r.m_pcount);
	reserve(l.vcount() * r.vcount());
	uint size = 0;
	for(uint li=0; li<l.m_size; li+=l.vsize()) {
		const double *lv = &l.m_data[li];
		for(uint ri=0; ri<r.m_size; ri+=r.vsize()) {
			const double *rv = &r.m_data[ri];
			for(uint lj=0; lj<l.m_vdim; lj++) {
				for(uint rj=0; rj<lj; rj++) m_data[size++] = lv[lj] * rv[rj] + lv[rj] * rv[lj];
				m_data[size++] = lv[lj] * rv[lj];
			}
			for(uint lj=0; lj<l.m_pcount; lj++) {
				const double *lp = &lv[m_vdim + m_pdim * lj];
				for(uint rj=0; rj<r.m_pcount; rj++) {
					const double *rp = &rv[m_vdim + m_pdim * rj];
					for(uint k=0; k<m_pdim; k++) m_data[size++] = lp[k] + rp[k];
				}
			}
		}
	}
}
void Quadrature::reserve(const uint vcount) {
	resize(vsize() * vcount);
}

Quadrature &Quadrature::relocate(const Vector4 &p) {
	for(uint i=0; i<m_size; ) {
		i += m_vdim;
		for(uint j=0; j<m_pcount; j++) {
			m_data[i++] += p.x;
			if(m_pdim >= 2) m_data[i++] += p.y;
			if(m_pdim >= 3) m_data[i++] += p.z;
			if(m_pdim >= 4) m_data[i++] += p.t;
		}
	}
	return *this;
}
Quadrature &Quadrature::invert(const Vector4 &p) {
	uint i, j;
	// sum up vectors
	Buffer<double> v(m_vdim, 0.0);
	for(i=0; i<m_size; ) {
		for(j=0; j<m_vdim; j++, i++) v[j] += m_data[i];
		for(j=0; j<m_pcount; j++) {
			m_data[i++] -= p.x;
			if(m_pdim >= 2) m_data[i++] -= p.y;
			if(m_pdim >= 3) m_data[i++] -= p.z;
			if(m_pdim >= 4) m_data[i++] -= p.t;
		}
	}
	// compute square of summed vector
	double sq = 0.0;
	for(i=0; i<m_vdim; i++) sq += v[i] * v[i];
	if(sq == 0.0) return *this;
	if(m_pcount > 0) sq *= m_pcount * m_pcount;
	// rescale vectors
	const double div = 1.0 / sq;
	for(i=0; i<m_size; i+=psize()) {
		for(j=0; j<m_vdim; j++, i++) m_data[i] *= div;
	}
	return *this;
}

void Quadrature::push(const double &v, uint &size) const {
	m_data[size++] = v;
}
void Quadrature::push(const Vector4 &v, uint &size) const {
	m_data[size++] = v.x;
	if(m_pdim == 1) return;
	m_data[size++] = v.y;
	if(m_pdim == 2) return;
	m_data[size++] = v.z;
	if(m_pdim == 3) return;
	m_data[size++] = v.t;
}
void Quadrature::push(const TwoVector4 &v, uint &size) const {
	m_data[size++] = v.xy;
	if(m_pdim == 2) return;
	m_data[size++] = v.xz;
	m_data[size++] = v.yz;
	if(m_pdim == 3) return;
	m_data[size++] = v.xt;
	m_data[size++] = v.yt;
	m_data[size++] = v.zt;
}
void Quadrature::push(const ThreeVector4 &v, uint &size) const {
	m_data[size++] = v.xyz;
	if(m_pdim == 3) return;
	m_data[size++] = v.xyt;
	m_data[size++] = v.xzt;
	m_data[size++] = v.yzt;
}
void Quadrature::push(const FourVector4 &v, uint &size) const {
	m_data[size++] = v.xyzt;
}
void Quadrature::push(const Buffer<double> &v, uint &size) const {
	for(uint i=0; i<v.size(); i++) m_data[size++] = v[i];
}
void Quadrature::push(const double *v, const uint vs, uint &size) const {
	for(uint i=0; i<vs; i++) m_data[size++] = v[i];
}

