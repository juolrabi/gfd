#include "BlockInterpolator.hpp"
#include <iostream>
using namespace std;
using namespace gfd;

BlockInterpolator::BlockInterpolator() {
	m_dim = 1;
	m_fields = 1;
	m_values = 0;
	m_order = 1;
	m_d = Vector4(1,1,1,1);
}

void BlockInterpolator::clear() {
	m_getter.clear();
	m_dim = 1;
	m_fields = 1;
	m_values = 0;
	m_order = 1;
	m_d = Vector4(1,1,1,1);
}

void BlockInterpolator::init(const FormGrade grade, const PartMesh &mesh, const Vector4 &d, const uint num) {
	clear();

	BlockIntegrator block;
	block.init(grade, mesh, num);

	m_d = d;
	m_dim = mesh.getDimension();
	m_fields = FormGradeVectorDimension(grade, m_dim);
	m_values = block.getLocals();
	const Buffer<bool> ignore = getBoundaries(mesh, grade);

	// init m_order and m_getter
	uint i, j, k, l;
	MatrixN A;
	Buffer<VectorN> v(m_values);
	Buffer<VectorN> w(m_values);
	m_order = getOrder(m_values / m_fields);
	for(i=0; i<m_values; i++) {
		// compute v[i] and w[i] (w[i] ~= weight * v[i])
		v[i].toVectorN(m_order * m_fields);
		w[i].toVectorN(m_order * m_fields);
		if(i < ignore.size() && ignore[i]) continue;

		const Quadrature &q = block.getQuadrature(i);
		for(uint qs=0; qs<q.size();) {
			const double *qv = &q[qs];
			qs += q.vdim();
			for(j=0; j<q.pcount(); j++) {
				Vector4 p(q[qs++],0,0,0);
				if(m_dim >= 2) p.y = q[qs++];
				if(m_dim >= 3) p.z = q[qs++];
				if(m_dim >= 4) p.t = q[qs++];
				const VectorN vj = getOrderFactors(p);
				const VectorN wj = getWeight(p) * vj;
				uint jj = 0;
				for(k=0; k<m_order; k++) {
					for(l=0; l<m_fields; l++, jj++) {
						v[i][jj] += qv[l] * vj[k];
						w[i][jj] += qv[l] * wj[k];
					}
				}
			}
		}
		A += w[i].outerProduct(v[i]);
	}

	while(true) {
		// initialize m_getter
		const uint ofs = m_order * m_fields;
		m_getter.resize(m_values * ofs);
		if(m_getter.size() == 0) break;

		const MatrixN invA(A.inverse());
		if(invA.size() == ofs) {
			// A is invertible
			// update values of m_getter
			for(i=0; i<m_values; i++) {
				const VectorN vi = invA * w[i];

				// check stability
				double sum = 0.0;
				for(j=0; j<ofs; j++) sum += fabs(vi[j] * v[i][j]);
				if(sum > 2.0) break; // non-stable polynomial

				// set m_getter
				for(j=0; j<ofs; j++) m_getter[i * ofs + j] = vi[j];
			}
			if(i == m_values) break; // found stable polynomial
		}

		// decrease m_order
		m_order = getOrder(m_order - 1);
		A.toMatrixN(m_order * m_fields);
	}
}

Buffer<bool> BlockInterpolator::getBoundaries(const Mesh &mesh, const FormGrade grade) const {
	Buffer<bool> boun;
	const uint gdim = FormGradeDimension(grade);
	if(FormGradeIsPrim(grade)) return boun;
	else if(gdim == 0) boun.resize(mesh.getNodeSize());
	else if(gdim == 1) boun.resize(mesh.getEdgeSize());
	else if(gdim == 2) boun.resize(mesh.getFaceSize());
	else if(gdim == 3) boun.resize(mesh.getBodySize());
	else return boun;
	boun.fill(false);

	if(mesh.getEdgeSize() == 0) return boun; // 0d
	uint i, j;
	if(gdim < 1) {
		for(i=0; i<mesh.getNodeSize(); i++) {
			if(mesh.getNodeEdges(i).size() < 2) boun[i] = true;
		}
	}
	if(mesh.getFaceSize() == 0) return boun; // 1d
	if(gdim < 2) {
		for(i=0; i<mesh.getEdgeSize(); i++) {
			if(mesh.getEdgeFaces(i).size() < 2) {
				Buffer<uint> par;
				if(gdim == 1) par = Buffer<uint>(1, i);
				else par = mesh.getEdgeNodes(i);
				for(j=0; j<par.size(); j++) boun[par[j]] = true;
			}
		}
	}
	if(mesh.getBodySize() == 0) return boun; // 2d
	if(gdim < 3) {
		for(i=0; i<mesh.getFaceSize(); i++) {
			if(mesh.getFaceBodies(i).size() < 2) {
				Buffer<uint> par;
				if(gdim == 2) par = Buffer<uint>(1, i);
				else if(gdim == 1) par = mesh.getFaceEdges(i);
				else par = mesh.getFaceNodes(i);
				for(j=0; j<par.size(); j++) boun[par[j]] = true;
			}
		}
	}
	if(mesh.getQuadSize() == 0) return boun; // 3d
	if(gdim < 4) {
		for(i=0; i<mesh.getBodySize(); i++) {
			if(mesh.getBodyQuads(i).size() < 2) {
				Buffer<uint> par;
				if(gdim == 3) par = Buffer<uint>(1, i);
				else if(gdim == 2) par = mesh.getBodyFaces(i);
				else if(gdim == 1) par = mesh.getBodyEdges(i);
				else par = mesh.getBodyNodes(i);
				for(j=0; j<par.size(); j++) boun[par[j]] = true;
			}
		}
	}
	return boun; // 4d
}

uint BlockInterpolator::getOrder(const uint order) const
{
	uint i = 0;

	// zeroth degree
	uint j = 1;
	if(j > order) return i;
	i = j;

	// first degree
	j += m_dim;
	if(j > order) return i;
	i = j;

	// second degree
	j += m_dim * (m_dim + 1) / 2;
	if(j > order) return i;
	i = j;

	// third degree
	j += m_dim * (m_dim + 1) * (m_dim + 2) / 6;
	if(j > order) return i;
	i = j;

	// fourth degree
	j += m_dim * (m_dim + 1) * (m_dim + 2) * (m_dim + 3) / 24;
	if(j > order) return i;
	i = j;

	// fifth degree
	j += m_dim * (m_dim + 1) * (m_dim + 2) * (m_dim + 3) * (m_dim + 4) / 120;
	if(j > order) return i;
	return j;
}
VectorN BlockInterpolator::getOrderFactors(const Vector4 &p) const
{
	uint i = 0;
	VectorN v;
	v.val.resize(m_order);
	if(i == m_order) return v;

	// zeroth degree
	v[i++] = 1.0;
	if(i == m_order) return v;

	// first degree
	if(m_dim >= 1) v[i++] = p.x;
	if(m_dim >= 2) v[i++] = p.y;
	if(m_dim >= 3) v[i++] = p.z;
	if(m_dim >= 4) v[i++] = p.t;
	if(i == m_order) return v;

	// second degree
	if(m_dim >= 1) v[i++] = p.x * p.x;
	if(m_dim >= 2) {
		v[i++] = p.x * p.y;
		v[i++] = p.y * p.y;
	}
	if(m_dim >= 3) {
		v[i++] = p.x * p.z;
		v[i++] = p.y * p.z;
		v[i++] = p.z * p.z;
	}
	if(m_dim >= 4) {
		v[i++] = p.x * p.t;
		v[i++] = p.y * p.t;
		v[i++] = p.z * p.t;
		v[i++] = p.t * p.t;
	}
	if(i == m_order) return v;

	// third degree
	if(m_dim >= 1) v[i++] = p.x * p.x * p.x;
	if(m_dim >= 2) {
		v[i++] = p.x * p.x * p.y;
		v[i++] = p.x * p.y * p.y;
		v[i++] = p.y * p.y * p.y;
	}
	if(m_dim >= 3) {
		v[i++] = p.x * p.x * p.z;
		v[i++] = p.x * p.y * p.z;
		v[i++] = p.x * p.z * p.z;
		v[i++] = p.y * p.y * p.z;
		v[i++] = p.y * p.z * p.z;
		v[i++] = p.z * p.z * p.z;
	}
	if(m_dim >= 4) {
		v[i++] = p.x * p.x * p.t;
		v[i++] = p.x * p.y * p.t;
		v[i++] = p.x * p.z * p.t;
		v[i++] = p.x * p.t * p.t;
		v[i++] = p.y * p.y * p.t;
		v[i++] = p.y * p.z * p.t;
		v[i++] = p.y * p.t * p.t;
		v[i++] = p.z * p.z * p.t;
		v[i++] = p.z * p.t * p.t;
		v[i++] = p.t * p.t * p.t;
	}
	if(i == m_order) return v;

	// fourth degree
	if(m_dim >= 1) v[i++] = p.x * p.x * p.x * p.x;
	if(m_dim >= 2) {
		v[i++] = p.x * p.x * p.x * p.y;
		v[i++] = p.x * p.x * p.y * p.y;
		v[i++] = p.x * p.y * p.y * p.y;
		v[i++] = p.y * p.y * p.y * p.y;
	}
	if(m_dim >= 3) {
		v[i++] = p.x * p.x * p.x * p.z;
		v[i++] = p.x * p.x * p.y * p.z;
		v[i++] = p.x * p.x * p.z * p.z;
		v[i++] = p.x * p.y * p.y * p.z;
		v[i++] = p.x * p.y * p.z * p.z;
		v[i++] = p.x * p.z * p.z * p.z;
		v[i++] = p.y * p.y * p.y * p.z;
		v[i++] = p.y * p.y * p.z * p.z;
		v[i++] = p.y * p.z * p.z * p.z;
		v[i++] = p.z * p.z * p.z * p.z;
	}
	if(m_dim >= 4) {
		v[i++] = p.x * p.x * p.x * p.t;
		v[i++] = p.x * p.x * p.y * p.t;
		v[i++] = p.x * p.x * p.z * p.t;
		v[i++] = p.x * p.x * p.t * p.t;
		v[i++] = p.x * p.y * p.y * p.t;
		v[i++] = p.x * p.y * p.z * p.t;
		v[i++] = p.x * p.y * p.t * p.t;
		v[i++] = p.x * p.z * p.z * p.t;
		v[i++] = p.x * p.z * p.t * p.t;
		v[i++] = p.x * p.t * p.t * p.t;
		v[i++] = p.y * p.y * p.y * p.t;
		v[i++] = p.y * p.y * p.z * p.t;
		v[i++] = p.y * p.y * p.t * p.t;
		v[i++] = p.y * p.z * p.z * p.t;
		v[i++] = p.y * p.z * p.t * p.t;
		v[i++] = p.y * p.t * p.t * p.t;
		v[i++] = p.z * p.z * p.z * p.t;
		v[i++] = p.z * p.z * p.t * p.t;
		v[i++] = p.z * p.t * p.t * p.t;
		v[i++] = p.t * p.t * p.t * p.t;
	}
	if(i == m_order) return v;

	// fifth degree
	if(m_dim >= 1) v[i++] = p.x * p.x * p.x * p.x * p.x;
	if(m_dim >= 2) {
		v[i++] = p.x * p.x * p.x * p.x * p.y;
		v[i++] = p.x * p.x * p.x * p.y * p.y;
		v[i++] = p.x * p.x * p.y * p.y * p.y;
		v[i++] = p.x * p.y * p.y * p.y * p.y;
		v[i++] = p.y * p.y * p.y * p.y * p.y;
	}
	if(m_dim >= 3) {
		v[i++] = p.x * p.x * p.x * p.x * p.z;
		v[i++] = p.x * p.x * p.x * p.y * p.z;
		v[i++] = p.x * p.x * p.x * p.z * p.z;
		v[i++] = p.x * p.x * p.y * p.y * p.z;
		v[i++] = p.x * p.x * p.y * p.z * p.z;
		v[i++] = p.x * p.x * p.z * p.z * p.z;
		v[i++] = p.x * p.y * p.y * p.y * p.z;
		v[i++] = p.x * p.y * p.y * p.z * p.z;
		v[i++] = p.x * p.y * p.z * p.z * p.z;
		v[i++] = p.x * p.z * p.z * p.z * p.z;
		v[i++] = p.y * p.y * p.y * p.y * p.z;
		v[i++] = p.y * p.y * p.y * p.z * p.z;
		v[i++] = p.y * p.y * p.z * p.z * p.z;
		v[i++] = p.y * p.z * p.z * p.z * p.z;
		v[i++] = p.z * p.z * p.z * p.z * p.z;
	}
	if(m_dim >= 4) {
		v[i++] = p.x * p.x * p.x * p.x * p.t;
		v[i++] = p.x * p.x * p.x * p.y * p.t;
		v[i++] = p.x * p.x * p.x * p.z * p.t;
		v[i++] = p.x * p.x * p.x * p.t * p.t;
		v[i++] = p.x * p.x * p.y * p.y * p.t;
		v[i++] = p.x * p.x * p.y * p.z * p.t;
		v[i++] = p.x * p.x * p.y * p.t * p.t;
		v[i++] = p.x * p.x * p.z * p.z * p.t;
		v[i++] = p.x * p.x * p.z * p.t * p.t;
		v[i++] = p.x * p.x * p.t * p.t * p.t;
		v[i++] = p.x * p.y * p.y * p.y * p.t;
		v[i++] = p.x * p.y * p.y * p.z * p.t;
		v[i++] = p.x * p.y * p.y * p.t * p.t;
		v[i++] = p.x * p.y * p.z * p.z * p.t;
		v[i++] = p.x * p.y * p.z * p.t * p.t;
		v[i++] = p.x * p.y * p.t * p.t * p.t;
		v[i++] = p.x * p.z * p.z * p.z * p.t;
		v[i++] = p.x * p.z * p.z * p.t * p.t;
		v[i++] = p.x * p.z * p.t * p.t * p.t;
		v[i++] = p.x * p.t * p.t * p.t * p.t;
		v[i++] = p.y * p.y * p.y * p.y * p.t;
		v[i++] = p.y * p.y * p.y * p.z * p.t;
		v[i++] = p.y * p.y * p.y * p.t * p.t;
		v[i++] = p.y * p.y * p.z * p.z * p.t;
		v[i++] = p.y * p.y * p.z * p.t * p.t;
		v[i++] = p.y * p.y * p.t * p.t * p.t;
		v[i++] = p.y * p.z * p.z * p.z * p.t;
		v[i++] = p.y * p.z * p.z * p.t * p.t;
		v[i++] = p.y * p.z * p.t * p.t * p.t;
		v[i++] = p.y * p.t * p.t * p.t * p.t;
		v[i++] = p.z * p.z * p.z * p.z * p.t;
		v[i++] = p.z * p.z * p.z * p.t * p.t;
		v[i++] = p.z * p.z * p.t * p.t * p.t;
		v[i++] = p.z * p.t * p.t * p.t * p.t;
		v[i++] = p.t * p.t * p.t * p.t * p.t;
	}
	return v;
}

double BlockInterpolator::getWeight(const Vector4 &p) const
{
	double wei = 1.0;

	if(m_dim < 1) return wei;
	const double px = fabs(p.x / m_d.x);
	if(px >= 1.0) return 0.0;
	wei *= 1.0 + px * px * (2.0 * px - 3.0);

	if(m_dim < 2) return wei;
	const double py = fabs(p.y / m_d.y);
	if(py >= 1.0) return 0.0;
	wei *= 1.0 + py * py * (2.0 * py - 3.0);

	if(m_dim < 3) return wei;
	const double pz = fabs(p.z / m_d.z);
	if(pz >= 1.0) return 0.0;
	wei *= 1.0 + pz * pz * (2.0 * pz - 3.0);

	if(m_dim < 4) return wei;
	const double pt = fabs(p.t / m_d.t);
	if(pt >= 1.0) return 0.0;
	wei *= 1.0 + pt * pt * (2.0 * pt - 3.0);

	return wei;
}


