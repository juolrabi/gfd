#include "MeshIntegrator.hpp"

using namespace gfd;

MeshIntegrator::MeshIntegrator(const PartMesh &mesh, const FormGrade grade, const uint num) : m_mesh(mesh) {
	m_grade = grade;
	m_num = num;
}

uint MeshIntegrator::getFields() const { 
	const uint gdim = FormGradeDimension(m_grade);
	if(gdim == 0 || gdim == 4) return 1;
	const uint dim = m_mesh.getDimension();
	if(gdim == 1) return dim;
	if(gdim == 2) return dim * (dim - 1) / 2;
	return dim * (dim - 1) * (dim - 2) / 6;
}

uint MeshIntegrator::getWedgeFields() const { 
	const uint fields = getFields();
	return fields * (fields + 1) / 2; 
}

uint MeshIntegrator::getLocals() const {
	const uint gdim = FormGradeDimension(m_grade);
	if(gdim == 0) return m_mesh.getNodeLocals();
	if(gdim == 1) return m_mesh.getEdgeLocals();
	if(gdim == 2) return m_mesh.getFaceLocals();
	if(gdim == 3) return m_mesh.getBodyLocals();
	return m_mesh.getQuadLocals();
}

const Buffer< pair<uint,uint> > &MeshIntegrator::getExternals() const {
	const uint gdim = FormGradeDimension(m_grade);
	if(gdim == 0) return m_mesh.getExternalNodes();
	if(gdim == 1) return m_mesh.getExternalEdges();
	if(gdim == 2) return m_mesh.getExternalFaces();
	if(gdim == 3) return m_mesh.getExternalBodies();
	return m_mesh.getExternalQuads();
}

void MeshIntegrator::gatherSetter0(const uint i, Buffer<double> &q, uint &qs) const {
	if(m_grade == fg_prim0) q.gather(1.0, qs);
	else if(m_grade == fg_dual0) gatherVector(m_mesh.getNodeDualVector(i).dualof(), q, qs);
	else if(m_grade == fg_prim1) gatherVector(m_mesh.getEdgeVector(i), q, qs);
	else if(m_grade == fg_dual1) gatherVector(m_mesh.getEdgeDualVector(i).dualof(), q, qs);
	else if(m_grade == fg_prim2) gatherVector(m_mesh.getFaceVector(i), q, qs);
	else if(m_grade == fg_dual2) gatherVector(m_mesh.getFaceDualVector(i).dualof(), q, qs);
	else if(m_grade == fg_prim3) gatherVector(m_mesh.getBodyVector(i), q, qs);
	else if(m_grade == fg_dual3) gatherVector(m_mesh.getBodyDualVector(i).dualof(), q, qs);
	else if(m_grade == fg_prim4) gatherVector(m_mesh.getQuadVector(i), q, qs);
	else if(m_grade == fg_dual4) q.gather(m_mesh.getQuadDualVector(i), qs);
}

void MeshIntegrator::gatherWedgeSetter0(const uint i, Buffer<double> &q, uint &qs) const {
	uint prims = 0;
	Buffer<double> prim;
	uint duals = 0;
	Buffer<double> dual;
	const uint gdim = FormGradeDimension(m_grade);
	if(gdim == 0) gatherVector(m_mesh.getNodeDualVector(i).dualof(), q, qs);
	else if(gdim == 1) {
		gatherVector(m_mesh.getEdgeVector(i), prim, prims);
		gatherVector(m_mesh.getEdgeDualVector(i).dualof(), dual, duals);
		gatherWedgeQuadrature(prim, prims, dual, duals, Vector4(0,0,0,0), q, qs);
	} 
	else if(gdim == 2) {
		gatherVector(m_mesh.getFaceVector(i), prim, prims);
		gatherVector(m_mesh.getFaceDualVector(i).dualof(), dual, duals);
		gatherWedgeQuadrature(prim, prims, dual, duals, Vector4(0,0,0,0), q, qs);
	}
	else if(gdim == 3) {
		gatherVector(m_mesh.getBodyVector(i), prim, prims);
		gatherVector(m_mesh.getBodyDualVector(i).dualof(), dual, duals);
		gatherWedgeQuadrature(prim, prims, dual, duals, Vector4(0,0,0,0), q, qs);
	}
	else if(gdim == 4) {
		FourVector4 v = m_mesh.getQuadVector(i);
		if(v.dual() < 0.0) v = -v; // divide with QuadDualVector
		gatherVector(v, q, qs);
	}
}

void MeshIntegrator::gatherSetter(const uint i, Buffer<double> &q, uint &qs) const {
	if(m_num == 0) {
		gatherSetter0(i, q, qs);
		return;
	}
	Buffer<Vector4> p;
	uint ps = 0;
	uint vs = 0;
	if(m_grade == fg_prim0) {
		q.gather(1.0, qs);
		gatherQuadrature(m_mesh.getNodePosition(i), 1.0, q, qs);
	}
	else if(m_grade == fg_dual0) {
		Buffer<FourVector4> v;
		m_mesh.gatherNodeDualSimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, true, q, qs);
	}
	else if(m_grade == fg_prim1) {
		Buffer<Vector4> v;
		m_mesh.gatherEdgeSimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, false, q, qs);
	}
	else if(m_grade == fg_dual1) {
		Buffer<ThreeVector4> v;
		m_mesh.gatherEdgeDualSimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, true, q, qs);
	}
	else if(m_grade == fg_prim2) {
		Buffer<TwoVector4> v;
		m_mesh.gatherFaceSimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, false, q, qs);
	}
	else if(m_grade == fg_dual2) {
		Buffer<TwoVector4> v;
		m_mesh.gatherFaceDualSimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, true, q, qs);
	}
	else if(m_grade == fg_prim3) {
		Buffer<ThreeVector4> v;
		m_mesh.gatherBodySimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, false, q, qs);
	}
	else if(m_grade == fg_dual3) {
		Buffer<Vector4> v;
		m_mesh.gatherBodyDualSimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, true, q, qs);
	}
	else if(m_grade == fg_prim4) {
		Buffer<FourVector4> v;
		m_mesh.gatherQuadSimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, false, q, qs);
	}
	else if(m_grade == fg_dual4) {
		q.gather(m_mesh.getQuadDualVector(i), qs);
		gatherQuadrature(m_mesh.getQuadPosition(i), 1.0, q, qs);
	}
}

void MeshIntegrator::gatherWedgeSetter(const uint i, Buffer<double> &q, uint &qs) const {
	if(m_num == 0) {
		gatherWedgeSetter0(i, q, qs);
		return;
	}
	
	uint ps = 0;
	uint vs = 0;
	Buffer<Vector4> p;
	uint dps = 0;
	uint dvs = 0;
	Buffer<Vector4> dp;
	uint prims = 0;
	Buffer<double> prim;
	uint duals = 0;
	Buffer<double> dual;
	const uint gdim = FormGradeDimension(m_grade);
	if(gdim == 0) {
		Buffer<FourVector4> v;
		m_mesh.gatherNodeDualSimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, true, q, qs);
	}
	else if(gdim == 1) {
		Buffer<Vector4> v;
		m_mesh.gatherEdgeSimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, false, prim, prims);

		Buffer<ThreeVector4> dv;
		m_mesh.gatherEdgeDualSimplices(i, dp, dps, dv, dvs);
		gatherQuadrature(dp, dps, dv, dvs, true, dual, duals);

		gatherWedgeQuadrature(prim, prims, dual, duals, m_mesh.getEdgePosition(i), q, qs);
	}
	else if(gdim == 2) {
		Buffer<TwoVector4> v;
		m_mesh.gatherFaceSimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, false, prim, prims);

		Buffer<TwoVector4> dv;
		m_mesh.gatherFaceDualSimplices(i, dp, dps, dv, dvs);
		gatherQuadrature(dp, dps, dv, dvs, true, dual, duals);

		gatherWedgeQuadrature(prim, prims, dual, duals, m_mesh.getFacePosition(i), q, qs);
	}
	else if(gdim == 3) {
		Buffer<ThreeVector4> v;
		m_mesh.gatherBodySimplices(i, p, ps, v, vs);
		gatherQuadrature(p, ps, v, vs, false, prim, prims);

		Buffer<Vector4> dv;
		m_mesh.gatherBodyDualSimplices(i, dp, dps, dv, dvs);
		gatherQuadrature(dp, dps, dv, dvs, true, dual, duals);

		gatherWedgeQuadrature(prim, prims, dual, duals, m_mesh.getBodyPosition(i), q, qs);
	}
	else {
		Buffer<FourVector4> v;
		m_mesh.gatherQuadSimplices(i, p, ps, v, vs);
		FourVector4 sum(0.0);
		for(uint i=0; i<vs; i++) sum += v[i];
		if(sum.dual() < 0.0) {
			for(uint i=0; i<vs; i++) v[i] = -v[i];
		}
		gatherQuadrature(p, ps, v, vs, false, q, qs);
	}
}

void MeshIntegrator::gatherQuadrature(const Vector4 &p, const double w, Buffer<double> &q, uint &qs) const {
	q.gather(w, qs);
	gatherVector(p, q, qs);
}
void MeshIntegrator::gatherQuadrature(const Buffer<Vector4> &p, const double w, Buffer<double> &q, uint &qs) const {
	if(p.size() == 1) gatherQuadrature(p[0], w, q, qs);
	else if(p.size() == 2) {
		uint n = uint(m_num * fabs(w) + 0.5);
		if(n < 2) n = 2;
		const double div = sqrt((n-1) * (n+1));
		const Vector4 v1 = (p[1] - p[0]) / div;
		const Vector4 d0 = p[0] + (div + 1.0 - n) / 2.0 * v1;
		const double dw = w / double(n);
		for(uint i=0; i<n; i++) gatherQuadrature(d0 + i * v1, dw, q, qs);
	}
	else if(p.size() == 3) {
		uint n = uint(m_num * sqrt(fabs(w)) + 0.5);
		if(n < 2) n = 2;
		const double div = sqrt((n-1) * (n+2));
		const Vector4 v1 = (p[1] - p[0]) / div;
		const Vector4 v2 = (p[2] - p[0]) / div;
		const Vector4 d0 = p[0] + (div + 1.0 - n) / 3.0 * (v1 + v2);
		const double dw = 2.0 * w / (double(n) * double(n + 1));
		for(uint i=0; i<n; i++) {
			const Vector4 d1 = d0 + i * v1;
			for(uint j=0; j<n-i; j++) gatherQuadrature(d1 + j * v2, dw, q, qs);
		}
	}
	else if(p.size() == 4) {
		uint n = uint(m_num * pow(fabs(w), 1.0 / 3.0) + 0.5);
		if(n < 2) n = 2;
		const double div = sqrt((n-1) * (n+3));
		const Vector4 v1 = (p[1] - p[0]) / div;
		const Vector4 v2 = (p[2] - p[0]) / div;
		const Vector4 v3 = (p[3] - p[0]) / div;
		const Vector4 d0 = p[0] + (div + 1.0 - n) / 4.0 * (v1 + v2 + v3);
		const double dw = 6.0 * w / (double(n) * double(n + 1) * double(n + 2));
		for(uint i=0; i<n; i++) {
			const Vector4 d1 = d0 + i * v1;
			for(uint j=0; j<n-i; j++) {
				const Vector4 d2 = d1 + j * v2;
				for(uint k=0; k<n-i-j; k++) gatherQuadrature(d2 + k * v3, dw, q, qs);
			}
		}
	}
	else if(p.size() == 5) {
		uint n = uint(m_num * pow(fabs(w), 0.25) + 0.5);
		if(n < 2) n = 2;
		const double div = sqrt((n-1) * (n+4));
		const Vector4 v1 = (p[1] - p[0]) / div;
		const Vector4 v2 = (p[2] - p[0]) / div;
		const Vector4 v3 = (p[3] - p[0]) / div;
		const Vector4 v4 = (p[4] - p[0]) / div;
		const Vector4 d0 = p[0] + (div + 1.0 - n) / 5.0 * (v1 + v2 + v3 + v4);
		const double dw = 24.0 * w / (double(n) * double(n + 1) * double(n + 2) * double(n + 3));
		for(uint i=0; i<n; i++) {
			const Vector4 d1 = d0 + i * v1;
			for(uint j=0; j<n-i; j++) {
				const Vector4 d2 = d1 + j * v2;
				for(uint k=0; k<n-i-j; k++) {
					const Vector4 d3 = d2 + k * v3;
					for(uint l=0; l<n-i-j-k; l++) gatherQuadrature(d3 + l * v4, dw, q, qs);
				}
			}
		}
	}
}
void MeshIntegrator::gatherVector(const Vector4 &v, Buffer<double> &q, uint &qs) const {
	const uint dim = m_mesh.getDimension();
	q.gather(v.x, qs);
	if(dim == 1) return;
	q.gather(v.y, qs);
	if(dim == 2) return;
	q.gather(v.z, qs);
	if(dim == 3) return;
	q.gather(v.t, qs);
}
void MeshIntegrator::gatherVector(const TwoVector4 &v, Buffer<double> &q, uint &qs) const {
	const uint dim = m_mesh.getDimension();
	q.gather(v.xy, qs);
	if(dim == 2) return;
	q.gather(v.xz, qs);
	q.gather(v.yz, qs);
	if(dim == 3) return;
	q.gather(v.xt, qs);
	q.gather(v.yt, qs);
	q.gather(v.zt, qs);
}
void MeshIntegrator::gatherVector(const ThreeVector4 &v, Buffer<double> &q, uint &qs) const {
	const uint dim = m_mesh.getDimension();
	q.gather(v.xyz, qs);
	if(dim == 3) return;
	q.gather(v.xyt, qs);
	q.gather(v.xzt, qs);
	q.gather(v.yzt, qs);
}
void MeshIntegrator::gatherVector(const FourVector4 &v, Buffer<double> &q, uint &qs) const {
	q.gather(v.xyzt, qs);
}
template<typename V> void MeshIntegrator::gatherQuadrature(const Buffer<Vector4> &p, const uint ps, const Buffer<V> &v, const uint vs, const bool dual, Buffer<double> &q, uint &qs) const {
	if(vs == 0) {
		const uint qsize = qs + getFields();
		if(qsize > q.size()) q.resize(qsize);
		while(qs < qsize) q[qs++] = 0.0;
		return;
	}
	
	// compute vector
	V sumv(v[0]);
	for(uint i=1; i<vs; i++) sumv += v[i];
	if(dual) gatherVector(sumv.dualof(), q, qs);
	else gatherVector(sumv, q, qs);
	sumv /= sumv.lensq();

	// compute weights and positions
	const uint psize = ps / vs;
	if(m_num == 1) {
		Vector4 sump(0,0,0,0);
		for(uint i=0; i<vs; i++) {
			const double fac = v[i].dot(sumv) / double(psize);
			for(uint j=0; j<psize; j++) sump += fac * p[psize * i + j];
		}
		gatherQuadrature(sump, 1.0, q, qs);
	}
	else if(m_num >= 2) {
		Buffer<Vector4> pp(psize);
		for(uint i=0; i<vs; i++) {
			for(uint j=0; j<psize; j++) pp[j] = p[psize * i + j];
			gatherQuadrature(pp, v[i].dot(sumv), q, qs);
		}
	}
}
void MeshIntegrator::gatherWedgeQuadrature(const Buffer<double> &prim, const uint prims, const Buffer<double> &dual, const uint duals, const Vector4 &p0, Buffer<double> &q, uint &qs) const {
	uint i, j, k;
	const uint fields = getFields();
	const uint dims = m_mesh.getDimension() + 1;
	const uint qsize = qs + getWedgeFields() + (prims - fields) * (duals - fields) / dims;
	if(qsize > q.size()) q.resize(qsize);
	for(i=0; i<fields; i++) {
		for(j=0; j<i; j++) q[qs++] = prim[i] * dual[j] + prim[j] * dual[i];
		q[qs++] = prim[i] * dual[i];
	}
	VectorN p(p0);
	for(i=fields; i<prims; i+=dims) {
		for(j=fields; j<duals; j+=dims) {
			q[qs++] = prim[i] * dual[j];
			for(k=1; k<dims; k++) q[qs++] = prim[i+k] + dual[j+k] - p[k-1];
		}
	}
}

