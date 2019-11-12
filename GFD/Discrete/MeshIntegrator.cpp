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

void MeshIntegrator::getSetter(const uint i, Buffer<double> &q) const {
	Buffer<Vector4> p;
	uint ps = 0;
	uint vs = 0;
	if(m_grade == fg_prim0) {
		uint qs = 0;
		q.gather(1.0, qs);
		if(m_num > 0) gatherQuadrature(m_mesh.getNodePosition(i), 1.0, q, qs);
		q.resize(qs);
	}
	else if(m_grade == fg_dual0) {
		Buffer<FourVector4> v;
		if(m_num == 0) v.gather(m_mesh.getNodeDualVector(i), vs);
		else m_mesh.gatherNodeDualSimplices(i, p, ps, v, vs);
		createQuadrature(q, p, ps, v, vs, true);
	}
	else if(m_grade == fg_prim1) {
		Buffer<Vector4> v;
		if(m_num == 0) v.gather(m_mesh.getEdgeVector(i), vs);
		else m_mesh.gatherEdgeSimplices(i, p, ps, v, vs);
		createQuadrature(q, p, ps, v, vs, false);
	}
	else if(m_grade == fg_dual1) {
		Buffer<ThreeVector4> v;
		if(m_num == 0) v.gather(m_mesh.getEdgeDualVector(i), vs);
		else m_mesh.gatherEdgeDualSimplices(i, p, ps, v, vs);
		createQuadrature(q, p, ps, v, vs, true);
	}
	else if(m_grade == fg_prim2) {
		Buffer<TwoVector4> v;
		if(m_num == 0) v.gather(m_mesh.getFaceVector(i), vs);
		else m_mesh.gatherFaceSimplices(i, p, ps, v, vs);
		createQuadrature(q, p, ps, v, vs, false);
	}
	else if(m_grade == fg_dual2) {
		Buffer<TwoVector4> v;
		if(m_num == 0) v.gather(m_mesh.getFaceDualVector(i), vs);
		else m_mesh.gatherFaceDualSimplices(i, p, ps, v, vs);
		createQuadrature(q, p, ps, v, vs, true);
	}
	else if(m_grade == fg_prim3) {
		Buffer<ThreeVector4> v;
		if(m_num == 0) v.gather(m_mesh.getBodyVector(i), vs);
		else m_mesh.gatherBodySimplices(i, p, ps, v, vs);
		createQuadrature(q, p, ps, v, vs, false);
	}
	else if(m_grade == fg_dual3) {
		Buffer<Vector4> v;
		if(m_num == 0) v.gather(m_mesh.getBodyDualVector(i), vs);
		else m_mesh.gatherBodyDualSimplices(i, p, ps, v, vs);
		createQuadrature(q, p, ps, v, vs, true);
	}
	else if(m_grade == fg_prim4) {
		Buffer<FourVector4> v;
		if(m_num == 0) v.gather(m_mesh.getQuadVector(i), vs);
		else m_mesh.gatherQuadSimplices(i, p, ps, v, vs);
		createQuadrature(q, p, ps, v, vs, false);
	}
	else if(m_grade == fg_dual4) {
		uint qs = 0;
		q.gather(m_mesh.getQuadDualVector(i), qs);
		if(m_num > 0) gatherQuadrature(m_mesh.getQuadPosition(i), 1.0, q, qs);
		q.resize(qs);
	}
}

void MeshIntegrator::getWedgeSetter(const uint i, Buffer<double> &q) const {
	Buffer<Vector4> p;
	Buffer<Vector4> dp;
	uint ps = 0;
	uint vs = 0;
	uint dps = 0;
	uint dvs = 0;
	Buffer<double> prim;
	Buffer<double> dual;
	const uint gdim = FormGradeDimension(m_grade);
	if(gdim == 0) {
		Buffer<FourVector4> v;
		if(m_num == 0) v.gather(m_mesh.getNodeDualVector(i), vs);
		else m_mesh.gatherNodeDualSimplices(i, p, ps, v, vs);
		createQuadrature(q, p, ps, v, vs, true);
	}
	else if(gdim == 1) {
		Buffer<Vector4> v;
		if(m_num == 0) v.gather(m_mesh.getEdgeVector(i), vs);
		else m_mesh.gatherEdgeSimplices(i, p, ps, v, vs);
		createQuadrature(prim, p, ps, v, vs, false);

		Buffer<ThreeVector4> dv;
		if(m_num == 0) dv.gather(m_mesh.getEdgeDualVector(i), dvs);
		else m_mesh.gatherEdgeDualSimplices(i, dp, dps, dv, dvs);
		createQuadrature(dual, dp, dps, dv, dvs, true);

		createWedgeQuadrature(q, prim, dual, m_mesh.getEdgePosition(i));
	}
	else if(gdim == 2) {
		Buffer<TwoVector4> v;
		if(m_num == 0) v.gather(m_mesh.getFaceVector(i), vs);
		else m_mesh.gatherFaceSimplices(i, p, ps, v, vs);
		createQuadrature(prim, p, ps, v, vs, false);

		Buffer<TwoVector4> dv;
		if(m_num == 0) dv.gather(m_mesh.getFaceDualVector(i), dvs);
		else m_mesh.gatherFaceDualSimplices(i, dp, dps, dv, dvs);
		createQuadrature(dual, dp, dps, dv, dvs, true);

		createWedgeQuadrature(q, prim, dual, m_mesh.getFacePosition(i));
	}
	else if(gdim == 3) {
		Buffer<ThreeVector4> v;
		if(m_num == 0) v.gather(m_mesh.getBodyVector(i), vs);
		else m_mesh.gatherBodySimplices(i, p, ps, v, vs);
		createQuadrature(prim, p, ps, v, vs, false);

		Buffer<Vector4> dv;
		if(m_num == 0) dv.gather(m_mesh.getBodyDualVector(i), dvs);
		else m_mesh.gatherBodyDualSimplices(i, dp, dps, dv, dvs);
		createQuadrature(dual, dp, dps, dv, dvs, true);

		createWedgeQuadrature(q, prim, dual, m_mesh.getBodyPosition(i));
	}
	else {
		Buffer<FourVector4> v;
		if(m_num == 0) v.gather(m_mesh.getQuadVector(i), vs);
		else m_mesh.gatherQuadSimplices(i, p, ps, v, vs);
		createQuadrature(q, p, ps, v, vs, true);
		if(q[0] < 0.0) q[0] = -q[0];
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
template<typename V> void MeshIntegrator::createQuadrature(Buffer<double> &q, const Buffer<Vector4> &p, const uint ps, const Buffer<V> &v, const uint vs, const bool dual) const {
	if(vs == 0) {
		q.resize(getFields());
		q.fill(0.0);
		return;
	}

	// compute vector
	uint qs = 0;
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
	q.resize(qs);
}
void MeshIntegrator::createWedgeQuadrature(Buffer<double> &q, const Buffer<double> &prim, const Buffer<double> &dual, const Vector4 &p0) const {
	uint i, j, k, l;

	const uint fields = getFields();
	const uint wfields = getWedgeFields();
	const uint dims = m_mesh.getDimension() + 1;
	const uint prims = (prim.size() - fields) / dims;
	const uint duals = (dual.size() - fields) / dims;
	q.resize(wfields + prims * duals * dims);
	for(i=0, k=0; i<fields; i++)
	{
		for(j=0; j<i; j++) q[k++] = prim[i] * dual[j] + prim[j] * dual[i];
		q[k++] = prim[i] * dual[i];
	}

	VectorN p(p0);
	for(i=0, k=0; i<prims; i++)
	{
		const uint ii = fields + dims * i;
		for(j=0; j<duals; j++, k++)
		{
			const uint jj = fields + dims * j;
			const uint kk = wfields + dims * k;
			q[kk] = prim[ii] * dual[jj];
			for(l=1; l<dims; l++) q[kk+l] = prim[ii+l] + dual[jj+l] - p[l-1];
		}
	}
}

