#include "BlockIntegrator.hpp"
#include <iostream>
using namespace std;
using namespace gfd;

BlockIntegrator::BlockIntegrator()
{
	m_dim = 1;
	m_fields = 1;
	m_values = 0;
}

void BlockIntegrator::clear()
{
	m_setter.clear();
	m_ext.clear();
	m_dim = 1;
	m_fields = 1;
	m_values = 0;
}

void BlockIntegrator::gatherQuadrature(const Vector4 &p, const double w, Buffer<double> &q, uint &qs) const
{
	q.gather(w, qs);
	gatherVector(p, q, qs);
}
void BlockIntegrator::gatherQuadrature(const Buffer<Vector4> &p, const double w, const uint num, Buffer<double> &q, uint &qs) const
{
	if(p.size() == 1) gatherQuadrature(p[0], w, q, qs);
	else if(p.size() == 2)
	{
		uint n = uint(num * fabs(w) + 0.5);
		if(n < 2) n = 2;
		const double div = sqrt((n-1) * (n+1));
		const Vector4 v1 = (p[1] - p[0]) / div;
		const Vector4 d0 = p[0] + (div + 1.0 - n) / 2.0 * v1;
		const double dw = w / double(n);
		for(uint i=0; i<n; i++) gatherQuadrature(d0 + i * v1, dw, q, qs);
	}
	else if(p.size() == 3)
	{
		uint n = uint(num * sqrt(fabs(w)) + 0.5);
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
	else if(p.size() == 4)
	{
		uint n = uint(num * pow(fabs(w), 1.0 / 3.0) + 0.5);
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
	else if(p.size() == 5)
	{
		uint n = uint(num * pow(fabs(w), 0.25) + 0.5);
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
void BlockIntegrator::gatherVector(const Vector4 &v, Buffer<double> &q, uint &qs) const
{
	q.gather(v.x, qs);
	if(m_dim == 1) return;
	q.gather(v.y, qs);
	if(m_dim == 2) return;
	q.gather(v.z, qs);
	if(m_dim == 3) return;
	q.gather(v.t, qs);
}
void BlockIntegrator::gatherVector(const TwoVector4 &v, Buffer<double> &q, uint &qs) const
{
	q.gather(v.xy, qs);
	if(m_dim == 2) return;
	q.gather(v.xz, qs);
	q.gather(v.yz, qs);
	if(m_dim == 3) return;
	q.gather(v.xt, qs);
	q.gather(v.yt, qs);
	q.gather(v.zt, qs);
}
void BlockIntegrator::gatherVector(const ThreeVector4 &v, Buffer<double> &q, uint &qs) const
{
	q.gather(v.xyz, qs);
	if(m_dim == 3) return;
	q.gather(v.xyt, qs);
	q.gather(v.xzt, qs);
	q.gather(v.yzt, qs);
}
void BlockIntegrator::gatherVector(const FourVector4 &v, Buffer<double> &q, uint &qs) const
{
	q.gather(v.xyzt, qs);
}
/*void BlockIntegrator::gatherDualVector(const Vector4 &v, Buffer<double> &q, uint &qs) const
{
	if(m_dim > 3) {
		q.gather(v.x, qs);
		q.gather(v.y, qs);
		q.gather(v.z, qs);
	}
	q.gather(v.t, qs);
}
void BlockIntegrator::gatherDualVector(const TwoVector4 &v, Buffer<double> &q, uint &qs) const
{
	if(m_dim > 3) {
		q.gather(v.xy, qs);
		q.gather(v.xz, qs);
		q.gather(v.yz, qs);
	}
	if(m_dim > 2) {
		q.gather(v.xt, qs);
		q.gather(v.yt, qs);
	}
	q.gather(v.zt, qs);
}
void BlockIntegrator::gatherDualVector(const ThreeVector4 &v, Buffer<double> &q, uint &qs) const
{
	if(m_dim > 3) q.gather(v.xyz, qs);
	if(m_dim > 2) q.gather(v.xyt, qs);
	if(m_dim > 1) q.gather(v.xzt, qs);
	q.gather(v.yzt, qs);
}
void BlockIntegrator::gatherDualVector(const FourVector4 &v, Buffer<double> &q, uint &qs) const
{
	q.gather(v.xyzt, qs);
}
*/
template<typename V> void BlockIntegrator::createQuadrature(Buffer<double> &q, const Buffer<Vector4> &p, const uint ps, const Buffer<V> &v, const uint vs, const uint num, const bool dual) const
{
	if(vs == 0)
	{
		q.resize(m_fields);
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
	if(num == 1)
	{
		Vector4 sump(0,0,0,0);
		for(uint i=0; i<vs; i++)
		{
			const double fac = v[i].dot(sumv) / double(psize);
			for(uint j=0; j<psize; j++) sump += fac * p[psize * i + j];
		}
		gatherQuadrature(sump, 1.0, q, qs);
	}
	else if(num >= 2)
	{
		Buffer<Vector4> pp(psize);
		for(uint i=0; i<vs; i++)
		{
			for(uint j=0; j<psize; j++) pp[j] = p[psize * i + j];
			gatherQuadrature(pp, v[i].dot(sumv), num, q, qs);
		}
	}

	q.resize(qs);
}
void BlockIntegrator::createProductQuadrature(Buffer<double> &q, const Buffer<double> &prim, const Buffer<double> &dual, const Vector4 &p0) const
{
	uint i, j, k, l;

	const uint fields = m_fields * (m_fields + 1) / 2;
	const uint dims = m_dim + 1;
	const uint prims = (prim.size() - m_fields) / dims;
	const uint duals = (dual.size() - m_fields) / dims;
	q.resize(fields + prims * duals * dims);
	for(i=0, k=0; i<m_fields; i++)
	{
		for(j=0; j<i; j++) q[k++] = prim[i] * dual[j] + prim[j] * dual[i];
		q[k++] = prim[i] * dual[i];
	}

	VectorN p(p0);
	for(i=0, k=0; i<prims; i++)
	{
		const uint ii = m_fields + dims * i;
		for(j=0; j<duals; j++, k++)
		{
			const uint jj = m_fields + dims * j;
			const uint kk = fields + dims * k;
			q[kk] = prim[ii] * dual[jj];
			for(l=1; l<dims; l++) q[kk+l] = prim[ii+l] + dual[jj+l] - p[l-1];
		}
	}

/*	const uint fields = uint(sqrt(2.0 * m_fields));
	const uint dims = m_dim + 1;
	const uint prims = (prim.size() - fields) / dims;
	const uint duals = (dual.size() - fields) / dims;
	q.resize(m_fields + prims * duals * dims);
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
			const uint kk = m_fields + dims * k;
			q[kk] = prim[ii] * dual[jj];
			for(l=1; l<dims; l++) q[kk+l] = prim[ii+l] + dual[jj+l] - p[l-1];
		}
	}
*/
}
void BlockIntegrator::init(const FormGrade grade, const PartMesh &mesh, const uint num)
{
	clear();
	uint i;
	m_dim = mesh.getDimension();
	m_fields = getFieldDimension(grade, m_dim);
	if(grade == fg_prim0)
	{
		m_values = mesh.getNodeLocals();
		m_setter.resize(m_values);
		for(i=0; i<m_setter.size(); i++)
		{
			uint qs = 0;
			Buffer<double> &q = m_setter[i];
			q.gather(1.0, qs);
			if(num > 0) gatherQuadrature(mesh.getNodePosition(i), 1.0, q, qs);
			q.resize(qs);
		}
	}
	else if(grade == fg_dual0)
	{
		m_values = mesh.getNodeLocals();
		m_setter.resize(mesh.getNodeSize());
		Buffer<Vector4> p;
		Buffer<FourVector4> v;
		for(i=0; i<m_setter.size(); i++)
		{
			uint ps = 0;
			uint vs = 0;
			if(num == 0) v.gather(mesh.getNodeDualVector(i), vs);
			else mesh.gatherNodeDualSimplices(i, p, ps, v, vs);
			createQuadrature(m_setter[i], p, ps, v, vs, num, true);
		}
		m_ext.resize(m_setter.size() - m_values);
		for(i=0; i<m_ext.size(); i++) m_ext[i] = pair<uint,uint>(mesh.getNodePart(m_values + i), mesh.getNodeLink(m_values + i));
	}
	else if(grade == fg_prim1)
	{
		m_values = mesh.getEdgeLocals();
		m_setter.resize(m_values);
		Buffer<Vector4> p;
		Buffer<Vector4> v;
		for(i=0; i<m_setter.size(); i++)
		{
			uint ps = 0;
			uint vs = 0;
			if(num == 0) v.gather(mesh.getEdgeVector(i), vs);
			else mesh.gatherEdgeSimplices(i, p, ps, v, vs);
			createQuadrature(m_setter[i], p, ps, v, vs, num, false);
		}
	}
	else if(grade == fg_dual1)
	{
		m_values = mesh.getEdgeLocals();
		m_setter.resize(mesh.getEdgeSize());
		Buffer<Vector4> p;
		Buffer<ThreeVector4> v;
		for(i=0; i<m_setter.size(); i++)
		{
			uint ps = 0;
			uint vs = 0;
			if(num == 0) v.gather(mesh.getEdgeDualVector(i), vs);
			else mesh.gatherEdgeDualSimplices(i, p, ps, v, vs);
			createQuadrature(m_setter[i], p, ps, v, vs, num, true);
		}
		m_ext.resize(m_setter.size() - m_values);
		for(i=0; i<m_ext.size(); i++) m_ext[i] = pair<uint,uint>(mesh.getEdgePart(m_values + i), mesh.getEdgeLink(m_values + i));
	}
	else if(grade == fg_prim2)
	{
		m_values = mesh.getFaceLocals();
		m_setter.resize(m_values);
		Buffer<Vector4> p;
		Buffer<TwoVector4> v;
		for(i=0; i<m_setter.size(); i++)
		{
			uint ps = 0;
			uint vs = 0;
			if(num == 0) v.gather(mesh.getFaceVector(i), vs);
			else mesh.gatherFaceSimplices(i, p, ps, v, vs);
			createQuadrature(m_setter[i], p, ps, v, vs, num, false);
		}
	}
	else if(grade == fg_dual2)
	{
		m_values = mesh.getFaceLocals();
		m_setter.resize(mesh.getFaceSize());
		Buffer<Vector4> p;
		Buffer<TwoVector4> v;
		for(i=0; i<m_setter.size(); i++)
		{
			uint ps = 0;
			uint vs = 0;
			if(num == 0) v.gather(mesh.getFaceDualVector(i), vs);
			else mesh.gatherFaceDualSimplices(i, p, ps, v, vs);
			createQuadrature(m_setter[i], p, ps, v, vs, num, true);
		}
		m_ext.resize(m_setter.size() - m_values);
		for(i=0; i<m_ext.size(); i++) m_ext[i] = pair<uint,uint>(mesh.getFacePart(m_values + i), mesh.getFaceLink(m_values + i));
	}
	else if(grade == fg_prim3)
	{
		m_values = mesh.getBodyLocals();
		m_setter.resize(m_values);
		Buffer<Vector4> p;
		Buffer<ThreeVector4> v;
		for(i=0; i<m_setter.size(); i++)
		{
			uint ps = 0;
			uint vs = 0;
			if(num == 0) v.gather(mesh.getBodyVector(i), vs);
			else mesh.gatherBodySimplices(i, p, ps, v, vs);
			createQuadrature(m_setter[i], p, ps, v, vs, num, false);
		}
	}
	else if(grade == fg_dual3)
	{
		m_values = mesh.getBodyLocals();
		m_setter.resize(mesh.getBodySize());
		Buffer<Vector4> p;
		Buffer<Vector4> v;
		for(i=0; i<m_setter.size(); i++)
		{
			uint ps = 0;
			uint vs = 0;
			if(num == 0) v.gather(mesh.getBodyDualVector(i), vs);
			else mesh.gatherBodyDualSimplices(i, p, ps, v, vs);
			createQuadrature(m_setter[i], p, ps, v, vs, num, true);
		}
		m_ext.resize(m_setter.size() - m_values);
		for(i=0; i<m_ext.size(); i++) m_ext[i] = pair<uint,uint>(mesh.getBodyPart(m_values + i), mesh.getBodyLink(m_values + i));
	}
	else if(grade == fg_prim4)
	{
		m_values = mesh.getQuadLocals();
		m_setter.resize(m_values);
		Buffer<Vector4> p;
		Buffer<FourVector4> v;
		for(i=0; i<m_setter.size(); i++)
		{
			uint ps = 0;
			uint vs = 0;
			if(num == 0) v.gather(mesh.getQuadVector(i), vs);
			else mesh.gatherQuadSimplices(i, p, ps, v, vs);
			createQuadrature(m_setter[i], p, ps, v, vs, num, false);
		}
	}
	else if(grade == fg_dual4)
	{
		m_values = mesh.getQuadLocals();
		m_setter.resize(m_values);
		for(i=0; i<m_setter.size(); i++)
		{
			uint qs = 0;
			Buffer<double> &q = m_setter[i];
			q.gather(mesh.getQuadDualVector(i), qs);
			if(num > 0) gatherQuadrature(mesh.getQuadPosition(i), 1.0, q, qs);
			q.resize(qs);
		}
	}
}

void BlockIntegrator::initWedge(const FormGrade grade, const PartMesh &mesh, const uint num)
{
	clear();
	m_dim = mesh.getDimension();
	m_fields = getFieldDimension(grade, m_dim);
	const uint gdim = FormGradeDimension(grade);
	if(gdim == 0) init(fg_dual0, mesh, num);
	else if(gdim == 1)
	{
		m_values = mesh.getEdgeLocals();
		m_setter.resize(mesh.getEdgeSize());
		Buffer<Vector4> p;
		Buffer<Vector4> v;
		Buffer<Vector4> dp;
		Buffer<ThreeVector4> dv;
		for(uint i=0; i<m_setter.size(); i++)
		{
			const Vector4 p0 = mesh.getEdgePosition(i);

			uint ps = 0;
			uint vs = 0;
			if(num == 0) v.gather(mesh.getEdgeVector(i), vs);
			else mesh.gatherEdgeSimplices(i, p, ps, v, vs);
			Buffer<double> prim;
			createQuadrature(prim, p, ps, v, vs, num, false);

			uint dps = 0;
			uint dvs = 0;
			if(num == 0) dv.gather(mesh.getEdgeDualVector(i), dvs);
			else mesh.gatherEdgeDualSimplices(i, dp, dps, dv, dvs);
			Buffer<double> dual;
			createQuadrature(dual, dp, dps, dv, dvs, num, true);

			createProductQuadrature(m_setter[i], prim, dual, p0);
		}
		m_ext.resize(m_setter.size() - m_values);
		for(uint i=0; i<m_ext.size(); i++) m_ext[i] = pair<uint,uint>(mesh.getEdgePart(m_values + i), mesh.getEdgeLink(m_values + i));
	}
	else if(gdim == 2)
	{
		m_values = mesh.getFaceLocals();
		m_setter.resize(mesh.getFaceSize());
		Buffer<Vector4> p;
		Buffer<TwoVector4> v;
		Buffer<Vector4> dp;
		Buffer<TwoVector4> dv;
		for(uint i=0; i<m_setter.size(); i++)
		{
			const Vector4 p0 = mesh.getFacePosition(i);

			uint ps = 0;
			uint vs = 0;
			if(num == 0) v.gather(mesh.getFaceVector(i), vs);
			else mesh.gatherFaceSimplices(i, p, ps, v, vs);
			Buffer<double> prim;
			createQuadrature(prim, p, ps, v, vs, num, false);

			uint dps = 0;
			uint dvs = 0;
			if(num == 0) dv.gather(mesh.getFaceDualVector(i), dvs);
			else mesh.gatherFaceDualSimplices(i, dp, dps, dv, dvs);
			Buffer<double> dual;
			createQuadrature(dual, dp, dps, dv, dvs, num, true);

			createProductQuadrature(m_setter[i], prim, dual, p0);
		}
		m_ext.resize(m_setter.size() - m_values);
		for(uint i=0; i<m_ext.size(); i++) m_ext[i] = pair<uint,uint>(mesh.getFacePart(m_values + i), mesh.getFaceLink(m_values + i));
	}
	else if(gdim == 3)
	{
		m_values = mesh.getBodyLocals();
		m_setter.resize(mesh.getBodySize());
		Buffer<Vector4> p;
		Buffer<ThreeVector4> v;
		Buffer<Vector4> dp;
		Buffer<Vector4> dv;
		for(uint i=0; i<m_setter.size(); i++)
		{
			const Vector4 p0 = mesh.getBodyPosition(i);

			uint ps = 0;
			uint vs = 0;
			if(num == 0) v.gather(mesh.getBodyVector(i), vs);
			else mesh.gatherBodySimplices(i, p, ps, v, vs);
			Buffer<double> prim;
			createQuadrature(prim, p, ps, v, vs, num, false);

			uint dps = 0;
			uint dvs = 0;
			if(num == 0) dv.gather(mesh.getBodyDualVector(i), dvs);
			else mesh.gatherBodyDualSimplices(i, dp, dps, dv, dvs);
			Buffer<double> dual;
			createQuadrature(dual, dp, dps, dv, dvs, num, true);

			createProductQuadrature(m_setter[i], prim, dual, p0);
		}
		m_ext.resize(m_setter.size() - m_values);
		for(uint i=0; i<m_ext.size(); i++) m_ext[i] = pair<uint,uint>(mesh.getBodyPart(m_values + i), mesh.getBodyLink(m_values + i));
	}
	else
	{
		init(fg_prim4, mesh, num);
		for(uint i=0; i<m_setter.size(); i++)
		{
			if(m_setter[i][0] < 0.0) m_setter[i][0] = -m_setter[i][0];
		}
	}
	m_fields = m_fields * (m_fields + 1) / 2;
}

void BlockIntegrator::getVectors(double *result, const Buffer<double *> &exterm) const
{
	for(uint i=0; i<m_setter.size(); i++)
	{
		double *ival = (i < m_values ? &result[m_fields * i] : exterm[i - m_values]);
		const Buffer<double> &setteri = m_setter[i];
		for(uint k=0; k<m_fields; k++) ival[k] += setteri[k];
	}
}


