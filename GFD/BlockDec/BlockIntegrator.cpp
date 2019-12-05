#include "BlockIntegrator.hpp"
#include "../Discrete/MeshIntegrator.hpp"
#include <iostream>
using namespace std;
using namespace gfd;

BlockIntegrator::BlockIntegrator() {
	m_dim = 1;
	m_fields = 1;
	m_values = 0;
}

void BlockIntegrator::clear() {
	m_setter.clear();
	m_ext.clear();
	m_dim = 1;
	m_fields = 1;
	m_values = 0;
}

void BlockIntegrator::init(const FormGrade grade, const PartMesh &mesh, const uint num) {
	clear();
	MeshIntegrator intg(mesh, grade, num, 0, mesh.getDimension());
	m_dim = mesh.getDimension();
	m_fields = intg.getFields();
	m_values = intg.getLocals();
	if(FormGradeIsDual(grade)) m_ext = intg.getExternals();
	m_setter.resize(m_values + m_ext.size());
	for(uint i=0; i<m_setter.size(); i++) intg.getQuadrature(i, m_setter[i]);
}

void BlockIntegrator::initWedge(const FormGrade grade, const PartMesh &mesh, const uint num) {
	clear();
	MeshIntegrator intg0(mesh, grade, num, 0, mesh.getDimension());
	MeshIntegrator intg1(mesh, FormGradeDual(grade), num, 0, mesh.getDimension());
	m_dim = mesh.getDimension();
	const uint fields = intg0.getFields();
	m_fields = fields * (fields + 1) / 2;
	m_values = intg0.getLocals();
	m_ext = intg0.getExternals();
	m_setter.resize(m_values + m_ext.size());
	const uint size = m_dim + fields;
	for(uint i=0; i<m_setter.size(); i++) {
		Buffer<double> q0;
		intg0.getQuadrature(i, q0);
		Buffer<double> q1;
		intg1.getQuadrature(i, q1);
		const Vector4 p = intg0.getPosition(i);
		Buffer<double> &q = m_setter[i];
		q.resize((q0.size() / size) * (q1.size() / size) * (m_dim + m_fields));
		uint n = 0;
		for(uint j0=0; j0<q0.size(); j0+=size) {
			const double *v0 = &q0[j0 + m_dim];
			for(uint j1=0; j1<q1.size(); j1+=size) {
				const double *v1 = &q1[j1 + m_dim];
				q[n++] = q0[j0] + q1[j1] - p.x;
				if(m_dim >= 2) q[n++] = q0[j0+1] + q1[j1+1] - p.y;
				if(m_dim >= 3) q[n++] = q0[j0+2] + q1[j1+2] - p.z;
				if(m_dim >= 4) q[n++] = q0[j0+3] + q1[j1+3] - p.t;
				for(uint k=0; k<fields; k++) {
					for(uint l=0; l<k; l++) q[n++] = v0[k] * v1[l] + v0[l] * v1[k];
					q[n++] = v0[k] * v1[k];
				}
			}
		}
	}
}

void BlockIntegrator::getVectors(double *result, const Buffer<double *> &exterm) const {
	for(uint i=0; i<m_setter.size(); i++) {
		double *ival = (i < m_values ? &result[m_fields * i] : exterm[i - m_values]);
		const Buffer<double> &setteri = m_setter[i];
		for(uint j=m_dim; j<setteri.size(); j+=m_dim+m_fields) {
			const double *setterij = &setteri[j];
			for(uint k=0; k<m_fields; k++) ival[k] += setterij[k];
		}
	}
}

