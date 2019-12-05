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
	for(uint i=0; i<m_setter.size(); i++) {
		uint setters = 0;
		intg.gatherSetter(i, m_setter[i], setters);
		m_setter[i].resize(setters);
	}
}

void BlockIntegrator::initWedge(const FormGrade grade, const PartMesh &mesh, const uint num) {
	clear();
	MeshIntegrator intg(mesh, grade, num, 0, mesh.getDimension());

	m_dim = mesh.getDimension();
	m_fields = intg.getWedgeFields();
	m_values = intg.getLocals();
	m_ext = intg.getExternals();
	m_setter.resize(m_values + m_ext.size());
	for(uint i=0; i<m_setter.size(); i++) {
		uint setters = 0;
		intg.gatherWedgeSetter(i, m_setter[i], setters);
		m_setter[i].resize(setters);
	}
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

