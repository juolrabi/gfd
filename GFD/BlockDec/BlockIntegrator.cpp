#include "BlockIntegrator.hpp"
#include "../Discrete/MeshIntegrator.hpp"
#include <iostream>
using namespace std;
using namespace gfd;

BlockIntegrator::BlockIntegrator() {
	m_values = 0;
}

void BlockIntegrator::clear() {
	m_q.clear();
	m_ext.clear();
	m_values = 0;
}

void BlockIntegrator::init(const FormGrade grade, const PartMesh &mesh, const int num) {
	clear();
	MeshIntegrator intg(mesh, grade, num, 0, mesh.getDimension());
	m_values = intg.getLocals();
	if(FormGradeIsDual(grade)) m_ext = intg.getExternals();
	m_q.resize(m_values + m_ext.size());
	for(uint i=0; i<m_q.size(); i++) {
		intg.getQuadrature(i, m_q[i]);
	}
}

void BlockIntegrator::initWedge(const FormGrade grade, const PartMesh &mesh, const int num) {
	clear();
	MeshIntegrator intg0(mesh, grade, num, 0, mesh.getDimension());
	MeshIntegrator intg1(mesh, FormGradeDual(grade), num, 0, mesh.getDimension());
	m_values = intg0.getLocals();
	m_ext = intg0.getExternals();
	m_q.resize(m_values + m_ext.size());
	for(uint i=0; i<m_q.size(); i++) {
		Quadrature q0;
		intg0.getQuadrature(i, q0);
		q0.relocate(-intg0.getPosition(i));
		Quadrature q1;
		intg1.getQuadrature(i, q1);		
		m_q[i].init(q0, q1);
	}

}


