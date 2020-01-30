#include "TimeIntegrator.hpp"
//#include "../Types/MpiEasy.hpp"

using namespace gfd;

double selectPositive(const double &r) { return (r > 0.0 ? r : 0.0); }
double selectNegative(const double &r) { return (r < 0.0 ? r : 0.0); }

TimeIntegrator::TimeIntegrator(const Buffer< Sparse<double> > &d, const Buffer< Diagonal<double> > &a, const Buffer< Column<double> > &f, const double dtime) {
	uint i, j;
	m_time = 0.0;

	// find number of forms to compute and initialize member buffers
	m_forms = d.size() / 2 + 1;
	m_d.resize(2 * (m_forms - 1));
	m_a.resize(m_forms);
	m_e.resize(m_forms);
	m_v.resize(m_forms);
	m_f.resize(m_forms);

	// check the dimensions of sparse matrices in Buffer d
//	uint sumheight = 0;
	Buffer<uint> height(m_forms);
	for(i=0; i<m_forms; i++) {
		if(i == 0) height[i] = d[0].m_height;
		else {
			const uint upd = 2 * i - 2;
			const uint leftd = 2 * i - 1;
			height[i] = d[leftd].m_height;
			if(d[upd].m_height != height[i-1] || d[upd].m_width != height[i]) {
				cout << "TimeIntegrator::TimeIntegrator -> dimensions of derivative d[" << upd << "] do not match." << endl;
				m_d[upd].setSparse(height[i], height[i-1]);
			}
			else m_d[upd] = d[upd];
			if(d[leftd].m_width != height[i-1]) {
				cout << "TimeIntegrator::TimeIntegrator -> width of derivative d[" << leftd << "] do not match." << endl;
				m_d[leftd].setSparse(height[i-1], height[i]);
			}
			else m_d[leftd] = d[leftd];
		}
//		sumheight += height[i];
	}

	// compute number of time step by CFL-condition
	m_steps = 1;
	for(i=0; i<m_forms; i++) {
		Sparse<double> product;
		if(i == 0) product.setTimes(m_d[0], m_d[1]);
		else {
			product.setTimes(m_d[2 * i - 1], m_d[2 * j - 2]);
			if(i + 1 < m_forms) product += m_d[2 * i] * m_d[2 * i + 1];
		}
		for(j=0; j<height[i]; j++) { 
			const uint steps = uint(dtime * sqrt(-0.5 * product.getValue(j, j)) + 1.0);
			if(m_steps < steps) m_steps = steps;
		}
	}
	maxMPI(&m_steps, 1);
	m_dtime = dtime / double(m_steps);


	// set integration terms
	for(i=0; i<m_forms; i++) {
		if(a.size() > i && a[i].m_height == height[i]) m_a[i].setNegation(a[i]);
		else m_a[i].setSparse(height[i]);
		m_e[i].setFunction(m_a[i], selectPositive);
		m_e[i].trimSparse();
		m_a[i].setFunction(m_a[i], selectNegative);
		m_a[i].trim();

		Diagonal<double> dt(height[i], 0.0);
		for(j=0; j<height[i]; j++) dt.m_val[j] = 2.0 / (2.0 / m_dtime - m_a[i].getValue(j));
		if(i > 0) m_d[2*i-1].setTimes(dt, m_d[2*i-1]);
		if(i + 1 < m_forms) m_d[2*i].setTimes(dt, m_d[2*i]);
		m_a[i].setTimes(dt, m_a[i]);
		m_e[i].setTimes(dt, m_e[i]);
		m_e[i].scale(0.5);
		if(f.size() > i && f[i].m_height == height[i]) m_f[i].setTimes(dt, f[i]);
		else m_f[i].setSparse(height[i]); 
		m_v[i].setSparse(height[i]); // initialize emission term as zero
	}
}

void TimeIntegrator::integratePeriod(Buffer< Column<double> > &v, Buffer<double (*)(const double)> &func) {
	if(v.size() < m_forms) {
		cout << "TimeIntegrator::integratePeriod -> argument v has only " << v.size() << " terms. It should have at least " << m_forms << " terms." << endl;
		return;
	}
	uint i, j, step;
	for(step=0; step<m_steps; step++) {
		for(j=0; j<2; j++) { // iterate even (0) and odd (1)
			for(i=j; i<m_forms; i+=2) {
				v[i] += m_a[i] * v[i]; // absorption
				v[i] += m_v[i]; // emission (1st of 2)
				m_v[i] -= m_e[i] * v[i]; // decrease emission term
				if(i > 0) v[i] += m_d[2 * i - 1] * v[i - 1]; // derivative of v[i-1]
				if(i + 1 < m_forms) v[i] += m_d[2 * i] * v[i + 1]; // derivative of v[i+1]
				if(func.size() > i) v[i] += scaled(m_f[i], func[i](m_time)); // source function
				m_v[i] += m_e[i] * v[i]; // increase emission term
				v[i] += m_v[i]; // emission (2nd of 2)
			}
			m_time += 0.5 * m_dtime;
		}
	}	
}
/*
TimeIntegrator::TimeIntegrator(const Sparse<double> &eDo, const Diagonal<double> &eAe, const Sparse<double> &oDe, const Diagonal<double> &oAo, const Column<double> &eF, const double dtime) {
	uint i;
	m_time = 0.0;
	
	// compute number of time step by CFL-condition
	m_steps = 1;
	Sparse<double> product;
	product.setTimes(eDo, oDe);
	for(i=0; i<product.m_height; i++) {
		const uint steps = uint(dtime * sqrt(-0.5 * product.getValue(i, i)) + 1.0);
		if(m_steps < steps) m_steps = steps;
	}
	product.setTimes(oDe, eDo);
	for(i=0; i<product.m_height; i++) {
		const uint steps = uint(dtime * sqrt(-0.5 * product.getValue(i, i)) + 1.0);
		if(m_steps < steps) m_steps = steps;
	}
	maxMPI(&m_steps, 1);
	m_dtime = dtime / double(m_steps);

	// set integration terms
	const uint eheight = eDo.m_height;
	if(eAe.m_height == eheight) m_eAe.setNegation(eAe);
	else m_eAe.setSparse(eheight);
	m_eEe.setFunction(m_eAe, selectPositive);
	m_eEe.trimSparse();
	m_eAe.setFunction(m_eAe, selectNegative);
	m_eAe.trim();

	Diagonal<double> dt(eheight, 0.0);
	for(i=0; i<eheight; i++) dt.m_val[i] = 2.0 / (2.0 / m_dtime - m_eAe.getValue(i));
	m_eDo.setTimes(dt, eDo);
	m_eAe.setTimes(dt, m_eAe);
	m_eEe.setTimes(dt, m_eEe);
	m_eEe.scale(0.5);
	m_eF.setTimes(dt, eF);
	m_eE.setSparse(eheight); // initialize emission term as zero

	const uint oheight = oDe.m_height;
	if(oAo.m_height == oheight) m_oAo.setNegation(oAo);
	else m_oAo.setSparse(oheight);
	m_oEo.setFunction(m_oAo, selectPositive);
	m_oEo.trimSparse();
	m_oAo.setFunction(m_oAo, selectNegative);
	m_oAo.trim();

	dt.setFullOfZeros(oheight);
	for(i=0; i<oheight; i++) dt.m_val[i] = 2.0 / (2.0 / m_dtime - m_oAo.getValue(i));
	m_oDe.setTimes(dt, oDe);
	m_oAo.setTimes(dt, m_oAo);
	m_oEo.setTimes(dt, m_oEo);
	m_oEo.scale(0.5);
	//m_oF.setTimes(dt, oF);
	m_oE.setSparse(oheight); // initialize emission term as zero
}
void TimeIntegrator::integratePeriod(Column<double> &e, Column<double> &o, double efunc(const double)) {
	Column<double> c;
	for(uint step=0; step<m_steps; step++) {
		o += m_oAo * o; // absorption
		o += m_oE; // emission
		c.setTimes(m_oDe, e); // derivative
		//c += scaled(m_oF, ofunc(m_time)); // source function
		o += c;
		m_oE += m_oEo * c;
		o += m_oE; // emission
		m_time += 0.5 * m_dtime;

		e += m_eAe * e; // absorption
		e += m_eE; // emission
		c.setTimes(m_eDo, o); // derivative
		c += scaled(m_eF, efunc(m_time)); // source function
		e += c;
		m_eE += m_eEe * c;
		e += m_eE; // emission
		m_time += 0.5 * m_dtime;
	}
}
*/