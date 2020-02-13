#include "TimeIntegrator.hpp"
//#include "../Types/MpiEasy.hpp"
#include "../Types/Random.hpp"

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
	Buffer<uint> height(m_forms);
	for(i=0; i<m_forms; i++) {
		if(i == 0) height[i] = d[0].m_height;
		else {
			const uint upd = 2 * i - 2;
			const uint leftd = 2 * i - 1;
			height[i] = d[leftd].m_height;
			if(d[upd].m_height != height[i-1] || d[upd].m_width != height[i]) {
				cout << "TimeIntegrator::TimeIntegrator -> dimensions of derivative d[" << upd << "] do not match." << endl;
				return;
			}
			if(d[leftd].m_width != height[i-1]) {
				cout << "TimeIntegrator::TimeIntegrator -> width of derivative d[" << leftd << "] do not match." << endl;
				return;
			}
		}
	}

	// compute number of time steps by CFL-condition
	m_steps = 1;
	for(i=0; i<m_forms; i++) {
		Sparse<double> product;
		if(i == 0) product.setTimes(d[0], d[1]);
		else {
			product.setTimes(d[2 * i - 1], d[2 * j - 2]);
			if(i + 1 < m_forms) product += d[2 * i] * d[2 * i + 1];
		}
		for(j=0; j<height[i]; j++) { 
			const uint steps = uint(dtime * sqrt(-0.5 * product.getValue(j, j)) + 1.0);
			if(m_steps < steps) m_steps = steps;
		}
	}
	maxMPI(&m_steps, 1);
	m_dtime = dtime / double(m_steps);

	// for testing: randomize categories for each term
	const uint num = 3;
	Random rnd(5345345 * getMPIrank());
	Buffer< Buffer<uint> > c(m_forms);
	for(i=0; i<c.size(); i++) {
		c[i].resize(height[i]);
		for(j=0; j<c[i].size(); j++) c[i][j] = rnd.getUint() % num;
	}
	Buffer<double> ctime(num, m_dtime);
	for(i=1; i<num; i++) ctime[i] = ctime[i-1] / 3.0;

	// set integration terms
	for(i=0; i<m_forms; i++) {
		Diagonal<double> ai(0.0);
		if(a.size() > i && a[i].m_height == height[i]) ai.setNegation(a[i]);
		else ai.setSparse(height[i]);
		Diagonal<double> ei(0.0);
		ei.setFunction(ai, selectPositive);
		ei.trimSparse();
		ei.scale(0.5);
		ai.setFunction(ai, selectNegative);
		ai.trim();

		Diagonal<double> dt(height[i], 0.0);
		for(j=0; j<height[i]; j++) dt.m_val[j] = 2.0 / (2.0 / ctime[c[i][j]] - ai.getValue(j));
		if(i > 0) m_d[2*i-1].init(c[i], dt * d[2*i-1]);
		if(i + 1 < m_forms) m_d[2*i].init(c[i], dt * d[2*i]);
		m_a[i].init(c[i], dt * ai);
		m_e[i].init(c[i], dt * ei);
		if(f.size() > i && f[i].m_height == height[i]) m_f[i].init(c[i], dt * f[i]);
		m_v[i].setSparse(height[i]); // initialize emission term as zero
	}

/*	// set integration terms
	for(i=0; i<m_forms; i++) {
		if(i > 0) m_d[2*i-1].init(c[i], d[2*i-1]);
		if(i + 1 < m_forms) m_d[2*i].init(c[i], d[2*i]);
		//if(i > 0) m_d[2*i-1].init(c[i], d[2*i-1], c[i-1]);
		//if(i + 1 < m_forms) m_d[2*i].init(c[i], d[2*i], c[i+1]);
		//if(i > 0) m_d[2*i-1].init(c[i], d[2*i-1]);
		//if(i + 1 < m_forms) m_d[2*i].init(d[2*i], c[i+1]);
		//if(i > 0) m_d[2*i-1].init(d[2*i-1],  m_d[2*i-2].getCategories());
		//if(i + 1 < m_forms) m_d[2*i].init(c[i], d[2*i], c[i+1]);
	}
	for(i=0; i<m_d.size(); i++) {
		const Buffer< Sparse<double> > &term = m_d[i].m_term;
		double dtime = m_dtime;
		for(j=0; j<term.size(); j++) {
			term[j].scale(dtime);
			dtime /= 3.0;
		}
	}
*/}

void TimeIntegrator::integratePeriod(Buffer< Column<double> > &v, Buffer<double (*)(const double)> &func) {
	if(v.size() < m_forms) {
		cout << "TimeIntegrator::integratePeriod -> argument v has only " << v.size() << " terms. It should have at least " << m_forms << " terms." << endl;
		return;
	}

	const uint exponent = 3;

	uint i, j, k, step;
	uint maxterms = 1;
	for(i=0; i<m_d.size(); i++) {
		const uint terms = m_d[i].m_term.size();
		if(maxterms < terms) maxterms = terms;
	}
	uint ticks = 2;
	for(i=1; i<maxterms; i++) ticks *= exponent;
	const double ticktime = m_dtime / double(ticks);

	for(step=0; step<m_steps; step++) {
		for(j=0; j<ticks; j++) {
			uint cat = 0;
			for(i=ticks/2; (j % i)!=0; i/=exponent) cat++;

			if(m_time+1e-13 < m_dtime) {
				cout << "cat " << cat << " " << j << " " << (j%2) << endl;
			}

/*
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
*/


		// this need repairing

			for(i=j%2; i<m_forms; i+=2) {
				const Buffer< Diagonal<double> > &ai = m_a[i].m_term;
				for(k=cat; k<ai.size(); k++) v[i] += ai[k] * v[i]; // absorption

				v[i] += m_v[i]; // emission (1st of 2)

				//const Buffer< Diagonal<double> > &ei = m_e[i].m_term;
				//for(k=cat; k<ei.size(); k++) m_v[i] -= ei[k] * v[i]; // decrease emission term
				
				if(i > 0) {
					const Buffer< Sparse<double> > &di = m_d[2 * i - 1].m_term;
					for(k=cat; k<di.size(); k++) v[i] += di[k] * v[i - 1]; // derivative of v[i-1]
				}
				if(i + 1 < m_forms) {
					const Buffer< Sparse<double> > &di = m_d[2 * i].m_term;
					for(k=cat; k<di.size(); k++) v[i] += di[k] * v[i + 1]; // derivative of v[i+1]
				}

				if(func.size() > i) {
					const Buffer< Column<double> > &fi = m_f[i].m_term;
					for(k=cat; k<fi.size(); k++) v[i] += scaled(fi[k], func[i](m_time)); // source function
				}

				//for(k=cat; k<ei.size(); k++) m_v[i] += ei[k] * v[i]; // increase emission term
				
				v[i] += m_v[i]; // emission (2nd of 2)
			}
			m_time += ticktime;
		}
	}
}
void TimeIntegrator::integratePeriod(Buffer< Column<double> > &v) {
	if(v.size() < m_forms) {
		cout << "TimeIntegrator::integratePeriod -> argument v has only " << v.size() << " terms. It should have at least " << m_forms << " terms." << endl;
		return;
	}

	const uint exponent = 3;

	uint i, j, k, step;
	uint maxterms = 1;
	for(i=0; i<m_d.size(); i++) {
		const uint terms = m_d[i].m_term.size();
		if(maxterms < terms) maxterms = terms;
	}
	uint ticks = 2;
	for(i=1; i<maxterms; i++) ticks *= exponent;
	const double ticktime = m_dtime / double(ticks);

/*	// odd terms half step backward
	for(i=0; i<m_forms; i+=2) {
		if(i > 0) {
			const Buffer< Sparse<double> > &term = m_d[2 * i - 1].m_term;
			for(k=0; k<term.size(); k++) v[i] += scaled(-0.5, term[k] * v[i - 1]); // derivative of v[i-1]
		}
		if(i + 1 < m_forms) {
			const Buffer< Sparse<double> > &term = m_d[2 * i].m_term;
			for(k=0; k<term.size(); k++) v[i] += scaled(-0.5, term[k] * v[i + 1]); // derivative of v[i+1]
		}
	}
*/

	for(step=0; step<m_steps; step++) {
		for(j=0; j<ticks; j++) {
			uint cat = 0;
			for(i=ticks/2; (j % i)!=0; i/=exponent) cat++;

			if(m_time+1e-13 < m_dtime) {
				cout << "cat " << cat << " " << j << " " << (j%2) << endl;
			}

			for(i=j%2; i<m_forms; i+=2) {
				if(i > 0) {
					const Buffer< Sparse<double> > &term = m_d[2 * i - 1].m_term;
					for(k=cat; k<term.size(); k++) v[i] += term[k] * v[i - 1]; // derivative of v[i-1]
				}
				if(i + 1 < m_forms) {
					const Buffer< Sparse<double> > &term = m_d[2 * i].m_term;
					for(k=cat; k<term.size(); k++) v[i] += term[k] * v[i + 1]; // derivative of v[i+1]
				}
			}
			m_time += ticktime;
		}
	}

/*	// odd terms half step forward
	for(i=0; i<m_forms; i+=2) {
		if(i > 0) {
			const Buffer< Sparse<double> > &term = m_d[2 * i - 1].m_term;
			for(k=0; k<term.size(); k++) v[i] += scaled(0.5, term[k] * v[i - 1]); // derivative of v[i-1]
		}
		if(i + 1 < m_forms) {
			const Buffer< Sparse<double> > &term = m_d[2 * i].m_term;
			for(k=0; k<term.size(); k++) v[i] += scaled(0.5, term[k] * v[i + 1]); // derivative of v[i+1]
		}
	}
*/
/*	for(step=0; step<m_steps; step++) {
		for(j=0; j<2; j++) { // iterate even (0) and odd (1)
			for(i=j; i<m_forms; i+=2) {
				if(i > 0) v[i] += m_d[2 * i - 1].m_term[0] * v[i - 1]; // derivative of v[i-1]
				if(i + 1 < m_forms) v[i] += m_d[2 * i].m_term[0] * v[i + 1]; // derivative of v[i+1]
			}
			m_time += 0.5 * m_dtime;
		}
	}	
*/}

/*TimeIntegrator::TimeIntegrator(const Buffer< Sparse<double> > &d, const Buffer< Diagonal<double> > &a, const Buffer< Column<double> > &f, const double dtime) {
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
*/