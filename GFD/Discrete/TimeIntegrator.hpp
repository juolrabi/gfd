/**
 * TimeIntegrator implements routines for solving time-dependent wave problems.
 * Integrating following pair of equation in time:
 *   \partial_t e + eAe e = eDo o + eF,
 *   \partial_t o + oAo h = oDe e.
 * Author: Jukka Räbinä, University of Jyväskylä, 2020.
 */

#ifndef _TIMEINTEGRATOR_HPP_INCLUDED_
#define _TIMEINTEGRATOR_HPP_INCLUDED_

#include "Form.hpp"

namespace gfd
{

class TimeIntegrator
{
public:
	TimeIntegrator(const Buffer< Sparse<double> > &d, const Buffer< Diagonal<double> > &a, const Buffer< Column<double> > &f, const double dtime);
//	TimeIntegrator(const Sparse<double> &eDo, const Diagonal<double> &eAe, const Sparse<double> &oDe, const Diagonal<double> &oAo, const Column<double> &eF, const double dtime);
	virtual ~TimeIntegrator() { }

	void integratePeriod(Buffer< Column<double> > &v, Buffer<double (*)(const double)> &func);
//	void integratePeriod(Column<double> &e, Column<double> &o, double efunc(const double));

protected:
	double m_time; // current integration time (starts from zero)
	double m_dtime; // time step duration
	double m_steps; // time steps per period
	uint m_forms; // number of forms
	Buffer< Sparse<double> > m_d; // derivatives
	Buffer< Diagonal<double> > m_a; // absorption terms
	Buffer< Diagonal<double> > m_e; // emission terms
	Buffer< Column<double> > m_v; // emission vector
	Buffer< Column<double> > m_f; // source terms

/*	Sparse<double> m_eDo;
	Diagonal<double> m_eAe;
	Column<double> m_eF;
	Sparse<double> m_oDe;
	Diagonal<double> m_oAo;

	// treatment of emittance (ie. negative absorption)
	Column<double> m_eE;
	Diagonal<double> m_eEe;
	Column<double> m_oE;
	Diagonal<double> m_oEo;
*/

};

}

#endif //_TIMEINTEGRATOR_HPP_INCLUDED_
