/**
 * TimeIntegrator implements routines for solving time-dependent wave problems.
 * Integrating following pair of equation in time:
 *   \partial_t e + eAe e = eDo o + eF,
 *   \partial_t o + oAo h = oDe e.
 * Author: Jukka Räbinä, University of Jyväskylä, 2020.
 */

#ifndef _TIMEINTEGRATOR_HPP_INCLUDED_
#define _TIMEINTEGRATOR_HPP_INCLUDED_

#include "Split.hpp"

namespace gfd
{

class TimeIntegrator
{
public:
//	TimeIntegrator(const Buffer< Sparse<double> > &d, const double dtime);
	TimeIntegrator(const Buffer< Sparse<double> > &d, const Buffer< Diagonal<double> > &a, const Buffer< Column<double> > &f, const double dtime);
//	TimeIntegrator(const Buffer< Sparse<double> > &d, const Buffer< Diagonal<double> > &a, const Buffer< Column<double> > &f, const double dtime);
	virtual ~TimeIntegrator() { }

	void integratePeriod(Buffer< Column<double> > &v);
	void integratePeriod(Buffer< Column<double> > &v, Buffer<double (*)(const double)> &func);

protected:
	double m_time; // current integration time (starts from zero)
	double m_dtime; // time step duration
	double m_steps; // time steps per period
	uint m_forms; // number of forms
	Buffer< SparseSplit<double> > m_d; // derivatives
	Buffer< DiagonalSplit<double> > m_a; // absorption terms
	Buffer< DiagonalSplit<double> > m_e; // emission terms
	Buffer< ColumnSplit<double> > m_v; // emission vector
	Buffer< ColumnSplit<double> > m_f; // source terms

/*	Buffer< Sparse<double> > m_d; // derivatives
	Buffer< Diagonal<double> > m_a; // absorption terms
	Buffer< Diagonal<double> > m_e; // emission terms
	Buffer< Column<double> > m_v; // emission vector
	Buffer< Column<double> > m_f; // source terms
*/
};

}

#endif //_TIMEINTEGRATOR_HPP_INCLUDED_
