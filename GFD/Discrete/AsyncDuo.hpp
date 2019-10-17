/*
AsyncDuo.hpp
*/

#ifndef _ASYNCDUO_HPP_INCLUDED_
#define _ASYNCDUO_HPP_INCLUDED_

#include "Dec.hpp"

namespace gfd
{

template<typename A, typename B>
class AsyncDuo
{
public:
	AsyncDuo() {

	}
	AsyncDuo &setDuo(const Sparse<A> &a, const Column<double> &even, const Sparse<B> &b, const Column<double> &odd) {

		return *this;
	}
	template<typename E, typename O> bool operate(Column<E> &even, Column<O> &odd) {
		odd += m_a * even;
		even += m_b * odd;
		return *this;
	}

protected:
	Sparse<A> m_a;
	Sparse<B> m_b;
};



};

#endif //_ASYNCDUO_HPP_INCLUDED_
