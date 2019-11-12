/**
 * Sequence.hpp implements class Sequence, which contains a sequence of Operations.
 * Operation is a general representation for any matrix multiplication kind of operation.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _SEQUENCE_HPP_INCLUDED_
#define _SEQUENCE_HPP_INCLUDED_

#include "Column.hpp"

namespace gfd
{

class Operation
{
public:
	virtual ~Operation() { }
	virtual Operation *clone() const = 0;
	virtual void operate() const = 0;
};

template<typename T, typename L, typename R>
class SparseOperation : public Operation
{
public:
	SparseOperation(Column<T> *t, const Sparse<L> *l, const Column<R> *r) {
		m_t = t;
		m_l = l;
		m_r = r;
	}
	virtual ~SparseOperation() { }
	virtual Operation *clone() const {
		return new SparseOperation<T,L,R>(m_t, m_l, m_r);
	}
	virtual void operate() const {
		(*m_t) += (*m_l) * (*m_r);
//		m_o->setTimes(*m_l, *m_r);
	}

protected:
	Column<T> *m_t;
	const Sparse<L> *m_l;
	const Column<R> *m_r;
};
template<typename T, typename L, typename R>
SparseOperation<T,L,R> makeSparseOperation(Column<T> *t, const Sparse<L> *l, const Column<R> *r) {
	return SparseOperation<T,L,R>(t,l,r);
}

template<typename T, typename L, typename R>
class DiagonalOperation : public Operation
{
public:
	DiagonalOperation(Column<T> *t, const Diagonal<L> *l, const Column<R> *r) {
		m_t = t;
		m_l = l;
		m_r = r;
	}
	virtual ~DiagonalOperation() { }
	virtual Operation *clone() const {
		return new DiagonalOperation<T,L,R>(m_t, m_l, m_r);
	}
	virtual void operate() const {
		(*m_t) += (*m_l) * (*m_r);
	}

protected:
	Column<T> *m_t;
	const Diagonal<L> *m_l;
	const Column<R> *m_r;
};
template<typename T, typename L, typename R>
DiagonalOperation<T,L,R> makeDiagonalOperation(Column<T> *t, const Diagonal<L> *l, const Column<R> *r) {
	return DiagonalOperation<T,L,R>(t,l,r);
}


class Sequence
{
public:
	Sequence() {
		m_opers = 0;
	}
	virtual ~Sequence() { clear(); }
	void clear() {
		for(uint i=0; i<m_opers; i++) delete m_oper[i];
		m_oper.clear();
		m_opers = 0;
	}
	void gather(const Operation &oper) {
		m_oper.gather(oper.clone(), m_opers);
	}
	void operate() const {
		for(uint i=0; i<m_opers; i++) m_oper[i]->operate();
	}

protected:
	uint m_opers;
	Buffer<Operation *> m_oper;
};



};

#endif //_SEQUENCE_HPP_INCLUDED_
