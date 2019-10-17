/*
Diagonal.hpp implements a diagonal matrix. The diagonal matrix can be sparse or full of values.
*/

#ifndef _DIAGONAL_HPP_INCLUDED_
#define _DIAGONAL_HPP_INCLUDED_

#include "Discrete.hpp"

using namespace std;

namespace gfd
{

template<typename T>
class Diagonal : public Discrete<T>
{
public:
	// constructors
	Diagonal(const T &zero = 0) : Discrete<T>::Discrete(zero) { } // initialize empty diagonal matrix
	template<typename R> Diagonal(const Diagonal<R> &r) : Discrete<T>::m_zero(r.m_zero) { setCopy(r); } // copy existing diagonal matrix
	Diagonal(const uint height, const T &zero) : Discrete<T>::m_zero(zero) { setFullOfZeros(height); } // initialize full diagonal matrix with given height
	Diagonal(const Buffer<T> &val, const T &zero) : Discrete<T>::m_zero(zero) { setFull(val); } // initialize full diagonal matrix
	Diagonal(const uint height, const Buffer< pair<uint, T> > &val, const T &zero) : Discrete<T>::m_zero(zero) { setSparse(height, val); } // initialize sparse diagonal matrix
	virtual ~Diagonal() { }

	// print functions
	void printData() const {
		cout << "Diagonal matrix:" << endl;
		Discrete<T>::printVectorData();
	}
	void printShape() const {
		uint i, j, k;
		Diagonal data(this->m_zero);
		const uint ranks = getMPIranks();
		uint col0 = 0;
		for(k=0; k<ranks; k++) {
			Discrete<T>::getVectorFromRank(k, 0, data);
			for(i=0; i<data.m_height; i++) {
				for(j=0; j<col0; j++) cout << "  ";
				uint row = i;
				if(!data.m_full && !data.searchIndex(row, data.m_row, 0, data.m_row.size())) {
					for(j=0; j<data.m_height; j++) cout << ". ";
				}
				else {
					for(j=0; j<i; j++) cout << ". ";
					cout << data.m_val[row];
					for(j=i+1; j<data.m_height; j++) cout << " .";
				}
				cout << endl;
			}
			col0 += data.m_height;
		}
	}

	// set functions
	Diagonal &setFullOfZeros(const uint height) { Discrete<T>::setVectorFullOfZeros(height); return *this; }
	Diagonal &setFull(const Buffer<T> &val) { Discrete<T>::setVectorFull(val); return *this; }
	Diagonal &setSparse(const uint height, const Buffer< pair<uint, T> > &val) { Discrete<T>::setVectorSparse(height, val); return *this; }
	template<typename R> Diagonal &setCopy(const Diagonal<R> &r) { Discrete<T>::setVectorCopy(r); return *this; }
	template<typename R> Diagonal &setNegation(const Diagonal<R> &r) { Discrete<T>::setVectorNegation(r); return *this; }
	template<typename L, typename R> Diagonal &setPlus(const Diagonal<L> &l, const Diagonal<R> &r) { Discrete<T>::setVectorPlus(l, r); return *this; }
	template<typename L, typename R> Diagonal &setMinus(const Diagonal<L> &l, const Diagonal<R> &r) { Discrete<T>::setVectorMinus(l, r); return *this; }
	template<typename L, typename R> Diagonal &setTimes(const Diagonal<L> &l, const Diagonal<R> &r) { Discrete<T>::setVectorTimes(l, r); return *this; }
	template<typename L, typename R> Diagonal &setScaleRight(const Diagonal<L> &l, const R &r) { Discrete<T>::setVectorScaleRight(l, r); return *this; }
	template<typename L, typename R> Diagonal &setScaleLeft(const L &l, const Diagonal<R> &r) { Discrete<T>::setVectorScaleLeft(l, r); return *this; }

	// modifier functions
	template<typename R> Diagonal &operator+=(const Diagonal<R> &r) { return setPlus(*this, r); }
	template<typename R> Diagonal &operator-=(const Diagonal<R> &r) { return setMinus(*this, r); }
	template<typename R> Diagonal &operator*=(const Diagonal<R> &r) { return setTimes(*this, r); }
	template<typename R> Diagonal &scaleRight(const R &r) { return setScaleRight(*this, r); }
	template<typename L> Diagonal &scaleLeft(const L &l) { return setScaleLeft(l, *this); }

	// trim functions
	Diagonal &trimFull() { Discrete<T>::trimVectorFull(); return *this; } // convert sparse to full
	Diagonal &trimSparse() { Discrete<T>::trimVectorSparse(); return *this; } // convert full to sparse
	Diagonal &trim() { Discrete<T>::trimVector(); return *this; } // remove all zero instances

	// get functions
	Buffer<T> getBuffer() const { return Discrete<T>::getVectorBuffer(); } // return diagonal values in the format of Buffer<T>
	const T &getValue(const uint i, const uint j) const {
		if(i == j) return Discrete<T>::getVectorValue(i);
		return this->m_zero;
	}

};

// operators
template<typename L, typename R, typename O = decltype(declval<L &>() + declval<R &>())> Diagonal<O> operator+(const Diagonal<L> &l, const Diagonal<R> &r) {
	Diagonal<O> o(l.m_zero + r.m_zero);
	return o.setPlus(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() - declval<R &>())> Diagonal<O> operator-(const Diagonal<L> &l, const Diagonal<R> &r) {
	Diagonal<O> o(l.m_zero - r.m_zero);
	return o.setMinus(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Diagonal<O> operator*(const Diagonal<L> &l, const Diagonal<R> &r) {
	Diagonal<O> o(l.m_zero * r.m_zero);
	return o.setTimes(l, r);
}
template<typename R> Diagonal<R> operator-(const Diagonal<R> &r) {
	Diagonal<R> o(r.m_zero);
	return o.setNegation(r);
}

};

#endif //_DIAGONAL_HPP_INCLUDED_
