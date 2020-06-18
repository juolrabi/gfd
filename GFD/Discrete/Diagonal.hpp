/**
 * Diagonal.hpp implements a diagonal matrix. The diagonal matrix can be sparse or full of values.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
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
	template<typename R> Diagonal(const Diagonal<R> &r) : Discrete<T>::Discrete(r.m_zero) { setCopy(r); } // copy existing diagonal matrix
	Diagonal(const uint height, const T &zero) : Discrete<T>::Discrete(zero) { setFullOfZeros(height); } // initialize full diagonal matrix with given height
	Diagonal(const Buffer<T> &val, const T &zero) : Discrete<T>::Discrete(zero) { setFull(val); } // initialize full diagonal matrix
	Diagonal(const uint height, const Buffer< pair<uint, T> > &val, const T &zero) : Discrete<T>::Discrete(zero) { setSparse(height, val); } // initialize sparse diagonal matrix
	virtual ~Diagonal() { }

	// print functions
	void printData() const {
		cout << "Diagonal matrix:" << endl;
		Discrete<T>::printData();
	}
	void printShape() const {
		uint i, j, k;
		Diagonal data(this->m_zero);
		const uint ranks = getMPIranks();
		uint col0 = 0;
		for(k=0; k<ranks; k++) {
			Discrete<T>::getFromRank(k, 0, data);
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
	Diagonal &setFullOfZeros(const uint height) { Discrete<T>::setFullOfZeros(height); return *this; }
	Diagonal &setFull(const Buffer<T> &val) { Discrete<T>::setFull(val); return *this; }
	Diagonal &setSparse(const uint height, const Buffer< pair<uint, T> > &val = Buffer< pair<uint, T> >()) { Discrete<T>::setSparse(height, val); return *this; }
	template<typename R> Diagonal &setCopy(const Diagonal<R> &r) { Discrete<T>::setCopy(r); return *this; }
	template<typename R> Diagonal &setFunction(const Diagonal<R> &r, T func(const R &)) { Discrete<T>::setFunction(r, func); return *this; }
	template<typename R> Diagonal &setNegation(const Diagonal<R> &r) { Discrete<T>::setNegation(r); return *this; }
	template<typename R> Diagonal &setInverse(const Diagonal<R> &r) { Discrete<T>::setInverse(r); return *this; }
	template<typename L, typename R> Diagonal &setUnion(const Diagonal<L> &l, const Diagonal<R> &r, T func(const L &, const R &)) { Discrete<T>::setUnion(l, r, func); return *this; }
	template<typename L, typename R> Diagonal &setIntersection(const Diagonal<L> &l, const Diagonal<R> &r, T func(const L &, const R &)) { Discrete<T>::setIntersection(l, r, func); return *this; }
	template<typename L, typename R> Diagonal &setPlus(const Diagonal<L> &l, const Diagonal<R> &r) { Discrete<T>::setPlus(l, r); return *this; }
	template<typename L, typename R> Diagonal &setMinus(const Diagonal<L> &l, const Diagonal<R> &r) { Discrete<T>::setMinus(l, r); return *this; }
	template<typename L, typename R> Diagonal &setTimes(const Diagonal<L> &l, const Diagonal<R> &r) { Discrete<T>::setTimes(l, r); return *this; }
	template<typename L, typename R> Diagonal &setScale(const Diagonal<L> &l, const R &r) { Discrete<T>::setScale(l, r); return *this; }
	template<typename L, typename R> Diagonal &setScale(const L &l, const Diagonal<R> &r) { Discrete<T>::setScale(l, r); return *this; }

	// modifier functions
	template<typename R> Diagonal &operator+=(const Diagonal<R> &r) { return setPlus(*this, r); }
	template<typename R> Diagonal &operator-=(const Diagonal<R> &r) { return setMinus(*this, r); }
	template<typename R> Diagonal &operator*=(const Diagonal<R> &r) { return setTimes(*this, r); }
	template<typename R> Diagonal &scale(const R &r) { return setScale(*this, r); }
	Diagonal &negate() { return setNegation(*this); }
	Diagonal &invert() { return setInverse(*this); }

	// trim functions
	Diagonal &trimFull() { Discrete<T>::trimFull(); return *this; } // convert sparse to full
	Diagonal &trimSparse() { Discrete<T>::trimSparse(); return *this; } // convert full to sparse
	Diagonal &trimOptimal(const double limit = 0.5) { Discrete<T>::trimOptimal(limit); return *this; } // convert to full if(number of non-empty rows > limit * m_height), otherwise convert to sparse
	Diagonal &trim() { Discrete<T>::trim(); return *this; } // remove all zero instances

	// get functions
	Buffer<T> getBuffer() const { return Discrete<T>::getBuffer(); } // return diagonal values in the format of Buffer<T>
	const T &getValue(const uint i) const { return Discrete<T>::getValue(i); }
	template<typename R> double getDot(const Diagonal<R> &r) const { return Discrete<T>::getDot(r); }
	template<typename R> double getDotDot(const Diagonal<R> &r) const { return Discrete<T>::getDotDot(r); }
	double getLensq() const { return Discrete<T>::getLensq(); }
	double getLensqDot() const { return Discrete<T>::getLensqDot(); }
	template<typename R> double getProduct(const Diagonal<R> &r, double func(const T &, const R &)) const { return Discrete<T>::getProduct(r, func); }

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
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Diagonal<O> scaled(const Diagonal<L> &l, const R &r) {
	Diagonal<O> o(l.m_zero * r);
	return o.setScale(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Diagonal<O> scaled(const L &l, const Diagonal<R> &r) {
	Diagonal<O> o(l * r.m_zero);
	return o.setScale(l, r);
}
template<typename R> Diagonal<R> inverse(const Diagonal<R> &r) {
	Diagonal<R> o(r.m_zero);
	return o.setInverse(r);
}
template<typename R> Diagonal<R> operator-(const Diagonal<R> &r) {
	Diagonal<R> o(r.m_zero);
	return o.setNegation(r);
}

}

#endif //_DIAGONAL_HPP_INCLUDED_
