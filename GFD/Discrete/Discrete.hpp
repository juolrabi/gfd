/**
 * Discrete.hpp implements a class, where to inherit discrete form, discrete Hodge, discrete exterior derivative etc.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _DISCRETE_HPP_INCLUDED_
#define _DISCRETE_HPP_INCLUDED_

#include "../Types/Types.hpp"
#include "../Types/MpiEasy.hpp"
#include <iostream>

using namespace std;

namespace gfd {

template<typename T,typename R> T functionNegation(const R &r) { return -r; }
template<typename T,typename R> T functionInverse(const R &r) { return 1.0 / r; }
template<typename T,typename L,typename R> T functionPlus(const L &l, const R &r) { return l + r; }
template<typename T,typename L,typename R> T functionMinus(const L &l, const R &r) { return l - r; }
template<typename T,typename L,typename R> T functionTimes(const L &l, const R &r) { return l * r; }
template<typename T,typename L,typename R> T functionDividedBy(const L &l, const R &r) { return l / r; }
template<typename T,typename L,typename R> T functionReversePlus(const R &r, const L &l) { return l + r; }
template<typename T,typename L,typename R> T functionReverseMinus(const R &r, const L &l) { return l - r; }
template<typename L,typename R> double functionDot(const L &l, const R &r) { return l.dot(r); }

template<typename T>
class Discrete
{
public:
	// constructors
	Discrete(const T &zero) : m_zero(zero) {
		m_height = 0;
		m_full = true;
	}
	virtual ~Discrete() { }

	// all variables are public
	T m_zero; // template zero value
	uint m_height; // matrix height
	bool m_full; // condition of full rows
	Buffer<uint> m_row; // row indices in the case of sparse rows
	Buffer<T> m_val; // term values

protected:
	void printData() const {
		uint i;
		cout << "m_height = " << m_height << endl;
		cout << "m_full = " << m_full << endl;
		cout << "m_row =";
		for(i=0; i<m_row.size(); i++) cout << " " << m_row[i];
		cout << endl;
		cout << "m_val =";
		for(i=0; i<m_val.size(); i++) cout << " " << m_val[i];
		cout << endl;
	}
	void getFromRank(const uint from, const uint to, Discrete &recv) const {
		if(getMPIrank() == from) {
			if(getMPIrank() == to) {
				recv.setCopy(*this);
				return;
			}
			Discrete send(*this);
			sendMPI(&send.m_zero, sizeof(T), to, 0);
			sendMPI(&send.m_height, sizeof(uint), to, 0);
			sendMPI(&send.m_full, sizeof(bool), to, 0);
			uint isize = send.m_row.size();
			sendMPI(&isize, sizeof(uint), to, 0);
			if(isize > 0) sendMPI(&send.m_row[0], isize * sizeof(uint), to, 0);
			isize = send.m_val.size();
			sendMPI(&isize, sizeof(uint), to, 0);
			if(isize > 0) sendMPI(&send.m_val[0], isize * sizeof(T), to, 0);
		}
		else if(getMPIrank() == to) {
			recvMPI(&recv.m_zero, sizeof(T), from, 0);
			recvMPI(&recv.m_height, sizeof(uint), from, 0);
			recvMPI(&recv.m_full, sizeof(bool), from, 0);
			uint isize;
			recvMPI(&isize, sizeof(uint), from, 0);
			recv.m_row.resize(isize);
			if(isize > 0) recvMPI(&recv.m_row[0], isize * sizeof(uint), from, 0);
			recvMPI(&isize, sizeof(uint), from, 0);
			recv.m_val.resize(isize);
			if(isize > 0) recvMPI(&recv.m_val[0], isize * sizeof(T), from, 0);
		}
	}
	Discrete &setEmpty() {
		m_height = 0;
		m_full = true;
		m_row.clear();
		m_val.clear();
		return *this;
	}
	Discrete &setFullOfZeros(const uint height) {
		m_height = height;
		m_full = true;
		m_row.clear();
		m_val.resize(m_height);
		m_val.fill(m_zero);
		return *this;
	}
	Discrete &setFull(const Buffer<T> &val) {
		m_height = val.size();
		m_full = true;
		m_row.clear();
		m_val = val;
		return *this;
	}
	Discrete &setSparse(const uint height, const Buffer< pair<uint, T> > &val) {
		m_height = height;
		m_full = false;
		const Buffer<uint> ord = bubbleSort(val);
		m_row.resize(ord.size());
		m_val.resize(ord.size());
		for(uint i=0; i<ord.size(); i++) {
			m_row[i] = val[ord[i]].first;
			m_val[i] = val[ord[i]].second;
		}
		return *this;
	}
	template<typename R> Discrete &setShape(const Discrete<R> &r) {
		m_height = r.m_height;
		m_full = r.m_full;
		m_row = r.m_row;
		m_val.resize(r.m_val.size());
		return *this;
	}
	template<typename R> Discrete &setCopy(const Discrete<R> &r) {
		setShape(r);
		for(uint i=0; i<m_val.size(); i++) m_val[i] = r.m_val[i];
		return *this;
	}
	template<typename L, typename R> Discrete &setUnion(const Discrete<L> &l, const Discrete<R> &r, T func(const L &, const R &)) {
		if(orMPI(l.m_height != r.m_height)) return setEmpty(); // the heights do not match
		uint i = 0;
		if(l.m_full) {
			if(r.m_full) { // both l and r are full column vectors
				if(this != (void*)&l && this != (void*)&r) setShape(l);
				for(i=0; i<l.m_height; i++) m_val[i] = func(l.m_val[i], r.m_val[i]);
				return *this;
			}
			// only l is a full column vector
			uint ri = 0;
			if(this == (void*)&r) {
				Buffer<T> val(l.m_height);
				while(ri < r.m_row.size()) {
					while(i < r.m_row[ri]) { val[i] = func(l.m_val[i], r.m_zero); i++; }
					val[i] = func(l.m_val[i], r.m_val[ri]); i++; ri++;
				}
				while(i < l.m_height) { val[i] = func(l.m_val[i], r.m_zero); i++; }
				m_val.swap(val);
				setShape(l);
			}
			else {
				if(this != (void*)&l) setShape(l);
				while(ri < r.m_row.size()) {
					while(i < r.m_row[ri]) { m_val[i] = func(l.m_val[i], r.m_zero); i++; }
					m_val[i] = func(l.m_val[i], r.m_val[ri]); i++; ri++;
				}
				while(i < l.m_height) { m_val[i] = func(l.m_val[i], r.m_zero); i++; }
			}
			return *this;
		}
		if(r.m_full) { // only r is a full column vector
			uint li = 0;
			if(this == (void*)&l) {
				Buffer<T> val(r.m_height);
				while(li < l.m_row.size()) {
					while(i < l.m_row[li]) { val[i] = func(l.m_zero, r.m_val[i]); i++; }
					val[i] = func(l.m_val[li], r.m_val[i]); i++; li++;
				}
				while(i < r.m_height) { val[i] = func(l.m_zero, r.m_val[i]); i++; }
				m_val.swap(val);
				setShape(r);
			}
			else {
				if(this != (void*)&r) setShape(r);
				while(li < l.m_row.size()) {
					while(i < l.m_row[li]) { m_val[i] = func(l.m_zero, r.m_val[i]); i++; }
					m_val[i] = func(l.m_val[li], r.m_val[i]); i++; li++;
				}
				while(i < r.m_height) { m_val[i] = func(l.m_zero, r.m_val[i]); i++; }
			}
			return *this;
		}
		// both l and r are sparse column vectors
		uint rowsize = l.m_row.size() + r.m_row.size();
		if(rowsize > l.m_height) rowsize = l.m_height;
		Buffer<uint> row(rowsize);
		Buffer<T> val(rowsize);
		uint li = 0;
		uint ri = 0;
		while(li < l.m_row.size() && ri < r.m_row.size()) {
			if(l.m_row[li] < r.m_row[ri]) {
				row[i] = l.m_row[li];
				val[i] = func(l.m_val[li], r.m_zero);
				i++; li++;
			}
			else if(r.m_row[ri] < l.m_row[li]) {
				row[i] = r.m_row[ri];
				val[i] = func(l.m_zero, r.m_val[ri]);
				i++; ri++;
			}
			else {
				row[i] = l.m_row[li];
				val[i] = func(l.m_val[li], r.m_val[ri]);
				i++; li++; ri++;
			}
		}
		while(li < l.m_row.size()) {
			row[i] = l.m_row[li];
			val[i] = func(l.m_val[li], r.m_zero);
			i++; li++;
		}
		while(ri < r.m_row.size()) {
			row[i] = r.m_row[ri];
			val[i] = func(l.m_zero, r.m_val[ri]);
			i++; ri++;
		}
		m_full = false;
		m_height = l.m_height;
		m_row.copy(row, i);
		m_val.copy(val, i);
		return *this;
	}
	template<typename L, typename R> Discrete &setPlus(const Discrete<L> &l, const Discrete<R> &r) {
		return setUnion(l, r, functionPlus);
	}
	template<typename L, typename R> Discrete &setMinus(const Discrete<L> &l, const Discrete<R> &r) {
		return setUnion(l, r, functionMinus);
	}
	template<typename L, typename R> Discrete &setIntersection(const Discrete<L> &l, const Discrete<R> &r, T func(const L &, const R &)) {
		if(orMPI(l.m_height != r.m_height)) return setEmpty(); // the heights do not match
		uint i, j, k;
		if(l.m_full) {
			if(this != (void*)&l && this != (void*)&r) setShape(r);
			if(r.m_full) { // both l and r are full column vectors
				for(i=0; i<m_height; i++) m_val[i] = func(l.m_val[i], r.m_val[i]);
				return *this;
			}
			// only l is a full column vector
			for(i=0; i<r.m_row.size(); i++) m_val[i] = func(l.m_val[r.m_row[i]], r.m_val[i]);
			if(this == (void*)&l) setShape(r);
			return *this;
		}
		if(this != (void*)&l && this != (void*)&r) setShape(l);
		if(r.m_full) { // only r is a full column vector
			for(i=0; i<l.m_row.size(); i++) m_val[i] = func(l.m_val[i], r.m_val[l.m_row[i]]);
			if(this == (void*)&r) setShape(l);
			return *this;
		}
		// both l and r are sparse column vectors
		for(i=0,j=0,k=0; i<l.m_row.size(); i++) {
			while(j < r.m_row.size() && r.m_row[j] < l.m_row[i]) j++;
			if(j < r.m_row.size() && r.m_row[j] == l.m_row[i]) {
				m_row[k] = l.m_row[i];
				m_val[k++] = func(l.m_val[i], r.m_val[j]);
			}
		}
		m_row.resize(k);
		m_val.resize(k);
		return *this;
	}
	template<typename L, typename R> Discrete &setTimes(const Discrete<L> &l, const Discrete<R> &r) {
		return setIntersection(l, r, functionTimes);
	}
	template<typename L, typename R> Discrete &setScale(const Discrete<L> &l, const R &r) {
		if(this != (void*)&l) setShape(l);
		for(uint i=0; i<m_val.size(); i++) m_val[i] = l.m_val[i] * r;
		return *this;
	}
	template<typename L, typename R> Discrete &setScale(const L &l, const Discrete<R> &r) {
		if(this != (void*)&r) setShape(r);
		for(uint i=0; i<m_val.size(); i++) m_val[i] = l * r.m_val[i];
		return *this;
	}
	template<typename R> Discrete &setFunction(const Discrete<R> &r, T func(const R &)) {
		if(this != (void*)&r) setShape(r);
		for(uint i=0; i<m_val.size(); i++) m_val[i] = func(r.m_val[i]);
		return *this;
	}
	template<typename R> Discrete &setNegation(const Discrete<R> &r) { return setFunction(r, functionNegation); }
	template<typename R> Discrete &setInverse(const Discrete<R> &r) { return setFunction(r, functionInverse); }

	Discrete &trimFull() { // convert sparse vector to full vector
		if(m_full) return *this; // already full
		m_full = true;
		Buffer<T> val(m_height, m_zero);
		for(uint i=0; i<m_row.size(); i++) val[m_row[i]] = m_val[i];
		m_val.swap(val);
		m_row.clear();
		return *this;
	}
	Discrete &trimSparse() { // convert full vector to sparse vector
		if(!m_full) return trim(); // already sparse
		m_full = false;
		m_row.resize(m_height);
		for(uint i=0; i<m_height; i++) m_row[i] = i;
		return trim();
	}
	Discrete &trimOptimal(const double limit) { // convert to full if(number of non-empty rows > limit * m_height), otherwise convert to sparse
		trimSparse();
		if(double(m_row.size()) > limit * m_height) trimFull();
		return *this;
	}
	Discrete &trim() { // remove all zero instances from a sparse column vector
		if(m_full) return *this; // not sparse
		uint i, j;
		for(i=0,j=0; i<m_row.size(); i++) {
			if(m_val[i] == m_zero) continue;
			if(j != i) {
				m_row[j] = m_row[i];
				m_val[j] = m_val[i];
			}
			j++;
		}
		m_row.resize(j);
		m_val.resize(j);
		return *this;
	}

	const T &getValue(uint i) const {
		if(m_full) {
			if(i >= m_val.size()) return m_zero;
		}
		else if(!searchIndex(i, m_row, 0, m_row.size())) return m_zero;
		return m_val[i];
	}
	template<typename R> double getDot(const Discrete<R> &r) const { return getProduct(r, functionTimes<double,T,R>); }
	template<typename R> double getDotDot(const Discrete<R> &r) const { return getProduct(r, functionDot); }
	double getLensq() const { return getProduct(*this, functionTimes<double,T,T>); }
	double getLensqDot() const { return getProduct(*this, functionDot); }
	template<typename R> double getProduct(const Discrete<R> &r, double func(const T &, const R &)) const {
		if(orMPI(m_height != r.m_height)) return 0.0; // the heights do not match
		uint i, j;
		double sum = 0.0;
		if(this == (void*)&r) { // self operation
			for(i=0; i<m_val.size(); i++) sum += func(m_val[i], m_val[i]);
			sumMPI(&sum, 1);
			return sum;
		}
		else if(m_full) {
			if(r.m_full) { // both this and r are full column vectors
				for(i=0; i<m_height; i++) sum += func(m_val[i], r.m_val[i]);
			}
			else { // only this is a full column vector
				for(i=0; i<r.m_row.size(); i++) sum += func(m_val[r.m_row[i]], r.m_val[i]);
			}
		}
		else {
			if(r.m_full) { // only r is a full column vector
				for(i=0; i<m_row.size(); i++) sum += func(m_val[i], r.m_val[m_row[i]]);
			}
			else { // both this and r are sparse column vectors
				for(i=0,j=0; i<m_row.size(); i++) {
					while(j < r.m_row.size() && r.m_row[j] < m_row[i]) j++;
					if(j >= r.m_row.size()) break;
					if(r.m_row[j] == m_row[i]) sum += func(m_val[i], r.m_val[j]);
				}
			}
		}
		sumMPI(&sum, 1);
		return sum;
	}
	Buffer<T> getBuffer() const {
		if(m_full) return m_val;
		Buffer<T> val(m_height, m_zero);
		for(uint i=0; i<m_row.size(); i++) val[m_row[i]] = m_val[i];
		return val;
	}
	template<typename R> Buffer<uint> bubbleSort(const Buffer< pair<uint, R> > &val) const { // val[i].first should be unique
		uint i, j;
		Buffer<uint> ord(val.size());
		for(i=0; i<val.size(); i++) {
			const uint first = val[i].first;
			for(j=i; j>0 && val[ord[j-1]].first>first; j--) ord[j] = ord[j-1];
			ord[j] = i;
		}
		return ord;
	}
	bool searchIndex(uint &i, const Buffer<uint> &inc, uint i0, uint i1) const { // find index with binary search algorithm
		const uint target = i;
		while(i0 < i1) {
			i = (i0 + i1) / 2;
			if(target < inc[i]) i1 = i;
			else if(inc[i] < target) i0 = i + 1;
			else return true;
		}
		return false;
	}

};

}

#endif //_DISCRETE_HPP_INCLUDED_
