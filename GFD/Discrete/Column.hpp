/**
 * Column.hpp implements a column vector. The column vector can be sparse or full of values.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _COLUMN_HPP_INCLUDED_
#define _COLUMN_HPP_INCLUDED_

#include "Sparse.hpp"

using namespace std;

namespace gfd
{

template<typename T>
class Column : public Diagonal<T>
{
public:
	// constructors
	Column(const T &zero = 0) : Diagonal<T>::Diagonal(zero) { } // initialize empty column vector
	template<typename R> Column(const Diagonal<R> &r) : Diagonal<T>::Diagonal(r) { } // copy existing column vector
	Column(const uint height, const T &zero) : Diagonal<T>::Diagonal(height, zero) { } // initialize full column vector with given height
	Column(const Buffer<T> &val, const T &zero) : Diagonal<T>::Diagonal(val, zero) { } // initialize full column vector
	Column(const uint height, const Buffer< pair<uint, T> > &val, const T &zero) : Diagonal<T>::Diagonal(height, val, zero) { } // initialize sparse column vector
	virtual ~Column() { }

	// print functions
	void printData() const {
		cout << "Column vector:" << endl;
		Discrete<T>::printData();
	}
	void printShape() const {
		uint i, k;
		Column data(this->m_zero);
		const uint ranks = getMPIranks();
		for(k=0; k<ranks; k++) {
			Discrete<T>::getFromRank(k, 0, data);
			for(i=0; i<data.m_height; i++) {
				uint row = i;
				if(!data.m_full && !data.searchIndex(row, data.m_row, 0, data.m_row.size())) cout << "." << endl;
				else cout << data.m_val[row] << endl;
			}
		}
	}

	// set functions
	Column &setFullOfZeros(const uint height) { Discrete<T>::setFullOfZeros(height); return *this; }
	Column &setFull(const Buffer<T> &val) { Discrete<T>::setFull(val); return *this; }
	Column &setSparse(const uint height, const Buffer< pair<uint, T> > &val = Buffer< pair<uint, T> >()) { Discrete<T>::setSparse(height, val); return *this; }
	template<typename R> Column &setCopy(const Diagonal<R> &r) { Discrete<T>::setCopy(r); return *this; }
	template<typename R> Column &setFunction(const Diagonal<R> &r, T func(const R &)) { Discrete<T>::setFunction(r, func); return *this; }
	template<typename R> Column &setNegation(const Diagonal<R> &r) { Discrete<T>::setNegation(r); return *this; }
	template<typename R> Column &setInverse(const Diagonal<R> &r) { Discrete<T>::setInverse(r); return *this; }
	template<typename L, typename R> Column &setUnion(const Diagonal<L> &l, const Diagonal<R> &r, T func(const L &, const R &)) { Discrete<T>::setUnion(l, r, func); return *this; }
	template<typename L, typename R> Column &setIntersection(const Diagonal<L> &l, const Diagonal<R> &r, T func(const L &, const R &)) { Discrete<T>::setIntersection(l, r, func); return *this; }
	template<typename L, typename R> Column &setPlus(const Diagonal<L> &l, const Diagonal<R> &r) { Discrete<T>::setPlus(l, r); return *this; }
	template<typename L, typename R> Column &setMinus(const Diagonal<L> &l, const Diagonal<R> &r) { Discrete<T>::setMinus(l, r); return *this; }
	template<typename L, typename R> Column &setTimes(const Diagonal<L> &l, const Diagonal<R> &r) { Discrete<T>::setTimes(l, r); return *this; }
	template<typename L, typename R> Column &setTimes(const Sparse<L> &l, const Diagonal<R> &r) {
		if(orMPI(l.m_width != r.m_height)) { Discrete<T>::setEmpty(); return *this; } // the dimensions do not match
		uint i, j, k, n;
		Buffer<T> val(l.m_beg.size(), this->m_zero);
		if(r.m_full) {
			// send terms
			for(i=0; i<l.m_send.size(); ) {
				const uint ranki = l.m_send[i++];
				Buffer<R> vali(j = l.m_send[i++]);
				while(j > 0) vali[--j] = r.m_val[l.m_send[i++]];
				sendMPI(&vali[0], vali.size() * sizeof(R), ranki, 0);
			}
			// local terms
			i = l.m_beg.size();
			j = l.m_col.size();
			while(i > 0) {
				T &sum = val[--i];
				while(j > l.m_beg[i]) { --j; sum += l.m_val[j] * r.m_val[l.m_col[j]]; }
			}
			// receive terms
			for(i=0,n=0; i<l.m_recv.size(); )	{
				const uint ranki = l.m_recv[i++];
				Buffer<R> vali(j = l.m_recv[i++]);
				recvMPI(&vali[0], vali.size() * sizeof(R), ranki, 0);
				while(j > 0) {
					const R &valj = vali[--j];
					for(k=l.m_recv[i++]; k>0; k--,i++,n++) val[l.m_recv[i]] += l.m_rval[n] * valj;
				}
			}
		}
		else {
			// send sparse data
			for(i=0; i<l.m_send.size(); ) {
				const uint ranki = l.m_send[i++];
				uint vals = 0;
				Buffer<pair<uint,R> > vali(l.m_send[i++]);
				for(j=0; j<vali.size(); j++) {
					uint row = l.m_send[i++];
					if(!Discrete<T>::searchIndex(row, r.m_row, 0, r.m_row.size())) continue;
					vali[vals++] = pair<uint,R>(j, r.m_val[row]);
				}
				sendMPI(&vals, sizeof(uint), ranki, 0);
				if(vals > 0) sendMPI(&vali[0], vals * sizeof(pair<uint,R>), ranki, 1);
			}
			// local terms
			i = l.m_beg.size();
			j = l.m_col.size();
			while(i > 0) {
				T &sum = val[--i];
				uint lastcol = r.m_row.size();
				while(j > l.m_beg[i]) {
					uint col = l.m_col[--j];
					if(!Discrete<T>::searchIndex(col, r.m_row, 0, lastcol)) continue;
					sum += l.m_val[j] * r.m_val[col];
					lastcol = col;
				}
			}
			// receive sparse data
			for(i=0,n=0; i<l.m_recv.size(); )	{
				const uint ranki = l.m_recv[i++];
				const uint sizei = l.m_recv[i++];
				uint vals;
				recvMPI(&vals, sizeof(uint), ranki, 0);
				if(vals > 0) {
					Buffer<pair<uint,R> > vali(vals);
					recvMPI(&vali[0], vals * sizeof(pair<uint,R>), ranki, 1);
					for(j=0,vals=0; j<vali.size(); j++,vals++) {
						while(vals < vali[j].first) { n += l.m_recv[i]; i += 1 + l.m_recv[i]; vals++; } // idle for the next
						for(k=l.m_recv[i++]; k>0; k--,i++,n++) val[l.m_recv[i]] += l.m_rval[n] * vali[j].second;
					}
				}
				while(vals < sizei) { n += l.m_recv[i]; i += 1 + l.m_recv[i]; vals++; } // idle for the next
			}
		}
		this->m_height = l.m_height;
		this->m_full = l.m_full;
		this->m_row = l.m_row;
		this->m_val.swap(val);
		return *this;
	}
	template<typename L, typename R> Column &setIntersection(const Sparse<L> &l, const Diagonal<R> &r, void func(T &, const L &, const R &)) {
		if(orMPI(l.m_width != r.m_height)) { Discrete<T>::setEmpty(); return *this; } // the dimensions do not match
		uint i, j, k, n;
		Buffer<T> val(l.m_beg.size(), this->m_zero);
		if(r.m_full) {
			// send terms
			for(i=0; i<l.m_send.size(); ) {
				const uint ranki = l.m_send[i++];
				Buffer<R> vali(j = l.m_send[i++]);
				while(j > 0) vali[--j] = r.m_val[l.m_send[i++]];
				sendMPI(&vali[0], vali.size() * sizeof(R), ranki, 0);
			}
			// local terms
			i = l.m_beg.size();
			j = l.m_col.size();
			while(i > 0) {
				T &sum = val[--i];
				while(j > l.m_beg[i]) { --j; func(sum, l.m_val[j], r.m_val[l.m_col[j]]); }
			}
			// receive terms
			for(i=0,n=0; i<l.m_recv.size(); )	{
				const uint ranki = l.m_recv[i++];
				Buffer<R> vali(j = l.m_recv[i++]);
				recvMPI(&vali[0], vali.size() * sizeof(R), ranki, 0);
				while(j > 0) {
					const R &valj = vali[--j];
					for(k=l.m_recv[i++]; k>0; k--,i++,n++) func(val[l.m_recv[i]], l.m_rval[n], valj);
				}
			}
		}
		else {
			// send sparse data
			for(i=0; i<l.m_send.size(); ) {
				const uint ranki = l.m_send[i++];
				uint vals = 0;
				Buffer<pair<uint,R> > vali(l.m_send[i++]);
				for(j=0; j<vali.size(); j++) {
					uint row = l.m_send[i++];
					if(!Discrete<T>::searchIndex(row, r.m_row, 0, r.m_row.size())) continue;
					vali[vals++] = pair<uint,R>(j, r.m_val[row]);
				}
				sendMPI(&vals, sizeof(uint), ranki, 0);
				if(vals > 0) sendMPI(&vali[0], vals * sizeof(pair<uint,R>), ranki, 1);
			}
			// local terms
			i = l.m_beg.size();
			j = l.m_col.size();
			while(i > 0) {
				T &sum = val[--i];
				uint lastcol = r.m_row.size();
				while(j > l.m_beg[i]) {
					uint col = l.m_col[--j];
					if(!Discrete<T>::searchIndex(col, r.m_row, 0, lastcol)) continue;
					func(sum, l.m_val[j], r.m_val[col]);
					lastcol = col;
				}
			}
			// receive sparse data
			for(i=0,n=0; i<l.m_recv.size(); )	{
				const uint ranki = l.m_recv[i++];
				const uint sizei = l.m_recv[i++];
				uint vals;
				recvMPI(&vals, sizeof(uint), ranki, 0);
				if(vals > 0) {
					Buffer<pair<uint,R> > vali(vals);
					recvMPI(&vali[0], vals * sizeof(pair<uint,R>), ranki, 1);
					for(j=0,vals=0; j<vali.size(); j++,vals++) {
						while(vals < vali[j].first) { n += l.m_recv[i]; i += 1 + l.m_recv[i]; vals++; } // idle for the next
						for(k=l.m_recv[i++]; k>0; k--,i++,n++) func(val[l.m_recv[i]], l.m_rval[n], vali[j].second);
					}
				}
				while(vals < sizei) { n += l.m_recv[i]; i += 1 + l.m_recv[i]; vals++; } // idle for the next
			}
		}
		this->m_height = l.m_height;
		this->m_full = l.m_full;
		this->m_row = l.m_row;
		this->m_val.swap(val);
		return *this;
	}
	template<typename R> Column &setRowLensq(const Sparse<R> &r) {
		uint i, j, k, n;
		Buffer<T> val(r.m_beg.size(), this->m_zero);
		// local terms
		i = r.m_beg.size();
		j = r.m_col.size();
		while(i > 0) {
			--i;
			while(j > r.m_beg[i]) { --j; val[i] += r.m_val[j] * r.m_val[j]; }
		}
		// receive terms
		for(i=0,n=0; i<r.m_recv.size(); ) {
			i++;
			j = r.m_recv[i++];
			while(j > 0) {
				--j;
				for(k=r.m_recv[i++]; k>0; k--,i++,n++) val[r.m_recv[i]] += r.m_rval[n] * r.m_rval[n];
			}
		}
		this->m_height = r.m_height;
		this->m_full = r.m_full;
		this->m_row = r.m_row;
		this->m_val.swap(val);
		return *this;
	}
	template<typename R> Column &setRowMaxAbs(const Sparse<R> &r) {
		uint i, j, k, n;
		Buffer<T> val(r.m_beg.size(), this->m_zero);
		// local terms
		i = r.m_beg.size();
		j = r.m_col.size();
		while(i > 0) {
			--i;
			while(j > r.m_beg[i]) { 
				--j; 
				if(val[i] < abs(r.m_val[j])) val[i] = abs(r.m_val[j]);
			}
		}
		// receive terms
		for(i=0,n=0; i<r.m_recv.size(); ) {
			i++;
			j = r.m_recv[i++];
			while(j > 0) {
				--j;
				for(k=r.m_recv[i++]; k>0; k--,i++,n++) {
					if(val[r.m_recv[i]] < abs(r.m_rval[n])) val[r.m_recv[i]] = abs(r.m_rval[n]);
				}
			}
		}
		this->m_height = r.m_height;
		this->m_full = r.m_full;
		this->m_row = r.m_row;
		this->m_val.swap(val);
		return *this;
	}
	template<typename L, typename R> Column &setScale(const Diagonal<L> &l, const R &r) { Discrete<T>::setScale(l, r); return *this; }
	template<typename L, typename R> Column &setScale(const L &l, const Diagonal<R> &r) { Discrete<T>::setScale(l, r); return *this; }

	// modifier functions
	template<typename R> Column &operator+=(const Diagonal<R> &r) { return setPlus(*this, r); }
	template<typename R> Column &operator-=(const Diagonal<R> &r) { return setMinus(*this, r); }
	template<typename R> Column &scale(const R &r) { return setScale(*this, r); }
	Column &negate() { return setNegation(*this); }

	// trim functions
	Column &trimFull() { Discrete<T>::trimFull(); return *this; } // convert sparse to full
	Column &trimSparse() { Discrete<T>::trimSparse(); return *this; } // convert full to sparse
	Column &trimOptimal(const double limit = 0.5) { Discrete<T>::trimOptimal(limit); return *this; } // convert to full if(number of non-empty rows > limit * m_height), otherwise convert to sparse
	Column &trim() { Discrete<T>::trim(); return *this; } // remove all zero instances

};

// operators
template<typename L, typename R, typename O = decltype(declval<L &>() + declval<R &>())> Column<O> operator+(const Column<L> &l, const Column<R> &r) {
	Column<O> o(l.m_zero + r.m_zero);
	return o.setPlus(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() - declval<R &>())> Column<O> operator-(const Column<L> &l, const Column<R> &r) {
	Column<O> o(l.m_zero - r.m_zero);
	return o.setMinus(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Column<O> operator*(const Diagonal<L> &l, const Column<R> &r) {
	Column<O> o(l.m_zero * r.m_zero);
	return o.setTimes(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Column<O> operator*(const Sparse<L> &l, const Column<R> &r) {
	Column<O> o(l.m_zero * r.m_zero);
	return o.setTimes(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Column<O> scaled(const Column<L> &l, const R &r) {
	Column<O> o(l.m_zero * r);
	return o.setScale(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Column<O> scaled(const L &l, const Column<R> &r) {
	Column<O> o(l * r.m_zero);
	return o.setScale(l, r);
}
template<typename R> Column<R> operator-(const Column<R> &r) {
	Column<R> o(r.m_zero);
	return o.setNegation(r);
}

}

#endif //_COLUMN_HPP_INCLUDED_
