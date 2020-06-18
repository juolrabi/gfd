/**
 * Sparse.hpp implements a sparse matrix operator that operates discrete forms (class Form)
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _SPARSE_HPP_INCLUDED_
#define _SPARSE_HPP_INCLUDED_

#include "Diagonal.hpp"
#include <map>

using namespace std;

namespace gfd
{

template<typename T>
class Sparse : public Discrete<T>
{
public:
	// constructors
	Sparse(const T &zero = 0) : Discrete<T>::Discrete(zero) { m_width = 0; } // initialize empty sparse matrix
	template<typename R> Sparse(const Sparse<R> &r) : Discrete<T>::Discrete(r.m_zero) { setCopy(r); } // copy existing sparse matrix
	template<typename R> Sparse(const Diagonal<R> &r) : Discrete<T>::Discrete(r.m_zero) { setDiagonal(r); } // copy diagonal matrix
	Sparse(const uint width, const Buffer< Buffer< pair<uint, T> > > &val, const T &zero) : Discrete<T>::Discrete(zero) { setFull(width, val); } // initialize full sparse matrix
	Sparse(const uint width, const uint height, const Buffer< pair<uint, Buffer< pair<uint, T> > > > &val, const T &zero) : Discrete<T>::Discrete(zero) { setSparse(width, height, val); } // initialize sparse sparse matrix
	virtual ~Sparse() { }

	// print functions
	void printData() const {
		cout << "Sparse matrix:" << endl;
		cout << "m_width = " << m_width << endl;
		Discrete<T>::printData();
		uint i;
		cout << "m_beg =";
		for(i=0; i<m_beg.size(); i++) cout << " " << m_beg[i];
		cout << endl;
		cout << "m_col =";
		for(i=0; i<m_col.size(); i++) cout << " " << m_col[i];
		cout << endl;
		cout << "m_send =";
		for(i=0; i<m_send.size(); i++) cout << " " << m_send[i];
		cout << endl;
		cout << "m_recv =";
		for(i=0; i<m_recv.size(); i++) cout << " " << m_recv[i];
		cout << endl;
		cout << "m_rval =";
		for(i=0; i<m_rval.size(); i++) cout << " " << m_rval[i];
		cout << endl;
	}
	void printShape() const {
		uint i, j, k, l;
		const uint irank = getMPIrank();
		const uint ranks = getMPIranks();
		uint col0 = 0;
		if(irank > 0) recvMPI(&col0, sizeof(uint), irank - 1, 0);
		uint col1 = col0 + m_width;
		sendMPI(&col1, sizeof(uint), (irank + 1) % ranks, 0);
		uint totalwidth = 0;
		if(irank == 0) recvMPI(&totalwidth, sizeof(uint), ranks - 1, 0);

		Buffer< Buffer< pair<uint,T> > > buf(this->m_height);
		Buffer<uint> row = this->m_row;
		if(this->m_full) {
			row.resize(this->m_height);
			for(i=0; i<row.size(); i++) row[i] = i;
		}
		i = m_beg.size();
		j = m_col.size();
		while(i > 0) {
			--i;
			while(j > m_beg[i]) {
				--j;
				buf[row[i]].push_back(pair<uint,T>(m_col[j] + col0, this->m_val[j]));
			}
		}
		for(i=0; i<m_send.size(); ) {
			const uint ranki = m_send[i++];
			Buffer<uint> col(m_send[i++]);
			for(j=0; j<col.size(); j++,i++) col[j] = m_send[i] + col0;
			sendMPI(&col[0], col.size() * sizeof(uint), ranki, 1);
		}
		for(i=0,l=0; i<m_recv.size(); )	{
			const uint ranki = m_recv[i++];
			Buffer<uint> col(m_recv[i++]);
			recvMPI(&col[0], col.size() * sizeof(uint), ranki, 1);
			for(j=0; j<col.size(); j++) {
				for(k=m_recv[i++]; k>0; k--,i++,l++) {
					buf[row[m_recv[i]]].push_back(pair<uint,T>(col[j], m_rval[l]));
				}
			}
		}
		if(irank > 0) {
			uint num = buf.size();
			sendMPI(&num, sizeof(uint), 0, 2);
			for(i=0; i<num; i++) {
				uint numm = buf[i].size();
				sendMPI(&numm, sizeof(uint), 0, 3);
				sendMPI(&buf[i][0], numm * sizeof(pair<uint,T>), 0, 4);
			}
		}
		else {
			for(i=0; i<buf.size(); i++) {
				for(j=0; j<totalwidth; j++) {
					for(k=0; k<buf[i].size() && buf[i][k].first != j; k++);
					if(k < buf[i].size()) cout << buf[i][k].second << " ";
					else cout << ". ";
				}
				cout << endl;
			}
			for(uint r=1; r<ranks; r++) {
				uint num;
				recvMPI(&num, sizeof(uint), r, 2);
				for(i=0; i<num; i++) {
					uint numm;
					recvMPI(&numm, sizeof(uint), r, 3);
					Buffer< pair<uint,T> > bufi(numm);
					recvMPI(&bufi[0], numm * sizeof(pair<uint,T>), r, 4);
					for(j=0; j<totalwidth; j++) {
						for(k=0; k<bufi.size() && bufi[k].first != j; k++);
						if(k < bufi.size()) cout << bufi[k].second << " ";
						else cout << ". ";
					}
					cout << endl;
				}

			}
		}
	}

	// set functions
	Sparse &setFull(const uint width, const Buffer< Buffer< pair<uint, T> > > &val, const Buffer< pair<uint,uint> > &ext = Buffer< pair<uint,uint> >()) {
		const uint thisRank = getMPIrank();
		Buffer<uint> srank(getMPIranks()-1);
		for(uint i=0; i<srank.size(); i++) srank[i] = (i < thisRank ? i : i + 1);
		return setFull(width, val, ext, srank, srank);
	}
	Sparse &setFull(const uint width, const Buffer< Buffer< pair<uint, T> > > &val, const Buffer< pair<uint,uint> > &ext, const Buffer<uint> &srank, const Buffer<uint> &rrank) {
		uint i, j, k;
		m_width = width;
		this->m_height = val.size();
		this->m_full = true;
		this->m_row.clear();
		m_beg.resize(this->m_height);
		// estimate number of terms
		for(i=0,k=0; i<this->m_height; i++) k += val[i].size();
		m_col.resize(k);
		this->m_val.resize(k);
		// set terms
		Buffer< Buffer< pair<uint, T> > > extv(ext.size());
		for(i=0,k=0; i<this->m_height; i++) {
			m_beg[i] = k;
			const Buffer<uint> ord = Discrete<T>::bubbleSort(val[i]);
			for(j=0; j<ord.size(); j++) {
				m_col[k] = val[i][ord[j]].first;
				this->m_val[k] = val[i][ord[j]].second;
				if(m_col[k] < m_width) k++;
				else extv[m_col[k] - m_width].push_back(pair<uint, T>(i, this->m_val[k]));
			}
		}
		m_col.resize(k);
		this->m_val.resize(k);
		return setCommunication(ext, extv, srank, rrank);
	}
	Sparse &setSparse(const uint width, const uint height, const Buffer< pair<uint, Buffer< pair<uint, T> > > > &val = Buffer< pair<uint, Buffer< pair<uint, T> > > >(), const Buffer< pair<uint,uint> > &ext = Buffer< pair<uint,uint> >()) {
		const uint thisRank = getMPIrank();
		Buffer<uint> srank(getMPIranks()-1);
		for(uint i=0; i<srank.size(); i++) srank[i] = (i < thisRank ? i : i + 1);
		return setSparse(width, height, val, ext, srank, srank);
	}
	Sparse &setSparse(const uint width, const uint height, const Buffer< pair<uint, Buffer< pair<uint, T> > > > &val, const Buffer< pair<uint,uint> > &ext, const Buffer<uint> &srank, const Buffer<uint> &rrank) {
		uint i, j, k, l;
		m_width = width;
		this->m_height = height;
		this->m_full = false;
		const Buffer<uint> rord = Discrete<T>::bubbleSort(val);
		this->m_row.resize(rord.size());
		m_beg.resize(rord.size());
		// estimate number of terms
		for(i=0,k=0; i<rord.size(); i++) k += val[rord[i]].second.size();
		m_col.resize(k);
		this->m_val.resize(k);
		// set terms
		Buffer< Buffer< pair<uint, T> > > extv(ext.size());
		for(i=0,j=0,k=0; i<rord.size(); i++) {
			const Buffer< pair<uint, T> > &vali = val[rord[i]].second;
			const Buffer<uint> ord = Discrete<T>::bubbleSort(vali);
			if(ord.empty()) continue;
			this->m_row[j] = val[rord[i]].first;
			m_beg[j] = k;
			for(l=0; l<ord.size(); l++) {
				m_col[k] = vali[ord[l]].first;
				this->m_val[k] = vali[ord[l]].second;
				if(m_col[k] < m_width) k++;
				else extv[m_col[k] - m_width].push_back(pair<uint, T>(j, this->m_val[k]));
			}
			j++;
		}
		this->m_row.resize(j);
		m_beg.resize(j);
		m_col.resize(k);
		this->m_val.resize(k);
		return setCommunication(ext, extv, srank, rrank);
	}
	template<typename R> Sparse &setCopy(const Sparse<R> &r) {
		if(this == (void*)&r) return *this;
		setShape(r);
		uint i;
		for(i=0; i<this->m_val.size(); i++) this->m_val[i] = r.m_val[i];
		for(i=0; i<m_rval.size(); i++) m_rval[i] = r.m_rval[i];
		return *this;
	}
	template<typename R> Sparse &setDiagonal(const Diagonal<R> &r) {
		setEmpty();
		Discrete<T>::setCopy(r);
		m_width = r.m_height;
		m_beg.resize(r.m_val.size());
		for(uint i=0; i<m_beg.size(); i++) m_beg[i] = i;
		if(r.m_full) m_col = m_beg;
		else m_col = r.m_row;
		return *this;
	}
	template<typename R> Sparse &setColumn(const Diagonal<R> &r) { // härö: this function is wrong
		setEmpty(); 
		Discrete<T>::setCopy(r);
		m_width = 1;
		m_beg.resize(r.m_val.size());
		for(uint i=0; i<m_beg.size(); i++) m_beg[i] = i;
		m_col.resize(m_beg.size());
		m_col.fill(0);
		return *this;
	}
	template<typename R> Sparse &setFunction(const Sparse<R> &r, T func(const R &)) {
		if(this != (void*)&r) setShape(r);
		uint i;
		for(i=0; i<this->m_val.size(); i++) this->m_val[i] = func(r.m_val[i]);
		for(i=0; i<m_rval.size(); i++) m_rval[i] = func(r.m_rval[i]);
		return *this;
	}
	template<typename R> Sparse &setNegation(const Sparse<R> &r) { return setFunction(r, functionNegation); }
	template<typename R> Sparse &setTranspose(const Sparse<R> &r) {
		uint i, j, k, m, n;
		// local terms
		Buffer<uint> cols(r.m_width, 0);
		for(i=0; i<r.m_col.size(); i++) cols[r.m_col[i]]++;
		Buffer<uint> beg(r.m_width);
		for(i=0, k=0; i<beg.size(); i++) {
			beg[i] = k;
			k += cols[i];
			cols[i] = beg[i];
		}
		Buffer<uint> col(r.m_col.size());
		Buffer<T> val(r.m_col.size());
		for(i=0; i<r.m_beg.size(); i++) {
			const uint j0 = r.m_beg[i];
			const uint j1 = (i+1 < r.m_beg.size() ? r.m_beg[i+1] : r.m_col.size());
			const uint coli = (r.m_full ? i : r.m_row[i]);
			for(j=j0; j<j1; j++) {
				const uint ii = cols[r.m_col[j]]++;
				col[ii] = coli;
				val[ii] = r.m_val[j];
			}
		}
		m_beg.swap(beg);
		beg.clear();
		m_col.swap(col);
		col.clear();
		this->m_val.swap(val);
		val.clear();

		// initialize communication
		uint sends = 0;
		Buffer<uint> send;
		Buffer<uint> comm;
		for(i=0,n=0; i<r.m_recv.size(); ) {
			const uint ranki = r.m_recv[i++];
			uint comms = 0;
			for(j=r.m_recv[i++]; j>0; j--) {
				uint num = r.m_recv[i++];
				sendMPI(&num, sizeof(uint), ranki, 0);
				Buffer< pair<uint,T> > buf(num);
				for(k=0; k<num; k++,i++,n++) {
					buf[k] = pair<uint,T>(r.m_recv[i], r.m_rval[n]);
					comm.gatherOnce(r.m_recv[i], comms);
				}
				sendMPI(&buf[0], num * sizeof(pair<uint,T>), ranki, 1);
			}
			// insert terms into send
			send.gather(ranki, sends); // gather rank
			send.gather(comms, sends); // gather number of communications
			for(m=0; m<comms; m++) send.gather((r.m_full ? comm[m] : r.m_row[comm[m]]), sends);
		}
		uint recvs = 0;
		Buffer<uint> recv;
		uint rvals = 0;
		Buffer<T> rval;
		for(i=0; i<r.m_send.size(); ) {
			const uint ranki = r.m_send[i++];
			uint comms = 0;
			Buffer< Buffer< pair<uint,T> > > buf(r.m_send[i++]);
			for(j=0; j<buf.size(); j++) {
				uint bufs;
				recvMPI(&bufs, sizeof(uint), ranki, 0);
				buf[j].resize(bufs);
				recvMPI(&buf[j][0], bufs * sizeof(pair<uint,T>), ranki, 1);
				for(k=0; k<bufs; k++) comm.gatherOnce(buf[j][k].first, comms);
			}
			// insert terms to recv
			recv.gather(ranki, recvs); // gather rank
			recv.gather(comms, recvs); // gather number of communications
			for(m=0; m<comms; m++) {
				const uint recvs0 = recvs;
				recv.gather(0, recvs); // gather number of copies (set to zero and update later)
				for(j=0; j<buf.size(); j++) {
					for(k=0; k<buf[j].size(); k++) {
						if(buf[j][k].first == comm[m]) {
							recv.gather(r.m_send[i + j], recvs); // gather row to copy to
							rval.gather(buf[j][k].second, rvals); // gather term value
						}
					}
				}
				recv[recvs0] = recvs - recvs0 - 1; // update number of copies
			}
			i += buf.size();
		}
		m_send.copy(send, sends);
		send.clear();
		m_recv.copy(recv, recvs);
		recv.clear();
		m_rval.copy(rval, rvals);
		rval.clear();
		comm.clear();
		// initialize matrix shape
		const bool rfull = r.m_full;
		const uint rwidth = r.m_width;
		m_width = r.m_height;
		this->m_height = rwidth;
		this->m_full = true;
		this->m_row.clear();
		if(rfull) return *this;
		return toSparse();
	}
	template<typename L, typename R> Sparse &setUnion(const Sparse<L> &l, const Sparse<R> &r, T func(const L &, const R &)) {
		if(orMPI(l.m_width != r.m_width || l.m_height != r.m_height)) return setEmpty(); // the dimensions do not match
		m_width = l.m_width;
		this->m_height = l.m_height;
		Buffer<uint> beg;
		Buffer<uint> row;
		uint cols = 0;
		Buffer<uint> col(l.m_col.size() + r.m_col.size());
		Buffer<T> val(col.size());
		Buffer<uint> lto(l.m_beg.size());
		Buffer<uint> rto(r.m_beg.size());
		if(l.m_full) {
			beg.resize(l.m_height);
			if(r.m_full) { // both l and r are full
				for(uint i=0; i<l.m_height; i++) {
					beg[i] = cols;
					lto[i] = i;
					rto[i] = i;
					mergeRow(l, i, r, i, val, col, cols, func);
				}
			}
			else { // only l is full
				uint i = 0, ri = 0;
				while(ri < r.m_row.size()) {
					while(i < r.m_row[ri]) {
						beg[i] = cols;
						lto[i] = i;
						mergeRow(l, i++, r, uint(-1), val, col, cols, func);
					}
					beg[i] = cols;
					lto[i] = i;
					rto[ri] = i;
					mergeRow(l, i++, r, ri++, val, col, cols, func);
				}
				while(i < l.m_height) {
					beg[i] = cols;
					lto[i] = i;
					mergeRow(l, i++, r, uint(-1), val, col, cols, func);
				}
			}
			this->m_full = true;
		}
		else if(r.m_full) { // only r is full
			beg.resize(r.m_height);
			uint i=0, li = 0;
			while(li < l.m_row.size()) {
				while(i < l.m_row[li])  {
					beg[i] = cols;
					rto[i] = i;
					mergeRow(l, uint(-1), r, i++, val, col, cols, func);
				}
				beg[i] = cols;
				lto[li] = i;
				rto[i] = i;
				mergeRow(l, li++, r, i++, val, col, cols, func);
			}
			while(i < r.m_height) {
				beg[i] = cols;
				rto[i] = i;
				mergeRow(l, uint(-1), r, i++, val, col, cols, func);
			}
			this->m_full = true;
		}
		else { // both l and r are sparse
			uint rowsize = l.m_row.size() + r.m_row.size();
			if(rowsize > l.m_height) rowsize = l.m_height;
			row.resize(rowsize);
			beg.resize(rowsize);
			uint i = 0, li = 0, ri = 0;
			while(li < l.m_row.size() && ri < r.m_row.size()) {
				beg[i] = cols;
				if(l.m_row[li] < r.m_row[ri]) {
					lto[li] = i;
					row[i++] = l.m_row[li];
					mergeRow(l, li++, r, uint(-1), val, col, cols, func);
				}
				else if(l.m_row[li] > r.m_row[ri]) {
					rto[ri] = i;
					row[i++] = r.m_row[ri];
					mergeRow(l, uint(-1), r, ri++, val, col, cols, func);
				}
				else {
					lto[li] = i;
					rto[ri] = i;
					row[i++] = l.m_row[li];
					mergeRow(l, li++, r, ri++, val, col, cols, func);
				}
			}
			while(li < l.m_row.size()) {
				beg[i] = cols;
				lto[li] = i;
				row[i++] = l.m_row[li];
				mergeRow(l, li++, r, uint(-1), val, col, cols, func);
			}
			while(ri < r.m_row.size()) {
				beg[i] = cols;
				rto[ri] = i;
				row[i++] = r.m_row[ri];
				mergeRow(l, uint(-1), r, ri++, val, col, cols, func);
			}
			beg.resize(i);
			row.resize(i);
			this->m_full = false;
		}
		val.resize(cols);
		col.resize(cols);
		this->m_row.swap(row);
		m_beg.swap(beg);
		m_col.swap(col);
		this->m_val.swap(val);
		row.clear();
		beg.clear();
		col.clear();
		val.clear();
		// initialize m_send
		uint i, j, k, m, n;
		uint iranks = 0;
		Buffer<uint> irank;
		for(i=0; i<l.m_send.size(); ) {
			const uint ranki = l.m_send[i++];
			irank.gather(ranki, iranks);
			sendMPI(&l.m_send[i+1], l.m_send[i] * sizeof(uint), ranki, 0);
			i += 1 + l.m_send[i];
		}
		for(i=0; i<r.m_send.size(); ) {
			const uint ranki = r.m_send[i++];
			irank.gatherOnce(ranki, iranks);
			sendMPI(&r.m_send[i+1], r.m_send[i] * sizeof(uint), ranki, 1);
			i += 1 + r.m_send[i];
		}
		Buffer<uint> sends(iranks, 0);
		Buffer< Buffer<uint> > send(iranks);
		for(i=0; i<l.m_send.size(); ) {
			const uint ranki = l.m_send[i++];
			const uint ii = irank.findFirst(ranki, iranks);
			for(j=l.m_send[i++]; j>0; j--) send[ii].gather(l.m_send[i++], sends[ii]);
		}
		for(i=0; i<r.m_send.size(); ) {
			const uint ranki = r.m_send[i++];
			const uint ii = irank.findFirst(ranki, iranks);
			for(j=r.m_send[i++]; j>0; j--) send[ii].gatherOnce(r.m_send[i++], sends[ii]);
		}
		for(i=0,m=0; i<iranks; i++) m += 2 + sends[i];
		m_send.clear();
		m_send.resize(m);
		for(i=0,m=0; i<iranks; i++) {
			m_send[m++] = irank[i];
			m_send[m++] = sends[i];
			for(j=0; j<sends[i]; j++) m_send[m++] = send[i][j];
		}
		send.clear();
		// initialize m_recv and m_rval
		iranks = 0;
		uint terms = 0;
		for(i=0; i<l.m_recv.size(); ) {
			const uint ranki = l.m_recv[i++];
			irank.gather(ranki, iranks);
			terms += l.m_recv[i];
			for(j=l.m_recv[i++]; j>0; j--) i += 1 + l.m_recv[i];
		}
		for(i=0; i<r.m_recv.size(); ) {
			const uint ranki = r.m_recv[i++];
			irank.gatherOnce(ranki, iranks);
			terms += r.m_recv[i];
			for(j=r.m_recv[i++]; j>0; j--) i += 1 + r.m_recv[i];
		}
		Buffer<uint> recvs(iranks, 0);
		Buffer< Buffer<pair<uint,uint> > > recv(iranks);
		Buffer< Buffer<pair<uint,pair<L,R> > > > term(terms);
		terms = 0;
		for(i=0,n=0; i<l.m_recv.size(); ) {
			const uint ranki = l.m_recv[i++];
			const uint ii = irank.findFirst(ranki, iranks);
			Buffer<uint> buf(l.m_recv[i++]);
			recvMPI(&buf[0], buf.size() * sizeof(uint), ranki, 0);
			for(j=0; j<buf.size(); j++) {
				recv[ii].gather(pair<uint,uint>(buf[j], terms), recvs[ii]);
				Buffer<pair<uint,pair<L,R> > > &termj = term[terms++];
				termj.resize(l.m_recv[i++]);
				for(k=0; k<termj.size(); k++,i++,n++) termj[k] = pair<uint,pair<L,R> >(lto[l.m_recv[i]], pair<L,R>(l.m_rval[n], r.m_zero));
			}
		}
		for(i=0,n=0; i<r.m_recv.size(); ) {
			const uint ranki = r.m_recv[i++];
			const uint ii = irank.findFirst(ranki, iranks);
			Buffer<uint> buf(r.m_recv[i++]);
			recvMPI(&buf[0], buf.size() * sizeof(uint), ranki, 1);
			for(j=0; j<buf.size(); j++) {
				for(k=0; k<recvs[ii]; k++) {
					if(recv[ii][k].first == buf[j]) break;
				}
				if(k < recvs[ii]) { // insert into existing term
					Buffer<pair<uint,pair<L,R> > > &termj = term[recv[ii][k].second];
					for(k=r.m_recv[i++]; k>0; k--,i++, n++) {
						const uint rtoi = rto[r.m_recv[i]];
						for(m=0; m<termj.size() && termj[m].first!=rtoi; m++);
						if(m < termj.size()) termj[m].second.second = r.m_rval[n];
						else termj.push_back(pair<uint,pair<L,R> >(rtoi, pair<L,R>(l.m_zero, r.m_rval[n])));
					}
				}
				else { // create new term
					recv[ii].gather(pair<uint,uint>(buf[j], terms), recvs[ii]);
					Buffer<pair<uint,pair<L,R> > > &termj = term[terms++];
					termj.resize(r.m_recv[i++]);
					for(k=0; k<termj.size(); k++,i++,n++) termj[k] = pair<uint,pair<L,R> >(rto[r.m_recv[i]], pair<L,R>(l.m_zero, r.m_rval[n]));
				}
			}
		}
		for(i=0,m=0,n=0; i<iranks; i++) {
			m += 2;
			for(j=0; j<recvs[i]; j++) {
				const uint sizej = term[recv[i][j].second].size();
				m += 1 + sizej;
				n += sizej;
			}
		}
		m_recv.clear();
		m_recv.resize(m);
		m_rval.clear();
		m_rval.resize(n);
		for(i=0,m=0,n=0; i<iranks; i++) {
			m_recv[m++] = irank[i];
			m_recv[m++] = recvs[i];
			for(j=0; j<recvs[i]; j++) {
				const Buffer<pair<uint,pair<L,R> > > &termj = term[recv[i][j].second];
				m_recv[m++] = termj.size();
				for(k=0; k<termj.size(); k++,m++,n++) {
					m_recv[m] = termj[k].first;
					m_rval[n] = func(termj[k].second.first, termj[k].second.second);
				}
			}
		}
		recv.clear();
		term.clear();
		return *this;
	}
	template<typename L, typename R> Sparse &setUnion(const Sparse<L> &l, const Diagonal<R> &r, T func(const L &, const R &)) {
		if(orMPI(l.m_width != r.m_height || l.m_height != r.m_height)) return setEmpty(); // the dimensions do not match
		m_width = l.m_width;
		this->m_height = l.m_height;
		Buffer<uint> beg;
		Buffer<uint> row;
		uint cols = 0;
		Buffer<uint> col(l.m_col.size() + r.m_val.size());
		Buffer<T> val(col.size());
		Buffer<uint> lto(l.m_beg.size());
		if(l.m_full) {
			beg.resize(l.m_height);
			if(r.m_full) { // both l and r are full
				for(uint i=0; i<l.m_height; i++) {
					beg[i] = cols;
					lto[i] = i;
					mergeRow(l, i, r, i, val, col, cols, func);
				}
			}
			else { // only l is full
				uint i = 0, ri = 0;
				while(ri < r.m_row.size()) {
					while(i < r.m_row[ri]) {
						beg[i] = cols;
						lto[i] = i;
						mergeRow(l, i++, r, uint(-1), val, col, cols, func);
					}
					beg[i] = cols;
					lto[i] = i;
					mergeRow(l, i++, r, ri++, val, col, cols, func);
				}
				while(i < l.m_height) {
					beg[i] = cols;
					lto[i] = i;
					mergeRow(l, i++, r, uint(-1), val, col, cols, func);
				}
			}
			this->m_full = true;
		}
		else if(r.m_full) { // only r is full
			beg.resize(r.m_height);
			uint i=0, li = 0;
			while(li < l.m_row.size()) {
				while(i < l.m_row[li])  {
					beg[i] = cols;
					mergeRow(l, uint(-1), r, i++, val, col, cols, func);
				}
				beg[i] = cols;
				lto[li] = i;
				mergeRow(l, li++, r, i++, val, col, cols, func);
			}
			while(i < r.m_height) {
				beg[i] = cols;
				mergeRow(l, uint(-1), r, i++, val, col, cols, func);
			}
			this->m_full = true;
		}
		else { // both l and r are sparse
			uint rowsize = l.m_row.size() + r.m_row.size();
			if(rowsize > l.m_height) rowsize = l.m_height;
			row.resize(rowsize);
			beg.resize(rowsize);
			uint i = 0, li = 0, ri = 0;
			while(li < l.m_row.size() && ri < r.m_row.size()) {
				beg[i] = cols;
				if(l.m_row[li] < r.m_row[ri]) {
					lto[li] = i;
					row[i++] = l.m_row[li];
					mergeRow(l, li++, r, uint(-1), val, col, cols, func);
				}
				else if(l.m_row[li] > r.m_row[ri]) {
					row[i++] = r.m_row[ri];
					mergeRow(l, uint(-1), r, ri++, val, col, cols, func);
				}
				else {
					lto[li] = i;
					row[i++] = l.m_row[li];
					mergeRow(l, li++, r, ri++, val, col, cols, func);
				}
			}
			while(li < l.m_row.size()) {
				beg[i] = cols;
				lto[li] = i;
				row[i++] = l.m_row[li];
				mergeRow(l, li++, r, uint(-1), val, col, cols, func);
			}
			while(ri < r.m_row.size()) {
				beg[i] = cols;
				row[i++] = r.m_row[ri];
				mergeRow(l, uint(-1), r, ri++, val, col, cols, func);
			}
			beg.resize(i);
			row.resize(i);
			this->m_full = false;
		}
		val.resize(cols);
		col.resize(cols);
		this->m_row.swap(row);
		m_beg.swap(beg);
		m_col.swap(col);
		this->m_val.swap(val);
		row.clear();
		beg.clear();
		col.clear();
		val.clear();
		// initialize send, recv and rval
		if(this != (void*)&l) {
			m_send = l.m_send;
			m_recv.resize(l.m_recv.size());
			m_rval.resize(l.m_rval.size());
		}
		uint i, j, k;
		for(i=0; i<m_recv.size(); )	{
			m_recv[i] = l.m_recv[i];
			i++;
			m_recv[i] = l.m_recv[i];
			for(j=m_recv[i++]; j>0; j--) {
				m_recv[i] = l.m_recv[i];
				for(k=m_recv[i++]; k>0; k--,i++) {
					m_recv[i] = lto[l.m_recv[i]];
				}
			}
		}
		for(i=0; i<m_rval.size(); i++) m_rval[i] = func(l.m_rval[i], r.m_zero);
		return *this;
	}
	template<typename L, typename R> Sparse &setPlus(const Sparse<L> &l, const Sparse<R> &r) { return setUnion(l, r, functionPlus); }
	template<typename L, typename R> Sparse &setMinus(const Sparse<L> &l, const Sparse<R> &r) { return setUnion(l, r, functionMinus); }
	template<typename L, typename R> Sparse &setPlus(const Sparse<L> &l, const Diagonal<R> &r) { return setUnion(l, r, functionPlus); }
	template<typename L, typename R> Sparse &setMinus(const Sparse<L> &l, const Diagonal<R> &r) { return setUnion(l, r, functionMinus); }
	template<typename L, typename R> Sparse &setPlus(const Diagonal<L> &l, const Sparse<R> &r) { return setUnion(r, l, functionReversePlus); }
	template<typename L, typename R> Sparse &setMinus(const Diagonal<L> &l, const Sparse<R> &r) { return setUnion(r, l, functionReverseMinus); }
	template<typename L, typename R> Sparse &setTimes(const Sparse<L> &l, const Sparse<R> &r) {
		if(orMPI(l.m_width != r.m_height)) return setEmpty(); // the dimensions do not match
		uint i, j, k, n, m;
		// send external links
		uint sranks = 0;
		Buffer<uint> srank;
		for(i=0; i<r.m_send.size(); ) {
			const uint ranki = r.m_send[i++];
			srank.gather(ranki, sranks);
			sendMPI(&r.m_send[i+1], r.m_send[i] * sizeof(uint), ranki, 0);
			i += 1 + r.m_send[i];
		}
		// receive external links
		uint rranks = 0;
		Buffer<uint> rrank;
		Buffer< Buffer< pair<pair<uint,uint>,R> > > rext(r.m_beg.size());
		for(i=0,n=0; i<r.m_recv.size(); ) {
			const uint ranki = r.m_recv[i++];
			rrank.gather(ranki, rranks);
			Buffer<uint> link(r.m_recv[i++]);
			recvMPI(&link[0], link.size() * sizeof(uint), ranki, 0);
			for(j=0; j<link.size(); j++) {
				for(k=r.m_recv[i++]; k>0; k--,i++,n++) {
					rext[r.m_recv[i]].push_back(pair<pair<uint,uint>,R>(pair<uint,uint>(ranki, link[j]), r.m_rval[n]));
				}
			}
		}
		// send external rows
		Buffer<uint> cranks(rranks, 0); // number of concatenated ranks
		Buffer< Buffer<uint> > crank(rranks); // concatenated ranks
		for(i=0; i<l.m_send.size(); ) {
			const uint ranki = l.m_send[i++];
			srank.gatherOnce(ranki, sranks);
			for(j=l.m_send[i++]; j>0; j--) {
				uint rrow = l.m_send[i++];
				if(!r.m_full && !Discrete<T>::searchIndex(rrow, r.m_row, 0, r.m_row.size())) {
					uint zero = 0;
					sendMPI(&zero, sizeof(uint), ranki, 0);
					continue;
				}
				const uint rbeg = r.m_beg[rrow];
				uint rsize = (rrow+1<r.m_beg.size() ? r.m_beg[rrow+1] : r.m_col.size()) - rbeg;
				const Buffer< pair<pair<uint,uint>,R> > &rrext = rext[rrow];
				uint tsize = rsize + rrext.size(); // total size = number of local terms + number of external terms
				sendMPI(&tsize, sizeof(uint), ranki, 0);
				if(tsize == 0) continue;
				sendMPI(&rsize, sizeof(uint), ranki, 1);
				Buffer< pair<uint,R> > rrloc(rsize);
				for(m=0; m<rsize; m++) rrloc[m] = pair<uint,R>(r.m_col[rbeg+m], r.m_val[rbeg+m]);
				sendMPI(&rrloc[0], rsize * sizeof(pair<uint,R>), ranki, 2);
				sendMPI(&rrext[0], rrext.size() * sizeof(pair<pair<uint,uint>,R>), ranki, 3);
				for(m=0; m<rrext.size(); m++) {
					const uint irank = rrext[m].first.first;
					if(irank == ranki) continue;
					const uint ii = rrank.findFirst(irank);
					crank[ii].gatherOnce(ranki, cranks[ii]);
				}
			}
		}
		// send concatenated neighbors
		for(i=0; i<cranks.size(); i++) {
			sendMPI(&cranks[i], sizeof(uint), rrank[i], 2);
			if(cranks[i] == 0) continue;
			sendMPI(&crank[i][0], cranks[i] * sizeof(uint), rrank[i], 3);
		}
		// local terms
		Buffer< Buffer< pair<uint,T> > > buf(l.m_beg.size());
		map<pair<uint,uint>,uint> ext;
		for(i=0; i<buf.size(); i++) {
			Buffer< pair<uint,T> > &bufi = buf[i];
			const uint jend = (i+1<l.m_beg.size() ? l.m_beg[i+1] : l.m_col.size());
			for(j=l.m_beg[i]; j<jend; j++) {
				uint rrow = l.m_col[j];
				if(!r.m_full && !Discrete<T>::searchIndex(rrow, r.m_row, 0, r.m_row.size())) continue;
				const L &lval = l.m_val[j];
				const uint rbeg = r.m_beg[rrow];
				const uint rsize = (rrow+1<r.m_beg.size() ? r.m_beg[rrow+1] : r.m_col.size()) - rbeg;
				for(k=0; k<rsize; k++) sumTerm(pair<uint,T>(r.m_col[rbeg+k], lval * r.m_val[rbeg+k]), bufi);

				const Buffer< pair<pair<uint,uint>,R> > &rrext = rext[rrow];
				for(k=0; k<rrext.size(); k++) sumTerm(getTerm(rrext[k].first, lval * rrext[k].second, r.m_width, ext), bufi);
			}
		}
		// receive external rows
		const uint thisRank = getMPIrank();
		for(i=0,n=0; i<l.m_recv.size(); ) {
			const uint ranki = l.m_recv[i++];
			rrank.gatherOnce(ranki, rranks);
			for(j=l.m_recv[i++]; j>0; j--) {
				uint tsize, rsize;
				recvMPI(&tsize, sizeof(uint), ranki, 0);
				if(tsize == 0) {
					for(k=l.m_recv[i++]; k>0; k--, i++, n++);
					continue;
				}
				recvMPI(&rsize, sizeof(uint), ranki, 1);
				Buffer< pair<uint,R> > rrloc(rsize);
				recvMPI(&rrloc[0], rsize * sizeof(pair<uint,R>), ranki, 2);
				Buffer< pair<pair<uint,uint>,R> > rrext(tsize - rsize);
				recvMPI(&rrext[0], rrext.size() * sizeof(pair<pair<uint,uint>,R>), ranki, 3);
				for(k=l.m_recv[i++]; k>0; k--, i++, n++) {
					const L &lval = l.m_rval[n];
					Buffer< pair<uint,T> > &bufi = buf[l.m_recv[i]];
					for(m=0; m<rsize; m++) sumTerm(getTerm(pair<uint,uint>(ranki,rrloc[m].first), lval * rrloc[m].second, r.m_width, ext), bufi);
					for(m=0; m<rrext.size(); m++) sumTerm(getTerm(rrext[m].first, lval * rrext[m].second, r.m_width, ext), bufi);
				}
				for(m=0; m<rrext.size(); m++) {
					const uint irank = rrext[m].first.first;
					if(irank == thisRank) continue;
					rrank.gatherOnce(irank, rranks);
				}
			}
		}
		Buffer<pair<uint,uint> > bext(ext.size());
		for(map<pair<uint,uint>,uint>::iterator it = ext.begin(); it != ext.end(); it++) bext[it->second - r.m_width] = it->first;
		ext.clear();
		// receive concatenated neighbors
		for(i=0; i<r.m_send.size(); ) {
			const uint ranki = r.m_send[i++];
			i += r.m_send[i] + 1;
			uint iranks;
			recvMPI(&iranks, sizeof(uint), ranki, 2);
			if(iranks == 0) continue;
			Buffer<uint> irank(iranks);
			recvMPI(&irank[0], iranks * sizeof(uint), ranki, 3);
			for(j=0; j<iranks; j++) srank.gatherOnce(irank[j], sranks);
		}
		srank.resize(sranks);
		rrank.resize(rranks);
		// initialize the product matrix
		if(l.m_full) setFull(r.m_width, buf, bext, srank, rrank);
		else {
			Buffer< Buffer< pair<uint,T> > > fbuf(l.m_height);
			for(i=0; i<buf.size(); i++) fbuf[l.m_row[i]].swap(buf[i]);
			setFull(r.m_width, fbuf, bext, srank, rrank);
		}
		if(l.m_full && r.m_full) return *this;
		return toSparse();
	}
	template<typename L, typename R> Sparse &setTimes(const Diagonal<L> &l, const Sparse<R> &r) {
		if(orMPI(l.m_height != r.m_height)) return setEmpty(); // the dimensions do not match
		uint i, j, k, n;
		if(l.m_full) { // l is a full diagonal matrix
			if(this != (void*)&r) setShape(r); // this and r are not equal -> initialize this
			i = m_beg.size();
			j = m_col.size();
			while(i > 0) {
				--i;
				const uint row = (r.m_full ? i : r.m_row[i]);
				while(j > r.m_beg[i]) { --j; this->m_val[j] = l.m_val[row] * r.m_val[j]; }
			}
			for(i=0,n=0; i<r.m_recv.size(); )	{
				i++;
				for(j=r.m_recv[i++]; j>0; j--) {
					for(k=r.m_recv[i++]; k>0; k--,i++,n++) {
						const uint row = (r.m_full ? r.m_recv[i] : r.m_row[r.m_recv[i]]);
						m_rval[n] = l.m_val[row] * r.m_rval[n];
					}
				}
			}
			return *this;
		}
		// l is a sparse diagonal matrix
		if(this != (void*)&r) { // this and r are not equal -> initialize this
			this->m_full = true;
			m_width = r.m_width;
			this->m_height = r.m_height;
			m_beg.resize(r.m_beg.size());
			m_col.resize(r.m_col.size());
			this->m_val.resize(r.m_val.size());
			m_send.resize(r.m_send.size());
			m_recv.resize(r.m_recv.size());
			m_rval.resize(r.m_rval.size());
		}
		uint rows = 0;
		uint cols = 0;
		if(this->m_full) this->m_row.resize(l.m_row.size());
		Buffer<uint> lto(r.m_beg.size(), uint(-1));
		Buffer<uint> rto(r.m_beg.size());
		for(i=0,j=0; i<l.m_row.size(); i++) {
			const uint row = l.m_row[i];
			if(r.m_full) j = row;
			else {
				while(j < r.m_row.size() && r.m_row[j] < row) j++;
				if(j >= r.m_row.size() || r.m_row[j] > row) continue;
			}
			lto[j] = i;
			rto[j] = rows;
			this->m_row[rows] = row;
			m_beg[rows] = cols;
			rows++;
			const uint k0 = r.m_beg[j];
			const uint k1 = (j+1<r.m_beg.size() ? r.m_beg[j+1] : r.m_col.size());
			for(k=k0; k<k1; k++) {
				m_col[cols] = r.m_col[k];
				this->m_val[cols] = l.m_val[i] * r.m_val[k];
				cols++;
			}
		}
		this->m_row.resize(rows);
		m_beg.resize(rows);
		m_col.resize(cols);
		this->m_val.resize(cols);
		// modify and trim m_recv and m_rval
		uint recvs = 0;
		uint rvals = 0;
		for(i=0,n=0; i<r.m_recv.size(); )	{
			const uint ranki = r.m_recv[i++];
			uint vals = 0;
			Buffer<uint> val(r.m_recv[i++]);
			const uint recvRanki = recvs++;
			const uint recvVals = recvs++;
			for(j=0; j<val.size(); j++) {
				uint copys = 0;
				const uint recvCopys = recvs++;
				for(k=r.m_recv[i++]; k>0; k--,i++,n++) {
					if(lto[r.m_recv[i]] == uint(-1)) continue;
					m_rval[rvals++] = l.m_val[lto[r.m_recv[i]]] * r.m_rval[n];
					m_recv[recvs++] = rto[r.m_recv[i]];
					copys++;
				}
				if(copys > 0) {
					m_recv[recvCopys] = copys;
					val[vals++] = j;
				}
				else recvs = recvCopys;
			}
			sendMPI(&vals, sizeof(uint), ranki, 0);
			if(vals > 0) {
				sendMPI(&val[0], vals * sizeof(uint), ranki, 1);
				m_recv[recvRanki] = ranki;
				m_recv[recvVals] = vals;
			}
			else recvs = recvRanki;
		}
		m_recv.resize(recvs);
		m_rval.resize(rvals);
		lto.clear();
		rto.clear();
		// trim m_send
		uint sends = 0;
		for(i=0; i<r.m_send.size(); ) {
			const uint ranki = r.m_send[i++];
			const uint sizei = r.m_send[i++];
			uint vals;
			recvMPI(&vals, sizeof(uint), ranki, 0);
			if(vals > 0) {
				m_send[sends++] = ranki;
				m_send[sends++] = vals;
				Buffer<uint> val(vals);
				recvMPI(&val[0], vals * sizeof(uint), ranki, 1);
				for(j=0; j<val.size(); j++) m_send[sends++] = r.m_send[i + val[j]];
			}
			i += sizei;
		}
		m_send.resize(sends);
		this->m_full = false;
		return *this;
	}
	template<typename L, typename R> Sparse &setTimes(const Sparse<L> &l, const Diagonal<R> &r) {
		if(orMPI(l.m_width != r.m_height)) return setEmpty(); // the dimensions do not match
		uint i, j, k, n;
		if(r.m_full) { // the diagonal matrix is full
			if(this != (void*)&l) setShape(l); // this and l are not equal -> initialize this
			// send full data
			for(i=0; i<m_send.size(); ) {
				const uint ranki = m_send[i++];
				Buffer<R> val(j = m_send[i++]);
				while(j > 0) val[--j] = r.m_val[m_send[i++]];
				sendMPI(&val[0], val.size() * sizeof(R), ranki, 0);
			}
			// local terms
			for(i=0; i<m_col.size(); i++) this->m_val[i] = l.m_val[i] * r.m_val[m_col[i]];
			// receive terms
			for(i=0,n=0; i<m_recv.size(); )	{
				const uint ranki = m_recv[i++];
				Buffer<R> val(j = m_recv[i++]);
				recvMPI(&val[0], val.size() * sizeof(R), ranki, 0);
				while(j > 0) {
					const R &valj = val[--j];
					for(k=m_recv[i++]; k>0; k--,i++,n++) m_rval[n] = l.m_rval[n] * valj;
				}
			}
			return *this;
		}
		// the diagonal matrix is sparse
		if(this != (void*)&l) { // this and l are not equal -> initialize this
			this->m_full = l.m_full;
			m_width = l.m_width;
			this->m_height = l.m_height;
			this->m_row = l.m_row;
			m_beg.resize(l.m_beg.size());
			m_col.resize(l.m_col.size());
			this->m_val.resize(l.m_val.size());
			m_send.resize(l.m_send.size());
			m_recv.resize(l.m_recv.size());
			m_rval.resize(l.m_rval.size());
		}
		// send sparse data and reshape m_send
		uint sends = 0;
		for(i=0,n=0; i<m_send.size(); ) {
			const uint ranki = l.m_send[i++];
			uint vals = 0;
			Buffer<pair<uint,R> > val(l.m_send[i++]);
			const uint sendRanki = sends++;
			const uint sendVals = sends++;
			for(j=0; j<val.size(); j++,i++) {
				uint row = l.m_send[i];
				if(!Discrete<T>::searchIndex(row, r.m_row, 0, r.m_row.size())) continue;
				val[vals++] = pair<uint,R>(j, r.m_val[row]);
				m_send[sends++] = l.m_send[i];
			}
			sendMPI(&vals, sizeof(uint), ranki, 0);
			if(vals > 0) {
				sendMPI(&val[0], vals * sizeof(pair<uint,R>), ranki, 1);
				m_send[sendRanki] = ranki;
				m_send[sendVals] = vals;
			}
			else sends = sendRanki; // undo insert of ranki
		}
		m_send.resize(sends);
		// local terms
		for(i=0,j=0,k=0; i<m_col.size(); i++) {
			while(k<m_beg.size() && l.m_beg[k]<=i) m_beg[k++] = j;
			uint col = l.m_col[i];
			if(!Discrete<T>::searchIndex(col, r.m_row, 0, r.m_row.size())) continue;
			m_col[j] = l.m_col[i];
			this->m_val[j++] = l.m_val[i] * r.m_val[col];
		}
		while(k<m_beg.size()) m_beg[k++] = j;
		m_col.resize(j);
		this->m_val.resize(j);
		// receive sparse data and reshape m_recv and m_rval
		uint recvs = 0;
		uint rvals = 0;
		for(i=0,n=0; i<m_recv.size(); )	{
			const uint ranki = l.m_recv[i++];
			const uint sizei = l.m_recv[i++];
			uint vals;
			recvMPI(&vals, sizeof(uint), ranki, 0);
			if(vals > 0) {
				m_recv[recvs++] = ranki;
				m_recv[recvs++] = vals;
				Buffer<pair<uint,R> > val(vals);
				recvMPI(&val[0], vals * sizeof(pair<uint,R>), ranki, 1);
				for(j=0,vals=0; j<val.size(); j++,vals++) {
					while(vals < val[j].first) { n += l.m_recv[i]; i += 1 + l.m_recv[i]; vals++; } // idle for the next
					m_recv[recvs++] = l.m_recv[i];
					for(k=l.m_recv[i++]; k>0; k--,i++,n++) {
						m_recv[recvs++] = l.m_recv[i];
						m_rval[rvals++] = l.m_rval[n] * val[j].second;
					}
				}
			}
			while(vals < sizei) { n += l.m_recv[i]; i += 1 + l.m_recv[i]; vals++; } // idle for the next
		}
		m_recv.resize(recvs);
		m_rval.resize(rvals);
		return toSparse();
	}
	template<typename L, typename R> Sparse &setOuter(const Diagonal<L> &l, const Diagonal<R> &r) {
		uint i, j, k;
		const uint irank = getMPIrank();
		const uint ranks = getMPIranks();
		// send r data to all neighbors
		for(i=0; i<ranks; i++) {
			if(i == irank) continue;
			j = r.m_val.size();
			sendMPI(&j, sizeof(uint), i, 0);
			if(j == 0) continue;
			sendMPI(&r.m_val[0], j * sizeof(R), i, 0);
			if(r.m_full) continue;
			sendMPI(&r.m_row[0], j * sizeof(uint), i, 0);
		}
		// receive r data from all neighbors
		Buffer< pair<uint,R> > rval(r.m_val.size());
		for(k=0; k<rval.size(); k++) rval[k] = pair<uint,R>((r.m_full ? k : r.m_row[k]), r.m_val[k]);
		Buffer< pair<uint,uint> > ext;
		for(i=0; i<ranks; i++) {
			if(i == irank) continue;
			recvMPI(&j, sizeof(uint), i, 0);
			if(j == 0) continue;
			const uint rvals = rval.size();
			const uint exts = ext.size();
			Buffer<R> ival(j);
			recvMPI(&ival[0], j * sizeof(R), i, 0);
			rval.resize(rvals + j);
			for(k=0; k<j; k++) rval[rvals + k] = pair<uint,R>(r.m_height + exts + k, ival[k]);
			ext.resize(exts + j);
			if(r.m_full) {
				for(k=0; k<j; k++) ext[exts + k] = pair<uint,uint>(i,k);
				continue;
			}
			Buffer<uint> iext(j);
			recvMPI(&iext[0], j * sizeof(uint), i, 0);
			for(k=0; k<j; k++) ext[exts + k] = pair<uint,uint>(i,iext[k]);
		}
		if(l.m_full) {
			Buffer< Buffer< pair<uint,T> > > val(l.m_height);
			for(i=0; i<val.size(); i++) {
				Buffer< pair<uint,T> > &ival = val[i];
				ival.resize(rval.size());
				for(j=0; j<rval.size(); j++) {
					ival[j] = pair<uint,T>(rval[j].first, l.m_val[i] * rval[j].second);
				}
			}
			return setFull(r.m_height, val, ext);
		}
		Buffer< pair<uint,Buffer< pair<uint,T> > > > val(l.m_row.size());
		for(i=0; i<val.size(); i++) {
			pair<uint,Buffer< pair<uint,T> > > &ival = val[i];
			ival.first = l.m_row[i];
			ival.second.resize(rval.size());
			for(j=0; j<rval.size(); j++) {
				ival.second[j] = pair<uint,T>(rval[j].first, l.m_val[i] * rval[j].second);
			}
		}
		return setSparse(r.m_height, l.m_height, val, ext);
	}
	template<typename L, typename R> Sparse &setScale(const Sparse<L> &l, const R &r) {
		if(this != (void*)&l) setShape(l);
		uint i;
		for(i=0; i<this->m_val.size(); i++) this->m_val[i] = l.m_val[i] * r;
		for(i=0; i<m_rval.size(); i++) m_rval[i] = l.m_rval[i] * r;
		return *this;
	}
	template<typename L, typename R> Sparse &setScale(const L &l, const Sparse<R> &r) {
		if(this != (void*)&r) setShape(r);
		uint i;
		for(i=0; i<this->m_val.size(); i++) this->m_val[i] = l * r.m_val[i];
		for(i=0; i<m_rval.size(); i++) m_rval[i] = l * r.m_rval[i];
		return *this;
	}

	// modifier functions
	template<typename R> Sparse &operator+=(const Sparse<R> &r) { return setPlus(*this, r); }
	template<typename R> Sparse &operator-=(const Sparse<R> &r) { return setMinus(*this, r); }
	template<typename R> Sparse &operator+=(const Diagonal<R>& r) { return setPlus(*this, r); }
	template<typename R> Sparse &operator-=(const Diagonal<R>& r) { return setMinus(*this, r); }
	template<typename R> Sparse &operator*=(const Sparse<R> &r) { return setTimes(*this, r); }
	template<typename R> Sparse &operator*=(const Diagonal<R> &r) { return setTimes(*this, r); } // multiply from right by a diagonal matrix
	template<typename R> Sparse &scale(const R &r) { return setScale(*this, r); } // multiply from right by a scalar
	Sparse &negate() { return setNegation(*this); }

	// trim functions
	Sparse &trimFull() { // convert sparse to full
		if(this->m_full) return trim(); // already full
		sparseToFull();
		return trim();
	}
	Sparse &trimSparse() { // convert full to sparse
		if(!this->m_full) return trim(); // already sparse
		this->m_full = false;
		this->m_row.resize(this->m_height);
		for(uint i=0; i<this->m_height; i++) this->m_row[i] = i;
		return trim();
	}
	Sparse &trimOptimal(const double limit) { // convert to full if(number of non-empty rows > limit * m_height), otherwise convert to sparse
		trimSparse();
		if(double(this->m_row.size()) > limit * this->m_height) sparseToFull();
		return *this;
	}
	Sparse &trim() { // remove all zero instances
		uint i, j, k, l;
		// trim data for MPI receive
		uint recvs = 0;
		uint rvals = 0;
		for(i=0,l=0; i<m_recv.size(); )	{
			const uint ranki = m_recv[i++];
			Buffer<bool> need(m_recv[i++]);
			m_recv[recvs++] = ranki;
			uint &needs = m_recv[recvs++];
			needs = 0;
			for(j=0; j<need.size(); j++) {
				const uint recvs0 = recvs++;
				for(k=m_recv[i++]; k>0; k--,i++,l++) {
					if(m_rval[l] == this->m_zero) continue;
					m_recv[recvs++] = m_recv[i];
					m_rval[rvals++] = m_rval[l];
				}
				m_recv[recvs0] = recvs - recvs0 - 1;
				need[j] = (m_recv[recvs0] > 0);
				if(need[j]) needs++;
				else recvs--; // remove communication of current term
			}
			sendMPI(&need[0], need.size() * sizeof(bool), ranki, 0);
			if(needs == 0) recvs -= 2; // remove communication with ranki
		}
		m_recv.resize(recvs);
		m_rval.resize(rvals);
		// trim data fo MPI send
		uint sends = 0;
		for(i=0; i<m_send.size(); ) {
			const uint ranki = m_send[i++];
			Buffer<bool> need(m_send[i++]);
			recvMPI(&need[0], need.size() * sizeof(bool), ranki, 0);
			m_send[sends++] = ranki;
			const uint sends0 = sends++;
			for(j=0; j<need.size(); j++,i++) {
				if(need[j]) m_send[sends++] = m_send[i];
			}
			m_send[sends0] = sends - sends0 - 1;
			if(m_send[sends0] == 0) sends -= 2; // remove communication with ranki
		}
		m_send.resize(sends);
		// trim columns
		for(i=0,j=0,k=0; i<m_col.size(); i++) {
			while(k<m_beg.size() && m_beg[k]<=i) m_beg[k++] = j;
			if(this->m_val[i] == this->m_zero) continue;
			if(j != i) {
				m_col[j] = m_col[i];
				this->m_val[j] = this->m_val[i];
			}
			j++;
		}
		while(k<m_beg.size()) m_beg[k++] = j;
		m_col.resize(j);
		this->m_val.resize(j);
		if(this->m_full) return *this;
		return toSparse();
	}

	// get functions
	const T &getValue(uint i, uint j) const {
		if(this->m_full) {
			if(i >= this->m_height) return this->m_zero;
		}
		else if(!this->searchIndex(i, this->m_row, 0, this->m_row.size())) return this->m_zero;
		if(!this->searchIndex(j, m_col, m_beg[i], (i+1<m_beg.size() ? m_beg[i+1] : m_col.size()))) return this->m_zero;
		return this->m_val[j];
	}

	// all variables are public
	uint m_width; // matrix width
	Buffer<uint> m_beg; // begin index for each row
	Buffer<uint> m_col; // column indices
	Buffer<uint> m_send; // coded data for MPI send
	Buffer<uint> m_recv; // coded data for MPI receive
	Buffer<T> m_rval; // multiplication terms for MPI receive

protected:
	Sparse &setEmpty() {
		Discrete<T>::setEmpty();
		m_width = 0;
		m_beg.clear();
		m_col.clear();
		m_send.clear();
		m_recv.clear();
		m_rval.clear();
		return *this;
	}
	Sparse &setCommunication(const Buffer< pair<uint,uint> > &ext, const Buffer< Buffer< pair<uint, T> > > &extv, const Buffer<uint> &srank, const Buffer<uint> &rrank) {
		uint i, j, k;
		// send information of required terms to each rank
		m_recv.clear();
		m_rval.clear();
		for(i=0; i<rrank.size(); i++) {
			uint links = 0;
			uint rvals = 0;
			for(j=0; j<ext.size(); j++) {
				if(ext[j].first != rrank[i]) continue;
				if(extv[j].empty()) continue;
				links++;
				rvals += extv[j].size();
			}
			sendMPI(&links, sizeof(uint), rrank[i], 0);
			if(links == 0) continue;
			Buffer<uint> link(links); // links to receive from rrank[i]
			Buffer<T> rval(rvals); // a buffer to combine into m_rval
			Buffer<uint> recv(2 + links + rvals); // a buffer to combine into m_recv
			links = 0;
			rvals = 0;
			uint recvs = 0;
			recv[recvs++] = rrank[i]; // rank to receive from
			recv[recvs++] = link.size(); // number of terms to receive
			for(j=0; j<ext.size(); j++) {
				if(ext[j].first != rrank[i]) continue;
				if(extv[j].empty()) continue;
				link[links++] = ext[j].second; // the term id in the neighbor rank
				recv[recvs++] = extv[j].size(); // number of copies taken from the term
				for(k=0; k<extv[j].size(); k++) {
					recv[recvs++] = extv[j][k].first; // a row to copy the term to
					rval[rvals++] = extv[j][k].second; // a multiplier for the term
				}
			}
			sendMPI(&link[0], link.size() * sizeof(uint), rrank[i], 1);
			m_recv.combine(recv);
			m_rval.combine(rval);
		}
		// receive information of required terms to be sent
		m_send.clear();
		for(i=0; i<srank.size(); i++) {
			uint sends;
			recvMPI(&sends, sizeof(uint), srank[i], 0);
			if(sends == 0) continue;
			Buffer<uint> send(2 + sends);
			send[0] = srank[i]; // rank to send to
			send[1] = sends; // number of terms to send
			recvMPI(&send[2], sends * sizeof(uint), srank[i], 1);
			m_send.combine(send);
		}
		return *this;
	}
	template<typename R> Sparse &setShape(const Sparse<R> &r) {
		Discrete<T>::setShape(r);
		m_width = r.m_width;
		m_beg = r.m_beg;
		m_col = r.m_col;
		m_send = r.m_send;
		m_recv = r.m_recv;
		m_rval.resize(r.m_rval.size());
		return *this;
	}
	template<typename R> pair<uint,R> getTerm(const pair<uint,uint> &ext, const R &val, const uint locs, map<pair<uint,uint>,uint> &exts) const {
		if(ext.first == getMPIrank()) return pair<uint,R>(ext.second, val);
		map<pair<uint,uint>,uint>::iterator it = exts.find(ext);
		if(it != exts.end()) return pair<uint,R>(it->second, val);
		const uint col = exts.size() + locs;
		exts.insert(make_pair(ext, col));
		return pair<uint,R>(col, val);
	}
	void sumTerm(const pair<uint,T> &val, Buffer< pair<uint,T> > &buf) {
		uint i0 = 0;
		uint i1 = buf.size();
		while(true) {
			if(i0 == i1) {
				buf.insert(val, i0);
				break;
			}
			const uint i = (i0 + i1) / 2;
			if(val.first < buf[i].first) i1 = i;
			else if(buf[i].first < val.first) i0 = i + 1;
			else {
				buf[i].second += val.second;
				break;
			}
		}
	}
	template<typename L, typename R> void mergeRow(const Sparse<L> &l, const uint li, const Diagonal<R> &r, const uint ri, Buffer<T> &val, Buffer<uint> &col, uint &cols, T func(const L &, const R &)) {
		uint lj = 0, l_end = 0;
		if(li < l.m_beg.size()) {
			lj = l.m_beg[li];
			l_end = (li+1<l.m_beg.size() ? l.m_beg[li+1] : l.m_col.size());
		}
		if(ri < r.m_val.size()) {
			const uint rcol = (r.m_full ? ri : r.m_row[ri]);
			while(lj < l_end && l.m_col[lj] < rcol) {
				val[cols] = func(l.m_val[lj], r.m_zero);
				col[cols] = l.m_col[lj];
				cols++; lj++;
			}
			if(lj < l_end && l.m_col[lj] == rcol) {
				val[cols] = func(l.m_val[lj], r.m_val[ri]);
				col[cols] = rcol;
				cols++; lj++;
			}
			else {
				val[cols] = func(l.m_zero, r.m_val[ri]);
				col[cols] = rcol;
				cols++;
			}
		}
		while(lj < l_end) {
			val[cols] = func(l.m_val[lj], r.m_zero);
			col[cols] = l.m_col[lj];
			cols++; lj++;
		}
	}
	template<typename L, typename R> void mergeRow(const Sparse<L> &l, const uint li, const Sparse<R> &r, const uint ri, Buffer<T> &val, Buffer<uint> &col, uint &cols, T func(const L &, const R &)) {
		uint lj = 0, rj = 0;
		uint l_end = 0, r_end = 0;
		if(li < l.m_beg.size()) {
			lj = l.m_beg[li];
			l_end = (li+1<l.m_beg.size() ? l.m_beg[li+1] : l.m_col.size());
		}
		if(ri < r.m_beg.size()) {
			rj = r.m_beg[ri];
			r_end = (ri+1<r.m_beg.size() ? r.m_beg[ri+1] : r.m_col.size());
		}
		while(lj < l_end && rj < r_end) {
			if(l.m_col[lj] < r.m_col[rj]) {
				val[cols] = func(l.m_val[lj], r.m_zero);
				col[cols] = l.m_col[lj];
				cols++; lj++;
			}
			else if(l.m_col[lj] > r.m_col[rj]) {
				val[cols] = func(l.m_zero, r.m_val[rj]);
				col[cols] = r.m_col[rj];
				cols++; rj++;
			}
			else {
				val[cols] = func(l.m_val[lj], r.m_val[rj]);
				col[cols] = l.m_col[lj];
				cols++; lj++; rj++;
			}
		}
		while(lj < l_end) {
			val[cols] = func(l.m_val[lj], r.m_zero);
			col[cols] = l.m_col[lj];
			cols++; lj++;
		}
		while(rj < r_end) {
			val[cols] = func(l.m_zero, r.m_val[rj]);
			col[cols] = r.m_col[rj];
			cols++; rj++;
		}
	}
	Sparse &toSparse() { // make sparse matrix with no empty rows
		uint i, j, k;
		if(this->m_full) {
			this->m_full = false;
			this->m_row.resize(this->m_height);
			for(i=0; i<this->m_height; i++) this->m_row[i] = i;
		}
		Buffer<uint> rrow(this->m_row.size(), uint(-1));
		for(i=0; i<m_recv.size(); i++) {
			for(j=m_recv[++i]; j>0; j--) {
				for(k=m_recv[++i]; k>0; k--) { ++i; rrow[m_recv[i]] = m_recv[i]; }
			}
		}
		for(i=0,j=0; i<this->m_row.size(); i++) {
			if((i+1<m_beg.size() ? m_beg[i+1] : m_col.size()) == m_beg[i] && rrow[i] != i) continue;
			if(i != j) {
				this->m_row[j] = this->m_row[i];
				m_beg[j] = m_beg[i];
				rrow[i] = j;
			}
			j++;
		}
		this->m_row.resize(j);
		m_beg.resize(j);
		for(i=0; i<m_recv.size(); i++) {
			for(j=m_recv[++i]; j>0; j--) {
				for(k=m_recv[++i]; k>0; k--) { ++i; m_recv[i] = rrow[m_recv[i]]; }
			}
		}
		return *this;
	}

	void sparseToFull() {
		this->m_full = true;
		uint i, j, k;
		for(i=0; i<m_recv.size(); )	{
			i++;
			for(j=m_recv[i++]; j>0; j--) {
				for(k=m_recv[i++]; k>0; k--,i++) m_recv[i] = this->m_row[m_recv[i]];
			}
		}
		Buffer<uint> beg(this->m_height);
		for(i=0,j=0; i<this->m_row.size(); i++) {
			while(j <= this->m_row[i]) beg[j++] = m_beg[i];
		}
		while(j < this->m_height) beg[j++] = m_col.size();
		m_beg.swap(beg);
		this->m_row.clear();
	}

};

// operators
template<typename L, typename R, typename O = decltype(declval<L &>() + declval<R &>())> Sparse<O> operator+(const Sparse<L> &l, const Sparse<R> &r) {
	Sparse<O> o(l.m_zero + r.m_zero);
	return o.setPlus(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() + declval<R &>())> Sparse<O> operator+(const Sparse<L> &l, const Diagonal<R> &r) {
	Sparse<O> o(l.m_zero + r.m_zero);
	return o.setPlus(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() + declval<R &>())> Sparse<O> operator+(const Diagonal<L> &l, const Sparse<R> &r) {
	Sparse<O> o(l.m_zero + r.m_zero);
	return o.setPlus(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() - declval<R &>())> Sparse<O> operator-(const Sparse<L> &l, const Sparse<R> &r) {
	Sparse<O> o(l.m_zero - r.m_zero);
	return o.setMinus(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() - declval<R &>())> Sparse<O> operator-(const Sparse<L> &l, const Diagonal<R> &r) {
	Sparse<O> o(l.m_zero + r.m_zero);
	return o.setMinus(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() - declval<R &>())> Sparse<O> operator-(const Diagonal<L> &l, const Sparse<R> &r) {
	Sparse<O> o(l.m_zero + r.m_zero);
	return o.setMinus(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Sparse<O> operator*(const Sparse<L> &l, const Sparse<R> &r) {
	Sparse<O> o(l.m_zero * r.m_zero);
	return o.setTimes(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Sparse<O> operator*(const Sparse<L> &l, const Diagonal<R> &r) {
	Sparse<O> o(l.m_zero * r.m_zero);
	return o.setTimes(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Sparse<O> operator*(const Diagonal<L> &l, const Sparse<R> &r) {
	Sparse<O> o(l.m_zero * r.m_zero);
	return o.setTimes(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Sparse<O> scaled(const Sparse<L> &l, const R &r) {
	Sparse<O> o(l.m_zero * r);
	return o.setScale(l, r);
}
template<typename L, typename R, typename O = decltype(declval<L &>() * declval<R &>())> Sparse<O> scaled(const L &l, const Sparse<R> &r) {
	Sparse<O> o(l * r.m_zero);
	return o.setScale(l, r);
}
template<typename R> Sparse<R> transpose(const Sparse<R>& r) {
	Sparse<R> o(r.m_zero);
	return o.setTranspose(r);
}
template<typename R> Sparse<R> operator-(const Sparse<R> &r) {
	Sparse<R> o(r.m_zero);
	return o.setNegation(r);
}


}

#endif //_SPARSE_HPP_INCLUDED_
