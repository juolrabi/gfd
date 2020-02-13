/**
 * Split.hpp implements split classes of Column, Diagonal, and Sparse.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _SPLIT_HPP_INCLUDED_
#define _SPLIT_HPP_INCLUDED_

#include "Column.hpp"

namespace gfd
{

class Split
{
protected:
	uint getNumberOfTerms(const Buffer<uint> &categ) {
		uint terms = 0;
		for(uint i=0; i<categ.size(); i++) {
			if(categ[i] >= terms) terms = categ[i] + 1;
		}
		return terms;
	}
	Diagonal<bool> getSplitter(const Buffer<uint> &cat, const uint term) {
		uint i;
		uint terms = 0;
		for(i=0; i<cat.size(); i++) {
			if(cat[i] == term) terms++;
		}
		if(2 * terms >= cat.size()) {
			Buffer<bool> choose(cat.size(), false);
			for(i=0; i<cat.size(); i++) {
				if(cat[i] == term) choose[i] = true;
			}
			return Diagonal<bool>(choose, false);
		}
		Buffer< pair<uint,bool> > choose(terms);
		terms = 0;
		for(i=0; i<cat.size(); i++) {
			if(cat[i] == term) choose[terms++] = pair<uint,bool>(i, true);
		}
		return Diagonal<bool>(cat.size(), choose, false);
	}
	Diagonal<bool> getSplitterUntil(const Buffer<uint> &cat, const uint term) {
		uint i;
		uint terms = 0;
		for(i=0; i<cat.size(); i++) {
			if(cat[i] <= term) terms++;
		}
		if(2 * terms >= cat.size()) {
			Buffer<bool> choose(cat.size(), false);
			for(i=0; i<cat.size(); i++) {
				if(cat[i] <= term) choose[i] = true;
			}
			return Diagonal<bool>(choose, false);
		}
		Buffer< pair<uint,bool> > choose(terms);
		terms = 0;
		for(i=0; i<cat.size(); i++) {
			if(cat[i] <= term) choose[terms++] = pair<uint,bool>(i, true);
		}
		return Diagonal<bool>(cat.size(), choose, false);
	}
	void getCategories(Buffer<uint> &cat, const Buffer<uint> &row, const uint term) {
		for(uint i=0; i<row.size(); i++) cat[row[i]] = term;
	}

};

template<typename T>
class ColumnSplit : public Split
{
public:
	ColumnSplit() { }
	template<typename R> ColumnSplit(const Column<R> &term) { init(term); }
	template<typename R> ColumnSplit(const Buffer<uint> &lcat, const Column<R> &term) { init(lcat, term); }
	template<typename R> void init(const Column<R> &term) {
		m_term.resize(1);
		m_term[0].setCopy(term);
	}
	template<typename R> void init(const Buffer<uint> &lcat, const Column<R> &term) {
		m_term.resize(Split::getNumberOfTerms(lcat));
		for(uint i=0; i<m_term.size(); i++) {
			m_term[i].m_zero = term.m_zero;
			m_term[i].setTimes(Split::getSplitter(lcat, i), term);
		}
	}
	Buffer<uint> getCategories() const {
		Buffer<uint> cat(m_term[0].m_height, 0);
		for(uint i=1; i<m_term.size(); i++) {
			if(m_term[i].m_full) {
				cat.fill(i);
				continue;
			}
			const Buffer<uint> &row = m_term[i].m_row;
			for(uint j=0; j<row.size(); j++) cat[row[j]] = i;
		}
		return cat;
	}
	Buffer< Column<T> > m_term;
};

template<typename T>
class DiagonalSplit : public Split
{
public:
	DiagonalSplit() { }
	template<typename R> DiagonalSplit(const Diagonal<R> &term) { init(term); }
	template<typename R> DiagonalSplit(const Buffer<uint> &lcat, const Diagonal<R> &term) { init(lcat, term); }
	template<typename R> DiagonalSplit(const Diagonal<R> &term, const Buffer<uint> &rcat) { init(term, rcat); }
	template<typename R> DiagonalSplit(const Buffer<uint> &lcat, const Diagonal<R> &term, const Buffer<uint> &rcat) { init(lcat, term, rcat); }
	template<typename R> void init(const Diagonal<R> &term) {
		m_term.resize(1);
		m_term[0].setCopy(term);
	}
	template<typename R> void init(const Buffer<uint> &lcat, const Diagonal<R> &term) {
		m_term.resize(Split::getNumberOfTerms(lcat));
		for(uint i=0; i<m_term.size(); i++) {
			m_term[i].m_zero = term.m_zero;
			m_term[i].setTimes(Split::getSplitter(lcat, i), term);
		}
	}
	template<typename R> void init(const Diagonal<R> &term, const Buffer<uint> &rcat) {
		m_term.resize(Split::getNumberOfTerms(rcat));
		for(uint i=0; i<m_term.size(); i++) {
			m_term[i].m_zero = term.m_zero;
			m_term[i].setTimes(term, Split::getSplitter(rcat, i));
		}
	}
	template<typename R> void init(const Buffer<uint> &lcat, const Diagonal<R> &term, const Buffer<uint> &rcat) {
		const uint lterms = Split::getNumberOfTerms(lcat);
		const uint rterms = Split::getNumberOfTerms(rcat);
		m_term.resize(lterms > rterms ? lterms : rterms);
		for(uint i=0; i<m_term.size(); i++) {
			m_term[i].m_zero = term.m_zero;
			m_term[i].setTimes(Split::getSplitter(lcat, i) * term, Split::getSplitterUntil(rcat, i));
			if(i > 0) m_term[i] += Split::getSplitterUntil(lcat, i-1) * term * Split::getSplitter(rcat, i);
		}
	}
	Buffer<uint> getCategories() const {
		Buffer<uint> cat(m_term[0].m_height, 0);
		for(uint i=1; i<m_term.size(); i++) {
			if(m_term[i].m_full) {
				cat.fill(i);
				continue;
			}
			const Buffer<uint> &row = m_term[i].m_row;
			for(uint j=0; j<row.size(); j++) cat[row[j]] = i;
		}
		return cat;
	}
	Buffer< Diagonal<T> > m_term;
};

template<typename T>
class SparseSplit : public Split
{
public:
	SparseSplit() { }
	template<typename R> SparseSplit(const Sparse<R> &term) { init(term); }
	template<typename R> SparseSplit(const Buffer<uint> &lcat, const Sparse<R> &term) { init(lcat, term); }
	template<typename R> SparseSplit(const Sparse<R> &term, const Buffer<uint> &rcat) { init(term, rcat); }
	template<typename R> SparseSplit(const Buffer<uint> &lcat, const Sparse<R> &term, const Buffer<uint> &rcat) { init(lcat, term, rcat); }
	template<typename R> void init(const Sparse<R> &term) {
		m_term.resize(1);
		m_term[0].setCopy(term);
	}
	template<typename R> void init(const Buffer<uint> &lcat, const Sparse<R> &term) {
		m_term.resize(Split::getNumberOfTerms(lcat));
	cout << "numberoftemrs" << m_term.size() << endl;
		for(uint i=0; i<m_term.size(); i++) {
			m_term[i].m_zero = term.m_zero;
			m_term[i].setTimes(Split::getSplitter(lcat, i), term);
		}
	}
	template<typename R> void init(const Sparse<R> &term, const Buffer<uint> &rcat) {
		m_term.resize(Split::getNumberOfTerms(rcat));
		for(uint i=0; i<m_term.size(); i++) {
			m_term[i].m_zero = term.m_zero;
			m_term[i].setTimes(term, Split::getSplitter(rcat, i));
		}
	}
	template<typename R> void init(const Buffer<uint> &lcat, const Sparse<R> &term, const Buffer<uint> &rcat) {
		const uint lterms = Split::getNumberOfTerms(lcat);
		const uint rterms = Split::getNumberOfTerms(rcat);
		m_term.resize(lterms > rterms ? lterms : rterms);
		for(uint i=0; i<m_term.size(); i++) {
			m_term[i].m_zero = term.m_zero;
			m_term[i].setTimes(Split::getSplitter(lcat, i) * term, Split::getSplitterUntil(rcat, i));
			if(i > 0) m_term[i] += Split::getSplitterUntil(lcat, i-1) * term * Split::getSplitter(rcat, i);
		}
	}
	Buffer<uint> getCategories() const {
		Buffer<uint> cat(m_term[0].m_height, 0);
		for(uint i=1; i<m_term.size(); i++) {
			if(m_term[i].m_full) {
				cat.fill(i);
				continue;
			}
			const Buffer<uint> &row = m_term[i].m_row;
			for(uint j=0; j<row.size(); j++) cat[row[j]] = i;
		}
		return cat;
	}
	Buffer< Sparse<T> > m_term;
};

}

#endif //_SPLIT_HPP_INCLUDED_
