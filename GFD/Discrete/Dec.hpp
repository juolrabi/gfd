/**
 * class Dec implements tools of discrete exterior calculus.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _DEC_HPP_INCLUDED_
#define _DEC_HPP_INCLUDED_

#include "MeshIntegrator.hpp"
#include "../Mesh/PartMesh.hpp"
#include "../Types/MpiEasy.hpp"

namespace gfd
{

class Dec
{
public:
	Dec(const PartMesh &mesh);
	Dec(const PartMesh &mesh, const uint lowdim, const uint highdim);
	virtual ~Dec() { }

	void setLowDimension(const uint lowdim);
	void setHighDimension(const uint highdim);

	Diagonal<bool> &integrateFlags(const FormGrade grade, const UintSet &flag, Diagonal<bool> &result) const;

	Sparse<sign> &integrateDerivative(const FormGrade grade, Sparse<sign> &result) const;
	Sparse<sign> &integrateDerivative(const FormGrade grade, const UintSet &flag, Sparse<sign> &result) const;
	Sparse<double> &integrateCurvatureDerivative(SymMatrix4 curv(const Vector4 &), const FormGrade grade, const UintSet &flag, Sparse<double> &result) const;

	template<typename T> Column<T> &integrateZeroForm(const FormGrade grade, Column<T> &result) const {
		MeshIntegrator intg(m_mesh, grade, 0, m_lowdim, m_highdim);
		const uint locs = intg.getLocals();
		initResult(locs, result);
		return result;
	}
	template<typename T> Column<T> &integrateForm(T func(const Buffer<double> &), const int num, const FormGrade grade, Column<T> &result) const {
		MeshIntegrator intg(m_mesh, grade, num, m_lowdim, m_highdim);
		const uint locs = intg.getLocals();
		initResult(locs, result);
		// integrate over mesh elements
		if(FormGradeIsPrim(grade)) {
			for(uint i=0; i<locs; i++) integrateElement(intg, i, func, result.m_val[i]);
			return result;
		}
		// dual integration from base up
		const Buffer< pair<uint,uint> > &ext = intg.getExternals();
		Buffer<T> extval(ext.size(), result.m_zero);
		const uint baselocs = intg.getBaseLocals();
		for(uint i=0; i<baselocs; i++) {
			integrateBaseElement(intg, i, func, result.m_val, extval);
		}
		shareExternals(ext, extval, result.m_val);
		return result;
	}
	template<typename T> Column<T> &integrateForm(T func(const Buffer<double> &), const int num, const FormGrade grade, const UintSet &flag, Column<T> &result) const {
		MeshIntegrator intg(m_mesh, grade, num, m_lowdim, m_highdim);
		const uint locs = intg.getLocals();
		initResult(locs, result);
		// integrate over mesh elements
		for(uint i=0; i<locs; i++) {
			if(flag.includes(intg.getFlag(i))) integrateElement(intg, i, func, result.m_val[i]);
		}
		if(FormGradeIsPrim(grade)) return result;
		// external integration
		const Buffer< pair<uint,uint> > &ext = intg.getExternals();
		Buffer<T> extval(ext.size(), result.m_zero);
		for(uint i=0; i<extval.size(); i++) {
			const uint ii = locs + i;
			if(flag.includes(intg.getFlag(ii))) integrateElement(intg, ii, func, extval[i]);
		}
		shareExternals(ext, extval, result.m_val);
		return result;
	}
	template<typename T> Column<T> &integrateForm(const UintSet &baseflag, T func(const Buffer<double> &), const int num, const FormGrade grade, Column<T> &result) const {
		MeshIntegrator intg(m_mesh, grade, num, m_lowdim, m_highdim);
		const uint locs = intg.getLocals();
		initResult(locs, result);
		// integrate over mesh elements
		if(FormGradeIsPrim(grade)) { // primal integration from base up
			const uint basesize = intg.getBaseSize();
			for(uint i=0; i<basesize; i++) {
				if(!baseflag.includes(intg.getBaseFlag(i))) continue;
				Buffer< pair<uint,Quadrature> > q;
				intg.getBaseQuadrature(i, q);
				for(uint j=0; j<q.size(); j++) {
					if(q[j].first >= locs) continue;
					q[j].second.integrate(func, result.m_val[q[j].first]);
				}
			}
			return result;
		}
		// dual integration from base up
		const Buffer< pair<uint,uint> > &ext = intg.getExternals();
		Buffer<T> extval(ext.size(), result.m_zero);
		const uint baselocs = intg.getBaseLocals();
		for(uint i=0; i<baselocs; i++) {
			if(baseflag.includes(intg.getBaseFlag(i))) integrateBaseElement(intg, i, func, result.m_val, extval);
		}
		shareExternals(ext, extval, result.m_val);
		return result;
	}

	template<typename T> Diagonal<T> &integrateZeroHodge(const FormGrade grade, Diagonal<T> &result) const {
		MeshIntegrator intg(m_mesh, grade, 0, m_lowdim, m_highdim);
		const uint locs = intg.getLocals();
		initResult(locs, result);
	}
	template<typename T> Diagonal<T> &integrateHodge(T func(const Buffer<double> &), const int num, const FormGrade grade, Diagonal<T> &result) const {
		const bool isprim = FormGradeIsPrim(grade);
		MeshIntegrator intg0(m_mesh, (isprim ? grade : FormGradeDual(grade)), num, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, (!isprim ? grade : FormGradeDual(grade)), num, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		// integrate dual quadrature
		const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
		Buffer<Quadrature> q1(locs + ext.size(), intg1.getEmptyQuadrature()); 
		const uint baselocs = intg1.getBaseLocals();
		for(uint i=0; i<baselocs; i++) {
			Buffer< pair<uint,Quadrature> > q;
			intg1.getBaseQuadrature(i, q);
			for(uint j=0; j<q.size(); j++) q1[q[j].first].combine(q[j].second);
		}
		combineExternals(ext, q1);
		for(uint i=0; i<locs; i++) {
			const Vector4 p = intg0.getPosition(i);
			Quadrature q0;
			intg0.getQuadrature(i, q0);
			if(isprim) q0.invert(p).integrateProduct(q1[i], func, result.m_val[i]);
			else q0.integrateProduct(q1[i].invert(p), func, result.m_val[i]);
			q1[i].clear();
		}
		return result;
	}
	template<typename T> Diagonal<T> &integrateHodge(T func(const Buffer<double> &), const int num, const FormGrade grade, const UintSet &flag, Diagonal<T> &result) const {
		MeshIntegrator intg0(m_mesh, grade, num, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, FormGradeDual(grade), num, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		if(FormGradeIsPrim(grade)) { // integrate dual quadrature divided by primal
			const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
			Buffer<T> extval(ext.size(), result.m_zero);
			const uint size = locs + ext.size();
			for(uint i=0; i<size; i++) {
				if(!flag.includes(intg0.getFlag(i))) continue;
				Quadrature q0;
				intg0.getQuadrature(i, q0);
				q0.invert(intg0.getPosition(i));
				Quadrature q1;
				intg1.getQuadrature(i, q1);
				T &val = (i < locs ? result.m_val[i] : extval[i - locs]);
				q0.integrateProduct(q1, func, val);
			}
			shareExternals(ext, extval, result.m_val);
			return result;
		}
		// integrate primal quadrature divided by dual
		const Buffer< pair<uint,uint> > &ext = intg0.getExternals();
		Buffer<Quadrature> q0(locs + ext.size(), intg0.getEmptyQuadrature()); 
		for(uint i=0; i<ext.size(); i++) {
			const uint ii = locs + i;
			if(flag.includes(intg0.getFlag(ii))) intg0.getQuadrature(ii, q0[ii]);
		}
		combineExternals(ext, q0);
		for(uint i=0; i<locs; i++) {
			if(!flag.includes(intg0.getFlag(i))) continue;
			Quadrature q;
			intg0.getQuadrature(i, q);
			q0[i].combine(q);
			q0[i].invert(intg0.getPosition(i));
			Quadrature q1;
			intg1.getQuadrature(i, q1);
			q0[i].integrateProduct(q1, func, result.m_val[i]);
			q0[i].clear();
		}
		return result;
	}
	template<typename T> Diagonal<T> &integrateHodge(const UintSet &baseflag, T func(const Buffer<double> &), const int num, const FormGrade grade, Diagonal<T> &result) const {
		MeshIntegrator intg0(m_mesh, grade, num, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, FormGradeDual(grade), num, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		if(FormGradeIsPrim(grade)) { // integrate dual quadrature divided by primal
			const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
			Buffer<Quadrature> q1(locs + ext.size(), intg1.getEmptyQuadrature()); 
			const uint baselocs = intg1.getBaseLocals();
			for(uint i=0; i<baselocs; i++) {
				if(!baseflag.includes(intg1.getBaseFlag(i))) continue;
				Buffer< pair<uint,Quadrature> > q;
				intg1.getBaseQuadrature(i, q);
				for(uint j=0; j<q.size(); j++) q1[q[j].first].combine(q[j].second);
			}
			combineExternals(ext, q1);
			for(uint i=0; i<locs; i++) {
				if(q1[i].empty()) continue;
				const Vector4 p = intg0.getPosition(i);
				Quadrature q0;
				intg0.getQuadrature(i, q0);
				q0.invert(p);
				q0.integrateProduct(q1[i], func, result.m_val[i]);
				q1[i].clear();
			}
			return result;
		}
		// integrate primal quadrature divided by dual
		const Buffer< pair<uint,uint> > &ext = intg0.getExternals();
		Buffer<Quadrature> q0(locs + ext.size(), intg0.getEmptyQuadrature()); 
		for(uint i=0; i<ext.size(); i++) { intg0.getQuadrature(locs + i, q0[locs + i]); }
		combineExternals(ext, q0);
		Buffer<bool> done(locs, false);
		const uint basesize = intg1.getBaseSize();
		for(uint i=0; i<basesize; i++) {
			if(!baseflag.includes(intg1.getBaseFlag(i))) continue;
			Buffer< pair<uint,Quadrature> > q1;
			intg1.getBaseQuadrature(i, q1);
			for(uint j=0; j<q1.size(); j++) {
				const uint jj = q1[j].first;
				if(jj >= locs) continue;
				if(!done[jj]) { // create quadrature q0[jj]
					Quadrature q;
					intg0.getQuadrature(jj, q);
					q0[jj].combine(q);
					q0[jj].invert(intg0.getPosition(jj));
					done[jj] = true;
				}
				q0[jj].integrateProduct(q1[j].second, func, result.m_val[jj]);
			}
		}
		return result;
	}

	template<typename T> Diagonal<T> &integrateWedge(T func(const Buffer<double> &), const int num, const FormGrade grade, Diagonal<T> &result) const {
		const bool isprim = FormGradeIsPrim(grade);
		MeshIntegrator intg0(m_mesh, (isprim ? grade : FormGradeDual(grade)), num, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, (!isprim ? grade : FormGradeDual(grade)), num, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
		Buffer<Quadrature> q1(locs + ext.size(), intg1.getEmptyQuadrature()); 
		const uint baselocs = intg1.getBaseLocals();
		for(uint i=0; i<baselocs; i++) {
			Buffer< pair<uint,Quadrature> > q;
			intg1.getBaseQuadrature(i, q);
			for(uint j=0; j<q.size(); j++) q1[q[j].first].combine(q[j].second);
		}
		combineExternals(ext, q1);
		for(uint i=0; i<locs; i++) {
			const Vector4 p = intg0.getPosition(i);
			Quadrature q0;
			intg0.getQuadrature(i, q0);
			q0.relocate(-p).integrateProduct(q1[i], func, result.m_val[i]);
			q1[i].clear();
		}
		return result;
	}
	template<typename T> Diagonal<T> &integrateWedge(T func(const Buffer<double> &), const int num, const FormGrade grade, const UintSet &flag, Diagonal<T> &result) const {
		const bool isprim = FormGradeIsPrim(grade);
		MeshIntegrator intg0(m_mesh, (isprim ? grade : FormGradeDual(grade)), num, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, (!isprim ? grade : FormGradeDual(grade)), num, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
		Buffer<T> extval(ext.size(), result.m_zero);
		const uint size = locs + ext.size();
		for(uint i=0; i<size; i++) {
			if(!flag.includes(intg0.getFlag(i))) continue;
			Quadrature q0;
			intg0.getQuadrature(i, q0);
			q0.relocate(-intg0.getPosition(i));
			Quadrature q1;
			intg1.getQuadrature(i, q1);
			T &val = (i < locs ? result.m_val[i] : extval[i - locs]);
			q0.integrateProduct(q1, func, val);
		}
		shareExternals(ext, extval, result.m_val);
		return result;
	}
	template<typename T> Diagonal<T> &integrateWedge(const UintSet &baseflag, T func(const Buffer<double> &), const int num, const FormGrade grade, Diagonal<T> &result) const {
		MeshIntegrator intg0(m_mesh, grade, num, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, FormGradeDual(grade), num, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		if(FormGradeIsPrim(grade)) { // integrate dual quadrature divided by primal
			const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
			Buffer<Quadrature> q1(locs + ext.size(), intg1.getEmptyQuadrature()); 
			const uint baselocs = intg1.getBaseLocals();
			for(uint i=0; i<baselocs; i++) {
				if(!baseflag.includes(intg1.getBaseFlag(i))) continue;
				Buffer< pair<uint,Quadrature> > q;
				intg1.getBaseQuadrature(i, q);
				for(uint j=0; j<q.size(); j++) q1[q[j].first].combine(q[j].second);
			}
			combineExternals(ext, q1);
			for(uint i=0; i<locs; i++) {
				if(q1[i].empty()) continue;
				const Vector4 p = intg0.getPosition(i);
				Quadrature q0;
				intg0.getQuadrature(i, q0);
				q0.relocate(-p);
				q0.integrateProduct(q1[i], func, result.m_val[i]);
				q1[i].clear();
			}
			return result;
		}
		// integrate primal quadrature divided by dual
		const Buffer< pair<uint,uint> > &ext = intg0.getExternals();
		Buffer<T> extval(ext.size(), result.m_zero);
		const uint basesize = intg1.getBaseSize();
		for(uint i=0; i<basesize; i++) {
			if(!baseflag.includes(intg1.getBaseFlag(i))) continue;
			Buffer< pair<uint,Quadrature> > q1;
			intg1.getBaseQuadrature(i, q1);
			for(uint j=0; j<q1.size(); j++) {
				const uint jj = q1[j].first;
				Quadrature q;
				intg0.getQuadrature(jj, q);
				q.relocate(-intg0.getPosition(jj));
				T &val = (jj < locs ? result.m_val[jj] : extval[jj - locs]);
				q.integrateProduct(q1[j].second, func, val);
			}
		}
		shareExternals(ext, extval, result.m_val);
		return result;
	}
	template<typename T> void getFullBuffer(const FormGrade grade, const Column<T> &val, Buffer<T> &result) {
		MeshIntegrator intg(m_mesh, grade, 0, m_lowdim, m_highdim);
		const Buffer< pair<uint,uint> > &ext = intg.getExternals();
		result = val.getBuffer();
		result.resize(val.m_height + ext.size());
		const Buffer< pair<uint,uint> > mext = getMyExternals(ext);
		const uint rank = getMPIrank();
		for(uint i=0; i<mext.size(); i++) {
			sendMPI(&result[mext[i].second], sizeof(T), mext[i].first, 0);
		}
		for(uint i=0; i<ext.size(); i++) {
			if(ext[i].first == rank) result[val.m_height + i] = result[ext[i].second];
			else recvMPI(&result[val.m_height + i], sizeof(T), ext[i].first, 0);
		}
	}

protected:
	const PartMesh &m_mesh;
	uint m_lowdim;
	uint m_highdim;

	template<typename C> void initResult(const uint locs, C &result) const {
		if(!result.m_full || result.m_height != locs) result.setFullOfZeros(locs);
	}
	void combineExternals(const Buffer< pair<uint,uint> > &ext, Buffer<Quadrature> &q) const;
	template<typename T> void shareExternals(const Buffer< pair<uint,uint> > &ext, const Buffer<T> &extval, Buffer<T> &locval) const {
		const Buffer< pair<uint,uint> > mext = getMyExternals(ext);
		const uint rank = getMPIrank();
		for(uint i=0; i<ext.size(); i++) {
			if(ext[i].first == rank) locval[ext[i].second] += extval[i];
			else sendMPI(&extval[i], sizeof(T), ext[i].first, 0);
		}
		T val;
		for(uint i=0; i<mext.size(); i++) {
			recvMPI(&val, sizeof(T), mext[i].first, 0);
			locval[mext[i].second] += val;
		}
	}
	template<typename T> void integrateElement(const MeshIntegrator &intg, const uint i, T func(const Buffer<double> &), T &val) const {
		Quadrature q;
		intg.getQuadrature(i, q);
		q.integrate(func, val);
	}
	template<typename T> void integrateBaseElement(const MeshIntegrator &intg, const uint i, T func(const Buffer<double> &), Buffer<T> &locval, Buffer<T> &extval) const {
		const uint locs = locval.size();
		Buffer< pair<uint,Quadrature> > q;
		intg.getBaseQuadrature(i, q);
		for(uint j=0; j<q.size(); j++) {
			const uint jj = q[j].first;
			T &val = (jj < locs ? locval[jj] : extval[jj - locs]);
			q[j].second.integrate(func, val);
		}
	}
	Buffer< pair<uint,uint> > getMyExternals(const Buffer< pair<uint,uint> > &ext) const;

};

}

#endif //_DEC_HPP_INCLUDED_
