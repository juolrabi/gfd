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

	Sparse<sign> &integrateDerivative(const FormGrade grade, Sparse<sign> &result) const;
	Sparse<sign> &integrateDerivative(const FormGrade grade, const UintSet &flag, Sparse<sign> &result) const;
	Sparse<double> &integrateCurvatureDerivative(SymMatrix4 curv(const Vector4 &), const FormGrade grade, const UintSet &flag, Sparse<double> &result) const;

	template<typename T> Column<T> &integrateZeroForm(const FormGrade grade, Column<T> &result) const {
		MeshIntegrator intg(m_mesh, grade, 0, m_lowdim, m_highdim);
		const uint locs = intg.getLocals();
		initResult(locs, result);
	}
	template<typename T> Column<T> &integrateConstantForm(const T *f, const FormGrade grade, Column<T> &result) const {
		MeshIntegrator intg(m_mesh, grade, 0, m_lowdim, m_highdim);
		const uint locs = intg.getLocals();
		initResult(locs, result);
		// integrate over mesh elements
		if(FormGradeIsPrim(grade)) {
			for(uint i=0; i<locs; i++) integrateConstantElement(intg, i, f, result.m_val[i]);
			return result;
		}
		// dual integration from base up
		const Buffer< pair<uint,uint> > &ext = intg.getExternals();
		Buffer<T> extval(ext.size(), result.m_zero);
		const uint baselocs = intg.getBaseLocals();
		for(uint i=0; i<baselocs; i++) {
			integrateConstantBaseElement(intg, i, f, result.m_val, extval);
		}
		shareExternals(ext, extval, result.m_val);
		return result;
	}
	template<typename T> Column<T> &integrateConstantForm(const T *f, const FormGrade grade, const UintSet &flag, Column<T> &result) const {
		MeshIntegrator intg(m_mesh, grade, 0, m_lowdim, m_highdim);
		const uint locs = intg.getLocals();
		initResult(locs, result);
		// integrate over mesh elements
		for(uint i=0; i<locs; i++) {
			if(flag.includes(intg.getFlag(i))) integrateConstantElement(intg, i, f, result.m_val[i]);
		}
		if(FormGradeIsPrim(grade)) return result;
		// external integration
		const Buffer< pair<uint,uint> > &ext = intg.getExternals();
		Buffer<T> extval(ext.size(), result.m_zero);
		for(uint i=0; i<extval.size(); i++) {
			const uint ii = locs + i;
			if(flag.includes(intg.getFlag(ii))) integrateConstantElement(intg, ii, f, extval[i]);
		}
		shareExternals(ext, extval, result.m_val);
		return result;
	}
	template<typename T> Column<T> &integrateConstantForm(const UintSet &baseflag, const T *f, const FormGrade grade, Column<T> &result) const {
		MeshIntegrator intg(m_mesh, grade, 0, m_lowdim, m_highdim);
		const uint locs = intg.getLocals();
		initResult(locs, result);
		// integrate over mesh elements
		if(FormGradeIsPrim(grade)) { // primal integration from base up
			const uint basesize = intg.getBaseSize();
			for(uint i=0; i<basesize; i++) {
				Buffer< pair<uint, Buffer<double> > > q;
				intg.getBaseVector(i, q);
				for(uint j=0; j<q.size(); j++) {
					T &val = result.m_val[q[j].first];
					const Buffer<double> &qj = q[j].second;
					for(uint k=0; k<qj.size(); k++) val += qj[k] * f[k];
				}
			}
			return result;
		}
		// dual integration from base up
		const Buffer< pair<uint,uint> > &ext = intg.getExternals();
		Buffer<T> extval(ext.size(), result.m_zero);
		const uint baselocs = intg.getBaseLocals();
		for(uint i=0; i<baselocs; i++) {
			if(baseflag.includes(intg.getBaseFlag(i))) integrateConstantBaseElement(intg, i, f, result.m_val, extval);
		}
		shareExternals(ext, extval, result.m_val);
		return result;
	}
	template<typename T> Column<T> &integrateForm(void func(const Vector4 &, T *), const uint num, const FormGrade grade, Column<T> &result) const {
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
	template<typename T> Column<T> &integrateForm(void func(const Vector4 &, T *), const uint num, const FormGrade grade, const UintSet &flag, Column<T> &result) const {
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
	template<typename T> Column<T> &integrateForm(const UintSet &baseflag, void func(const Vector4 &, T *), const uint num, const FormGrade grade, Column<T> &result) const {
		MeshIntegrator intg(m_mesh, grade, num, m_lowdim, m_highdim);
		const uint locs = intg.getLocals();
		initResult(locs, result);
		// integrate over mesh elements
		if(FormGradeIsPrim(grade)) { // primal integration from base up
			const uint basesize = intg.getBaseSize();
			for(uint i=0; i<basesize; i++) {
				if(!baseflag.includes(intg.getBaseFlag(i))) continue;
				Buffer< pair<uint, Buffer<double> > > q;
				intg.getBaseQuadrature(i, q);
				for(uint j=0; j<q.size(); j++) {
					intg.integrateQuadrature(q[j].second, func, result.m_val[q[j].first]);
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
	template<typename T> Diagonal<T> &integrateUnitHodge(const FormGrade grade, Diagonal<T> &result) const {
		const bool isprim = FormGradeIsPrim(grade);
		MeshIntegrator intg0(m_mesh, (isprim ? grade : FormGradeDual(grade)), 0, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, (!isprim ? grade : FormGradeDual(grade)), 0, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		// integrate dual quadrature
		const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
		Buffer< Buffer<double> > q1(locs + ext.size()); 
		const uint baselocs = intg1.getBaseLocals();
		for(uint i=0; i<baselocs; i++) {
			Buffer< pair<uint, Buffer<double> > > q;
			intg1.getBaseVector(i, q);
			for(uint j=0; j<q.size(); j++) q1[q[j].first].mergesum(q[j].second);
		}
		mergesumExternals(ext, q1);
		for(uint i=0; i<locs; i++) {
			Buffer<double> q0;
			intg0.getVector(i, q0);
			if(isprim) intg0.integrateUnitProductVector(intg0.invertVector(q0), q1[i], result.m_val[i]);
			else intg1.integrateUnitProductVector(q0, intg1.invertVector(q1[i]), result.m_val[i]);
			q1[i].clear();
		}
		return result;
	}
	template<typename T> Diagonal<T> &integrateUnitHodge(const FormGrade grade, const UintSet &flag, Diagonal<T> &result) const {
		MeshIntegrator intg0(m_mesh, grade, 0, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, FormGradeDual(grade), 0, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		if(FormGradeIsPrim(grade)) { // integrate dual quadrature divided by primal
			const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
			Buffer<T> extval(ext.size(), result.m_zero);
			const uint size = locs + ext.size();
			for(uint i=0; i<size; i++) {
				if(!flag.includes(intg0.getFlag(i))) continue;
				Buffer<double> q0;
				intg0.getVector(i, q0);
				intg0.invertVector(q0);
				Buffer<double> q1;
				intg1.getVector(i, q1);
				T &val = (i < locs ? result.m_val[i] : extval[i - locs]);
				intg0.integrateUnitProductVector(q0, q1, val);
			}
			shareExternals(ext, extval, result.m_val);
			return result;
		}
		// integrate primal quadrature divided by dual
		const Buffer< pair<uint,uint> > &ext = intg0.getExternals();
		Buffer< Buffer<double> > q0(locs + ext.size()); 
		for(uint i=0; i<ext.size(); i++) {
			const uint ii = locs + i;
			if(flag.includes(intg0.getFlag(ii))) intg0.getVector(ii, q0[ii]);
		}
		mergesumExternals(ext, q0);
		for(uint i=0; i<locs; i++) {
			if(!flag.includes(intg0.getFlag(i))) continue;
			Buffer<double> q;
			q0[i].mergesum(intg0.getVector(i, q));
			intg0.invertVector(q0[i]);
			Buffer<double> q1;
			intg1.getVector(i, q1);
			intg0.integrateUnitProductVector(q0[i], q1, result.m_val[i]);
			q0[i].clear();
		}
		return result;
	}
	template<typename T> Diagonal<T> &integrateUnitHodge(const UintSet &baseflag, const FormGrade grade, Diagonal<T> &result) const {
		MeshIntegrator intg0(m_mesh, grade, 0, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, FormGradeDual(grade), 0, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		if(FormGradeIsPrim(grade)) { // integrate dual quadrature divided by primal
			const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
			Buffer< Buffer<double> > q1(locs + ext.size()); 
			const uint baselocs = intg1.getBaseLocals();
			for(uint i=0; i<baselocs; i++) {
				if(!baseflag.includes(intg1.getBaseFlag(i))) continue;
				Buffer< pair<uint, Buffer<double> > > q;
				intg1.getBaseVector(i, q);
				for(uint j=0; j<q.size(); j++) q1[q[j].first].mergesum(q[j].second);
			}
			mergesumExternals(ext, q1);
			for(uint i=0; i<locs; i++) {
				if(q1[i].empty()) continue;
				Buffer<double> q0;
				intg0.getVector(i, q0);
				intg0.invertVector(q0);
				intg0.integrateUnitProductVector(q0, q1[i], result.m_val[i]);
				q1[i].clear();
			}
			return result;
		}
		// integrate primal quadrature divided by dual
		const Buffer< pair<uint,uint> > &ext = intg0.getExternals();
		Buffer< Buffer<double> > q0(locs + ext.size()); 
		for(uint i=0; i<ext.size(); i++) { intg0.getVector(locs + i, q0[locs + i]); }
		mergesumExternals(ext, q0);
		Buffer<bool> done(locs, false);
		const uint basesize = intg1.getBaseSize();
		for(uint i=0; i<basesize; i++) {
			if(!baseflag.includes(intg1.getBaseFlag(i))) continue;
			Buffer< pair<uint, Buffer<double> > > q1;
			intg1.getBaseVector(i, q1);
			for(uint j=0; j<q1.size(); j++) {
				const uint jj = q1[j].first;
				if(!done[jj]) { // create quadrature q0[jj]
					Buffer<double> q;
					intg0.getVector(jj, q);
					q0[jj].mergesum(q);
					intg0.invertVector(q0[jj]);
					done[jj] = true;
				}
				intg0.integrateUnitProductVector(q0[jj], q1[j].second, result.m_val[jj]);
			}
		}
		return result;
	}
	template<typename T> Diagonal<T> &integrateConstantHodge(const T *f, const FormGrade grade, Diagonal<T> &result) const {
		const bool isprim = FormGradeIsPrim(grade);
		MeshIntegrator intg0(m_mesh, (isprim ? grade : FormGradeDual(grade)), 0, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, (!isprim ? grade : FormGradeDual(grade)), 0, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		// integrate dual quadrature
		const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
		Buffer< Buffer<double> > q1(locs + ext.size()); 
		const uint baselocs = intg1.getBaseLocals();
		for(uint i=0; i<baselocs; i++) {
			Buffer< pair<uint, Buffer<double> > > q;
			intg1.getBaseVector(i, q);
			for(uint j=0; j<q.size(); j++) q1[q[j].first].mergesum(q[j].second);
		}
		mergesumExternals(ext, q1);
		for(uint i=0; i<locs; i++) {
			Buffer<double> q0;
			intg0.getVector(i, q0);
			if(isprim) intg0.integrateProductVector(intg0.invertVector(q0), q1[i], f, result.m_val[i]);
			else intg1.integrateProductVector(q0, intg1.invertVector(q1[i]), f, result.m_val[i]);
			q1[i].clear();
		}
		return result;
	}
	template<typename T> Diagonal<T> &integrateConstantHodge(const T *f, const FormGrade grade, const UintSet &flag, Diagonal<T> &result) const {
		MeshIntegrator intg0(m_mesh, grade, 0, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, FormGradeDual(grade), 0, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		if(FormGradeIsPrim(grade)) { // integrate dual quadrature divided by primal
			const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
			Buffer<T> extval(ext.size(), result.m_zero);
			const uint size = locs + ext.size();
			for(uint i=0; i<size; i++) {
				if(!flag.includes(intg0.getFlag(i))) continue;
				Buffer<double> q0;
				intg0.getVector(i, q0);
				intg0.invertVector(q0);
				Buffer<double> q1;
				intg1.getVector(i, q1);
				T &val = (i < locs ? result.m_val[i] : extval[i - locs]);
				intg0.integrateProductVector(q0, q1, f, val);
			}
			shareExternals(ext, extval, result.m_val);
			return result;
		}
		// integrate primal quadrature divided by dual
		const Buffer< pair<uint,uint> > &ext = intg0.getExternals();
		Buffer< Buffer<double> > q0(locs + ext.size()); 
		for(uint i=0; i<ext.size(); i++) {
			const uint ii = locs + i;
			if(flag.includes(intg0.getFlag(ii))) intg0.getVector(ii, q0[ii]);
		}
		mergesumExternals(ext, q0);
		for(uint i=0; i<locs; i++) {
			if(!flag.includes(intg0.getFlag(i))) continue;
			Buffer<double> q;
			q0[i].mergesum(intg0.getVector(i, q));
			intg0.invertVector(q0[i]);
			Buffer<double> q1;
			intg1.getVector(i, q1);
			intg0.integrateProductVector(q0[i], q1, f, result.m_val[i]);
			q0[i].clear();
		}
		return result;
	}
	template<typename T> Diagonal<T> &integrateConstantHodge(const UintSet &baseflag, const T *f, const FormGrade grade, Diagonal<T> &result) const {
		MeshIntegrator intg0(m_mesh, grade, 0, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, FormGradeDual(grade), 0, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		if(FormGradeIsPrim(grade)) { // integrate dual quadrature divided by primal
			const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
			Buffer< Buffer<double> > q1(locs + ext.size()); 
			const uint baselocs = intg1.getBaseLocals();
			for(uint i=0; i<baselocs; i++) {
				if(!baseflag.includes(intg1.getBaseFlag(i))) continue;
				Buffer< pair<uint, Buffer<double> > > q;
				intg1.getBaseVector(i, q);
				for(uint j=0; j<q.size(); j++) q1[q[j].first].mergesum(q[j].second);
			}
			mergesumExternals(ext, q1);
			for(uint i=0; i<locs; i++) {
				if(q1[i].empty()) continue;
				Buffer<double> q0;
				intg0.getVector(i, q0);
				intg0.invertVector(q0);
				intg0.integrateProductVector(q0, q1[i], f, result.m_val[i]);
				q1[i].clear();
			}
			return result;
		}
		// integrate primal quadrature divided by dual
		const Buffer< pair<uint,uint> > &ext = intg0.getExternals();
		Buffer< Buffer<double> > q0(locs + ext.size()); 
		for(uint i=0; i<ext.size(); i++) { intg0.getVector(locs + i, q0[locs + i]); }
		mergesumExternals(ext, q0);
		Buffer<bool> done(locs, false);
		const uint basesize = intg1.getBaseSize();
		for(uint i=0; i<basesize; i++) {
			if(!baseflag.includes(intg1.getBaseFlag(i))) continue;
			Buffer< pair<uint, Buffer<double> > > q1;
			intg1.getBaseVector(i, q1);
			for(uint j=0; j<q1.size(); j++) {
				const uint jj = q1[j].first;
				if(!done[jj]) { // create quadrature q0[jj]
					Buffer<double> q;
					intg0.getVector(jj, q);
					q0[jj].mergesum(q);
					intg0.invertVector(q0[jj]);
					done[jj] = true;
				}
				intg0.integrateProductVector(q0[jj], q1[j].second, f, result.m_val[jj]);
			}
		}
		return result;
	}
	template<typename T> Diagonal<T> &integrateHodge(void func(const Vector4 &, T *), const uint num, const FormGrade grade, Diagonal<T> &result) const {
		const bool isprim = FormGradeIsPrim(grade);
		MeshIntegrator intg0(m_mesh, (isprim ? grade : FormGradeDual(grade)), num, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, (!isprim ? grade : FormGradeDual(grade)), num, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		// integrate dual quadrature
		const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
		Buffer< Buffer<double> > q1(locs + ext.size()); 
		const uint baselocs = intg1.getBaseLocals();
		for(uint i=0; i<baselocs; i++) {
			Buffer< pair<uint, Buffer<double> > > q;
			intg1.getBaseQuadrature(i, q);
			for(uint j=0; j<q.size(); j++) q1[q[j].first].combine(q[j].second);
		}
		combineExternals(ext, q1);
		for(uint i=0; i<locs; i++) {
			const Vector4 p = intg0.getPosition(i);
			Buffer<double> q0;
			intg0.getQuadrature(i, q0);
			if(isprim) intg0.integrateProductQuadrature(intg0.invertQuadrature(p, q0), q1[i], func, result.m_val[i]);
			else intg1.integrateProductQuadrature(q0, intg1.invertQuadrature(p, q1[i]), func, result.m_val[i]);
			q1[i].clear();
		}
		return result;
	}
	template<typename T> Diagonal<T> &integrateHodge(void func(const Vector4 &, T *), const uint num, const FormGrade grade, const UintSet &flag, Diagonal<T> &result) const {
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
				Buffer<double> q0;
				intg0.getQuadrature(i, q0);
				intg0.invertQuadrature(intg0.getPosition(i), q0);
				Buffer<double> q1;
				intg1.getQuadrature(i, q1);
				T &val = (i < locs ? result.m_val[i] : extval[i - locs]);
				intg0.integrateProductQuadrature(q0, q1, func, val);
			}
			shareExternals(ext, extval, result.m_val);
			return result;
		}
		// integrate primal quadrature divided by dual
		const Buffer< pair<uint,uint> > &ext = intg0.getExternals();
		Buffer< Buffer<double> > q0(locs + ext.size()); 
		for(uint i=0; i<ext.size(); i++) {
			const uint ii = locs + i;
			if(flag.includes(intg0.getFlag(ii))) intg0.getQuadrature(ii, q0[ii]);
		}
		combineExternals(ext, q0);
		for(uint i=0; i<locs; i++) {
			if(!flag.includes(intg0.getFlag(i))) continue;
			Buffer<double> q;
			intg0.getQuadrature(i, q);
			q0[i].combine(q);
			intg0.invertQuadrature(intg0.getPosition(i), q0[i]);
			Buffer<double> q1;
			intg1.getQuadrature(i, q1);
			intg0.integrateProductQuadrature(q0[i], q1, func, result.m_val[i]);
			q0[i].clear();
		}
		return result;
	}
	template<typename T> Diagonal<T> &integrateHodge(const UintSet &baseflag, void func(const Vector4 &, T *), const uint num, const FormGrade grade, Diagonal<T> &result) const {
		MeshIntegrator intg0(m_mesh, grade, num, m_lowdim, m_highdim);
		MeshIntegrator intg1(m_mesh, FormGradeDual(grade), num, m_lowdim, m_highdim);
		const uint locs = intg0.getLocals();
		initResult(locs, result);
		if(FormGradeIsPrim(grade)) { // integrate dual quadrature divided by primal
			const Buffer< pair<uint,uint> > &ext = intg1.getExternals();
			Buffer< Buffer<double> > q1(locs + ext.size()); 
			const uint baselocs = intg1.getBaseLocals();
			for(uint i=0; i<baselocs; i++) {
				if(!baseflag.includes(intg1.getBaseFlag(i))) continue;
				Buffer< pair<uint, Buffer<double> > > q;
				intg1.getBaseQuadrature(i, q);
				for(uint j=0; j<q.size(); j++) q1[q[j].first].combine(q[j].second);
			}
			combineExternals(ext, q1);
			for(uint i=0; i<locs; i++) {
				if(q1[i].empty()) continue;
				const Vector4 p = intg0.getPosition(i);
				Buffer<double> q0;
				intg0.getQuadrature(i, q0);
				intg0.invertQuadrature(p, q0);
				intg0.integrateProductQuadrature(q0, q1[i], func, result.m_val[i]);
				q1[i].clear();
			}
			return result;
		}
		// integrate primal quadrature divided by dual
		const Buffer< pair<uint,uint> > &ext = intg0.getExternals();
		Buffer< Buffer<double> > q0(locs + ext.size()); 
		for(uint i=0; i<ext.size(); i++) { intg0.getQuadrature(locs + i, q0[locs + i]); }
		combineExternals(ext, q0);
		Buffer<bool> done(locs, false);
		const uint basesize = intg1.getBaseSize();
		for(uint i=0; i<basesize; i++) {
			if(!baseflag.includes(intg1.getBaseFlag(i))) continue;
			Buffer< pair<uint, Buffer<double> > > q1;
			intg1.getBaseQuadrature(i, q1);
			for(uint j=0; j<q1.size(); j++) {
				const uint jj = q1[j].first;
				if(!done[jj]) { // create quadrature q0[jj]
					Buffer<double> q;
					intg0.getQuadrature(jj, q);
					q0[jj].combine(q);
					intg0.invertQuadrature(intg0.getPosition(jj), q0[jj]);
					done[jj] = true;
				}
				intg0.integrateProductQuadrature(q0[jj], q1[j].second, func, result.m_val[jj]);
			}
		}
		return result;
	}

protected:
	const PartMesh &m_mesh;
	uint m_lowdim;
	uint m_highdim;

	template<typename C> void initResult(const uint locs, C &result) const {
		if(!result.m_full || result.m_height != locs) result.setFullOfZeros(locs);
	}
	void mergesumExternals(const Buffer< pair<uint,uint> > &ext, Buffer< Buffer<double> > &q) const;
	void combineExternals(const Buffer< pair<uint,uint> > &ext, Buffer< Buffer<double> > &q) const;
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
	template<typename T> void integrateElement(const MeshIntegrator &intg, const uint i, void func(const Vector4 &, T *), T &val) const {
		Buffer<double> q;
		intg.getQuadrature(i, q);
		intg.integrateQuadrature(q, func, val);
	}
	template<typename T> void integrateConstantElement(const MeshIntegrator &intg, const uint i, const T *f, T &val) const {
		Buffer<double> q;
		intg.getVector(i, q);
		intg.integrateVector(q, f, val);
	}
	template<typename T> void integrateBaseElement(const MeshIntegrator &intg, const uint i, void func(const Vector4 &, T *), Buffer<T> &locval, Buffer<T> &extval) const {
		const uint locs = locval.size();
		Buffer< pair<uint, Buffer<double> > > q;
		intg.getBaseQuadrature(i, q);
		for(uint j=0; j<q.size(); j++) {
			const uint jj = q[j].first;
			T &val = (jj < locs ? locval[jj] : extval[jj - locs]);
			intg.integrateQuadrature(q[j].second, func, val);
		}
	}
	template<typename T> void integrateConstantBaseElement(const MeshIntegrator &intg, const uint i, const T *f, Buffer<T> &locval, Buffer<T> &extval) const {
		const uint locs = locval.size();
		Buffer< pair<uint, Buffer<double> > > q;
		intg.getBaseVector(i, q);
		for(uint j=0; j<q.size(); j++) {
			const uint jj = q[j].first;
			T &val = (jj < locs ? locval[jj] : extval[jj - locs]);
			intg.integrateVector(q[j].second, f, val);
		}
	}
	Buffer< pair<uint,uint> > getMyExternals(const Buffer< pair<uint,uint> > &ext) const;

};

}

#endif //_DEC_HPP_INCLUDED_
