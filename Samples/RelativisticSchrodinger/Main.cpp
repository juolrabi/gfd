/**
 * A program for simulating wave propagation.
 * The code is under construction.
 * Author: Jukka Räbinä, University of Jyväskylä, 2020.
 */

#include "../../GFD/BlockDec/BlockMesh.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include "../../GFD/Types/Types.hpp"
#include "../../GFD/Types/Random.hpp"
#include "../../GFD/Types/Buffer.hpp"
#include "../../GFD/Types/MpiEasy.hpp"
#include "../../GFD/Discrete/Split.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include "../../GFD/Discrete/TimeIntegrator.hpp"
#include <iostream>
#include <chrono>

using namespace std;
using namespace gfd;

double get0Form(const Buffer<double> &q) {
	const Vector2 p(q[1], q[2]);
	const double plen = 3.0 * p.len();
	if(plen < 1.0) return 3.0 * (1.0 + cos(PI * plen)) * q[0] * p.x;
	return 0.0;

/*	const Vector2 p(q[1], q[2]);
	const double plen = 3.0 * (p - Vector2(0.1,0.2)).len();
	if(plen < 1.0) return 3.0 * (1.0 + cos(PI * plen)) * q[0];
	return 0.0;
*/}
double get1Form(const Buffer<double> &q) {
	const Vector2 p(q[2], q[3]);
	const double plen = 3.0 * p.len();
	if(plen < 1.0) return 3.0 * (1.0 + cos(PI * plen)) * q[0];
	return q[1];
}
double functionPotential(const Buffer<double> &q) {
	const Vector2 p(q[1], q[2]);
	return -20.0 / p.len();
}

double sinWave(const double time) {
	return 50.0 * sin(PIx2 * time);
}
double sourceVector(const Buffer<double> &q) {
	const Vector2 p(q[1], q[2]);
	if((p - Vector2(1.00001,1.00001)).lensq() < 0.01) return q[0];
	return 0.0;
}
const double PML_MAGNITUDE = 50.0;
Couple<double> getPML1(const Buffer<double> &q) {
	//return Couple<double>(PML_MAGNITUDE, 1.0);
	const Vector2 p(q[1], q[2]);
	double sum = 0.0;
	if(p.x > 1.5) sum += p.x - 1.5;
	else if(-p.x > 0.5) sum += -p.x - 0.5;
	if(p.y > 1.5) sum += p.y - 1.5;
	else if(-p.y > 0.5) sum += -p.y - 0.5;
	if(sum == 0.0) return Couple<double>(0.0, 0.0);
	return Couple<double>(PML_MAGNITUDE * sum * q[0], q[0]);
}
double functionDivision(const Couple<double> &v) {
	if(v.a == 0.0) return 0.0;
	return v.a / v.b;
}
Couple<Vector2> getPML2(const Buffer<double> &q) {
	//return Couple<Vector2>(-PML_MAGNITUDE * Vector2(1,1), Vector2(1,1));
	const Vector2 p(q[2], q[3]);
	Vector2 sum(0,0);
	if(p.x > 1.5) sum += (p.x - 1.5) * Vector2(q[0], -q[1]);
	else if(-p.x > 0.5) sum += (-p.x - 0.5) * Vector2(q[0], -q[1]);
	if(p.y > 1.5) sum += (p.y - 1.5) * Vector2(-q[0], q[1]);
	else if(-p.y > 0.5) sum += (-p.y - 0.5) * Vector2(-q[0], q[1]);
	if(sum == Vector2(0,0)) return Couple<Vector2>(sum, sum);
	return Couple<Vector2>(PML_MAGNITUDE * sum, Vector2(q[0], q[1]));
}
double functionDivision(const Couple<Vector2> &v) {
	if(v.a == Vector2(0,0)) return 0.0;
	return v.a.dot(v.b) / v.b.lensq();
}
double getMaximumDiagonal(const Sparse<double> &d0, const double c0, const Sparse<double> &d1, const double c1) {
	double maxval = 0.0;
	Sparse<double> test;
	test.setTimes(d0, d1);
	for(uint i=0; i<test.m_height; i++) {
		const double val = c0 + test.getValue(i, i);
		if(maxval < val) maxval = val;
	}
	test.setTimes(d1, d0);
	for(uint i=0; i<test.m_height; i++) {
		const double val = c1 + test.getValue(i, i);
		if(maxval < val) maxval = val;
	}
	return maxval;
}
double getEigenMode(Column<double> &v, const Sparse<double> &dd, const Diagonal<double> &pot, const double sigma = 0.0) {
	Sparse<double> ddpot;
	ddpot.setPlus(dd, pot);
	const double dtime = sqrt(2.0 / getMaximumDiagonal(ddpot, 0.0, ddpot, 0.0));
	cout << "dtime = " << dtime << endl;

	// search for initial mode
	double E = 0.0;
	uint i;
	for(i=0; i<3000; i++) {
		Form<double> grad(0.0);
		grad.setTimes(ddpot, v);
		grad -= (pot * v).scale(sigma * E * (i < 100 ? 0.01 * i : 1.0));
		double sumdot = 0.0;
		double sumsq = 0.0;
		for(uint j=0; j<v.m_height; j++) {
			sumdot -= v.m_val[j] * grad.m_val[j];
			sumsq += v.m_val[j] * v.m_val[j];
		}
		sumMPI(&sumdot, 1);
		sumMPI(&sumsq, 1);
		const double newE = sumdot / (sumsq * (1.0 - sigma * E));
		cout << "E = " << newE << endl;
		if(i > 100 && fabs(E - newE) < 1e-4) break;
		E = newE;

		// update and normalize v
		grad.scale(dtime);
		v -= grad;
		const double div = 0.5 * sqrt(double(v.m_height) / sumsq);
		for(uint j=0; j<v.m_height; j++) v.m_val[j] *= div;
	}
	cout << "ended at i = " << i << endl;
	return E;
}

int main()
{
	initMPI();

	// 2d
	const uint w = 41;
	const double h = (w>1 ? 2.0 / double(w - 1) : 2.0);
	Buffer<uchar> data(w * w, 3);
	BlockMesh bm;
	bm.init(Vector2(-1,-1), Vector2(h,h), w, w, &data[0]);

/*	// 3d
	const uint w = 11;
	const double h = (w>1 ? 3.0 / double(w - 1) : 3.0);
	Buffer<uchar> data(w * w * w, 3);
	for(uint i=0; i<w; i++) {
		const uint ii = (i-5) * (i-5);
		for(uint j=0; j<w; j++) {
			const uint jj = (j-5) * (j-5);
			for(uint k=0; k<w; k++) {
				const uint kk = (k-5) * (k-5);
				if(ii + jj + kk < 7) data[(k * w + j) * w + i] = 4;
			}
		}
	}
	BlockMesh bm;
	bm.init(Vector3(-1,-1, -1), Vector3(h,h,h), w, w, w, &data[0]);
*/

	// relativistic Schrödinger
	const double sigma = 0.001;
	Derivative d0;
	bm.integrateDerivative(fg_prim1, d0);
	Hodge<double> h1;
	bm.integrateHodge(HodgeUnit1, 0, fg_prim2, h1);
	Hodge<double> dh0;
	bm.integrateHodge(HodgeUnit2, 0, fg_dual1, dh0);

	// initialize forms and system matrices
	Buffer< Form<double> > v(4, Form<double>(0.0));
	bm.integrateForm(get0Form, 1, fg_dual2, v[0]);
	//v[0].setFullOfZeros(dh0.m_height);
	v[1].setFullOfZeros(h1.m_height);
	v[2].setFullOfZeros(dh0.m_height);
	v[3].setFullOfZeros(dh0.m_height);

	Diagonal<double> pot(0.0);
	bm.integrateForm(functionPotential, 1, fg_dual2, pot);
	double maxsq = 0.0;
	for(uint i=0; i<pot.m_height; i++) {
		const double sq = pot.m_val[i] * pot.m_val[i];
		if(maxsq < sq) maxsq = sq;
	}
	cout << "maxsq = " << maxsq << endl;


	Buffer< Sparse<double> > d(2, Sparse<double>(0.0));
	d[0].setTranspose(d0);
	d[1].setTimes(d0, dh0);
	d[1].setTimes(h1, d[1]);
	const double E = getEigenMode(v[0], d[1] * d[0], pot, sigma);
	d[0].scale(1.0 / sigma);
	const double dtime = sqrt(2.0 / getMaximumDiagonal(d[0], 1.0 / (sigma * sigma), d[1], maxsq));
	cout << "dtime = " << dtime << endl;
	const double dtime_sigma = dtime / sigma;
	d[0].scale(dtime);
	d[1].scale(dtime);
	pot.scale(dtime);

	const uint steps = uint(0.1 * PIx2 / (E * dtime));
	cout << "total time to simulate = " << 10 * steps * dtime << endl;
	v[2].setTimes(d[0], v[0]);
	v[2].scale(1.0 / (dtime * E - dtime_sigma));

//	const uint steps = uint(0.02 / dtime + 0.5);

	for(uint k=0; k<10; k++) {
		Picture pic0(500,500);
		Buffer<double> valr;
		Buffer<double> vali;
		for(uint j=0; j<pic0.getHeight(); j++) {
			for(uint i=0; i<pic0.getWidth(); i++) {
				const Vector4 p(0.004 * i - 1.0, 0.004 * j - 1.0,0,0);
				bm.interpolateForm(fg_dual2, v[0], p, valr);
				sumMPI(&valr[0], valr.size());
				bm.interpolateForm(fg_dual2, v[1], p, vali);
				sumMPI(&vali[0], vali.size());
				const double sq = sqrt(valr[0] * valr[0] + vali[0] * vali[0]);
				//const Vector4 col(sq, sq, sq, 1.0);
				const Vector4 col(valr[0], vali[0], sq, 1.0);
				pic0.setColor(i, j, col);
			}
		}
		if(getMPIrank() == 0)
		{
			Text path;
			path << "kuva" << sigma << "_" << k << ".bmp" << endl;
			pic0.save(path.str(), false);
		}

		// iterate 
		for(uint i=0; i<steps; i++) {
			// odd (imaginary terms)
			v[1] -= pot * v[0];
			v[1] += d[1] * v[2];
			v[3] += scaled(v[2], dtime_sigma);
			v[3] += d[0] * v[0];
			//even (real terms)
			v[0] += pot * v[1];
			v[0] -= d[1] * v[3];
			v[2] -= scaled(v[3], dtime_sigma);
			v[2] -= d[0] * v[1];
		}
	}

/*
	// non-relativistic Schrödinger
	Derivative d0;
	bm.integrateDerivative(fg_prim1, d0);
	Hodge<double> h1;
	bm.integrateHodge(HodgeUnit1, 0, fg_prim2, h1);
	Hodge<double> dh0;
	bm.integrateHodge(HodgeUnit2, 0, fg_dual1, dh0);

	// initialize forms and system matrices
	Buffer< Form<double> > v(2, Form<double>(0.0));
	bm.integrateForm(get0Form, 1, fg_dual2, v[0]);
	//v[0].setFullOfZeros(dh0.m_height);
	v[1].setFullOfZeros(h1.m_height);

	Diagonal<double> pot(0.0);
	bm.integrateForm(functionPotential, 1, fg_dual2, pot);

	Buffer< Sparse<double> > d(1, Sparse<double>(0.0));
	d[0].setTranspose(d0);
	d[0].setTimes(dh0, d[0]);
	d[0].setTimes(d0, d[0]);
	d[0].setTimes(h1, d[0]);
	const double E = getEigenMode(v[0], d[0], pot, 0.0);
//	const double E = 10.0 * PIx2;
	d[0].setPlus(d[0], pot);
	const double dtime = sqrt(2.0 / getMaximumDiagonal(d[0], 0.0, d[0], 0.0));
	cout << "dtime = " << dtime << endl;
	d[0].scale(dtime);

	const uint steps = uint(0.1 * PIx2 / (fabs(E) * dtime) + 0.5);
	cout << "total time to simulate = " << 10 * steps * dtime << endl;
	v[1].setFullOfZeros(v[1].m_height);

	for(uint k=0; k<10; k++) {

		Picture pic0(500,500);
		Buffer<double> valr;
		Buffer<double> vali;
		for(uint j=0; j<pic0.getHeight(); j++) {
			for(uint i=0; i<pic0.getWidth(); i++) {
				const Vector4 p(0.004 * i - 1.0, 0.004 * j - 1.0,0,0);
				bm.interpolateForm(fg_dual2, v[0], p, valr);
				sumMPI(&valr[0], valr.size());
				bm.interpolateForm(fg_dual2, v[1], p, vali);
				sumMPI(&vali[0], vali.size());
				const double sq = sqrt(valr[0] * valr[0] + vali[0] * vali[0]);
				//const Vector4 col(sq, sq, sq, 1.0);
				const Vector4 col(valr[0], vali[0], sq, 1.0);
				pic0.setColor(i, j, col);
			}
		}
		if(getMPIrank() == 0)
		{
			Text path;
			path << "kuva0_" << k << ".bmp" << endl;
			pic0.save(path.str(), false);
		}

		// iterate 
		for(uint i=0; i<steps; i++) { // schrödinger
			// odd (imaginary terms)
			v[1] -= d[0] * v[0];
			//even (real terms)
			v[0] += d[0] * v[1];
		}
	}
*/

	finalizeMPI();

    return 0;
}
