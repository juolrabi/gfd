/**
 * A program for simulating wave propagation.
 * The code is under construction.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
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
	if(plen < 1.0) return 3.0 * (1.0 + cos(PI * plen)) * q[0];
	return 0.0;
}
double get1Form(const Buffer<double> &q) {
	const Vector2 p(q[2], q[3]);
	const double plen = 3.0 * p.len();
	if(plen < 1.0) return 3.0 * (1.0 + cos(PI * plen)) * q[0];
	return q[1];
}

double sinWave(const double time) {
//	return 200.0 * (time > 3.25 ? 1.0 : sin(PIx2 * time));
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
int main()
{
	initMPI();
	auto starttime = chrono::system_clock::now();

	// 2d
	const uint w = 21;
	const double h = (w>1 ? 3.0 / double(w - 1) : 3.0);
	Buffer<uchar> data(w * w, 3);
//	data[4] = 11;
	for(uint i=0; i<w; i++) {
		const uint ii = (i-10) * (i-10);
		for(uint j=0; j<w; j++) {
			const uint jj = (j-10) * (j-10);
			if(ii + jj < 7) data[j * w + i] = 5;
			//else if(ii + jj < 25) data[j * w + i] = 4;
		}
	}
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
	// create a PartMesh from BlockMesh
	PartMesh mesh(getMPIrank(), getMPIranks(), 2);
	bm.toMesh(mesh);
/*
{
	FormGrade grade = fg_dual1; 
	auto func = HodgeUnit2;
	const int num = 1;
	Dec dec(mesh, 0, mesh.getDimension());

	Hodge<double> f0;
	auto starttime = chrono::system_clock::now();
	bm.integrateWedge(func, num, grade, f0);
	cout << "f0 time: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;

	Hodge<double> f1;
	starttime = chrono::system_clock::now();
	dec.integrateWedge(func, num, grade, f1);
	cout << "f1 time: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;

	Hodge<double> f2;
	starttime = chrono::system_clock::now();
	dec.integrateWedge(func, num, grade, UINTSETALL, f2);
	cout << "f2 time: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;

	Hodge<double> f3;
	starttime = chrono::system_clock::now();
	dec.integrateWedge(UINTSETALL, func, num, grade, f3);
	cout << "f3 time: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;

	double summ = 0.0;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double sum3 = 0.0;
	for(uint j=0; j<f0.m_height; j++) {
		const double val = f0.getValue(j,j);
		summ += val * val;
		const double dif1 = val - f1.getValue(j,j);
		sum1 += dif1 * dif1;
		const double dif2 = val - f2.getValue(j,j);
		sum2 += dif2 * dif2;
		const double dif3 = val - f3.getValue(j,j);
		sum3 += dif3 * dif3;
	}
	sumMPI(&summ, 1);
	sumMPI(&sum1, 1);
	sumMPI(&sum2, 1);
	sumMPI(&sum3, 1);

	if(getMPIrank() == 0) {
		cout << "summ = " << summ << endl;
		cout << "sum1 = " << sum1 << endl;
		cout << "sum2 = " << sum2 << endl;
		cout << "sum3 = " << sum3 << endl;
	}

}

	finalizeMPI();

    return 0;
*/

	// draw the mesh
	MeshDrawer drawer;
	const Vector3 vo(0,0,0);
	const Vector3 vp(0,0,10);
	const Vector3 vx = TwoVector3(Vector3(0,1,0), vp-vo).dual().unit() / 4.5;
	const Vector3 vy = TwoVector3(vp-vo, vx).dual().unit() / 4.5;
	drawer.initPosition(Vector4(vp,0), Vector4(vo,0), Vector4(vx,0), Vector4(vy,0));
	drawer.initSvg(500, 500);
	drawer.drawBoundaryFaces(mesh, Vector3(0,0.5,1));
	Text name;
	name << "mesh" << getMPIrank() << ".svg";
	drawer.saveSvg(name.str());

	Text text;
	mesh.writeStatistics(text);
	text.save("stat.txt");

	// initialize derivatives and integrate Hodge operators


/*	Derivative d0;
	bm.integrateDerivative(fg_prim0, d0);
	Hodge<double> h1;
	bm.integrateHodge(HodgeUnit2, 0, fg_prim1, h1);
	Hodge<double> dh0;
	bm.integrateHodge(HodgeUnit1, 0, fg_dual0, dh0);

	// initialize forms and system matrices
	Form<double> e(dh0.m_height, 0.0);
	//bm.integrateForm(get0Form, 2, fg_prim0, e);
	Form<double> o(h1.m_height, 0.0);

	Sparse<double> oDe;
	oDe.setNegation(d0);
	Sparse<double> eDo;
	eDo.setTranspose(d0);
	eDo.setTimes(dh0, eDo);
	eDo.setTimes(eDo, h1);
	
	Column<Couple<double> > pml1(Couple<double>(0.0, 0.0));
	bm.integrateForm(getPML1, 1, fg_prim0, pml1);
	Diagonal<double> eAe(0.0);
	eAe.setFunction(pml1, functionDivision);
	eAe.trimSparse();
	Column<Couple<Vector2> > pml2(Couple<Vector2>(Vector2(0,0), Vector2(0,0)));
	bm.integrateForm(getPML2, 1, fg_prim1, pml2);
	Diagonal<double> oAo(0.0);
	oAo.setFunction(pml2, functionDivision);
	oAo.trimSparse();

	//Diagonal<double> eAe(e.m_height, Buffer< pair<uint,double> >(), 0.0);
	//Diagonal<double> oAo(o.m_height, Buffer< pair<uint,double> >(), 0.0);
	Column<double> eF;//(e.m_height, Buffer< pair<uint,double> >(), 0.0);
	bm.integrateForm(sourceVector, 1, fg_prim0, eF);
	eF.trimSparse();


	cout << "koko " << eAe.m_row.size() << " " << eAe.m_val.size() << " " << eAe.m_height << endl;
	TimeIntegrator intg(eDo, eAe, oDe, oAo, eF, 1.0);
*/

	Derivative d0;
	bm.integrateDerivative(fg_prim0, d0);
	Hodge<double> h1;
	bm.integrateHodge(HodgeUnit2, 0, fg_prim1, h1);
	Hodge<double> dh0;
	bm.integrateHodge(HodgeUnit1, 0, fg_dual0, dh0);

	// initialize forms and system matrices
	Buffer< Form<double> > v(2, Form<double>(0.0));
	//v[0].setFullOfZeros(dh0.m_height);
	bm.integrateForm(get0Form, 1, fg_prim0, v[0]);
	v[1].setFullOfZeros(h1.m_height);

	Buffer< Sparse<double> > d(2, Sparse<double>(0.0));
	d[0].setTranspose(d0);
	d[0].setTimes(dh0, d[0]);
	d[0].setTimes(d[0], h1);
	d[1].setNegation(d0);
	
	Buffer< Diagonal<double> > a(2, Diagonal<double>(0.0));
	Column<Couple<double> > pml1(Couple<double>(0.0, 0.0));
	bm.integrateForm(getPML1, 1, fg_prim0, pml1);
	a[0].setFunction(pml1, functionDivision);
	a[0].trimSparse();
	Column<Couple<Vector2> > pml2(Couple<Vector2>(Vector2(0,0), Vector2(0,0)));
	bm.integrateForm(getPML2, 1, fg_prim1, pml2);
	a[1].setFunction(pml2, functionDivision);
	a[1].trimSparse();

	Buffer< Form<double> > f(1, Form<double>(0.0));
	bm.integrateForm(sourceVector, 1, fg_prim0, f[0]);
	f[0].trimSparse();

	Buffer< double (*)(const double) > func(1, sinWave);

	TimeIntegrator intg(d, a, f, 1.0);
//	TimeIntegrator intg(d, 1.0);
	

	for(uint k=0; k<=10; k++) {
		intg.integratePeriod(v, func);
//		intg.integratePeriod(v);

		Picture pic0(500,500);
		Buffer<double> val;
		for(uint j=0; j<pic0.getHeight(); j++) {
			for(uint i=0; i<pic0.getWidth(); i++) {
				const Vector4 p(0.006 * i - 1.0, 0.006 * j - 1.0,0,0);
				bm.interpolateForm(fg_prim0, v[0], p, val);
				sumMPI(&val[0], val.size());
				const Vector4 col((val.size() > 0 ? val[0] : 0.0), (val.size() > 1 ? val[1] : 0.0), (val.size() > 2 ? val[2] : 0.0), 1.0);
				pic0.setColor(i, j, col);
			}
		}
		if(getMPIrank() == 0)
		{
			Text path;
			path << "kuva" << k << ".bmp" << endl;
			pic0.save(path.str(), true);
		}
	}

	// compute check value
	double check = 0.0;
	Buffer<double> val(3, 0.0);
	for(uint i=0; i<v[0].m_val.size(); i++) check += v[0].m_val[i];
	for(uint i=0; i<v[1].m_val.size(); i++) check += v[1].m_val[i];
	sumMPI(&check, 1);
	if(getMPIrank() == 0) {
		cout << "check = " << check << " = 246.182" << endl;
	}




/*	
	Form<double> f0;
	bm.integrateForm(get0Form, 2, fg_prim0, f0);
	Derivative d0;
	bm.integrateDerivative(fg_prim0, d0);
	Derivative d0T;
	d0T.setTranspose(d0);
	Hodge<double> h1;
	bm.integrateHodge(HodgeUnit2, 0, fg_prim1, h1);
	Hodge<double> dh0;
	bm.integrateHodge(HodgeUnit1, 0, fg_dual0, dh0);
	Sparse<double> delta0 = dh0 * d0T * h1;
	Form<double> f1;
	f1.setFullOfZeros(d0.m_height);

	uint i, j;
	Buffer<uint> c0(f0.m_height, 0);
	Buffer<uint> c1(f1.m_height, 0);
	uint steps0 = 1;
{
	ullong sum = 0;
	uint num = 0;
	Sparse<double> test;
	test.setTimes(d0, delta0);
	for(i=0; i<test.m_height; i++) {
		c1[i] = uint(sqrt(0.5 * test.getValue(i, i)) + 1.0);
		sum += c1[i];
		num++;
	}
	test.setTimes(delta0, d0);
	for(i=0; i<test.m_height; i++) {
		c0[i] = uint(sqrt(0.5 * test.getValue(i, i)) + 1.0);
		sum += c0[i];
		num++;
	}
	sumMPI(&sum, 1);
	sumMPI(&num, 1);
	cout << "average " << double(sum) / double(num) << endl;
	steps0 = uint(sum / num);

	for(i=0; i<c0.size(); i++) {
		uint steps = steps0;
		for(j=0; steps<c0[i]; j++) steps *= 2;
		c0[i] = j;
	}
	for(i=0; i<c1.size(); i++) {
		uint steps = steps0;
		for(j=0; steps<c1[i]; j++) steps *= 2;
		c1[i] = j;
	}
//	c0[i] /= uint(sqrt(mindiag));
//	for(i=0; i<c1.size(); i++) c1[i] /= uint(sqrt(mindiag));
//	delta0.scaleRight(-2.0 / mindiag);
}


//	Random rnd(5345345 * getMPIrank());
//	Buffer<uint> c0(f0.m_height, 0);
//	for(i=0; i<c0.size(); i++) c0[i] = rnd.getUint() % 2;
//	Buffer<uint> c1(f1.m_height, 0);
//	for(i=0; i<c1.size(); i++) c1[i] = rnd.getUint() % 2;


	SparseSplit<double> sd0(c1, d0, c0);
	SparseSplit<double> sdelta0(delta0, sd0.getCategories());
	cout << sdelta0.m_term.size() << " " << sd0.m_term.size() << endl;
	uint maxstep = 1;
	for(i=0; i<sdelta0.m_term.size() || i<sd0.m_term.size(); i++) {
		if(i<sd0.m_term.size()) sd0.m_term[i].scaleRight(1.0 / double(steps0 * maxstep));
		maxstep *= 2;
		if(i<sdelta0.m_term.size()) sdelta0.m_term[i].scaleRight(-1.0 / double(steps0 * maxstep));
	}
	delta0.scaleRight(-1.0 / double(steps0 * maxstep));
*/
/*	const uint maxpow = 2;
	uint maxstep = 1;
	for(i=1; i<maxpow; i++) maxstep *= 3;

	Random rnd(5345345 * getMPIrank());
	Buffer<uint> c0(f0.m_height, 0);
	for(i=0; i<c0.size(); i++) c0[i] = rnd.getUint() % maxpow;
	Buffer<uint> c1(f1.m_height, 0);
	for(i=0; i<c1.size(); i++) c1[i] = rnd.getUint() % maxpow;
	SparseSplit<double> sd0(c1, d0);
	SparseSplit<double> sdelta0(c0, delta0);

	double fac = 1.0;
	for(i=1; i<sd0.m_term.size(); i++) {
		fac /= 3.0;
		sd0.m_term[i].scaleRight(fac);
	}
	fac = 1.0;
	for(i=1; i<sdelta0.m_term.size(); i++) {
		fac /= 3.0;
		sdelta0.m_term[i].scaleRight(fac);
	}
*/
/*	// iterate
	for(uint k=0; k<=10*steps0; k++)
	{
		for(i=0; i<maxstep; i++) {
			uint step = maxstep;
			for(j=0; j<sd0.m_term.size(); j++) {
				if((i % step) == (step / 2)) {
					f1 += sd0.m_term[j] * f0;
					if(k == 0) cout << "odd " << j << endl;
				}
				step /= 2;
			}
			step = maxstep;
			for(j=0; j<sdelta0.m_term.size(); j++) {
				step /= 2;
				if((i % step) == (step / 2)) {
					f0 += sdelta0.m_term[j] * f1;
					if(k == 0) cout << "even " << j + 1 << endl;
				}
			}

//			f0 += delta0 * f1;
//			if(k == 0) cout << "even" << endl;
		}
*/
/*		for(i=0; i<2*maxstep; i++) {
			uint step = 2*maxstep;
			for(j=0; j<maxpow; j++) {
				if((i % step) == 0) {
					if(j < sdelta0.m_term.size()) f0 += sdelta0.m_term[j] * f1;
					if(k == 0) cout << "even " << j << endl;
				}
				else if(((i + maxstep) % step) == 0) {
					if(j < sd0.m_term.size()) f1 += sd0.m_term[j] * f0;
					if(k == 0) cout << "odd " << j << endl;
				}
				step /= 3;
			}
		}
*/
/*
		if((k % steps0) != 0) continue;

		Picture pic0(500,500);
		Buffer<double> val;
		for(j=0; j<pic0.getHeight(); j++)
		{
			for(i=0; i<pic0.getWidth(); i++)
			{
				const Vector4 p(0.006 * i - 1.0, 0.006 * j - 1.0,0,0);
				bm.interpolateForm(fg_prim0, f0, p, val);
				sumMPI(&val[0], val.size());
				const Vector4 col((val.size() > 0 ? val[0] : 0.0), (val.size() > 1 ? val[1] : 0.0), (val.size() > 2 ? val[2] : 0.0), 1.0);
				pic0.setColor(i, j, col);
			}
		}
		if(getMPIrank() == 0)
		{
			Text path;
			path << "kuva" << k / steps0 << ".bmp" << endl;
			pic0.save(path.str(), true);
		}
	}

*/
/*	// draw f0 at the end
	double check = 0.0;
	Picture pic0(500,500);
	Buffer<double> val(3, 0.0);
	for(uint j=0; j<pic0.getHeight(); j++) {
		for(uint i=0; i<pic0.getWidth(); i++) {
			const Vector4 p(0.006 * i - 1.0, 0.006 * j - 1.0,0,0);
			bm.interpolateForm(fg_prim0, f0, p, val);
			sumMPI(&val[0], val.size());
			const Vector4 col((val.size() > 0 ? val[0] : 0.0), (val.size() > 1 ? val[1] : 0.0), (val.size() > 2 ? val[2] : 0.0), 1.0);
			check += col.toVector3().len();
			pic0.setColor(i, j, col);
		}
	}
	if(getMPIrank() == 0) {
		cout << "check = " << check << endl;
		pic0.save("kuva.bmp", true);
	}
*/
	if(getMPIrank() == 0) cout << "Elapsed time: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
	finalizeMPI();


    return 0;
}
