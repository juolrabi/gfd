/**
 * A program for simulating DC currents in superconductor.
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
	if(plen < 1.0) return 3.0 * (1.0 + cos(PI * plen)) * q[0];
	return 0.0;
}
double get1Form(const Buffer<double> &q) {
	const Vector2 p(q[2], q[3]);
	const double plen = 3.0 * p.len();
	if(plen < 1.0) return 3.0 * (1.0 + cos(PI * plen)) * q[0];
	return q[1];
}
uint functionCut(const Buffer<double> &q) {
	const Vector2 p(q[1], q[2]);
	if(fabs(p.y - 0.5) <= 0.300001 && fabs(p.x - 1.5) > 0.300001) return 1;
	return 0;
}

double sinWave(const double time) {
//	return 200.0 * (time > 3.25 ? 1.0 : sin(PIx2 * time));
	return 50.0 * sin(PIx2 * time);
}
double sourceVector(const Buffer<double> &q) {
	const Vector2 p(q[1], q[2]);
	if((p - Vector2(1.00001,1.00001)).lensq() < 0.01) return -50.0 * q[0];
	if((p - Vector2(1.00001,0.00001)).lensq() < 0.01) return 50.0 * q[0];
	return 0.0;
}
uint functionAbsolute(const sign &s) {
	if(s == 0) return 0;
	return 1;
}
uint functionSelectZeros(const uint &s) {
	if(s == 0) return 1;
	return 0;
}
int main()
{
	initMPI();

	// 2d
	const uint w = 21;
	const double h = (w>1 ? 3.0 / double(w - 1) : 3.0);
	Buffer<uchar> data(w * w, 3);
//	data[4] = 11;
	for(uint i=0; i<w; i++) {
		const uint ii = (i-10) * (i-10);
		for(uint j=0; j<w; j++) {
			const uint jj = (j-10) * (j-10);
			if(ii + jj < 7) data[j * w + i] = 11;
			//else if(ii + jj < 25) data[j * w + i] = 4;
		}
	}
	BlockMesh bm;
	bm.init(Vector2(-1,-1), Vector2(h,h), w, w, &data[0]);

	// create a PartMesh from BlockMesh
	PartMesh mesh(getMPIrank(), getMPIranks(), 2);
	bm.toMesh(mesh);

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
	Derivative d0;
	bm.integrateDerivative(fg_prim0, d0);
	Hodge<double> h1;
	bm.integrateHodge(HodgeUnit2, 0, fg_prim1, h1);
	Hodge<double> dh0;
	bm.integrateHodge(HodgeUnit1, 0, fg_dual0, dh0);

	// initialize forms and system matrices
	Buffer< Form<double> > v(2, Form<double>(0.0));
	Form<double> current(0.0);
	v[0].setFullOfZeros(dh0.m_height); // div current
	v[1].setFullOfZeros(h1.m_height); // current

	Buffer< Sparse<double> > d(2, Sparse<double>(0.0));
	d[0].setTranspose(d0);
	d[0].setTimes(dh0, d[0]);
	d[0].setTimes(d[0], h1);
	d[1].setNegation(d0);

/*	Column<uint> cut(0);
	bm.integrateForm(functionCut, 1, fg_prim0, cut);
	cut.trimSparse();
	Sparse<uint> absd0(0);
	absd0.setFunction(d0, functionAbsolute);
	cut.setTimes(absd0, cut);
	cut.setFunction(cut, functionSelectZeros);
	cut.trim();
	d[1].setTimes(cut, d[1]);
*/	


	Buffer< Form<double> > f(1, Form<double>(0.0));
	bm.integrateForm(sourceVector, 1, fg_prim0, f[0]);
	f[0].trimSparse();

	Sparse<double> test;
	test.setTimes(d[0], d[1]);
	double maxdiag = 0.0;
	for(uint i=0; i<test.m_height; i++) {
		const double diag = test.getValue(i,i);
		if(maxdiag > diag) maxdiag = diag;
	}
	cout << "maxdiag = " << maxdiag << endl;
	d[1].scale(-1.0 / maxdiag);

	for(uint k=0; k<=10; k++) {
		for(uint j=0; j<5000; j++) {
			v[0] = d[0] * v[1] + f[0];
			v[1] += d[1] * v[0];
		}
		double sum = 0.0;
		for(uint j=0; j<v[0].m_height; j++) {
			sum += v[0].m_val[j] * v[0].m_val[j];
		}
		sumMPI(&sum, 1);
		cout << "div sum = " << sum << endl;

		Picture pic0(500,500);
		Buffer<double> val;
		for(uint j=0; j<pic0.getHeight(); j++) {
			for(uint i=0; i<pic0.getWidth(); i++) {
				const Vector4 p(0.01 * i - 2.0, 0.01 * j - 2.0,0,0);
				bm.interpolateForm(fg_prim1, v[1], p, val);
				sumMPI(&val[0], val.size());
//				double sq = 0.0;
//				for(uint x=0; x<val.size(); x++) sq += val[x] * val[x];
//				sq = sqrt(sq);
//				const Vector4 col(sq,sq,sq, 1.0);
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

	finalizeMPI();

    return 0;
}
