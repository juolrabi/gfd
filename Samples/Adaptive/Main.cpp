/**
 * A program for testing recursive adaptation method.
 * The code is under construction.
 * Author: Jukka Räbinä, University of Jyväskylä, 2020.
 */

#include "Adaptive.hpp"
#include "../../GFD/Output/MeshDrawer.hpp"
#include "../../GFD/Discrete/Dec.hpp"
#include "../../GFD/Discrete/Split.hpp"
#include <iostream>
#include <chrono>

using namespace std;
using namespace gfd;

double get1Form(const Buffer<double> &q) {
	const Vector2 p(q[2], q[3]);
	const double plen = 3.0 * p.len();
	if(plen < 1.0) return 3.0 * (1.0 + cos(PI * plen)) * q[0];
	return q[1];
}


Matrix4 interpolate1Form(const uint node, const PartMesh &mesh, const Buffer<double> &val) {
	if(node >= mesh.getNodeSize()) return Matrix4(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1);
	Buffer<uint> par = mesh.getNodeEdges(node);
	uint pars = par.size();
	const Buffer<uint> vol = mesh.getNodeFaces(node);
	for(uint i=0; i<vol.size(); i++) {
		const Buffer<uint> &e = mesh.getFaceEdges(vol[i]);
		for(uint j=0; j<e.size(); j++) par.gatherOnce(e[j], pars);
	}
	par.resize(pars);
//	const Buffer<uint> &par = mesh.getNodeEdges(node);
	if(par.size() >= 12) {
		MatrixN A;
		VectorN b;
		for(uint i=0; i<par.size(); i++) {
			const Vector3 p = mesh.getEdgeAverage(par[i]).toVector3();
			const Vector3 v = mesh.getEdgeVector(par[i]).toVector3();
			VectorN w(12);
			w[0] = p.x * v.x; w[1] = p.y * v.x; w[2] = p.z * v.x; w[3] = v.x;
			w[4] = p.x * v.y; w[5] = p.y * v.y; w[6] = p.z * v.y; w[7] = v.y;
			w[8] = p.x * v.z; w[9] = p.y * v.z; w[10] = p.z * v.z; w[11] = v.z;
			A += w.outerProduct(w);
			b += w * val[par[i]];
		}
		const VectorN x = A.inverse() * b;
		if(x.size() == 12) return Matrix4(x[0],x[1],x[2],x[3], x[4],x[5],x[6],x[7], x[8],x[9],x[10],x[11], 0,0,0,1);
	}
	if(par.size() >= 3) {
		SymMatrix3 A(ZEROSYMMATRIX3);
		Vector3 b(0,0,0);
		for(uint i=0; i<par.size(); i++) {
			const Vector3 v = mesh.getEdgeVector(par[i]).toVector3();
			A += v.outerProduct();
			b += v * val[par[i]];
		}
		const Vector3 x = A.inverse() * b;
		return Matrix4(0,0,0,x.x, 0,0,0,x.y, 0,0,0,x.z, 0,0,0,1);
	}
	return Matrix4(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1);
}

int main()
{
	initMPI();

//	Term<bool> cubic(false);
//	Term<bool> bcc(true);
	const double adaptr = 1.0;
	const double adapth = 0.1;
	Adaptive adapt(Vector2(-adaptr, -adaptr), Vector2(adaptr, adaptr), adapth);
	adapt.refineSphereBoundary(Vector3(0,0,0), 0.7, 4);
//	abcc.fillSphere(Vector3(0,0,0), 0.7, &bcc);
	adapt.splitPixel(10);

	uint i, j;
	const double h = 0.004;
	Picture pic0(500,500);
	for(j=0; j<pic0.getHeight(); j++) {
		for(i=0; i<pic0.getWidth(); i++) {
			const Vector3 p(h * i - 1.0, h * j - 1.0, 0);
			const double level = 0.17 * adapt.findLevel(p);
			pic0.setColor(i, j, Vector4(level,level,level,1));
		}
	}
	pic0.save("kuva.bmp", false);

/*	const uint index = adapt.findIndex(Vector2(0.687, 0.067));

	Buffer<uint> trace = adapt.findTrace(Vector2(0.687, 0.067));
	cout << "trace";
	for(i=0; i<trace.size(); i++) cout << " " << trace[i];
	cout << endl;

	trace = adapt.getTrace(index);
	cout << "trace";
	for(i=0; i<trace.size(); i++) cout << " " << trace[i];
	cout << endl;
*/
	auto starttime = chrono::system_clock::now();
	BuilderMesh mesh(2);
	adapt.createMesh(mesh);
	cout << "mesh created: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;


	// draw the mesh
	MeshDrawer drawer;
	const Vector3 vo(0,0,0);
	const Vector3 vp(2.5,1.5,10);
	const Vector3 vx = TwoVector3(Vector3(0,1,0), vp-vo).dual().unit() / 2.5;
	const Vector3 vy = TwoVector3(vp-vo, vx).dual().unit() / 2.5;
	drawer.initPosition(Vector4(vp,0), Vector4(vo,0), Vector4(vx,0), Vector4(vy,0));
	drawer.initSvg(500, 500);
	drawer.drawBoundaryFaces(mesh, Vector3(0,0.5,1));
	drawer.saveSvg("mesh.svg");

	cout << "mesh drawn: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;

	// save mesh statistics
	Text text;
	mesh.writeStatistics(text);
	text.save("stat.txt");


/*
	cout << "adapt generated: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;

	// create mesh
	BuilderMesh mesh(3);
	adapt.createMesh(mesh);
	cout << "mesh generated: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;

	PartMesh pmesh(0, 1, 3);
	pmesh.createPartFromFlags(mesh);
	mesh.clear();
	cout << "pmesh generated: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;

	Dec dec(pmesh, 0, pmesh.getDimension());

	// initialize forms and integrate initial solution
	Form<double> f1;
	dec.integrateForm(get1Form, 2, fg_prim1, f1);
	cout << "form integrated: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;


	// initialize derivatives and integrate Hodge operators
	Derivative d1;
	dec.integrateDerivative(fg_prim1, d1);
	Derivative d1T;
	d1T.setTranspose(d1);
	Hodge<double> h2;
	dec.integrateHodge(HodgeUnit1, 2, fg_prim2, h2);
	Hodge<double> dh1;
	dec.integrateHodge(HodgeUnit2, 2, fg_dual1, dh1);
	Sparse<double> delta1 = dh1 * d1T * h2;
	Form<double> f2;
	f2.setFullOfZeros(d1.m_height);
	cout << "operators integrated: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;

	uint steps0 = 1;
	Buffer<uint> c1(f1.m_height, 0);
	Buffer<uint> c2(f2.m_height, 0);
{
	ullong sum = 0;
	uint num = 0;
	Sparse<double> test;
	test.setTimes(d1, delta1);
	for(i=0; i<test.m_height; i++) {
		c2[i] = uint(sqrt(0.5 * test.getValue(i, i)) + 1.0);
		sum += c2[i];
		num++;
	}
	test.setTimes(delta1, d1);
	for(i=0; i<test.m_height; i++) {
		c1[i] = uint(sqrt(0.5 * test.getValue(i, i)) + 1.0);
		sum += c1[i];
		num++;
	}
	sumMPI(&sum, 1);
	sumMPI(&num, 1);
	cout << "average " << double(sum) / double(num) << endl;
	steps0 = uint(sum / num);

	for(i=0; i<c1.size(); i++) {
		uint steps = steps0;
		for(j=0; steps<c1[i]; j++) steps *= 2;
		c1[i] = j;
	}
	for(i=0; i<c2.size(); i++) {
		uint steps = steps0;
		for(j=0; steps<c2[i]; j++) steps *= 2;
		c2[i] = j;
	}
//	c0[i] /= uint(sqrt(mindiag));
//	for(i=0; i<c1.size(); i++) c1[i] /= uint(sqrt(mindiag));
//	delta0.scaleRight(-2.0 / mindiag);
}

	SparseSplit<double> sd1(c2, d1, c1);
	SparseSplit<double> sdelta1(delta1, sd1.getCategories());
	cout << sdelta1.m_term.size() << " " << sd1.m_term.size() << endl;
	uint maxstep = 1;
	for(i=0; i<sdelta1.m_term.size() || i<sd1.m_term.size(); i++) {
		if(i<sd1.m_term.size()) sd1.m_term[i].scaleRight(1.0 / double(steps0 * maxstep));
		maxstep *= 2;
		if(i<sdelta1.m_term.size()) sdelta1.m_term[i].scaleRight(-1.0 / double(steps0 * maxstep));
	}
	delta1.scaleRight(-1.0 / double(steps0 * maxstep));

	// iterate
	for(uint k=0; k<=10*steps0; k++)
	{
//		f2 += d1 * f1;
//		f1 -= (delta1 * f2).scaleRight(0.0005);
		for(i=0; i<maxstep; i++) {
			uint step = maxstep;
			for(j=0; j<sd1.m_term.size(); j++) {
				if((i % step) == (step / 2)) {
					f2 += sd1.m_term[j] * f1;
					if(k == 0) cout << "odd " << j << endl;
				}
				step /= 2;
			}
			step = maxstep;
			for(j=0; j<sdelta1.m_term.size(); j++) {
				step /= 2;
				if((i % step) == (step / 2)) {
					f1 += sdelta1.m_term[j] * f2;
					if(k == 0) cout << "even " << j + 1 << endl;
				}
			}
		}

		if((k % steps0) != 0) continue;

		uint node = NONE;
		Matrix4 val(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1);
		for(j=0; j<pic0.getHeight(); j++) {
			for(i=0; i<pic0.getWidth(); i++) {
				const Vector3 p(h * i - 1.0, h * j - 1.0, 0.0);
				const uint newnode = pmesh.findNode(Vector4(p, 0), abcch * abcch, node, false);
				if(newnode != node) {
					node = newnode;
					val = interpolate1Form(node, pmesh, f1.m_val);
				}
				pic0.setColor(i, j, val * Vector4(p,1));
			}
		}
		Text path;
		path << "field" << k / steps0 << ".bmp" << endl;
		pic0.save(path.str(), true);
		cout << "field drawn: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
	}

	// draw the mesh
	MeshDrawer drawer;
	const Vector3 vo(0,0,0);
	const Vector3 vp(2.5,1.5,-10);
	const Vector3 vx = TwoVector3(Vector3(0,1,0), vp-vo).dual().unit() / 2.5;
	const Vector3 vy = TwoVector3(vp-vo, vx).dual().unit() / 2.5;
	drawer.initPosition(Vector4(vp,0), Vector4(vo,0), Vector4(vx,0), Vector4(vy,0));
	drawer.initSvg(500, 500);
	drawer.drawBoundaryFaces(pmesh, Vector3(0,0.5,1));
	drawer.saveSvg("mesh.svg");

	cout << "mesh drawn: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;

	// save mesh statistics
	Text text;
	pmesh.writeStatistics(text);
	text.save("stat.txt");
*/
	cout << "Elapsed time: " << chrono::duration<double>(chrono::system_clock::now() - starttime).count() << " seconds" << endl;
	finalizeMPI();

    return 0;
}
