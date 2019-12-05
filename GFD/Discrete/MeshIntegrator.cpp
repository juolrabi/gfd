#include "MeshIntegrator.hpp"

using namespace gfd;

MeshIntegrator::MeshIntegrator(const PartMesh &mesh, const FormGrade grade, const uint num, const uint lowdim, const uint highdim)
 : m_mesh(mesh) {
	m_grade = grade;
	m_num = num;
	m_lowdim = (lowdim < 4 ? lowdim : 4);
	m_highdim = (highdim < 4 ? highdim : 4);
}

uint MeshIntegrator::getFields() const { 
	const uint dim = m_mesh.getDimension();
	switch (FormGradeDimension(m_grade)) {
	case 1: return dim;
	case 2: return dim * (dim - 1) / 2;
	case 3: return dim * (dim - 1) * (dim - 2) / 6;
	default: return 1;
	}
}
uint MeshIntegrator::getLocals(const uint gdim) const {
	switch (gdim) {
	case 0: return m_mesh.getNodeLocals();
	case 1: return m_mesh.getEdgeLocals();
	case 2: return m_mesh.getFaceLocals();
	case 3: return m_mesh.getBodyLocals();
	default: return m_mesh.getQuadLocals();
	}
}
uint MeshIntegrator::getSize(const uint gdim) const {
	switch (gdim) {
	case 0: return m_mesh.getNodeSize();
	case 1: return m_mesh.getEdgeSize();
	case 2: return m_mesh.getFaceSize();
	case 3: return m_mesh.getBodySize();
	default: return m_mesh.getQuadSize();
	}
}
const Buffer< pair<uint,uint> > &MeshIntegrator::getExternals(const uint gdim) const {
	switch (gdim) {
	case 0: return m_mesh.getExternalNodes();
	case 1: return m_mesh.getExternalEdges();
	case 2: return m_mesh.getExternalFaces();
	case 3: return m_mesh.getExternalBodies();
	default: return m_mesh.getExternalQuads();
	}
}
uint MeshIntegrator::getFlag(const uint i, const uint gdim) const {
	switch (gdim) {
	case 0: return m_mesh.getNodeFlag(i);
	case 1: return m_mesh.getEdgeFlag(i);
	case 2: return m_mesh.getFaceFlag(i);
	case 3: return m_mesh.getBodyFlag(i);
	default: return m_mesh.getQuadFlag(i);
	}
}

Vector4 MeshIntegrator::getPosition(const uint i, const uint gdim) const {
	switch (gdim) {
	case 0: return m_mesh.getNodePosition(i);
	case 1: return m_mesh.getEdgePosition(i);
	case 2: return m_mesh.getFacePosition(i);
	case 3: return m_mesh.getBodyPosition(i);
	default: return m_mesh.getQuadPosition(i);
	}
}

Vector4 unitEdge(const PartMesh &m, const uint j) { return m.getEdgeVector(j).unit(); }
TwoVector4 unitFace(const PartMesh &m, const uint j) { return m.getFaceVector(j).unit(); }
ThreeVector4 unitBody(const PartMesh &m, const uint j) { return m.getBodyVector(j).unit(); }
FourVector4 unitQuad(const PartMesh &m, const uint j) { return m.getQuadVector(j).unit(); }
double vectorNodeEdge(const Vector4 *p, const Vector4 &b) { return FourVector4(p[1] - p[0], b.dual()).dualof(); }
double vectorNodeFace(const Vector4 *p, const TwoVector4 &b) { return FourVector4(p[1] - p[0], p[2] - p[0], b.dual()).dualof() / 2.0; }
double vectorNodeBody(const Vector4 *p, const ThreeVector4 &b) { return FourVector4(p[1] - p[0], p[2] - p[0], p[3] - p[0], b.dual()).dualof() / 6.0; }
double vectorNodeQuad(const Vector4 *p, const FourVector4 &b) { return FourVector4(p[1] - p[0], p[2] - p[0], p[3] - p[0], p[4] - p[0]).dualof() * b.dual() / 24.0; }
Vector4 vectorEdgeNode(const Vector4 *p, const double &b) { return (p[1] - p[0]) * b; }
Vector4 vectorEdgeFace(const Vector4 *p, const TwoVector4 &b) { return ThreeVector4(p[1] - p[0], b.dual()).dualof(); }
Vector4 vectorEdgeBody(const Vector4 *p, const ThreeVector4 &b) { return ThreeVector4(p[1] - p[0], p[2] - p[0], b.dual()).dualof() / 2.0; }
Vector4 vectorEdgeQuad(const Vector4 *p, const FourVector4 &b) { return ThreeVector4(p[1] - p[0], p[2] - p[0], p[3] - p[0]).dualof() * (b.dual() / 6.0); }
TwoVector4 vectorFaceNode(const Vector4 *p, const double &b) { return TwoVector4(p[1] - p[0], p[2] - p[0]) * (b / 2.0); }
TwoVector4 vectorFaceEdge(const Vector4 *p, const Vector4 &b) { return TwoVector4(p[1] - p[0], b); }
TwoVector4 vectorFaceBody(const Vector4 *p, const ThreeVector4 &b) { return TwoVector4(p[1] - p[0], b.dual()).dualof(); }
TwoVector4 vectorFaceQuad(const Vector4 *p, const FourVector4 &b) { return TwoVector4(p[1] - p[0], p[2] - p[0]).dualof() * (b.dual() / 2.0); }
ThreeVector4 vectorBodyNode(const Vector4 *p, const double &b) { return ThreeVector4(p[1] - p[0], p[2] - p[0], p[3] - p[0]) * (b / 6.0); }
ThreeVector4 vectorBodyEdge(const Vector4 *p, const Vector4 &b) { return ThreeVector4(p[1] - p[0], p[2] - p[0], b) / 2.0; }
ThreeVector4 vectorBodyFace(const Vector4 *p, const TwoVector4 &b) { return ThreeVector4(p[1] - p[0], b); }
ThreeVector4 vectorBodyQuad(const Vector4 *p, const FourVector4 &b) { return (p[1] - p[0]).dualof() * b.dual(); }
FourVector4 vectorQuadNode(const Vector4 *p, const double &b) { return FourVector4(p[1] - p[0], p[2] - p[0], p[3] - p[0], p[4] - p[0]) * (b / 24.0); }
FourVector4 vectorQuadEdge(const Vector4 *p, const Vector4 &b) { return FourVector4(p[1] - p[0], p[2] - p[0], p[3] - p[0], b) / 6.0; }
FourVector4 vectorQuadFace(const Vector4 *p, const TwoVector4 &b) { return FourVector4(p[1] - p[0], p[2] - p[0], b) / 2.0; }
FourVector4 vectorQuadBody(const Vector4 *p, const ThreeVector4 &b) { return FourVector4(p[1] - p[0], b); }

double vectorNodeEdge(const Vector4 &v, const Vector4 &b) { return FourVector4(v, b.dual()).dualof(); }
double vectorNodeFace(const TwoVector4 &v, const TwoVector4 &b) { return FourVector4(v, b.dual()).dualof(); }
double vectorNodeBody(const ThreeVector4 &v, const ThreeVector4 &b) { return FourVector4(v, b.dual()).dualof(); }
double vectorNodeQuad(const FourVector4 &v, const FourVector4 &b) { return v.dualof() * b.dual(); }
Vector4 vectorEdgeNode(const Vector4 &v, const double &b) { return v * b; }
Vector4 vectorEdgeFace(const Vector4 &v, const TwoVector4 &b) { return ThreeVector4(v, b.dual()).dualof(); }
Vector4 vectorEdgeBody(const TwoVector4 &v, const ThreeVector4 &b) { return ThreeVector4(v, b.dual()).dualof(); }
Vector4 vectorEdgeQuad(const ThreeVector4 &v, const FourVector4 &b) { return v.dualof() * b.dual(); }
TwoVector4 vectorFaceNode(const TwoVector4 &v, const double &b) { return v * b; }
TwoVector4 vectorFaceEdge(const Vector4 &v, const Vector4 &b) { return TwoVector4(v, b); }
TwoVector4 vectorFaceBody(const Vector4 &v, const ThreeVector4 &b) { return TwoVector4(v, b.dual()).dualof(); }
TwoVector4 vectorFaceQuad(const TwoVector4 &v, const FourVector4 &b) { return v.dualof() * b.dual(); }
ThreeVector4 vectorBodyNode(const ThreeVector4 &v, const double &b) { return v * b; }
ThreeVector4 vectorBodyEdge(const TwoVector4 &v, const Vector4 &b) { return ThreeVector4(v, b); }
ThreeVector4 vectorBodyFace(const Vector4 &v, const TwoVector4 &b) { return ThreeVector4(v, b); }
ThreeVector4 vectorBodyQuad(const Vector4 &v, const FourVector4 &b) { return v.dualof() * b.dual(); }
FourVector4 vectorQuadNode(const FourVector4 &v, const double &b) { return v * b; }
FourVector4 vectorQuadEdge(const ThreeVector4 &v, const Vector4 &b) { return FourVector4(v, b); }
FourVector4 vectorQuadFace(const TwoVector4 &v, const TwoVector4 &b) { return FourVector4(v, b); }
FourVector4 vectorQuadBody(const Vector4 &v, const ThreeVector4 &b) { return FourVector4(v, b); }

Vector4 vectorEdge(const Vector4 *p) { return p[1] - p[0]; }
TwoVector4 vectorFace(const Vector4 *p) { return TwoVector4(p[1] - p[0], (p[2] - p[0]) / 2.0); }
ThreeVector4 vectorBody(const Vector4 *p) { return ThreeVector4(p[1] - p[0], p[2] - p[0], p[3] - p[0]) / 6.0; }
FourVector4 vectorQuad(const Vector4 *p) { return FourVector4(p[1] - p[0], p[2] - p[0], p[3] - p[0], p[4] - p[0]) / 24.0; }
void MeshIntegrator::removeExternalSimplices(const uint locs, Buffer<uint> &j, Buffer<Vector4> &p) const {
	uint elems = j.size();
	if(elems == 0) return;
	const uint pdim = p.size() / elems;
	for(uint i=elems; i>0; ) {
		if(j[--i] < locs) continue;
		j[i] = j[--elems];
		for(uint k=0; k<pdim; k++) p[pdim * i + k] = p[pdim * elems + k];
		j.resize(elems); // resize tables here because this is very rare event
		p.resize(pdim * elems);
	}
}
template<typename V> void MeshIntegrator::removeExternalVectors(const uint locs, Buffer<uint> &j, Buffer<V> &v) const {
	for(uint i=j.size(); i>0; ) {
		if(j[--i] < locs) continue;
		j[i] = j.back();
		v[i] = v.back();
		j.pop_back();
		v.pop_back();
	}
}
template<typename V> Buffer<double> &MeshIntegrator::createEntryQuadrature(const uint pdim, const Buffer<Vector4> &p, V vector(const Vector4 *), Buffer<double> &q) const {
	const uint elems = p.size() / pdim;
	if(elems == 0) return q;
	uint vs = 0;
	Buffer<double> v(elems * getFields());
	for(uint i=0; i<elems; i++) setVector(vector(&p[pdim * i]), v, vs);
	return createMultiQuadrature(elems, p, v, q);
}
template<typename V, typename B> Buffer<double> &MeshIntegrator::createQuadrature(const Buffer<uint> &j, const Buffer<Vector4> &p, V vector(const Vector4 *, const B &), B base(const PartMesh &, const uint), Buffer<double> &q) const {
	const uint elems = j.size();
	if(elems == 0) return q;
	const uint pdim = p.size() / elems;
	// compute base vectors
	uint bs = 0;
	Buffer<B> b(elems);
	Buffer<uint> link(elems, elems);
	for(uint i=0; i<elems; i++) {
		if(link[i] != elems) continue;
		b[bs] = base(m_mesh, j[i]);
		link[i] = bs;
		for(uint k=i+1; k<elems; k++) {
			if(j[k] == j[i]) link[k] = bs;
		}
		bs++;
	}
	// compute quadrature
	uint vs = 0;
	Buffer<double> v(elems * getFields());
	for(uint i=0; i<elems; i++) setVector(vector(&p[pdim * i], b[link[i]]), v, vs);
	return createMultiQuadrature(elems, p, v, q);
}
template<typename V, typename B> Buffer< pair<uint, Buffer<double> > > &MeshIntegrator::createBaseQuadrature(const Buffer<uint> &j, const Buffer<Vector4> &p, V vector(const Vector4 *, const B &), const B &b, Buffer< pair<uint, Buffer<double> > > &q) const {
	const uint elems = j.size();
	if(elems == 0) return q;
	const uint pdim = p.size() / elems;
	// create groups by table j
	uint qjs = 0;
	Buffer<uint> qj(elems);
	for(uint i=0; i<elems; i++) qj.gatherOnce(j[i], qjs);
	// compute quadratures for each group
	const uint fields = getFields();
	q.resize(qjs);
	for(uint i=0; i<qjs; i++) {
		uint size = 0;
		for(uint k=0; k<elems; k++) {
			if(j[k] == qj[i]) size++;
		}
		uint ps = 0;
		Buffer<Vector4> pi(size * pdim);
		uint vs = 0;
		Buffer<double> vi(size * fields);
		for(uint k=0; k<elems; k++) {
			if(j[k] != qj[i]) continue;
			const Vector4 *pk = &p[pdim * k];
			for(uint l=0; l<pdim; l++) pi[ps++] = pk[l];
			setVector(vector(pk, b), vi, vs);
		}
		q[i].first = qj[i];
		createMultiQuadrature(size, pi, vi, q[i].second);
	}
	return q;
}
void initZero(double &v) { v = 0.0; }
void initZero(Vector4 &v) { v = Vector4(0,0,0,0); }
void initZero(TwoVector4 &v) { v = TwoVector4(0,0,0,0,0,0); }
void initZero(ThreeVector4 &v) { v = ThreeVector4(0,0,0,0); }
void initZero(FourVector4 &v) { v = FourVector4(0); }
template<typename V, typename B, typename R> Buffer<double> &MeshIntegrator::createVector(const Buffer<uint> &j, const Buffer<V> &v, R vector(const V &, const B &), B base(const PartMesh &, const uint), Buffer<double> &q) const {
	R res;
	initZero(res);
	for(uint i=0; i<j.size(); i++) res += vector(v[i], base(m_mesh, j[i])); 
	return createSingleVector(res, q);
}
Buffer<double> &MeshIntegrator::getVector(const uint i, Buffer<double> &q) const {
	Buffer<uint> j;
	q.clear();
	switch (m_grade) {
	case fg_prim0: {
		switch (m_lowdim) {
		case 0: return createSingleVector(1.0, q);
		default: return q;
		}
	}
	case fg_dual0: {
		switch (m_highdim) {
		case 0: return createSingleVector(1.0, q);
		case 1: {
			Buffer<Vector4> v;
			m_mesh.getNodeEdgeVectors(i, j, v);
			removeExternalVectors(m_mesh.getEdgeLocals(), j, v);
			return createVector(j, v, vectorNodeEdge, unitEdge, q);
		}
		case 2: {
			Buffer<TwoVector4> v;
			m_mesh.getNodeFaceVectors(i, j, v);
			removeExternalVectors(m_mesh.getFaceLocals(), j, v);
			return createVector(j, v, vectorNodeFace, unitFace, q);
		}
		case 3: {
			Buffer<ThreeVector4> v;
			m_mesh.getNodeBodyVectors(i, j, v);
			removeExternalVectors(m_mesh.getBodyLocals(), j, v);
			return createVector(j, v, vectorNodeBody, unitBody, q);
		}
		case 4: {
			Buffer<FourVector4> v;
			m_mesh.getNodeQuadVectors(i, j, v);
			removeExternalVectors(m_mesh.getQuadLocals(), j, v);
			return createVector(j, v, vectorNodeQuad, unitQuad, q);
		}
		default: return q;
		}
	}
	case fg_prim1: {
		switch (m_lowdim) {
		case 0: return createSingleVector(m_mesh.getEdgeVector(i), q);
		case 1: return createSingleVector(unitEdge(m_mesh, i), q);
		default: return q;
		}
	}
	case fg_dual1: {
		switch (m_highdim) {
		case 1: return createSingleVector(unitEdge(m_mesh, i), q);
		case 2: {
			Buffer<Vector4> v;
			m_mesh.getEdgeFaceVectors(i, j, v);
			removeExternalVectors(m_mesh.getFaceLocals(), j, v);
			return createVector(j, v, vectorEdgeFace, unitFace, q);
		}
		case 3: {
			Buffer<TwoVector4> v;
			m_mesh.getEdgeBodyVectors(i, j, v);
			removeExternalVectors(m_mesh.getBodyLocals(), j, v);
			return createVector(j, v, vectorEdgeBody, unitBody, q);
		}
		case 4: {
			Buffer<ThreeVector4> v;
			m_mesh.getEdgeQuadVectors(i, j, v);
			removeExternalVectors(m_mesh.getQuadLocals(), j, v);
			return createVector(j, v, vectorEdgeQuad, unitQuad, q);
		}
		default: return q;
		}
	}
	case fg_prim2: {
		switch (m_lowdim) {
		case 0: return createSingleVector(m_mesh.getFaceVector(i), q);
		case 1: {
			Buffer<Vector4> v;
			m_mesh.getFaceEdgeVectors(i, j, v);
			return createVector(j, v, vectorFaceEdge, unitEdge, q);
		}
		case 2: return createSingleVector(unitFace(m_mesh, i), q);
		default: return q;
		}
	}
	case fg_dual2: {
		switch (m_highdim) {
		case 2: return createSingleVector(unitFace(m_mesh, i), q);
		case 3: {
			Buffer<Vector4> v;
			m_mesh.getFaceBodyVectors(i, j, v);
			removeExternalVectors(m_mesh.getBodyLocals(), j, v);
			return createVector(j, v, vectorFaceBody, unitBody, q);
		}
		case 4: {
			Buffer<TwoVector4> v;
			m_mesh.getFaceQuadVectors(i, j, v);
			removeExternalVectors(m_mesh.getQuadLocals(), j, v);
			return createVector(j, v, vectorFaceQuad, unitQuad, q);
		}
		default: return q;
		}
	}
	case fg_prim3: {
		switch (m_lowdim) {
		case 0: return createSingleVector(m_mesh.getBodyVector(i), q);
		case 1: {
			Buffer<TwoVector4> v;
			m_mesh.getBodyEdgeVectors(i, j, v);
			return createVector(j, v, vectorBodyEdge, unitEdge, q);
		}
		case 2: {
			Buffer<Vector4> v;
			m_mesh.getBodyFaceVectors(i, j, v);
			return createVector(j, v, vectorBodyFace, unitFace, q);
		}
		case 3: return createSingleVector(unitBody(m_mesh, i), q);
		default: return q;
		}
	}
	case fg_dual3: {
		switch (m_highdim) {
		case 3: return createSingleVector(unitBody(m_mesh, i), q);
		case 4: {
			Buffer<Vector4> v;
			m_mesh.getBodyQuadVectors(i, j, v);
			removeExternalVectors(m_mesh.getQuadLocals(), j, v);
			return createVector(j, v, vectorBodyQuad, unitQuad, q);
		}
		default: return q;
		}
	}
	case fg_prim4: {
		switch (m_lowdim) {
		case 0: return createSingleVector(m_mesh.getQuadVector(i), q);
		case 1: {
			Buffer<ThreeVector4> v;
			m_mesh.getQuadEdgeVectors(i, j, v);
			return createVector(j, v, vectorQuadEdge, unitEdge, q);
		}
		case 2: {
			Buffer<TwoVector4> v;
			m_mesh.getQuadFaceVectors(i, j, v);
			return createVector(j, v, vectorQuadFace, unitFace, q);
		}
		case 3: {
			Buffer<Vector4> v;
			m_mesh.getQuadBodyVectors(i, j, v);
			return createVector(j, v, vectorQuadBody, unitBody, q);
		}
		case 4: return createSingleVector(unitQuad(m_mesh, i), q);
		default: return q;
		}
	}
	case fg_dual4: {
		switch (m_highdim) {
		case 4: return createSingleVector(unitQuad(m_mesh, i), q);
		default: return q;
		}
	}
	default: return q;
	}
}
template<typename V, typename B, typename R> Buffer< pair<uint, Buffer<double> > > &MeshIntegrator::createBaseVector(const Buffer<uint> &j, const Buffer<V> &v, R vector(const V &, const B &), const B &b, Buffer< pair<uint, Buffer<double> > > &q) const {
	q.resize(j.size());
	for(uint i=0; i<q.size(); i++) {
		q[i].first = j[i];
		createSingleVector(vector(v[i], b), q[i].second); 
	}
	return q;
}
Buffer< pair<uint, Buffer<double> > > &MeshIntegrator::getBaseVector(const uint i, Buffer< pair<uint, Buffer<double> > > &q) const {
	Buffer<uint> j;
	q.clear();
	switch (m_grade) {
	case fg_prim0: {
		switch (m_lowdim) {
		case 0: {
			if(i >= m_mesh.getNodeLocals()) return q;
			return createSingleBaseVector(i, 1.0, q);
		}
		default: return q;
		}
	}
	case fg_dual0: {
		switch (m_highdim) {
		case 0: return createSingleBaseVector(i, 1.0, q);
		case 1: {
			Buffer<Vector4> v;
			m_mesh.getEdgeNodeVectors(i, j, v);
			return createBaseVector(j, v, vectorNodeEdge, unitEdge(m_mesh, i), q);
		}
		case 2: {
			Buffer<TwoVector4> v;
			m_mesh.getFaceNodeVectors(i, j, v);
			return createBaseVector(j, v, vectorNodeFace, unitFace(m_mesh, i), q);
		}
		case 3: {
			Buffer<ThreeVector4> v;
			m_mesh.getBodyNodeVectors(i, j, v);
			return createBaseVector(j, v, vectorNodeBody, unitBody(m_mesh, i), q);
		}
		case 4: {
			Buffer<FourVector4> v;
			m_mesh.getQuadNodeVectors(i, j, v);
			return createBaseVector(j, v, vectorNodeQuad, unitQuad(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_prim1: {
		switch (m_lowdim) {
		case 0: {
			Buffer<Vector4> v;
			m_mesh.getNodeEdgeVectors(i, j, v);
			removeExternalVectors(m_mesh.getEdgeLocals(), j, v);
			return createBaseVector(j, v, vectorEdgeNode, 1.0, q);
		}
		case 1: {
			if(i >= m_mesh.getEdgeLocals()) return q;
			return createSingleBaseVector(i, unitEdge(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_dual1: {
		switch (m_highdim) {
		case 1: return createSingleBaseVector(i, unitEdge(m_mesh, i), q);
		case 2: {
			Buffer<Vector4> v;
			m_mesh.getFaceEdgeVectors(i, j, v);
			return createBaseVector(j, v, vectorEdgeFace, unitFace(m_mesh, i), q);
		}
		case 3: {
			Buffer<TwoVector4> v;
			m_mesh.getBodyEdgeVectors(i, j, v);
			return createBaseVector(j, v, vectorEdgeBody, unitBody(m_mesh, i), q);
		}
		case 4: {
			Buffer<ThreeVector4> v;
			m_mesh.getQuadEdgeVectors(i, j, v);
			return createBaseVector(j, v, vectorEdgeQuad, unitQuad(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_prim2: {
		switch (m_lowdim) {
		case 0: {
			Buffer<TwoVector4> v;
			m_mesh.getNodeFaceVectors(i, j, v);
			removeExternalVectors(m_mesh.getFaceLocals(), j, v);
			return createBaseVector(j, v, vectorFaceNode, 1.0, q);
		}
		case 1: {
			Buffer<Vector4> v;
			m_mesh.getEdgeFaceVectors(i, j, v);
			removeExternalVectors(m_mesh.getFaceLocals(), j, v);
			return createBaseVector(j, v, vectorFaceEdge, unitEdge(m_mesh, i), q);
		}
		case 2: {
			if(i >= m_mesh.getFaceLocals()) return q;
			return createSingleBaseVector(i, unitFace(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_dual2: {
		switch (m_highdim) {
		case 2: return createSingleBaseVector(i, unitFace(m_mesh, i), q);
		case 3: {
			Buffer<Vector4> v;
			m_mesh.getBodyFaceVectors(i, j, v);
			return createBaseVector(j, v, vectorFaceBody, unitBody(m_mesh, i), q);
		}
		case 4: {
			Buffer<TwoVector4> v;
			m_mesh.getQuadFaceVectors(i, j, v);
			return createBaseVector(j, v, vectorFaceQuad, unitQuad(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_prim3: {
		switch (m_lowdim) {
		case 0: {
			Buffer<ThreeVector4> v;
			m_mesh.getNodeBodyVectors(i, j, v);
			removeExternalVectors(m_mesh.getBodyLocals(), j, v);
			return createBaseVector(j, v, vectorBodyNode, 1.0, q);
		}
		case 1: {
			Buffer<TwoVector4> v;
			m_mesh.getEdgeBodyVectors(i, j, v);
			removeExternalVectors(m_mesh.getBodyLocals(), j, v);
			return createBaseVector(j, v, vectorBodyEdge, unitEdge(m_mesh, i), q);
		}
		case 2: {
			Buffer<Vector4> v;
			m_mesh.getFaceBodyVectors(i, j, v);
			removeExternalVectors(m_mesh.getBodyLocals(), j, v);
			return createBaseVector(j, v, vectorBodyFace, unitFace(m_mesh, i), q);
		}
		case 3: {
			if(i >= m_mesh.getBodyLocals()) return q;
			return createSingleBaseVector(i, unitBody(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_dual3: {
		switch (m_highdim) {
		case 3: return createSingleBaseVector(i, unitBody(m_mesh, i), q);
		case 4: {
			Buffer<Vector4> v;
			m_mesh.getQuadBodyVectors(i, j, v);
			return createBaseVector(j, v, vectorBodyQuad, unitQuad(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_prim4: {
		switch (m_lowdim) {
		case 0: {
			Buffer<FourVector4> v;
			m_mesh.getNodeQuadVectors(i, j, v);
			removeExternalVectors(m_mesh.getQuadLocals(), j, v);
			return createBaseVector(j, v, vectorQuadNode, 1.0, q);
		}
		case 1: {
			Buffer<ThreeVector4> v;
			m_mesh.getEdgeQuadVectors(i, j, v);
			removeExternalVectors(m_mesh.getQuadLocals(), j, v);
			return createBaseVector(j, v, vectorQuadEdge, unitEdge(m_mesh, i), q);
		}
		case 2: {
			Buffer<TwoVector4> v;
			m_mesh.getFaceQuadVectors(i, j, v);
			removeExternalVectors(m_mesh.getQuadLocals(), j, v);
			return createBaseVector(j, v, vectorQuadFace, unitFace(m_mesh, i), q);
		}
		case 3: {
			Buffer<Vector4> v;
			m_mesh.getBodyQuadVectors(i, j, v);
			removeExternalVectors(m_mesh.getQuadLocals(), j, v);
			return createBaseVector(j, v, vectorQuadBody, unitBody(m_mesh, i), q);
		}
		case 4: {
			if(i >= m_mesh.getQuadLocals()) return q;
			return createSingleBaseVector(i, unitQuad(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_dual4: {
		switch (m_highdim) {
		case 4: return createSingleBaseVector(i, unitQuad(m_mesh, i), q);
		default: return q;
		}
	}
	default: return q;
	}
}
Buffer<double> &MeshIntegrator::getQuadrature(const uint i, Buffer<double> &q) const {
	Buffer<Vector4> p;
	Buffer<uint> j;
	q.clear();
	switch (m_grade) {
	case fg_prim0: {
		switch (m_lowdim) {
		case 0: return createSingleQuadrature(m_mesh.getNodePosition(i), 1.0, q);
		default: return q;
		}
	}
	case fg_dual0: {
		switch (m_highdim) {
		case 0: {
			if(i >= m_mesh.getNodeLocals()) return q;
			return createSingleQuadrature(m_mesh.getNodePosition(i),1.0, q);
		}
		case 1: {
			m_mesh.getNodeEdgeSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getEdgeLocals(), j, p);
			return createQuadrature(j, p, vectorNodeEdge, unitEdge, q);
		}
		case 2: {
			m_mesh.getNodeFaceSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getFaceLocals(), j, p);
			return createQuadrature(j, p, vectorNodeFace, unitFace, q);
		}
		case 3: {
			m_mesh.getNodeBodySimplices(i, j, p);
			removeExternalSimplices(m_mesh.getBodyLocals(), j, p);
			return createQuadrature(j, p, vectorNodeBody, unitBody, q);
		}
		case 4: {
			m_mesh.getNodeQuadSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getQuadLocals(), j, p);
			return createQuadrature(j, p, vectorNodeQuad, unitQuad, q);
		}
		default: return q;
		}
	}
	case fg_prim1: {
		switch (m_lowdim) {
		case 0: {
			m_mesh.getEdgeSimplices(i, p);
			return createEntryQuadrature(2, p, vectorEdge, q);
		}
		case 1: return createSingleQuadrature(m_mesh.getEdgePosition(i), unitEdge(m_mesh, i), q);
		default: return q;
		}
	}
	case fg_dual1: {
		switch (m_highdim) {
		case 1: {
			if(i >= m_mesh.getEdgeLocals()) return q;
			return createSingleQuadrature(m_mesh.getEdgePosition(i), unitEdge(m_mesh, i), q);
		}
		case 2: {
			m_mesh.getEdgeFaceSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getEdgeLocals(), j, p);
			return createQuadrature(j, p, vectorEdgeFace, unitFace, q);
		}
		case 3: {
			m_mesh.getEdgeBodySimplices(i, j, p);
			removeExternalSimplices(m_mesh.getBodyLocals(), j, p);
			return createQuadrature(j, p, vectorEdgeBody, unitBody, q);
		}
		case 4: {
			m_mesh.getEdgeQuadSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getQuadLocals(), j, p);
			return createQuadrature(j, p, vectorEdgeQuad, unitQuad, q);
		}
		default: return q;
		}
	}
	case fg_prim2: {
		switch (m_lowdim) {
		case 0: {
			m_mesh.getFaceSimplices(i, p);
			return createEntryQuadrature(3, p, vectorFace, q);
		}
		case 1: {
			m_mesh.getFaceEdgeSimplices(i, j, p);
			return createQuadrature(j, p, vectorFaceEdge, unitEdge, q);
		}
		case 2: return createSingleQuadrature(m_mesh.getFacePosition(i), unitFace(m_mesh, i), q);
		default: return q;
		}
	}
	case fg_dual2: {
		switch (m_highdim) {
		case 2: {
			if(i >= m_mesh.getFaceLocals()) return q;
			return createSingleQuadrature(m_mesh.getFacePosition(i), unitFace(m_mesh, i), q);
		}
		case 3: {
			m_mesh.getFaceBodySimplices(i, j, p);
			removeExternalSimplices(m_mesh.getBodyLocals(), j, p);
			return createQuadrature(j, p, vectorFaceBody, unitBody, q);
		}
		case 4: {
			m_mesh.getFaceQuadSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getQuadLocals(), j, p);
			return createQuadrature(j, p, vectorFaceQuad, unitQuad, q);
		}
		default: return q;
		}
	}
	case fg_prim3: {
		switch (m_lowdim) {
		case 0: {
			m_mesh.getBodySimplices(i, p);
			return createEntryQuadrature(4, p, vectorBody, q);
		}
		case 1: {
			m_mesh.getBodyEdgeSimplices(i, j, p);
			return createQuadrature(j, p, vectorBodyEdge, unitEdge, q);
		}
		case 2: {
			m_mesh.getBodyFaceSimplices(i, j, p);
			return createQuadrature(j, p, vectorBodyFace, unitFace, q);
		}
		case 3: return createSingleQuadrature(m_mesh.getBodyPosition(i), unitBody(m_mesh, i), q);
		default: return q;
		}
	}
	case fg_dual3: {
		switch (m_highdim) {
		case 3: {
			if(i >= m_mesh.getBodyLocals()) return q;
			return createSingleQuadrature(m_mesh.getBodyPosition(i), unitBody(m_mesh, i), q);
		}
		case 4: {
			m_mesh.getBodyQuadSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getQuadLocals(), j, p);
			return createQuadrature(j, p, vectorFaceQuad, unitQuad, q);
		}
		default: return q;
		}
	}
	case fg_prim4: {
		switch (m_lowdim) {
		case 0: {
			m_mesh.getQuadSimplices(i, p);
			return createEntryQuadrature(5, p, vectorQuad, q);
		}
		case 1: {
			m_mesh.getQuadEdgeSimplices(i, j, p);
			return createQuadrature(j, p, vectorQuadEdge, unitEdge, q);
		}
		case 2: {
			m_mesh.getQuadFaceSimplices(i, j, p);
			return createQuadrature(j, p, vectorQuadFace, unitFace, q);
		}
		case 3: {
			m_mesh.getQuadBodySimplices(i, j, p);
			return createQuadrature(j, p, vectorQuadBody, unitBody, q);
		}
		case 4: return createSingleQuadrature(m_mesh.getQuadPosition(i), unitQuad(m_mesh, i), q);
		default: return q;
		}
	}
	case fg_dual4: {
		switch (m_highdim) {
		case 4: {
			if(i >= m_mesh.getQuadLocals()) return q;
			return createSingleQuadrature(m_mesh.getQuadPosition(i), unitQuad(m_mesh, i), q);
		}
		default: return q;
		}
	}
	default: return q;
	}
}

Buffer< pair<uint, Buffer<double> > > &MeshIntegrator::getBaseQuadrature(const uint i, Buffer< pair<uint, Buffer<double> > > &q) const {
	Buffer<Vector4> p;
	Buffer<uint> j;
	q.clear();
	switch (m_grade) {
	case fg_prim0: {
		switch (m_lowdim) {
		case 0: {
			if(i >= m_mesh.getNodeLocals()) return q;
			return createSingleBaseQuadrature(i, m_mesh.getNodePosition(i), 1.0, q);
		}
		default: return q;
		}
	}
	case fg_dual0: {
		switch (m_highdim) {
		case 0: return createSingleBaseQuadrature(i, m_mesh.getNodePosition(i),1.0, q);
		case 1: {
			m_mesh.getEdgeNodeSimplices(i, j, p);
			return createBaseQuadrature(j, p, vectorNodeEdge, unitEdge(m_mesh, i), q);
		}
		case 2: {
			m_mesh.getFaceNodeSimplices(i, j, p);
			return createBaseQuadrature(j, p, vectorNodeFace, unitFace(m_mesh, i), q);
		}
		case 3: {
			m_mesh.getBodyNodeSimplices(i, j, p);
			return createBaseQuadrature(j, p, vectorNodeBody, unitBody(m_mesh, i), q);
		}
		case 4: {
			m_mesh.getQuadNodeSimplices(i, j, p);
			return createBaseQuadrature(j, p, vectorNodeQuad, unitQuad(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_prim1: {
		switch (m_lowdim) {
		case 0: {
			m_mesh.getNodeEdgeSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getEdgeLocals(), j, p);
			return createBaseQuadrature(j, p, vectorEdgeNode, 1.0, q);
		}
		case 1: {
			if(i >= m_mesh.getEdgeLocals()) return q;
			return createSingleBaseQuadrature(i, m_mesh.getEdgePosition(i), unitEdge(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_dual1: {
		switch (m_highdim) {
		case 1: return createSingleBaseQuadrature(i, m_mesh.getEdgePosition(i), unitEdge(m_mesh, i), q);
		case 2: {
			m_mesh.getFaceEdgeSimplices(i, j, p);
			createBaseQuadrature(j, p, vectorEdgeFace, unitFace(m_mesh, i), q);
			return q;
		}
		case 3: {
			m_mesh.getBodyEdgeSimplices(i, j, p);
			return createBaseQuadrature(j, p, vectorEdgeBody, unitBody(m_mesh, i), q);
		}
		case 4: {
			m_mesh.getQuadEdgeSimplices(i, j, p);
			return createBaseQuadrature(j, p, vectorEdgeQuad, unitQuad(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_prim2: {
		switch (m_lowdim) {
		case 0: {
			m_mesh.getNodeFaceSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getFaceLocals(), j, p);
			return createBaseQuadrature(j, p, vectorFaceNode, 1.0, q);
		}
		case 1: {
			m_mesh.getEdgeFaceSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getFaceLocals(), j, p);
			return createBaseQuadrature(j, p, vectorFaceEdge, unitEdge(m_mesh, i), q);
		}
		case 2: {
			if(i >= m_mesh.getFaceLocals()) return q;
			return createSingleBaseQuadrature(i, m_mesh.getFacePosition(i), unitFace(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_dual2: {
		switch (m_highdim) {
		case 2: return createSingleBaseQuadrature(i, m_mesh.getFacePosition(i), unitFace(m_mesh, i), q);
		case 3: {
			m_mesh.getBodyFaceSimplices(i, j, p);
			return createBaseQuadrature(j, p, vectorFaceBody, unitBody(m_mesh, i), q);
		}
		case 4: {
			m_mesh.getQuadFaceSimplices(i, j, p);
			return createBaseQuadrature(j, p, vectorFaceQuad, unitQuad(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_prim3: {
		switch (m_lowdim) {
		case 0: {
			m_mesh.getNodeBodySimplices(i, j, p);
			removeExternalSimplices(m_mesh.getBodyLocals(), j, p);
			return createBaseQuadrature(j, p, vectorBodyNode, 1.0, q);
		}
		case 1: {
			m_mesh.getEdgeBodySimplices(i, j, p);
			removeExternalSimplices(m_mesh.getBodyLocals(), j, p);
			return createBaseQuadrature(j, p, vectorBodyEdge, unitEdge(m_mesh, i), q);
		}
		case 2: {
			m_mesh.getFaceBodySimplices(i, j, p);
			removeExternalSimplices(m_mesh.getBodyLocals(), j, p);
			return createBaseQuadrature(j, p, vectorBodyFace, unitFace(m_mesh, i), q);
		}
		case 3: {
			if(i >= m_mesh.getBodyLocals()) return q;
			return createSingleBaseQuadrature(i, m_mesh.getBodyPosition(i), unitBody(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_dual3: {
		switch (m_highdim) {
		case 3: return createSingleBaseQuadrature(i, m_mesh.getBodyPosition(i), unitBody(m_mesh, i), q);
		case 4: {
			m_mesh.getQuadFaceSimplices(i, j, p);
			return createBaseQuadrature(j, p, vectorFaceQuad, unitQuad(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_prim4: {
		switch (m_lowdim) {
		case 0: {
			m_mesh.getNodeQuadSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getQuadLocals(), j, p);
			return createBaseQuadrature(j, p, vectorQuadNode, 1.0, q);
		}
		case 1: {
			m_mesh.getEdgeQuadSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getQuadLocals(), j, p);
			return createBaseQuadrature(j, p, vectorQuadEdge, unitEdge(m_mesh, i), q);
		}
		case 2: {
			m_mesh.getFaceQuadSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getQuadLocals(), j, p);
			return createBaseQuadrature(j, p, vectorQuadFace, unitFace(m_mesh, i), q);
		}
		case 3: {
			m_mesh.getBodyQuadSimplices(i, j, p);
			removeExternalSimplices(m_mesh.getQuadLocals(), j, p);
			return createBaseQuadrature(j, p, vectorQuadBody, unitBody(m_mesh, i), q);
		}
		case 4: {
			if(i >= m_mesh.getQuadLocals()) return q;
			return createSingleBaseQuadrature(i, m_mesh.getQuadPosition(i), unitQuad(m_mesh, i), q);
		}
		default: return q;
		}
	}
	case fg_dual4: {
		switch (m_highdim) {
		case 4: return createSingleBaseQuadrature(i, m_mesh.getQuadPosition(i), unitQuad(m_mesh, i), q);
		default: return q;
		}
	}
	default: return q;
	}
}

Buffer<double> &MeshIntegrator::invertVector(Buffer<double> &q) const {
	uint i;
	// compute square of the vector
	double sq = 0.0;
	for(i=0; i<q.size(); i++) sq += q[i] * q[i];
	if(sq == 0.0) return q;
	// rescale vector 
	const double div = 1.0 / sq;
	for(i=0; i<q.size(); i++) q[i] *= div;
	return q;
}
Buffer<double> &MeshIntegrator::invertQuadrature(const Vector4 &p, Buffer<double> &q) const {
	uint i, j;
	const uint dim = m_mesh.getDimension();
	const uint fields = getFields();
	const uint size = dim + fields;
	// move positions and sum up vectors
	Buffer<double> v(fields, 0.0);
	for(i=0; i<q.size(); i+=size) {
		q[i] -= p.x;
		if(dim >= 2) q[i+1] -= p.y;
		if(dim >= 3) q[i+2] -= p.z;
		if(dim >= 4) q[i+3] -= p.t;
		const double *qi = &q[i + dim];
		for(j=0; j<fields; j++) v[j] += qi[j];
	}
	// compute square of summed vector
	double sq = 0.0;
	for(j=0; j<fields; j++) sq += v[j] * v[j];
	if(sq == 0.0) return q;
	// rescale vectors
	const double div = 1.0 / sq;
	for(i=dim; i<q.size(); i+=size) {
		double *qi = &q[i];
		for(j=0; j<fields; j++) qi[j] *= div;
	}
	return q;
}
template<typename V> Buffer<double> &MeshIntegrator::createSingleVector(const V &v, Buffer<double> &q) const {
	uint qs = 0;
	q.resize(getFields());
	setVector(v, q, qs);
	return q;
}
template<typename V> Buffer< pair<uint, Buffer<double> > > &MeshIntegrator::createSingleBaseVector(const uint i, const V &v, Buffer< pair<uint, Buffer<double> > > &q) const {
	q.resize(1);
	q[0].first = i;
	createSingleVector(v, q[0].second);
	return q;
}
template<typename V> Buffer<double> &MeshIntegrator::createSingleQuadrature(const Vector4 &p, const V &v, Buffer<double> &q) const {
	uint qs = 0;
	q.resize(m_mesh.getDimension() + getFields());
	setVector(p, q, qs);
	setVector(v, q, qs);
	return q;
}
template<typename V> Buffer< pair<uint, Buffer<double> > > &MeshIntegrator::createSingleBaseQuadrature(const uint i, const Vector4 &p, const V &v, Buffer< pair<uint, Buffer<double> > > &q) const {
	q.resize(1);
	q[0].first = i;
	createSingleQuadrature(p, v, q[0].second);
	return q;
}

Buffer<double> &MeshIntegrator::createMultiQuadrature(const uint elems, const Buffer<Vector4> &p, const Buffer<double> &v, Buffer<double> &q) const {
	uint i, j, n;
	uint qs = 0;
	const uint pdim = p.size() / elems;
	const uint vdim = v.size() / elems;
	const uint qdim = m_mesh.getDimension() + vdim;
	switch (m_num) {
	case 0: { // compress simplices into one average point
		Buffer<double> sumv(vdim, 0.0);
		for(i=0; i<elems; i++) {
			const double *vi = &v[vdim * i];
			for(n=0; n<vdim; n++) sumv[n] += vi[n];
		}
		Vector4 sump(0,0,0,0);
		for(i=0; i<elems; i++) {
			double dot = 0.0;
			const double *vi = &v[vdim * i];
			for(n=0; n<vdim; n++) dot += sumv[n] * vi[n];
			const Vector4 *pi = &p[pdim * i];
			for(n=0; n<pdim; n++) sump += dot * pi[n];
		}
		double sumsq = 0.0;
		for(n=0; n<vdim; n++) sumsq += sumv[n] * sumv[n];
		sump /= double(pdim) * sumsq;
		q.resize(qdim);
		setVector(sump, q, qs);
		for(n=0; n<vdim; n++) q[qs++] = sumv[n];
		return q;
	}
	case 1: { // use one average point per simplex
		q.resize(qdim * elems);
		for(i=0; i<elems; i++) {
			const Vector4 *pi = &p[pdim * i];
			Vector4 sump(0,0,0,0);
			for(n=0; n<pdim; n++) sump += pi[n];
			sump /= double(pdim);
			setVector(sump, q, qs);
			const double *vi = &v[vdim * i];
			for(n=0; n<vdim; n++) q[qs++] = vi[n];
		}
		return q;
	}
	default: { // use several average points for the quadrature
		const double a = 1.0 / sqrt((m_num - 1) * (m_num + pdim - 1));
		const double b = (1.0 / a - double(m_num - 1)) / double(pdim);
		uint num = 1;
		for(n=1; n<pdim; n++) num = (num * (m_num + n - 1)) / n;
		const double c = 1.0 / double(num);
		Buffer<double> cvi(vdim);
		Buffer<Vector4> di(pdim - 1);
		q.resize(qdim * elems * num);
		for(i=0; i<elems; i++) {
			const Vector4 *pi = &p[pdim * i];
			Vector4 p0 = pi[0];
			for(n=1; n<pdim; n++) {
				di[n-1] = a * (pi[n] - pi[0]);
				p0 += b * di[n-1];
			}
			const double *vi = &v[vdim * i];
			for(n=0; n<vdim; n++) cvi[n] = c * vi[n];
			if(pdim == 1) {
				setVector(p0, q, qs);
				for(n=0; n<vdim; n++) q[qs++] = cvi[n];
				continue;
			}
			for(j=0; j<m_num; j++) {
				const Vector4 p1 = p0 + j * di[0];
				if(pdim == 2) {
					setVector(p1, q, qs);
					for(n=0; n<vdim; n++) q[qs++] = cvi[n];
					continue;
				}
				for(uint k=j; k<m_num; k++) {
					const Vector4 p2 = p1 + (k - j) * di[1];
					if(pdim == 3) {
						setVector(p2, q, qs);
						for(n=0; n<vdim; n++) q[qs++] = cvi[n];
						continue;
					}
					for(uint l=k; l<m_num; l++) {
						const Vector4 p3 = p2 + (l - k) * di[2];
						if(pdim == 4) {
							setVector(p3, q, qs);
							for(n=0; n<vdim; n++) q[qs++] = cvi[n];
							continue;
						}
						for(uint m=l; m<m_num; m++) {
							setVector(p3 + (m - l) * di[3], q, qs);
							for(n=0; n<vdim; n++) q[qs++] = cvi[n];
						}
					}
				}
			}
		}
		return q;
	}
	}
}

void MeshIntegrator::setVector(const double &v, Buffer<double> &q, uint &qs) const {
	q[qs++] = v;
}
void MeshIntegrator::setVector(const Vector4 &v, Buffer<double> &q, uint &qs) const {
	const uint dim = m_mesh.getDimension();
	q[qs++] = v.x;
	if(dim == 1) return;
	q[qs++] = v.y;
	if(dim == 2) return;
	q[qs++] = v.z;
	if(dim == 3) return;
	q[qs++] = v.t;
}
void MeshIntegrator::setVector(const TwoVector4 &v, Buffer<double> &q, uint &qs) const {
	const uint dim = m_mesh.getDimension();
	q[qs++] = v.xy;
	if(dim == 2) return;
	q[qs++] = v.xz;
	q[qs++] = v.yz;
	if(dim == 3) return;
	q[qs++] = v.xt;
	q[qs++] = v.yt;
	q[qs++] = v.zt;
}
void MeshIntegrator::setVector(const ThreeVector4 &v, Buffer<double> &q, uint &qs) const {
	q[qs++] = v.xyz;
	if(m_mesh.getDimension() == 3) return;
	q[qs++] = v.xyt;
	q[qs++] = v.xzt;
	q[qs++] = v.yzt;
}
void MeshIntegrator::setVector(const FourVector4 &v, Buffer<double> &q, uint &qs) const {
	q[qs++] = v.xyzt;
}

