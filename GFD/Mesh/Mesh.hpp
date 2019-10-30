/*
Mesh is a class for partitioning a 1--4 dimensional domain with polyhedral cells.
We use the following naming: Node = 0-cell, Edge = 1-cell, Face = 2-cell, Body = 3-cell, and Quad = 4-cell.
Cells are constructed recursively and are linked with their boundary cells and parent cells.
Each cell can be assigned with a flag (unsigned int)
*/

#ifndef _MESH_HPP_INCLUDED_
#define _MESH_HPP_INCLUDED_

#include "../Types/Types.hpp"
#include "../Types/Buffer.hpp"
#include "../Types/Text.hpp"
#include "../Types/UintSet.hpp"

namespace gfd
{

struct Node
{
	Buffer<uint> e;
};

struct Edge
{
	Buffer<uint> n;
	Buffer<uint> f;
};

struct Face
{
	Buffer<uint> e;
	Buffer<uint> b;
};

struct Body
{
	Buffer<uint> f;
	Buffer<uint> q;
};

struct Quad
{
	Buffer<uint> b;
};

class Mesh
{
public:
	Mesh(const uint dim = 4);
	virtual ~Mesh() { clear(); }
	virtual void clear();
	void swap(Mesh &mesh);

	// load and save mesh
	//bool loadMesh(const std::string &path);
	//bool saveMesh(const std::string &path) const;
	bool loadJRMesh(const std::string &path);
	bool saveJRMesh(const std::string &path) const;

	// statistics
	void writeStatistics(Text &text, const UintSet &flag = UINTSETALL) const;

	// vector space dimension
	uint getDimension() const { return m_dim; }

	// number of cells
	uint getNodeSize() const { return m_nsize; }
	uint getEdgeSize() const { return m_esize; }
	uint getFaceSize() const { return m_fsize; }
	uint getBodySize() const { return m_bsize; }
	uint getQuadSize() const { return m_qsize; }

	// intersections
	uint getEdgeIntersection(const uint e0, const uint e1) const { return m_e[e0].n.getFirstIntersection(m_e[e1].n, NONE); }
	uint getFaceIntersection(const uint f0, const uint f1) const { return m_f[f0].e.getFirstIntersection(m_f[f1].e, NONE); }
	uint getBodyIntersection(const uint b0, const uint b1) const { return m_b[b0].f.getFirstIntersection(m_b[b1].f, NONE); }
	uint getQuadIntersection(const uint q0, const uint q1) const { return m_q[q0].b.getFirstIntersection(m_q[q1].b, NONE); }

	// neighbors
	const Buffer<uint> &getNodeEdges(const uint n) const { return m_n[n].e; }
	Buffer<uint> getNodeFaces(const uint n) const;
	Buffer<uint> getNodeBodies(const uint n) const;
	Buffer<uint> getNodeQuads(const uint n) const;
	const Buffer<uint> &getEdgeNodes(const uint e) const { return m_e[e].n; }
	const Buffer<uint> &getEdgeFaces(const uint e) const { return m_e[e].f; }
	Buffer<uint> getEdgeBodies(const uint e) const;
	Buffer<uint> getEdgeQuads(const uint e) const;
	Buffer<uint> getFaceNodes(const uint f) const;
	const Buffer<uint> &getFaceEdges(const uint f) const { return m_f[f].e; }
	const Buffer<uint> &getFaceBodies(const uint f) const { return m_f[f].b; }
	Buffer<uint> getFaceQuads(const uint f) const;
	Buffer<uint> getBodyNodes(const uint b) const;
	Buffer<uint> getBodyEdges(const uint b) const;
	const Buffer<uint> &getBodyFaces(const uint b) const { return m_b[b].f; }
	const Buffer<uint> &getBodyQuads(const uint b) const { return m_b[b].q; }
	Buffer<uint> getQuadNodes(const uint q) const;
	Buffer<uint> getQuadEdges(const uint q) const;
	Buffer<uint> getQuadFaces(const uint q) const;
	const Buffer<uint> &getQuadBodies(const uint q) const { return m_q[q].b; }

	// get single neighbor
	uint getNodeAnyEdge(const uint n) const { return m_n[n].e[0]; }
	uint getNodeAnyFace(const uint n) const { return m_e[m_n[n].e[0]].f[0]; }
	uint getNodeAnyBody(const uint n) const { return m_f[m_e[m_n[n].e[0]].f[0]].b[0]; }
	uint getNodeAnyQuad(const uint n) const { return m_b[m_f[m_e[m_n[n].e[0]].f[0]].b[0]].q[0]; }
	uint getEdgeAnyNode(const uint e) const { return m_e[e].n[0]; }
	uint getEdgeAnyFace(const uint e) const { return m_e[e].f[0]; }
	uint getEdgeAnyBody(const uint e) const { return m_f[m_e[e].f[0]].b[0]; }
	uint getEdgeAnyQuad(const uint e) const { return m_b[m_f[m_e[e].f[0]].b[0]].q[0]; }
	uint getFaceAnyNode(const uint f) const { return m_e[m_f[f].e[0]].n[0]; }
	uint getFaceAnyEdge(const uint f) const { return m_f[f].e[0]; }
	uint getFaceAnyBody(const uint f) const { return m_f[f].b[0]; }
	uint getFaceAnyQuad(const uint f) const { return m_b[m_f[f].b[0]].q[0]; }
	uint getBodyAnyNode(const uint b) const { return m_e[m_f[m_b[b].f[0]].e[0]].n[0]; }
	uint getBodyAnyEdge(const uint b) const { return m_f[m_b[b].f[0]].e[0]; }
	uint getBodyAnyFace(const uint b) const { return m_b[b].f[0]; }
	uint getBodyAnyQuad(const uint b) const { return m_b[b].q[0]; }
	uint getQuadAnyNode(const uint q) const { return m_e[m_f[m_b[m_q[q].b[0]].f[0]].e[0]].n[0]; }
	uint getQuadAnyEdge(const uint q) const { return m_f[m_b[m_q[q].b[0]].f[0]].e[0]; }
	uint getQuadAnyFace(const uint q) const { return m_b[m_q[q].b[0]].f[0]; }
	uint getQuadAnyBody(const uint q) const { return m_q[q].b[0]; }

	// more neighbors
	uint getEdgeOtherNode(const uint e, const uint n) const;

	// incidence
	sign getEdgeIncidence(const uint e, const uint n) const;
	sign getFaceIncidence(const uint f, const uint e) const;
	sign getBodyIncidence(const uint b, const uint f) const;
	sign getQuadIncidence(const uint q, const uint b) const;

	// circumcenter positions
	double getNodePosition1(const uint n) const { return m_p[n * m_dim]; }
	Vector2 getNodePosition2(const uint n) const { const uint i = n * m_dim; return Vector2(m_p[i], m_p[i+1]); }
	Vector3 getNodePosition3(const uint n) const { const uint i = n * m_dim; return Vector3(m_p[i], m_p[i+1], m_p[i+2]); }
	Vector4 getNodePosition4(const uint n) const { const uint i = n * m_dim; return Vector4(m_p[i], m_p[i+1], m_p[i+2], m_p[i+3]); }
	Vector4 getNodePosition(const uint n) const;
	double getEdgePosition1(const uint e) const;
	Vector2 getEdgePosition2(const uint e) const;
	Vector3 getEdgePosition3(const uint e) const;
	Vector4 getEdgePosition4(const uint e) const;
	Vector4 getEdgePosition(const uint e) const;
	Vector2 getFacePosition2(const uint f) const;
	Vector3 getFacePosition3(const uint f) const;
	Vector4 getFacePosition4(const uint f) const;
	Vector4 getFacePosition(const uint f) const;
	Vector3 getBodyPosition3(const uint b) const;
	Vector4 getBodyPosition4(const uint b) const;
	Vector4 getBodyPosition(const uint b) const;
	Vector4 getQuadPosition(const uint q) const;

	// average positions
	double getEdgeAverage1(const uint e) const;
	Vector2 getEdgeAverage2(const uint e) const;
	Vector3 getEdgeAverage3(const uint e) const;
	Vector4 getEdgeAverage4(const uint e) const;
	Vector4 getEdgeAverage(const uint e) const;
	Vector2 getFaceAverage2(const uint f) const;
	Vector3 getFaceAverage3(const uint f) const;
	Vector4 getFaceAverage4(const uint f) const;
	Vector4 getFaceAverage(const uint f) const;
	Vector3 getBodyAverage3(const uint b) const;
	Vector4 getBodyAverage4(const uint b) const;
	Vector4 getBodyAverage(const uint b) const;
	Vector4 getQuadAverage(const uint q) const;

	// circumcenter determinants
	double getNodeWeight(const uint n) const { if(n < m_w.size()) return m_w[n]; return 0.0; }
	SymMatrix4 getMetric() const;
	Vector4 getTransformed(const Vector4 &r) const;

	// dual average positions
	double getNodeDualAverage1(const uint n) const;
	Vector2 getNodeDualAverage2(const uint n) const;
	Vector3 getNodeDualAverage3(const uint n) const;
	Vector4 getNodeDualAverage4(const uint n) const;
	Vector4 getNodeDualAverage(const uint n) const;
	Vector2 getEdgeDualAverage2(const uint e) const;
	Vector3 getEdgeDualAverage3(const uint e) const;
	Vector4 getEdgeDualAverage4(const uint e) const;
	Vector4 getEdgeDualAverage(const uint e) const;
	Vector3 getFaceDualAverage3(const uint f) const;
	Vector4 getFaceDualAverage4(const uint f) const;
	Vector4 getFaceDualAverage(const uint f) const;
	Vector4 getBodyDualAverage4(const uint b) const;
	Vector4 getBodyDualAverage(const uint b) const;

	// primal volume vectors
	double getEdgeVector1(const uint e) const;
	Vector2 getEdgeVector2(const uint e) const;
	Vector3 getEdgeVector3(const uint e) const;
	Vector4 getEdgeVector4(const uint e) const;
	Vector4 getEdgeVector(const uint e) const;
	TwoVector2 getFaceVector2(const uint f) const;
	TwoVector3 getFaceVector3(const uint f) const;
	TwoVector4 getFaceVector4(const uint f) const;
	TwoVector4 getFaceVector(const uint f) const;
	ThreeVector3 getBodyVector3(const uint b) const;
	ThreeVector4 getBodyVector4(const uint b) const;
	ThreeVector4 getBodyVector(const uint b) const;
	FourVector4 getQuadVector(const uint q) const;

	// dual volume vectors
	double getNodeDualVector1(const uint n) const;
	TwoVector2 getNodeDualVector2(const uint n) const;
	ThreeVector3 getNodeDualVector3(const uint n) const;
	FourVector4 getNodeDualVector4(const uint n) const;
	FourVector4 getNodeDualVector(const uint n) const;
	double getEdgeDualVector1(const uint e) const;
	Vector2 getEdgeDualVector2(const uint e) const;
	TwoVector3 getEdgeDualVector3(const uint e) const;
	ThreeVector4 getEdgeDualVector4(const uint e) const;
	ThreeVector4 getEdgeDualVector(const uint e) const;
	double getFaceDualVector2(const uint f) const;
	Vector3 getFaceDualVector3(const uint f) const;
	TwoVector4 getFaceDualVector4(const uint f) const;
	TwoVector4 getFaceDualVector(const uint f) const;
	double getBodyDualVector3(const uint b) const;
	Vector4 getBodyDualVector4(const uint b) const;
	Vector4 getBodyDualVector(const uint b) const;
	double getQuadDualVector(const uint q) const;

	// diagonal Hodge terms
	double getNodeHodge(const uint n) const;
	double getNodeHodge(const uint n, const double &metric) const;
	double getEdgeHodge(const uint e) const;
	double getEdgeHodge(const uint e, const SymMatrix4 &metric) const;
	double getFaceHodge(const uint f) const;
	double getFaceHodge(const uint f, const SymTwoMatrix4 &metric) const;
	double getBodyHodge(const uint b) const;
	double getBodyHodge(const uint b, const SymThreeMatrix4 &metric) const;
	double getQuadHodge(const uint q) const;
	double getQuadHodge(const uint q, const SymFourMatrix4 &metric) const;

	// get deviation vector from the cell plane to position p (independent of metric)
	Vector4 getEdgeDeviation(const uint e, const Vector4 &p) const;
	Vector4 getFaceDeviation(const uint f, const Vector4 &p) const;
	Vector4 getBodyDeviation(const uint b, const Vector4 &p) const;

	// get projection of vector d that is orthogonal to the cell (depend on the metric)
	Vector4 getEdgeOrthogonal(const uint e, const Vector4 &d) const;
	Vector4 getFaceOrthogonal(const uint f, const Vector4 &d) const;
	Vector4 getBodyOrthogonal(const uint b, const Vector4 &d) const;

	// find elements
	uint findNode(const Vector4 &p, const double zerolensq = 1e-13, uint curr = 0, const bool assured = true) const;
	uint findEdge(const uint n0, const uint n1) const;
	uint findFace(const Buffer<uint> &e) const;
	uint findBody(const Buffer<uint> &f) const;
	uint findQuad(const Buffer<uint> &b) const;

	// flags
	uint getNodeFlag(const uint i) const { if(i < m_nflag.size()) return m_nflag[i]; return 0; }
	uint getEdgeFlag(const uint i) const { if(i < m_eflag.size()) return m_eflag[i]; return 0; }
	uint getFaceFlag(const uint i) const { if(i < m_fflag.size()) return m_fflag[i]; return 0; }
	uint getBodyFlag(const uint i) const { if(i < m_bflag.size()) return m_bflag[i]; return 0; }
	uint getQuadFlag(const uint i) const { if(i < m_qflag.size()) return m_qflag[i]; return 0; }

	// add and remove cells
	uint addNode(const Vector4 &p);
	uint addEdge(const uint n0, const uint n1);
	uint addFace(const Buffer<uint> &e);
	uint addBody(const Buffer<uint> &f);
	uint addQuad(const Buffer<uint> &b);
	void removeNode(const uint n);
	void removeEdge(const uint e);
	void removeFace(const uint f);
	void removeBody(const uint b);
	void removeQuad(const uint q);

	// transformation and relocation
	void setNodePosition(const uint n, const Vector4 &p);
	void transform(const Matrix4 &mat);
	void move(const Vector4 &vec);
	void setNodeWeight(const uint n, const double w);
	void setMetric(const SymMatrix4 &m);

	// modify flags
	void setNodeFlag(const uint n, const uint flag);
	void setEdgeFlag(const uint e, const uint flag);
	void setFaceFlag(const uint f, const uint flag);
	void setBodyFlag(const uint b, const uint flag);
	void setQuadFlag(const uint q, const uint flag);

	// resizing the element buffers (use this only for optimization, if you know the element sizes in advance)
	void resizeNodeBuffer(const uint size);
	void resizeEdgeBuffer(const uint size);
	void resizeFaceBuffer(const uint size);
	void resizeBodyBuffer(const uint size);
	void resizeQuadBuffer(const uint size);

	void gatherNodeSimplices(const uint n, Buffer<Vector4> &p, uint &ps, Buffer<double> &v, uint &vs) const;
	void gatherEdgeSimplices(const uint e, Buffer<Vector4> &p, uint &ps, Buffer<Vector4> &v, uint &vs) const;
	void gatherFaceSimplices(const uint f, Buffer<Vector4> &p, uint &ps, Buffer<TwoVector4> &v, uint &vs) const;
	void gatherBodySimplices(const uint b, Buffer<Vector4> &p, uint &ps, Buffer<ThreeVector4> &v, uint &vs) const;
	void gatherQuadSimplices(const uint q, Buffer<Vector4> &p, uint &ps, Buffer<FourVector4> &v, uint &vs) const;
	void gatherNodeDualSimplices(const uint n, Buffer<Vector4> &p, uint &ps, Buffer<FourVector4> &v, uint &vs) const;
	void gatherEdgeDualSimplices(const uint e, Buffer<Vector4> &p, uint &ps, Buffer<ThreeVector4> &v, uint &vs) const;
	void gatherFaceDualSimplices(const uint f, Buffer<Vector4> &p, uint &ps, Buffer<TwoVector4> &v, uint &vs) const;
	void gatherBodyDualSimplices(const uint b, Buffer<Vector4> &p, uint &ps, Buffer<Vector4> &v, uint &vs) const;
	void gatherQuadDualSimplices(const uint q, Buffer<Vector4> &p, uint &ps, Buffer<double> &v, uint &vs) const;

/*	Buffer<double> getNodeQuadrature(const uint n) const;
	Buffer<double> getEdgeQuadrature(const uint e, const double h) const;
	Buffer<double> getFaceQuadrature(const uint f, const double h) const;
	Buffer<double> getBodyQuadrature(const uint b, const double h) const;
	Buffer<double> getQuadQuadrature(const uint q, const double h) const;
	Buffer<double> getNodeDualQuadrature(const uint n, const double h) const;
	Buffer<double> getEdgeDualQuadrature(const uint e, const double h) const;
	Buffer<double> getFaceDualQuadrature(const uint f, const double h) const;
	Buffer<double> getBodyDualQuadrature(const uint b, const double h) const;
	Buffer<double> getQuadDualQuadrature(const uint q) const;
*/
protected:

	uint m_dim; // vector space dimension
	Buffer<double> m_p; // vector coordinates for each node

	// mesh elements
	uint m_nsize;
	Buffer<Node> m_n; // nodes
	uint m_esize;
	Buffer<Edge> m_e; // edges
	uint m_fsize;
	Buffer<Face> m_f; // faces
	uint m_bsize;
	Buffer<Body> m_b; // bodies
	uint m_qsize;
	Buffer<Quad> m_q; // quads

	// flags
	Buffer<uint> m_nflag; // node flags (optional)
	Buffer<uint> m_eflag; // edge flags (optional)
	Buffer<uint> m_fflag; // face flags (optional)
	Buffer<uint> m_bflag; // body flags (optional)
	Buffer<uint> m_qflag; // quad flags (optional)

	// circumcenter computation (squared distance of v is v.dot(m * v) + w)
	Buffer<double> m_m; // symmetric matrix to determine dot product (optional)
	Buffer<double> m_w; // node weights (optional)

	double getMetric1() const { return m_m[0]; }
	SymMatrix2 getMetric2() const { return SymMatrix2(m_m[0], m_m[1], m_m[2]); }
	SymMatrix3 getMetric3() const { return SymMatrix3(m_m[0], m_m[1], m_m[2], m_m[3], m_m[4], m_m[5]); }
	SymMatrix4 getMetric4() const { return SymMatrix4(m_m[0], m_m[1], m_m[2], m_m[3], m_m[4], m_m[5], m_m[6], m_m[7], m_m[8], m_m[9]); }
	double getTransformed1(const double &r) const { return (m_m.empty() ? r : getMetric1() * r); }
	Vector2 getTransformed2(const Vector2 &r) const { return (m_m.empty() ? r : getMetric2() * r); }
	Vector3 getTransformed3(const Vector3 &r) const { return (m_m.empty() ? r : getMetric3() * r); }
	Vector4 getTransformed4(const Vector4 &r) const { return (m_m.empty() ? r : getMetric4() * r); }

public:
	void orderFaceEdges(const uint f);
	void orderBodyFaces(const uint b);
	void orderQuadBodies(const uint q);
/*
	void insertQuadrature(const Vector4 &p0, const double w, Buffer<double> &q, uint &qs) const;
	void insertQuadrature(const Vector4 &p0, const Vector4 &p1, const Vector4 &w, const double h, Buffer<double> &q, uint &qs) const;
	void insertQuadrature(const Vector4 &p0, const Vector4 &p1, const Vector4 &p2, const TwoVector4 &w, const double h, Buffer<double> &q, uint &qs) const;
	void insertQuadrature(const Vector4 &p0, const Vector4 &p1, const Vector4 &p2, const Vector4 &p3, const ThreeVector4 &w, const double h, Buffer<double> &q, uint &qs) const;
	void insertQuadrature(const Vector4 &p0, const Vector4 &p1, const Vector4 &p2, const Vector4 &p3, const Vector4 &p4, const FourVector4 &w, const double h, Buffer<double> &q, uint &qs) const;
*/
};

}

#endif //_MESH_HPP_INCLUDED_
