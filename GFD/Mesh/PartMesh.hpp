/**
 * PartMesh is a partitioned Mesh.
 * Author: Jukka Räbinä, University of Jyväskylä, 2019.
 */

#ifndef _PARTMESH_HPP_INCLUDED_
#define _PARTMESH_HPP_INCLUDED_

#include "Mesh.hpp"

using namespace std;

namespace gfd
{

class PartMesh : public Mesh
{
public:
	PartMesh(const uint part, const uint parts, const uint dim = 4);
	virtual ~PartMesh() { clear(); }
	void clear();

	void createPartFromFlags(const Mesh &mesh);
	void createPart(const Mesh &mesh, const Buffer<uint> &npart, 
		const Buffer<uint> &epart = Buffer<uint>(), const Buffer<uint> &fpart = Buffer<uint>(), 
		const Buffer<uint> &bpart = Buffer<uint>(), const Buffer<uint> &qpart = Buffer<uint>());
//	void createPart(const Mesh &mesh, const Buffer< pair<uint,uint> > &extn, const Buffer<uint> &repeat = Buffer<uint>()); 
	void createCombined(Buffer<const PartMesh *> &mesh);

	uint getPart() const { return m_part; }
	uint getNumberOfParts() const { return m_parts; }

	uint getNodeLocals() const { return m_nsize - m_extn.size(); }
	uint getEdgeLocals() const { return m_esize - m_exte.size(); }
	uint getFaceLocals() const { return m_fsize - m_extf.size(); }
	uint getBodyLocals() const { return m_bsize - m_extb.size(); }
	uint getQuadLocals() const { return m_qsize - m_extq.size(); }

	void setNodePart(const uint i, const uint part) { const uint j = i - getNodeLocals(); if(j < m_extn.size()) m_extn[j].first = part; }
	void setEdgePart(const uint i, const uint part) { const uint j = i - getEdgeLocals(); if(j < m_exte.size()) m_exte[j].first = part; }
	void setFacePart(const uint i, const uint part) { const uint j = i - getFaceLocals(); if(j < m_extf.size()) m_extf[j].first = part; }
	void setBodyPart(const uint i, const uint part) { const uint j = i - getBodyLocals(); if(j < m_extb.size()) m_extb[j].first = part; }
	void setQuadPart(const uint i, const uint part) { const uint j = i - getQuadLocals(); if(j < m_extq.size()) m_extq[j].first = part; }

	void setNodeLink(const uint i, const uint link) { const uint j = i - getNodeLocals(); if(j < m_extn.size()) m_extn[j].second = link; }
	void setEdgeLink(const uint i, const uint link) { const uint j = i - getEdgeLocals(); if(j < m_exte.size()) m_exte[j].second = link; }
	void setFaceLink(const uint i, const uint link) { const uint j = i - getFaceLocals(); if(j < m_extf.size()) m_extf[j].second = link; }
	void setBodyLink(const uint i, const uint link) { const uint j = i - getBodyLocals(); if(j < m_extb.size()) m_extb[j].second = link; }
	void setQuadLink(const uint i, const uint link) { const uint j = i - getQuadLocals(); if(j < m_extq.size()) m_extq[j].second = link; }

	uint getNodePart(const uint i) const { const uint j = i - getNodeLocals(); if(j < m_extn.size()) return m_extn[j].first; return m_part; }
	uint getEdgePart(const uint i) const { const uint j = i - getEdgeLocals(); if(j < m_exte.size()) return m_exte[j].first; return m_part; }
	uint getFacePart(const uint i) const { const uint j = i - getFaceLocals(); if(j < m_extf.size()) return m_extf[j].first; return m_part; }
	uint getBodyPart(const uint i) const { const uint j = i - getBodyLocals(); if(j < m_extb.size()) return m_extb[j].first; return m_part; }
	uint getQuadPart(const uint i) const { const uint j = i - getQuadLocals(); if(j < m_extq.size()) return m_extq[j].first; return m_part; }

	uint getNodeLink(const uint i) const { const uint j = i - getNodeLocals(); if(j < m_extn.size()) return m_extn[j].second; return i; }
	uint getEdgeLink(const uint i) const { const uint j = i - getEdgeLocals(); if(j < m_exte.size()) return m_exte[j].second; return i; }
	uint getFaceLink(const uint i) const { const uint j = i - getFaceLocals(); if(j < m_extf.size()) return m_extf[j].second; return i; }
	uint getBodyLink(const uint i) const { const uint j = i - getBodyLocals(); if(j < m_extb.size()) return m_extb[j].second; return i; }
	uint getQuadLink(const uint i) const { const uint j = i - getQuadLocals(); if(j < m_extq.size()) return m_extq[j].second; return i; }

	void setExternalNodes(const Buffer< pair<uint,uint> > &ext) { m_extn = ext; }
	void setExternalEdges(const Buffer< pair<uint,uint> > &ext) { m_exte = ext; }
	void setExternalFaces(const Buffer< pair<uint,uint> > &ext) { m_extf = ext; }
	void setExternalBodies(const Buffer< pair<uint,uint> > &ext){ m_extb = ext; }
	void setExternalQuads(const Buffer< pair<uint,uint> > &ext) { m_extq = ext; }

	const Buffer< pair<uint,uint> > &getExternalNodes() const { return m_extn; }
	const Buffer< pair<uint,uint> > &getExternalEdges() const { return m_exte; }
	const Buffer< pair<uint,uint> > &getExternalFaces() const { return m_extf; }
	const Buffer< pair<uint,uint> > &getExternalBodies() const { return m_extb; }
	const Buffer< pair<uint,uint> > &getExternalQuads() const { return m_extq; }

protected:

	uint m_part; // current part
	uint m_parts; // number of parts

	// external elements (which are the last elements of m_n, m_e, m_f, m_b and m_q)
	Buffer< pair<uint,uint> > m_extn;
	Buffer< pair<uint,uint> > m_exte;
	Buffer< pair<uint,uint> > m_extf;
	Buffer< pair<uint,uint> > m_extb;
	Buffer< pair<uint,uint> > m_extq;

	virtual const string getJRMeshType() const { return "JRMP"; }
	virtual bool loadJRMeshMore(std::ifstream &fs);
	virtual bool saveJRMeshMore(std::ofstream &fs) const;

};

}

#endif //_PARTMESH_HPP_INCLUDED_
