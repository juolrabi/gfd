#include "PartMesh.hpp"
#include "../Types/Text.hpp"
#include <fstream>
#include <iostream>

using namespace gfd;

PartMesh::PartMesh(const uint part, const uint parts, const uint dim)
: Mesh(dim)
{
	m_part = part;
	m_parts = parts;
}

void PartMesh::clear()
{
	Mesh::clear();

	m_extn.clear();
	m_exte.clear();
	m_extf.clear();
	m_extb.clear();
	m_extq.clear();
}

void PartMesh::createPartFromFlags(const Mesh &mesh)
{
	uint i, j;

	const uint nsize = mesh.getNodeSize();
	resizeNodeBuffer(nsize);
	Buffer<uint> n(nsize);
	for(i=0; i<nsize; i++)
	{
		const uint ipart = mesh.getNodeFlag(i) % m_parts;
		if(ipart != m_part) continue;
		n[i] = addNode(mesh.getNodePosition(i));
	}
	const uint nlocs = getNodeSize();
	m_extn.resize(nsize - nlocs);
	for(i=0; i<nsize; i++)
	{
		const uint ipart = mesh.getNodeFlag(i) % m_parts;
		if(ipart == m_part) continue;
		n[i] = addNode(mesh.getNodePosition(i));
		m_extn[n[i] - nlocs] = pair<uint,uint>(ipart, mesh.getNodeFlag(i) / m_parts);
	}
	const uint esize = mesh.getEdgeSize();
	resizeEdgeBuffer(esize);
	Buffer<uint> e(esize);
	for(i=0; i<esize; i++)
	{
		const uint ipart = mesh.getEdgeFlag(i) % m_parts;
		if(ipart != m_part) continue;
		const Buffer<uint> &par = mesh.getEdgeNodes(i);
		e[i] = addEdge(n[par[0]], n[par[1]]);
	}
	const uint elocs = getEdgeSize();
	m_exte.resize(esize - elocs);
	for(i=0; i<esize; i++)
	{
		const uint ipart = mesh.getEdgeFlag(i) % m_parts;
		if(ipart == m_part) continue;
		const Buffer<uint> &par = mesh.getEdgeNodes(i);
		e[i] = addEdge(n[par[0]], n[par[1]]);
		m_exte[e[i] - elocs] = pair<uint,uint>(ipart, mesh.getEdgeFlag(i) / m_parts);
	}
	const uint fsize = mesh.getFaceSize();
	resizeFaceBuffer(fsize);
	Buffer<uint> f(fsize);
	for(i=0; i<fsize; i++)
	{
		const uint ipart = mesh.getFaceFlag(i) % m_parts;
		if(ipart != m_part) continue;
		Buffer<uint> par = mesh.getFaceEdges(i);
		for(j=0; j<par.size(); j++) par[j] = e[par[j]];
		f[i] = addFace(par);
	}
	const uint flocs = getFaceSize();
	m_extf.resize(fsize - flocs);
	for(i=0; i<fsize; i++)
	{
		const uint ipart = mesh.getFaceFlag(i) % m_parts;
		if(ipart == m_part) continue;
		Buffer<uint> par = mesh.getFaceEdges(i);
		for(j=0; j<par.size(); j++) par[j] = e[par[j]];
		f[i] = addFace(par);
		m_extf[f[i] - flocs] = pair<uint,uint>(ipart, mesh.getFaceFlag(i) / m_parts);
	}
	const uint bsize = mesh.getBodySize();
	resizeBodyBuffer(bsize);
	Buffer<uint> b(bsize);
	for(i=0; i<bsize; i++)
	{
		const uint ipart = mesh.getBodyFlag(i) % m_parts;
		if(ipart != m_part) continue;
		Buffer<uint> par = mesh.getBodyFaces(i);
		for(j=0; j<par.size(); j++) par[j] = f[par[j]];
		b[i] = addBody(par);
	}
	const uint blocs = getBodySize();
	m_extb.resize(bsize - blocs);
	for(i=0; i<bsize; i++)
	{
		const uint ipart = mesh.getBodyFlag(i) % m_parts;
		if(ipart == m_part) continue;
		Buffer<uint> par = mesh.getBodyFaces(i);
		for(j=0; j<par.size(); j++) par[j] = f[par[j]];
		b[i] = addBody(par);
		m_extb[b[i] - blocs] = pair<uint,uint>(ipart, mesh.getBodyFlag(i) / m_parts);
	}
	const uint qsize = mesh.getQuadSize();
	resizeQuadBuffer(qsize);
	for(i=0; i<qsize; i++)
	{
		const uint ipart = mesh.getQuadFlag(i) % m_parts;
		if(ipart != m_part) continue;
		Buffer<uint> par = mesh.getQuadBodies(i);
		for(j=0; j<par.size(); j++) par[j] = b[par[j]];
		addQuad(par);
	}
	const uint qlocs = getQuadSize();
	m_extq.resize(qsize - qlocs);
	for(i=0; i<qsize; i++)
	{
		const uint ipart = mesh.getQuadFlag(i) % m_parts;
		if(ipart == m_part) continue;
		Buffer<uint> par = mesh.getQuadBodies(i);
		for(j=0; j<par.size(); j++) par[j] = b[par[j]];
		const uint q = addQuad(par);
		m_extq[q - qlocs] = pair<uint,uint>(ipart, mesh.getQuadFlag(i) / m_parts);
	}
}
/*
void PartMesh::createPart(const Mesh &mesh, const Buffer<uint> &npart, 
		const Buffer<uint> &epart, const Buffer<uint> &fpart, 
		const Buffer<uint> &bpart, const Buffer<uint> &qpart) {
	uint i, j;
	clear();
	
	const uint nsize = mesh.getNodeSize();
	const uint esize = mesh.getEdgeSize();
	const uint fsize = mesh.getFaceSize();
	const uint bsize = mesh.getBodySize();
	const uint qsize = mesh.getQuadSize();

	if(npart.size() < nsize) return;
	if(epart.size() < esize) return;
	if(fpart.size() < fsize) return;
	if(bpart.size() < bsize) return;
	if(qpart.size() < qsize) return;

	Buffer<uint> nlink(nsize);
	Buffer<uint> elink(esize);
	Buffer<uint> flink(fsize);
	Buffer<uint> blink(bsize);

	Buffer<uint> n;
	Buffer<uint> e;
	Buffer<uint> f;
	Buffer<uint> b;
	Buffer<uint> q;

	uint ns = 0;
	uint es = 0;
	uint fs = 0;
	uint bs = 0;
	uint qs = 0;

	Buffer<uint> ni(nsize, NONE);
	Buffer<uint> ei(esize, NONE);
	Buffer<uint> fi(fsize, NONE);
	Buffer<uint> bi(bsize, NONE);

	// find local elements
	Buffer<uint> locs(m_parts, 0);
	for(i=0; i<nsize; i++) {
		if(npart[i] == m_part) {
			ni[i] = ns;
			n.gather(i, ns);
		}
		else nlink[i] = locs[npart[i]]++;
	}
	locs.fill(0);
	for(i=0; i<esize; i++) {
		if(epart[i] == m_part) {
			ei[i] = es;
			e.gather(i, es);
		}
		else elink[i] = locs[epart[i]]++;
	}
	locs.fill(0);
	for(i=0; i<fsize; i++) {
		if(fpart[i] == m_part) {
			fi[i] = fs;
			f.gather(i, fs);
		}
		else flink[i] = locs[fpart[i]]++;
	}
	locs.fill(0);
	for(i=0; i<bsize; i++) {
		if(bpart[i] == m_part) {
			bi[i] = bs;
			b.gather(i, bs);
		}
		else blink[i] = locs[bpart[i]]++;
	}
	locs.clear();
	for(i=0; i<qsize; i++) {
		if(qpart[i] == m_part) q.gather(i, qs);
	}

	// numbers of local elements
	const uint n0 = ns;
	const uint e0 = es;
	const uint f0 = fs;
	const uint b0 = bs;

	// external elements
	for(i=0; i<qs; i++) {
		const Buffer<uint> &ele = mesh.getQuadBodies(q[i]);
		for(j=0; j<ele.size(); j++) {
			if(bi[ele[j]] != NONE) continue;
			bi[ele[j]] = bs;
			b.gather(ele[j], bs);
		}
	}
	for(i=0; i<bs; i++) {
		const Buffer<uint> &ele = mesh.getBodyFaces(b[i]);
		for(j=0; j<ele.size(); j++) {
			if(fi[ele[j]] != NONE) continue;
			fi[ele[j]] = fs;
			f.gather(ele[j], fs);
		}
	}
	for(i=0; i<fs; i++) {
		const Buffer<uint> &ele = mesh.getFaceEdges(f[i]);
		for(j=0; j<ele.size(); j++) {
			if(ei[ele[j]] != NONE) continue;
			ei[ele[j]] = es;
			e.gather(ele[j], es);
		}
	}
	for(i=0; i<es; i++) {
		const Buffer<uint> &ele = mesh.getEdgeNodes(e[i]);
		for(j=0; j<ele.size(); j++) {
			if(ni[ele[j]] != NONE) continue;
			ni[ele[j]] = ns;
			n.gather(ele[j], ns);
		}
	}

	// create mesh
	setMetric(mesh.getMetric());
	resizeNodeBuffer(ns);
	m_extn.resize(ns - n0);
	for(i=0; i<ns; i++) {
		addNode(mesh.getNodePosition(n[i]));
		setNodeWeight(i, mesh.getNodeWeight(n[i]));
		setNodeFlag(i, mesh.getNodeFlag(n[i]));
		if(i >= n0) m_extn[i - n0] = pair<uint,uint>(npart[n[i]], nlink[n[i]]);
	}
	nlink.clear();
	n.clear();
	resizeEdgeBuffer(es);
	m_exte.resize(es - e0);
	for(i=0; i<es; i++) {
		const Buffer<uint> &ele = mesh.getEdgeNodes(e[i]);
		addEdge(ni[ele[0]], ni[ele[1]]);
		setEdgeFlag(i, mesh.getEdgeFlag(e[i]));
		if(i >= e0) m_exte[i - e0] = pair<uint,uint>(epart[e[i]], elink[e[i]]);
	}
	ni.clear();
	elink.clear();
	e.clear();
	resizeFaceBuffer(fs);
	m_extf.resize(fs - f0);
	for(i=0; i<fs; i++) {
		Buffer<uint> ele = mesh.getFaceEdges(f[i]);
		for(j=0; j<ele.size(); j++) ele[j] = ei[ele[j]];
		addFace(ele);
		setFaceFlag(i, mesh.getFaceFlag(f[i]));
		if(i >= f0) m_extf[i - f0] = pair<uint,uint>(fpart[f[i]], flink[f[i]]);
	}
	ei.clear();
	flink.clear();
	f.clear();
	resizeBodyBuffer(bs);
	m_extb.resize(bs - b0);
	for(i=0; i<bs; i++) {
		Buffer<uint> ele = mesh.getBodyFaces(b[i]);
		for(j=0; j<ele.size(); j++) ele[j] = fi[ele[j]];
		addBody(ele);
		setBodyFlag(i, mesh.getBodyFlag(b[i]));
		if(i >= b0) m_extb[i - b0] = pair<uint,uint>(bpart[b[i]], blink[b[i]]);
	}
	fi.clear();
	blink.clear();
	b.clear();
	resizeQuadBuffer(qs);
	for(i=0; i<qs; i++) {
		Buffer<uint> ele = mesh.getQuadBodies(q[i]);
		for(j=0; j<ele.size(); j++) ele[j] = bi[ele[j]];
		addQuad(ele);
		setQuadFlag(i, mesh.getQuadFlag(q[i]));
	}
}
*/
void PartMesh::createPart(const Mesh &mesh, const Buffer<uint> &npart, 
		const Buffer<uint> &epart, const Buffer<uint> &fpart, 
		const Buffer<uint> &bpart, const Buffer<uint> &qpart) {
	uint i, j;
	clear();
	
	const uint nsize = mesh.getNodeSize();
	const uint esize = mesh.getEdgeSize();
	const uint fsize = mesh.getFaceSize();
	const uint bsize = mesh.getBodySize();
	const uint qsize = mesh.getQuadSize();

	if(npart.size() < nsize) return;
	if(epart.size() < esize) return;
	if(fpart.size() < fsize) return;
	if(bpart.size() < bsize) return;
	if(qpart.size() < qsize) return;

	Buffer<uint> nlink(nsize);
	Buffer<uint> elink(esize);
	Buffer<uint> flink(fsize);
	Buffer<uint> blink(bsize);

	Buffer<uint> n;
	Buffer<uint> e;
	Buffer<uint> f;
	Buffer<uint> b;
	Buffer<uint> q;

	uint ns = 0;
	uint es = 0;
	uint fs = 0;
	uint bs = 0;
	uint qs = 0;

	Buffer<uint> ni(nsize, NONE);
	Buffer<uint> ei(esize, NONE);
	Buffer<uint> fi(fsize, NONE);
	Buffer<uint> bi(bsize, NONE);

	// find local elements
	Buffer<uint> locs(m_parts, 0);
	for(i=0; i<nsize; i++) {
		if(npart[i] == m_part) {
			ni[i] = ns;
			n.gather(i, ns);
		}
		else if(npart[i] < m_parts) nlink[i] = locs[npart[i]]++;
	}
	locs.fill(0);
	for(i=0; i<esize; i++) {
		if(epart[i] == m_part) {
			ei[i] = es;
			e.gather(i, es);
		}
		else if(epart[i] < m_parts) elink[i] = locs[epart[i]]++;
	}
	locs.fill(0);
	for(i=0; i<fsize; i++) {
		if(fpart[i] == m_part) {
			fi[i] = fs;
			f.gather(i, fs);
		}
		else if(fpart[i] < m_parts) flink[i] = locs[fpart[i]]++;
	}
	locs.fill(0);
	for(i=0; i<bsize; i++) {
		if(bpart[i] == m_part) {
			bi[i] = bs;
			b.gather(i, bs);
		}
		else if(bpart[i] < m_parts) blink[i] = locs[bpart[i]]++;
	}
	locs.clear();
	for(i=0; i<qsize; i++) {
		if(qpart[i] == m_part) q.gather(i, qs);
	}

	// numbers of local elements
	const uint n0 = ns;
	const uint e0 = es;
	const uint f0 = fs;
	const uint b0 = bs;

	// external elements
	for(i=0; i<qs; i++) {
		const Buffer<uint> &ele = mesh.getQuadBodies(q[i]);
		for(j=0; j<ele.size(); j++) {
			if(bi[ele[j]] != NONE) continue;
			bi[ele[j]] = bs;
			b.gather(ele[j], bs);
		}
	}
	for(i=0; i<bs; i++) {
		const Buffer<uint> &ele = mesh.getBodyFaces(b[i]);
		for(j=0; j<ele.size(); j++) {
			if(fi[ele[j]] != NONE) continue;
			fi[ele[j]] = fs;
			f.gather(ele[j], fs);
		}
	}
	for(i=0; i<fs; i++) {
		const Buffer<uint> &ele = mesh.getFaceEdges(f[i]);
		for(j=0; j<ele.size(); j++) {
			if(ei[ele[j]] != NONE) continue;
			ei[ele[j]] = es;
			e.gather(ele[j], es);
		}
	}
	for(i=0; i<es; i++) {
		const Buffer<uint> &ele = mesh.getEdgeNodes(e[i]);
		for(j=0; j<ele.size(); j++) {
			if(ni[ele[j]] != NONE) continue;
			ni[ele[j]] = ns;
			n.gather(ele[j], ns);
		}
	}

	// create mesh
	setMetric(mesh.getMetric());
	resizeNodeBuffer(ns);
	m_extn.resize(ns - n0);
	for(i=0; i<ns; i++) {
		uint ii = n[i];
		addNode(mesh.getNodePosition(ii));
		setNodeWeight(i, mesh.getNodeWeight(ii));
		setNodeFlag(i, mesh.getNodeFlag(ii));
		if(i < n0) continue;
		while(npart[ii] >= m_parts) ii = npart[ii] - m_parts;
		m_extn[i - n0] = pair<uint,uint>(npart[ii], nlink[ii]);

	}
	nlink.clear();
	n.clear();
	resizeEdgeBuffer(es);
	m_exte.resize(es - e0);
	for(i=0; i<es; i++) {
		uint ii = e[i];
		const Buffer<uint> &ele = mesh.getEdgeNodes(ii);
		addEdge(ni[ele[0]], ni[ele[1]]);
		setEdgeFlag(i, mesh.getEdgeFlag(ii));
		if(i < e0) continue;
		while(epart[ii] >= m_parts) ii = epart[ii] - m_parts;
		m_exte[i - e0] = pair<uint,uint>(epart[ii], elink[ii]);
	}
	ni.clear();
	elink.clear();
	e.clear();
	resizeFaceBuffer(fs);
	m_extf.resize(fs - f0);
	for(i=0; i<fs; i++) {
		uint ii = f[i];
		Buffer<uint> ele = mesh.getFaceEdges(ii);
		for(j=0; j<ele.size(); j++) ele[j] = ei[ele[j]];
		addFace(ele);
		setFaceFlag(i, mesh.getFaceFlag(ii));
		if(i < f0) continue;
		while(fpart[ii] >= m_parts) ii = fpart[ii] - m_parts;
		m_extf[i - f0] = pair<uint,uint>(fpart[ii], flink[ii]);
	}
	ei.clear();
	flink.clear();
	f.clear();
	resizeBodyBuffer(bs);
	m_extb.resize(bs - b0);
	for(i=0; i<bs; i++) {
		uint ii = b[i];
		Buffer<uint> ele = mesh.getBodyFaces(ii);
		for(j=0; j<ele.size(); j++) ele[j] = fi[ele[j]];
		addBody(ele);
		setBodyFlag(i, mesh.getBodyFlag(ii));
		if(i < b0) continue;
		while(bpart[ii] >= m_parts) ii = bpart[ii] - m_parts;
		m_extb[i - b0] = pair<uint,uint>(bpart[ii], blink[ii]);
	}
	fi.clear();
	blink.clear();
	b.clear();
	resizeQuadBuffer(qs);
	for(i=0; i<qs; i++) {
		const uint ii = q[i];
		Buffer<uint> ele = mesh.getQuadBodies(ii);
		for(j=0; j<ele.size(); j++) ele[j] = bi[ele[j]];
		addQuad(ele);
		setQuadFlag(i, mesh.getQuadFlag(ii));
	}
}


void PartMesh::createCombined(Buffer<const PartMesh *> &mesh)
{
	uint i, j, k, exts;
	clear();

	// nodes
	Buffer< Buffer<uint> > n(m_parts);
	for(i=0; i<mesh.size(); i++)
	{
		Buffer<uint> &elem = n[mesh[i]->getPart()];
		elem.resize(mesh[i]->getNodeLocals());
		for(j=0; j<elem.size(); j++)
		{
			elem[j] = m_nsize;
			addNode(mesh[i]->getNodePosition(j));
			setNodeWeight(elem[j], mesh[i]->getNodeWeight(j));
		}
	}
	exts = 0;
	for(i=0; i<mesh.size(); i++)
	{
		for(j=mesh[i]->getNodeLocals(); j<mesh[i]->getNodeSize(); j++)
		{
			const uint part = mesh[i]->getNodePart(j);
			const uint link = mesh[i]->getNodeLink(j);
			Buffer<uint> &elem = n[part];
			if(link >= elem.size())
			{
				const uint osize = elem.size();
				elem.resize(link + 1);
				for(k=osize; k<=link; k++) elem[k] = NONE;
			}
			if(elem[link] != NONE) continue;
			elem[link] = m_nsize;
			addNode(mesh[i]->getNodePosition(j));
			setNodeWeight(elem[link], mesh[i]->getNodeWeight(j));
			m_extn.gather(pair<uint,uint>(part, link), exts);
		}
	}
	m_extn.resize(exts);

	// edges
	Buffer< Buffer<uint> > e(m_parts);
	for(i=0; i<mesh.size(); i++)
	{
		Buffer<uint> &elem = e[mesh[i]->getPart()];
		elem.resize(mesh[i]->getEdgeLocals());
		for(j=0; j<elem.size(); j++)
		{
			elem[j] = m_esize;
			Buffer<uint> par = mesh[i]->getEdgeNodes(j);
			for(k=0; k<par.size(); k++) par[k] = n[mesh[i]->getNodePart(par[k])][mesh[i]->getNodeLink(par[k])];
			addEdge(par[0], par[1]);
		}
	}
	exts = 0;
	for(i=0; i<mesh.size(); i++)
	{
		for(j=mesh[i]->getEdgeLocals(); j<mesh[i]->getEdgeSize(); j++)
		{
			const uint part = mesh[i]->getEdgePart(j);
			const uint link = mesh[i]->getEdgeLink(j);
			Buffer<uint> &elem = e[part];
			if(link >= elem.size())
			{
				const uint osize = elem.size();
				elem.resize(link + 1);
				for(k=osize; k<=link; k++) elem[k] = NONE;
			}
			if(elem[link] != NONE) continue;
			elem[link] = m_esize;
			Buffer<uint> par = mesh[i]->getEdgeNodes(j);
			for(k=0; k<par.size(); k++) par[k] = n[mesh[i]->getNodePart(par[k])][mesh[i]->getNodeLink(par[k])];
			addEdge(par[0], par[1]);
			m_exte.gather(pair<uint,uint>(part, link), exts);
		}
	}
	m_exte.resize(exts);
	n.clear();

	// faces
	Buffer< Buffer<uint> > f(m_parts);
	for(i=0; i<mesh.size(); i++)
	{
		Buffer<uint> &elem = f[mesh[i]->getPart()];
		elem.resize(mesh[i]->getFaceLocals());
		for(j=0; j<elem.size(); j++)
		{
			elem[j] = m_fsize;
			Buffer<uint> par = mesh[i]->getFaceEdges(j);
			for(k=0; k<par.size(); k++) par[k] = e[mesh[i]->getEdgePart(par[k])][mesh[i]->getEdgeLink(par[k])];
			addFace(par);
		}
	}
	exts = 0;
	for(i=0; i<mesh.size(); i++)
	{
		for(j=mesh[i]->getFaceLocals(); j<mesh[i]->getFaceSize(); j++)
		{
			const uint part = mesh[i]->getFacePart(j);
			const uint link = mesh[i]->getFaceLink(j);
			Buffer<uint> &elem = f[part];
			if(link >= elem.size())
			{
				const uint osize = elem.size();
				elem.resize(link + 1);
				for(k=osize; k<=link; k++) elem[k] = NONE;
			}
			if(elem[link] != NONE) continue;
			elem[link] = m_fsize;
			Buffer<uint> par = mesh[i]->getFaceEdges(j);
			for(k=0; k<par.size(); k++) par[k] = e[mesh[i]->getEdgePart(par[k])][mesh[i]->getEdgeLink(par[k])];
			addFace(par);
			m_extf.gather(pair<uint,uint>(part, link), exts);
		}
	}
	m_extf.resize(exts);
	e.clear();

	// bodies
	Buffer< Buffer<uint> > b(m_parts);
	for(i=0; i<mesh.size(); i++)
	{
		Buffer<uint> &elem = b[mesh[i]->getPart()];
		elem.resize(mesh[i]->getBodyLocals());
		for(j=0; j<elem.size(); j++)
		{
			elem[j] = m_bsize;
			Buffer<uint> par = mesh[i]->getBodyFaces(j);
			for(k=0; k<par.size(); k++) par[k] = f[mesh[i]->getFacePart(par[k])][mesh[i]->getFaceLink(par[k])];
			addBody(par);
		}
	}
	exts = 0;
	for(i=0; i<mesh.size(); i++)
	{
		for(j=mesh[i]->getBodyLocals(); j<mesh[i]->getBodySize(); j++)
		{
			const uint part = mesh[i]->getBodyPart(j);
			const uint link = mesh[i]->getBodyLink(j);
			Buffer<uint> &elem = b[part];
			if(link >= elem.size())
			{
				const uint osize = elem.size();
				elem.resize(link + 1);
				for(k=osize; k<=link; k++) elem[k] = NONE;
			}
			if(elem[link] != NONE) continue;
			elem[link] = m_bsize;
			Buffer<uint> par = mesh[i]->getBodyFaces(j);
			for(k=0; k<par.size(); k++) par[k] = f[mesh[i]->getFacePart(par[k])][mesh[i]->getFaceLink(par[k])];
			addBody(par);
			m_extb.gather(pair<uint,uint>(part, link), exts);
		}
	}
	m_extb.resize(exts);
	f.clear();

	// quads
	Buffer< Buffer<uint> > q(m_parts);
	for(i=0; i<mesh.size(); i++)
	{
		Buffer<uint> &elem = q[mesh[i]->getPart()];
		elem.resize(mesh[i]->getQuadLocals());
		for(j=0; j<elem.size(); j++)
		{
			elem[j] = m_qsize;
			Buffer<uint> par = mesh[i]->getQuadBodies(j);
			for(k=0; k<par.size(); k++) par[k] = b[mesh[i]->getBodyPart(par[k])][mesh[i]->getBodyLink(par[k])];
			addQuad(par);
		}
	}
	exts = 0;
	for(i=0; i<mesh.size(); i++)
	{
		for(j=mesh[i]->getQuadLocals(); j<mesh[i]->getQuadSize(); j++)
		{
			const uint part = mesh[i]->getQuadPart(j);
			const uint link = mesh[i]->getQuadLink(j);
			Buffer<uint> &elem = q[part];
			if(link >= elem.size())
			{
				const uint osize = elem.size();
				elem.resize(link + 1);
				for(k=osize; k<=link; k++) elem[k] = NONE;
			}
			if(elem[link] != NONE) continue;
			elem[link] = m_qsize;
			Buffer<uint> par = mesh[i]->getQuadBodies(j);
			for(k=0; k<par.size(); k++) par[k] = b[mesh[i]->getBodyPart(par[k])][mesh[i]->getBodyLink(par[k])];
			addQuad(par);
			m_extq.gather(pair<uint,uint>(part, link), exts);
		}
	}
	m_extq.resize(exts);
}

bool PartMesh::loadJRMesh(const std::string &path)
{
	if(m_parts == 1) return Mesh::loadJRMesh(path); // full load

	// load part
	std::ifstream fs(path.c_str(), std::ios::binary | std::ios::in);
	if(fs.fail()) return false;

	// header
	char id[4];
	uint parts;
	fs.read(id, 4);
	fs.read((char*)&parts, 4);
	if(std::string(id).compare(0, 4, "JRMP") != 0 || parts != m_parts)
	{
		fs.close();
		return false;
	}

	clear();
	uint i, size, index;

	// dimension
	fs.read((char*)&m_dim, sizeof(uint));

	// node positions
	fs.read((char*)&m_nsize, sizeof(uint));
	resizeNodeBuffer(m_nsize);
	fs.read((char*)&m_p[0], m_dim * m_nsize * sizeof(double));

	// edges
	fs.read((char*)&m_esize, sizeof(uint));
	resizeEdgeBuffer(m_esize);
	for(i=0; i<m_esize; i++)
	{
		fs.read((char*)&index, sizeof(uint));
		const uint n0 = index;
		fs.read((char*)&index, sizeof(uint));
		const uint n1 = index;
		m_e[i].n.resize(2);
		m_e[i].n[0] = n0;
		m_n[n0].e.push_back(i);
		m_e[i].n[1] = n1;
		m_n[n1].e.push_back(i);
	}

	// faces
	if(m_dim > 1)
	{
		fs.read((char*)&m_fsize, sizeof(uint));
		resizeFaceBuffer(m_fsize);
		for(i=0; i<m_fsize; i++)
		{
			fs.read((char*)&size, sizeof(uint));
			m_f[i].e.resize(uint(size));
			for(uint j=0; j<m_f[i].e.size(); j++)
			{
				fs.read((char*)&index, sizeof(uint));
				m_f[i].e[j] = index;
				m_e[index].f.push_back(i);
			}
			orderFaceEdges(i);
		}
	}

	// bodies
	if(m_dim > 2)
	{
		fs.read((char*)&m_bsize, sizeof(uint));
		resizeBodyBuffer(m_bsize);
		for(i=0; i<m_bsize; i++)
		{
			fs.read((char*)&size, sizeof(uint));
			m_b[i].f.resize(uint(size));
			for(uint j=0; j<m_b[i].f.size(); j++)
			{
				fs.read((char*)&index, sizeof(uint));
				m_b[i].f[j] = index;
				m_f[index].b.push_back(i);
			}
			orderBodyFaces(i);
		}
	}

	// quads
	if(m_dim > 3)
	{
		fs.read((char*)&m_qsize, sizeof(uint));
		resizeQuadBuffer(m_qsize);
		for(i=0; i<m_qsize; i++)
		{
			fs.read((char*)&size, sizeof(uint));
			m_q[i].b.resize(uint(size));
			for(uint j=0; j<m_q[i].b.size(); j++)
			{
				fs.read((char*)&index, sizeof(uint));
				m_q[i].b[j] = index;
				m_b[index].q.push_back(i);
			}
			orderQuadBodies(i);
		}
	}

	// flags
	fs.read((char*)&size, sizeof(uint));
	m_nflag.resize(size);
	if(size > 0) fs.read((char*)&m_nflag[0], size * sizeof(uint));
	fs.read((char*)&size, sizeof(uint));
	m_eflag.resize(size);
	if(size > 0) fs.read((char*)&m_eflag[0], size * sizeof(uint));
	fs.read((char*)&size, sizeof(uint));
	m_fflag.resize(size);
	if(size > 0) fs.read((char*)&m_fflag[0], size * sizeof(uint));
	fs.read((char*)&size, sizeof(uint));
	m_bflag.resize(size);
	if(size > 0) fs.read((char*)&m_bflag[0], size * sizeof(uint));
	fs.read((char*)&size, sizeof(uint));
	m_qflag.resize(size);
	if(size > 0) fs.read((char*)&m_qflag[0], size * sizeof(uint));

	// circumcenter computation
	fs.read((char*)&size, sizeof(uint));
	m_m.resize(size);
	if(size > 0) fs.read((char*)&m_m[0], size * sizeof(double));
	fs.read((char*)&size, sizeof(uint));
	m_w.resize(size);
	if(size > 0) fs.read((char*)&m_w[0], size * sizeof(double));

	// external elements
	fs.read((char*)&size, sizeof(uint));
	m_extn.resize(size);
	if(size > 0) fs.read((char*)&m_extn[0], size * sizeof(pair<uint,uint>));
	fs.read((char*)&size, sizeof(uint));
	m_exte.resize(size);
	if(size > 0) fs.read((char*)&m_exte[0], size * sizeof(pair<uint,uint>));
	fs.read((char*)&size, sizeof(uint));
	m_extf.resize(size);
	if(size > 0) fs.read((char*)&m_extf[0], size * sizeof(pair<uint,uint>));
	fs.read((char*)&size, sizeof(uint));
	m_extb.resize(size);
	if(size > 0) fs.read((char*)&m_extb[0], size * sizeof(pair<uint,uint>));
	fs.read((char*)&size, sizeof(uint));
	m_extq.resize(size);
	if(size > 0) fs.read((char*)&m_extq[0], size * sizeof(pair<uint,uint>));

	fs.close();
	return true;
}

bool PartMesh::saveJRMesh(const std::string &path) const
{
	if(m_parts == 1) return Mesh::saveJRMesh(path); // save full

	// save part
	std::ofstream fs(path.c_str(), std::ios_base::binary | std::ios::trunc);
	if(fs.fail()) return false;

	// header
	const char *id = "JRMP";
	fs.write(id, 4);
	fs.write((char*)&m_parts, sizeof(uint));

	uint i, j, size;

	// dimension
	fs.write((char*)&m_dim, sizeof(uint));

	// node positions
	fs.write((char*)&m_nsize, sizeof(uint));
	fs.write((char*)&m_p[0], m_dim * m_nsize * sizeof(double));

	// edges
	fs.write((char*)&m_esize, sizeof(uint));
	for(i=0; i<m_esize; i++)
	{
		const Buffer<uint> &n = getEdgeNodes(i);
		fs.write((char*)&n[0], sizeof(uint));
		fs.write((char*)&n[1], sizeof(uint));
	}

	// faces
	if(m_dim > 1)
	{
		fs.write((char*)&m_fsize, sizeof(uint));
		for(i=0; i<m_fsize; i++)
		{
			const Buffer<uint> &e = getFaceEdges(i);
			size = e.size();
			fs.write((char*)&size, sizeof(uint));
			for(j=0; j<size; j++) fs.write((char*)&e[j], sizeof(uint));
		}
	}

	// bodies
	if(m_dim > 2)
	{
		fs.write((char*)&m_bsize, sizeof(uint));
		for(i=0; i<m_bsize; i++)
		{
			const Buffer<uint> &f = getBodyFaces(i);
			size = f.size();
			fs.write((char*)&size, sizeof(uint));
			for(j=0; j<size; j++) fs.write((char*)&f[j], sizeof(uint));
		}
	}

	// quads
	if(m_dim > 3)
	{
		fs.write((char*)&m_qsize, sizeof(uint));
		for(i=0; i<m_qsize; i++)
		{
			const Buffer<uint> &b = getQuadBodies(i);
			size = b.size();
			fs.write((char*)&size, sizeof(uint));
			for(j=0; j<size; j++) fs.write((char*)&b[j], sizeof(uint));
		}
	}

	// flags
	size = m_nflag.size();
	while(size > 0 && m_nflag[size - 1] == 0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_nflag[0], size * sizeof(uint));
	size = m_eflag.size();
	while(size > 0 && m_eflag[size - 1] == 0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_eflag[0], size * sizeof(uint));
	size = m_fflag.size();
	while(size > 0 && m_fflag[size - 1] == 0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_fflag[0], size * sizeof(uint));
	size = m_bflag.size();
	while(size > 0 && m_bflag[size - 1] == 0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_bflag[0], size * sizeof(uint));
	size = m_qflag.size();
	while(size > 0 && m_qflag[size - 1] == 0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_qflag[0], size * sizeof(uint));

	// circumcenter computation
	size = m_m.size();
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_m[0], size * sizeof(double));
	size = m_w.size();
	while(size > 0 && m_w[size - 1] == 0.0) size--;
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_w[0], size * sizeof(double));

	// parts and links
	size = m_extn.size();
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_extn[0], size * sizeof(double));
	size = m_exte.size();
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_exte[0], size * sizeof(double));
	size = m_extf.size();
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_extf[0], size * sizeof(double));
	size = m_extb.size();
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_extb[0], size * sizeof(double));
	size = m_extq.size();
	fs.write((char*)&size, sizeof(uint));
	if(size > 0) fs.write((char*)&m_extq[0], size * sizeof(double));

	fs.close();
	return true;
}
