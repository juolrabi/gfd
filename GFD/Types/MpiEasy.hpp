#ifndef _MPIEASY_HPP_INCLUDED_
#define _MPIEASY_HPP_INCLUDED_

#include "Types.hpp"

namespace gfd
{

/*struct External
{
	uint part;
	uint link;
	External() {}
	External(const uint apart, const uint alink) {
		part = apart;
		link = alink;
	}
};
*/
void initMPI();
void finalizeMPI();
uint getMPIrank();
uint getMPIranks();
void sendMPI(void *data, const uint size, const uint rank, const uint tag);
void recvMPI(void *data, const uint size, const uint rank, const uint tag);
bool orMPI(bool condition);
bool andMPI(bool condition);
void sumMPI(uint *sum, const uint size);
void sumMPI(ullong *sum, const uint size);
void sumMPI(double *sum, const uint size);
void minMPI(uint *sum, const uint size);
void minMPI(ullong *sum, const uint size);
void minMPI(double *sum, const uint size);
void maxMPI(uint *sum, const uint size);
void maxMPI(ullong *sum, const uint size);
void maxMPI(double *sum, const uint size);
void barrierMPI();

};

#endif //_MPIEASY_HPP_INCLUDED_
