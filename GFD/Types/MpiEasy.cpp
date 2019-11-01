#include "MpiEasy.hpp"
#include <mpi.h>

namespace gfd
{

const uint MAXMPISIZE = 4000;

void initMPI()
{
	MPI_Init(NULL, NULL);
}

void finalizeMPI()
{
	MPI_Finalize();
}

uint getMPIrank()
{
	uint rank;
	MPI_Comm_rank(MPI_COMM_WORLD, (int*)&rank);
	return rank;
}

uint getMPIranks()
{
	uint ranks;
	MPI_Comm_size(MPI_COMM_WORLD, (int*)&ranks);
	return ranks;
}

void sendMPI(void *data, const uint size, const uint rank, const uint tag)
{
	uint j;
	for(j=0; j+MAXMPISIZE<size; j+=MAXMPISIZE) MPI_Send(&((char *)data)[j], MAXMPISIZE, MPI_CHAR, rank, tag, MPI_COMM_WORLD);
	MPI_Send(&((char *)data)[j], size - j, MPI_CHAR, rank, tag, MPI_COMM_WORLD);
}

void recvMPI(void *data, const uint size, const uint rank, const uint tag)
{
	uint j;
	for(j=0; j+MAXMPISIZE<size; j+=MAXMPISIZE) MPI_Recv(&((char *)data)[j], MAXMPISIZE, MPI_CHAR, rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&((char *)data)[j], size - j, MPI_CHAR, rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

bool orMPI(const bool condition) {
	uint val = (condition ? 1 : 0);
	maxMPI(&val, 1);
	return val == 1;
}
bool andMPI(const bool condition) {
	uint val = (condition ? 1 : 0);
	minMPI(&val, 1);
	return val == 1;
}

void sumMPI(uint *sum, const uint size)
{
	uint *send = new uint[size];
	memcpy(send, sum, size * sizeof(uint));
	MPI_Reduce(send, sum, size, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(sum, size, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	delete[] send;
}
void sumMPI(ullong *sum, const uint size)
{
	ullong *send = new ullong[size];
	memcpy(send, sum, size * sizeof(ullong));
	MPI_Reduce(send, sum, size, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(sum, size, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	delete[] send;
}
void sumMPI(double *sum, const uint size)
{
	double *send = new double[size];
	memcpy(send, sum, size * sizeof(double));
	MPI_Reduce(send, sum, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(sum, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	delete[] send;
}

void minMPI(uint *sum, const uint size)
{
	uint *send = new uint[size];
	memcpy(send, sum, size * sizeof(uint));
	MPI_Reduce(send, sum, size, MPI_UNSIGNED, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Bcast(sum, size, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	delete[] send;
}
void minMPI(ullong *sum, const uint size)
{
	ullong *send = new ullong[size];
	memcpy(send, sum, size * sizeof(ullong));
	MPI_Reduce(send, sum, size, MPI_UNSIGNED_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Bcast(sum, size, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	delete[] send;
}
void minMPI(double *sum, const uint size)
{
	double *send = new double[size];
	memcpy(send, sum, size * sizeof(double));
	MPI_Reduce(send, sum, size, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Bcast(sum, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	delete[] send;
}

void maxMPI(uint *sum, const uint size)
{
	uint *send = new uint[size];
	memcpy(send, sum, size * sizeof(uint));
	MPI_Reduce(send, sum, size, MPI_UNSIGNED, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(sum, size, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	delete[] send;
}
void maxMPI(ullong *sum, const uint size)
{
	ullong *send = new ullong[size];
	memcpy(send, sum, size * sizeof(ullong));
	MPI_Reduce(send, sum, size, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(sum, size, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
	delete[] send;
}
void maxMPI(double *sum, const uint size)
{
	double *send = new double[size];
	memcpy(send, sum, size * sizeof(double));
	MPI_Reduce(send, sum, size, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Bcast(sum, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	delete[] send;
}

void barrierMPI()
{
	MPI_Barrier(MPI_COMM_WORLD);
}

}
