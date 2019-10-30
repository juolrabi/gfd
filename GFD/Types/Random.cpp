#include "Random.hpp"

using namespace gfd;

const uint MT_M = 397;
const uint LAST_31BITS = 2147483647;
const uint THE_32ND_BIT = 2147483648;
const uint EVEN_ODD[2] = {0, 2567483615};
const double UINT2DOUBLE = 1.0 / 4294967296.0;
const double INT2DOUBLE = 1.0 / 2147483648.0;

uint Random::ms_seed = 0;

void Random::init(const uint seed) {
    m_mt[0] = seed;
    for(m_i=1; m_i<MT_N; m_i++) {
    	//m_mt[m_i] = (69069 * m_mt[m_i-1]) & 0xffffffff;
		m_mt[m_i] = (1812433253 * (m_mt[m_i-1] ^ (m_mt[m_i-1] >> 30)) + m_i);
    }
}

uint Random::getUint() {
	if(m_i >= MT_N) {
		uint i = 0;
		while(i < MT_N - MT_M) {
			const uint y = (THE_32ND_BIT & m_mt[i]) | (LAST_31BITS & m_mt[i+1]);
			m_mt[i] = (m_mt[i + MT_M] ^ (y >> 1) ^ EVEN_ODD[y & 1]);
			i++;
		}
		while(i < MT_N - 1) {
			const uint y = (THE_32ND_BIT & m_mt[i]) | (LAST_31BITS & m_mt[i+1]);
			m_mt[i] = (m_mt[i + MT_M - MT_N] ^ (y >> 1) ^ EVEN_ODD[y & 1]);
			i++;
		}
		const uint y = (THE_32ND_BIT & m_mt[i]) | (LAST_31BITS & m_mt[0]);
		m_mt[i] = m_mt[MT_M-1] ^ (y >> 1) ^ EVEN_ODD[y & 1];
		m_i = 0;
	}

	uint y = m_mt[m_i++];
	y ^=  y >> 11;
	y ^= (y << 7) & 2636928640;
	y ^= (y << 15) & 4022730752;
	y ^=  y >> 18;
	return y;
}

double Random::getUniform() {
  return UINT2DOUBLE * getUint();
}

double Random::getExponential() {
  return -log(1.0 - UINT2DOUBLE * getUint());
}

double Random::getGaussian() {
  return sqrt(2.0 * getExponential()) * cos(PIx2 * getUniform());
}

Vector2 Random::getUniformSphere2() {
  Vector2 v;
  do {
    v.x = INT2DOUBLE * int(getUint());
    v.y = INT2DOUBLE * int(getUint());
  } while(v.lensq() > 1.0);
  return v;
}

Vector3 Random::getUniformSphere3() {
  Vector3 v;
  do {
    v.x = INT2DOUBLE * int(getUint());
    v.y = INT2DOUBLE * int(getUint());
    v.z = INT2DOUBLE * int(getUint());
  } while(v.lensq() > 1.0);
  return v;
}

Vector4 Random::getUniformSphere4() {
  Vector4 v;
  do {
    v.x = INT2DOUBLE * int(getUint());
    v.y = INT2DOUBLE * int(getUint());
    v.z = INT2DOUBLE * int(getUint());
    v.t = INT2DOUBLE * int(getUint());
  } while(v.lensq() > 1.0);
  return v;
}

