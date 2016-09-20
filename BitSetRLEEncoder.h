/*
 * BitSetRLEEncoder.h
 *
 *  Created on: 02/nov/2014
 *      Author: fernando
 */

#ifndef BITSETRLEENCODER_H_
#define BITSETRLEENCODER_H_
#include <boost/dynamic_bitset.hpp>

namespace std {

class BitSetRLEEncoder {
public:
	static boost::dynamic_bitset<> *compress(boost::dynamic_bitset<> *bitSet);
	static boost::dynamic_bitset<> *decompress(boost::dynamic_bitset<> *bitSet);
};

} /* namespace std */

#endif /* BITSETRLEENCODER_H_ */
