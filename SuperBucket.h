/*
 * SuperBucket.h
 *
 *  Created on: 23/dic/2014
 *      Author: fernando
 */

#ifndef SUPERBUCKET_H_
#define SUPERBUCKET_H_

#include <boost/dynamic_bitset.hpp>
#include <sdsl/int_vector.hpp>
#include "BitReader.h"

namespace std {

class EFMIndex;


class SuperBucket {
public:
	SuperBucket(EFMIndex *index,uint64_t superBucketNumber,uint32_t compactAlphabetSize);
	virtual ~SuperBucket();
	void load(BitReader *bitReader);
	int64_t getRemappedCode(uint32_t compactAlphabetCode);
	friend class EFMIndex;
	friend class Bucket;
private:
	EFMIndex *index;
	uint64_t superBucketNumber;
	boost::dynamic_bitset<> *charactersBitSet; //identifies the compact alphabet's symbols occurring in this super-bucket
	sdsl::int_vector<0> charactersRemappings;
	uint32_t *decode_map;
	uint64_t *previousSuperbucketsOccurrences;  //number of occurrences of characters in previous super buckets
	uint32_t alphaSize;
	uint64_t maxOccurrence; //maximum value stored in previousSuperbucketsOccurrences array
	bool loaded;
};

} /* namespace std */

#endif /* SUPERBUCKET_H_ */
