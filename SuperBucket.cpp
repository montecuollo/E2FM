/*
 * SuperBucket.cpp
 *
 *  Created on: 23/dic/2014
 *      Author: fernando
 */

#include "SuperBucket.h"
#include "EFMIndex.h"
#include "BitSetRLEEncoder.h"


namespace std {

SuperBucket::SuperBucket(EFMIndex *index,uint64_t superBucketNumber,uint32_t compactAlphabetSize) {
	this->index=index;
	this->superBucketNumber=superBucketNumber;
	charactersBitSet=new boost::dynamic_bitset<>();
	charactersBitSet->resize(compactAlphabetSize);
	previousSuperbucketsOccurrences=new uint64_t[compactAlphabetSize];
	for (uint32_t i=0;i<compactAlphabetSize;i++)
		previousSuperbucketsOccurrences[i]=0;
	alphaSize=0;
	maxOccurrence=0;
	decode_map=NULL;
	loaded=false;
}

SuperBucket::~SuperBucket() {
	charactersBitSet->resize(0);
	delete charactersBitSet;
	delete[] previousSuperbucketsOccurrences;
	if (decode_map!=NULL)
		delete[] decode_map;
}

void SuperBucket::load(BitReader *bitReader){
		uint64_t *superBucketStarts=index->superBucketsStarts;
		bitReader->gotoBookmark(new Bookmark(superBucketStarts[superBucketNumber],0));

		uint32_t compactAlphabetSize=index->compactAlphabetSize;
		bool bitmapStoredInCompressedFormat=bitReader->read(1)==1;
		if (bitmapStoredInCompressedFormat){
			uint64_t compressedBitSetSize=bitReader->getInt64();
			boost::dynamic_bitset<> *compressed=new boost::dynamic_bitset<>();
			compressed->resize(compressedBitSetSize);
			for (uint64_t k=0;k<compressedBitSetSize;k++)
				compressed->set(k,(bitReader->read(1)==1)?true:false);
			charactersBitSet=BitSetRLEEncoder::decompress(compressed);
		}
		else{
			uint64_t bitSetSize=bitReader->getInt64();
			charactersBitSet=new boost::dynamic_bitset<>();
			charactersBitSet->resize(bitSetSize);
			for (uint64_t k=0;k<bitSetSize;k++)
				charactersBitSet->set(k,(bitReader->read(1)==1)?true:false);
		}

		alphaSize=charactersBitSet->count();

		//compute decoding map for the super-bucket
		decode_map=new uint32_t[compactAlphabetSize];
		uint32_t code=0;
		for(uint32_t j=0; j<compactAlphabetSize; j++)
			if((*charactersBitSet)[j]){
				decode_map[code]=j;
			  code++;
		}


		if (superBucketNumber>0){
			int8_t neededBits=bitReader->getUByte();
			for (uint32_t k=0;k<compactAlphabetSize;k++)
				previousSuperbucketsOccurrences[k]=bitReader->read(neededBits);
		}

		loaded=true;

}


int64_t SuperBucket::getRemappedCode(uint32_t compactAlphabetCode) {
		if ((*charactersBitSet)[compactAlphabetCode]==0)
			return -1;
		uint32_t remappedCode=0;
		for (uint32_t i=0;i<compactAlphabetCode;i++)
			if ((*charactersBitSet)[i])
				remappedCode++;
		return remappedCode;
}



} /* namespace std */
