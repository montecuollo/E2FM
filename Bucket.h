/*
 * Bucket.h
 *
 *  Created on: 23/dic/2014
 *      Author: fernando
 */

#ifndef BUCKET_H_
#define BUCKET_H_
#include <boost/dynamic_bitset.hpp>
#include <sdsl/vectors.hpp>
#include "SuperBucket.h"
#include "Bookmark.h"

//#if BALANCED_TREE_MTF_DECODE == 0
#include "ArrayMTFList.h"
//#else
#include "BalancedTreeMTFList.h"
//#endif
//#include "ArrayMTFList.h"

#include "BitWriter.h"
#include "BitReader.h"
#include "RandomGenerator.h"


namespace std {

class EFMIndex;
class BucketEncoder;

typedef struct CryptoParams {

	CryptoParams(uint32_t *substitutionKey){
		this->substitutionKey=substitutionKey;
		this->offset=0;
	}

	CryptoParams(uint32_t *substitutionKey,uint32_t offset){
		this->substitutionKey=substitutionKey;
		this->offset=offset;
	}

   	uint32_t *substitutionKey;
	uint32_t offset;

} CryptoParams;


typedef struct EncodingResult {

	EncodingResult(uint32_t *sequence,uint64_t realLength){
		this->sequence=sequence;
		this->realLength=realLength;
	}

	uint32_t *sequence;
	uint64_t realLength;

} EncodingResult;


enum OffsetInformationsType {
	CountOnly,         //counts only the occurrences until a certain bucket offset
	RetrieveCharacter  //retrieve the character, reporting the occurrences count too
};

typedef struct OffsetInformations {
	uint64_t occurrencesCount;  //occurrences (whithin the bucket) of a character until a certain offset
	int64_t superBucketCode;   //super bucket code of the character located at a certain offset

	OffsetInformations(uint64_t occurrencesCount,int64_t superBucketCode){
		this->occurrencesCount=occurrencesCount;
		this->superBucketCode=superBucketCode;
	}

} OffsetInformations;




class Bucket {
public:
	Bucket(EFMIndex *index,SuperBucket *superBucket, uint64_t bucketNumber,uint32_t  *sbwt,uint64_t startPosition, uint64_t length,uint32_t superBucketAlphaSize);
	void buildCharactersMap();
	virtual ~Bucket();
	void write(BitWriter *bitWriter);
	void getCryptoParams();
	EncodingResult *encodeText();
	void encodeText(BucketEncoder *encoder);
	void dumpBucketBwtCodes();
	friend class EFMIndex;
	friend class BucketEncoder;
private:
	void writeText(BitWriter *bitWriter);
	void readHeader(BitReader *bitReader);
	void readText(BitReader *bitReader,RandomGenerator *bucketEncryptionRandomgenerator);
	uint64_t expandSingleCharacterUntilTo(uint64_t offset);
	uint32_t *computeBucketKeyStream(RandomGenerator *bucketEncryptionRandomGenerator,uint64_t bucketNumber,uint64_t startOffset,uint64_t endOffset,uint32_t bucketAlphaSize);
	uint64_t readTextUntilTo(BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator,uint64_t offset);
	uint32_t invertPoly(uint32_t compressedCharacter,uint32_t *substitutionKey,uint64_t substitutionKeyOffset);
	//OffsetInformations *getOffsetInformationsOld(BitReader *bitReader,int superBucketCode,int logicalOffset,OffsetInformationsType informationType);
	OffsetInformations *getOffsetInformations(BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator,uint32_t superBucketCode,uint64_t logicalOffset,OffsetInformationsType informationType);
	int64_t getRemappedCode(uint32_t superBucketCode);
	uint64_t occ(uint32_t bucketCode,uint64_t physicalOffset);
	int64_t getMarkedTextPosition(BitReader *bitReader,uint64_t bo);
	uint64_t getMarkedBwtPosition(BitReader *bitReader,uint64_t textPosition);
	void readMarkedRows(BitReader *bitReader);


	bool loaded; 					//true if the bucket header has been loaded
	bool loadingInProgress;	    //true if loading is in progress
	uint64_t bucketNumber; 		        //absolute number of this bucket
	bool isFirstInSuperBucket;    //true if this bucket is the first in its superbucket
	bool isOdd;
	uint32_t *sbwt;  					//reference to sbwt of entire super-text (used in compression phase)
	uint32_t *bucketText; 				//bucket portion of sbwt (used in decompression phase)

	uint64_t startPosition; 				//start position of bucket within sbwt
	uint64_t length;   					//length of the bucket (bucket ends at position startPosition+length-1)
	uint32_t superBucketAlphaSize;		//size of super-bucket alphabet
	uint32_t bucketAlphaSize; 			//size of local (bucket-relative) alphabet
	uint32_t bucketPolyAlphaSize;			//size of local (bucket-relative) alphabet after the poly-alphabetic substitution

	boost::dynamic_bitset<> *charactersBitSet;	   		//charactersBitSet[i]=1 if i_th super bucket character appears in this bucket


	sdsl::int_vector<0> charactersRemappings; //charactersRemappings[i] = code of the i_th super bucket character in this bucket

	uint32_t *decodeMap;
	uint64_t  *previousBucketsOccurrences;  //occurrences of all the superbucket's characters in previous buckets, used in decompression phase
	uint64_t previousBucketsOccurrencesMaximum;  //maximum occurrence of the superbucket's characters in previous buckets

	vector<uint64_t> *charactersPositions;  //array of vectors: charactersPositions[c][j] contains the position of the j^ + 1 occurrence of c inside the bucket
	SuperBucket *superBucket;			//super-bucket containing this bucket
	map<uint64_t,uint64_t> markedRows;  //marked rows whithin this bucket
	                          //keys are position of the suffix array rows from this bucket start
	                          //values are positions in text, divided by markingRate
	EFMIndex *index;               //index containing this bucket

	CryptoParams *cryptoParams;   //cryptographic parameters (substitution key and relative start offset)




	int64_t actualPhysicalOffset;			   //actual offset: the characters until
		                                       //the actualOffset are yet available (they have been loaded previously)
		                                       //this offset is the logical one (no matter of odd and even buckets)
	uint64_t actualCompressedOffset;  //stores the next position of the compressed
		                                       //from which is necessary to read characters
		                                       //in order to advance the actualOffset
	Bookmark *compressedTextStartBookmark;   //bookmark marking the start of the compressed text
	Bookmark *markedRowsStartBookmark;      //bookmark marking the start of the marked rows
	uint64_t compressedLength;		//compressed length of the bucket text
	int bitsPerTextCharacter;  //bits needed to encode a text character

	boost::dynamic_bitset<> *markedRowsBitmap;     //bitmap used to verify if a certain row is marked
	uint64_t *markedRowsPositions;  //marked rows positions: array of length equal to the number of the bucket's marked rows
	    					   //It contains in the i^th position, the text position of the i^th marked row of this bucket
	map<uint64_t,uint64_t> markedRowsCache;
	map<uint64_t,uint64_t> markedRowsReverseCache;

	#if BALANCED_TREE_MTF_DECODE == 0
	ArrayMTFList *mtfList;		   //Move To Front characters list
	#else
	BalancedTreeMTFList *mtfList;
	int timestamp; //current timestamp (used by BalancedTreeMTFList)
	#endif


	mutex textReadingMutex;
	mutex textReadingStatusMutex;
	int64_t textReadingUntilTo;//contains a value > -1 (the untilTo position) only if a text reading is in progress

	mutex markedRowsMutex;   //protects markedRowsCache,markedRowsReverseCache and allMarkedRowsLoaded
	bool allMarkedRowsLoaded; //true if all the marked rows information have been loaded yet


};

} /* namespace std */

#endif /* BUCKET_H_ */
