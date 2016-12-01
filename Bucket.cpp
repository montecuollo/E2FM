/*
 * Bucket.cpp
 *
 *  Created on: 23/dic/2014
 *      Author: fernando
 */

#include "Common.h"
#include "Bucket.h"
#include "Utils.h"
#include "EFMIndex.h"
#include <boost/dynamic_bitset.hpp>
#include "BitSetRLEEncoder.h"
#include "Bookmark.h"
#include "RandomGenerator.h"



namespace std {

Bucket::Bucket(EFMIndex *index,SuperBucket *superBucket, uint64_t bucketNumber,uint32_t  *sbwt,uint64_t startPosition, uint64_t length,uint32_t superBucketAlphaSize){
	this->index=index;
	this->superBucket=superBucket;
	this->bucketNumber=bucketNumber;
	this->sbwt=sbwt;
	this->startPosition=startPosition;
	this->length=length;
	this->superBucketAlphaSize=superBucketAlphaSize;
	this->loaded=false;
	this->loadingInProgress=false;

	this->actualPhysicalOffset=-1;
	this->actualCompressedOffset=0;
	this->bucketAlphaSize=0;

	this->charactersBitSet=NULL;
	this->decodeMap=NULL;
	this->previousBucketsOccurrences=NULL;
	this->charactersPositions=NULL;
	this->cryptoParams=NULL;
	this->mtfList=NULL;
	this->compressedTextStartBookmark=NULL;
	this->markedRowsStartBookmark=NULL;
	this->markedRowsPositions=NULL;
	this->markedRowsBitmap=NULL;
	this->isOdd=false;

	this->bucketText=NULL;
	this->textReadingUntilTo=-1;
	this->allMarkedRowsLoaded=false;
	this->compressedLength=0;
}


Bucket::~Bucket() {
	if (charactersBitSet!=NULL){
		charactersBitSet->resize(0);
		delete charactersBitSet;
	}
	if (decodeMap!=NULL)
		delete[] decodeMap;
	if (previousBucketsOccurrences!=NULL)
		delete[] previousBucketsOccurrences;
	if (charactersPositions!=NULL)
		delete[] charactersPositions;
	if (cryptoParams!=NULL)
		delete cryptoParams;
	if (compressedTextStartBookmark!=NULL)
		delete compressedTextStartBookmark;
	if (markedRowsStartBookmark!=NULL)
		delete markedRowsStartBookmark;
	if (markedRowsPositions!=NULL)
		delete[] markedRowsPositions;
	if (mtfList!=NULL)
		delete mtfList;
}


void Bucket::buildCharactersMap(){
	charactersRemappings.resize(superBucketAlphaSize);  //number of elements
	charactersRemappings.width(Utils::int_log2(bucketAlphaSize-1)); //bits for each bucket character code
	// compute characters remappings map and decoding map
	bucketAlphaSize=0;
	for(uint32_t j=0; j<superBucketAlphaSize; j++)
		if((*charactersBitSet)[j]){
		  uint32_t code=bucketAlphaSize;
		  charactersRemappings[j]= code;
		  if (decodeMap !=NULL)
			  decodeMap[code]= j;
		  bucketAlphaSize++;
		}
}


/**
 * Writes the bucket on resilient storage
 * @param bitWriter
 * @throws Exception
 */
void Bucket::write(BitWriter *bitWriter){

	//starting position of the super bucket within the BWT
	uint64_t sbStart=superBucket->superBucketNumber * index->superBucketSize;
	//starting position of the bucket within the BWT
	uint64_t start=startPosition;

	//start of the bucket in the compressed file
	index->bucketsStarts[bucketNumber]=bitWriter->getBookmark()->getPosition();
	bitWriter->writeInt64(0);  //leave space for the marked rows array start position
    bitWriter->writeInt64(0);  //leave space for the compressed text start position

    if(start != sbStart)   // if this isn't the first bucket of this super bucket
		if(!isOdd){ // and if it is not odd
			//writes a table containing, for each super bucket's character, the number
			//of its occurrences in the previous buckets of this superbucket
			uint8_t neededBits=Utils::int_log2(previousBucketsOccurrencesMaximum);
			bitWriter->writeUByte(neededBits);
			for(uint32_t k=0; k<superBucket->alphaSize; k++)
				bitWriter->integerEncode(previousBucketsOccurrences[k], neededBits);  //TODO: verificare che codifichi a 64 bit
	}


	// write characters bit set to the file, testing if the RLE compressed format
	// produces a gain in space over the uncompressed format
    boost::dynamic_bitset<> *bitsetToStore=charactersBitSet;
    bool storeInCompressedFormat=false;
    if (index->tryToStoreCompressedBitmaps){
    	boost::dynamic_bitset<> *compressed=BitSetRLEEncoder::compress(charactersBitSet);
    	storeInCompressedFormat=compressed->size() < charactersBitSet->size();
    	if (storeInCompressedFormat)
    		bitsetToStore=compressed;
    }
    bitWriter->write(1,(storeInCompressedFormat==true)?1:0);
	bitWriter->writeInt64(bitsetToStore->size());
	for (uint64_t i=0;i<bitsetToStore->size();i++)
			bitWriter->write(1, (*bitsetToStore)[i]);


	if (index->computeStatistics){
		index->statistics->bucketBitmapsCumulativeSize+=bitsetToStore->size();
	}


	//if necessary, writes a bit set identifying the marked row positions in this bucket
	if (index->markedRowsPercentage>0){
		  uint64_t bl=length;
		  boost::dynamic_bitset<> *bitSet=new boost::dynamic_bitset<>();
		  bitSet->resize(bl);
		  uint64_t row=0;
		  markedRowsPositions=new uint64_t[markedRows.size()];
		  for (auto it=markedRows.begin();it!=markedRows.end();++it){
				bitSet->set(it->first);
				markedRowsPositions[row]=it->second;
				row++;
		  }
		  boost::dynamic_bitset<> *compressedBitSet=BitSetRLEEncoder::compress(bitSet);
		  //Size of the marked rows bitmap
		  bitWriter->writeInt64(compressedBitSet->size());
		  //bitmap
		  for (uint64_t i=0;i<compressedBitSet->size();i++)
			bitWriter->write(1, (*compressedBitSet)[i]);
	}



	//write text only if this isn't a single-character bucket
	bitWriter->flush();
	if (bucketAlphaSize==1)
		compressedTextStartBookmark=bitWriter->getBookmark();
	else {
		writeText(bitWriter);
	}
	bitWriter->flush();

	//if necessary, writes an array containing, for each marked row,
	//the corresponding text position (divided by markingRate)
	uint64_t markedRowsArrayStartPosition=0;
	if (index->markedRowsPercentage>0){
	   markedRowsArrayStartPosition=bitWriter->getBookmark()->getPosition();
	   //array of the marked rows positions
	   uint32_t markedRowsValuesNeededBits=index->markedRowsValuesNeededBits;
	   if (index->computeStatistics){
	   		index->statistics->bucketAverageNumberOfMarkedRows+=markedRows.size();
	   		index->statistics->bucketMarkedRowsArrayAverageSize+=markedRows.size()*markedRowsValuesNeededBits;
	   	}

	   for (uint64_t i=0;i<markedRows.size();i++){
		  bitWriter->write(markedRowsValuesNeededBits,markedRowsPositions[i]);
	   }
	   bitWriter->flush();
	}

	//write at the bucket header's start the suspended information about the start of the compressed text
	Bookmark *eobBookmark=bitWriter->getBookmark();						  //bookmark pointing to the end of bucket
	Bookmark *sobBookmark=new Bookmark(index->bucketsStarts[bucketNumber],0);  //bookmark pointing to the start of bucket
	bitWriter->gotoBookmark(sobBookmark);
	bitWriter->writeInt64(markedRowsArrayStartPosition);
	bitWriter->writeInt64(compressedTextStartBookmark->getPosition());
	bitWriter->gotoBookmark(eobBookmark);  //restart to write from the end of the bucket
}


void Bucket::writeText(BitWriter *bitWriter){
	uint32_t *sequence;
	uint64_t sequenceLength;
	if (bucketText!=NULL){
		sequence=bucketText;
		sequenceLength=compressedLength;
	} else {
		//std::chrono::time_point<std::chrono::system_clock> startTime,endTime;
		//startTime=std::chrono::system_clock::now();
		if (cryptoParams==NULL)
			getCryptoParams();
		//endTime=std::chrono::system_clock::now();
		//double totalTime=(double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()/(double)1000;
		//cout << "CryptoParams: " <<totalTime << endl;

		//startTime=std::chrono::system_clock::now();
		EncodingResult *er=encodeText();
		//endTime=std::chrono::system_clock::now();
		//totalTime=(double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()/(double)1000;
		//cout << "EncodeText: " <<totalTime << endl;

		sequence=er->sequence;
		sequenceLength=er->realLength;
	}
	//startTime=std::chrono::system_clock::now();
	uint8_t neededBits = Utils::int_log2(bucketAlphaSize);  //+1 cause of RLE0
	bitWriter->writeInt64(sequenceLength);	 //it's part of the header at loading time (TODO: handle as suspended information)
	compressedTextStartBookmark=bitWriter->getBookmark(); //here really starts the bucket's text
	//if (bucketNumber==12799)
	//	cout << "absoluteBitPosition " << compressedTextStartBookmark->getAbsoluteBitPosition() << endl;
	for(uint64_t j=0;j<sequenceLength;j++){
		//if (bucketNumber==12799)
		//	cout << "sequence[" << j << "]=" << sequence[j] << endl;
		//if (bucketNumber==12799 && (j==511 || j==512 || j==513)){
		//	cout << "abs " <<bitWriter->getBookmark()->getAbsoluteBitPosition() <<endl;
		//	cout << "pos " <<bitWriter->getBookmark()->getPosition() <<endl;
		//	cout << "bb "<<bitWriter->getBookmark()->getBitsBuffer() <<endl;
		//	cout << "lb "<< (int64_t)bitWriter->getBookmark()->getLiveBits() <<endl;

		//}
		bitWriter->write(neededBits,sequence[j]);
	}
	//endTime=std::chrono::system_clock::now();
	//totalTime=(double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()/(double)1000;
	//cout << "EffectiveWrite: " <<totalTime << endl;
}

void Bucket::readText(BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator){
	if (bucketAlphaSize>1)
		//readTextUntilTo(bitReader,index.getBucketSize());
		actualPhysicalOffset=readTextUntilTo(bitReader,bucketEncryptionRandomGenerator,index->bucketSize-1);
	else
		//expandSingleCharacterUntilTo(index.getBucketSize());
		expandSingleCharacterUntilTo(index->bucketSize-1);


	//for (int64_t i=0;i<this->length;i++)
	//	cout <<   superBucket->decode_map[ decodeMap[bucketText[i]]] << endl;

}



uint32_t *Bucket::computeBucketKeyStream(RandomGenerator *bucketEncryptionRandomGenerator,uint64_t bucketNumber,uint64_t startOffset,uint64_t endOffset,uint32_t bucketAlphaSize){
		uint64_t keyFragmentLength=endOffset-startOffset+1;
		uint32_t* key;
		try{
			key=new uint32_t[keyFragmentLength];
		} catch(std::bad_alloc&){
			index->coutMutex.lock();
			cout << "keyFragmentLength =" << keyFragmentLength << endl;
			index->coutMutex.unlock();
			throw;
		}

		//Initialize a pseudo-casual number generator
		uint32_t maxValue=bucketAlphaSize;  //this upper bound is included in the range of generated numbers

		bucketEncryptionRandomGenerator->populateKeyStreamBuffer(bucketNumber,maxValue,startOffset,endOffset);
		for (uint64_t i=0;i<keyFragmentLength;i++)
			key[i]=bucketEncryptionRandomGenerator->nextInt();
		return key;
}


uint64_t Bucket::readTextUntilTo(BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator,uint64_t offset){
		//start from the character following the last available until now
		uint64_t skipBits=actualCompressedOffset*bitsPerTextCharacter;
		uint64_t start=compressedTextStartBookmark->getAbsoluteBitPosition() +skipBits;
		Bookmark *startBookmark=new Bookmark(start);
		Bookmark *endBookmark=new Bookmark(compressedTextStartBookmark->getAbsoluteBitPosition()+
						offset*bitsPerTextCharacter);
		if (!index->loadWholeFileInMemory)
			bitReader->setBlockSize((endBookmark->getPosition()-startBookmark->getPosition()+1));
		bitReader->gotoBookmark(startBookmark);

		uint32_t *substitutionKeyFragment;
		uint64_t fragmentStartOffset;
		uint64_t fragmentEndOffset;

		fragmentStartOffset=actualCompressedOffset;
		fragmentEndOffset=min(offset+10, compressedLength-1);
		//substitutionKeyFragment=index->encryptionManager->computeBucketKeyStream(bucketNumber,
		//		 fragmentStartOffset, fragmentEndOffset, bucketAlphaSize);
		substitutionKeyFragment=computeBucketKeyStream(bucketEncryptionRandomGenerator,bucketNumber,
				 fragmentStartOffset, fragmentEndOffset, bucketAlphaSize);


		//if not yet done, initialize the MTF List for this bucket
		if (mtfList==NULL)
		#if BALANCED_TREE_MTF_DECODE == 0
			mtfList=new ArrayMTFList(bucketAlphaSize);
		#else
			{mtfList=new BalancedTreeMTFList(bucketAlphaSize);
			timestamp=-1;
			}
		Pair *pair;
		BalancedTreeNode *node;
		#endif


		int64_t p=actualPhysicalOffset+1;
		int64_t i=actualCompressedOffset;
		int64_t compressedChar=-1;
		while (p<=offset  && i<compressedLength){
			//read the next character
			if (compressedChar==-1){
				Bookmark *b=bitReader->getBookmark();  //TODO: rimuovere (solo per debug)
				compressedChar=bitReader->read(bitsPerTextCharacter);
			}
			//invert the polyalphabetic substitution
			int64_t	h=i-fragmentStartOffset;   	 //the  element in position 0 whithin the key fragment is relative to the bucket's element of position fragmentStartOffset
			int64_t c=(compressedChar - substitutionKeyFragment[h])%(bucketAlphaSize+1); //+1 cause of RLE0
			while (c<0)
			   c=c+(bucketAlphaSize+1); //+1 cause of RLE0
			//invert RLE0 and MTF
			if (c>1){
				   int64_t position=c-1;  //inversion of RLE0
				   #if BALANCED_TREE_MTF_DECODE == 0
				   uint32_t bucketCode=mtfList->moveToFront(position);
				   #else
				   uint32_t bucketCode=mtfList->getCharacter(position,&node,&pair);
				   mtfList->moveToFront(bucketCode,timestamp,node,pair);
				   timestamp--;
				   #endif

				   if (isOdd)
						bucketText[length-p-1]=bucketCode;
				   else
						bucketText[p]=bucketCode;
				   charactersPositions[bucketCode].push_back(p);
				   p++;
				   i++;
				   compressedChar=-1;
				}
				else{
					string s="";
					int64_t x=0;    //number of consequent 0s and 1s
					int64_t inverted;

					while (i+x<compressedLength && (inverted=invertPoly(compressedChar,substitutionKeyFragment,i+x-fragmentStartOffset))<2){
						s+=to_string(inverted);
						x++;
						compressedChar=bitReader->read(bitsPerTextCharacter);
					};

					char *cs=(char *)s.c_str();
					char *endptr;
					uint64_t N = strtol(cs, &endptr, 2);
					if (endptr == cs) {
						throw runtime_error("Error in readTextUtilTo: No digits were found");
						exit(EXIT_FAILURE);
					}

					uint64_t m=pow(2.0, (double)x) + N-1;
					for (uint64_t j=0;j<m;j++){
						#if BALANCED_TREE_MTF_DECODE == 0
						uint32_t bucketCode=mtfList->get(0);	//bring the element in position 0
						#else
						uint32_t bucketCode=mtfList->getFirstCharacter();	//bring the element in position 0
						//if (bucketCodeNew!=bucketCode)
						//   cout << "Problem at coffset " << i << endl;
						#endif

						if (isOdd)
							bucketText[length-p-1]=bucketCode;
						else
							bucketText[p]=bucketCode;
						charactersPositions[bucketCode].push_back(p);
						p++;
					}
					i+=x;
				}
		}

		actualCompressedOffset=i;  //the next to be read

		//actualPhysicalOffset=p-1;  //the last read
		return p-1;

}


uint32_t Bucket::invertPoly(uint32_t compressedCharacter,uint32_t *substitutionKey,uint64_t substitutionKeyOffset){
	int64_t c=((int64_t)compressedCharacter - substitutionKey[substitutionKeyOffset])%(bucketAlphaSize+1); //+1 cause of RLE0
	while (c<0)
	   c=c+(bucketAlphaSize+1); //+1 cause of RLE0
	return c;
}


inline uint64_t Bucket::expandSingleCharacterUntilTo(uint64_t offset){
	//start from the character following the last available until now
	uint64_t p=actualPhysicalOffset+1;
	uint32_t bucketCode=0; //the char is unique in this bucket and so it's remapped on code 0
	while (p<=offset){
		if (isOdd)
			bucketText[length-p-1]=bucketCode;
		else
			bucketText[p]=bucketCode;
		charactersPositions[bucketCode].push_back(p);
		p++;
	}
	//actualPhysicalOffset=p-1;  //the last read
	return p-1;
}



void Bucket::getCryptoParams(){
	cryptoParams=index->getCryptoParams(bucketNumber,length,bucketAlphaSize);
}


//TODO: verificare bene l'equivalenza con quella scritta in Java
/**
 * Performs a single-pass Move to front, Run length encoding and Poly-alphabetic substitution
 * @return
 */
EncodingResult *Bucket::encodeText(){
	uint32_t *bucketText=new uint32_t[length];
	// ROVESCIA i bucket dispari, disponendo i caratteri dall'ultimo al primo
	if(isOdd)
		for(uint64_t j=0; j<length; j++)
			bucketText[length-j-1] = sbwt[startPosition+j];
	else
		memcpy(bucketText,sbwt+startPosition,length*sizeof(int));

	if (cryptoParams==NULL)
		getCryptoParams();
	uint32_t *polyAlphabeticKey=cryptoParams->substitutionKey;
	uint64_t keyLength=length;  //some final elements will not be used, cause of RLE0

	uint64_t keyStartOffset=0;

	#if BALANCED_TREE_MTF_ENCODE == 0
	ArrayMTFList *mtfList=new ArrayMTFList(bucketAlphaSize);
	#else
	BalancedTreeMTFList *mtfList=new BalancedTreeMTFList(bucketAlphaSize);
	Pair *pair;
	BalancedTreeNode *node;
	int64_t timestamp=-1;
	#endif

	uint32_t *output=new uint32_t[length];
	uint64_t outputSize=0;
	uint64_t i=0;
	uint64_t h=0;

	while (i<length){
			uint32_t c=bucketText[i];
			#if BALANCED_TREE_MTF_ENCODE == 0
			int64_t position=mtfList->indexOf(c);   //Move to front
			#else
			int64_t position=mtfList->indexOf(c,&node,&pair);   //Move to front
			//if (bucketNumber==12799 || bucketNumber==6399){
			//	index->coutMutex.lock();
			//	cout << "bucketText["<<i <<"]="<<c << "  ";
			//    cout << "position=" << position << "  ";
			//    index->coutMutex.unlock();
			//}
			#endif

			if (position>0){
				#if BALANCED_TREE_MTF_ENCODE == 0
				mtfList->moveToFront(position);
				#else
				mtfList->moveToFront(c,timestamp,node,pair);
				timestamp--;
				#endif

				h=(outputSize+keyStartOffset)%keyLength; //h is the index of the polyAlphabeticKey's element to use
				output[outputSize]=( position+1 + polyAlphabeticKey[h])%(bucketAlphaSize+1);  //position+1 is the RLE0 encoding of the MTF value position
				//if (bucketNumber==12799 || bucketNumber==6399){
				//	index->coutMutex.lock();
				//	cout << "timestamp=" << timestamp << "  ";
				//	cout << "h=" << h << "  ";
				//	cout << "keyStartOffset=" << keyStartOffset << "  ";
				//	cout << "keyLength=" << keyLength << "  ";
				//	cout << "poliAlphabeticKey["<<h <<"]="<<polyAlphabeticKey[h] << "  ";
				//	cout << "output["<<outputSize <<"]="<<output[outputSize] << endl;
				//	index->coutMutex.unlock();
				//}
				//index->coutMutex.unlock();
				outputSize++;
				i++;
			}else{
				/* RLE0 of a run of 0s */
				uint64_t m=1;
				#if BALANCED_TREE_MTF_ENCODE == 0
				while (i+m<length && mtfList->indexOf(bucketText[i+m])==0)
				#else
				while (i+m<length && mtfList->indexOf(bucketText[i+m],&node,&pair)==0)
				#endif
					m++;
				//At this point there is a run of 0s to encode with RLE0, whose length is m
				uint64_t nb= floor(log2(m+1));
				//compute the number N to encode in binary
				uint64_t N= m+1 - (uint64_t)pow(2, nb);
				boost::dynamic_bitset<>  binaryDigits(nb,N);
				for (int64_t j=binaryDigits.size()-1;j>=0;j--){
					h=(outputSize+keyStartOffset)%keyLength;
					if (binaryDigits[j]==0)	  //0a becomes 0
						output[outputSize]=(polyAlphabeticKey[h])%(bucketAlphaSize+1);  //0 is the result of RLE0 encoding
					else //0b becomes 10
						output[outputSize]=( 1 + polyAlphabeticKey[h])%(bucketAlphaSize+1);  //1 is the result of RLE0 encoding
					outputSize++;
				}
				i+=m;
			}
	}

	return new EncodingResult(output,outputSize);
}


void Bucket::dumpBucketBwtCodes(){
	uint64_t bucketBwtStartIndex=bucketNumber*index->bucketSize;
	for (int32_t i=0;i<=actualPhysicalOffset;i++ ){
		cout<< "bwt[" << bucketBwtStartIndex << "]=" << superBucket->decode_map[decodeMap[bucketText[i]]]<< endl;
				  //"\t(bucketText["<< i <<"])"<<endl  ;
		bucketBwtStartIndex++;
	}
}


void Bucket::readHeader(BitReader *bitReader){
	  uint64_t *bucketStarts=index->bucketsStarts;
	  if (!index->loadWholeFileInMemory)
		  bitReader->setBlockSize(16);
	  bitReader->gotoBookmark(new Bookmark(bucketStarts[bucketNumber],0));

	  //retrieve the marked rows array start position and create the corresponding bookmark
	  uint64_t markedRowsStartPosition=bitReader->getInt64();
	  markedRowsStartBookmark=new Bookmark(markedRowsStartPosition,0);
	  //retrieve the compressed text starting position and create the corresponding bookmark
	  uint64_t compressedTextStartPosition=bitReader->getInt64();
	  compressedTextStartBookmark=new Bookmark(compressedTextStartPosition,0);
	  if (!index->loadWholeFileInMemory)
		  bitReader->setBlockSize((uint32_t)(compressedTextStartBookmark->getPosition()-bucketStarts[bucketNumber]));
	  if( bucketNumber % (index->superBucketSize/index->bucketSize) !=0){      // if not the lastCharacterPosition bucket read occurrences from file
			Bookmark *thisBucketStart=NULL;
			if (isOdd){
				thisBucketStart=bitReader->getBookmark();
				bitReader->gotoBookmark(new Bookmark(bucketStarts[bucketNumber+1]+16,0));
					//previousBucketOccurrences starts at offset 16 of each bucket
			}
			uint8_t neededBits=bitReader->getUByte();
			previousBucketsOccurrences=new uint64_t[superBucketAlphaSize];
			for(uint32_t k=0; k<superBucketAlphaSize; k++)
				previousBucketsOccurrences[k]= bitReader->integerDecode(neededBits);

			if (isOdd)
				bitReader->gotoBookmark(thisBucketStart);

	  }


	  bool bitmapStoredInCompressedFormat=bitReader->read(1)==1;
	  if (bitmapStoredInCompressedFormat){
				uint64_t compressedBitSetSize=bitReader->getInt64();
				boost::dynamic_bitset<> *compressed=new boost::dynamic_bitset<>();
				compressed->resize(compressedBitSetSize);
				for (uint64_t k=0;k<compressedBitSetSize;k++)
					compressed->set(k,(bitReader->read(1)==1)?true:false);
				charactersBitSet=BitSetRLEEncoder::decompress(compressed);
	  } else{
				uint64_t bitSetSize=bitReader->getInt64();
				charactersBitSet=new boost::dynamic_bitset<>();
				charactersBitSet->resize(bitSetSize);
				for (uint64_t k=0;k<bitSetSize;k++)
					charactersBitSet->set(k,(bitReader->read(1)==1)?true:false);
	  }

	  //build decode map and compute bucketAlphaSize
	  decodeMap=new uint32_t[superBucketAlphaSize];
	  bucketAlphaSize=0;
	  for(uint32_t j=0; j<superBucketAlphaSize; j++)
		if((*charactersBitSet)[j]){
		  uint32_t code=bucketAlphaSize;
		  //decodeMap.put(code, j);
		  decodeMap[code]= j;
		  bucketAlphaSize++;
	  }


	  //if necessary, reads the bit set identifying the marked row positions in this bucket
	  if (index->markedRowsPercentage>0){
		  uint64_t compressedBitSetSize=bitReader->getInt64();
		  boost::dynamic_bitset<> *compressed=new boost::dynamic_bitset<>();
		  compressed->resize(compressedBitSetSize);
		  for (uint64_t k=0;k<compressedBitSetSize;k++)
		  		compressed->set(k,(bitReader->read(1)==1)?true:false);
		  markedRowsBitmap=BitSetRLEEncoder::decompress(compressed);
		  markedRowsBitmap->resize(length);
	  }


	  bitReader->gotoNextByteStart();
	  Bookmark *b=bitReader->getBookmark();
	  compressedLength=bitReader->getInt64();

	  //initialize the bucket text
	  bucketText=new uint32_t[length];
	  charactersPositions=new vector<uint64_t>[bucketAlphaSize];
	  bitsPerTextCharacter = Utils::int_log2(bucketAlphaSize);  //+1 cause of RLE0

	  //The bucket's header has been loaded
	  loaded=true;
	  //The bucket's text will be loaded after (when it's needed)
	  this->actualPhysicalOffset=-1;
	  this->actualCompressedOffset=0;

}



/**
 * Get the informations about a certain offset
 * @param bitReader
 * @param logicalOffset (no matter if the bucket is odd or not)
 * @param informationType
 * @return
 * @throws Exception
 */
/**
 * @param bitReader
 * @param superBucketCode
 * @param logicalOffset
 * @param informationType
 * @return
 * @throws Exception
 */
OffsetInformations *Bucket::getOffsetInformations(BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator,uint32_t superBucketCode,uint64_t logicalOffset,OffsetInformationsType informationType){
	if (informationType==OffsetInformationsType::CountOnly && (*charactersBitSet)[superBucketCode]==0)
		return new OffsetInformations(isFirstInSuperBucket?0:this->previousBucketsOccurrences[superBucketCode],-1);
	int64_t physicalOffset=isOdd?length-logicalOffset-1:logicalOffset;


	//Retrieves the actual physical offset (if a current text reading operation is in progress
	//afo is that before the starting of the current reading operation)
	int64_t afo;
	int64_t utfo;

	bool actAsReader=false;
	textReadingStatusMutex.lock();
	afo=actualPhysicalOffset;
	utfo=textReadingUntilTo;
	if (physicalOffset>afo && (utfo==-1 || physicalOffset>utfo))
		actAsReader=true;
	textReadingStatusMutex.unlock();


	if (actAsReader){   //ACT AS READER, BECAUSE no reading operation is in progress or the current reading operation doesn't arrive to the desired offset
			//acquire reading lock: wait until a previous reader ends
			textReadingMutex.lock();
			textReadingStatusMutex.lock();
			afo=actualPhysicalOffset;
			textReadingStatusMutex.unlock();

			if (physicalOffset > afo){
				//updates the utfo information, for the benefit of all the other threads
				textReadingStatusMutex.lock();
				textReadingUntilTo=physicalOffset;
				textReadingStatusMutex.unlock();

				if (bucketAlphaSize>1)
					afo=readTextUntilTo(bitReader,bucketEncryptionRandomGenerator, physicalOffset);
				else
					afo=expandSingleCharacterUntilTo(physicalOffset);

				//updates the utfo information, for the benefit of all the other threads, signaling that no reading is currently in progress
				textReadingStatusMutex.lock();
				textReadingUntilTo=-1;
				actualPhysicalOffset=afo;
				textReadingStatusMutex.unlock();
			}

			//release reading lock
			textReadingMutex.unlock();
	} else{
			while (utfo!=-1){  //waits the ending of the current reading operation (if there is a current reading operation)
				textReadingStatusMutex.lock();
				utfo=textReadingUntilTo;
				textReadingStatusMutex.unlock();
			}
	}


	uint32_t sbCode;
	int64_t bCode;
	if (informationType==OffsetInformationsType::CountOnly){
		sbCode=superBucketCode;
		bCode=getRemappedCode(sbCode);
	}
	else{
		bCode=bucketText[logicalOffset];
		sbCode=decodeMap[bCode];
	}

	//if (bucketNumber==12799 || bucketNumber==6399)
	//	dumpBucketBwtCodes();

	return new OffsetInformations(occ(bCode,physicalOffset),
			informationType==OffsetInformationsType::CountOnly?-1:sbCode);
}



/**
 * Return the number of occurrences of the character c till position q
 * This function is used by exactMatch and getRowPosition
 * @param c	the character
 * @param q	the position
 * @return
 */
uint64_t Bucket::occ(uint32_t bucketCode,uint64_t physicalOffset){
	uint64_t result=0;
	uint64_t p=0;
	//List<Integer> characterPositions=charactersPositions[bucketCode];
	uint64_t size=charactersPositions[bucketCode].size();

	uint32_t superBucketCode=decodeMap[bucketCode];
	if (isOdd){
		//determino tutte le occorrenze che hanno una posizione logica maggiore del logicalOffset (
		//(equivalentemente una posizione fisica minore del physicalOffset)
		//perchè le devo sottrarre al numero di occorrenze contenute per il carattere in questione nel super-blocco successivo
		bool foundGreaterOrEqual=false;
		while (p<size && !foundGreaterOrEqual){  //il ciclo si arresta appena trova una posizione fisica <= del physicalOffset
			if (charactersPositions[bucketCode][p]>=physicalOffset){
				foundGreaterOrEqual=true;
				result=p;
			}
			else p++;
		}
		if (!foundGreaterOrEqual)
			result=size;

		return this->previousBucketsOccurrences[superBucketCode] - result;  //it's not previous, it's the next Bucket occurrences
	}
	else{
		//determino tutte le occorrenze che hanno una posizione logica minore o uguale del logicalOffset (
		//(equivalentemente una posizione fisica minore o uguale del physicalOffset)
		//perchè le devo aggiungere al numero di occorrenze contenute per il carattere in questione nell'header del blocco
		bool foundGreater=false;
		while (p<size && !foundGreater){
			if (charactersPositions[bucketCode][p]>physicalOffset){
				foundGreater=true;
				result=p;
			}
			else p++;
		}
		if (!foundGreater)
			result=size;
		return (isFirstInSuperBucket?0:this->previousBucketsOccurrences[superBucketCode]) + result;
	}
}



int64_t Bucket::getRemappedCode(uint32_t superBucketCode) {
	if ((*charactersBitSet)[superBucketCode]==0)
		return -1;
	uint32_t remappedCode=0;
	for (uint32_t i=0;i<superBucketCode;i++)
		if ((*charactersBitSet)[i])
			remappedCode++;
	return remappedCode;
}


uint64_t Bucket::getMarkedBwtPosition(BitReader *bitReader,uint64_t textPosition){
	int64_t bo;

	markedRowsMutex.lock();
	auto it=markedRowsReverseCache.find(textPosition);
	if (it==markedRowsReverseCache.end())
		bo=-1;
	else
		bo=it->second;
	markedRowsMutex.unlock();

	if (bo!=-1)
		return bucketNumber*index->bucketSize+ bo;
	markedRowsMutex.lock();
	//load in memory all the marked rows informations
	if (!allMarkedRowsLoaded){
		readMarkedRows(bitReader);
		allMarkedRowsLoaded=true;
	}
	uint64_t retVal=bucketNumber*index->bucketSize+markedRowsReverseCache[textPosition];
	markedRowsMutex.unlock();
	return retVal;
}

/**
 * Return the marked text position
 * @param bo bucket offset
 * @return
 */
int64_t Bucket::getMarkedTextPosition(BitReader *bitReader,uint64_t bo){
	if (!(*markedRowsBitmap)[bo])
		return -1;
	int64_t markedRowPosition;

	markedRowsMutex.lock();
	auto it=markedRowsCache.find(bo);
	if (it==markedRowsCache.end())
		markedRowPosition=-1;
	else
		markedRowPosition=it->second;

	if (markedRowPosition==-1){
		//find the position of the marked row:
		//this is equal to the number of the 1 bits until the bo position
		int64_t pos=-1;
		for (uint64_t i=0;i<=bo;i++)
			if ((*markedRowsBitmap)[i])
			  pos++;

		int64_t skipBits=pos*index->markedRowsValuesNeededBits;
		int64_t markedRowAbsolutePosition=markedRowsStartBookmark->getAbsoluteBitPosition() +skipBits;
		Bookmark *markedRowBookmark=new Bookmark(markedRowAbsolutePosition);
		bitReader->gotoBookmark(markedRowBookmark);

		markedRowPosition=bitReader->read(index->markedRowsValuesNeededBits);
		markedRowsCache[bo]=markedRowPosition;
		markedRowsReverseCache[markedRowPosition]=bo;
	}
	markedRowsMutex.unlock();
	return markedRowPosition;
}


/**
 * Loads all the marked rows
 * @throws Exception
 */
void Bucket::readMarkedRows(BitReader *bitReader){
	bitReader->gotoBookmark(markedRowsStartBookmark);
	uint32_t actualBlockSize=bitReader->getBlockSize();

	//set the optimal block size to obtain with a single "logical" I/O operation
	//the entire set of the marked rows
	uint64_t allBits=markedRowsBitmap->count()*index->markedRowsValuesNeededBits;
	uint64_t endOfMarkedRowsAbsolutePosition=markedRowsStartBookmark->getAbsoluteBitPosition() +allBits;
	Bookmark *markedRowsEndBookmark=new Bookmark(endOfMarkedRowsAbsolutePosition);
	if (!index->loadWholeFileInMemory)
		bitReader->setBlockSize(markedRowsEndBookmark->getPosition()-markedRowsStartBookmark->getPosition()+1);

	for (uint64_t i=0;i<length;i++){ //iterate over the true bits
		if ((*markedRowsBitmap)[i]){
			uint64_t markedRowPosition=bitReader->read(index->markedRowsValuesNeededBits);
			markedRowsReverseCache[markedRowPosition]=i;
			markedRowsCache[i]=markedRowPosition;
		}
	}
	if (!index->loadWholeFileInMemory)
		bitReader->setBlockSize(actualBlockSize);
}






} /* namespace std */
