/*
 * EFMIndex.h
 *
 *  Created on: 19/dic/2014
 *      Author: fernando
 */

#ifndef EFMINDEX_H_
#define EFMINDEX_H_
#include <vector>
#include "Statistics.h"
#include "ScrambledSuperAlphabet.h"
#include "ecrypt-portable.h"
#include "EncryptionManager.h"
#include "FastBWTransformer.h"
#include <boost/dynamic_bitset.hpp>
#include <sdsl/int_vector.hpp>
#include "SuperBucket.h"
#include "Bucket.h"
#include "BitWriter.h"
#include "BitReader.h"
#include "FileBitReader.h"
#include "RandomGenerator.h"
#include <condition_variable>
#include "Common.h"

namespace std {

typedef RandomGenerator* PRG;
typedef BitReader* PBR;


class LocateNotSupportedException: public exception
{
  virtual const char* what() const throw()
  {
    return "Locate not supported by this index";
  }
};





/**
 * Interval of integer numbers.
 * This class is used
 * 1) to contain the result of an exact Match or one of an inexact match results;
 * 2) to express the scrambled alphabet characters intervals corresponding to a variable
 *    super pattern character.
 *
 * @author fernando
 *
 */
typedef struct Interval {
	Interval(uint64_t sp,uint64_t ep){
		this->sp=sp;
		this->ep=ep;
	}

	Interval(){}

	uint64_t sp;  //the interval starts at this position
	uint64_t ep;  //the interval ends at this position

	/**
	 * @return Length of the interval
	 */
	uint64_t getLength(){
		return ep-sp+1;
	}

	friend ostream& operator<<(ostream& os, const Interval& i);

} Interval;


struct SuperPattern;


typedef struct CompactAlphabetInterval:Interval {
	CompactAlphabetInterval(SuperPattern *superPattern,uint64_t sp,uint64_t ep){
		this->superPattern=superPattern;
		this->sp=sp;
		this->ep=ep;
	}

	CompactAlphabetInterval(){}

	SuperPattern *superPattern;  //super-pattern from which this interval stems



	friend ostream& operator<<(ostream& os, const Interval& i);

} CompactAlphabetInterval;






enum SuperPatternCharacterType {
	/**
	 * single character
	 */
	fixedCharacter,

	/**
	 * character with a prefix and/or a suffix of wildcard symbols ?
	 * (a ? stands for "any original alphabet character")
	 */
	variableCharacter
};


typedef struct SuperPatternCharacter {
	/**
	 * text of this supercharacter respect to the original alphabet
	 */
	string text;
	/**
	 * fixed character (without prefix and suffix) or variable character (with a prefix or a suffix or both)
	 */
	SuperPatternCharacterType  type;
	/**
	 * length of the wildcard ? (any original alphabet character) prefix
	 */
	uint32_t prefixLength;
	/**
	 * length of the wildcard ? (any original alphabet character) suffix
	 */
	uint32_t suffixLength;
	/**
	 * all possibile values of the prefix
	 */
	uint32_t *prefixValues;
	uint32_t prefixValuesSize; //real size of the prefixValues array

	/**
	 * Last value of the range of values corresponding to the prefix
	 * (suffixRangeStart is not stored because it's always 0)
	 */
	uint32_t suffixRangeEnd;

	/**
	 * start of the fixed part within the pattern. This can be <0 if this
	 * super-character has a prefix and in this case the prefix will contain
	 * exactly (-fixedPartStartInPattern) wildcard ? symbols
	 */
	uint32_t fixedPartStartInPattern;

	/**
	 * end of the fixed part within the pattern.
	 */
	uint32_t fixedPartEndInPattern;

	/**
	 * Fixed value of super-character code.
	 */
	uint32_t fixedValue;

	/**
	 * Intervals of super alphabet containing the real super character codes
	 */
	vector<CompactAlphabetInterval*> compactAlphabetIntervals;


	friend ostream& operator<<(ostream& os, const SuperPatternCharacter& sc);

} SuperPatternCharacter;


typedef struct SuperPattern{

	SuperPattern(){
		length=-1;
		shiftAmount=-1;
		isHat=false;
	}

	/** Array of symbols of two possible kinds:
	 *   1) super-characters;
	 *   2) wildcard ? characters, encoded as superAlphabetSize
	 */
	SuperPatternCharacter **text;
	//number of super-characters in this super-pattern
	uint32_t length;

	/**
	 * Number of right shift positions from which this super-pattern
	 * has been obtained
	 */
	uint32_t shiftAmount;

	/**
	 * True if it's a "hat super-pattern" (not including the last variable super-character of the whole super-pattern)
	 */
	bool isHat;

	friend ostream& operator<<(ostream& os, const SuperPattern& sp);


} SuperPattern;


typedef struct Match {
	Match();
	virtual ~Match();

	SuperPattern *superPattern;
	Interval *suffixArrayInterval;
	vector<uint64_t> *occurrences;

	friend ostream& operator<<(ostream& os, const Match& m);

} Match;



class OriginalSymbolsFinder {
public:
	uint64_t firstTextPosition;
	uint64_t lastTextPosition;
	string &text;
	set<char> chars;

	OriginalSymbolsFinder(string &text,uint64_t firstTextPosition,uint64_t lastTextPosition):text(text){
		this->firstTextPosition=firstTextPosition;
		this->lastTextPosition=lastTextPosition;
	}

	void run(){
		for (uint64_t i=firstTextPosition;i<=lastTextPosition;i++)
				chars.insert(text[i]);
	}

};


class EFMIndex{
public:
	EFMIndex();
	virtual ~EFMIndex();
	uint64_t getBucketSize() const;
	void setBucketSize(uint64_t bucketSize);
	uint64_t getBucketsNumber();
	uint64_t getSuperBucketsNumber();
	uint64_t getNumberOfLoadedBuckets();
	uint64_t getNumberOfLoadedSuperBuckets();
	bool isComputeStatistics() const;
	void setComputeStatistics(bool computeStatistics);
	uint8_t getMarkedRowsPercentage();
	void setMarkedRowsPercentage(uint8_t markedRowsPercentage);
	vector<char> getOriginalSymbols();
	void setOriginalSymbols(const vector<char> originalSymbols);
	Statistics* getStatistics();
	ScrambledSuperAlphabet* getSuperAlphabet();
	void setSuperAlphabet(ScrambledSuperAlphabet* superAlphabet);
	uint8_t getSuperAlphabetOrder() const;
	void setSuperAlphabetOrder(uint8_t superAlphabetOrder);
	uint64_t getSuperBucketSize();
	void setSuperBucketSize(uint64_t superBucketSize);
	uint64_t getSuperbucketNumber(uint64_t bucketAbsoluteNumber);
	u8 *getEncryptionKey();
    void setEncryptionKey(u8* encryptionKey);
    void setEncryptionKeyFilePath(string encryptionKeyFilePath);
    string getEncryptionKeyFilePath();

    void setNumberOfThreads(uint32_t numberOfThreads);
    uint32_t getNumberOfThreads();

    void build(string &originalText,string fileName);
    void load(string fileName);
    void readText();
    void close();


    void exactMatch(string pattern);
    //vector<Match*> *exactMatchOld(string pattern);
    //vector<Match*> *exactMatch(SuperPattern *superPattern);
    vector<uint64_t> *locateOccurrences(vector<Match*> *matches);
    uint64_t countOccurrences(vector<Match*> *matches);
    string extractSnippet(uint64_t from, uint64_t to);
    friend class Bucket;
    friend class SuperBucket;
    friend class BucketEncoder;

protected:
    void buildInMemory(string &originalText);
    void initializeStatistics(uint64_t originalTextLength);
    void finalizeStatistics();

    string trimInitialAndFinalNs(string &originalText);
    void computeOriginalSymbols(string &originalText);
    void buildSuperAlphabet();
    void buildSuperText(string &originalText);
    void remapBWT();
    void buildSuperBuckets();
    void buildBuckets();
    void encodeBuckets();
    void buildBuckets(uint64_t superBucketNumber);
    Bucket *allocateBucket(uint64_t currentBucket,uint32_t *sbwt);
    void writeIndex();
    virtual void writeHeader();
    virtual void readHeader();

    void writeBuckets();
    virtual void createSpaceForSuspendedInfos();
    virtual void writeSuspendedInformations();
    virtual void readSuspendedInformations();

    inline void checkSuperBucketLoadedOld(BitReader *bitReader,uint64_t currentSuperBucket);
    inline void checkBucketLoadedOld(BitReader *bitReader,uint64_t currentBucket);
    inline void checkSuperBucketLoaded(BitReader *bitReader,uint64_t sn);
    inline void checkBucketLoaded(BitReader *bitReader,uint64_t bn);
    CryptoParams *getCryptoParams(uint64_t bucketNumber,uint64_t bucketLength,uint32_t bucketAlphaSize);
    inline vector<SuperPattern*> computeSuperPatterns(string pattern);
    inline vector<Match*> *exactMatch(SuperPattern *superPattern,Interval *interval,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator);
    void exactMatch(CompactAlphabetInterval *interval,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator);
    inline void findWholePatternMatches(Match *match,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator);

    //vector<Match*> *optimizedExactMatch(SuperPattern *superPattern,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator);
    string extractSuperCharacter(uint64_t position,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator);
    inline void getBwtPositionCoordinates(uint64_t bwtPosition,uint64_t *superBucketAbsoluteNumber,uint64_t *bucketAbsoluteNumber,uint64_t *bucketOffset);
    Interval *backwardStep(uint32_t c,Interval *si,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator);
    Interval *backwardStepIfMatch(uint64_t saIndex,string pattern,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator);
    uint64_t getRowPosition(uint64_t i,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator);
    string extractSnippet(uint64_t stFrom,uint8_t stFromOffset,uint64_t stTo,uint8_t stToOffset);
    string extractAll();
    vector<uint64_t> *locateOccurrences(Match *match);

    void reconstructMarkedRows();

	vector<char> originalSymbols;
	uint8_t superAlphabetOrder;
	ScrambledSuperAlphabet *superAlphabet;
	uint8_t markedRowsPercentage; //percentage of marked chars in BWT (chars whose position is linked to their relative text position)
	uint64_t bucketSize;
	uint64_t superBucketSize;
	uint64_t bucketsNumber;        //number of buckets
	uint64_t superBucketsNumber;	  //number of super-buckets
	SuperBucket **superBuckets;       //array of super buckets
	uint64_t *superBucketsStarts;        //starting positions of super buckets in the compressed text
	Bucket **buckets; 				//array of buckets
	uint64_t *bucketsStarts;       		//starting positions of buckets in the compressed text

	Statistics *statistics;
	bool computeStatistics;

	u8 *encryptionKey;
	string encryptionKeyFilePath;
	EncryptionManager *encryptionManager;


	/**
	 * Chromosomes reference sequences start with a very high number of N characters and
	 * also end with an high number of N characters: the suppression of these characters from BWT leads
	 * to two important advantages:
	 * 1) a better compression rate;
	 * 2) an higher BWT computing speed;
	 * Instead of effectively storing these characters, the index will take them into account
	 * simply storing the number of initial and final N characters
	 */
	uint64_t initialNCharacters;
	uint64_t finalNCharacters;


	uint32_t *text;		 //original text in the super-alphabet space
	uint32_t *debugText; //original super text, used only for debugging purposes

	uint64_t textLength;  //length of text in terms of super-alphabet characters
	uint32_t compactAlphabetSize;    //size of compact super-alphabet, i.e. number of super-characters really appearing in the super-text
	boost::dynamic_bitset<> *charactersBitSet; //identifies the compact alphabet's symbols
	uint64_t *precedingCharactersOccurrences; //number of the characters preceding each alphabet's character within the text

	sdsl::int_vector<0> charactersRemappings;
	unordered_map<uint32_t,uint32_t> decodeMap;


	FastBWTransformer *bwTransformer;
	uint32_t *sbwt;
	uint64_t lastCharacterPosition;

	uint32_t markingRate; //marking rate of suffix array (the markedRows array will contain an element for each position which is multiple of markingRate)
	uint64_t markedRowsNumber;  //number of marked rows
	uint32_t markedRowsKeysNeededBits;  //bits needed to store a marked row key
	uint32_t markedRowsValuesNeededBits; //bits needed to store a marked row value
	map<uint64_t,uint64_t>  *markedRows;  // contains, for each marked suffix array row (key),
		                       // its position within the text (value)
	uint64_t *markedRowsBuckets;   //contains, for each marked position, the bucket storing the marked row
	  	  	  	  	  	  	  //it's used in the extract operation
    						  //of the suffix corresponding to the i^th marked row

	uint32_t maximumBucketAlphaSize;
	BitWriter *bitWriter;				//bit writer
	BitReader *bitReader; 				//bit reader
	string fileName;  //absolute path of the index file on disk

	Bookmark *startSuspendedInfos;  //start byte of suspended infos within the header
	Bookmark *superBucketsPositionsStart;//start byte of super buckets positions in the compressed file
	Bookmark *bucketsPositionsStart;//start byte of buckets positions list in the compressed file
	Bookmark *characterOccurrencesStart; //start byte of character occurrences list
	Bookmark *markedRowsBucketsStart; //start byte of the marked rows buckets table

	bool tryToStoreCompressedBitmaps=false;
	bool encryptCumulativeFrequencies=true;
	uint64_t numberOfLoadedSuperBuckets=0;
	uint64_t numberOfLoadedBuckets=0;
	bool *bucketHeaderLoadingInProgress;
	int64_t *bucketTextLoadingInProgress;  //its elements are -1 if text loading isn't in progress,
									   //otherwise the untilTo position of the current loading operation
	bool *superBucketLoadingInProgress;
	//MUTEX used to build super-text in parallel
	mutex superTextBuildMutex;

	//MUTEXES used to implement a multi-threaded access to index structures
	mutex superBucketLoadingMutex;
	mutex bucketLoadingMutex;
	mutex bucketLoadingInProgressMutex;
	mutex getCryptoParamsMutex;
	mutex compactAlphabetIntervalsMutex; //protects nextCompactAlphabetIntervalToProcess
	mutex allSuperPatterMatchesMutex; //protects allSuperPatternMatches


	mutex readingProtectionMutex; //PROVA

	vector<CompactAlphabetInterval*> allCompactAlphabetIntervals;
	uint32_t nextCompactAlphabetIntervalToProcess;

	vector<Match*> allSuperPatternMatches;

	RandomGenerator *bucketEncryptionRandomGenerator;
	PRG *randomGenerators;
	PBR *bitReaders;


	uint32_t numberOfThreads;

	bool loadWholeFileInMemory;
	char *wholeFileBuffer;
	uint64_t wholeFileSize;

	class SuperTextBuilder {
	public:

		EFMIndex&  index;
		uint64_t firstTextPosition;
		uint64_t lastTextPosition;
		string &originalText;


		SuperTextBuilder(EFMIndex &index,string &originalText,uint64_t firstTextPosition,uint64_t lastTextPosition):index(index),originalText(originalText){
			this->firstTextPosition=firstTextPosition;
			this->lastTextPosition=lastTextPosition;
		}

		void run(){
			uint8_t k=index.superAlphabetOrder;

			for (uint64_t i=firstTextPosition;i<=lastTextPosition;i++){
				uint32_t symbolOrder=index.superAlphabet->getCode(originalText,i*k);
				index.text[i]=symbolOrder;
				/*
				index.superTextBuildMutex.lock();
				if ((*(index.charactersBitSet))[symbolOrder]==0){
					index.compactAlphabetSize++;
					index.charactersBitSet->set(symbolOrder);
				};
				index.precedingCharactersOccurrences[symbolOrder]++;
				index.superTextBuildMutex.unlock();
				*/
			}
		}

	};


	class BucketEncoder {
	public:
		uint64_t firstBucket;
		uint64_t lastBucket;
		RandomGenerator *bucketEncryptionRandomGenerator;
		EFMIndex *index;
		u8 *polySeed;

		BucketEncoder(EFMIndex *index,uint64_t firstBucket,uint64_t lastBucket,u8 *polySeed){
			this->index=index;
			this->firstBucket=firstBucket;
			this->lastBucket=lastBucket;
			this->polySeed=polySeed;
			bucketEncryptionRandomGenerator=NULL;
		}

		~BucketEncoder(){
			if (bucketEncryptionRandomGenerator!=NULL)
				delete bucketEncryptionRandomGenerator;
		}

		void run(){
			std::chrono::time_point<std::chrono::system_clock> startTime,endTime;
			startTime=std::chrono::system_clock::now();
			bucketEncryptionRandomGenerator=new RandomGenerator(131072,polySeed);
			for (int64_t bucketNumber=firstBucket;bucketNumber<=lastBucket;bucketNumber++){
				Bucket *bucket=index->buckets[bucketNumber];

				//GET CRYPTOPARAMS
				uint32_t* bucketKeyStream=new uint32_t[bucket->length];
				//Initialize a pseudo-casual number generator
				uint32_t maxValue=bucket->bucketAlphaSize;  //this upper bound is included in the range of generated numbers
				bucketEncryptionRandomGenerator->populateKeyStreamBuffer(bucketNumber,maxValue,0,bucket->length-1);
				for (uint64_t i=0;i<bucket->length;i++)
					bucketKeyStream[i]=bucketEncryptionRandomGenerator->nextInt();
				bucket->cryptoParams=new CryptoParams(bucketKeyStream);

				//ENCODE BUCKET TEXT

				EncodingResult *er=bucket->encodeText();
				bucket->bucketText=er->sequence;
				bucket->compressedLength=er->realLength;
			}

			endTime=std::chrono::system_clock::now();
			double totalTime=(double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()/(double)1000;
			#if INDEX_DEBUG_LEVEL > 0
			index->coutMutex.lock();
			cout << "\t [" << firstBucket << "," << lastBucket << "]" << "\t totalTime: " << totalTime  << endl;
			index->coutMutex.unlock();
			#endif
		}


	};


	//bool executePatternMatching;
	//bool indexClosingInProgress;
	//bool patternMatchingCompleted;
	mutex coutMutex;
	//vector<bool> terminatedPatternMatching;
	//vector<bool> notifiedCompletion;
	std::condition_variable matcherThreadsWakeUpCV;
	std::mutex matcherThreadsWakeUpMutex;
	std::vector<std::thread*> matcherThreads;
	std::condition_variable mainThreadWakeUpCV;
	std::mutex mainThreadWakeUpMutex;


	bool allMatchingThreadsFinished(){
		uint32_t i=0;
		uint32_t nt=matcherThreads.size();
		bool finished=true;
		while (i<nt && finished){
			finished = finished && matchers[i]->finishedPatternMatching;
			i++;
		}
		return finished;
	}

	/*
	bool allMatchingThreadsNotified(){
		uint32_t i=0;
		uint32_t nt=matcherThreads.size();
		bool allNotified=true;
		while (i<nt && allNotified){
			allNotified = allNotified && notifiedCompletion[i];
			i++;
		}
		return allNotified;
	}

	*/

	class IntervalMatcher {
		public:
			EFMIndex&  index;
			u8 *polySeed;
			string filename;//filename
			RandomGenerator *bucketEncryptionRandomGenerator;
			BitReader *bitReader;
			uint32_t threadNumber;
			bool doTerminate;
			bool doPatternMatching;
			bool finishedPatternMatching;

			/*IntervalMatcher(EFMIndex &index, u8 *polySeed, string filename):index(index){
				this->polySeed=polySeed;
				this->filename=filename;
			};*/

			IntervalMatcher(EFMIndex &index, uint32_t threadNumber,BitReader *bitReader, RandomGenerator *bucketEncryptionRandomGenerator):index(index){
				this->bitReader=bitReader;
				this->bucketEncryptionRandomGenerator=bucketEncryptionRandomGenerator;
				this->threadNumber=threadNumber;
				polySeed=NULL;
				doTerminate=false;
				doPatternMatching=false;
			};

			void run(){
				do{
					{
						//acquire work to do
						std::unique_lock<std::mutex> lk(index.matcherThreadsWakeUpMutex);
						index.matcherThreadsWakeUpCV.wait(lk, [&]{
								return doPatternMatching || doTerminate;
					     	 });
					}
					if (doPatternMatching){
						//DO PATTERN MATCHING UNTIL THERE ARE COMPACT ALPHABET INTERVAL TO ANALYZE
						uint32_t n=index.allCompactAlphabetIntervals.size();
						index.compactAlphabetIntervalsMutex.lock();
						uint32_t ti=index.nextCompactAlphabetIntervalToProcess;
						index.nextCompactAlphabetIntervalToProcess++;
						index.compactAlphabetIntervalsMutex.unlock();

						while (ti < n)	{
							//gets from allCompactAlphabetIntervals the next compact alphabet interval
							index.exactMatch(index.allCompactAlphabetIntervals[ti],bitReader,bucketEncryptionRandomGenerator);

							index.compactAlphabetIntervalsMutex.lock();
							ti=index.nextCompactAlphabetIntervalToProcess;
							index.nextCompactAlphabetIntervalToProcess++;
							index.compactAlphabetIntervalsMutex.unlock();
						}

						doPatternMatching=false;
						index.mainThreadWakeUpMutex.lock();
						finishedPatternMatching=true;
						index.mainThreadWakeUpMutex.unlock();
						index.mainThreadWakeUpCV.notify_one();
					}

				} while (!doTerminate);

			}
		};


	typedef IntervalMatcher* PIM;

	PIM *matchers;

};





} /* namespace std */

#endif /* EFMINDEX_H_ */
