/*
 * EFMIndex.cpp
 *
 *  Created on: 19/dic/2014
 *      Author: fernando
 */
#include <chrono>
#include <set>
#include <cmath>
#include <thread>
#include <algorithm>
#include "EFMIndex.h"
#include "EncryptionManager.h"
#include "Utils.h"
#include "FileBitWriter.h"
#include "FileBitReader.h"
#include "MemoryBitReader.h"
#include "BitSetRLEEncoder.h"
#include "Bookmark.h"
#include "Common.h"
#include "RandomGenerator.h"


namespace std { /* namespace std */

typedef pair<uint64_t,uint64_t> IntIntPair;

typedef pair<uint64_t,uint64_t> UIntIntPair;

ostream& operator<<(ostream& os, const Interval& i){
	os << "[" <<i.sp << "," << i.ep << "]";
	return os;
}

ostream& operator<<(ostream& os, const CompactAlphabetInterval& i){
	os << "[" <<i.sp << "," << i.ep << "]";
	return os;
}


ostream& operator<<(ostream& os, const SuperPatternCharacter& sc) {
	os << sc.text;
	return os;
}

ostream& operator<<(ostream& os, const SuperPattern& sp){
	for (uint32_t sc=0;sc<sp.length;sc++){
		os << sp.text[sc];
		if (sc<sp.length-1)
			os << "-";
		return os;
	}
}

ostream& operator<<(ostream& os, const Match& m){
	os << *m.suffixArrayInterval;
	return os;
}

Match::Match(){
	occurrences=NULL;
}

Match::~Match() {
	if (occurrences!=NULL)
		delete occurrences;
}



EFMIndex::EFMIndex() {
	bucketHeaderLoadingInProgress=NULL;
	bitWriter=NULL;
	bitReader=NULL;
	superBucketsStarts=NULL;
	bucketsStarts=NULL;

	encryptionKey=NULL;
	encryptionManager=NULL;
	bucketEncryptionRandomGenerator=NULL;
	text=NULL;
	debugText=NULL;
	charactersBitSet=NULL;
	precedingCharactersOccurrences=NULL;
	bwTransformer=NULL;
	sbwt=NULL;
	startSuspendedInfos=NULL;
	superBucketsPositionsStart=NULL;
	bucketsPositionsStart=NULL;
	characterOccurrencesStart=NULL;
	markedRowsBucketsStart=NULL;
	buckets=NULL;
	superBuckets=NULL;

	numberOfThreads = std::thread::hardware_concurrency();  //default


	statistics=new Statistics();
	loadWholeFileInMemory=true;
	wholeFileBuffer=NULL;

}

EFMIndex::~EFMIndex() {
	if (superBucketsStarts!=NULL)
		delete[] superBucketsStarts;
	if (bucketsStarts!=NULL)
		delete[] bucketsStarts;
	if (statistics!=NULL)
		delete statistics;
	if (encryptionManager!=NULL)
		delete encryptionManager;
	if (text!=NULL)
		delete text;
	if (charactersBitSet!=NULL)
		delete charactersBitSet;
	if (precedingCharactersOccurrences!=NULL)
		delete[] precedingCharactersOccurrences;
	if (bwTransformer!=NULL)
		delete bwTransformer;
	if (sbwt!=NULL)
		delete[] sbwt;
	if (markedRowsBuckets!=NULL)
		delete[] markedRowsBuckets;
	if (bitReader!=NULL)
		delete bitReader;
	if (bitWriter!=NULL)
		delete bitWriter;
	if (bucketHeaderLoadingInProgress!=NULL)
		delete[] bucketHeaderLoadingInProgress;
	if (startSuspendedInfos!=NULL)
		delete startSuspendedInfos;
	if (superBucketsPositionsStart!=NULL)
		delete superBucketsPositionsStart;
	if 	(bucketsPositionsStart!=NULL)
		delete bucketsPositionsStart;
	if 	(characterOccurrencesStart!=NULL)
		delete characterOccurrencesStart;
	if	(markedRowsBucketsStart!=NULL)
		delete markedRowsBucketsStart;
	if (buckets!=NULL){
		for (uint64_t i=0;i<bucketsNumber;i++){
			delete buckets[i];
			buckets[i]=NULL;
		}
		delete[] buckets;
	}
	if (superBuckets!=NULL){
		for (uint64_t i=0;i<superBucketsNumber;i++){
			delete superBuckets[i];
			superBuckets[i]=NULL;
		}
		delete[] superBuckets;
	}

	if (wholeFileBuffer!=NULL)
		delete[] wholeFileBuffer;

}

uint64_t EFMIndex::getBucketSize() const {
	return bucketSize;
}

uint64_t EFMIndex::getSuperbucketNumber(uint64_t bucketAbsoluteNumber){
		return bucketAbsoluteNumber/(superBucketSize/bucketSize);
}

inline void EFMIndex::getBwtPositionCoordinates(uint64_t bwtPosition,uint64_t *superBucketAbsoluteNumber,uint64_t *bucketAbsoluteNumber,uint64_t *bucketOffset){
	*bucketAbsoluteNumber= bwtPosition/bucketSize;
	*bucketOffset= bwtPosition%bucketSize;
	*superBucketAbsoluteNumber=*bucketAbsoluteNumber/(superBucketSize/bucketSize);
}

void EFMIndex::setBucketSize(uint64_t bucketSize) {
	this->bucketSize = bucketSize;
}

uint64_t EFMIndex::getBucketsNumber() {
	return bucketsNumber;
}

uint64_t EFMIndex::getSuperBucketsNumber() {
	return superBucketsNumber;
}

uint64_t EFMIndex::getNumberOfLoadedBuckets() {
	return numberOfLoadedBuckets;
}

uint64_t EFMIndex::getNumberOfLoadedSuperBuckets() {
	return numberOfLoadedSuperBuckets;
}

bool EFMIndex::isComputeStatistics() const {
	return computeStatistics;
}

void EFMIndex::setComputeStatistics(bool computeStatistics) {
	this->computeStatistics = computeStatistics;
}

uint8_t EFMIndex::getMarkedRowsPercentage() {
	return markedRowsPercentage;
}

void EFMIndex::setMarkedRowsPercentage(uint8_t markedRowsPercentage) {
	this->markedRowsPercentage = markedRowsPercentage;
}

vector<char> EFMIndex::getOriginalSymbols(){
	return originalSymbols;
}

void EFMIndex::setOriginalSymbols(const vector<char> originalSymbols) {
	this->originalSymbols = originalSymbols;
}

Statistics* EFMIndex::getStatistics() {
	return statistics;
}

ScrambledSuperAlphabet* EFMIndex::getSuperAlphabet() {
	return superAlphabet;
}

void EFMIndex::setSuperAlphabet(ScrambledSuperAlphabet* superAlphabet) {
	this->superAlphabet = superAlphabet;
}

uint8_t EFMIndex::getSuperAlphabetOrder() const {
	return superAlphabetOrder;
}

void EFMIndex::setSuperAlphabetOrder(uint8_t superAlphabetOrder) {
	this->superAlphabetOrder = superAlphabetOrder;
}

uint64_t EFMIndex::getSuperBucketSize(){
	return superBucketSize;
}

void EFMIndex::setSuperBucketSize(uint64_t superBucketSize) {
	this->superBucketSize = superBucketSize;
}


u8 *EFMIndex::getEncryptionKey(){
	return encryptionKey;
}

void EFMIndex::setEncryptionKey(u8* encryptionKey){
	this->encryptionKey=encryptionKey;
}

void EFMIndex::setEncryptionKeyFilePath(string encryptionKeyFilePath){
	this->encryptionKeyFilePath=encryptionKeyFilePath;
}

string EFMIndex::getEncryptionKeyFilePath(){
	return encryptionKeyFilePath;
}

void EFMIndex::setNumberOfThreads(uint32_t numberOfThreads){
	uint32_t nt = std::thread::hardware_concurrency();
	if (numberOfThreads<nt)
		this->numberOfThreads=numberOfThreads;
	else{
		this->numberOfThreads=nt;
	}
}

uint32_t EFMIndex::getNumberOfThreads(){
	return this->numberOfThreads;
}

void EFMIndex::build(string &originalText,string fileName){
		buildInMemory(originalText);

		std::chrono::time_point<std::chrono::system_clock> startTime,endTime;

		if (computeStatistics)
			startTime=std::chrono::system_clock::now();


		this->fileName=fileName;
		bitWriter=new FileBitWriter(fileName);
		bitWriter->setBlockSize(2*bucketSize);


		writeIndex();

		if (computeStatistics){
			endTime=std::chrono::system_clock::now();
			double totalTime=(double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()/(double)1000;
		    statistics->indexSaveTime =(int)(totalTime);
		}
}


/**
 * Builds the index structure in memory
 * @param originalText
 * @throws Exception
 */
void EFMIndex::buildInMemory(string &originalText) {
	#if INDEX_DEBUG_LEVEL > 0
	cout << "Original text length: " << originalText.length() <<endl;
	#endif

	if (computeStatistics)
		initializeStatistics(originalText.length());

	initialNCharacters=0;
	finalNCharacters=0;
	//originalText=trimInitialAndFinalNs(originalText);

	if (originalSymbols.size()==0){
		std::chrono::time_point<std::chrono::system_clock> currentTime=std::chrono::system_clock::now();
		std:time_t printableTime= std::chrono::system_clock::to_time_t(currentTime);
		char mbstr[80];
		std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
		cout << "Computing original symbols, starting at "<< mbstr << endl;
		computeOriginalSymbols(originalText);
	}

	//Initializes the encryption manager
	//if the encryptionKey is not assigned then try to load it from the keyfile
	if (encryptionKey!=NULL)
		encryptionManager=new EncryptionManager(encryptionKey);
	else
		encryptionManager=new EncryptionManager(encryptionKeyFilePath);


	std::chrono::time_point<std::chrono::system_clock> currentTime=std::chrono::system_clock::now();
	std::time_t printableTime= std::chrono::system_clock::to_time_t(currentTime);
	char mbstr[80];
	std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
	cout << "Building super-alphabet, starting at " << mbstr<< endl;
	buildSuperAlphabet();

	currentTime=std::chrono::system_clock::now();
	printableTime= std::chrono::system_clock::to_time_t(currentTime);
	std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
	cout << "Building super-text, starting at " << mbstr << endl;
	buildSuperText(originalText);

	debugText=text;

	cout << "Computing scrambled BWT" << endl;
	//uint32_t cores = 24;


	uint32_t numberOfElementaryRanges=1000;
	if (superAlphabet->order > 3){
		if (originalText.length()>1000000000)
			numberOfElementaryRanges=1000;
	} else if (superAlphabet->order==3)
		numberOfElementaryRanges=100;
	else if (superAlphabet->order==2)
		numberOfElementaryRanges=10;
	#if INDEX_DEBUG_LEVEL > 0
	cout << "Number of elementary ranges: " << numberOfElementaryRanges <<endl;
	#endif
	originalText.resize(0);
	originalText.shrink_to_fit();

	bwTransformer=new FastBWTransformer(numberOfThreads,superAlphabet, numberOfElementaryRanges);
	if (markedRowsPercentage>0)
		markingRate=100/markedRowsPercentage;
	else
		markingRate=0;
	bwTransformer->computeBWT(text,textLength, markingRate);

	#if INDEX_DEBUG_LEVEL > 0
	  bwTransformer->verifyBWT();
	#endif

	sbwt=bwTransformer->getBWT();
	lastCharacterPosition=bwTransformer->getLastCharacterPosition();

	if (computeStatistics){
		statistics->sortTime=(uint64_t)bwTransformer->getSortingTime();
		statistics->bwtComputationTime=(uint64_t)bwTransformer->getBwtComputationTime();
	}

	currentTime=std::chrono::system_clock::now();
	printableTime= std::chrono::system_clock::to_time_t(currentTime);
	std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
	cout << "Building index structure,  starting at " << mbstr <<endl;

	//Computes the number of buckets and super-buckets
	superBucketsNumber=textLength/superBucketSize;
	bucketsNumber=(textLength + bucketSize-1)/bucketSize;

	if (markedRowsPercentage>0){
		markedRows=bwTransformer->getMarkedRows();
		markedRowsKeysNeededBits=Utils::int_log2(bucketSize-1);
		markedRowsValuesNeededBits=Utils::int_log2(textLength/markingRate);
		markedRowsBuckets=new uint64_t[textLength/markingRate + 1];
		for (auto it=markedRows->begin(); it!=markedRows->end(); ++it){
			 markedRowsBuckets[it->second]=it->first/bucketSize;
			 //cout << it->first << "," << it->second*markingRate << endl;
		}
	}

	remapBWT();

	//for (int64_t i=52424704; i<=52427866; i++)
	//	cout << "bwt[" << i << "]="<<sbwt[i] << endl;

	std::chrono::time_point<std::chrono::system_clock> startTime,endTime;
	if (computeStatistics)
	startTime=std::chrono::system_clock::now();


	buildSuperBuckets();

	buildBuckets();


	if (computeStatistics){
		endTime=std::chrono::system_clock::now();
	    double totalTime=std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
	   statistics->indexStructureBuildTime =(int)(totalTime);
	}


	#if MULTI_THREADED_BUCKET_ENCODING == 1
	currentTime=std::chrono::system_clock::now();
	printableTime= std::chrono::system_clock::to_time_t(currentTime);
	std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
	cout << "Encoding and encrypting buckets (second encryption stage), starting at "  << mbstr <<endl;
	if (computeStatistics)
		startTime=std::chrono::system_clock::now();
	encodeBuckets();
	if (computeStatistics){
		endTime=std::chrono::system_clock::now();
		double totalTime=std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
	   statistics->bucketEncodingTime =(int)(totalTime);
	}
	#endif

}


string EFMIndex::trimInitialAndFinalNs(string &originalText){
		bool foundFinalNonNCharacter=false;
		int64_t otl=originalText.length();
		int64_t i=otl-1;
		while (i>=0 && !foundFinalNonNCharacter){
			if (originalText[i]!='N')
				foundFinalNonNCharacter=true;
			else
				i--;
		}
		finalNCharacters=originalText.length()-(i+1);

		bool foundInitialNonNCharacter=false;
		i=0;
		while (i<otl && !foundInitialNonNCharacter){
			if (!foundInitialNonNCharacter && originalText[i]!='N'){
				foundInitialNonNCharacter=true;
				initialNCharacters=i;
			}
			else
				i++;
		}
		return originalText.substr(initialNCharacters,otl-finalNCharacters-initialNCharacters);
}



void EFMIndex::computeOriginalSymbols(string &originalText){
	set<char> chars;
	int64_t otl=originalText.length();
	//for (int64_t i=0;i<otl;i++)
	//	chars.insert(originalText[i]);

	uint32_t nt = numberOfThreads; //Create threads
	//uint32_t nt=24;
	int64_t perThreadCharacters = otl / nt;
	std::vector<std::thread*> threads;
	std::vector<OriginalSymbolsFinder*> finders;
	for (uint32_t t = 0; t < nt; t++) {
		int64_t firstTextPosition = t * perThreadCharacters;
		int64_t lastTextPosition = (t + 1) * perThreadCharacters - 1;
		if (t == nt - 1)
			lastTextPosition = otl - 1;
		OriginalSymbolsFinder* finder=new OriginalSymbolsFinder(originalText,firstTextPosition, lastTextPosition);

		finders.push_back(finder);
		threads.push_back(new std::thread(&OriginalSymbolsFinder::run, finder));
	}
	for (uint32_t i = 0; i < nt; i++){
		threads[i]->join();
		chars.insert(finders[i]->chars.begin(),finders[i]->chars.end());
	}

	originalSymbols.clear();
	if (chars.find('$')==chars.end())
		originalSymbols.push_back('$');

	originalSymbols.insert(originalSymbols.end(),chars.begin(),chars.end());
	#if INDEX_DEBUG_LEVEL > 0
	for (auto c:originalSymbols){
		cout << "\t" << c << endl;
	}
	#endif

}


void EFMIndex::buildSuperAlphabet(){
	std::chrono::time_point<std::chrono::system_clock> startTime,endTime;
	if (computeStatistics)
		startTime=std::chrono::system_clock::now();

	uint32_t superAlphabetCardinality=pow(originalSymbols.size(),superAlphabetOrder);

	superAlphabet=new ScrambledSuperAlphabet(originalSymbols,superAlphabetOrder,
								   encryptionManager->computeScramblingKey(superAlphabetCardinality));
	uint32_t allNCode=superAlphabet->getCode("NNNN");


	if (computeStatistics){
	   endTime=std::chrono::system_clock::now();
	   double totalTime=std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
	   statistics->superAlphabetBuildTime =(int)(totalTime);

	}

	#if INDEX_DEBUG_LEVEL > 0
	cout << "Super-alphabet size: " << superAlphabet->getSize() <<endl;
	#endif

}


/**
	 * Builds superText, i.e. the representation of the originalText in terms of the super-alphabet
	 * symbols, storing it in text variable.
	 * The absolute frequencies of super-alphabet symbols are stored in precedingCharactersOccurrences
	 * ,boolean value indicating that i_th super-character appears in text is stored in
	 * charactersBitmap[i] and real super-alphabet size is stored in compactSuperAlphabetSize
	 * Remaps super-alphabet to a compact one, filling the array characterRemappings, of size compactSuperAlphabetSize,
	 * containing remappings of super-alphabet symbols to the compact super-alphabet ones.
	 */
void EFMIndex::buildSuperText(string &originalText){
	std::chrono::time_point<std::chrono::system_clock> startTime,endTime;
	if (computeStatistics)
		startTime=std::chrono::system_clock::now();

	uint8_t k=superAlphabetOrder;

	uint64_t otl=originalText.length();
	//int64_t originalTextLength=originalText.length();

	//Transforms text into an array of Super-characters
	//computes the array length
	textLength=ceil((double)otl/(double)k) + 1;
	text=new uint32_t[textLength];

	compactAlphabetSize=0;

	charactersBitSet=new boost::dynamic_bitset<>();
	charactersBitSet->clear();
	charactersBitSet->resize(superAlphabet->getSize());


	uint32_t saSize=superAlphabet->getSize();
	precedingCharactersOccurrences=new uint64_t[saSize];
	for (uint32_t i=0;i<saSize;i++)
		precedingCharactersOccurrences[i]=0;


	uint32_t nt = numberOfThreads; //Create threads
	uint64_t stc=otl/k;  //super-text characters to split between several threads
	uint64_t perThreadSuperCharacters = stc / nt; //super-characters computed by each thread
	std::vector<std::thread*> threads;
	std::vector<SuperTextBuilder*> builders;
	for (uint32_t t = 0; t < nt; t++) {
		uint64_t firstTextPosition = t * perThreadSuperCharacters;
		uint64_t lastTextPosition = (t + 1) * perThreadSuperCharacters - 1;
		if (t == nt - 1)
			lastTextPosition = stc - 1;
		SuperTextBuilder* builder=new SuperTextBuilder(*this,originalText,firstTextPosition, lastTextPosition);
		builders.push_back(builder);
		threads.push_back(new std::thread(&SuperTextBuilder::run, builder));
	}
	for (uint32_t i = 0; i < nt; i++){
		threads[i]->join();
	}


	for (uint64_t i=0;i<stc;i++){
			uint32_t symbolOrder=text[i];
			if ((*(charactersBitSet))[symbolOrder]==0){
				compactAlphabetSize++;
				charactersBitSet->set(symbolOrder);
			};
			precedingCharactersOccurrences[symbolOrder]++;

	}


	uint8_t remainder=otl%k;
	if (remainder>0){
		string paddedSuperCharacterValue=originalText.substr((otl/k)*k);
		for (uint8_t i=remainder;i<k;i++)
			paddedSuperCharacterValue+="$";
		uint32_t fillingSuperCharacter=superAlphabet->getCode(paddedSuperCharacterValue);
		text[textLength-2]=fillingSuperCharacter;
		compactAlphabetSize++;
		//charactersBitmap[fillingSuperCharacter]=true;
		charactersBitSet->set(fillingSuperCharacter);
		precedingCharactersOccurrences[fillingSuperCharacter]++;
	}


	uint32_t allDollarsSymbolOrder=0;
	text[textLength-1]=allDollarsSymbolOrder;
	compactAlphabetSize++;
	//charactersBitmap[allDollarsSymbolOrder]=true;
	charactersBitSet->set(allDollarsSymbolOrder);
	precedingCharactersOccurrences[allDollarsSymbolOrder]++;


	//Remaps super-alphabet on a compact one, filling charactersRemappings array
	charactersRemappings.resize(superAlphabet->getSize()); //number of elements
	charactersRemappings.width(Utils::int_log2(compactAlphabetSize-1)); //bits for each compact alphabet code

	//charactersRemappings=new CompactUnsignedArray(superAlphabet.getSize(), compactAlphabetSize-1);
	uint32_t firstFreeSpace=0;
	for (uint32_t i=0;i<superAlphabet->getSize();i++){
		if ((*charactersBitSet)[i]){
			charactersRemappings[i]=firstFreeSpace;
			precedingCharactersOccurrences[firstFreeSpace]=precedingCharactersOccurrences[i];
			firstFreeSpace++;
		}
	}

	//At this point only fist alphaSize elements of precedingCharactersOccurrences are significant
	//And now let precedingCharactersOccurrences[i] contain the number of characters in super-text
	//which are alphabetically smaller than character of order i
	uint64_t appo=0,sum=0;
	for (uint32_t i=0;i<=compactAlphabetSize;i++){
		appo=precedingCharactersOccurrences[i];
		precedingCharactersOccurrences[i]=sum;
		sum +=appo;
	}

	//Collects statistics
	if (computeStatistics){
		endTime=std::chrono::system_clock::now();
		double totalTime=std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
		statistics->textLength=textLength;
		statistics->superTextBuildTime =(int)totalTime;
	}

	#if INDEX_DEBUG_LEVEL > 0
	cout << "Super-text length: " << textLength <<endl;
	#endif

}


void EFMIndex::remapBWT(){
	for (uint64_t i=0;i<textLength;i++)
		sbwt[i]=charactersRemappings[sbwt[i]];
}


void EFMIndex::buildSuperBuckets(){
	//allocates super buckets data structures
	if (textLength%superBucketSize>0)
		superBucketsNumber++;
	superBucketsStarts=new uint64_t[superBucketsNumber];
	superBuckets=new SuperBucket*[superBucketsNumber];  //array of object pointers
	for (uint64_t i=0;i<superBucketsNumber;i++){
		SuperBucket *sb=new SuperBucket(this,i,compactAlphabetSize);
		superBuckets[i]=sb;
	}

	//initializes super buckets data structures
	uint64_t currentSuperBucketIndex=0;
	uint64_t currentSuperBucketUpperLimit=superBucketSize;
	SuperBucket *currentSuperBucket=superBuckets[0];

	for (uint64_t i=0;i<textLength;i++){
		if (i==currentSuperBucketUpperLimit){
			currentSuperBucketIndex++;
			currentSuperBucketUpperLimit+=superBucketSize;
			currentSuperBucket=superBuckets[currentSuperBucketIndex];
		}

		//initializes previousSuperbucketsOccurrences of current super bucket and
		// characters map only if RemappingLevel > none
		uint32_t currentChar=sbwt[i];

		if ((*currentSuperBucket->charactersBitSet)[currentChar]==false){
			(*currentSuperBucket->charactersBitSet).set(currentChar);
			//currentSuperBucket.occurrences[currentChar]=0;  non necessario
			currentSuperBucket->alphaSize++;
		}

		currentSuperBucket->previousSuperbucketsOccurrences[currentChar]++;
	}

	//for each compact alphabet character computes, for each super bucket, the number
	// of occurrences in previous super buckets
	uint64_t occ;
	for (uint32_t k=0;k<compactAlphabetSize;k++){
		occ=0;
		for (uint64_t i=0;i<superBucketsNumber;i++){
			uint64_t temp=superBuckets[i]->previousSuperbucketsOccurrences[k];
			superBuckets[i]->previousSuperbucketsOccurrences[k]=occ;
			if (occ > superBuckets[i]->maxOccurrence)
				superBuckets[i]->maxOccurrence=occ;
			occ+=temp;
		}
	}

	//builds the character map for each super bucket
	for (uint64_t i=0;i<superBucketsNumber;i++){
		SuperBucket *sb=superBuckets[i];
		sb->charactersRemappings.resize(compactAlphabetSize);  //number of elements
		sb->charactersRemappings.width(Utils::int_log2(sb->alphaSize-1)); //bits for each super-bucket character code
		uint32_t temp=0;
		for (uint32_t k=0;k<compactAlphabetSize;k++)
			if ((*sb->charactersBitSet)[k])
				sb->charactersRemappings[k]=temp++;
	}


}


void EFMIndex::buildBuckets(){
	//allocates buckets data structures
	bucketsStarts=new uint64_t[bucketsNumber];
	buckets=new Bucket*[bucketsNumber];

	for (uint64_t num=0;num<superBucketsNumber;num++)
		buildBuckets(num);
}


void EFMIndex::encodeBuckets(){
	uint32_t nt = std::thread::hardware_concurrency(); //Allocate as many threads as the number of the microprocessor's cores to use

	int64_t perThreadBuckets = bucketsNumber / nt; //buckets encoded by each thread
	std::vector<std::thread*> threads;
	std::vector<BucketEncoder*> encoders;
	for (uint32_t t = 0; t < nt; t++) {
		int64_t firstBucket = t * perThreadBuckets;
		int64_t lastBucket = (t + 1) * perThreadBuckets - 1;
		if (t == nt - 1)
			lastBucket = bucketsNumber - 1;
		BucketEncoder* encoder=new BucketEncoder(this,firstBucket, lastBucket,encryptionManager->polySeed);
		encoders.push_back(encoder);
		threads.push_back(new std::thread(&BucketEncoder::run, encoder));
		//encoder->run();
	}
	for (uint32_t i = 0; i < nt; i++){
		threads[i]->join();
	}
	for (uint32_t i = 0; i < nt; i++){
		delete encoders[i];
	}
}



/**
 * Builds the buckets of a superbucket
 * @param num
 * @throws IOException
 */
void EFMIndex::buildBuckets(uint64_t superBucketNumber){
	uint64_t j;
	uint32_t c;
	uint64_t len;
	uint32_t neededBits;
	uint64_t *bucketOccurrences;
	uint64_t maxOccurrence;

	SuperBucket *sb=superBuckets[superBucketNumber];

	//starting position of the super bucket
	uint64_t sbStart=superBucketNumber*superBucketSize;

	//size of the super bucket's alphabet (needed to encode symbols in this super bucket
	//using a number of bits as low as possible)
	uint32_t sbAlphaSize=sb->alphaSize;
	uint64_t sbEnd=min(sbStart+ superBucketSize, textLength);
	uint64_t currentBucket=sbStart/bucketSize;

	//initializes the occurrences vector
	bucketOccurrences=new uint64_t[sbAlphaSize];
	for (uint32_t i=0;i<sbAlphaSize;i++)
		bucketOccurrences[i]=0;
	//maximum number of occurrences,
	//from which the number of bits needed to encode the occurrences will be calculated
	maxOccurrence=0;


	//in the following cycle the "start" variable contains the start index of the
	//current bucket within the BWT
    for(uint64_t start=sbStart; start<sbEnd; start+=bucketSize, currentBucket++) {

        //length of the bucket
	    len = min(bucketSize, sbEnd-start);

		bool isOdd = false;
		if((currentBucket%2 ==0) && (start != sbStart) && (currentBucket != bucketsNumber-1))
			isOdd = true; // decide quali bucket da rovesciare: per essi non verranno memorizzate le occorrenze dei caratteri nei bucket precedenti

		Bucket *bucket=allocateBucket(currentBucket,sbwt); //new Bucket(this,sb,currentBucket,sbwt, start, len, sbAlphaSize,isOdd);
	    buckets[currentBucket]=bucket;

	    if(start != sbStart)   // if this isn't the first bucket of this super bucket
			if(!isOdd){ // and if it is not odd
				//bits needed for each entry of the occurrences table (containing, for each super bucket's character, the number
				//of its occurrences in the previous buckets of this superbucket)
				neededBits=ceil(Utils::int_log2(maxOccurrence));
				if (computeStatistics)
					statistics->bucketOccurrencesTablesCumulativeSize+=(neededBits*sb->alphaSize);
				//clone the bucket occurrences vector, assigning the clone to previousBucketsOccurrences
				bucket->previousBucketsOccurrences=new uint64_t[sbAlphaSize];
				copy(bucketOccurrences,bucketOccurrences+sbAlphaSize,bucket->previousBucketsOccurrences);
				bucket->previousBucketsOccurrencesMaximum=maxOccurrence;
		}

	    for(j=0; j<len; j++) {      // update previousBucketsOccurrences and remap
	    		// do the remapping at super bucket level
	    		sbwt[start + j] = sb->charactersRemappings[sbwt[start+j]];
	    		bucketOccurrences[sbwt[start + j] ]++;   // used in the next bucket
	    		if (bucketOccurrences[sbwt[start + j]]>maxOccurrence)
	    			maxOccurrence=bucketOccurrences[sbwt[start + j]];
	    }

		//there is a remapping at bucket level and so we must compute
		//characters bit set and characters remapping
	    bucket->charactersBitSet=new boost::dynamic_bitset<>();
		bucket->charactersBitSet->resize(bucket->superBucketAlphaSize);
		for(j=0;j<bucket->length;j++) {              // compute local boolean map
			    c = sbwt[bucket->startPosition+j];   // remapped char
			    bucket->charactersBitSet->set(c);
		}

	    bucket->bucketAlphaSize=bucket->charactersBitSet->count();
	    bucket->buildCharactersMap();

	    for(j=0;j<bucket->length;j++) // remap bucket section of sbwt to use bucket remapped codes
			sbwt[bucket->startPosition+j]=bucket->charactersRemappings[sbwt[bucket->startPosition+j]];

		if (computeStatistics)
				statistics->bucketAlphabetAverageSize+=bucket->bucketAlphaSize;


		if (bucket->bucketAlphaSize > maximumBucketAlphaSize)
			maximumBucketAlphaSize=bucket->bucketAlphaSize;


		//if necessary, computes marking informations, storing for each marked row:
		//row position in bucket
		//corresponding text position (divided by markingRate)
		if (markedRowsPercentage>0){
			for(j=0;j<bucket->length;j++){
				auto it=markedRows->find(bucket->startPosition+j);
			    if (it!=markedRows->end())
				  bucket->markedRows[j]= it->second;
			}

			if (computeStatistics)
				statistics->bucketMarkedRowsTablesCumulativeSize+=(bucket->markedRows.size()*(markedRowsKeysNeededBits+markedRowsValuesNeededBits));
		}

	}
}


Bucket *EFMIndex::allocateBucket(uint64_t currentBucket,uint32_t *sbwt){
		uint64_t start=(currentBucket)*bucketSize;  //start character of bwt in this bucket
		uint64_t end=min((currentBucket+1)*bucketSize-1,textLength-1);
		uint64_t length=end-start+1 ;   //number of characters in bucket
		//find superblock
		uint64_t t=getSuperbucketNumber(currentBucket);
		SuperBucket *sb=superBuckets[t];
		//coutMutex.lock();
		//cout << pthread_self() << " sbAlphaSize " << sb->alphaSize <<  " " << endl;
		//coutMutex.unlock();

		Bucket *b=new Bucket(this,sb,currentBucket,sbwt, start, length, sb->alphaSize);
		b->isOdd = false;
		b->isFirstInSuperBucket= (b->bucketNumber % (superBucketSize/bucketSize))==0;
		if((b->bucketNumber%2 ==0) &&  !b->isFirstInSuperBucket && (b->bucketNumber != bucketsNumber-1))
			b->isOdd = true; // i bucket odd sono rovesciati e non riportano le occorrenze
		buckets[currentBucket]=b;
		return b;
}



/**
 * Writes the entire index on a file or on a blob field using the bitWriter object
 * @throws Exception
 */
void EFMIndex::writeIndex() {
	std::chrono::time_point<std::chrono::system_clock> currentTime=std::chrono::system_clock::now();
	std::time_t printableTime= std::chrono::system_clock::to_time_t(currentTime);
	char mbstr[100];
	std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
	cout << "Writing index to disk, starting at " << mbstr <<endl;

	bitWriter->open();

	writeHeader();

	writeBuckets();

	//set file pointer to start byte of super buckets info
	Bookmark *eofBookmark=bitWriter->getBookmark();
	//write the suspended informations
	writeSuspendedInformations();
	//return to write at the end of the file
	bitWriter->gotoBookmark(eofBookmark);

	bitWriter->flush();
	uint64_t compressedSize=bitWriter->getBookmark()->getPosition();
	bitWriter->close();


	if (computeStatistics){
		statistics->compressedSize=compressedSize;
		finalizeStatistics();
	}

}


void EFMIndex::writeHeader(){
	uint8_t neededBits;

	//Writes the file type (Encrypted FM-index)
	bitWriter->writeUByte('E');  //Encrypted and Compressed
	bitWriter->writeUByte('G');  //Genomic index
	//Writes the original symbols (symbols of the primigenious alphabet)
	bitWriter->writeUByte(originalSymbols.size());
	for (uint32_t i=0;i<originalSymbols.size();i++)
		bitWriter->writeUByte(originalSymbols[i]);
	//Super-alphabet order
	bitWriter->writeUByte(superAlphabet->order);
	//Super-text length
	bitWriter->writeInt64(textLength);
	//Number of initial N characters
	bitWriter->writeInt64(initialNCharacters);
	//Number of final N characters
	bitWriter->writeInt64(finalNCharacters);
	//Position in BWT of the last character of the text
	bitWriter->writeInt64(lastCharacterPosition);
	//Superbuckets size
	bitWriter->writeInt64(superBucketSize);
	//Buckets size
	bitWriter->writeInt64(bucketSize);

	//marked rows percentage
	bitWriter->writeUByte(markedRowsPercentage);

	bitWriter->writeInt64(maximumBucketAlphaSize);

	startSuspendedInfos=bitWriter->getBookmark();

	// create free space for suspended Infos
	createSpaceForSuspendedInfos();


	boost::dynamic_bitset<> *bitsetToStore=charactersBitSet;
	bool storeInCompressedFormat=false;
	if (tryToStoreCompressedBitmaps){
		boost::dynamic_bitset<> *compressed=BitSetRLEEncoder::compress(charactersBitSet);
		storeInCompressedFormat=compressed->size() < charactersBitSet->size();
		if (storeInCompressedFormat)
				bitsetToStore=compressed;
	 }
	bitWriter->write(1,(storeInCompressedFormat==true)?1:0);
	//Size of characters bitmap
	bitWriter->writeInt64(bitsetToStore->size());
	//Characters Bitmap
	for (uint64_t i=0;i<bitsetToStore->size();i++)
		bitWriter->write(1, (*bitsetToStore)[i]);
	bitWriter->flush();

	characterOccurrencesStart=bitWriter->getBookmark();
	//For each character of compact alphabet, writes the occurrences of characters preceding it
	//computes number of needed bits for each occurrence (occurrences are cumulative and so this
	//number can be computed looking at the last character)
	neededBits=ceil(Utils::int_log2(precedingCharactersOccurrences[compactAlphabetSize]));
	bitWriter->writeUByte(neededBits);

	if (encryptCumulativeFrequencies){
		//Poly-alphabetic substitution of the cumulative occurrences, followed by a their transposition,
		//before writing them to the index file. This is necessary
		//to avoid that an attacker can obtain informations on the frequence distribution of the
		//scrambled alphabet characters
		uint64_t *encryptedOccurrences=new uint64_t[compactAlphabetSize+1];

		//POLYALPHABETIC SUBSTITUTION (respect to this, maximumFrequency+1 is the alphabet size: the alphabet, in fact,
		//is that of cumulative frequencies)
		//estimates the rounded up maximum frequency from needed bits
		//int maximumFrequency=(int)(pow(2.0, (double)neededBits))-1;
		uint64_t maximumFrequency=(int64_t)(pow(2.0, (double)neededBits))-1;
		uint64_t *polyAlphabeticKey=encryptionManager->computeCumulativeFreqsKeyStream(compactAlphabetSize+1, maximumFrequency);
		for (uint32_t i=0;i<=compactAlphabetSize;i++){
			encryptedOccurrences[i]=(precedingCharactersOccurrences[i]+polyAlphabeticKey[i])%(maximumFrequency+1);
			if (encryptedOccurrences[i]<0)
				cout << "overflow problem at" << i << endl;
		}

		//TRANSPOSITION
		uint32_t *permutation=encryptionManager->computeFrequenciesPermutation(compactAlphabetSize+1,maximumFrequency);

		for (uint32_t i=0;i<=compactAlphabetSize;i++)
			bitWriter->write(neededBits,encryptedOccurrences[permutation[i]]);
	}
	else
		for (uint32_t i=0;i<=compactAlphabetSize;i++)
			bitWriter->write(neededBits,precedingCharactersOccurrences[i]);
	bitWriter->flush();

	markedRowsBucketsStart=bitWriter->getBookmark();
	//if it's necessary, write a table storing, for each marked position, the bucket containing the marked row
	if (markedRowsPercentage>0){
		neededBits=(int)ceil(Utils::int_log2(bucketsNumber-1));
		for (uint64_t i=0;i<markedRows->size();i++)
			bitWriter->write(neededBits,markedRowsBuckets[i]);
		bitWriter->flush();
	}

	if (computeStatistics){
		statistics->alphabetSize=(uint32_t) pow((double)originalSymbols.size(), (double)superAlphabet->order);
		statistics->compactAlphabetSize=compactAlphabetSize;
		statistics->headerBitmapSize=bitsetToStore->size();
		statistics->headerOccurrencesTableSize=neededBits*compactAlphabetSize;
		statistics->buckets=bucketsNumber;;
		statistics->superbuckets=superBucketsNumber;
	}


	superBucketsPositionsStart=bitWriter->getBookmark();
	//Leaves space to accomodate super buckets informations (starting positions)
	for (uint64_t i=0;i<superBucketsNumber;i++)
		bitWriter->writeInt64(0);

	//Writes Super buckets informations and collects statistics (if necessary)
	int64_t alternativeFreqsCumulativeSize=0;

	for (uint64_t num=0;num<superBucketsNumber;num++){
		superBucketsStarts[num]=bitWriter->getBookmark()->getPosition();
		map<uint64_t,uint64_t> frequenciesMultiplicities;
		SuperBucket *sb=superBuckets[num];

		//writes the super bucket alphabet size and the bit set to file

		//writes the size of the super bucket's alphabet size
		//bitWriter->writeInt((int)sb.alphaSize);   //can be obtained directly from bitset cardinality

		//writes the super bucket's characters bitmap
		boost::dynamic_bitset<> *bitsetToStore=sb->charactersBitSet;
		storeInCompressedFormat=false;
		if (tryToStoreCompressedBitmaps){
			boost::dynamic_bitset<> *compressed=BitSetRLEEncoder::compress(sb->charactersBitSet);
			storeInCompressedFormat=compressed->size() < charactersBitSet->size();
			if (storeInCompressedFormat)
					bitsetToStore=compressed;
		}
		bitWriter->write(1,(storeInCompressedFormat==true)?1:0);
		//Size of characters bitmap
		bitWriter->writeInt64(bitsetToStore->size());
		//Characters Bitmap
		for (uint64_t i=0;i<bitsetToStore->size();i++)
			bitWriter->write(1, (*bitsetToStore)[i]);

		if (computeStatistics){
			statistics->superbucketAlphabetAverageSize+=sb->alphaSize;
			statistics->superbucketBitmapsCumulativeSize+=bitsetToStore->size();
		}


		if (computeStatistics){
			for (int32_t i=0;i<compactAlphabetSize;i++){
				int64_t mult=0;
				if (frequenciesMultiplicities.find(sb->previousSuperbucketsOccurrences[i]) != frequenciesMultiplicities.end() )
					mult=frequenciesMultiplicities[sb->previousSuperbucketsOccurrences[i]];
				mult++;
				frequenciesMultiplicities.insert(IntIntPair(sb->previousSuperbucketsOccurrences[i],mult));
			}
			//if (num%100==0)
			//	System.out.println(frequenciesMultiplicities.size());

		}


		//if the super bucket is not the first, writes a table containing, for each character of the
		//compact alphabet, the number of its occurrences in previous super buckets.
		if (num>0){
			neededBits=Utils::int_log2(sb->maxOccurrence);
			bitWriter->writeUByte(neededBits);
			for (uint32_t k=0;k<compactAlphabetSize;k++)
				bitWriter->write(neededBits,sb->previousSuperbucketsOccurrences[k]);

			if (computeStatistics){
				statistics->superbucketOccurrencesTablesCumulativeSize+=neededBits*compactAlphabetSize;
				int64_t maxInterval=0;
				int64_t last=-1;
				for ( auto it = frequenciesMultiplicities.begin(), end = frequenciesMultiplicities.end(); it != end; ++it ){
					int64_t f=it->first;
					int64_t interval=f-last;
					if (interval>maxInterval)
						maxInterval=interval;
					last=f;
				}

				if (frequenciesMultiplicities.size()*Utils::int_log2(maxInterval) + compactAlphabetSize*Utils::int_log2(frequenciesMultiplicities.size())<neededBits*compactAlphabetSize)
					alternativeFreqsCumulativeSize+=(frequenciesMultiplicities.size()*Utils::int_log2(maxInterval) +
							compactAlphabetSize*Utils::int_log2(frequenciesMultiplicities.size()));
				else
					alternativeFreqsCumulativeSize+=neededBits*compactAlphabetSize;

			}
		}
		bitWriter->flush();
	}
	//System.out.println("standard: "+statistics.superbucketOccurrencesTablesCumulativeSize/8);
	//System.out.println("\t\talternative: "+alternativeFreqsCumulativeSize/8);


	bucketsPositionsStart=bitWriter->getBookmark();
	//Leaves space to accomodate buckets informations (starting positions)
	for (uint64_t i=0;i<bucketsNumber;i++)
		bitWriter->writeInt64(0);

}


void EFMIndex::readHeader(){
	//Read the first two bytes and checks file type
	char firstByte=bitReader->getUByte();
	char secondByte=bitReader->getUByte();
	if (firstByte!='E' || secondByte!='G')
		throw runtime_error("File type is not Encrypted and Compressed Genomic Index");

	//Read the original symbols (symbols of the primigenious alphabet)

	uint8_t numberOfOriginalSymbols=bitReader->getUByte();
	originalSymbols.resize(numberOfOriginalSymbols);
	for (uint8_t i=0;i<numberOfOriginalSymbols;i++)
		originalSymbols[i]= bitReader->getUByte();

	//read super-alphabet order
	superAlphabetOrder=bitReader->getUByte();

	//build super-alphabet in memory
	buildSuperAlphabet();

	//super-text length
	textLength=bitReader->getInt64();

	//Number of initial N characters
	initialNCharacters=bitReader->getInt64();
	//Number of final N characters
	finalNCharacters=bitReader->getInt64();

	//Position in BWT of the last character of the text
	lastCharacterPosition=bitReader->getInt64();
	//Writes superbuckets size
	superBucketSize=bitReader->getInt64();
	//Buckets size
	bucketSize=bitReader->getInt64();

	//compute the number of buckets
	bucketsNumber = textLength/bucketSize;
	if (textLength%bucketSize>0)
		bucketsNumber++;

	//marked rows percentage
	markedRowsPercentage=bitReader->getUByte();

	maximumBucketAlphaSize=bitReader->getInt64();

	//suspended informations
	readSuspendedInformations();

	Bookmark *b=bitReader->getBookmark();

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


	compactAlphabetSize=charactersBitSet->count();


	//fill charactersRemappings array, needed to remap the whole super-alphabet on the compact one
	charactersRemappings.resize(superAlphabet->getSize()); //number of elements
	charactersRemappings.width(Utils::int_log2(compactAlphabetSize-1)); //bits for each compact alphabet code

	uint32_t firstFreeSpace=0;
	for (uint32_t i=0;i<superAlphabet->getSize();i++){
		if ((*charactersBitSet)[i]){
			charactersRemappings[i]=firstFreeSpace;
			firstFreeSpace++;
		}
	}

	//compute decoding map (needed to remap the compact super-alphabet backward to the whole super-alphabet)
	uint32_t code=0;
	for(uint32_t j=0; j<superAlphabet->getSize(); j++)
		if((*charactersBitSet)[j]){
		  decodeMap[code]= j;
		  code++;
		}

	b=bitReader->getBookmark();
	bitReader->gotoNextByteStart();

	//reads the number of bits used to encode the characters occurrences
	uint8_t neededBits=bitReader->getUByte();

	//For each character of compact alphabet, reads the occurrences of characters preceding it
	precedingCharactersOccurrences=new uint64_t[compactAlphabetSize+1];
	if (encryptCumulativeFrequencies){
		uint64_t *encryptedOccurrences=new uint64_t[compactAlphabetSize+1];
		//estimates the rounded up maximum frequency from needed bits
		//int maximumFrequency=(int)(pow(2.0, (double)neededBits))-1;
		int64_t maximumFrequency=(int64_t)(pow(2.0, (double)neededBits))-1;
		uint32_t *permutation=encryptionManager->computeFrequenciesPermutation(compactAlphabetSize+1,maximumFrequency);
		//invert the TRANSPOSITION
		for (uint32_t i=0;i<=compactAlphabetSize;i++)
			encryptedOccurrences[permutation[i]]=bitReader->read(neededBits);
		//Invert the POLYALPHABETIC SUBSTITUTION (respect to this, maximumFrequency+1 is the alphabet size: the alphabet
		//is that of cumulative frequencies)
		uint64_t *polyAlphabeticKey=encryptionManager->computeCumulativeFreqsKeyStream(compactAlphabetSize+1, maximumFrequency);
		int64_t pco;
		for (uint32_t i=0;i<=compactAlphabetSize;i++){
			pco=(encryptedOccurrences[i]-polyAlphabeticKey[i])%(maximumFrequency+1);
			if (pco<0)
				pco=pco + (maximumFrequency+1);
			precedingCharactersOccurrences[i]=pco;
		}
	}
	else{
		for (uint32_t i=0;i<=compactAlphabetSize;i++)
			precedingCharactersOccurrences[i]=bitReader->read(neededBits);
	}

	bitReader->gotoNextByteStart();

	if (markedRowsPercentage>0){
		markingRate=100/markedRowsPercentage;
		markedRowsNumber=textLength/markingRate + 1;
		markedRowsKeysNeededBits=Utils::int_log2(bucketSize-1);
		markedRowsValuesNeededBits=Utils::int_log2(textLength/markingRate-1);
		markedRowsBuckets=new uint64_t[markedRowsNumber];

		neededBits=ceil(Utils::int_log2(bucketsNumber-1));
		for (uint64_t i=0;i<markedRowsNumber;i++)
			markedRowsBuckets[i]=bitReader->read(neededBits);
		bitReader->gotoNextByteStart();
	}


	//Super buckets informations
	//bitReader->gotoBookmark(superBucketsInfoStart);
	superBucketsNumber=textLength/superBucketSize;
	if (textLength%superBucketSize>0)
		superBucketsNumber++;


	//retrieve super buckets starting positions
	superBucketsStarts=new uint64_t[superBucketsNumber];
	bitReader->gotoBookmark(superBucketsPositionsStart);
	for (uint64_t i=0;i<superBucketsNumber;i++)
		superBucketsStarts[i]=bitReader->getInt64();


	superBuckets=new SuperBucket*[superBucketsNumber];
	superBucketLoadingInProgress=new bool[superBucketsNumber];
	for (uint64_t i=0;i<superBucketsNumber;i++){
		superBuckets[i]=NULL;
		superBucketLoadingInProgress[i]=false;
	}

	//retrieve buckets starting positions
	bucketsStarts=new uint64_t[bucketsNumber];
	bitReader->gotoBookmark(bucketsPositionsStart);
	for (uint64_t i=0;i<bucketsNumber;i++)
		bucketsStarts[i]=bitReader->getInt64();
	//allocates buckets vector and sets number of loaded buckets to 0
	buckets=new Bucket*[bucketsNumber];
	bucketHeaderLoadingInProgress=new bool[bucketsNumber];
	bucketTextLoadingInProgress=new int64_t[bucketsNumber];
	for (uint64_t i=0;i<bucketsNumber;i++){
		buckets[i]=NULL;
		bucketHeaderLoadingInProgress[i]=false;
		bucketTextLoadingInProgress[i]=-1;
	}

	numberOfLoadedBuckets=0;
	numberOfLoadedSuperBuckets=0;

}


void EFMIndex::reconstructMarkedRows(){
	ofstream ofs("/tmp/markedRows.txt");

	markedRows=new map<uint64_t,uint64_t>();
	for (uint64_t bn=0;bn<bucketsNumber;bn++){
		for (auto& kv : buckets[bn]->markedRowsCache){
			uint64_t row=kv.first +bn*bucketSize;
			uint64_t textPosition= kv.second;
			markedRows->insert(UIntIntPair(row,textPosition));
			ofs << row << "," << textPosition*markingRate << endl;
		}
	}

	ofs.close();
}

void EFMIndex::readText(){
	for (uint64_t currentBucket=0;currentBucket<bucketsNumber;currentBucket++){
		checkBucketLoaded(bitReader,currentBucket);
		buckets[currentBucket]->readText(bitReader,bucketEncryptionRandomGenerator);
		//cout << currentBucket << "->" <<buckets[currentBucket]->actualPhysicalOffset << endl;
		buckets[currentBucket]->readMarkedRows(bitReader);
	}
	delete bucketEncryptionRandomGenerator;
}


inline void EFMIndex::checkSuperBucketLoadedOld(BitReader *bitReader,uint64_t currentSuperBucket){
	std::lock_guard<std::mutex> guard(superBucketLoadingMutex);
	if (superBuckets[currentSuperBucket]==NULL){
		SuperBucket *sb=new SuperBucket(this,currentSuperBucket,compactAlphabetSize);
		superBuckets[currentSuperBucket]=sb;
		sb->load(bitReader);
		numberOfLoadedSuperBuckets++;
	}
}


inline void EFMIndex::checkBucketLoadedOld(BitReader *bitReader,uint64_t currentBucket){
	std::lock_guard<std::mutex> guard(bucketLoadingMutex);
	//setBucketLoadingInProgress(currentBucket, true);
	checkSuperBucketLoaded(bitReader,getSuperbucketNumber(currentBucket));
	if (buckets[currentBucket]==NULL)
		allocateBucket(currentBucket,NULL);
	Bucket *b=buckets[currentBucket];
	if (!b->loaded){
		b->readHeader(bitReader);
		numberOfLoadedBuckets++;
	}
	//setBucketLoadingInProgress(currentBucket, false);
}

inline void EFMIndex::checkSuperBucketLoaded(BitReader *bitReader,uint64_t sn){
	bool actAsLoader=false;
	bool loadingInProgress=false;
	superBucketLoadingMutex.lock();

	bool loaded=superBuckets[sn]!=NULL && superBuckets[sn]->loaded;
	if (!loaded){
		loadingInProgress=superBucketLoadingInProgress[sn];
		if (!loadingInProgress){
			actAsLoader=true;
			superBucketLoadingInProgress[sn]=true;
		}
	}
	superBucketLoadingMutex.unlock();

	if (!loaded){
		if (actAsLoader){
			if (superBuckets[sn]==NULL){
				SuperBucket *sb=new SuperBucket(this,sn,compactAlphabetSize);
				superBuckets[sn]=sb;
				sb->load(bitReader);
				numberOfLoadedSuperBuckets++;

				superBucketLoadingMutex.lock();
				superBucketLoadingInProgress[sn]=false;
				loaded=true;
				superBucketLoadingMutex.unlock();
			}
		} else{
			while (loadingInProgress){
				superBucketLoadingMutex.lock();
				loadingInProgress=superBucketLoadingInProgress[sn];
				superBucketLoadingMutex.unlock();
			}
		}
	}

}

inline void EFMIndex::checkBucketLoaded(BitReader *bitReader,uint64_t bn){
	checkSuperBucketLoaded(bitReader,getSuperbucketNumber(bn));

	//coutMutex.lock();
	//std::cout << pthread_self() << " check bucket loaded " <<  bn << " " << endl;
	//coutMutex.unlock();

	bool actAsLoader=false;
	bool loadingInProgress=false;
	bucketLoadingMutex.lock();
	bool loaded=buckets[bn]!=NULL && buckets[bn]->loaded;
	if (!loaded){
		loadingInProgress=bucketHeaderLoadingInProgress[bn];
		if (!loadingInProgress){
			actAsLoader=true;
			bucketHeaderLoadingInProgress[bn]=true;
		}
	}
	bucketLoadingMutex.unlock();

	if (!loaded){
		if (actAsLoader){
			//coutMutex.lock();
			//std::cout << pthread_self() << " loading bucket header " <<  bn << " " << endl;
			//coutMutex.unlock();

			if (buckets[bn]==NULL){
				allocateBucket(bn,NULL);
			Bucket *b=buckets[bn];
			b->readHeader(bitReader);
			numberOfLoadedBuckets++;


			bucketLoadingMutex.lock();
			bucketHeaderLoadingInProgress[bn]=false;
			loaded=true;
			bucketLoadingMutex.unlock();

			//coutMutex.lock();
			//std::cout << pthread_self() << " finished bucket header loading " <<  bn << " " << endl;
			//coutMutex.unlock();

			}
		} else{
			//coutMutex.lock();
			//std::cout << pthread_self() << " waiting for bucket header loading " <<  bn << " " << endl;
			//coutMutex.unlock();

			while (loadingInProgress){
				bucketLoadingMutex.lock();
				loadingInProgress=bucketHeaderLoadingInProgress[bn];
				bucketLoadingMutex.unlock();
			}
		}
	}


}




void EFMIndex::writeBuckets(){
		for (uint64_t currentBucket=0;currentBucket<bucketsNumber;currentBucket++)
			buckets[currentBucket]->write(bitWriter);
}


void EFMIndex::createSpaceForSuspendedInfos(){
	//Free space for character occurrences list start byte within this header
	bitWriter->writeInt64(0);
	//Free space for marked rows buckets table start byte
	bitWriter->writeInt64(0);
	//Free space for super buckets positions list start byte
	bitWriter->writeInt64(0);
	//Free space for buckets positions list start byte
	bitWriter->writeInt64(0);
}

void EFMIndex::writeSuspendedInformations() {
	bitWriter->gotoBookmark(superBucketsPositionsStart);
	for (uint64_t i=0;i<superBucketsNumber;i++)
		bitWriter->writeInt64(superBucketsStarts[i]);

	bitWriter->gotoBookmark(bucketsPositionsStart);
	for (uint64_t i=0;i<bucketsNumber;i++)
		bitWriter->writeInt64(bucketsStarts[i]);

	bitWriter->gotoBookmark(startSuspendedInfos);
	bitWriter->writeInt64(characterOccurrencesStart->getPosition());
	bitWriter->writeInt64(markedRowsBucketsStart->getPosition());
	bitWriter->writeInt64(superBucketsPositionsStart->getPosition());
	bitWriter->writeInt64(bucketsPositionsStart->getPosition());
}

void EFMIndex::readSuspendedInformations(){
	//suspended informations
	//Character occurrences list offset within this header
	characterOccurrencesStart=new Bookmark(bitReader->getInt64(),0);
	//Marked rows buckets table start byte
	markedRowsBucketsStart=new Bookmark(bitReader->getInt64(),0);
	//Super buckets positions list start byte
	superBucketsPositionsStart=new Bookmark(bitReader->getInt64(),0);
	//Buckets positions list start byte
	bucketsPositionsStart=new Bookmark(bitReader->getInt64(),0);
}




void EFMIndex::load(string fileName){
	this->fileName=fileName;

	if (encryptionKey!=NULL)
		encryptionManager=new EncryptionManager(encryptionKey);
	else
		encryptionManager=new EncryptionManager(encryptionKeyFilePath);

	//create all objects need for multi-threading pattern match and start matcher threads


	bucketEncryptionRandomGenerator=new RandomGenerator(131072,encryptionManager->polySeed);

	if (loadWholeFileInMemory){
		wholeFileBuffer=MemoryBitReader::loadDataFromFile(fileName,wholeFileSize);
		bitReader=new MemoryBitReader(wholeFileBuffer,wholeFileSize);

	} else{
		bitReader=new FileBitReader(fileName);
		bitReader->setBlockSize(1048576);
	}
	bitReader->open();
	readHeader();

	uint32_t nt;
	if (numberOfThreads>superAlphabetOrder)
	  nt=superAlphabetOrder;
	else
	  nt=numberOfThreads;

	cout << "Number of pattern matching threads " << nt << endl;

	if (nt>1){
		randomGenerators = new PRG[nt];
		bitReaders=new PBR[nt];
		matchers=new PIM[nt];

		for (uint32_t i=0;i<nt;i++){
			if (i>0){
				randomGenerators[i]=new RandomGenerator(131072,encryptionManager->polySeed);
				if (loadWholeFileInMemory){
					bitReaders[i]=new MemoryBitReader(wholeFileBuffer,wholeFileSize);

				} else{
					bitReaders[i]=new FileBitReader(fileName);
					bitReaders[i]->setBlockSize(1048576);
				}
				bitReaders[i]->open();
			} else{
				bitReaders[i]=bitReader;
				randomGenerators[i]=bucketEncryptionRandomGenerator;
			}

			IntervalMatcher* matcher=new IntervalMatcher(*this,i,bitReaders[i],randomGenerators[i]);
			matchers[i]=matcher;
		}

		for (uint32_t i=0;i<nt;i++)
			matcherThreads.push_back(new std::thread(&IntervalMatcher::run, matchers[i]));
	}


}



void EFMIndex::close(){
	uint32_t nt;
	if (numberOfThreads>superAlphabetOrder)
	  nt=superAlphabetOrder;
	else
	  nt=numberOfThreads;



	// wait until all worker threads terminate
	if (nt>1){
		//signals the worker threads to terminate (index is closing)
		{
			std::lock_guard<std::mutex> lk(matcherThreadsWakeUpMutex);
			for (uint32_t i = 0; i < nt; i++)
				matchers[i]->doTerminate=true;
		}
		matcherThreadsWakeUpCV.notify_all();


		for (uint32_t i = 0; i < nt; i++)
		   matcherThreads[i]->join();


		for (uint32_t i=0;i<nt;i++){
			if (randomGenerators[i]!=NULL){
				delete randomGenerators[i];
				randomGenerators[i]=NULL;
			}
			if (bitReaders[i]!=NULL){
				bitReaders[i]->close();
				delete bitReaders[i];
				bitReaders[i]=NULL;
			}
			if (matcherThreads[i]!=NULL){
				delete matcherThreads[i];
				matcherThreads[i]=NULL;
			}
			if (matchers[i]!=NULL){
				delete matchers[i];
				matchers[i]=NULL;
			}
		}
	} else{
		bitReader->close();
		delete bitReader;

		delete bucketEncryptionRandomGenerator;
	}
	bitReader=NULL;
	bucketEncryptionRandomGenerator=NULL;
	if (bitWriter!=NULL)
		bitWriter->close();

	if (encryptionManager!=NULL){
		delete encryptionManager;
		encryptionManager=NULL;
	}

	if (bitReaders!=NULL)
		delete[] bitReaders;
	if (matchers!=NULL)
		delete[] matchers;
}


CryptoParams *EFMIndex::getCryptoParams(uint64_t bucketNumber,uint64_t bucketLength,uint32_t bucketAlphaSize){
	//std::lock_guard<std::mutex> guard(getCryptoParamsMutex);
	uint32_t *bucketKeyStream;
	bucketKeyStream=encryptionManager->computeBucketKeyStream(bucketNumber,0, bucketLength-1, bucketAlphaSize);
	CryptoParams *cryptoParams=new CryptoParams(bucketKeyStream);
	return cryptoParams;
}

/*
vector<Match*> *EFMIndex::exactMatchOld(string pattern){
	vector<SuperPattern*> superPatterns=computeSuperPatterns(pattern);
	vector<Match*> *matches=new vector<Match*>();
	for (SuperPattern *superPattern:superPatterns){
		   vector<Match*> *superPatternMatches;
		   if (markedRowsPercentage>0)
		      superPatternMatches=optimizedExactMatch(superPattern);
		   else
			  superPatternMatches=exactMatch(superPattern);
		   if (superPatternMatches->size()>0)
			   matches->insert(matches->end(), superPatternMatches->begin(), superPatternMatches->end());
	}
	return matches;
}*/



void EFMIndex::exactMatch(string pattern){

	for (uint64_t i=0;i<allSuperPatternMatches.size();i++)
					delete allSuperPatternMatches[i];
	allSuperPatternMatches.clear();

	//Controlla che tutti i caratteri del pattern appartengano dell'alfabeto primigenio
	//Se c'è anche un solo carattere al di fuori, allora il pattern non ha occorrenze
	bool validPattern=true;
	uint64_t i=0;
	while (validPattern && i<pattern.length()){
		uint8_t s=originalSymbols.size();
		char c=pattern[i];
		bool found=false;
		uint8_t j=0;
		while (j<s &&!found){
			if (c==originalSymbols[j])
				found=true;
			else
			  j++;
		}
		validPattern=found;
		i++;
	}

	if (validPattern){
		vector<SuperPattern*> superPatterns=computeSuperPatterns(pattern);
		//Determina tutti gli intervalli dell'alfabeto compatto che corrisponderanno
		//a diverse esecuzioni in parallelo (multi-threaded) dell'algoritmo di backward search
		for (uint64_t i=0;i<allCompactAlphabetIntervals.size();i++)
			delete allCompactAlphabetIntervals[i];
		allCompactAlphabetIntervals.clear();
		nextCompactAlphabetIntervalToProcess=0;

		for (SuperPattern *superPattern:superPatterns){
			uint32_t p=superPattern->length-1;
			SuperPatternCharacter *lastCharacter=superPattern->text[p];
			CompactAlphabetInterval *initialInterval;
			if (lastCharacter->type==SuperPatternCharacterType::fixedCharacter){
				initialInterval=new CompactAlphabetInterval(superPattern,lastCharacter->fixedValue,lastCharacter->fixedValue);
				allCompactAlphabetIntervals.push_back(initialInterval);
			}
			else{
				SuperPatternCharacter *beforeLast;
				if (superPattern->length > 1 &&
						(beforeLast=superPattern->text[p-1])->type==SuperPatternCharacterType::fixedCharacter){
					//the hat super-pattern is not build explicitly, but only setting a isHat flag on the pattern
					superPattern->isHat=true;
					initialInterval=new CompactAlphabetInterval(superPattern,beforeLast->fixedValue,beforeLast->fixedValue);
					allCompactAlphabetIntervals.push_back(initialInterval);
				} else //it happens only for very small patterns (they could have made of only two variable super-characters)
					for (CompactAlphabetInterval *interval: lastCharacter->compactAlphabetIntervals)
						allCompactAlphabetIntervals.push_back(interval);
			}
		}

		//cout << " " << allCompactAlphabetIntervals.size() << " ";

		//execute pattern matching with multiple threads
		uint32_t nt;
		if (numberOfThreads>superAlphabetOrder)
		  nt=superAlphabetOrder;
		else
		  nt=numberOfThreads;

		if (nt>1){
			uint32_t workFinishedThreads=0;
			//terminatedPatternMatching.resize(nt);
		    //notifiedCompletion.resize(nt);
			//signals the worker threads to start pattern matching
			{
				std::lock_guard<std::mutex> lk(matcherThreadsWakeUpMutex);
				for (uint32_t i=0;i<nt;i++){
					//terminatedPatternMatching[i]=false;
					//notifiedCompletion[i]=false;
					matchers[i]->doPatternMatching=true;
					matchers[i]->doTerminate=false;
					matchers[i]->finishedPatternMatching=false;
				}
				matcherThreadsWakeUpCV.notify_all();
			}

			//wait until all worker threads finish their work
			//do
			//{
			{
				std::unique_lock<std::mutex> lk(mainThreadWakeUpMutex);
				mainThreadWakeUpCV.wait(lk, [&]{return allMatchingThreadsFinished();});
			}
			//	workFinishedThreads++;
			//} while (workFinishedThreads < nt);
		}
		else{//SINGLE-THREAD
			for (CompactAlphabetInterval *interval:allCompactAlphabetIntervals){
				exactMatch(interval,bitReader,bucketEncryptionRandomGenerator);
			}
		}

	}

}






inline vector<Match*> *EFMIndex::exactMatch(SuperPattern *superPattern,Interval *interval,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator){
		int64_t p,sp,ep,i;

		uint32_t l;   //super-pattern length
		if (superPattern->isHat)
			l=superPattern->length-1;
		else
			l=superPattern->length;
		p=l-1;
		SuperPatternCharacter *sc=superPattern->text[p];
		vector<Match*> *intervalMatches=new vector<Match*>();
		sp=precedingCharactersOccurrences[interval->sp];
		ep=precedingCharactersOccurrences[interval->ep+1]-1;
		i=p-1;
		while (sp<=ep && i>=0){
			sc=superPattern->text[i];
			if (sc->type==SuperPatternCharacterType::fixedCharacter){  //vale sicuramente per tutti i caratteri "interni" al super-pattern
				uint32_t c=sc->fixedValue;
				Interval *ei=backwardStep(c, new Interval(sp,ep),bitReader,bucketEncryptionRandomGenerator);
				sp=ei->sp;
				ep=ei->ep;
			}
			else{  //non può che essere il primo carattere del super-pattern (siamo, dunque, al passo finale)
				Interval *si=new Interval(sp,ep);
				for (int64_t saIndex=si->sp;saIndex<=si->ep;saIndex++){
					Interval *ei=backwardStepIfMatch(saIndex, sc->text,bitReader,bucketEncryptionRandomGenerator);
					if (ei!=NULL){
						Match *match=new Match();
						match->superPattern=superPattern;
						match->suffixArrayInterval=new Interval(ei->sp,ei->ep);
						intervalMatches->push_back(match);
					}

				}

			}
			i--;
		}
		/**
		 * Nel caso in cui il primo super-carattere del pattern (l'ultimo analizzato durante
		 * la backward search) sia a lunghezza variabile, ma non sia l'unico carattere del pattern,
		 * gli eventuali match sono già stati aggiunti ed esso va aggiunto solo nel caso in cui l'ultimo
		 * super-carattere sia a lunghezza fissa
		 */
		if ((sc->type==SuperPatternCharacterType::fixedCharacter || l==1) &&
			     ep >= sp){
				Match *match=new Match();
				match->superPattern=superPattern;
				match->suffixArrayInterval=new Interval(sp,ep);
				intervalMatches->push_back(match);
			}

		return intervalMatches;
	}


//The found match (given in input) is the suffix array range corresponding to the occurrences of the
//sub-pattern hatSuperPattern, composed of the first m-1 characters P_0,P_1,..., P_(m-2)
//of the real superPattern.
//Among them we must find the ones for which the following character matches with the last character
//P_(m-1) of the original pattern.
inline void EFMIndex::findWholePatternMatches(Match *match,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator){
		SuperPattern *superPattern=match->superPattern;
		SuperPatternCharacter *lastCharacter=superPattern->text[superPattern->length-1];
		Interval *sai=match->suffixArrayInterval;
		for (int64_t i=sai->sp;i<=sai->ep;i++){
			//find the suffix array position of this suffix
			int64_t superTextPosition=getRowPosition(i,bitReader,bucketEncryptionRandomGenerator);
			//extract the character in position i+m-1 (where m is the whole super pattern's length)
			uint32_t m=match->superPattern->length;
			string foundCharacter=extractSuperCharacter(superTextPosition+m-1,bitReader,bucketEncryptionRandomGenerator);
			string lastCharacterAsString=lastCharacter->text;
			bool matching=true;
			uint32_t j=0;
			while (j<lastCharacterAsString.size() && matching){
				if (lastCharacterAsString[j]!='?' && lastCharacterAsString[j]!=foundCharacter[j])
					matching=false;
				j++;
			}
			if (matching){
				Match *wholePatternMatch=new Match();
				wholePatternMatch->suffixArrayInterval=new CompactAlphabetInterval(match->superPattern,i, i);
				wholePatternMatch->superPattern=match->superPattern;
				wholePatternMatch->occurrences=new vector<uint64_t>();

				//computes the original text position of the occurrence
				int64_t originalTextPosition=superTextPosition * superAlphabetOrder +
						match->superPattern->shiftAmount;
				//Adds to the found position the number of initial N characters which
				//are not included in sBWT
				originalTextPosition+=initialNCharacters;
				wholePatternMatch->occurrences->push_back(originalTextPosition);

				allSuperPatterMatchesMutex.lock();
				allSuperPatternMatches.push_back(wholePatternMatch);
				allSuperPatterMatchesMutex.unlock();
			}
	}
}

void EFMIndex::exactMatch(CompactAlphabetInterval *interval,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator){
		int64_t p,sp,ep,i;
		SuperPattern *superPattern=interval->superPattern;
		uint32_t l;   //super-pattern length
		if (superPattern->isHat)
			l=superPattern->length-1;
		else
			l=superPattern->length;
		p=l-1;
		SuperPatternCharacter *sc=superPattern->text[p];
		sp=precedingCharactersOccurrences[interval->sp];
		ep=precedingCharactersOccurrences[interval->ep+1]-1;
		i=p-1;
		while (sp<=ep && i>=0){
			sc=superPattern->text[i];
			if (sc->type==SuperPatternCharacterType::fixedCharacter){  //vale sicuramente per tutti i caratteri "interni" al pattern
				uint32_t c=sc->fixedValue;
				Interval *ei=backwardStep(c, new Interval(sp,ep),bitReader,bucketEncryptionRandomGenerator);
				sp=ei->sp;
				ep=ei->ep;
			}
			else{  //non può che essere il carattere del super-pattern (siamo, dunque, al passo finale)
				Interval *si=new Interval(sp,ep);
				for (int64_t saIndex=si->sp;saIndex<=si->ep;saIndex++){
					Interval *ei=backwardStepIfMatch(saIndex, sc->text,bitReader,bucketEncryptionRandomGenerator);
					if (ei!=NULL){
						Match *match=new Match();
						match->superPattern=superPattern;
						match->suffixArrayInterval=new Interval(ei->sp,ei->ep);
						if (superPattern->isHat)
							findWholePatternMatches(match,bitReader,bucketEncryptionRandomGenerator);
						else{
							allSuperPatterMatchesMutex.lock();
							allSuperPatternMatches.push_back(match);
							allSuperPatterMatchesMutex.unlock();
						}
					}
				}

			}
			i--;
		}
		/**
		 * Nel caso in cui il primo super-carattere del pattern (l'ultimo analizzato durante
		 * la backward search) sia a lunghezza variabile, ma non sia l'unico carattere del pattern,
		 * gli eventuali match sono già stati aggiunti ed esso va aggiunto solo nel caso in cui l'ultimo
		 * super-carattere sia a lunghezza fissa
		 */
		if ((sc->type==SuperPatternCharacterType::fixedCharacter || l==1) &&
			     ep >= sp){
				Match *match=new Match();
				match->superPattern=superPattern;
				match->suffixArrayInterval=new Interval(sp,ep);
				if (superPattern->isHat)
					findWholePatternMatches(match,bitReader,bucketEncryptionRandomGenerator);
				else{
					allSuperPatterMatchesMutex.lock();
					allSuperPatternMatches.push_back(match);
					allSuperPatterMatchesMutex.unlock();
				}
		}

	}


Interval *EFMIndex::backwardStep(uint32_t c,Interval *si,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator){
	//std::lock_guard<std::mutex> guard(readingProtectionMutex);

	Interval *ei=new Interval();
	uint64_t sn,bn,bo;
	getBwtPositionCoordinates(si->sp-1,&sn,&bn,&bo);
	//computes the occurrences of the compact alphabet character c in this bucket
	//and in previous buckets of this superbucket
	uint64_t thisSuperBucketOccurrences=0;
	checkSuperBucketLoaded(bitReader,sn);
	int64_t superBucketCode=superBuckets[sn]->getRemappedCode(c);
	if (superBucketCode!=-1){
		checkBucketLoaded(bitReader,bn);
		Bucket *b=buckets[bn];
		thisSuperBucketOccurrences=b->getOffsetInformations(bitReader,bucketEncryptionRandomGenerator, superBucketCode, bo,
				   OffsetInformationsType::CountOnly)->occurrencesCount;
	}
	ei->sp=precedingCharactersOccurrences[c] +  superBuckets[sn]->previousSuperbucketsOccurrences[c] +
			thisSuperBucketOccurrences;

	getBwtPositionCoordinates(si->ep,&sn,&bn,&bo);

	superBucketCode=-1;
	thisSuperBucketOccurrences=0;
	checkSuperBucketLoaded(bitReader,sn);
	superBucketCode=superBuckets[sn]->getRemappedCode(c);
	if (superBucketCode!=-1){
		checkBucketLoaded(bitReader,bn);
		Bucket *b=buckets[bn];
		thisSuperBucketOccurrences=b->getOffsetInformations(bitReader,bucketEncryptionRandomGenerator, superBucketCode, bo,
				   OffsetInformationsType::CountOnly)->occurrencesCount;
	}
	ei->ep=precedingCharactersOccurrences[c] + superBuckets[sn]->previousSuperbucketsOccurrences[c] +
			thisSuperBucketOccurrences-1;
	return ei;
}


/**
 * Verifies if the BWT character in position saIndex matches with a pattern and, if so, does a backward step
 * @param pattern of the variable character to find (with ? initial wildcards)
 * @param saIndex	index in suffix array
 * @return
 * @throws Exception
 */
Interval *EFMIndex::backwardStepIfMatch(uint64_t saIndex,string pattern,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator){
	//std::lock_guard<std::mutex> guard(readingProtectionMutex);

	Interval *result=NULL;
	uint64_t sn,bn,bo;

	getBwtPositionCoordinates(saIndex,&sn,&bn,&bo);
	checkBucketLoaded(bitReader,bn);
	Bucket *b=buckets[bn];
	OffsetInformations *offsetInformations=b->getOffsetInformations(bitReader,bucketEncryptionRandomGenerator, -1, bo,
			   OffsetInformationsType::RetrieveCharacter);
	uint32_t superBucketCode=offsetInformations->superBucketCode;
	//int compactAlphabetCode=superBuckets[sn].decode_map.get(superBucketCode);
	uint32_t compactAlphabetCode=superBuckets[sn]->decode_map[superBucketCode];
	uint32_t scrambledAlphabetCode=decodeMap[compactAlphabetCode];
	string character=superAlphabet->getSymbol(scrambledAlphabetCode);
	bool matching=true;
	uint32_t i=0;
	while (i<pattern.length() && matching){
		if (pattern[i]!='?' && pattern[i]!=character[i])
			matching=false;
		i++;
	}
	if (matching){
		//compute LF-mapping
		uint64_t first=precedingCharactersOccurrences[compactAlphabetCode] +  superBuckets[sn]->previousSuperbucketsOccurrences[compactAlphabetCode] +
			offsetInformations->occurrencesCount -1;  //-1 because array is 0-based
		result=new Interval(first,first);
	}
	return result;
}



/**
	 * Builds all the super-patterns needed to query the index, build in terms
	 * of super-alphabet symbols, for an original alphabet pattern.
	 * Super-patterns will be constructed as follows:
	 *   1) compute superPatternLength as ceil(pattern.length/superAlphabetOrder)
	 *   2) compute a super-pattern for each possible right shift of the
	 *      pattern on  the superPatternLength positions
	 * @param pattern
	 * @return
	 */
vector<SuperPattern*> EFMIndex::computeSuperPatterns(string pattern) {
	vector<SuperPattern*> superPatterns;
	uint32_t m=pattern.length();
	uint8_t k=superAlphabetOrder;
	uint32_t no=superAlphabet->originalAlphabetSize;

	 //for each shift amount, compute all super-patterns corresponding to it
	for (uint8_t sa=0;sa<superAlphabetOrder;sa++){
		//compute super-pattern length as (m+sa)/k,
		//because m+sa represents the position of the last character of the original alphabet
		//within the shifted pattern
		uint32_t superPatternLength=ceil((double)(m+sa)/(double)k);
		SuperPattern *sp=new SuperPattern();
		sp->shiftAmount=sa;
		sp->text=new SuperPatternCharacter*[superPatternLength];
		sp->length=superPatternLength;
		uint32_t c=0;
		bool foundSuperCharacterOutsideOfCompactAlphabet=false;
		while (c<superPatternLength && !foundSuperCharacterOutsideOfCompactAlphabet){
			SuperPatternCharacter *sc=new SuperPatternCharacter();

			//compute the initial position of this super-character within the original pattern
			int32_t superCharacterInitialPositionInOriginalPattern=c*k-sa;
			if (superCharacterInitialPositionInOriginalPattern<0){
				sc->prefixLength=-superCharacterInitialPositionInOriginalPattern;
				sc->fixedPartStartInPattern=0;
			} else{
				sc->prefixLength=0;
				sc->fixedPartStartInPattern=superCharacterInitialPositionInOriginalPattern;
			}

			//compute the ending position of this super-character within the original pattern
			int32_t superCharacterEndingPositionInOriginalPattern=(c+1)*k-sa-1;
			if (superCharacterEndingPositionInOriginalPattern > m-1){
				sc->suffixLength=superCharacterEndingPositionInOriginalPattern-m+1;
				sc->fixedPartEndInPattern=m-1;
			} else{
				sc->suffixLength=0;
				sc->fixedPartEndInPattern=superCharacterEndingPositionInOriginalPattern;
			}

			//if needed, compute all possible prefix values
			if (sc->prefixLength>0){
				uint32_t npv=pow(no, sc->prefixLength);
				sc->prefixValuesSize=npv;
				sc->prefixValues=new uint[npv];
				for (uint32_t i=0;i<npv;i++){
					sc->prefixValues[i]=i*(uint32_t)pow((double) no, (double)(k - sc->prefixLength));
				}
			}

			//if needed, compute suffix range end (suffix range start is always 0)
			if (sc->suffixLength>0)
				sc->suffixRangeEnd = (uint)pow((double)no, (double)sc->suffixLength)-1;

			//compute fixed part value
			uint32_t exponent=sc->suffixLength;
			int64_t i=sc->fixedPartEndInPattern;
			sc->fixedValue=0;
			while (i>=sc->fixedPartStartInPattern){
				sc->fixedValue+=superAlphabet->getOriginalAlphabetCode(pattern[i])* (uint32_t)pow((double)no,(double) exponent);
				i--;
				exponent++;
			}

			//for an easier debug
			if (sc->prefixLength>0 || sc->suffixLength>0){
				sc->type=SuperPatternCharacterType::variableCharacter;

				sc->text="";
				for (i=0;i<sc->prefixLength;i++)
					sc->text+="?";
				for (i=sc->fixedPartStartInPattern;i<=sc->fixedPartEndInPattern;i++)
					sc->text+=pattern[i];
				for (i=0;i<sc->suffixLength;i++)
					sc->text+="?";

				if (c>1){ //this super-character is the last one and there is at least an internal super-character
					//compute all compact alphabet codes sub-intervals corresponding to this variable character
					vector<uint32_t> values;

					if (sc->suffixLength>0 && sc->prefixLength>0){
						//for each prefix value
						//  for each suffix value
						//    computes a single character code summing prefixValue,suffixValue and fixedValue
						for (i=0;i<sc->prefixValuesSize;i++){
							for (uint32_t j=0;j<=sc->suffixRangeEnd;j++){
								uint32_t entireAlphabetCode=superAlphabet->getCode(sc->prefixValues[i] + sc->fixedValue +j);
								if ((*charactersBitSet)[entireAlphabetCode])
									values.push_back(charactersRemappings[entireAlphabetCode]);
							}
						}
					}
					else if (sc->suffixLength>0){
						for (uint32_t j=0;j<=sc->suffixRangeEnd;j++){
								uint32_t entireAlphabetCode=superAlphabet->getCode(sc->fixedValue +j);
								if ((*charactersBitSet)[entireAlphabetCode])
									values.push_back(charactersRemappings[entireAlphabetCode]);
						}
					}else{  //only prefix length is > 0
						for (i=0;i<sc->prefixValuesSize;i++){
							uint32_t entireAlphabetCode=superAlphabet->getCode(sc->prefixValues[i] + sc->fixedValue);
							if ((*charactersBitSet)[entireAlphabetCode])
								values.push_back(charactersRemappings[entireAlphabetCode]);
						}
					}


					if (values.size()>0){
						sort(values.begin(), values.end());
						CompactAlphabetInterval *currentInterval=new CompactAlphabetInterval();
						currentInterval->superPattern=sp;
						currentInterval->sp=values[0];
						uint32_t h;
						for (h=1;h<values.size();h++){
							if (values[h]>values[h-1]+1){
								currentInterval->ep=values[h-1];
								sc->compactAlphabetIntervals.push_back(currentInterval);
								currentInterval=new CompactAlphabetInterval();
								currentInterval->superPattern=sp;
								currentInterval->sp=values[h];
							}
						}
						currentInterval->ep=values[h-1];
						sc->compactAlphabetIntervals.push_back(currentInterval);
					}
				}
			}
			else{
				sc->type=SuperPatternCharacterType::fixedCharacter;
				//transform character natural code into the scrambled super-alphabet code and remaps it vs. the compact one
				uint32_t entireAlphabetCode=superAlphabet->getCode(sc->fixedValue);
				uint32_t size=(*charactersBitSet).size();
				if ((*charactersBitSet)[entireAlphabetCode])
					sc->fixedValue=charactersRemappings[superAlphabet->getCode(sc->fixedValue)];
				else
					sc->fixedValue=-1;
			}

			if (sc->type==SuperPatternCharacterType::fixedCharacter && sc->fixedValue == -1
				//||sc->type==SuperPatternCharacterType::variableCharacter && sc->compactAlphabetIntervals.size() == 0
				)
				  foundSuperCharacterOutsideOfCompactAlphabet=true;
			else{
				sp->text[c]=sc;
				c++;
			}
		}
		if (!foundSuperCharacterOutsideOfCompactAlphabet)
			superPatterns.push_back(sp);
	}
	return superPatterns;
}


/**
	 * Retrieve the super-character that occupies a certain position in the super-text
	 * @param position
	 * @return
	 * @throws Exception
	 */
string EFMIndex::extractSuperCharacter(uint64_t position,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator){
	//std::lock_guard<std::mutex> guard(readingProtectionMutex);
	uint32_t superCharacter;
	//computes the start position of the backward search and the number of characters to skip
	//the start position is the marked row
	uint64_t i;
	uint64_t skip;
	//computes the next marked row progressive, assuming that the text is disposed on a circular tape;
	//the virtual position can be above the end of the text
	uint64_t nextMarkedRowVirtualProgressive=(position/markingRate)+1;
	//the real progressive corresponding to a virtual progressive above the end of the text is 0,
	//corresponding to the first marked row
	uint64_t markedRowRealProgressive=nextMarkedRowVirtualProgressive%markedRowsNumber;

	uint64_t bn=markedRowsBuckets[markedRowRealProgressive];
	checkBucketLoaded(bitReader,bn);
	i = buckets[bn]->getMarkedBwtPosition(bitReader, markedRowRealProgressive);

	uint64_t nextMarkedRowAbsolutePosition=nextMarkedRowVirtualProgressive*markingRate;
	//se la prossima riga marcata cade dopo la fine del testo, allora essa è la riga 0, che,
	//se si considera il testo come disposto su un nastro circolare, cade nella posizione
	//textLength;
	if (nextMarkedRowAbsolutePosition>textLength-1)
		nextMarkedRowAbsolutePosition=textLength;
	skip=nextMarkedRowAbsolutePosition-position-1;

	uint64_t sn,bo;
	uint32_t superBucketCode;
	uint32_t compactAlphabetCode;
	//cout << "skip" <<endl;
	for (uint64_t s=0;s<skip;s++){
		getBwtPositionCoordinates(i,&sn,&bn,&bo);
		checkBucketLoaded(bitReader,bn);
		Bucket *b=buckets[bn];
		OffsetInformations *offsetInformations=b->getOffsetInformations(bitReader,bucketEncryptionRandomGenerator, -1, bo, OffsetInformationsType::RetrieveCharacter);
		superBucketCode=offsetInformations->superBucketCode;
		compactAlphabetCode=b->superBucket->decode_map[superBucketCode];
		//a backward step is necessary
		i=precedingCharactersOccurrences[compactAlphabetCode] +  superBuckets[sn]->previousSuperbucketsOccurrences[compactAlphabetCode] +
					offsetInformations->occurrencesCount -1;  //-1 because array is 0-based
		//cout << i <<endl;

	}
	//there are no more characters to skip: extract the super-character in the position i of the BWT
	getBwtPositionCoordinates(i,&sn,&bn,&bo);
	checkBucketLoaded(bitReader,bn);
	Bucket *b=buckets[bn];
	OffsetInformations *offsetInformations=b->getOffsetInformations(bitReader,bucketEncryptionRandomGenerator, -1, bo, OffsetInformationsType::RetrieveCharacter);
	superBucketCode=offsetInformations->superBucketCode;
	compactAlphabetCode=b->superBucket->decode_map[superBucketCode];
	superCharacter=decodeMap[compactAlphabetCode];
	return superAlphabet->getSymbol(superCharacter);

}



uint64_t EFMIndex::countOccurrences(vector<Match*> *matches){
	uint64_t result=0;
	for (Match *match:*matches){
		result+=match->suffixArrayInterval->getLength();
	}
	return result;
}

vector<uint64_t> *EFMIndex::locateOccurrences(vector<Match*> *matches){
		vector<uint64_t> *occs=new vector<uint64_t>();
		for (Match *match:*matches){
			if (match->occurrences!=NULL)
				occs->insert(occs->end(),match->occurrences->begin(),match->occurrences->end());
			else{
				vector<uint64_t> *matchOccurrences=locateOccurrences(match);
				occs->insert(occs->end(),matchOccurrences->begin(),matchOccurrences->end());
			}
		}
		if  (occs->size()>1)
			sort(occs->begin(),occs->end());

		//TODO: capire perchè a volte vengono ritornate due posizioni uguali
		//Il codice seguente è stato aggiunto per fixare la PCR prima della presentazione del
		//progetto di ricerca 2013
		for (int64_t i=occs->size()-1;i>0;i--)
			if ((*occs)[i]==(*occs)[i-1])
				occs->erase(occs->begin()+i);
		return occs;
}


vector<uint64_t> *EFMIndex::locateOccurrences(Match *match){
	if (markedRowsPercentage==0)
		throw runtime_error("The index doesn't contain marked rows: locate operation not supported");

	vector<uint64_t> *occs=new vector<uint64_t>();
	Interval *sai=match->suffixArrayInterval;
	for (uint64_t i=sai->sp;i<=sai->ep;i++){
		uint64_t superTextPosition=getRowPosition(i,bitReader,bucketEncryptionRandomGenerator);
		//Computes the original text position adding shift amount to
		//the position in the original text of the super-character position.
		uint64_t originalTextPosition=superTextPosition * superAlphabetOrder +
				match->superPattern->shiftAmount;
		//Adds to the found position the number of initial N characters which
		//are not included in sBWT
		originalTextPosition+=initialNCharacters;

		occs->push_back(originalTextPosition);
	}
	return occs;
}

string EFMIndex::extractSnippet(uint64_t from, uint64_t to) {
	if (markedRowsPercentage==0)
		throw runtime_error("The index doesn't contain marked rows: the extraction operation is not supported");

	//initial and final N characters are not represented in the index
	from= from - initialNCharacters;
	to= to - initialNCharacters;

	if (from > to)
		return "";

	//computes from in terms of super text and from offset
	uint64_t stFrom=from/superAlphabetOrder;
	uint64_t stFromOffset=from%superAlphabetOrder;
	//computes to in terms of super text and to offset
	uint64_t stTo=to/superAlphabetOrder;
	uint64_t stToOffset=to%superAlphabetOrder;

	if (stTo>textLength-2)
		stTo=textLength-2;
	if (stFrom > stTo || (stFrom==stTo && stFromOffset > stToOffset))
			return "";

	if (stFrom==0 && stTo==textLength-1)
		return extractAll();
	else {
		return extractSnippet(stFrom,stFromOffset,stTo,stToOffset);
	}
}


/**
	 * Extracts the original alphabet snippet between the stFromOffset of the stFrom^ super-character
	 * to the stToOffset of the stTo^ super character
	 * @param stFrom       starting super-alphabet character
	 * @param stFromOffset offset of original alphabet starting character whithin the lastCharacterPosition super-character
	 * @param stTo		   ending super-alphabet character
	 * @param stToOffset   offset of original alphabet ending character whithin the last super-character
	 * @return
	 * @throws Exception
	 */

string EFMIndex::extractSnippet(uint64_t stFrom,uint8_t stFromOffset,uint64_t stTo,uint8_t stToOffset){
	/*
	string  sgood=extractSuperCharacter(337381048,bitReader,bucketEncryptionRandomGenerator);
	sgood=extractSuperCharacter(337381099,bitReader,bucketEncryptionRandomGenerator);
	sgood=extractSuperCharacter(337381149,bitReader,bucketEncryptionRandomGenerator);
	string sbad=extractSuperCharacter(337381150,bitReader,bucketEncryptionRandomGenerator);*/

	//allocates an array which will contain all the snippet's super-characters
	uint64_t remainingCharacters=stTo-stFrom+1;
	uint64_t snippetCharsLength=remainingCharacters;
	uint32_t *snippetChars=new uint32_t[snippetCharsLength];
	//computes the start position of the backward search and the number of characters to skip
	//the start position is the marked row
	uint64_t i;
	uint64_t skip;
	//computes the next marked row progressive, assuming that the text is disposed on a circular tape;
	//the virtual position can be above the end of the text
	uint64_t nextMarkedRowVirtualProgressive=(stTo/markingRate)+1;
	//the real progressive corresponding to a virtual progressive above the end of the text is 0,
	//corresponding to the first marked row
	uint64_t markedRowRealProgressive=nextMarkedRowVirtualProgressive%markedRowsNumber;

	uint64_t bn=markedRowsBuckets[markedRowRealProgressive];
	checkBucketLoaded(bitReader,bn);
	i = buckets[bn]->getMarkedBwtPosition(bitReader, markedRowRealProgressive);

	uint64_t nextMarkedRowAbsolutePosition=nextMarkedRowVirtualProgressive*markingRate;
	//se la prossima riga marcata cade dopo la fine del testo, allora essa è la riga 0, che,
	//se si considera il testo come disposto su un nastro circolare, cade nella posizione
	//textLength;
	if (nextMarkedRowAbsolutePosition>textLength-1)
		nextMarkedRowAbsolutePosition=textLength;
	skip=nextMarkedRowAbsolutePosition-stTo-1;


	uint64_t sn;
	uint64_t bo;
	uint32_t superBucketCode;
	uint32_t compactAlphabetCode;
	for (uint64_t s=0;s<skip;s++){
		getBwtPositionCoordinates(i,&sn,&bn,&bo);

		checkBucketLoaded(bitReader,bn);
		Bucket *b=buckets[bn];

		OffsetInformations *offsetInformations=b->getOffsetInformations(bitReader,bucketEncryptionRandomGenerator, -1, bo, OffsetInformationsType::RetrieveCharacter);
		superBucketCode=offsetInformations->superBucketCode;
		compactAlphabetCode=b->superBucket->decode_map[superBucketCode];
		//a backward step is necessary
		i=precedingCharactersOccurrences[compactAlphabetCode] +  superBuckets[sn]->previousSuperbucketsOccurrences[compactAlphabetCode] +
					offsetInformations->occurrencesCount -1;  //-1 because array is 0-based

	}

	while (remainingCharacters>0){
		getBwtPositionCoordinates(i,&sn,&bn,&bo);

		checkBucketLoaded(bitReader,bn);
		Bucket *b=buckets[bn];

		OffsetInformations *offsetInformations=b->getOffsetInformations(bitReader,bucketEncryptionRandomGenerator, -1, bo, OffsetInformationsType::RetrieveCharacter);
		superBucketCode=offsetInformations->superBucketCode;
		compactAlphabetCode=b->superBucket->decode_map[superBucketCode];

		snippetChars[remainingCharacters-1]=decodeMap[compactAlphabetCode];
		remainingCharacters--;

		if (remainingCharacters>0)  //a backward step is necessary
			i=precedingCharactersOccurrences[compactAlphabetCode] +  superBuckets[sn]->previousSuperbucketsOccurrences[compactAlphabetCode]  +
					offsetInformations->occurrencesCount -1;  //-1 because array is 0-based

	}

	if (snippetCharsLength==1)
		return superAlphabet->getSymbol(snippetChars[0]).substr(stFromOffset,stToOffset+1);
	else{
		string result="";
		if (snippetCharsLength>0){
			result+=superAlphabet->getSymbol(snippetChars[0]).substr(stFromOffset);
			for (uint64_t c=1;c<snippetCharsLength-1;c++)
				result+=superAlphabet->getSymbol(snippetChars[c]) ;
			result+=superAlphabet->getSymbol(snippetChars[snippetCharsLength-1]).substr(0, stToOffset+1);
		}
		return result;
	}

}


string EFMIndex::extractAll(){
	throw "Not yet implemented";
}




/**
 * Returns the position in text corresponding to a certain suffix array row
 * @param i the suffix array row
 * @return
 */
uint64_t EFMIndex::getRowPosition(uint64_t i,BitReader *bitReader,RandomGenerator *bucketEncryptionRandomGenerator){
	//std::lock_guard<std::mutex> guard(readingProtectionMutex);
	int64_t t=0;
	int64_t markedTextPosition=-1;
	uint64_t sn,bn,bo;
	//while (i%markingRate!=0){
	while (markedTextPosition==-1){
		getBwtPositionCoordinates(i,&sn,&bn,&bo);
		checkBucketLoaded(bitReader,bn);
		Bucket *b=buckets[bn];
		//markedTextPosition=b.markedRows.get(bo);
		markedTextPosition=b->getMarkedTextPosition(bitReader, bo);

		if (markedTextPosition==-1){  //a backward step is necessary
			//computes LF-mapping for this row
			OffsetInformations *offsetInformations=b->getOffsetInformations(bitReader,bucketEncryptionRandomGenerator, -1, bo,
				   OffsetInformationsType::RetrieveCharacter);
			uint32_t superBucketCode=offsetInformations->superBucketCode;
			//int compactAlphabetCode=b.getSuperBucket().decode_map.get(superBucketCode);
			uint32_t compactAlphabetCode=b->superBucket->decode_map[superBucketCode];
			i=precedingCharactersOccurrences[compactAlphabetCode] +
					superBuckets[sn]->previousSuperbucketsOccurrences[compactAlphabetCode] +
					offsetInformations->occurrencesCount -1;  //-1 because array is 0-based
			t++;
		}
	}
	//System.out.println("marking rate="+ markingRate+ " t="+t);
	uint64_t pos=(markedTextPosition*markingRate +t)%textLength;
	return  pos;
}


void EFMIndex::initializeStatistics(uint64_t originalTextLength) {
	if (statistics!=NULL){
		delete statistics;
	}
	statistics=new Statistics();
	statistics->originalTextLength=originalTextLength;
	statistics->superbucketAlphabetAverageSize=0;
	statistics->superbucketBitmapsCumulativeSize=0;
	statistics->superbucketOccurrencesTablesCumulativeSize=0;
	statistics->bucketAlphabetAverageSize=0;
	statistics->bucketBitmapsCumulativeSize=0;
	statistics->bucketOccurrencesTablesCumulativeSize=0;
	statistics->bucketMarkedRowsTablesCumulativeSize=0;
}


void EFMIndex::finalizeStatistics() {
	statistics->compressionRatio=(double)(statistics->compressedSize)/(double)(statistics->originalTextLength + initialNCharacters + finalNCharacters);
	statistics->realCompressionRatio=(double)(statistics->compressedSize)/(double)(superAlphabet->order*statistics->textLength);
	statistics->superbucketAlphabetAverageSize/=superBucketsNumber;
	statistics->superbucketBitmapAverageSize=statistics->superbucketBitmapsCumulativeSize/superBucketsNumber;
	statistics->superbucketOccurrencesTableAverageSize=statistics->superbucketOccurrencesTablesCumulativeSize/superBucketsNumber;
	statistics->bucketAlphabetAverageSize/=bucketsNumber;
	statistics->bucketBitmapAverageSize=statistics->bucketBitmapsCumulativeSize/bucketsNumber;
	statistics->bucketOccurrencesTableAverageSize=statistics->bucketOccurrencesTablesCumulativeSize/bucketsNumber;
	statistics->bucketMarkedRowsTableAverageSize=statistics->bucketMarkedRowsTablesCumulativeSize/bucketsNumber;

}






} /* namespace std */
