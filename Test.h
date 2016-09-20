/*
 * Test.h
 *
 *  Created on: 31/dic/2014
 *      Author: fernando
 */

#ifndef TEST_H_
#define TEST_H_

#include <vector>
#include <string>
#include "ecrypt-portable.h"
#include "EFMCollection.h"

namespace std {

typedef struct PatternSearchInfo{
	double loadingTime;
	double loadingTimesStdDev;
	double countTime;
	double countTimesStdDev;
	double countMeanTime;
	double countMeanTimesStdDev;
	bool locateSupported;
	double locateTime;
	double locateTimesStdDev;
	double locateMeanTime;
	double locateMeanTimesStdDev;
	double searchTime;
	double searchTimesStdDev;
	double searchMeanTime;
	double searchMeanTimesStdDev;
	double numberOfLoadedBuckets;
	double numberOfBuckets;       //total number of buckets in the index
	double numberOfLoadedSuperBuckets;
	double numberOfSuperBuckets;  //total number of super-buckets in the index
	double numberOfOccurrences;
	double meanUsedMemory;
	double usedMemoryDevStd;
	double usedMemory;
	double medianUsedMemory;
	double maxUsedMemory;
	double minUsedMemory;
	vector<uint64_t> *occurrences;
} PatternSearchInfo;


typedef struct PatternSet{
	int length;
	vector<string> patterns;
} PatternSet;

typedef struct TestGroup{
	string inputFileType;
	string databaseRoot;
	string sequencesDirectory;  //for tests on collections
	string referenceIdentifier;//for tests on collections
	string inputFileName;  //for tests on single chromosomes
	string indexFileName; //for tests on collections
	string indexParentDirectory;
	bool buildIndex;
	bool patternSearch;
	bool naiveSearch;
	bool loadReferenceInMemory;
	bool loadIndexInMemory;
	int blockSize;
	int userId;
	string userPrivateKey;
	string userPassword;
	vector<PatternSet> patternSets;
	vector<int> superAlphabetOrders;
	vector<int> bucketSizes;
	vector<int> markedRowsPercentages;
	int keyIndex;
} TestGroup;



class Test {
public:
	Test(string xmlFilePath);
	virtual ~Test();
	vector<TestGroup> testGroups;
	static PatternSearchInfo *executePatternSearch(u8 *encryptionKey,
				string pattern,string indexFileName,bool loadEntireIndexBeforeSearch);
	static void executePatternSearchInCollection (string indexFileName,u8* encryptionKey, string pattern,bool retrieveSequenceDescription,bool loadEntireIndexBeforeSearch);
private:
	void loadParameters(string xmlFilePath);
};


class TestSingle {
public:
	TestSingle(string xmlFilePath);
	virtual ~TestSingle();
	vector<TestGroup> testGroups;
private:
	void loadParameters(string xmlFilePath);
};





} /* namespace std */

#endif /* TEST_H_ */
