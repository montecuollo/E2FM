/*
 * EncryptionManager.cpp
 *
 *  Created on: 20/dic/2014
 *      Author: fernando
 */

#include "EncryptionManager.h"
#include <algorithm>
#include "RandomGenerator.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
namespace std {


EncryptionManager::EncryptionManager(string &keyFilePath) {
	std::ifstream is (keyFilePath, std::ifstream::binary);
	if (is) {
		// opens the key file and verifies it
		is.seekg (0, is.end);
		uint64_t length = is.tellg();
		is.seekg (0, is.beg);

		if (length !=64)
				throw runtime_error("File " + keyFilePath + "is not a valid key file");

		u8 * key = new u8[length];
		is.read ((char *)key,length);
		is.close();

		initSALSA20Keys(key);

		delete[] key;
	} else
		throw runtime_error ("Error opening key file "+ keyFilePath);
}

EncryptionManager::EncryptionManager(u8 *key) {
	initSALSA20Keys(key);
}

/**
 * Splits the encryption key into two SALSA20 512 bit keys: the first will be used for scrambling, while the second one
 * will be used to cipher the index blocks
 */
void EncryptionManager::initSALSA20Keys(u8 *key){
	bucketEncryptionRandomGenerator=NULL;
	monoSeed=new u8[32];
	polySeed=new u8[32];
	//the first part of the key is used to scramble the extended alphabet
	std::copy(key, key+32, monoSeed);
	//the second part of the key is used in bucket texts encodings
	std::copy(key+32, key+64, polySeed);
}

EncryptionManager::~EncryptionManager() {
	delete[] monoSeed;
	delete[] polySeed;
	//if (bucketEncryptionRandomGenerator!=NULL)
	//	delete bucketEncryptionRandomGenerator; //TODO: capire perchè produce Segmentation fault
}


uint32_t *EncryptionManager::computeBucketKeyStream(uint64_t bucketNumber,uint64_t startOffset,uint64_t endOffset,uint32_t bucketAlphaSize){
		uint64_t keyFragmentLength=endOffset-startOffset+1;
		uint32_t* key=new uint32_t[keyFragmentLength];

		//Initialize a pseudo-casual number generator
		uint32_t maxValue=bucketAlphaSize;  //this upper bound is included in the range of generated numbers
		if (bucketEncryptionRandomGenerator==NULL)
				bucketEncryptionRandomGenerator=new RandomGenerator(131072,polySeed);

		bucketEncryptionRandomGenerator->populateKeyStreamBuffer(bucketNumber,maxValue,startOffset,endOffset);
		for (uint64_t i=0;i<keyFragmentLength;i++)
			key[i]=bucketEncryptionRandomGenerator->nextInt();
		return key;
}


/**
 * Computes the scrambling key using the Fisher-Yates algorithm
 * (Knuth shuffle) similar to that provided by standard Collections.shuffle.
 * The algorithm proceeds as follows:
 * 1) it creates an array of dimension equal to the super alphabet size, initialized
 *    with the numbers from 0 to superAlphabetSize -1
 * 2) it shuffles the array's elements by the Fisher-Yates algorithm (Knuth shuffle),
 *    excluding its first element.
 * @throws Exception
 *
 */
vector<uint32_t> *EncryptionManager::computeScramblingKey(uint32_t alphabetCardinality){
	vector<uint32_t> *scramblingKey=new vector<uint32_t>();
	scramblingKey->resize(alphabetCardinality);
	for (uint32_t i=0;i<alphabetCardinality;i++)
		(*scramblingKey)[i]=i;
	//Initialize a pseudo-casual number generator using key as seed
	RandomGenerator *random=new RandomGenerator(16384,monoSeed,0);
	//Do a Knuth shuffle in order to obtain a random scrambling key
	shuffleExcludingFirstElement(scramblingKey,alphabetCardinality,random);
	delete random;
	return scramblingKey;
}



void EncryptionManager::shuffleExcludingFirstElement(vector<uint32_t> *array,uint32_t arraySize,RandomGenerator *rnd){
	// Shuffles the array excluding the first element
	for (uint32_t i=arraySize; i>2; i--){
		//computes casually the index of the element to swap with the i^th element of the array
		uint32_t toSwapWith;
		do
		  toSwapWith=rnd->nextInt(i);
		while (toSwapWith==0);
		swap(array, i-1, toSwapWith);
	}
}

void EncryptionManager::shuffle(uint32_t *array, uint32_t arraySize,RandomGenerator *rnd){
        // Shuffles the array excluding the first element
        for (uint32_t i=arraySize; i>1; i--){
        	//computes casually the index of the element to swap with the i^th element of the array
        	uint32_t toSwapWith;
        	do
        	  toSwapWith=rnd->nextInt(i);
        	while (toSwapWith==0);
            swap(array, i-1, toSwapWith);
        }
}



inline void EncryptionManager::swap(vector<uint32_t> *array, uint32_t i, uint32_t j) {
        uint32_t tmp = (*array)[i];
        (*array)[i] = (*array)[j];
        (*array)[j] = tmp;
}

inline void EncryptionManager::swap(uint32_t *array, uint32_t i, uint32_t j) {
        uint32_t tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
}

/*
	 * @param keyLength  Key length
	 * @param maxFrequence  maximum cumulative frequency of the alphabet's characters
	 * @return
	 */
	uint64_t *EncryptionManager::computeCumulativeFreqsKeyStream(uint32_t keyLength,uint64_t maxFrequence){
		uint64_t *key=new uint64_t[keyLength];

		//Initialize a pseudo-casual number generator, adding the alphabet cardinality
		//to the polySeed: this tries to add to the seed a dependency from the data
		//RandomGenerator randomGenerator=new Salsa20NativeRNG(salsa20KeyStreamStandardBufferSize,polySeed,2);
		RandomGenerator *random=new RandomGenerator(16384,polySeed,2);

		for (uint32_t i=0;i<keyLength;i++)
			key[i]=random->nextInt(maxFrequence+1);
		return key;
	}


	/**
	 *
	 * @param numberOfElements  number of elements of the cumulative frequencies vector
	 * @param maximumFrequency  value of the greatest cumulative frequency
	 * @return
	 */
	uint32_t *EncryptionManager::computeFrequenciesPermutation(uint32_t numberOfElements,uint64_t maximumFrequency){
		uint32_t *permutation=new uint32_t[numberOfElements];
		for (uint32_t i=0;i<numberOfElements;i++)
			permutation[i]=i;
		//Initialize a pseudo-casual number generator using key as seed
		RandomGenerator *random=new RandomGenerator(16384, monoSeed,1);
		//Do a Knuth shuffle in order to generate a random monoAlphabeticKey
		shuffle(permutation,numberOfElements,random);
		return permutation;
	}



	void EncryptionManager::generateKeyFile(string &keyFilePath) {
		std::random_device generator;

		u8* key=new u8[64];
		std::uniform_int_distribution<u8> distribution(0,255);

		for (uint32_t i=0;i<64;i++)
			key[i]=distribution(generator);

		ofstream keyFile(keyFilePath, ios::out | ios::binary);
		keyFile.write ((char*)key, 64);
		keyFile.close();
		return;
	}


} /* namespace std */
