/*
 * EncryptionManager.h
 *
 *  Created on: 20/dic/2014
 *      Author: fernando
 */

#ifndef ENCRYPTIONMANAGER_H_
#define ENCRYPTIONMANAGER_H_
#include "ecrypt-portable.h"
#include "RandomGenerator.h"
#include <vector>
#include <string>

namespace std {

class EncryptionManager {
public:
	EncryptionManager(string &keyFilePath);
	EncryptionManager(u8 *key);
	virtual ~EncryptionManager();
	void initSALSA20Keys(u8 *key);
	uint32_t *computeBucketKeyStream(uint64_t bucketNumber,uint64_t startOffset,uint64_t endOffset,uint32_t bucketAlphaSize);
	vector<uint32_t> *computeScramblingKey(uint32_t alphabetCardinality);
	uint64_t *computeCumulativeFreqsKeyStream(uint32_t keyLength,uint64_t maxFrequence);
	uint32_t *computeFrequenciesPermutation(uint32_t numberOfElements,uint64_t maximumFrequency);
	static void generateKeyFile(string &keyFilePath); //generate a random key of 64 bytes

	friend class EFMIndex;
private:
	void shuffleExcludingFirstElement(vector<uint32_t> *array, uint32_t arraySize,RandomGenerator *rnd);
	void shuffle(uint32_t *array, uint32_t arraySize,RandomGenerator *rnd);
	void swap(vector<uint32_t> *array, uint32_t i, uint32_t j);
	void swap(uint32_t *array, uint32_t i, uint32_t j);
	RandomGenerator *bucketEncryptionRandomGenerator;
	u8 *monoSeed;
	u8 *polySeed;
};

} /* namespace std */

#endif /* ENCRYPTIONMANAGER_H_ */
