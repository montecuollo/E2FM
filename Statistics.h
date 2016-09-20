/*
 * Statistics.h
 *
 *  Created on: 19/dic/2014
 *      Author: fernando
 */

#ifndef STATISTICS_H_
#define STATISTICS_H_
#include <stdint.h>

namespace std {

class Statistics {
public:
	Statistics();
	virtual ~Statistics();

	uint64_t originalTextLength;
	uint64_t textLength;
	uint32_t alphabetSize;
	uint32_t compactAlphabetSize;
	uint64_t compressedSize;
	double  compressionRatio;
	double  realCompressionRatio;
	uint64_t headerBitmapSize;
	uint64_t headerOccurrencesTableSize;
	uint64_t headerSize;
	uint64_t superbuckets;
	uint64_t buckets;
	uint64_t superbucketAlphabetAverageSize;
	uint64_t superbucketBitmapAverageSize;
	uint64_t superbucketBitmapsCumulativeSize;
	uint64_t superbucketOccurrencesTableAverageSize;
	uint64_t superbucketOccurrencesTablesCumulativeSize;
	uint64_t bucketAlphabetAverageSize;
	uint64_t bucketBitmapAverageSize;
	uint64_t bucketBitmapsCumulativeSize;
	uint64_t bucketOccurrencesTableAverageSize;
	uint64_t bucketOccurrencesTablesCumulativeSize;
	uint64_t bucketMarkedRowsTableAverageSize;
	uint64_t bucketMarkedRowsTablesCumulativeSize;
	//Timings
	uint64_t superAlphabetBuildTime;  //time needed to create the scrambled super alphabet
	uint64_t superTextBuildTime;      //time needed to create the super-text
	uint64_t sortTime;         		//time needed to sort the text suffixes during the computing of BWT
	uint64_t bwtComputationTime;
	uint64_t indexStructureBuildTime;	 //time needed to create the index structure (header, superbuckets and buckets)
	uint64_t bucketEncodingTime;
	uint64_t indexSaveTime;    //time needed to save the index (including Polyalphabetic substitution, MTF and RLE0 for each block)
	double superAlphabetBuildTimesStdDev;  //standard deviations between multiple executions times used to compute superAlphabetBuildTime, which is a mean
	double superTextBuildTimesStdDev;
	double sortTimesStdDev;
	double indexStructureBuildTimesStdDev;
	double indexSaveTimesStdDev;
};

} /* namespace std */

#endif /* STATISTICS_H_ */
