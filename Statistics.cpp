/*
 * Statistics.cpp
 *
 *  Created on: 19/dic/2014
 *      Author: fernando
 */

#include "Statistics.h"

namespace std {

Statistics::Statistics() {
	superAlphabetBuildTime=0;
	superTextBuildTime=0;
	sortTime=0;
	indexStructureBuildTime=0;
	indexSaveTime=0;
	bucketEncodingTime=0;
	bwtComputationTime=0;

}

Statistics::~Statistics() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
