/*
 * EFMCollection.h
 *
 *  Created on: 08/gen/2015
 *      Author: fernando
 */

#ifndef EFMCOLLECTION_H_
#define EFMCOLLECTION_H_

#include <string>
#include "EFMIndex.h"
#include "Bookmark.h"
#include "FMCollection.h"
#include <algorithm>
#include <chrono>
#include "Utils.h"
#include <boost/algorithm/string.hpp>

namespace std {

class EFMCollection : public EFMIndex, public FMCollection{
public:
	EFMCollection();
	virtual ~EFMCollection();

	void loadFasta(string fastaFileName, uint8_t k);
	virtual void build(string &indexFileName);
	virtual void build(string &fastaFilePath,string &indexFileName);
	virtual void load(string fileName);
	virtual string getSequence(uint32_t sequenceIndex);
	virtual string extract(uint32_t sequenceIndex,uint64_t from,uint64_t to);
	virtual void loadFastaHeaders();
	virtual SearchResult *search (string pattern,bool retrieveSequenceDescriptions);
	virtual void setBucketSize(uint64_t bucketSize);
	virtual void setMarkedRowsPercentage(uint8_t markedRowsPercentage);

private:
	virtual void writeSuspendedInformations();
	virtual void readSuspendedInformations();
	void createSpaceForSuspendedInfos();
	virtual void writeHeader();
	virtual void readHeader();
	void loadAddedSequences(uint8_t k);
	void loadSequencesFromFasta(string &fastaFilePath,uint8_t k);

	Bookmark *fastaHeadersStart;

};

} /* namespace std */

#endif /* EFMCOLLECTION_H_ */
