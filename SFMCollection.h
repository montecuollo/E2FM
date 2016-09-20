/*
 * SFMCollection.h
 *
 *  Created on: 20/gen/2015
 *      Author: fernando
 */

#ifndef SFMCOLLECTION_H_
#define SFMCOLLECTION_H_

#include "FMCollection.h"
#include <string>
#include <sdsl/suffix_trees.hpp>

using namespace sdsl;

namespace std {

class SFMCollection : public FMCollection{
public:
	SFMCollection();
	virtual ~SFMCollection();
	virtual void setBucketSize(uint64_t bucketSize);
	virtual void setMarkedRowsPercentage(uint8_t markedRowsPercentage);
	virtual void build(string &indexFilePath);
	virtual void build(string &fastaFilePath,string &indexFilePath);
	virtual void load(string fileName);
	virtual string getSequence(uint32_t sequenceIndex);
	string extract(uint32_t sequenceIndex,uint64_t from,uint64_t to);
	virtual SearchResult *search (string pattern,bool retrieveSequenceDescriptions);
	void close();
private:
	void loadSequences();
	void loadSequencesFromFasta(string &fastaFilePath);
	void retrieveTerminatorsPositions();
	virtual void loadFastaHeaders();

	uint64_t bucketSize;
	uint8_t markedRowsPercentage;
	csa_wt<wt_huff<rrr_vector<127> >,50 , 1024> fm_index;
};

} /* namespace std */

#endif /* SFMCOLLECTION_H_ */
