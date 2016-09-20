/*
 * FMCollection.h
 *
 *  Created on: 20/gen/2015
 *      Author: fernando
 */

#ifndef FMCOLLECTION_H_
#define FMCOLLECTION_H_

#include <string>
#include <vector>

#define fastaHeadersTerminator "\n";

namespace std {

/**
 *
 *
 * Occurrence of a pattern within a multi-Fasta file, containing a collection of sequences
 */
typedef struct  Occurrence {
	/**
	 *
	 * @param sequenceIndex index of the sequence within the FASTA file
	 * @param sequenceDescription	description of the sequence as stored in the relative FASTA header line
	 * @param patternPosition	position (0-based) of the searched pattern within the sequence
	 */
	Occurrence(uint32_t sequenceIndex,string sequenceDescription,uint64_t patternPosition){
		this->sequenceIndex=sequenceIndex;
		this->sequenceDescription=sequenceDescription;
		this->patternPosition=patternPosition;
	}

	friend ostream& operator<<(ostream& os, const Occurrence& o);


	uint32_t sequenceIndex;
	string sequenceDescription;
	uint64_t patternPosition;


} Occurrence;


typedef struct SearchResult {
	SearchResult(){

	}

	SearchResult(string pattern,double elapsedTime,vector<Occurrence> occurrences){
		this->pattern=pattern;
		this->elapsedTime=elapsedTime;
		this->occurrences=occurrences;
	}

	string pattern;
	double elapsedTime;
	vector<Occurrence> occurrences;
} SearchResult;



class FMCollection {
public:
	FMCollection();
	virtual ~FMCollection();
	void loadFasta(string fastaFileName, int k);
	virtual void build(string &fastaFilePath,string &indexFilePath)=0;
	virtual void build(string &indexFileName)=0;
	uint32_t getSize();
	void addSequence(string sequenceDescription,string sequenceFileName);
	virtual string getSequence(uint32_t sequenceIndex)=0;
	virtual string extract(uint32_t sequenceIndex,uint64_t from,uint64_t to)=0;
	string getSequence(string sequenceDescription);
	string getDescription(uint32_t sequenceIndex);
	virtual void loadFastaHeaders()=0;
	virtual SearchResult *search (string pattern,bool retrieveSequenceDescriptions)=0;
	void saveToFasta(string fastaFileName);
	//virtual void setMarkedRowsPercentage(uint8_t markedRowsPercentage)=0;
	virtual void setBucketSize(uint64_t bucketSize)=0;
	virtual void load(string fileName)=0;
protected:
	vector<int64_t> * terminatorsPositions;
	vector<string> *fastaHeaders;
	uint32_t size;
	vector<string> sequenceFileNames;
	vector<string> sequenceDescriptions;

	string headersBuffer;
	string buffer;


};

} /* namespace std */

#endif /* FMCOLLECTION_H_ */
