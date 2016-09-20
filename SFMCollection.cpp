/*
 * SFMCollection.cpp
 *
 *  Created on: 20/gen/2015
 *      Author: fernando
 */

#include "SFMCollection.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "Utils.h"
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io/sequence_file.h>



using namespace sdsl;

namespace std {

SFMCollection::SFMCollection() : FMCollection(){

}

SFMCollection::~SFMCollection() {

}

/**
 * Builds an index on the collection of sequences (DNA, RNA or protein sequences) yet added
 * with a series of addSequence invocations
 *
 * @param fileName, output file
 * @throws Exception
 */
void SFMCollection::build(string &fastaFilePath,string &indexFilePath){
	loadSequencesFromFasta(fastaFilePath);
	string memoryBufferFileName = "@build.buffer" + indexFilePath;
	store_to_file((const char*)buffer.c_str(), memoryBufferFileName);
	construct(fm_index, memoryBufferFileName, 1);
	store_to_file(fm_index, indexFilePath);
	//EFMIndex::build(buffer,fileName);
}

void SFMCollection::load(string fileName){
	load_from_file(fm_index, fileName);
	retrieveTerminatorsPositions();
}


void SFMCollection::build(string &indexFilePath){
	//not implemented
}


SearchResult *SFMCollection::search (string pattern,bool retrieveSequenceDescriptions){
	std::chrono::time_point<std::chrono::system_clock> startTime,endTime;
	startTime=std::chrono::system_clock::now();
	SearchResult *searchResult=new SearchResult();
	searchResult->pattern=pattern;
	uint32_t m=pattern.size();
	auto locations = locate(fm_index, pattern.begin(), pattern.begin()+m);
	sort(locations.begin(), locations.end());
	uint64_t nOccs=locations.size();
	for (uint64_t i=0;i<nOccs;i++){
		uint64_t occ=locations[i];
		int64_t sequenceIndex=-1;
		uint32_t j=0;
		uint32_t l=terminatorsPositions->size();
		while (j<l && sequenceIndex==-1){
			if (terminatorsPositions->at(j)>occ)
				sequenceIndex=j;
			else
				j++;
		}

		string sequenceDescription;
		if (retrieveSequenceDescriptions)
			sequenceDescription=getDescription(sequenceIndex);
		else
			sequenceDescription="";

		int64_t startingPosition;
		if (sequenceIndex>0)
			startingPosition=terminatorsPositions->at(sequenceIndex-1)+1;
		else
			startingPosition=0;

		searchResult->occurrences.push_back(Occurrence(sequenceIndex,sequenceDescription, occ-startingPosition));
	}

	endTime=std::chrono::system_clock::now();
	searchResult->elapsedTime=(double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()/(double)1000;
	return searchResult;
}



void SFMCollection::loadSequencesFromFasta(string &fastaFilePath){
	cout << "Loading sequences in memory" <<endl;
	stringstream ssText;
	stringstream ssHeaders;
	int64_t pos=0;
	int64_t missingToKAlignment;
	string terminator( 1,'&');
	string sequence;
	terminatorsPositions=new vector<int64_t>();

	seqan::SeqFileIn seqFileIn(fastaFilePath.c_str());
	seqan::StringSet<seqan::CharString> ids;
	seqan::StringSet<seqan::IupacString> seqs;
	seqan::readRecords(ids, seqs, seqFileIn);

	size=length(seqs);
	for (unsigned i = 0; i < length(seqs); ++i){
		 string sequenceDescription(seqan::toCString(ids[i]));
		 sequenceDescriptions.push_back(sequenceDescription);

		 ssHeaders << sequenceDescription << fastaHeadersTerminator;
		 ssText << seqs[i];
		 pos+=length(seqs[i]);
		 terminatorsPositions->push_back(pos);
		 ssText << terminator;
		 pos+=1;
	}
	cout << "\tloaded " << size << " sequences" <<endl;
	buffer=ssText.str();
	//cout << buffer << endl;
	cout << "\tconcatenated size:" << buffer.size() <<endl;
	headersBuffer=ssHeaders.str();

}




void SFMCollection::loadSequences(){
	cout << "Loading sequences in memory" <<endl;
	stringstream ssText;
	stringstream ssHeaders;
	int64_t pos=0;
	int64_t missingToKAlignment;
	string terminator( 1,'&');
	string sequence;
	terminatorsPositions=new vector<int64_t>();
	for (uint32_t i=0;i<size;i++){
		 ssHeaders << sequenceDescriptions[i] << fastaHeadersTerminator;
		 sequence=Utils::loadFasta(sequenceFileNames[i]);
		 ssText << sequence;
		 pos+=sequence.size();
		 terminatorsPositions->push_back(pos);
		 ssText << terminator;
		 pos+=1;
	}
	cout << "\tloaded " << size << " sequences" <<endl;
	buffer=ssText.str();
	cout << "\tconcatenated size:" << buffer.size() <<endl;
	headersBuffer=ssHeaders.str();
}

/**
 * Returns the sequence corresponding to the index sequenceIndex
 * @param sequenceIndex
 * @return
 * @throws Exception
 */
string SFMCollection::getSequence(uint32_t sequenceIndex){
	uint64_t from;
	if (sequenceIndex>0)
		from=terminatorsPositions->at(sequenceIndex-1)+1;
	else
		from=0;
	uint64_t to=terminatorsPositions->at(sequenceIndex);
	//TODo: implementare l'estrazione dello snippet
	//string snippet=extractSnippet(from, to);
	//int64_t p=snippet.find('&');
	//return snippet.substr(0,p);
	return "";
}

/**
 * Returns the snippet in the  positions [from,to] of the sequence corresponding to sequenceIndex
 * @param sequenceIndex
 * @return
 * @throws Exception
 */
string SFMCollection::extract(uint32_t sequenceIndex,uint64_t from,uint64_t to){


	uint64_t wholeSequenceFrom;  //absolute: relative to the full EFM-index text
	if (sequenceIndex>0)
		wholeSequenceFrom=terminatorsPositions->at(sequenceIndex-1)+1;
	else
		wholeSequenceFrom=0;
	uint64_t wholeSequenceTo=terminatorsPositions->at(sequenceIndex);

	uint64_t sequenceMaximumLength=wholeSequenceTo-wholeSequenceFrom+1;  //maximum length can be greater of effective length at most of k-1 characters
	if (from>sequenceMaximumLength-1)
		throw runtime_error("start index out of range");
	if (to>sequenceMaximumLength-1)
			throw runtime_error("end index out of range");


	uint64_t absoluteFrom=wholeSequenceFrom+from;
	uint64_t absoluteTo=wholeSequenceFrom+to;


	//TODO: implementare l'estrazione dello snippet
	//string snippet=extractSnippet(absoluteFrom,sequenceRelativeFromOffset,absoluteTo,sequenceRelativeToOffset);
	//int64_t p=snippet.find('&');
	//if (p>-1)
	//	return snippet.substr(0,p);
	//else
	//	return snippet;
}




void SFMCollection::setBucketSize(uint64_t bucketSize){
	this->bucketSize = bucketSize;
}

void SFMCollection::setMarkedRowsPercentage(uint8_t markedRowsPercentage) {
	this->markedRowsPercentage = markedRowsPercentage;
}

void SFMCollection::loadFastaHeaders(){
	//TODO: implementare
	/*
	bitReader->gotoBookmark(fastaHeadersStart);
	//reads the FASTA headers from disk
	uint32_t fastaHeadersSize=bitReader->getInt();
	headersBuffer.resize(fastaHeadersSize);
	for (uint32_t i=0;i<fastaHeadersSize;i++)
		headersBuffer[i]= bitReader->getUByte();
	*/
}

void SFMCollection::retrieveTerminatorsPositions(){
		string terminator="&";
		terminatorsPositions=new vector<int64_t>();
		auto locations = locate(fm_index, terminator.begin(), terminator.begin()+1);
		sort(locations.begin(), locations.end());
		uint32_t nOccs=locations.size();
		for (uint32_t i=0;i<nOccs;i++){
			int64_t occ=locations[i];
			terminatorsPositions->push_back(occ);
		}
}

void SFMCollection::close(){

}



} /* namespace std */
