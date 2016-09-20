/*
 * EFMCollection.cpp
 *
 *  Created on: 08/gen/2015
 *      Author: fernando
 */
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include "EFMCollection.h"


namespace std {



EFMCollection::EFMCollection():FMCollection() {
	fastaHeadersStart=NULL;
}

EFMCollection::~EFMCollection() {
	delete fastaHeaders;
	delete terminatorsPositions;
	delete fastaHeadersStart;
}


void EFMCollection::setBucketSize(uint64_t bucketSize){
	this->bucketSize = bucketSize;
}


void EFMCollection::setMarkedRowsPercentage(uint8_t markedRowsPercentage) {
	this->markedRowsPercentage = markedRowsPercentage;
}


/**
 * Returns the sequence corresponding to the index sequenceIndex
 * @param sequenceIndex
 * @return
 * @throws Exception
 */
string EFMCollection::getSequence(uint32_t sequenceIndex){
	uint64_t from;
	if (sequenceIndex>0)
		from=terminatorsPositions->at(sequenceIndex-1)+superAlphabetOrder;
	else
		from=0;
	uint64_t to=terminatorsPositions->at(sequenceIndex);
	string snippet=extractSnippet(from, to);
	uint64_t p=snippet.find('&');
	if (p>-1)
		return snippet.substr(0,p);
	else
		return snippet;
}




/**
 * Returns the snippet in the  positions [from,to] of the sequence corresponding to sequenceIndex
 * @param sequenceIndex
 * @return
 * @throws Exception
 */
string EFMCollection::extract(uint32_t sequenceIndex,uint64_t from,uint64_t to){
	int64_t wholeSequenceFrom;  //absolute: relative to the full EFM-index text
	if (sequenceIndex>0)
		wholeSequenceFrom=terminatorsPositions->at(sequenceIndex-1)+superAlphabetOrder;
	else
		wholeSequenceFrom=0;
	int64_t wholeSequenceTo=terminatorsPositions->at(sequenceIndex);

	int64_t sequenceMaximumLength=wholeSequenceTo-wholeSequenceFrom+1;  //maximum length can be greater of effective length at most of k-1 characters
	if (from>sequenceMaximumLength-1)
		throw runtime_error("start index out of range");
	if (to>sequenceMaximumLength-1)
			throw runtime_error("end index out of range");



	//computes from in terms of super text and from offset
	int64_t f=wholeSequenceFrom+from;
	int64_t t=wholeSequenceFrom + to;

	string snippet=extractSnippet(f,t);
	int64_t p=snippet.find('&');
	if (p>-1)
		return snippet.substr(0,p);
	else
		return snippet;
}



void EFMCollection::loadFastaHeaders(){
	bitReader->gotoBookmark(fastaHeadersStart);
	//reads the FASTA headers from disk
	uint32_t fastaHeadersSize=bitReader->getInt();
	headersBuffer.resize(fastaHeadersSize);
	for (uint32_t i=0;i<fastaHeadersSize;i++)
		headersBuffer[i]= bitReader->getUByte();
}



/**
 * Searchs for a pattern within the multi-FASTA file
 * @param pattern
 * @param retrieveSequenceDescriptions
 * @return the search result
 * @throws Exception
 */
SearchResult *EFMCollection::search (string pattern,bool retrieveSequenceDescriptions){

	//uint64_t rp=getRowPosition(54009,bitReader,bucketEncryptionRandomGenerator);

	std::chrono::time_point<std::chrono::system_clock> startTime,endTime;
	startTime=std::chrono::system_clock::now();
	SearchResult *searchResult=new SearchResult();
	searchResult->pattern=pattern;
	//if (pattern=="AATGGCGAATAAAAAGCTGC")
	//	cout <<endl;
	exactMatch(pattern);
	vector<uint64_t> *occs=locateOccurrences(&allSuperPatternMatches);
	for (uint64_t occ:*occs){
		int64_t sequenceIndex=-1;
		uint32_t i=0;
		uint32_t l=terminatorsPositions->size();
		while (i<l && sequenceIndex==-1){
			if (terminatorsPositions->at(i)>occ)
				sequenceIndex=i;
			else
				i++;
		}

		string sequenceDescription;
		if (retrieveSequenceDescriptions)
			sequenceDescription=getDescription(sequenceIndex);
		else
			sequenceDescription="";

		uint64_t startingPosition;
		if (sequenceIndex>0)
			startingPosition=terminatorsPositions->at(sequenceIndex-1)+superAlphabetOrder;
		else
			startingPosition=0;

		searchResult->occurrences.push_back(Occurrence(sequenceIndex,sequenceDescription, occ-startingPosition));
	}

	endTime=std::chrono::system_clock::now();
	searchResult->elapsedTime=(double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()/(double)1000;
	return searchResult;
}

void EFMCollection::loadAddedSequences(uint8_t k){
	cout << "Loading sequences in memory" <<endl;
	stringstream ssText;
	stringstream ssHeaders;
	int64_t pos=0;
	int64_t missingToKAlignment;
	string terminator( k,'&');
	string sequence;
	terminatorsPositions=new vector<int64_t>();
	for (uint32_t i=0;i<size;i++){
		 ssHeaders << sequenceDescriptions[i] << fastaHeadersTerminator;
		 sequence=Utils::loadFasta(sequenceFileNames[i]);
		 ssText << sequence;
		 pos+=sequence.size();

		 missingToKAlignment=pos%k;
		 if (missingToKAlignment>0)
			 missingToKAlignment=k-missingToKAlignment;
		 pos+=missingToKAlignment;
		 terminatorsPositions->push_back(pos);
		 string missingCharacters(missingToKAlignment,'&');
		 ssText << missingCharacters << terminator;
		 pos+=k;
	}
	cout << "\tloaded " << size << " sequences" <<endl;
	buffer=ssText.str();
	cout << "\tconcatenated size:" << buffer.size() <<endl;
	headersBuffer=ssHeaders.str();

}


void EFMCollection::loadSequencesFromFasta(string &fastaFilePath,uint8_t k){
	std::chrono::time_point<std::chrono::system_clock> encodingStartTime=std::chrono::system_clock::now();
	std:time_t printableTime= std::chrono::system_clock::to_time_t(encodingStartTime);
	char mbstr[80];
	std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
	cout << "Loading sequences in memory, starting at " << mbstr <<endl;

	vector<string> fileNames;
	stringstream ssText;
	stringstream ssHeaders;
	int64_t pos=0;
	int64_t missingToKAlignment;
	string terminator( k,'&');
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

		 missingToKAlignment=pos%k;
		 if (missingToKAlignment>0)
			 missingToKAlignment=k-missingToKAlignment;
		 pos+=missingToKAlignment;
		 terminatorsPositions->push_back(pos);
		 string missingCharacters(missingToKAlignment,'&');
		 ssText << missingCharacters << terminator;
		 pos+=k;
	}
	cout << "\tloaded " << size << " sequences" <<endl;
	buffer=ssText.str();
	//cout << buffer << endl;
	cout << "\tconcatenated size:" << buffer.size() <<endl;
	headersBuffer=ssHeaders.str();

}


void EFMCollection::writeSuspendedInformations(){
	EFMIndex::writeSuspendedInformations();
	bitWriter->writeInt64(fastaHeadersStart->getPosition());
}


void EFMCollection::readSuspendedInformations(){
	EFMIndex::readSuspendedInformations();
	fastaHeadersStart=new Bookmark(bitReader->getInt64(),0);

}

void EFMCollection::createSpaceForSuspendedInfos(){
	EFMIndex::createSpaceForSuspendedInfos();
	//Free space for the start byte of the compressed BZIP2 stream containing the FASTA headers
	bitWriter->writeInt64(0);
}

void EFMCollection::writeHeader(){
	EFMIndex::writeHeader();
	//writes the terminators positions to disk
	vector<int64_t> diffs;
	uint64_t maxDiff=0;
	uint64_t ptp=0;
	uint64_t tp;
	uint64_t diff;
	for (uint32_t i=0;i<terminatorsPositions->size();i++){
		tp=terminatorsPositions->at(i);
		diff=tp-ptp;
		ptp=tp;
		if (diff>maxDiff)
			maxDiff=diff;
		diffs.push_back(diff);
	}
	//Writes to disk the number of the collection elements (size of the diffs list)
	Bookmark* b=bitWriter->getBookmark();
	bitWriter->writeInt(diffs.size());
	uint8_t neededBits=Utils::int_log2(maxDiff); // (int)Math.ceil((Math.log(Math.ceil(maxDiff +1))/Math.log(2)));
	bitWriter->writeUByte(neededBits);
	for (uint32_t i=0;i<diffs.size();i++)
		bitWriter->write(neededBits, diffs[i]);

	bitWriter->flush();
	fastaHeadersStart=bitWriter->getBookmark();

	//writes FASTA headers to disk
	uint32_t fastaHeadersSize=headersBuffer.size();
	bitWriter->writeInt(fastaHeadersSize);
	for (uint32_t i=0;i<fastaHeadersSize;i++)
			bitWriter->writeUByte(headersBuffer[i]);

}

void EFMCollection::readHeader(){
	EFMIndex::readHeader();
	Bookmark *b=bitReader->getBookmark();

	// reads the terminators positions from disk
	terminatorsPositions=new vector<int64_t>();
	size=bitReader->getInt();
	uint8_t neededBits=bitReader->getUByte();
	int64_t tp=0;
	int64_t diff=0;
	for (uint32_t i=0;i<size;i++){
		diff=bitReader->read(neededBits);
		tp=diff+tp;
		terminatorsPositions->push_back(tp);
	}


	//readText();
	//reconstructMarkedRows();
}


/**
 * Builds an index on the collection of sequences (DNA, RNA or protein sequences) yet added
 * with a series of addSequence invocations
 *
 * @param fileName, output file
 */
void EFMCollection::build(string &indexFileName){
	loadAddedSequences(superAlphabetOrder);
	EFMIndex::build(buffer,fileName);
}

/**
 * Builds an index on the collection of sequences (DNA, RNA or protein sequences) contained in a FASTA file
 *
 * @param fileName, output file
 */
void EFMCollection::build(string &fastaFilePath,string &indexFileName){
	loadSequencesFromFasta(fastaFilePath,superAlphabetOrder);
	EFMIndex::build(buffer,indexFileName);
}


void EFMCollection::load(string fileName){
	EFMIndex::load(fileName);
}




} /* namespace std */
