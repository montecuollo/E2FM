/*
 * FMCollection.cpp
 *
 *  Created on: 20/gen/2015
 *      Author: fernando
 */

#include "FMCollection.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace std {

ostream& operator<<(ostream& os, const Occurrence& o){
	os << "[" << o.sequenceIndex << "," << o.patternPosition << "]" <<  o.sequenceDescription;
	return os;
}


FMCollection::FMCollection() {
	size=0;
	fastaHeaders=NULL;
	terminatorsPositions=NULL;
}

FMCollection::~FMCollection() {
	// TODO Auto-generated destructor stub
}





/**
 *
 * @return the number of sequences in the original multi-FASTA file (null if the collection is empty or not yet loaded)
 */
uint32_t FMCollection::getSize(){
	return size;
}

void FMCollection::addSequence(string sequenceDescription,string sequenceFileName){
	sequenceDescriptions.push_back(sequenceDescription);
	sequenceFileNames.push_back(sequenceFileName);
	size++;
}


string FMCollection::getSequence(string sequenceDescription){
	if (fastaHeaders==NULL)
		loadFastaHeaders();
	auto it = find(fastaHeaders->begin(), fastaHeaders->end(),sequenceDescription);
	if (it == fastaHeaders->end())
		throw runtime_error("The index doesn't contain a sequence named \""+ sequenceDescription+"\"");
	uint32_t sequenceIndex=it-fastaHeaders->begin();
	return getSequence(sequenceIndex);
}


/**
 * Returns the (FASTA header) description of a sequence identified by its index
 * @param sequenceIndex
 * @return
 * @throws Exception
 */
string FMCollection::getDescription(uint32_t sequenceIndex){
	if (fastaHeaders==NULL)
		loadFastaHeaders();
	fastaHeaders=new vector<string>();
	boost::split(*fastaHeaders, headersBuffer, boost::is_any_of("\n"));
	return fastaHeaders->at(sequenceIndex);
}


void FMCollection::saveToFasta(string fastaFileName){
	ofstream out(fastaFileName, std::ios::out | std::ios::binary);
	for (uint32_t i=0;i<size;i++){
		string sequence=getSequence(i);
		string sequenceDescription=getDescription(i);
		out << ">";
		out << sequenceDescription;
		out << "\n";
		out << sequence;
		out <<"\n";
	}
	out.close();
}






} /* namespace std */
