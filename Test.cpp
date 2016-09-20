/*
 * Test.cpp
 *
 *  Created on: 31/dic/2014
 *      Author: fernando
 */

#include "Test.h"
#include <cmath>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <iostream>

namespace std {

Test::Test(string xmlFilePath) {
	loadParameters(xmlFilePath);
}

Test::~Test() {
	// TODO Auto-generated destructor stub
}


void Test::loadParameters(string xmlFilePath){
	boost::property_tree::ptree pTree;
	boost::property_tree::read_xml(xmlFilePath,pTree);
	//read the test parameters from the XML file
	//Using boost::property_tree
	for( auto const& elementNode : pTree.get_child("tns:parameters") ) {
		//string s=elementNode.first;
		boost::property_tree::ptree groupTree = elementNode.second;
		TestGroup group;
		group.inputFileType = groupTree.get<string>("tns:inputFileType");
		group.sequencesDirectory = groupTree.get<string>("tns:sequencesDirectory");
		group.referenceIdentifier = groupTree.get<string>("tns:referenceIdentifier");
		group.indexFileName = groupTree.get<string>("tns:indexFileName");
		group.naiveSearch = groupTree.get<bool>("tns:naiveSearch");
		group.loadIndexInMemory = groupTree.get<bool>("tns:loadIndexInMemory");
		for( auto const& groupChildNode : groupTree.get_child("") ) {
			if (groupChildNode.first=="tns:patternSets"){
				PatternSet patternSet;
				boost::property_tree::ptree patternSetSubTree = groupChildNode.second;
				patternSet.length= patternSetSubTree.get<int>("tns:length");
				for( auto const& patternSetChildNode : patternSetSubTree.get_child("") ) {
					if (patternSetChildNode.first=="tns:patterns"){
						//string s=patternSetChildNode.first;
						string p=patternSetChildNode.second.data();
						patternSet.patterns.push_back(p);
					}
				}
				group.patternSets.push_back(patternSet);
			}
		}
		testGroups.push_back(group);
	}
}

//RIMUOVERE: NON FUNZIONA PIU' (RICERCHE SUPPORTATE SOLO SU COLLEZIONI DI SEQUENZE)
PatternSearchInfo *Test::executePatternSearch(u8* encryptionKey,
			string pattern,string indexFileName,bool loadEntireIndexBeforeSearch){
	PatternSearchInfo *info=new PatternSearchInfo();
	try{
		EFMIndex *indexHandler=new EFMIndex();
		indexHandler->setEncryptionKey(encryptionKey);
		std::chrono::time_point<std::chrono::system_clock> startTime,endTime;
		startTime=std::chrono::system_clock::now();
		indexHandler->load(indexFileName);
		if (loadEntireIndexBeforeSearch)
			indexHandler->readText();
		endTime=std::chrono::system_clock::now();
		info->loadingTime= std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
		info->numberOfBuckets=(double)indexHandler->getBucketsNumber();
		info->numberOfSuperBuckets=(double)indexHandler->getSuperBucketsNumber();
		startTime = std::chrono::system_clock::now();
		//vector<Match*> *matches=indexHandler->exactMatchOld(pattern);
		vector<Match*> *matches;
		info->numberOfOccurrences=indexHandler->countOccurrences(matches);
		endTime   = std::chrono::system_clock::now();
		info->countTime=(double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()/(double)1000;
 		info->countMeanTime=info->countTime/(double)info->numberOfOccurrences;
		info->numberOfLoadedBuckets=indexHandler->getNumberOfLoadedBuckets();
		info->numberOfLoadedSuperBuckets=indexHandler->getNumberOfLoadedSuperBuckets();
		try{
			startTime = std::chrono::system_clock::now();
			vector<uint64_t> *occs=  indexHandler->locateOccurrences(matches);
			//indexHandler->extractSnippet(389,399);
			endTime   = std::chrono::system_clock::now();
			info->occurrences=occs;
			info->numberOfLoadedBuckets=indexHandler->getNumberOfLoadedBuckets();
			info->numberOfLoadedSuperBuckets=indexHandler->getNumberOfLoadedSuperBuckets();
			info->locateSupported=true;
			info->locateTime=(double)std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count()/(double)1000;
			info->searchTime=info->countTime+info->locateTime;
			info->locateMeanTime=info->locateTime/info->numberOfOccurrences;
			info->searchMeanTime=info->searchTime/info->numberOfOccurrences;
		} catch(LocateNotSupportedException ex){
			info->locateSupported=false;
		}

		indexHandler->close();
		indexHandler=NULL;

	}
	catch(exception ex){
		cout << "Errore nella ricerca del pattern: " << ex.what() << endl;
	}
	return info;
}


void Test::executePatternSearchInCollection (string indexFileName,u8* encryptionKey, string pattern,bool retrieveSequenceDescription,bool loadEntireIndexBeforeSearch){
	cout << "Search of pattern " << pattern;
	EFMCollection *indexHandler=new EFMCollection();
	indexHandler->setEncryptionKey(encryptionKey);
	indexHandler->load(indexFileName);
	if (loadEntireIndexBeforeSearch)
		indexHandler->readText();
	SearchResult *results=indexHandler->search(pattern,retrieveSequenceDescription);
	cout << "\t Locate time: " << results->elapsedTime << "ms" << "\t" << "occurrences: " << results->occurrences.size()  << endl;
	for (Occurrence r:results->occurrences)
		cout <<"\t" << r << endl;

	//cout << "Loaded buckets: " << indexHandler->getNumberOfLoadedBuckets() << "/" << indexHandler->getBucketsNumber() << endl;
	//cout << "Elapsed time: " << results->elapsedTime << endl;

}


TestSingle::TestSingle(string xmlFilePath) {
	loadParameters(xmlFilePath);
}

TestSingle::~TestSingle() {
	// TODO Auto-generated destructor stub
}


void TestSingle::loadParameters(string xmlFilePath){
	boost::property_tree::ptree pTree;
	boost::property_tree::read_xml(xmlFilePath,pTree);
	//read the test parameters from the XML file
	//Using boost::property_tree
	for( auto const& elementNode : pTree.get_child("tns:parameters") ) {
		//string s=elementNode.first;
		boost::property_tree::ptree groupTree = elementNode.second;
		TestGroup group;
		group.inputFileType = groupTree.get<string>("tns:inputFileType");
		group.inputFileName = groupTree.get<string>("tns:inputFileName");
		group.indexParentDirectory = groupTree.get<string>("tns:indexParentDirectory");
		group.buildIndex = groupTree.get<bool>("tns:buildIndex");
		group.patternSearch = groupTree.get<bool>("tns:patternSearch");
		group.naiveSearch = groupTree.get<bool>("tns:naiveSearch");
		group.loadIndexInMemory = groupTree.get<bool>("tns:loadIndexInMemory");
		group.keyIndex=groupTree.get<int>("tns:keyIndex");
		for( auto const& groupChildNode : groupTree.get_child("") ) {
			if (groupChildNode.first=="tns:patternSets"){
				PatternSet patternSet;
				boost::property_tree::ptree patternSetSubTree = groupChildNode.second;
				patternSet.length= patternSetSubTree.get<int>("tns:length");
				for( auto const& patternSetChildNode : patternSetSubTree.get_child("") ) {
					if (patternSetChildNode.first=="tns:patterns"){
						//string s=patternSetChildNode.first;
						string p=patternSetChildNode.second.data();
						patternSet.patterns.push_back(p);
					}
				}
				group.patternSets.push_back(patternSet);
			} else if (groupChildNode.first=="tns:superAlphabetOrders"){
				string s = groupChildNode.second.data();
				int64_t order=stoi(s);
				group.superAlphabetOrders.push_back(order);
			} else if (groupChildNode.first=="tns:bucketSizes"){
						string s = groupChildNode.second.data();
						int32_t bs=stoi(s);
						group.bucketSizes.push_back(bs);
			} else if (groupChildNode.first=="tns:markedRowsPercentages"){
				string s = groupChildNode.second.data();
				int8_t mrp=stoi(s);
				group.markedRowsPercentages.push_back(mrp);
			}
		}
		testGroups.push_back(group);
	}
}




} /* namespace std */
