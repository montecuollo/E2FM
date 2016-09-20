/*
 * ScrambledSuperAlphabet.cpp
 *
 *  Created on: 19/dic/2014
 *      Author: fernando
 */

#include "ScrambledSuperAlphabet.h"
#include <algorithm>
#include <iostream>
#include <cmath>

namespace std {


ScrambledSuperAlphabet::ScrambledSuperAlphabet(vector<char> originalSymbols,uint8_t order,vector<uint32_t> *permutation) {
	this->order=order;
	this->originalSymbols=originalSymbols;
	this->originalAlphabetSize=originalSymbols.size();
	this->size = pow(originalAlphabetSize, order);
	this->permutation=permutation;
	this->allNSymbol=std::string(order, 'N');
	if (isOriginalSymbol('N')){
		allNCode=permutation->size()-1;
		uint32_t repeatedSymbolNaturalCode=getSymbolNaturalCode(allNSymbol);
		uint32_t repeatedSymbolPermutationIndex=repeatedSymbolNaturalCode;//+1;
		uint32_t repeatedSymbolCode=(*permutation)[repeatedSymbolPermutationIndex];
		uint32_t posOfRepeatedSymbolScrambledCode=permutationIndexOf(allNCode);
		(*permutation)[posOfRepeatedSymbolScrambledCode]=repeatedSymbolCode;
		(*permutation)[repeatedSymbolPermutationIndex]=allNCode;
	}
	uint32_t ps=permutation->size();
	uint32_t *inversePerm=new uint32_t[ps];
	for (uint32_t i=0;i<ps;i++){
		inversePerm[(*permutation)[i]]=i;
	}
	inverseMap=new string[permutation->size()];
	for (uint32_t ord=0;ord<ps;ord++){
		uint32_t symbolNaturalCode=inversePerm[ord];
		inverseMap[ord]=getSymbolByNaturalCode(symbolNaturalCode);
	}
	delete[] inversePerm;
}

ScrambledSuperAlphabet::~ScrambledSuperAlphabet() {
	delete[] inverseMap;
}

/**
* Returns the position of a certain value within the monoAlphabeticKey array
* @return
 */
inline uint32_t ScrambledSuperAlphabet::permutationIndexOf(uint32_t value){
	vector<uint32_t>::iterator it=find(permutation->begin(),permutation->end(),value);
	return it - permutation->begin();
}

inline uint32_t ScrambledSuperAlphabet::originalSymbolsIndexOf(char c) {
	vector<char>::iterator it=find(originalSymbols.begin(),originalSymbols.end(),c);
	return it - originalSymbols.begin();
}

inline bool ScrambledSuperAlphabet::isOriginalSymbol(char c) {
	vector<char>::iterator it=find(originalSymbols.begin(),originalSymbols.end(),c);
	return it != originalSymbols.end();
}


//returns the positional notation code of a given symbol
inline uint32_t ScrambledSuperAlphabet::getSymbolNaturalCode(string symbol) {
		uint32_t positionalNotationCode = 0;
		uint32_t b = originalAlphabetSize;
		for (uint8_t i = 0; i < order; i++)
			positionalNotationCode += originalSymbolsIndexOf(symbol[i]) * pow(b, order - i - 1);

		return positionalNotationCode;
}

//returns the positional notation code of a given symbol
inline uint32_t ScrambledSuperAlphabet::getSymbolNaturalCode(string &text,uint64_t symbolStartIndex) {
		uint32_t positionalNotationCode = 0;
		uint32_t b = originalAlphabetSize;
		for (uint8_t i = 0; i < order; i++)
			positionalNotationCode += originalSymbolsIndexOf(text[symbolStartIndex+i]) * pow(b, order - i - 1);

		return positionalNotationCode;
}

inline string ScrambledSuperAlphabet::getSymbolByNaturalCode(uint32_t naturalCode) {
	string retVal = "";
	if (naturalCode > 0) {
		uint32_t b = originalAlphabetSize;
		uint32_t r;
		uint32_t remainders[this->order];
		int64_t pos = -1;
		do {
			pos++;
			r = naturalCode % b;
			naturalCode = naturalCode / b;
			remainders[pos] = r;

		} while (naturalCode != 0);
		for (int64_t i = pos; i >= 0; i--)
			retVal += originalSymbols[remainders[i]];
	}
	uint32_t l = retVal.length();
	for (uint32_t i = l; i < this->order; i++)
		retVal = originalSymbols[0] + retVal;
	return retVal;
}



uint32_t ScrambledSuperAlphabet::getCode(uint32_t symbolNaturalCode){
	return (*permutation)[symbolNaturalCode];
}

string ScrambledSuperAlphabet::getSymbol(uint32_t superAlphabetCode){
	return inverseMap[superAlphabetCode];
}


uint32_t ScrambledSuperAlphabet::getSize(){
	return size;
}

uint8_t ScrambledSuperAlphabet::getOrder(){
	return order;
}

uint32_t ScrambledSuperAlphabet::getOriginalAlphabetSize(){
	return originalAlphabetSize;
}

vector<char> ScrambledSuperAlphabet::getOriginalSymbols(){
	return originalSymbols;
}


uint32_t ScrambledSuperAlphabet::getCode(string symbol){
	return (*permutation)[getSymbolNaturalCode(symbol)];
}

uint32_t ScrambledSuperAlphabet::getCode(string &text,uint64_t symbolStartIndex){
	return (*permutation)[getSymbolNaturalCode(text,symbolStartIndex)];
}

uint32_t  ScrambledSuperAlphabet::getOriginalAlphabetCode(char c){
	return originalSymbolsIndexOf(c);
}


} /* namespace std */
