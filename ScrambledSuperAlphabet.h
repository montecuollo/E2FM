/*
 * ScrambledSuperAlphabet.h
 *
 *  Created on: 19/dic/2014
 *      Author: fernando
 */

#ifndef SCRAMBLEDSUPERALPHABET_H_
#define SCRAMBLEDSUPERALPHABET_H_
#include <string>
#include <vector>

namespace std {

class ScrambledSuperAlphabet {
public:
	ScrambledSuperAlphabet(vector<char> originalSymbols,uint8_t order,vector<uint32_t> *permutation);
	virtual ~ScrambledSuperAlphabet();

	string getSymbol(uint32_t superAlphabetCode);

	uint32_t getCode(string symbol);
	uint32_t getCode(string &text,uint64_t symbolStartIndex); //used during super-text building
	uint32_t getCode(uint32_t symbolNaturalCode);
	uint32_t getSize();
	uint8_t getOrder();

	string toString(uint32_t text[]);

	string toString(uint32_t text[],uint64_t from, uint64_t to);

	uint32_t getOriginalAlphabetCode(char c);
	uint32_t getOriginalAlphabetSize();
	vector<char> getOriginalSymbols();

	friend class EFMIndex;
private:
	uint32_t getSymbolNaturalCode(string symbol);
	uint32_t getSymbolNaturalCode(string &text,uint64_t symbolStartIndex);
	string getSymbolByNaturalCode(uint32_t naturalCode);

	uint32_t originalSymbolsIndexOf(char c);
	bool isOriginalSymbol(char c);
	uint32_t permutationIndexOf(uint32_t value);

	vector<char> originalSymbols;  	//symbols of the primigenious (order 1) alphabet
	uint32_t originalAlphabetSize;  //number of the original alphabet symbols
	uint8_t order;   				//order (1 for the original alphabet, >1 if this is a super-alphabet)
	uint32_t size;                  //cardinality (number of symbols) of the super alphabet
	string *inverseMap;				//maps the super-alphabet code to their relative k-mers
	uint32_t allNCode;					//scrambled code of the repeated symbol (in DNA sequences it's the allN symbol)
	string allNSymbol;				//string composed of a number of N equal to the superalphabet's order
	vector<uint32_t> *permutation;

};






} /* namespace std */

#endif /* SCRAMBLEDSUPERALPHABET_H_ */
