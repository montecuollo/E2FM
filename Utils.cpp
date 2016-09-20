/*
 * Utils.cpp
 *
 *  Created on: 24/nov/2014
 *      Author: fernando
 */

#include "Utils.h"
#include <string>
#include <iostream>
#include <fstream>
namespace std {

Utils::Utils() {
	// TODO Auto-generated constructor stub

}

Utils::~Utils() {
	// TODO Auto-generated destructor stub
}

/**
 * Loads in a buffer of size sequenceLength + k-1 the sequence contained in the fastaFileName,
 * so that the buffer contains all the characters to compute the initial k-mers of the last text suffixes
 * @param fastaFileName
 * @throws Exception
 */
string Utils::loadFasta(string fastaFileName) {
	std::ifstream in(fastaFileName, std::ios::in);
	if (in) {
		std::string sequence,line;
		std::getline( in, line );
		//while( std::getline( in, line ).good() ){
		while( std::getline( in, line )){
		        sequence += line;
		}
		return sequence;
	}
}

/*
int Utils::int_log2(int u)    // compute # bits to represent u
{
  int i = 1;
  int r = 1;

  while((i<=32) && (r<u)){
    r=2*r+1;
    i = i+1;
  }

  return i;
}

*/

uint8_t Utils::int_log2(uint64_t u)    // compute # bits to represent u
{
  uint8_t i = 1;
  uint64_t r = 1;

  while((i<=64) && (r<u)){
    r=2*r+1;
    i = i+1;
  }

  return i;
}



string Utils::getFileName(string filePath){
	// Remove directory if present.
	// Do this before extension removal incase directory has a period character.
	const size_t last_slash_idx = filePath.find_last_of("\\/");
	if (std::string::npos != last_slash_idx)
	{
	    filePath.erase(0, last_slash_idx + 1);
	}

	// Remove extension if present.
	const size_t period_idx = filePath.rfind('.');
	if (std::string::npos != period_idx)
	{
	    filePath.erase(period_idx);
	}
	return filePath;

}

} /* namespace std */
