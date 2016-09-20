/*
 * Utils.h
 *
 *  Created on: 24/nov/2014
 *      Author: fernando
 */

#ifndef UTILS_H_
#define UTILS_H_
#include <string>
namespace std {

class Utils {
public:
	Utils();
	virtual ~Utils();
	//static int int_log2(int u);
	static uint8_t int_log2(uint64_t u);  //TODO: testare che sia completamente equivalente alla precedente a 32 bit
	static string loadFasta(string fastaFileName);
	static string getFileName(string filePath);

};

} /* namespace std */

#endif /* UTILS_H_ */
