/*
 * ArrayMTFList.cpp
 *
 *  Created on: 24/dic/2014
 *      Author: fernando
 */

#include "ArrayMTFList.h"
#include <algorithm>
#include <cstring>
#include <iostream>

namespace std {

ArrayMTFList::ArrayMTFList(uint32_t alphabetSize){
	this->alphabetSize=alphabetSize;
	elementData=new uint32_t[alphabetSize];
	for (uint32_t i=0;i<alphabetSize;i++)
		elementData[i]=i;
}

ArrayMTFList::~ArrayMTFList() {
	delete[] elementData;
}


/**
 * Moves to the front of the list the item actually at itemPosition
 * @param itemPosition position of the item to move to the front of the list
 */

uint32_t ArrayMTFList::moveToFront(uint32_t position){
	uint32_t elementToMove=elementData[position];
	//int64_t soi=sizeof(uint32_t); it's 4
	if (position >0){
		memmove(elementData+1, elementData, 4 * position);
		elementData[0]=elementToMove;
	}
	return elementToMove;
}

int64_t ArrayMTFList::indexOf(uint32_t item){
	for (uint32_t i = 0; i < alphabetSize; i++)
		if (elementData[i]==item)
			return i;
	return -1;
}


uint32_t ArrayMTFList::get(uint32_t position){
	return elementData[position];
}


void ArrayMTFList::dumpList(){
	cout << "[" ;
	for (uint32_t i=0;i<alphabetSize;i++){
		cout <<  elementData[i];
		if (i<alphabetSize -1)
			cout <<",";
	}
	cout << "]" << endl ;
}





} /* namespace std */
