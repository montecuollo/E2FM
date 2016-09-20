/*
 * MemoryBitReader.cpp
 *
 *  Created on: 19/nov/2015
 *      Author: fernando
 */

#include "MemoryBitReader.h"
#include <fstream>
#include <stdexcept>
#include <cmath>
namespace std {

MemoryBitReader::MemoryBitReader() {
	// TODO Auto-generated constructor stub

}

MemoryBitReader::MemoryBitReader(char *data,uint64_t size){
	blockBuffer=data;
	blockBufferEffectiveLength=size+7;
}

MemoryBitReader::~MemoryBitReader() {
	// TODO Auto-generated destructor stub
}

void MemoryBitReader::open(){
	BitReader::open();
	position=0;
	liveBytes=blockBufferEffectiveLength-7;
}

uint8_t MemoryBitReader::getByteFromBlockBuffer(bool *errorFlag){
	if (liveBytes == 0)
		*errorFlag=true;
	else{
		uint8_t thech=blockBuffer[position]&0xFF;
		position++;
		liveBytes--;
		*errorFlag=false;
		return thech;
	}
}


int32_t MemoryBitReader::getInt() {
    return read(32);
}

int64_t MemoryBitReader::getInt64() {
    return (read(56)<<8)|read(8);
}

/*
void MemoryBitReader::gotoBookmark(Bookmark *bookmark){
	liveBits=0;
	if (bookmark->liveBits==0){
	  position=bookmark->position;
	  liveBytes=blockBufferEffectiveLength-position;
	}
	else{
		uint8_t thech = blockBuffer[bookmark->position-1]&0xFF;
		//accoda i bit del carattere appena letto all'attuale contenuto del bit buffer
		bitsBuffer =  thech;
		liveBits=bookmark->liveBits;
		position=bookmark->position;
		liveBytes=blockBufferEffectiveLength-position;
	}
}*/

void MemoryBitReader::gotoBookmark(Bookmark *bookmark){
  position=bookmark->position;
  liveBytes=blockBufferEffectiveLength-position;
  liveBits=bookmark->liveBits;
}

void MemoryBitReader::gotoNextByteStart(){
	if (liveBits!=0){
		position++;
		liveBits=0;
	}
}

Bookmark *MemoryBitReader::getBookmark(){
  	return new Bookmark(position,bitsBuffer,liveBits);
}

void MemoryBitReader::close(){
	liveBytes=0;
}

char *MemoryBitReader::loadDataFromFile(string filePath,uint64_t &size){
	 char* data = NULL;
	 ifstream inFile;
	 inFile.open(filePath.c_str(), ios::in|ios::binary|ios::ate);
	 size = inFile.tellg();
	 inFile.seekg(0, ios::beg);
	 data = new char[7+size];
	 inFile.read(data, size);
	 for (uint64_t i=size;i<size+7;i++)
		 data[i]=0;

	 inFile.close();
	 return data;
}




union QuadWord{
	uint64_t value;
	uint8_t array[8];
};

union DoubleWord{
	uint32_t value;
	uint8_t array[4];
};

union Word{
	uint16_t value;
	uint8_t array[2];
};

int64_t MemoryBitReader::read(uint8_t n){

	char *pqwv=(blockBuffer+position); //pointer to quadword value
	uint8_t headBits;
	uint8_t tailBits;
	int64_t retVal;
	if (n<9){
		Word w;
		w.array[0]=pqwv[1];
		w.array[1]=pqwv[0];
		headBits=(8-liveBits)&7;
		tailBits=16-headBits-n;
		w.value=w.value >> tailBits;
		int16_t mask=(1u<<n)-1;
		retVal=w.value & mask;
	} else if (n<25){
		DoubleWord dw;
		dw.array[0]=pqwv[3];
		dw.array[1]=pqwv[2];
		dw.array[2]=pqwv[1];
		dw.array[3]=pqwv[0];
		headBits=(8-liveBits)&7;
		tailBits=32-headBits-n;
		dw.value=dw.value >> tailBits;
		int32_t mask=(1u<<n)-1;
		retVal=dw.value & mask;
	} else {
		QuadWord qw;
		qw.array[0]=pqwv[7];
		qw.array[1]=pqwv[6];
		qw.array[2]=pqwv[5];
		qw.array[3]=pqwv[4];
		qw.array[4]=pqwv[3];
		qw.array[5]=pqwv[2];
		qw.array[6]=pqwv[1];
		qw.array[7]=pqwv[0];
		headBits=(8-liveBits) & 7;  //=(8-liveBits)%8
		tailBits=64-headBits-n;
		qw.value=qw.value >> tailBits;
		uint64_t mask=(1ul<<n)-1;
		retVal=qw.value & mask;
	}

	//if (lbp > liveBytes-1)
	//	throw runtime_error("unexpected end of stream");
    liveBits=tailBits&7;  //=tailBits%8
    position+=(headBits+n)/8;
    return retVal;
}








} /* namespace std */
