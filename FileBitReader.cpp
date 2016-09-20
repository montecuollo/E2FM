/*
 * FileBitReader.cpp
 *
 *  Created on: 08/ott/2014
 *      Author: fernando
 */

#include "FileBitReader.h"
#include "BitReader.h"
#include "Bookmark.h"
#include <iostream>
#include <stdexcept>
namespace std {

FileBitReader::FileBitReader(string fileName){
		this->fileName=fileName;
}

FileBitReader::FileBitReader(){
}

FileBitReader::~FileBitReader() {
	// TODO Auto-generated destructor stub
}

void FileBitReader::open(){
	in.open(fileName, std::ifstream::binary);
	if (!in)
			throw runtime_error("Error in opening file " + fileName);
	if (blockSize==0)
	        	blockSize=512;
    blockBuffer=new char[blockSize];
    readBlocks=0;
    blockBufferEffectiveLength=0;
    liveBytes=0;

	BitReader::open();
}

void FileBitReader::open(string fileName){
	this->fileName=fileName;
	open();
}

void FileBitReader::close(){
	in.close();
	BitReader::close();
}

uint32_t FileBitReader::readBlockBuffer(){
	//clearBlockBuffer();
	in.read(blockBuffer,blockSize);
	blockBufferEffectiveLength=in.gcount();
	liveBytes=blockBufferEffectiveLength;
	readBlocks++;
	return liveBytes;
}


uint8_t FileBitReader::getByteFromBlockBuffer(bool *errorFlag){
	    	if (liveBytes == 0){
	    		try{
	    			liveBytes=readBlockBuffer();  //
	    		}
	    		catch(exception ex){
	    			*errorFlag=true;
	    			return 0;
	    		}
	    	}
	    	if (liveBytes < 0){
	    		return 0;
	    		*errorFlag=true;
	    	}
	    	else{
	    		uint8_t thech=blockBuffer[blockBufferEffectiveLength -liveBytes]&0xFF;
	    		liveBytes--;
	    		*errorFlag=false;
	    		return thech;
	    	}
}




 void FileBitReader::gotoBookmark(Bookmark *bookmark){
	liveBytes=0;
	liveBits=0;

	if (bookmark->getLiveBits()==0){
	  ifstream::pos_type p=bookmark->position;
	  in.clear();
	  in.seekg (p, in.beg);
	  //in.seekg (p);
	  bitsBuffer=bookmark->bitsBuffer;
	}
	else{
		//ifstream::pos_type p=bookmark->position-1;
		ifstream::pos_type p=bookmark->position;
		in.clear();
		in.seekg(p,in.beg);
		//in.seekg(p);
		read(8-bookmark->liveBits);
	}
}

Bookmark *FileBitReader::getBookmark(){
	ifstream::pos_type p=in.tellg();
   	return new Bookmark(p - liveBytes,bitsBuffer,liveBits);
}

void FileBitReader::gotoNextByteStart(){
    	if (liveBits!=0){
    		bitsBuffer=0;
    		liveBits=0;
    	}
}


int32_t FileBitReader::getInt(){
	int32_t acc=read(8);
	acc = (acc << 8)| read(8);
	acc = (acc << 8)| read(8);
	acc = (acc << 8)| read(8);
	return acc;
}

int64_t FileBitReader::getInt64() {
   	int64_t acc=read(8);
   	acc = (acc << 8)| read(8);
   	acc = (acc << 8)| read(8);
   	acc = (acc << 8)| read(8);
   	acc = (acc << 8)| read(8);
   	acc = (acc << 8)| read(8);
   	acc = (acc << 8)| read(8);
   	acc = (acc << 8)| read(8);

       return acc;
}


} /* namespace std */
