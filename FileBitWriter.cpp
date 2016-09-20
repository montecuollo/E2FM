/*
 * FileBitWriter.cpp
 *
 *  Created on: 08/ott/2014
 *      Author: fernando
 */

#include "FileBitWriter.h"
#include "Bookmark.h"
#include <fstream>

namespace std {

FileBitWriter::FileBitWriter(string fileName){
		this->fileName=fileName;
}


FileBitWriter::~FileBitWriter() {
	// TODO Auto-generated destructor stub
}


Bookmark *FileBitWriter::getBookmark(){
	uint64_t p=out.tellp();
   	return new Bookmark(p+liveBytes,bitsBuffer,liveBits);
}


/**
 * Goto bookmark, restoring file position and bit buffer state
 * @param bookmark
 * @throws IOException
 */

void FileBitWriter::gotoBookmark(Bookmark *bookmark){
	finishedWithStream();
	out.seekp(bookmark->getPosition());
	bitsBuffer=bookmark->getBitsBuffer();
	liveBits=bookmark->getLiveBits();
	clearBlockBuffer();
	liveBytes=0;
}


void FileBitWriter::open(){
	out.open(fileName,ios::out | ios::binary | ios::trunc);
	init();
}

void FileBitWriter::close(){
	finishedWithStream();
	out.close();
}



void FileBitWriter::writeBlockBuffer(){
	out.write(blockBuffer, liveBytes);
}





} /* namespace std */
