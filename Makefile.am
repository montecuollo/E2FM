# what flags you want to pass to the C compiler & linker
AM_CFLAGS=-ftrapv -march=native -fmessage-length=0 -pthread -Wchar-subscripts -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free 
AM_LDFLAGS=-ftrapv -march=native -pthread -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free 
AM_CPPFLAGS=-I/usr/local/boost_1_58_0 -I../../src/seqan-include
AM_CXXFLAGS=-ftrapv -march=native -fmessage-length=0 -pthread -std=c++11 -Wchar-subscripts -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS=efmc 
efmc_SOURCES=salsa20.S ecrypt-portable.h ecrypt-config.h ecrypt-machine.h ecrypt-sync.h Main.cpp ArrayMTFList.cpp ArrayMTFList.h BalancedTreeMTFList.cpp BalancedTreeMTFList.h BitReader.cpp BitReader.h BitSetRLEEncoder.cpp BitSetRLEEncoder.h BitWriter.cpp BitWriter.h Bookmark.cpp Bookmark.h Bucket.cpp Bucket.h Common.h EFMCollection.cpp EFMCollection.h EFMIndex.cpp EFMIndex.h EncryptionManager.cpp EncryptionManager.h FastBWTransformer.cpp FastBWTransformer.h MemoryBitReader.cpp MemoryBitReader.h FileBitReader.cpp FileBitReader.h FileBitWriter.cpp FileBitWriter.h FMCollection.cpp FMCollection.h RandomGenerator.cpp RandomGenerator.h ScrambledSuperAlphabet.cpp ScrambledSuperAlphabet.h SFMCollection.cpp SFMCollection.h Statistics.cpp Statistics.h SuperBucket.cpp SuperBucket.h Test.cpp Test.h tinydir.h Utils.cpp Utils.h 
