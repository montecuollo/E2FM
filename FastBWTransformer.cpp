/*
 * FastBWTransformer.cpp
 *
 *  Created on: 22/dic/2014
 *      Author: fernando
 */

#include "FastBWTransformer.h"
#include <iostream>
#include "Common.h"

namespace std {

std::mutex coutMutex;

typedef std::chrono::high_resolution_clock myclock;

#ifndef min2
#define min2(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

#define swap(x,a, b) {uint64_t t=x[a]; x[a]=x[b]; x[b]=t;}

#define rotChar(x,n,ri,p)({ ri + p < n ? x[ri + p]:x[ri + p -n];})

/**
 * Compares two rotations identified by rotationIndex1 and rotationIndex2,
 * starting from character depth; code is mutuated from String.compareTo
 * (as implemented in Open JDK 7)
 * @param text whole text
 * @param n whole text length
 * @param rotationIndex1
 * @param rotationIndex2
 * @param depth
 * @return
 */

//text length



int32_t compareRotations(uint32_t *text,uint64_t n,uint64_t rotationIndex1,uint64_t rotationIndex2,uint64_t depth){
	  uint64_t i=rotationIndex1+depth;
	  if (i >= n)
			   i = i - n;
	  uint64_t j=rotationIndex2+depth;
	  if (j >= n)
		   j = j - n;
	  int64_t m=n;
	  while (m-- != 0) {
			   if ( text[i] !=  text[j])
					 return text[i] - text[j];
			   else{
				   if (i<n-1)
					   i=i+1;
				   else
					   i=0;

				   if (j<n-1)
					   j=j+1;
				   else
					   j=0;
			   }
	  }
	  return 0;
}


/**
	 * Insertion sort of the subarray of rotations
	 * The algorithm is very easy to understand: it simulates at each step the insertion of a new element into
	 * an ordered array, inserting it at the end of the array and, starting from the end, swapping the new element
	 * with all the preceding elements a[j] until a[j]<=newElement (each step terminates where the new element has been
	 * put in the right place so that the extended array is ordered)
	 * @param sub-array to sort
	 * @param ns length of the sub-array
	 * @param depth	position of the rotations from which is necessary
	 *              to start the comparison of symbols
	 * @throws Exception
	 */
	void inssort(uint32_t *text,uint64_t n,uint64_t *a, uint64_t ns, uint64_t depth)
	{   uint64_t i, j;
		for (i = 1; i < ns; i++)
		  for (j = i; j > 0; j--) {
			 if (compareRotations(text,n,a[j-1], a[j],depth) <= 0)
				 break;
			 swap(a,j-1,j)
		  }
	}



/**
 * Sorts a sub-array of rotations: each subarray element is an integer value representing
 * the starting index of a certain rotation.
 * @param text	whole text
 * @param n	whole text length
 * @param randomGenerator random generator used to swap the element in position 0 with a random element
 * @param sa	sub-array to sort
 * @param ns    length of the sub-array to sort
 * @param depth	position of the rotations from which is necessary
 *              to start the comparison of symbols
 * @param rl	recursion level
 * @throws Exception
 */
void ssort(uint32_t *text,uint64_t n,minstd_rand0* randomGenerator,uint64_t* sa, uint64_t ns,uint64_t depth,uint64_t rl){

	  int64_t a, b, c, d, r,v;
	  if (ns <=  50) {
			inssort(text,n,sa, ns, depth);
			return;
	  }

	  //Scelta dell'elemento Pivot e swap di esso con l'elemento di posto 0
	  int64_t toSwapElement = (*randomGenerator)() % ns;
	  swap(sa,0,toSwapElement);



	  v = rotChar(text,n,sa[0],depth);

	  a = b = 1;
	  c = d = ns-1;

	  for (;;) {
		  while (b <= c && (r=rotChar(text,n,sa[b],depth) - v)<=0){
			 if (r==0){  if (a!=b){swap(sa,a, b);}
						 a++;
					  }
			 b++;
		  }
		  while(b <= c && (r=rotChar(text,n,sa[c],depth) -v)>=0){
			 if (r == 0){ if (c!=d){swap(sa,c, d);}
						 d--;
						}
			 c--;
		  }

		  if (b > c) break;
		  swap(sa, b, c);
		  b++;
		  c--;

	 }
	 if (a<b-a)
		 r=a;
	 else
		 r=b-a;


	 //vecswap(sa, 0, b-r, r);
	 uint64_t conta=r;
	 uint64_t* v1=sa;
	 uint64_t* v2=sa + b-r;
	 while (conta-- > 0){
		int64_t t = *v1;
		*v1=*v2;
		*v2=t;
		v1++;
		v2++;
	 }



	 //r = min2(d-c, ns-d-1);
	 if (d-c < ns-d-1)
		 r=d-c;
	 else
		 r=ns-d-1;
	 //vecswap(sa, b, ns-r, r);
	 conta=r;
	 v1=sa+b;
	 v2=sa + ns -r;
	 while (conta-- > 0){
		int64_t t = *v1;
		*v1=*v2;
		*v2=t;
		v1++;
		v2++;
	 }

	 r=b-a;
	 if (r>0)
		 ssort(text,n,randomGenerator,sa, r, depth,rl+1);  //orders left part

	 if (depth<n-1){  //Order central part in base of components from depth+1 to n
		 if (a+ns-d-1>0)
			 ssort(text,n,randomGenerator,sa+r, a + ns-d-1, depth+1,rl+1);
	 }
	 r=d-c;
	 if (r>0)
		 ssort(text,n,randomGenerator,sa+ns-r,r,depth,rl+1); //orders right part
	 return;
}


typedef pair<uint64_t,uint64_t> UIntIntPair;







ostream& operator<<(ostream& os, const ArrayInterval& i) {
	os << "[" << i.left << "," << i.right << "] (" << i.weight << ")";
	return os;
}

ostream& operator<<(ostream& os, const AlphabetRange& r) {
		os << "[" << r.firstCharacter << "," << r.lastCharacter << "] (" << r.numberOfSuffixes << " suffixes)";
		return os;
}

ostream& operator<<(ostream& os, const Repetition& r) {
	os << "[" << r.start << "," << r.length << "]";
	return os;
}

ostream& operator<<(ostream& os, const RepetitionInformations& r) {
	os << "mrl=" << r.maximumRepetitionLength << " {";
	for (uint64_t i = 0; i < r.repetitions.size(); i++) {
		os << r.repetitions[i];
		if (i < r.repetitions.size() - 1)
			os << ",";
	}
	os << "}";
	return os;
}






Distributor::Distributor(uint64_t firstTextPosition,uint64_t lastTextPosition,uint32_t *text,ScrambledSuperAlphabet *superAlphabet,
		uint32_t er,uint32_t elementaryRangesLength,
		AlphabetRange **elementaryRanges) {
			this->firstTextPosition = firstTextPosition;
			this->lastTextPosition = lastTextPosition;
			this->text=text;
			this->superAlphabet=superAlphabet;
			this->er=er;
			this->elementaryRangesLength=elementaryRangesLength;
			this->elementaryRanges=elementaryRanges;
}


void Distributor::run() {
	distribute();
}

void Distributor::distribute() {

	//create a map containing, for each possible repeated supercharacter, the related repetition informations
	uint32_t b = superAlphabet->getOriginalAlphabetSize();
	uint8_t k = superAlphabet->getOrder();
	vector<char>  originalSymbols= superAlphabet->getOriginalSymbols();
	for (uint32_t i = 1; i < b; i++){  //the first original character is always '$'
		string repeatedSuperCharacter(k, originalSymbols[i]);
		uint32_t repeatedSuperCharacterCode=superAlphabet->getCode(repeatedSuperCharacter);
		repetitions.insert(
				SymbolRepInfoPair(repeatedSuperCharacterCode,
						RepetitionInformations()));
	}


	int64_t pc = -1;  //previous character
	int64_t actualRepetitionLength = 1;

	for (int64_t i = firstTextPosition; i <= lastTextPosition; i++) {

		uint32_t c = text[i];
		uint32_t r = c / elementaryRangesLength;
		if (r > er - 1)
			r = er - 1;
		elementaryRanges[r]->addSuffix(i);


		//Handle the repetitions
		//Store only the repetitions that are long at least MINIMUM_REPETITION_LENGTH characters
		//and the repetitions starting from the first character (they could be the remaining part of a repetition started
		//in a text position handled by preceding threads)
		if (c != pc || i == lastTextPosition) {
			if (c != pc) {
				//Registra le ripetizioni di previousCharacter lunghe almeno MINIMUM_REPETITIONS_LENGTH;
				//Se il thread non sta gestendo il primo range di posizioni del testo registra anche quelle che iniziano sul primo carattere,
				//indipendentemente dalla loro lunghezza, perchè esse potrebbero costituire il prolungamento di una ripetizione iniziata nel
				//range precedente
				if (actualRepetitionLength >= MINIMUM_REPETITION_LENGTH || (firstTextPosition > 0 && i - actualRepetitionLength== firstTextPosition)) {
					//le ripetizioni che iniziano sul primo carattere non devono contribuire all'aggiornamento della repInfo.maximumRepetitionLength
					//se non sono di lunghezza sufficiente
					if (actualRepetitionLength>= MINIMUM_REPETITION_LENGTH && actualRepetitionLength > repetitions[pc].maximumRepetitionLength)
						repetitions[pc].maximumRepetitionLength = actualRepetitionLength;
					repetitions[pc].repetitions.push_back(Repetition(i - actualRepetitionLength,actualRepetitionLength));
				}
				actualRepetitionLength = 1;
				//registra la potenziale "ripetizione", per ora di lunghezza 1, che inizia sull'ultimo carattere, perchè essa potrebbe continuare
				//nel range successivo, non gestito da questo thread
				if (i == lastTextPosition)
							this->repetitions[c].repetitions.push_back(Repetition(i, actualRepetitionLength));

			} else {
				//Registra, indipendentemente dalla sua lunghezza, la ripetizione che termina sull'ultimo carattere del range, perchè essa potrebbe continuare
				//nel range successivo, non gestito da questo thread
				actualRepetitionLength++;
				//RepetitionInformations* repInfo = &this->repetitions[c];
				if (actualRepetitionLength >= MINIMUM_REPETITION_LENGTH && actualRepetitionLength> repetitions[c].maximumRepetitionLength)
					repetitions[c].maximumRepetitionLength =
							actualRepetitionLength;
				repetitions[c].repetitions.push_back(Repetition(i + 1 - actualRepetitionLength,actualRepetitionLength));
			}
		} else
			actualRepetitionLength++;

		pc = c;
	}

}


NoLongRepetitionsSorter::NoLongRepetitionsSorter(uint32_t threadNumber,AlphabetRange **ranges,uint32_t firstRange,
		uint32_t lastRange,uint32_t *text,uint64_t n,ScrambledSuperAlphabet* superAlphabet,
		unordered_map<uint32_t, RepetitionInformations> repetitions,uint32_t numberOfThreads) {
				this->ranges=ranges;
				this->firstRange = firstRange;
				this->lastRange = lastRange;
				this->text=text;
				this->n=n;
				this->superAlphabet=superAlphabet;
				this->repetitions=repetitions;
				this->numberOfThreads=numberOfThreads;

				//Initialize the random numbers generation engine
				myclock::time_point beginning = myclock::now();
				myclock::duration d = myclock::now() - beginning;
				unsigned seed = d.count();
				this->randomGenerator=new minstd_rand0(seed);
				this->threadNumber=threadNumber;
}

NoLongRepetitionsSorter::~NoLongRepetitionsSorter(){
	delete randomGenerator;
}

void NoLongRepetitionsSorter::sort(){
#if INDEX_DEBUG_LEVEL > 0
coutMutex.lock();
cout << "\tSorting thread " << threadNumber << ": starting sort" << endl;
cout << "\t\t\t Ranges: " << firstRange << "-" << lastRange << endl;
coutMutex.unlock();
#endif
for (uint32_t r=firstRange;r<=lastRange;r++){
	AlphabetRange* range=ranges[r];
	//double rangeEffectiveSortTime;
	if (!range->containsLongRepetitionsSuffixes && range->numberOfSuffixes>1){
			#if INDEX_DEBUG_LEVEL > 0
			coutMutex.lock();
			cout << "\t\tSorting thread " << threadNumber << ": sorting range " << r <<endl;
			coutMutex.unlock();
			#endif
			//it's not a range containing a long repetition character and so it can be orderered as a whole
			ssort(text,n,randomGenerator,range->suffixes,range->numberOfSuffixes,0,1);
	}
}
#if INDEX_DEBUG_LEVEL > 0
coutMutex.lock();
cout << "\tSorting thread " << threadNumber << ": sort finished" << endl;
coutMutex.unlock();
#endif
}


void NoLongRepetitionsSorter::run(){
	sort();
}


LongRepetitionsSorter::LongRepetitionsSorter(uint32_t threadNumber,uint32_t *text,uint64_t n) {
				this->text=text;
				this->n=n;
				this->threadNumber=threadNumber;
				//Initialize the random numbers generation engine
				myclock::time_point beginning = myclock::now();
				myclock::duration d = myclock::now() - beginning;
				unsigned seed = d.count();
				this->randomGenerator=new minstd_rand0(seed);
}

LongRepetitionsSorter::~LongRepetitionsSorter(){
	delete randomGenerator;
}



void LongRepetitionsSorter::sortSubRanges(AlphabetRange *range,uint32_t firstSubRange,uint32_t lastSubRange,uint32_t longRepetitionsSubRangesCount){
	#if INDEX_DEBUG_LEVEL > 0
	coutMutex.lock();
	cout << "\t\tSorting thread " << threadNumber << ": starting sort" << endl;
	cout << "\t\t\t\t sub-Ranges: " << firstSubRange << "-" << lastSubRange << endl;
	coutMutex.unlock();
	#endif
	for (uint32_t sri=firstSubRange;sri<=lastSubRange;sri++){
		#if INDEX_DEBUG_LEVEL > 0
		coutMutex.lock();
		cout << "\t\t\tSorting thread " << threadNumber << ": sorting sub-range " << sri <<endl;
		coutMutex.unlock();
		#endif

		LongRepetitionsRange* subRange=&range->longRepetitionsSubRanges[sri];
		ssort(text,n,randomGenerator,subRange->suffixes,subRange->numberOfSuffixes,sri*SUBINTERVAL_LENGTH+1,1);

		subRange->insertionPoint=-1;
		//if this sub-range is not the latter then compute the insertion point
		if (sri<longRepetitionsSubRangesCount-1){
			uint64_t i=0;
			while (i<subRange->numberOfSuffixes && subRange->insertionPoint<0){
				if (subRange->firstNonRepeatedGreaterThanRepeated[subRange->suffixes[i]])
					subRange->insertionPoint=i;
				else
					i++;
			}
			//Se non è stato trovato alcun suffisso maggiore di tutti quelli del subrange successivo,
			//allora i suffissi del subrange successivo vanno riportati successivamente a quelli del subrange corrente
			if (subRange->insertionPoint==-1)
				subRange->insertionPoint=subRange->numberOfSuffixes;
		}
	}
	#if INDEX_DEBUG_LEVEL > 0
	coutMutex.lock();
	cout << "\t\tSorting thread " << threadNumber << ": sort finished" << endl;
	coutMutex.unlock();
	#endif
}






void FastBWTransformer::split() {
		std::chrono::time_point<std::chrono::system_clock> currentTime;
		std::time_t printableTime;
		currentTime=std::chrono::system_clock::now();
		printableTime= std::chrono::system_clock::to_time_t(currentTime);
		char mbstr[100];
		std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
		cout << "\tSplitting suffixes into " << r << " ranges intervals,  starting at " << mbstr <<endl;


		//allocates all the elementary ranges, linking each of them to a range of the superalphabet's characters
		elementaryRanges = new AlphabetRange*[er];
		nonEmptyRangesNumber = 0;
		//ranges = new AlphabetRange[er];
		uint32_t superAlphabetSize = superAlphabet->getSize();
		vector<char> alphabet=superAlphabet->getOriginalSymbols();

		suffixesMaximumNumberPerRange =  ceil(
				((double) (n) / er) * (1 + risePercentage));
		elementaryRangesLength = superAlphabetSize / er;
		uint32_t rn = 0;
		for (uint32_t i = 0; i < superAlphabetSize && rn < er;i = i + elementaryRangesLength) {
			elementaryRanges[rn]=new AlphabetRange();
			elementaryRanges[rn]->firstCharacter = i;
			elementaryRanges[rn]->lastCharacter = i + elementaryRangesLength - 1;
			rn++;
		}
		//the last elementary range must correspond to a superalphabet's range ending with the last superalphabet's character
		//(the following instruction is needed because the above code doesn't guarantee this property)
		elementaryRanges[er - 1]->lastCharacter = superAlphabetSize - 1;

		//distributes the suffixes among the ranges
		std::chrono::time_point<std::chrono::system_clock> startTime, endTime;
		startTime = std::chrono::system_clock::now();

		//split suffixes into ranges (multithreaded version) and find the characters repetitions
		uint32_t nt = r; //Allocate as many threads as the number of the microprocessor's cores to use
		uint64_t perThreadSuffixes = n / nt;
		std::vector<std::thread*> threads;
		std::vector<Distributor*> distributors;
		for (uint32_t t = 0; t < nt; t++) {
			uint64_t firstTextPosition = t * perThreadSuffixes;
			uint64_t lastTextPosition = (t + 1) * perThreadSuffixes - 1;
			if (t == nt - 1)
				lastTextPosition = n - 1;
			Distributor* distributor=new Distributor(firstTextPosition, lastTextPosition,
					text, superAlphabet,er,elementaryRangesLength,elementaryRanges);

			distributors.push_back(distributor);
			threads.push_back(new std::thread(&Distributor::run, distributor));
		}
		for (uint32_t i = 0; i < nt; i++){
			threads[i]->join();
		}

		endTime = std::chrono::system_clock::now();
	    distributionTime = std::chrono::duration_cast<std::chrono::milliseconds>(
									endTime - startTime).count();

		//merge the repetitions
		//create a map containing, for each possibly repeated super-character, the related repetition informations
		uint32_t b = superAlphabet->getOriginalAlphabetSize();
		uint8_t k = superAlphabet->getOrder();
		vector<char>  originalSymbols= superAlphabet->getOriginalSymbols();
		for (uint32_t i = 1; i < b; i++){  //the first original character is always '$' and cannot be repeated
			string repeatedSuperCharacter(k, originalSymbols[i]);
			uint32_t repeatedSuperCharacterCode=superAlphabet->getCode(repeatedSuperCharacter);
			repetitions.insert(
					SymbolRepInfoPair(repeatedSuperCharacterCode,
							RepetitionInformations()));
		}

		for (auto it=repetitions.begin(); it!=repetitions.end(); ++it){
			RepetitionInformations *globalRepInfo = &it->second;
			for (uint32_t t = 0; t < nt; t++) {
				RepetitionInformations threadRepInfo =
						distributors[t]->repetitions[it->first];

				int64_t lastGlobalRepetitionIndex = (globalRepInfo->repetitions.size() > 0 ? globalRepInfo->repetitions.size() - 1 : -1);
				if (threadRepInfo.repetitions.size() > 0) {
					Repetition rep = threadRepInfo.repetitions[0];
					//verifica se è necessario effettuare il merge
					if (rep.start == distributors[t]->firstTextPosition && lastGlobalRepetitionIndex > -1
							&& globalRepInfo->repetitions[lastGlobalRepetitionIndex].start + globalRepInfo->repetitions[lastGlobalRepetitionIndex].length == rep.start) {
						uint64_t mergedRepetitionLength = globalRepInfo->repetitions[lastGlobalRepetitionIndex].length + rep.length;
						globalRepInfo->repetitions[lastGlobalRepetitionIndex].length = mergedRepetitionLength;
						if (mergedRepetitionLength >= MINIMUM_REPETITION_LENGTH) {
							if (mergedRepetitionLength > globalRepInfo->maximumRepetitionLength)
								globalRepInfo->maximumRepetitionLength =mergedRepetitionLength;
						}
					} else { //se non è necessario effettuare il merge
							 //indipendentemente dalla loro lunghezza, anche le ripetizioni che terminano alla fine dell'intervallo vanno aggiunte
						if (rep.length >= MINIMUM_REPETITION_LENGTH || (rep.start + rep.length - 1 == distributors[t]->lastTextPosition)) {
							globalRepInfo->repetitions.push_back(rep);
							if (rep.length >= MINIMUM_REPETITION_LENGTH) {
								if (rep.length > globalRepInfo->maximumRepetitionLength)
									globalRepInfo->maximumRepetitionLength = rep.length;
							}
						}
					}
					//se la lastGlobalRepetition, eventualmente unita alla prima del thread sotto esame, non supera la lunghezza minima,
					//è necessario eliminarla; l'eliminazione viene fatta solo se la ripetizione non si estende fino alla fino dell'intervallo di testo
					//esaminato dal thread, perchè in caso contrario la decisione va effettuata al passo successivo
					if (lastGlobalRepetitionIndex > -1 &&
							globalRepInfo->repetitions[lastGlobalRepetitionIndex].start+globalRepInfo->repetitions[lastGlobalRepetitionIndex].length < distributors[t]->lastTextPosition &&
							globalRepInfo->repetitions[lastGlobalRepetitionIndex].length < MINIMUM_REPETITION_LENGTH)
						globalRepInfo->repetitions.erase(globalRepInfo->repetitions.begin() + lastGlobalRepetitionIndex);

					//aggiunge le ripetizioni successive alla prima
					for (uint64_t j = 1;
							j < threadRepInfo.repetitions.size(); j++) {
						rep = threadRepInfo.repetitions[j];
						globalRepInfo->repetitions.push_back(rep);
						if (rep.length >= MINIMUM_REPETITION_LENGTH && rep.length > globalRepInfo->maximumRepetitionLength) //l'ultima ripetizione del thread potrebbe avere lunghezza inferiore a MINIMUM_REPETITION_LENGTH
																				 //ed in tal caso verrà eliminata al passo successivo, per cui non deve contribuire all'aggiornamento della
																				 //globalRepInfo.maximumRepetitionLength
							globalRepInfo->maximumRepetitionLength = rep.length;
					}

				}
			}

			int64_t lastGlobalRepetitionIndex = (globalRepInfo->repetitions.size() > 0 ?globalRepInfo->repetitions.size() - 1 : -1);
			Repetition *lastGlobalRepetition;
			if (lastGlobalRepetitionIndex > -1)
				lastGlobalRepetition = &globalRepInfo->repetitions[lastGlobalRepetitionIndex];

			if (globalRepInfo->repetitions.size() > 1) {
				int64_t firstGlobalRepetitionIndex = (globalRepInfo->repetitions.size() > 0 ? 0 : -1);
				Repetition *firstGlobalRepetition;
				if (firstGlobalRepetitionIndex > -1)
					firstGlobalRepetition =&globalRepInfo->repetitions[firstGlobalRepetitionIndex];

				//Siccome si andranno ad ordinare le rotazioni del testo, l'ultima ripetizione va fusa con la prima se
				//l'ultima termina alla fine del testo e la prima inizia dalla posizione 0; per far ciò si aggiunge semplicemente alla lunghezza dell'ultima
				//quella della prima.
				if (lastGlobalRepetitionIndex > -1
						&& firstGlobalRepetitionIndex > -1
						&& lastGlobalRepetition->start + lastGlobalRepetition->length == n
						&& firstGlobalRepetition->start == 0)
					lastGlobalRepetition->length = lastGlobalRepetition->length + firstGlobalRepetition->length;
			}
			//se la lastGlobalRepetition non supera la lunghezza minima è necessario eliminarla
			if (lastGlobalRepetitionIndex > -1 && lastGlobalRepetition->length < MINIMUM_REPETITION_LENGTH)
				globalRepInfo->repetitions.erase(
						globalRepInfo->repetitions.begin()
								+ lastGlobalRepetitionIndex);

			repetitions.insert(SymbolRepInfoPair(it->first, *globalRepInfo));

		}

		//Leave only the really repeated characters in the repetitions map
		for (auto it=repetitions.begin(); it!=repetitions.end();){
			if (it->second.maximumRepetitionLength==0)
				repetitions.erase(it++);
			else
				++it;
		}

		//cout << "Global repetitions" <<endl;
		//for (auto it=repetitions.begin(); it!=repetitions.end(); ++it)
		//	cout<< superAlphabet->getSymbol(it->first) <<":"<< it->second << endl;

		for (uint32_t i = 0; i < nt; i++)
			delete distributors[i];
		for (uint32_t i = 0; i < nt; i++)
			delete threads[i];


		for (uint32_t i = 0; i < er; i++)
			if (elementaryRanges[i]->numberOfSuffixes > 0)
				nonEmptyRangesNumber++;

		int64_t array[nonEmptyRangesNumber];
		ranges = new AlphabetRange*[nonEmptyRangesNumber];
		rn = 0;
		for (uint32_t i = 0; i < er; i++) {
			if (elementaryRanges[i]->numberOfSuffixes > 0) {
				ranges[rn] = elementaryRanges[i];
				array[rn] = elementaryRanges[i]->numberOfSuffixes;
				rn++;
			}
		}

		ArraySplitter *arraySplitter=new ArraySplitter();
		arraySplitter->setArray(array, nonEmptyRangesNumber);
		arraySplitter->setMaxSteps(100);
		arraySplitter->setMinPercentualGain(0.001);
		arraySplitter->setNumberOfIntervals(r);
		arraySplitter->split();
		meanRangesRelativeError = 0;
		maxRangesRelativeError = 0;
		optimalWeight = 0;

		endTime = std::chrono::system_clock::now();
		splittingTime = std::chrono::duration_cast<std::chrono::milliseconds>(
				endTime - startTime).count();
		splittingSteps = arraySplitter->getDoneSteps();
		intervals = arraySplitter->getIntervals();

		delete arraySplitter;

	}


void FastBWTransformer::markLongRepetitionsRanges(){
	std::chrono::time_point<std::chrono::system_clock> currentTime;
	std::time_t printableTime;
	currentTime=std::chrono::system_clock::now();
	printableTime= std::chrono::system_clock::to_time_t(currentTime);
	char mbstr[100];
	std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
	cout << "\tMarking of long repetitions ranges, starting at " << mbstr << endl;

	//The long repetitions have to be managed only if at least a repeated superCharacter (k times a repeated character)
	//is contained in this range: in this case the variable handleLongRepetition will take the value true.
	//In the following code will be assumed that at most a single repeated character can exist in a range
	//(This assumption is correct if the chosen number of elementary ranges is at least equal to the number
	//of the superalphabet's characters: this constraint must be enforced (TODO))
	for (uint32_t rn=0;rn<nonEmptyRangesNumber;rn++){
		AlphabetRange* range=ranges[rn];
		int64_t repeatedSuperCharacterCode=-1;
		auto it = (repetitions).begin();
		while ( it != (repetitions).end() && !range->containsLongRepetitionsSuffixes ){
			repeatedSuperCharacterCode=it->first;
			if (repeatedSuperCharacterCode>=range->firstCharacter && repeatedSuperCharacterCode<=range->lastCharacter){
				range->containsLongRepetitionsSuffixes=true;
				range->longRepeatedCharacter=repeatedSuperCharacterCode;
			}
			it++;
		}
	}

}


void FastBWTransformer::mergeSortResultsAndComputeBWT(){
		std::chrono::time_point<std::chrono::system_clock> startTime,endTime;
		startTime=std::chrono::system_clock::now();
		std::time_t printableTime= std::chrono::system_clock::to_time_t(startTime);
		char mbstr[100];
		std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
		cout << "\tMerging sort results and computing BWT, starting at " << mbstr <<endl;

		bwt=new uint32_t[n];
		markedRows=new map<uint64_t,uint64_t>();
		sp=0;
		for (uint32_t rn=0;rn<nonEmptyRangesNumber;rn++){
			AlphabetRange* range=ranges[rn];
			if (range->containsLongRepetitionsSuffixes){
				//The range has been ordered as made of three separate parts:
				//the precedingSubRange, the set of the subranges containing long repetitions and
				//the following subrange
				//The logic below is repeated to achieve better performance

				//PRECEDING SUBRANGE
				//int[] suffixes=range->precedingSubRange.suffixes;
				int64_t sn=range->precedingSubRange.numberOfSuffixes;
				for (int64_t i=0;i<sn;i++){
					uint64_t suffixStartingPosition=range->precedingSubRange.suffixes[i];
					if (markingRate > 0 && suffixStartingPosition%markingRate==0)
						markedRows->insert(UIntIntPair(sp,suffixStartingPosition/markingRate));
					if (suffixStartingPosition==0)
						first=sp;
					uint64_t tp;
					if (suffixStartingPosition>0)
						tp = suffixStartingPosition-1;
					else
						tp = suffixStartingPosition + n -1;  //the super-text is seen as a circular tape
					bwt[sp]=text[tp];
					sp++;
				}

				//process the SUBRANGES containing the long repetitions suffixes
				//this method is called directly on the range 0 and will be recursively called for the other subranges

				//processLongRepetitionsSubRange(range,0);
				//ITERATIVE UNROLLMENT OF THE RECURSIVE CALL processLongRepetitionsSubRange(range,0);
				int64_t srn=range->numberOfLongRepetitionsSubranges;
				for (int64_t sri=0;sri<srn-1;sri++){
					LongRepetitionsRange* subRange=&range->longRepetitionsSubRanges[sri];
					//aggiunge i suffissi precedenti al punto di inserzione
					sn=subRange->numberOfSuffixes;
					for (int64_t i=0;i<subRange->insertionPoint;i++){
						uint64_t suffixStartingPosition=subRange->suffixes[i];
						if (markingRate > 0 && suffixStartingPosition%markingRate==0)
							 markedRows->insert(UIntIntPair(sp,suffixStartingPosition/markingRate));
						if (suffixStartingPosition==0)
							first=sp;
						uint64_t tp;
						if (suffixStartingPosition>0)
							tp = suffixStartingPosition-1;
						else
							tp = n -1;  //the super-text is seen as a circular tape
						bwt[sp]=text[tp];
						sp++;
					}
				}

				LongRepetitionsRange* subRange=&range->longRepetitionsSubRanges[srn-1];
				sn=subRange->numberOfSuffixes;
				for (int64_t i=0;i<sn;i++){
					uint64_t suffixStartingPosition=subRange->suffixes[i];
					if (markingRate > 0 && suffixStartingPosition%markingRate==0)
						 markedRows->insert(UIntIntPair(sp,suffixStartingPosition/markingRate));
					if (suffixStartingPosition==0)
						first=sp;
					uint64_t tp;
					if (suffixStartingPosition>0)
						tp = suffixStartingPosition-1;
					else
						tp = n -1;  //the super-text is seen as a circular tape
					bwt[sp]=text[tp];
					sp++;
				}

				for (int64_t sri=srn-2;sri>=0;sri--){
					LongRepetitionsRange* subRange=&range->longRepetitionsSubRanges[sri];
					//aggiunge i suffissi precedenti al punto di inserzione
					sn=subRange->numberOfSuffixes;
					for (int64_t i=subRange->insertionPoint;i<sn;i++){
						uint64_t suffixStartingPosition=subRange->suffixes[i];
						if (markingRate > 0 && suffixStartingPosition%markingRate==0)
							 markedRows->insert(UIntIntPair(sp,suffixStartingPosition/markingRate));
						if (suffixStartingPosition==0)
							first=sp;
						uint64_t tp;
						if (suffixStartingPosition>0)
							tp = suffixStartingPosition-1;
						else
							tp = n -1;  //the super-text is seen as a circular tape
						bwt[sp]=text[tp];
						sp++;
					}
				}


				//FOLLOWING SUBRANGE
				//Subrange containing the suffixes that follow the long repetitions ones

				sn=range->followingSubRange.numberOfSuffixes;
				for (int64_t i=0;i<sn;i++){
					uint64_t suffixStartingPosition=range->followingSubRange.suffixes[i];
					if (markingRate > 0 && suffixStartingPosition%markingRate==0)
						markedRows->insert(UIntIntPair(sp,suffixStartingPosition/markingRate));
					if (suffixStartingPosition==0)
						first=sp;
					uint64_t tp;
					if (suffixStartingPosition>0)
						tp = suffixStartingPosition-1;
					else
						tp = n -1;  //the super-text is seen as a circular tape
					bwt[sp]=text[tp];
					sp++;
				}
			}
			else{
				//The range has been ordered as a whole
				int64_t sn=range->numberOfSuffixes;
				for (int64_t i=0;i<sn;i++){
					uint64_t suffixStartingPosition=range->suffixes[i];
					if (markingRate > 0 && suffixStartingPosition%markingRate==0)
						markedRows->insert(UIntIntPair(sp,suffixStartingPosition/markingRate));
					if (suffixStartingPosition==0)
						first=sp;

					uint64_t tp;
					if (suffixStartingPosition>0)
						tp = suffixStartingPosition-1;
					else
						tp = n -1;  //the super-text is seen as a circular tape

					bwt[sp]=text[tp];
					sp++;
				}

			}
		}
		endTime=std::chrono::system_clock::now();
		bwtComputationTime=std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

		#if INDEX_DEBUG_LEVEL > 0
		cout << "\t\tNumber of marked rows: " << markedRows->size() <<endl;
		#endif
	}


void FastBWTransformer::sortNoLongRepetitionsRanges(){
		std::chrono::time_point<std::chrono::system_clock> currentTime=std::chrono::system_clock::now();
		std::time_t printableTime= std::chrono::system_clock::to_time_t(currentTime);
		char mbstr[100];
		std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
		cout << "\tSorting ranges not containing long repetitions, starting at " << mbstr <<endl;

		uint32_t nt = r; //Allocate as many threads as the number of the microprocessor's cores to use

		std::vector<std::thread*> threads;
		std::vector<NoLongRepetitionsSorter*> sorters;

		for (uint32_t t = 0; t < nt; t++) {
			NoLongRepetitionsSorter* sorter=new NoLongRepetitionsSorter(t,ranges,intervals[t].getLeft(),intervals[t].getRight(),text,n,superAlphabet,repetitions,r);
			sorters.push_back(sorter);
			threads.push_back(new std::thread(&NoLongRepetitionsSorter::run, sorter));
		}

		for (uint32_t i = 0; i < nt; i++)
			threads[i]->join();
		for (uint32_t i = 0; i < nt; i++)
			delete threads[i];
		for (uint32_t i = 0; i < nt; i++)
			delete sorters[i];

	}



	void FastBWTransformer::sortLongRepetitionsRanges(){
		std::chrono::time_point<std::chrono::system_clock> currentTime=std::chrono::system_clock::now();
		std::time_t printableTime= std::chrono::system_clock::to_time_t(currentTime);
		char mbstr[100];
		std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
		cout << "\tSorting ranges containing long repetitions, starting at " << mbstr <<endl;

		//Initialize the random numbers generation engine
		myclock::time_point beginning = myclock::now();
		myclock::duration d = myclock::now() - beginning;
		unsigned seed = d.count();
		minstd_rand0* randomGenerator=new minstd_rand0(seed);

		for (uint32_t r=0;r<nonEmptyRangesNumber;r++){
			AlphabetRange* range=ranges[r];
			//double rangeEffectiveSortTime;
			if (range->containsLongRepetitionsSuffixes && range->numberOfSuffixes>1){
					#if INDEX_DEBUG_LEVEL > 0
					cout << "\t\tSorting range " << r << endl;
					#endif
					//it's a range containing a long repetition character and so it can't be orderered as a whole
					//it's necessary to split the range in several subranges:
					//1) the preceding range, containing the suffixes starting with supercharacters preceding the repeated supercharacter
					//     according to the natural superalphabet order;
					//2) the succeedingRange, containing the suffixes starting with supercharacters succeeding the repeated supercharacter;
					//3) a set of subranges , containing the suffixes starting with the repeated supercharacter, so that each of them contain
					//     only the suffixes starting with a number of repetead values contained in an interval ]i,j] of length subIntervalLength.
					//   This assures that at most subIntervalLength comparisons are required to order a couple of suffixes

					/*
					cout << "Sorting long rep. range " << r << endl;
					cout << "\tRange"<<r<<": repeated character: " << repeatedCharacter << endl;
					cout << "\tRange"<<r<<": suffixes number: " << range->numberOfSuffixes << endl;
					cout << "N repetitions" <<endl;
					cout << repetitions['N'] << endl;*/
					uint32_t repeatedSuperCharacterCode=range->longRepeatedCharacter;
					uint32_t mrl=repetitions[repeatedSuperCharacterCode].maximumRepetitionLength;
					uint32_t quotient=mrl/SUBINTERVAL_LENGTH;
					uint32_t remainder=mrl-quotient;
					//int longRepetitionsSubRangesCount=mrl/subIntervalLength + (mrl%subIntervalLength==0?0:1);
					uint32_t longRepetitionsSubRangesCount=quotient + (remainder==0?0:1);
					range->numberOfLongRepetitionsSubranges=longRepetitionsSubRangesCount;
					range->longRepetitionsSubRanges=new LongRepetitionsRange[longRepetitionsSubRangesCount];
					for (uint32_t sri=0;sri<longRepetitionsSubRangesCount;sri++){
						range->longRepetitionsSubRanges[sri].subRangeStart=sri*SUBINTERVAL_LENGTH;
					}


					/*Assign a suffix s to a subRange sri if s contains:
					  - more than subRanges[sri].subRangeStart initial repetead values;
					  - but no more than subRanges[sri+1].subRangeStart repeated values;
					*/
					//Create an array containing all the repetitions to have better performance in the following code
					uint32_t repeatedCharacterRepetitionsSize=repetitions[repeatedSuperCharacterCode].repetitions.size();  //size non
					Repetition repeatedCharacterRepetitions[repeatedCharacterRepetitionsSize];
					for (uint32_t h=0;h<repeatedCharacterRepetitionsSize;h++)
						repeatedCharacterRepetitions[h]=repetitions[repeatedSuperCharacterCode].repetitions[h];


					for (uint64_t si=0;si<range->numberOfSuffixes;si++){
						uint64_t s=range->suffixes[si];

						uint32_t firstSuperCharacterCode=text[s];
						if (firstSuperCharacterCode < repeatedSuperCharacterCode){
							range->precedingSubRange.addSuffix(s);
						} else if (firstSuperCharacterCode > repeatedSuperCharacterCode){
							range->followingSubRange.addSuffix(s);
						} else{
							//search for a repetition which includes the suffix s starting position
							Repetition *rep=NULL;
							bool beforeNextRepetition=false;
							uint32_t h=0;
							while (h<repeatedCharacterRepetitionsSize && !rep && !beforeNextRepetition){
								Repetition tryRep=repeatedCharacterRepetitions[h];
								if (s<tryRep.start)
									beforeNextRepetition=true;  //il suffisso viene prima della prossima ripetizione
																//e quindi la sua parte iniziale non contiene una ripetizione
								else{
									if (s<tryRep.start+tryRep.length)
										rep=&tryRep;
									else
										h++;
								}
							}
							int32_t sri;
							uint64_t numberOfLongRepeatedChars;
							if (!rep){
								sri=0;   //i suffissi che hanno il primo carattere non rientrante in una long repetition
										 //vanno assegnati al subrange 0, ossia a quello contenente i suffissi che
										 //iniziano con il più piccolo numero di caratteri con codice repeatedCharacterCode (allN)
								numberOfLongRepeatedChars=0;
							}
							else{
								//calcola il numero dei caratteri fino alla fine della ripetizione (incluso quello
								//da cui inizia il suffisso), ossia il numero di caratteri ripetuti
								//con cui inizia il suffisso ed è quello che serve per assegnarlo ad un subrange
								numberOfLongRepeatedChars=(rep->start+rep->length) - s;
								int64_t p=0; int64_t u=longRepetitionsSubRangesCount-1; sri=-1;
								while (sri==-1 && p<=u){
									int64_t m=(p+u)/2;
									if (numberOfLongRepeatedChars <= range->longRepetitionsSubRanges[m].subRangeStart)
											u=m-1;
									else if (m==(longRepetitionsSubRangesCount-1) ||
												 numberOfLongRepeatedChars <= range->longRepetitionsSubRanges[m+1].subRangeStart){
											sri=m;
									}
									else
											p=m+1;
								}
								if (sri==-1)
									sri=p;
							}
							LongRepetitionsRange *subRange=&range->longRepetitionsSubRanges[sri];
							subRange->addSuffix(s);

							//L'informazione seguente serve a determinare l'insertion point, cioè il punto in cui dovranno essere inseriti i suffissi
							//del range successivo: l'insertion point coincide con la posizione del primo suffisso del range sri che segue tutti quelli del range sri+1
							uint64_t tp=s+numberOfLongRepeatedChars;
							if (tp >= n)
								   tp = tp - n;
							subRange->firstNonRepeatedGreaterThanRepeated.insert(IntBoolPair(s,text[tp] > range->longRepeatedCharacter));
						}
					}

					//dispose the suffixes array

					delete[] range->suffixes;
					range->suffixes=NULL;

					if (range->precedingSubRange.numberOfSuffixes>1){
						//cout << "Range " << r << ": sorting preceding subrange" << "("<< range->precedingSubRange.numberOfSuffixes  <<")"<< endl;
							ssort(text,n,randomGenerator,range->precedingSubRange.suffixes,range->precedingSubRange.numberOfSuffixes,0,1);
					}

					//sort the sub-ranges (multi-threaded version)
					uint32_t nt=numberOfThreads;  //Allocate as many threads as the number of cores to use
					if (longRepetitionsSubRangesCount < nt)
						nt=longRepetitionsSubRangesCount;
					uint32_t perThreadSubRanges=longRepetitionsSubRangesCount/nt;
					std::vector<std::thread*> threads;
					std::vector<LongRepetitionsSorter*> sorters;
					for (uint32_t t = 0; t < nt; t++){
						//cout << (t*perThreadSubRanges) << "," << (t==nt-1?longRepetitionsSubRangesCount-1:(t+1)*perThreadSubRanges-1) << endl;

						LongRepetitionsSorter* sorter=new LongRepetitionsSorter(t,text,n);
						sorters.push_back(sorter);
						threads.push_back(new std::thread(&LongRepetitionsSorter::sortSubRanges,
															sorter,range,
															t*perThreadSubRanges,
															t==nt-1?longRepetitionsSubRangesCount-1:(t+1)*perThreadSubRanges-1,
															longRepetitionsSubRangesCount));
					}
					for (uint32_t i = 0; i < nt; i++)
						threads[i]->join();
					for (uint32_t i = 0; i < nt; i++)
						delete threads[i];
					for (uint32_t i = 0; i < nt; i++)
						delete sorters[i];

					if (range->followingSubRange.numberOfSuffixes>1){
						//cout << "Range " << r << ": sorting following subrange" << "("<< range->followingSubRange.numberOfSuffixes  <<")"<< endl;
						ssort(text,n,randomGenerator,range->followingSubRange.suffixes,range->followingSubRange.numberOfSuffixes,0,1);
					}

				}
			 }
	}


	void FastBWTransformer::computeBWT(uint32_t *text,uint64_t n,uint8_t markingRate){
			this->text=text;
			this->n=n;
			this->markingRate=markingRate;

			split();
			markLongRepetitionsRanges();

			std::chrono::time_point<std::chrono::system_clock> startTime,endTime;
			startTime=std::chrono::system_clock::now();

			sortNoLongRepetitionsRanges();
			sortLongRepetitionsRanges();

			endTime=std::chrono::system_clock::now();
			sortingTime=std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

			mergeSortResultsAndComputeBWT();

	}



	/**
	 *
	 * @param range			the alphabet range containing the subrange to process
	 * @param sri			the index of the subrange to process
	 * @returns the position of the next bwt character (the new starting position sp)
	 * @throws Exception
	 */

	/*
	void FastBWTransformer::processLongRepetitionsSubRange(AlphabetRange* range,int sri){
		LongRepetitionsRange* subRange=&range->longRepetitionsSubRanges[sri];
		int srn=range->numberOfLongRepetitionsSubranges;
		if (sri<srn-1){
			//aggiunge i suffissi precedenti al punto di inserzione
			int sn=subRange->numberOfSuffixes;
			for (int i=0;i<subRange->insertionPoint;i++){
				int suffixStartingPosition=subRange->suffixes[i];
				if (markingRate > 0 && suffixStartingPosition%markingRate==0)
					 markedRows->insert(IntIntPair(sp,suffixStartingPosition/markingRate));
				if (suffixStartingPosition==0)
					first=sp;
				int tp=suffixStartingPosition + n-1;
				if (tp >= n)
					tp = tp - n;
				bwt[sp]=text[tp];
				sp++;
			}


			//invoca, se esiste un subrange successivo, processSubRange(sri+1)
			if (sri<srn-1)
				processLongRepetitionsSubRange(range,sri+1);

			//aggiunge i suffissi che si susseguono dal punto di inserzione in poi
			for (int i=subRange->insertionPoint;i<sn;i++){
				int suffixStartingPosition=subRange->suffixes[i];
				if (markingRate > 0 && suffixStartingPosition%markingRate==0)
					 markedRows->insert(IntIntPair(sp,suffixStartingPosition/markingRate));
				if (suffixStartingPosition==0)
					first=sp;
				int tp=suffixStartingPosition + n-1;
				if (tp >= n)
					tp = tp - n;
				bwt[sp]=text[tp];
				sp++;
			}

		} else{  //per l'ultimo subRange non esiste alcun punto di inserzione
			int sn=subRange->numberOfSuffixes;
			for (int i=0;i<sn;i++){
				int suffixStartingPosition=subRange->suffixes[i];
				if (markingRate > 0 && suffixStartingPosition%markingRate==0)
					 markedRows->insert(IntIntPair(sp,suffixStartingPosition/markingRate));
				if (suffixStartingPosition==0)
					first=sp;
				int tp=suffixStartingPosition + n-1;
				if (tp >= n)
					tp = tp - n;
				bwt[sp]=text[tp];
				sp++;
			}
		}
	}

    */




} /* namespace std */
