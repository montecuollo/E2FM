/*
 * FastBWTransformer.h
 *
 *  Created on: 22/dic/2014
 *      Author: fernando
 */

#ifndef FASTBWTRANSFORMER_H_
#define FASTBWTRANSFORMER_H_

#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <map>
#include <cmath>
#include <thread>
#include <mutex>
#include <fstream>
#include <cerrno>
#include <unordered_map>
#include <vector>
#include <random>
#include <libgen.h>
#include <cstring>
#include "ScrambledSuperAlphabet.h"


namespace std {


#ifndef  MINIMUM_REPETITION_LENGTH
#define MINIMUM_REPETITION_LENGTH 1000
#endif

#ifndef  SUBINTERVAL_LENGTH
#define SUBINTERVAL_LENGTH 1000
#endif


typedef pair<string, long> StringLongPair;

typedef pair<uint64_t, bool> IntBoolPair;

class ArrayInterval {
private:
	uint32_t left;
	uint32_t right;
	double weight;
	double relativeError;

public:
	ArrayInterval(uint32_t left, uint32_t right, double weight, double relativeError) {
		this->left = left;
		this->right = right;
		this->weight = weight;
		this->relativeError = relativeError;
	}

	ArrayInterval(){}

	uint32_t getLeft() {
		return left;
	}
	void setLeft(uint32_t left) {
		this->left = left;
	}
	int getRight() {
		return right;
	}
	void setRight(uint32_t right) {
		this->right = right;
	}

	double getWeight() {
		return weight;
	}

	void setWeight(double weight) {
		this->weight = weight;
	}

	double getRelativeError() {
		return relativeError;
	}

	void setRelativeError(double relativeError) {
		this->relativeError = relativeError;
	}

	friend ostream& operator<<(ostream& os, const ArrayInterval& dt);

	friend class  ArraySplitter;
};



/**
 * A range of suffixes
 * @author fernando
 *
 */
class Range {
public:
	uint64_t numberOfSuffixes; //number of the suffixes really contained in this range
	uint64_t *suffixes;  //suffixes contained in this range: its dimension is greater
					     //than the number of the suffixes really contained in this range
	uint64_t suffixesSize;
	mutex *thisRangeMutex;

	void addSuffix(uint64_t startPosition) {
		thisRangeMutex->lock();
		uint64_t length = suffixesSize;
		if (numberOfSuffixes > suffixesSize - 1) {
			suffixesSize = numberOfSuffixes * 2;
			uint64_t *newArray = new uint64_t[suffixesSize];
			std::copy(suffixes, suffixes + length, newArray);
			delete[] suffixes;
			suffixes = newArray;
		}
		suffixes[numberOfSuffixes] = startPosition;
		numberOfSuffixes++;
		thisRangeMutex->unlock();
	}

	void addSuffixes(int64_t suffixesToAdd[], int64_t numberOfSuffixesToAdd) {
		std::copy(suffixesToAdd, suffixesToAdd + numberOfSuffixesToAdd,suffixes);
		numberOfSuffixes += numberOfSuffixesToAdd;
	}

	Range() {
		suffixesSize = 1000;
		suffixes = new uint64_t[suffixesSize];
		numberOfSuffixes = 0;
		thisRangeMutex = new mutex();
	}

	~Range(){
		if (suffixes)
			delete[] suffixes;
		delete thisRangeMutex;
	}
};

class LongRepetitionsRange: public Range {
public:
	/**
		 * L'HashMap seguente contiene, per ogni suffisso s, la seguente informazione booleana: il primo carattere non ripetuto segue nell'ordinamento
		 * il carattere che costituisce le ripetizioni lunghe?
		 */
	unordered_map<uint64_t, bool> firstNonRepeatedGreaterThanRepeated;
	uint64_t subRangeStart;
	int64_t insertionPoint; //-1 if not assigned

	/*void add(){
	 firstNonRepeatedGreaterThanRepeated.insert(IntBoolPair(1,true));
	 }*/

};

class AlphabetRange: public Range {
public:
	uint32_t firstCharacter;
	string firstCharacterSymbol;
	uint32_t lastCharacter;
	string lastCharacterSymbol;

	/**
	 * Contiene suffissi con lunghe ripetizioni
	 * **/
	bool containsLongRepetitionsSuffixes;

	/**
	 * Carattere ripetuto
	 */
	uint32_t longRepeatedCharacter;

	/**
	 * Sub-range costituito dai soli suffissi che precedono quelli con una lunga ripetizione
	 */
	Range precedingSubRange;

	/**
	 * Array dei sub-range contenenti ripetizioni lunghe
	 */
	LongRepetitionsRange* longRepetitionsSubRanges;
	uint32_t numberOfLongRepetitionsSubranges;

	/**
	 * Sub-range costituito dai soli suffissi che seguono quelli con una lunga ripetizione
	 */
	Range followingSubRange;

	AlphabetRange() :
			Range() {
		containsLongRepetitionsSuffixes = false;
		numberOfLongRepetitionsSubranges = 0;
	}

	friend ostream& operator<<(ostream& os, const AlphabetRange& r);

};


class ArraySplitter {

private:
	int64_t *array;			//array to split
	int64_t arrayLength;	//length of the array
	int numberOfIntervals;     			//number of intervals
	int *splitters;  			//positions of the numberOfIntervals-1 splitters
								//the splitter between two intervals is positioned on the first element
								//of the second interval
	int maxSteps;	  				//maximum number of steps
	double minPercentualGain;   	//minimum objective function gain
	double optimalWeight;			//optimal weight
	int steps;				//steps done by the splitting algorithm



	double computeIntervalWeight(int interval, int splitters[]) {
		int startPosition, endPosition;
		if (interval == 0)
			startPosition = 0;
		else
			startPosition = splitters[interval - 1];

		if (interval == numberOfIntervals-1)
			endPosition = arrayLength - 1;
		else
			endPosition = splitters[interval] - 1;
		double sum = 0;
		for (int i = startPosition; i <= endPosition; i++)
			sum = sum + array[i];
		return sum;
	}

	/**
	 * This objective function is the square quadratic error between all the interval weigths and the optimal weight,
	 * corresponding to an uniform weights distribution
	 * @param splitters
	 * @param optimalValue
	 * @return
	 */
	double computeObjectiveFunction(int splitters[], double optimalValue) {
		double result = 0;
		for (int interval = 0; interval < numberOfIntervals; interval++)
			result += pow(
					computeIntervalWeight(interval, splitters) - optimalValue,
					2);
		return result;
	}

	/**
	 * This objective is the sum of the absolute values of the relative errors between all the interval weigths and the optimal weight,
	 * corresponding to an uniform weights distribution
	 * @param splitters
	 * @param optimalValue
	 * @return
	 */
	/*private double computeObjectiveFunction(int[] splitters,double optimalValue){
	 double result=0;
	 for (int interval=0;interval<numberOfIntervals;interval++)
	 result+=Math.abs(computeIntervalWeight(interval, splitters) - optimalValue);
	 return result;
	 }*/

public:
	ArraySplitter(){
		array=NULL;
		splitters=NULL;
	}

	~ArraySplitter(){
		if (splitters)
			delete[] splitters;
	}

	void setNumberOfIntervals(int numberOfIntervals) {
		this->numberOfIntervals = numberOfIntervals;
	}

	int getNumberOfIntervals() {
		return numberOfIntervals;
	}

	void split() {
		//computes the sum of all the weights
		int n = 0;
		for (int i = 0; i < arrayLength; i++)
			n += array[i];
		int pn = numberOfIntervals;  					//intervals number
		int sn = pn - 1;									//splitters number
		//initial splitters assignment
		splitters = new int[sn];
		double intervalLength = ((double) arrayLength) / pn;
		for (int j = 0; j < sn; j++)
			splitters[j] = (int) round((j + 1) * intervalLength);
		optimalWeight = (double) n / pn;

		steps = 0;
		double gain = 1;
		while (steps < maxSteps && gain >= minPercentualGain) {
			double objectiveFunctionInitialValue = computeObjectiveFunction(
					splitters, optimalWeight);

			double maxGain = 0;
			int *maxGainSplitters = new int[sn];
			std::copy(splitters, splitters + sn, maxGainSplitters);
			for (int j = 0; j < sn; j++) {  //for each splitter
				//try to move left
				int leftmostPosition; //so that at least one element remains in the previous interval
				if (j == 0)
					leftmostPosition = 1;
				else
					leftmostPosition = splitters[j - 1] + 1;
				for (int trialPosition = splitters[j] - 1;
						trialPosition >= leftmostPosition; trialPosition--) {
					int *trialSplitters = new int[sn];
					std::copy(splitters, splitters + sn, trialSplitters);
					trialSplitters[j] = trialPosition;
					double objectiveFunctionTrialValue =
							computeObjectiveFunction(trialSplitters,
									optimalWeight);
					gain = objectiveFunctionInitialValue
							- objectiveFunctionTrialValue;
					if (gain > maxGain) {
						delete[] maxGainSplitters;
						maxGainSplitters = trialSplitters;
						maxGain = gain;
					} else
						delete[] trialSplitters;
				}
				//try to move right
				int rightmostPosition; //so that at least one element remains in the next interval
				if (j == sn - 1)
					rightmostPosition = arrayLength - 1;
				else
					rightmostPosition = splitters[j + 1] - 1;
				for (int trialPosition = splitters[j] + 1;
						trialPosition <= rightmostPosition; trialPosition++) {
					int *trialSplitters = new int[sn];
					std::copy(splitters, splitters + sn, trialSplitters);
					trialSplitters[j] = trialPosition;
					double objectiveFunctionTrialValue =
							computeObjectiveFunction(trialSplitters,
									optimalWeight);
					gain = objectiveFunctionInitialValue
							- objectiveFunctionTrialValue;
					if (gain > maxGain) {
						delete[] maxGainSplitters;
						maxGainSplitters = trialSplitters;
						maxGain = gain;
					} else
						delete[] trialSplitters;
				}
			}
			if (maxGain > 0) {
				delete[] splitters;
				splitters = maxGainSplitters;
				gain = maxGain / objectiveFunctionInitialValue;
			}
			else{
				delete[] maxGainSplitters;
			}
			steps++;
		}
	}

	int* getSplitters() {
		return splitters;
	}

	int getMaxSteps() {
		return maxSteps;
	}

	void setMaxSteps(int maxSteps) {
		this->maxSteps = maxSteps;
	}

	double getMinPercentualGain() {
		return minPercentualGain;
	}

	void setMinPercentualGain(double minPercentualGain) {
		this->minPercentualGain = minPercentualGain;
	}

	void setArray(int64_t *array, int64_t arrayLength) {
		this->array = array;
		this->arrayLength = arrayLength;
	}

	int64_t* getArray() {
		return array;
	}

	double getOptimalWeight() {
		return optimalWeight;
	}

	int getDoneSteps() {
		return steps;
	}

	ArrayInterval* getIntervals() {
		ArrayInterval* intervals = new ArrayInterval[numberOfIntervals];
		int left,right;
		for (int interval = 0; interval < numberOfIntervals; interval++){
			if (interval == 0)
				left = 0;
			else
				left = splitters[interval - 1];

			if (interval == numberOfIntervals - 1)
				right = arrayLength - 1;
			else
				right = splitters[interval] - 1;
			double weight = 0;
			for (int i = left; i <= right; i++)
				weight += array[i];
			double relativeError = abs((weight - optimalWeight) / optimalWeight);
			intervals[interval].left=left;
			intervals[interval].right=right;
			intervals[interval].weight=weight;
			intervals[interval].relativeError=relativeError;
		}
		return intervals;

	}

	ArrayInterval getInterval(int interval) {
		int left, right;
		if (interval == 0)
			left = 0;
		else
			left = splitters[interval - 1];

		if (interval == numberOfIntervals - 1)
			right = arrayLength - 1;
		else
			right = splitters[interval] - 1;
		double weight = 0;
		for (int i = left; i <= right; i++)
			weight += array[i];

		double relativeError = abs((weight - optimalWeight) / optimalWeight);
		ArrayInterval retVal(left, right, weight, relativeError);
		return retVal;
	}
};

class Repetition {
private:

public:
	int64_t start;  //index in text from which repetition starts
	int64_t length; //length of repetition (number of repeated characters)

	Repetition(int64_t start, int64_t length) {
		this->start = start;
		this->length = length;
	}

	Repetition(){
	}
	friend ostream& operator<<(ostream& os, const Repetition& dt);

};


class RepetitionInformations {
public:
	uint64_t maximumRepetitionLength;
	vector<Repetition> repetitions;

	RepetitionInformations() {
		maximumRepetitionLength = 0;
	}

	friend ostream& operator<<(ostream& os, const RepetitionInformations& dt);

};



typedef pair<uint32_t, RepetitionInformations> SymbolRepInfoPair;

class Distributor {
	public:
		uint64_t firstTextPosition;
		uint64_t lastTextPosition;
		uint32_t *text;
		ScrambledSuperAlphabet *superAlphabet;
		uint32_t er;
		uint32_t elementaryRangesLength;  //TODO: verificare che non servano 64 bit
		AlphabetRange** elementaryRanges;

		unordered_map<uint32_t, RepetitionInformations> repetitions;

		Distributor(uint64_t firstTextPosition,uint64_t lastTextPosition,uint32_t *text,ScrambledSuperAlphabet *superAlphabet,
						uint32_t er,uint32_t elementaryRangesLength,
						AlphabetRange **elementaryRanges);


		void run();

		void distribute();

};

class NoLongRepetitionsSorter {
		public:
			int threadNumber;
			AlphabetRange **ranges;
			uint32_t firstRange;
			uint32_t lastRange;
			uint32_t *text;
			uint64_t n; //text length
			ScrambledSuperAlphabet* superAlphabet;
			unordered_map<uint32_t, RepetitionInformations> repetitions;
			uint32_t numberOfThreads;
			minstd_rand0* randomGenerator;

			NoLongRepetitionsSorter(uint32_t threadNumber,AlphabetRange **ranges,uint32_t firstRange,
								uint32_t lastRange,uint32_t *text,uint64_t n,ScrambledSuperAlphabet* superAlphabet,
								unordered_map<uint32_t, RepetitionInformations> repetitions,uint32_t numberOfThreads);
			~NoLongRepetitionsSorter();
			void sort();
			void run();
};

class LongRepetitionsSorter {
		public:
			int threadNumber;
			uint32_t *text;
			uint64_t n; //text length
			minstd_rand0* randomGenerator;

			LongRepetitionsSorter(uint32_t threadNumber,uint32_t *text,uint64_t n);
			~LongRepetitionsSorter();
			void sortSubRanges(AlphabetRange *range,uint32_t firstSubRange,uint32_t lastSubRange,uint32_t longRepetitionsSubRangesCount);
			void sort();
			void run();
};




class FastBWTransformer {
private:
	uint32_t numberOfThreads;
	string fastaFileDirectory;
	string fastaFileNameWithoutExtension;
	uint32_t *text;
	uint32_t *bwt;
	uint8_t markingRate;
	uint64_t first; //pointer to the suffix array's first element (i.e. to the last text character);
	uint64_t sp; //used during the BWT computation
	ScrambledSuperAlphabet* superAlphabet;

	uint64_t n;	 //text length
	uint8_t k;   //order of the alphabet extension

	uint32_t r; //number of ranges of the suffix array whose suffixes can be ordered separately
	uint32_t er; //number of elementary ranges of the suffix array, composed of suffixes starting with super-characters contained in a range of the superalphabet
			//these ranges are aggregated to form the r desired ranges
	uint32_t elementaryRangesLength; //number of the extended alphabet characters falling into an elementary range
	uint64_t suffixesMaximumNumberPerRange; //maximum number of suffixes falling in each elementary range
	double risePercentage;

	AlphabetRange** elementaryRanges;
	uint32_t nonEmptyRangesNumber; //number of the elementary ranges containing at least one suffix
	mutex nonEmptyRangesNumberMutex; //mutex that protects the shared access to the above member variable
	AlphabetRange** ranges;//Ranges of the suffix array whose suffixes will be ordered separately: they are obtained
						  //aggregating consecutive elementary range so that the "weigth" of each range,
						  //defined as the number of its suffixes, is almost the same.

	ArrayInterval* intervals;//Intervals of the ranges array that will be sorted by separate threads

	unordered_map<uint32_t, RepetitionInformations> repetitions;
	map<uint64_t,uint64_t>  *markedRows;

	bool computeStatistics;
	double splittingTime;
	double distributionTime;
	double loadingTime;
	double sortingTime;
	double bwtComputationTime;
	double savingTime;
	double optimalWeight;
	double meanRangesRelativeError;
	double maxRangesRelativeError;
	double meanElementaryRangesRelativeError;
	double maxElementaryRangesRelativeError;
	uint32_t splittingSteps;




public:

	/**
	 *
	 * @param numberOfThreads	Number of parallel threads to use
	 * @param superAlphabetOrder Order of the alphabet extension used to represent the original sequence in order to speed up
	 *                           the BWT computation
	 * @param numberOfElementaryRanges	Number of the elementary ranges; each of them is made only of suffixes starting with super-characters contained
	 *                                  in a corresponding range of the superalphabet's characters
	 *                					The elementary ranges are finally aggregated to form the r ranges that will be separately orderer
	 * @throws Exception
	 */
	FastBWTransformer(uint32_t numberOfThreads, ScrambledSuperAlphabet *superAlphabet,
			uint32_t numberOfElementaryRanges) {
		this->numberOfThreads = numberOfThreads;
		r = numberOfThreads;
		this->superAlphabet = superAlphabet;
		er = numberOfElementaryRanges;
		k = superAlphabet->getOrder();
		risePercentage = 0.50;

		if (superAlphabet->getSize() < er)
			throw logic_error(
					"The extended alphabet's size must be greater than the number of elementary ranges");
		computeStatistics = false;
		splittingTime=0;
		distributionTime=0;
		loadingTime=0;
		sortingTime=0;
		bwtComputationTime=0;
		savingTime=0;
		optimalWeight=0;
		meanRangesRelativeError=0;
		maxRangesRelativeError=0;
		meanElementaryRangesRelativeError=0;
		maxElementaryRangesRelativeError=0;


		elementaryRanges=NULL;
		ranges=NULL;
		intervals=NULL;
		nonEmptyRangesNumber=0;

	}

	~FastBWTransformer(){
		if (elementaryRanges){
			for (uint32_t i=0;i<er;i++){
				AlphabetRange* range=elementaryRanges[i];
				if (range)
					delete range;
			}
			delete[] elementaryRanges;
		}

		if (ranges)
			delete[] ranges;

		if (intervals)
			delete[] intervals;

	}

	void mergeSortResultsAndComputeBWT();


	void computeBWT(uint32_t *text,uint64_t n,uint8_t markingRate);

	uint32_t *getBWT(){
		return bwt;
	}

	uint64_t getLastCharacterPosition(){
		return first;
	}

	map<uint64_t,uint64_t> *getMarkedRows(){
		return markedRows;
	}

	void vecswap(uint32_t *x,uint64_t i, uint64_t j, uint64_t n){
		while (n-- > 0) {
	        //Expansion of swap(x,i, j);
			uint32_t t=x[i];
			x[i]=x[j];
			x[j]=t;
	        i++;
	        j++;
	    }
	}


	void split();

	void markLongRepetitionsRanges();

	void sortLongRepetitionsRanges();

	void sortNoLongRepetitionsRanges();


	/**
		 *
		 * @param range			the alphabet range containing the subrange to process
		 * @param sri			the index of the subrange to process
		 * @returns the position of the next bwt character (the new starting position sp)
		 * @throws Exception
		 */

	void processLongRepetitionsSubRange(AlphabetRange* range,int sri);

	uint32_t* invertBWT(){
		uint64_t lFirst=first;  //it's the pointer to the first suffix of the suffix array and then it points to the last
						   //character of the text
		uint64_t *translate = new uint64_t[n];
		for (uint64_t i=0;i<n;i++)
			translate[i]=0;
		uint32_t s=superAlphabet->getSize();

		int64_t *sums = new int64_t[s];
		for (uint64_t i=0;i<s;i++)
			sums[i]=0;
		for(uint64_t k=0; k < n; k++){
			sums[bwt[k]+1]++;
		}

		for(uint64_t k=1; k < s; k++){
			sums[k] += sums[k-1];
		}
		for(uint64_t k=0; k < n; k++){
			uint64_t index = sums[bwt[k]]++;
			translate[index] = k;
		}
		//lFirst = translate[lFirst];
		uint32_t *reconstruct = new uint32_t[n];
		for (uint64_t k = 0; k <n; k++) {
			lFirst = translate[lFirst];
			reconstruct[k] = bwt[lFirst];
		}

		delete[] sums;
		delete[] translate;

		return reconstruct;
	}

	void verifyBWT(){
		std::chrono::time_point<std::chrono::system_clock> currentTime;
		std::time_t printableTime;
		currentTime=std::chrono::system_clock::now();
		printableTime= std::chrono::system_clock::to_time_t(currentTime);
		char mbstr[100];
		std::strftime(mbstr, sizeof(mbstr), "%H:%M:%S", std::localtime(&printableTime));
		cout << "Verifying BWT, starting at " << mbstr << endl;
		uint32_t* inverse=invertBWT();
		uint64_t i=0;
		int64_t mismatchAtPosition=-1;
		while(i<n && mismatchAtPosition==-1)
			if (!(text[i]==inverse[i]))
				mismatchAtPosition=i;
			else
				i++;
		if (mismatchAtPosition==-1)
			cout << "\tBWT is OK" << endl;
		else
			cout << "\tBWT inversion ERROR: Mismatch at position " << mismatchAtPosition << endl;
		delete[] inverse;
	}


		double getLoadingTime() {
			return loadingTime;
		}

		double getSplittingTime() {
			return splittingTime;
		}

		double getDistributionTime() {
			return distributionTime;
		}

		double getSortingTime() {
			return sortingTime;
		}

		double getBwtComputationTime() {
			return bwtComputationTime;
		}

		double getSavingTime() {
			return savingTime;
		}

		uint32_t getNonEmptyRangesNumber() {
			return nonEmptyRangesNumber;
	    }

		uint32_t getSplittingSteps() {
			return splittingSteps;
		}


		double getOptimalWeight() {
			return optimalWeight;
		}

		double getMeanRangesRelativeError() {
			return meanRangesRelativeError;
		}

		double getMaxRangesRelativeError() {
			return maxRangesRelativeError;
		}

		double getMeanElementaryRangesRelativeError() {
			return meanElementaryRangesRelativeError;
		}

		double getMaxElementaryRangesRelativeError() {
			return maxElementaryRangesRelativeError;
		}


		/**
		 * Compares two rotations identified by rotationIndex1 and rotationIndex2,
		 * starting from character depth; code is mutuated from String.compareTo
		 * (as implemented in Open JDK 7)
		 * @param rotationIndex1
		 * @param rotationIndex2
		 * @param depth
		 * @return
		 */
		int64_t compareRotations(uint64_t rotationIndex1,uint64_t rotationIndex2,uint64_t depth){
			  uint64_t i=rotationIndex1+depth;
			  if (i >= n)
				   i = i - n;
			  uint64_t j=rotationIndex2+depth;
			  if (j >= n)
				   j = j - n;
			  uint64_t m=n;  //TODO: verificare il ciclo tramite un breakpoint
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


		void verifySort(){
			bool sortOK=true;
			for (uint32_t rn=0;rn<nonEmptyRangesNumber;rn++){
				Range* r=ranges[rn];
				for (uint64_t i=0;i<r->numberOfSuffixes-1;i++){
					if (compareRotations(r->suffixes[i],r->suffixes[i+1],0)>0){
						cout << "ERROR at range "<< rn << " (suffix "<< i << "/" << r->numberOfSuffixes << ")"<< endl;
						sortOK=false;
						break;
					}
				}
			}
			if (sortOK)
				cout << "All ranges have been sorted correctly" << endl;
		}
};




} /* namespace std */

#endif /* FASTBWTRANSFORMER_H_ */
