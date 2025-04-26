// my_predictor.h
// This file contains my_predictor class.

#include "utils.h"
#include <bitset>
#include <time.h>
#include <cstdlib>

#define GHR_SIZE 131
#define BIMODAL_CTR_MAX 3
#define BIMODAL_CTR_INIT 2
#define TAGEPRED_CTR_MAX 7
#define TAGEPRED_CTR_INIT 0
#define BIMODALLOG 14
#define TAGELOG 13
#define NUMTAGETABLES 4
#define USEFUL_BITS_MAX 3

/////////////// STORAGE BUDGET JUSTIFICATION ////////////////
// Total PHT counters: 2^14 (Bimodal)
// Total PHT size = 2^14 * 2 bits/counter = 2^15 bits = 32KB
// GHR size: 131 bits
// Table size: 4 * 2^13 * (2 + 9 + 3) bits = 4 * 2^13 * 14 bits = 458KB
// Total Size = PHT size + GHR size + Table size
/////////////////////////////////////////////////////////////


//each entry within a TAGE table contains a 3-bit counter, a 8 or 9-bit tag, and 2-bit usefulness counter
struct TageEntry{
	INT32 ctr; //0-7 counter which dictates prediction (4-7 taken, 0-3 not taken),
	UINT32 tag; //tag, hashing function described below
	INT32 usefulBits; // 0-3 counter dictating usefulness of entry (0-3, 0 is least useful, 3 is most useful)
};

struct CompressedHist{
	UINT32 geomLength; //oriiginal length of the history, dictated via geometric sequence 
	UINT32 targetLength; //the desired target length for indexing the TAGE tables
	UINT32 compHist; //the compressed history itself

	void updateCompHist(bitset<GHR_SIZE> GHR){
		//Update the compressed history. The basic idea can be find in the PPM paper

		
		//create masks
		int mask = (1 << targetLength) - 1; //ensure that final result is only targetLength bits long
		int oldGHRMask = GHR[geomLength] << (geomLength % targetLength); //position the oldest GHR bit at position geomLength % targetLength
		int mask2 = 1 << targetLength; //extract the oldest bit in the compressed history

		//update history itself
		compHist = (compHist << 1) + GHR[0]; // this pushes the history to the left, placing the oldest bit out of range and placing the newest GHR bit in the LSB position (newest)
		compHist ^= ((compHist & mask2) >> targetLength); //extract the oldest bit in the comphistory, and create a new mask that has that bit in LSB to XOR with the new GHR bit
		compHist ^= oldGHRMask; //XOR the oldest GHR bit with the bit at position geomLength % targetLength
		compHist &= mask; //ensure that the final result is only targetLength bits long
	}
};

class my_update : public branch_update {
public:
	unsigned int index;
};

class my_predictor : public branch_predictor {
private:
	bitset<GHR_SIZE> GHR; //global history register
	int PHR; //16-bit Path history

	//Base Bimodal predictor
	UINT32* bimodal;
	UINT32 historyLength; //length of the bimodal history
	UINT32 numBimodalEntries; //number of entries in the bimodal predictor
	UINT32 bimodalLog;

	//TAGE tables
	TageEntry* tageTables[NUMTAGETABLES]; //array of 4 pointers to the 4 TAGE tables. Because the entries are of type Entry, this is a Entry* datatype
	UINT32 numTageTableEntries; //number of entries in a TAGE table
	UINT32 tageLog; //log2 of the number of entries in the TAGE table
	UINT32 geometric[NUMTAGETABLES]; //geometric sequence of history length for the TAGE tables

	//Compressed histories, useful for indexing
	CompressedHist compressedIndex[NUMTAGETABLES]; //compressed index for each TAGE table
	CompressedHist compressedTag[2][NUMTAGETABLES]; //use in compressed tag for each TAGE table

	//Prediction related variables
	bool primePred;
	bool altPred;
	int primeTable; //the table that is used for the prime prediction
	int altTable; //the table that is used for the alternative prediction

	//identifying data
	UINT32 index[NUMTAGETABLES]; //index for each TAGE table
	UINT32 tag[NUMTAGETABLES]; //tag for each TAGE table

	//clock for graceful reset over time 
	UINT32 clock;
	int clock_flip;

	INT32 altBetterCount; //counter to determine wheter altPred is more suitable to use. Intuition: newly allocated less likely to be useful

public:
	my_update u;
	branch_info bi; 


	my_predictor (void){ 
		//initialize bimodal predictor
		bimodalLog = BIMODALLOG;
		numBimodalEntries = 1 << bimodalLog;
		bimodal = new UINT32[numBimodalEntries];
		for (int i = 0; i < numBimodalEntries; i++) {
			bimodal[i] = BIMODAL_CTR_INIT; //initialize all entries to 2
		}
		//initialize TAGE tables
		tageLog = TAGELOG;
		numTageTableEntries = 1 << tageLog;
		for (int i = 0; i < NUMTAGETABLES; i++){
			tageTables[i] = new TageEntry[numTageTableEntries];//create each TAGE table with numTageTableEntries entries
			for (int j = 0; j < numTageTableEntries; j++) { //for each entry in the TAGE table
				tageTables[i][j].ctr = TAGEPRED_CTR_INIT; //initialize the counter to 0
				tageTables[i][j].tag = 0; //initialize the tag to 0
				tageTables[i][j].usefulBits = 0; //initialize the useful bits to 0	
			}
		}
		//define geometric sequence of history length used for the TAGE tables
		//geometric[0] is length of history used for rightmost table indexed using most history
		//geometric[3] is length of history used for leftmost table indexed using least history
		geometric[0] = 130;
		geometric[1] = 44;
		geometric[2] = 15;
		geometric[3] = 5;

		//initialize compressed index buffers
		for (int i = 0; i < NUMTAGETABLES; i++) {
			compressedIndex[i].compHist = 0; //initialize compressed history to 0
			compressedIndex[i].geomLength = geometric[i]; //set starting length to corresponding length from geometric series
			compressedIndex[i].targetLength = tageLog;//set target legnth to the log2 of the number of entries in the TAGE table
		}

		//initialize compressed tag buffers
		for (int j = 0; j < 2; j++){ //why 2? see PPM paper, CSR1 and CSR2 are both used in final tag value
			for (int i = 0; i < NUMTAGETABLES; i++){
				compressedTag[j][i].compHist = 0;
				compressedTag[j][i].geomLength = geometric[i];
				if (j == 0){
					compressedTag[j][i].targetLength = 9; //CSR1 is 9 bits long
				} else {
					compressedTag[j][i].targetLength = 8; //CSR2 is 8 bits long
				}
			}
		}

		//initialize prediction related variables
		primePred = -1;
		altPred = -1;
		primeTable = NUMTAGETABLES;
		altTable = NUMTAGETABLES;

		//initialize index and tag
		for (int i = 0; i < NUMTAGETABLES; i++){
			index[i] = 0;
			tag[i] = 0;
		}
		//initialize clock
		clock = 0;
		clock_flip = 0;

		//initialize altBetterCount
		altBetterCount = 8; //start with a value indicating that the altpred is likely to be useful, as the newly allocated entries are less likely to be useful

		//initialize GHR and PHR
		GHR.reset(); //initialize GHR to all 0s
		PHR = 0; //initialize PHR to 0
	}

	branch_update *predict (branch_info & b) {
		bi = b;
		if (b.br_flags & BR_CONDITIONAL) {
			//cout<<"predict"<<endl;
			//define default prediction just in case the predictor fails
			bool final_prediction = NOT_TAKEN;
			UINT32 pc = b.address;

			//base prediciton (bimodal predictor)
			bool basePrediction;
			UINT32 bimodalIndex = pc % numBimodalEntries;
			UINT32 bimodalCounter = bimodal[bimodalIndex];
			if (bimodalCounter > BIMODAL_CTR_MAX / 2) {
				basePrediction = TAKEN; //if counter is greater than 1, predict taken
			} else {
				basePrediction = NOT_TAKEN; //otherwise predict not taken
			}

			//TAGE predictor
			//create tag
			for (int i = 0; i < NUMTAGETABLES; i++) {
				tag[i] = pc ^ compressedTag[0][i].compHist ^ compressedTag[1][i].compHist; //XOR the PC with the compressed histories
				//formula from PPM paper (PC ^ CSR1 ^ CSR2)
				tag[i] &= ((1<<9) - 1); //ensure that the tag is 9 bits long
			}
			//create index
			index[0] = pc ^ (pc >> TAGELOG) ^ compressedIndex[0].compHist ^ PHR ^ (PHR >> TAGELOG); //XOR the PC, PHR, and compressed history
			index[1] = pc ^ (pc >> TAGELOG - 1) ^ compressedIndex[1].compHist ^ PHR; 
			index[2] = pc ^ (pc >> TAGELOG - 2) ^ compressedIndex[2].compHist ^ (PHR & 63); //XOR the PC and compressed history
			index[3] = pc ^ (pc >> TAGELOG - 3) ^ compressedIndex[2].compHist ^ (PHR & 7); //XOR the PC
			//create ensure that index is within the proper size
			for (int i = 0; i < NUMTAGETABLES; i++) {
				index[i] = index[i] & ((1<<TAGELOG) - 1); //ensure that the index is within the bounds of the TAGE table
			}
			//search for the best prediction
			primePred = -1;
			altPred = -1;
			primeTable = NUMTAGETABLES; //initialize prime table to invalid value
			altTable = NUMTAGETABLES; //initialize alt table to invalid value
			//find the table that contains primePred
			for (int i = 0; i < NUMTAGETABLES; i++){
				if (tageTables[i][index[i]].tag == tag[i]){
					primeTable = i;
					break;
				}
			}
			//find the table that contains altPred, this is always at a higher (smaller history used) than primeTable
			for (int i = primeTable + 1; i < NUMTAGETABLES; i++){
				if (tageTables[i][index[i]].tag == tag[i]){
					altTable = i;
					break;
				}
			}
			
			//determine final prediction based on the primeTable and altTable
			if (primeTable < NUMTAGETABLES){ //if a tagged example exists outside of base predictor
				//find altPred
				if (altTable == NUMTAGETABLES){ //no altPred found in TAGE tables, use base predictor
					altPred = basePrediction;
				} else { //altPred found in "smaller" TAGE table
					if (tageTables[altTable][index[altTable]].ctr >= TAGEPRED_CTR_MAX){ //altPred's counter is above taken threshold
						altPred = TAKEN;
					} else {
						altPred = NOT_TAKEN;
					}
				}
				//find primePred
				if ((tageTables[primeTable][index[primeTable]].ctr != 3) || 
					(tageTables[primeTable][index[primeTable]].ctr != 4) || 
					(tageTables[primeTable][index[primeTable]].usefulBits != 0) ||
					(altBetterCount < 8)) { // if the PrimeTable entry is not newly allocated, is useful, and alt is worse
					if (tageTables[primeTable][index[primeTable]].ctr >= TAGEPRED_CTR_MAX/2){ //primePred's counter is above taken threshold
						primePred = TAKEN;
					} else {
						primePred= NOT_TAKEN;
					}
					final_prediction = primePred;
				} else { //above conditions for choosing primePred are not met, so go with altPred instead
					final_prediction = altPred;
				}
			} else { //none found
				altPred = basePrediction;
				final_prediction = altPred;
			}
			u.direction_prediction (final_prediction);
		} else {
			u.direction_prediction (true);
		}
		u.target_prediction (0);
		return &u;
	}

	void update (branch_update *u, bool taken, unsigned int target) {
		if (bi.br_flags & BR_CONDITIONAL) {
			//cout << "update" << endl;
			//rename inputs for clarity
			UINT32 pc = bi.address;
			bool actualDir = taken; //actual direction of the branch
			bool predDir = u->direction_prediction();
			UINT32 branchTarget = target;
			 
			bool strong_old_present = false;
			bool new_entry = 0;
			
			//update TAGE tables
			if (primeTable < NUMTAGETABLES){
				// As per update policy
				// 1st update the useful counter
				if ((predDir != altPred)){ //if prediction is not alt we should update prime
					if (predDir == actualDir){ // if prime was actually correct
						tageTables[primeTable][index[primeTable]].usefulBits = SatIncrement(tageTables[primeTable][index[primeTable]].usefulBits, BIMODAL_CTR_MAX);
					} else { //prime was not correct
						tageTables[primeTable][index[primeTable]].usefulBits = SatDecrement(tageTables[primeTable][index[primeTable]].usefulBits);
					}
				}    
			 	// 2nd update the counters which provided the prediction  
				if(actualDir) { //actual direction taken
					tageTables[primeTable][index[primeTable]].ctr = SatIncrement(tageTables[primeTable][index[primeTable]].ctr, TAGEPRED_CTR_MAX);
				} else { //actual direction not taken
					tageTables[primeTable][index[primeTable]].ctr = SatDecrement(tageTables[primeTable][index[primeTable]].ctr);
				}
			} else { //update base bimodal because prediction did not come from TAGE tables
				UINT32 bimodalIndex  = (pc) % (numBimodalEntries);
				if(actualDir) {
					bimodal[bimodalIndex] = SatIncrement(bimodal[bimodalIndex], BIMODAL_CTR_MAX);
				} else {
					bimodal[bimodalIndex] = SatDecrement(bimodal[bimodalIndex]);
				}
			}
			// Check if the current Entry which gave the prediction is a newly allocated entry or not.
			if (primeTable < NUMTAGETABLES){ 
				if((tageTables[primeTable][index[primeTable]].usefulBits == 0) &&
				((tageTables[primeTable][index[primeTable]].ctr  == 3) || (tageTables[primeTable][index[primeTable]].ctr  == 4))){
					new_entry = true;
					//update altBetterCount based on altPred's relation to primePred and the actual direction
					if (primePred != altPred){
						if (altPred == actualDir){
							// Alternate prediction more useful is a counter to be of 4 bits
							if (altBetterCount < 15){  
								altBetterCount++;
							}    
						} else if (altBetterCount > 0){
								altBetterCount--;
						}
					}
				}
			}
		
		
			//Proceeding to allocation of the entry.
			//This only occurs if:
			//1. overall prediction is incorrec t
			//2. primeTable is not the table using the longest history
			// allocate entry on a table using a longer history
			if((!new_entry) || (new_entry && (primePred != actualDir))) {    //this type of allocation should not occur on correct new entries
				if (((predDir != actualDir) & (primeTable > 0))) {			
					for (int i = 0; i < primeTable; i++){
						if (tageTables[i][index[i]].usefulBits == 0);
								strong_old_present = true; //identified candidate for replacement
				  	}
					//if no suitable candidates found, decrement the usefulness of the corresponding entry in all tables using longer history
					if (strong_old_present == false){
						for (int i = primeTable - 1; i >= 0; i--){
							tageTables[i][index[i]].usefulBits--;
						}
					} else { //pick which table to use with allocation
						int count = 0;
						int bank_store[NUMTAGETABLES - 1] = {-1, -1, -1}; //store which tables have available slot
						int matchBank = 0;
						for (int i = 0; i < primeTable; i++) {
							if (tageTables[i][index[i]].usefulBits == 0){
								count++;
								bank_store[i] = i;
							}
						}  
						if(count == 1) { // only one table has available slot
							matchBank = bank_store[0];
						} else if(count > 1){ //allocate to the one with longer history 
							matchBank = bank_store[(count-1)]; 
						}
						//perform actual allocation
						for (int i = matchBank; i > -1; i--) {
							if ((tageTables[i][index[i]].usefulBits == 0)) {
								if(actualDir){    
									tageTables[i][index[i]].ctr = 4;
								} else {
									tageTables[i][index[i]].ctr = 3;
								}    
								tageTables[i][index[i]].tag = tag[i];
								tageTables[i][index[i]].usefulBits = 0;
								break; // Only 1 entry allocated
							}
						}
					}
				}
			}    
		
		
			// Periodic Useful bit Reset Logic ( Important so as to optimize compared to PPM paper)
			clock++;
			//for every 256 K instruction 1st MSB than LSB
			if(clock == (256*1024)) {
				// reset clock
				clock = 0;
				if(clock_flip == 1) {
					// this is the 1st time
					clock_flip = 0;
				} else {
					clock_flip = 1;
				}
				if(clock_flip == 1){// MSB turn
					for (int jj = 0; jj < NUMTAGETABLES; jj++)
					{    
						for (UINT32 ii = 0; ii < numTageTableEntries; ii++)
						{
							tageTables[jj][ii].usefulBits = tageTables[jj][ii].usefulBits & 1;
						}
					}
				} else {// LSB turn
					for (int jj = 0; jj < NUMTAGETABLES; jj++) {    
						for (UINT32 ii = 0; ii < numTageTableEntries; ii++) {
							tageTables[jj][ii].usefulBits = tageTables[jj][ii].usefulBits & 2;
						}
					}
				}
			}
		
			
			// update the GHR
			GHR = (GHR << 1);
			
			if(actualDir == TAKEN){
				GHR.set(0,1); 
			}
			//update compressed index and tag
			for (int i = 0; i < NUMTAGETABLES; i++) {			
					compressedIndex[i].updateCompHist(GHR);
					compressedTag[0][i].updateCompHist(GHR);
					compressedTag[1][i].updateCompHist(GHR);
					
			}
		  	// update PHR
			PHR = (PHR << 1); 
			if(pc & 1) {
				PHR = PHR + 1;
			}
			PHR = (PHR & ((1 << 16) - 1)); //ensure PHR is 16 bits

		}
	}
};
