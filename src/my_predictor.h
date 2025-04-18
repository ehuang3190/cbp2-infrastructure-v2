// my_predictor.h
// This file contains a sample my_predictor class.
// It is a simple 32,768-entry gshare with a history length of 15.
// Note that this predictor doesn't use the whole 32 kilobytes available
// for the CBP-2 contest; it is just an example.

#include "tage_predictor.h"

class my_update : public branch_update {
public:
	unsigned int index;
};

class my_predictor : public branch_predictor {
public:

	TAGE_PREDICTOR tage;
	my_update u;
	branch_info bi;


	my_predictor (void) { 
		// Initialize the predictor
		tage = TAGE_PREDICTOR();
		//cout << "TAGE predictor initialized" << endl;
	}

	branch_update *predict (branch_info & b) {
		bi = b;
		if (b.br_flags & BR_CONDITIONAL) {
			u.direction_prediction (tage.GetPrediction (b.address));
		} else {
			u.direction_prediction (true);
		}
		u.target_prediction (0);
		return &u;
	}

	void update (branch_update *u, bool taken, unsigned int target) {
		if (bi.br_flags & BR_CONDITIONAL) {
			tage.UpdatePredictor (bi.address, taken, u->direction_prediction(), target);
		}
	}
};
