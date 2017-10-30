#include "../interface/runDivide_class.hh"
#include "TTreeFormula.h"
#define DEBUG
#include <cassert>
#ifdef DEBUG
#include <TStopwatch.h>
#endif

#define MAX_TIME_DIFF 3600

runDivide_class::runDivide_class(void)
{
}

runDivide_class::runDivide_class(TCut region, std::set<TString> activeBranchList) :
	_region(region),
	_activeBranchList(activeBranchList)
{
}

runDivide_class::~runDivide_class(void)
{

}


void runDivide_class::ReadRunRangeLimits(std::string fileName)
{

	//no file provided for limits: no limits
	if(fileName == "") {
		std::cout << "[WARNING] No run range limits file provided" << std::endl;
		return;
	}

	std::ifstream file(fileName);
	if(!file.good()) {
		std::cerr << "[ERROR] File " << fileName << "to opened" << std::endl;
		return;
	}

	//filling the limits
	run_t limit_value;
	while(file.peek() != EOF && file.good()) {
		if(file.peek() == 10) { //new line
			file.get();
			continue;
		}
		while(file.peek() == 32) { // space
			file.get();
			continue;
		}
		if(file.peek() == 35) { // comment, the rest of the line can be skipped
			file.ignore(1000, 10);
			continue;
		}

		file >> limit_value;
#ifdef DEBUG
		std::cout << "[DEBUG] " << limit_value << std::endl;
#endif

		limits.insert(limit_value);
		auto keyVal = key(limit_value, 0); 
		auto itr = _runMap.find(keyVal);
		if( itr == _runMap.end()) itr = _runMap.insert(itr, std::make_pair(keyVal, line()));
		itr->second._limit=true;
	}

	return;
}


/** this method reads the ntuple and saves the number of events per run in a map */
void runDivide_class::LoadRunEventNumbers(TChain *tree, std::string runNumber_branchName, std::string lumiBlock_branchName, std::string runTime_branchName)
{
#ifdef DEBUG
	TStopwatch w;
	w.Start();
#endif
	run_t		runNumber;
	lumiBlock_t lumiBlock;
	time_t		runTime;

	tree->SetBranchStatus("*", 0);
	tree->GetBranch(runNumber_branchName.c_str())->ResetBit(kDoNotProcess);
	tree->GetBranch(lumiBlock_branchName.c_str())->ResetBit(kDoNotProcess);
	tree->GetBranch(runTime_branchName.c_str())->ResetBit(kDoNotProcess);
//	tree->SetBranchStatus(runNumber_branchName, 1);
//	tree->SetBranchStatus(runTime_branchName, 1);

	tree->SetBranchAddress(runNumber_branchName.c_str(), &runNumber);
	tree->SetBranchAddress(lumiBlock_branchName.c_str(), &lumiBlock);
	tree->SetBranchAddress(runTime_branchName.c_str(), &runTime);

	//loop over tree and count the events per run
	tree->GetEntries();

	for(std::set<TString>::const_iterator itr = _activeBranchList.begin();
	        itr != _activeBranchList.end();
	        itr++) {
		std::cout << "[STATUS] Enabling branch: " << *itr << std::endl;
		tree->SetBranchStatus(*itr, 1);
	}

	TTreeFormula * selector = (_region != "") ? new TTreeFormula("region", _region, tree) : NULL;
	std::cout << "[STATUS] selecting events in region " << _region << std::endl;

	Long64_t treenumber = -1;
	for(Long64_t ientry = 0; ientry < tree->GetEntriesFast(); ++ientry) {
		tree->GetEntry(ientry);
		if(tree->GetTreeNumber() != treenumber) {
			treenumber = tree->GetTreeNumber();
			selector->UpdateFormulaLeaves();
		}
		if(selector != NULL && selector->EvalInstance() == true) {
			auto keyVal = key(runNumber, lumiBlock);
			auto itr = _runMap.find(keyVal);
			if( itr == _runMap.end()) {
				_runMap.insert(itr, std::make_pair(keyVal, line(runTime, runNumber, lumiBlock)));
			} else {
				itr->second.update(runTime);
			}
		}
	}
	tree->ResetBranchAddresses();
#ifdef DEBUG
	w.Stop();
	w.Print();
#endif
	return;
}



// read the run range output file 
void runDivide_class::FillRunLimits(std::string fileName)
{
	//no file provided for limits: no limits
	if(fileName == "") {
		std::cout << "[WARNING] No run range file provided" << std::endl;
		return;
	}

	std::ifstream file(fileName);
	if(!file.good()) {
		std::cerr << "[ERROR] File " << fileName << "to opened" << std::endl;
		return;
	}

	//filling the limits
	run_t		run_min,  run_max;
	lumiBlock_t	lumi_min, lumi_max;
	time_t		time_min, time_max;
	Long64_t    nEvents;
	while(file.peek() != EOF && file.good()) {
		if(file.peek() == 10) { //new line
			file.get();
			continue;
		}
		if(file.peek() == 35) { // comment, the rest of the line can be skipped
			file.ignore(1000, 10);
			continue;
		}

		file >> run_min>> lumi_min>> time_min>> run_max>> lumi_max>> time_max >> nEvents;
		
		auto keyVal = key(run_min, lumi_min); 
		_runMap.insert( _runMap.end(), 
						std::make_pair(keyVal, line(time_min, run_max, lumi_max, time_max, nEvents)));

	}

	return;
}


// loop over the map containing the number of events per run-lumi and divide it into bins
void runDivide_class::FillRunLimits(unsigned int nEvents_min, float nEventsFrac_min)
{

#ifdef DEBUG
	std::cout << "[DEBUG] nEvents_min: " << nEvents_min << std::endl;
#endif
	
	// look over the map and start merging only LSs in the same run
    // start from the second, since the first has already been added
	// the current implementation is not optimal: the last LS can have  twice the required events
	
	auto key_itr = _runMap.begin(); 
	while(key_itr != _runMap.end()) { 
		auto key_next_itr = key_itr;
		key_next_itr++;

		if( // do not add events if already reached the thresold
			(key_itr->second._nEvents < nEvents_min) &&
			run(key_next_itr->first) == run(key_itr->first) && // do not merge with following run
            //add events if with the new LS the total number of events is closer to the threshold
			abs(key_next_itr->second._nEvents + key_itr->second._nEvents - nEvents_min) < abs(key_next_itr->second._nEvents  - nEvents_min)
			)
		{
			key_itr->second.add(key_next_itr->second); // add to the last
			_runMap.erase(key_next_itr);
		} else key_itr++;
	}

	key_itr = _runMap.begin();
	while(key_itr != _runMap.end()) { 
		auto key_next_itr = key_itr;
		key_next_itr++;
		
		if( key_next_itr!= _runMap.end() && 
			! key_next_itr->second._limit && 

            // do not add events from runs with time gap > 1h
			key_next_itr->second._timeMin -  key_itr->second._timeMax < MAX_TIME_DIFF
			&& 
            //add events if with the new run the total number of events is closer to the threshold
			( abs(key_next_itr->second._nEvents + key_itr->second._nEvents - nEvents_min) < abs(key_itr->second._nEvents  - nEvents_min) 
			  || abs(key_next_itr->second._nEvents + key_itr->second._nEvents - nEvents_min) < abs(key_next_itr->second._nEvents  - nEvents_min) ) 
			)
		{
			key_itr->second.add(key_next_itr->second); // add to the last
			_runMap.erase(key_next_itr);
		}else key_itr++;
	}

	//now remove empty runs
	key_itr = _runMap.begin();
	while(key_itr != _runMap.end()) { 
		if(key_itr->second._nEvents==0){
			key_itr = _runMap.erase(key_itr);
		}else key_itr++;
	}

	


	return;
}



//==============================
void runDivide_class::Divide(TChain *tree, std::string fileName,
                             unsigned int nEvents_min,
                             std::string runNumber_branchName,
							 std::string lumiBlock_branchName,
                             std::string time_branchName)
{

#ifdef DEBUG
	std::cout << "[DEBUG] Reading run range limits from file: " << fileName << std::endl;
#endif
	ReadRunRangeLimits(fileName);


#ifdef DEBUG
	std::cout << "[DEBUG] Loading events per run from tree" << std::endl;
#endif
	LoadRunEventNumbers(tree, runNumber_branchName, lumiBlock_branchName, time_branchName);

#ifdef DEBUG
	std::cout << "[DEBUG] Filling run limits" << std::endl;
#endif
	FillRunLimits(nEvents_min);
	return;

}



std::ostream& operator<<(std::ostream& os, runDivide_class& r)
{
	for(auto runM = r._runMap.begin(); runM != r._runMap.end(); runM++) {
		auto runN = runM;
		runN++;

		if(runN != r._runMap.end()) {
			os << r.printline(runM, runN) << std::endl;
		}else{
			os << r.printline(runM, runM) << std::endl;
		}
	}
	return os;
}


std::vector<TString> runDivide_class::GetRunRanges(void){
		std::vector<TString> v;

		for(auto key_itr = _runMap.begin(); key_itr != _runMap.end(); key_itr++) {
			auto key_next_itr = key_itr;
			key_next_itr++;
			if(key_next_itr == _runMap.end()) key_next_itr=key_itr;

			char range[60];
			sprintf(range, "%u-%u",
					run(key_itr->first), 
					run(key_next_itr->second._keyMax)
					);
			v.push_back(range);		
		}
		return v;
}

std::vector<std::string> runDivide_class::GetRunRangesSelectionString(bool onlyRuns){
	
	std::set<std::string> v;
	std::vector<std::string> vs;

	for(auto key_itr = _runMap.begin(); key_itr != _runMap.end(); key_itr++) {
		auto key_next_itr = key_itr;
		key_next_itr++;
		if(key_next_itr == _runMap.end()) key_next_itr=key_itr;
		
		char range[60];
		if(onlyRuns){
			sprintf(range, "runNumber_%u_%u",
					run(key_itr->first), 
					run(key_itr->second._keyMax)
				);
		}else{
			sprintf(range, "runNumber_%u_%u-lumiBlock_%u_%u",
					run(key_itr->first),
					run(key_next_itr->second._keyMax),
					lumi(key_itr->first),
					lumi(key_next_itr->second._keyMax)
				);
		}			
			
		v.insert(range);
	}
	for(auto vi : v) vs.push_back(vi); // this is to remove duplicated ranges
	return vs;
}
