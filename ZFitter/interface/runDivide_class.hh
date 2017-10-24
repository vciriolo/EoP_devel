#ifndef runDivide_class_hh
#define runDivide_class_hh

#include <TChain.h>
#include <TCut.h>

#include <iostream>
#include <fstream>
#include <ostream>

#include <map>

/** class runDivide_class
 * \author Shervin Nourbakhsh
 * \brief Define run-range bins containing a fix number of events
 *
 * \details
 *  - take a list of runs (limits) that have to be the first run of a run range
 *  - For each run:
 *    - counts the number of events
 *    - keep in memory the minTime and maxTime of the events in that run
 *  - agreegate runs respecting the limits
 *
 */
class runDivide_class
{
	// here define the types
	typedef ULong64_t key_t;
	typedef UInt_t run_t;
	typedef UInt_t time_t;
	typedef UShort_t lumiBlock_t;

public:
	runDivide_class(void);
	runDivide_class(TCut region, std::set<TString> activeBranchList);
	~runDivide_class(void);

	/*==================== define a sub-class to better handle some operations*/
	class line
	{
	public:
		Long64_t		_nEvents;
		time_t			_timeMin;
		time_t			_timeMax;
		key_t           _keyMax;
		// run_t			_runMax;
		// lumiBlock_t		_lumiMax;
		bool            _limit;
		line():
			_nEvents(0), _timeMin(99999999), _timeMax(0), _keyMax(0),  _limit(false)
		{
		}

		line(time_t t, run_t run, lumiBlock_t lumi):
			_nEvents(1), _timeMin(t), _timeMax(t), _limit(false)
		{
			_keyMax=runDivide_class::key(run, lumi);
		};

		/** increment the counters, to be called for each event with that runNumber
		 * param t Event time
		 * this method udpates the min and max time of the run and increment the number of events by 1
		 */
		void update(time_t t)   // increment
		{
			_nEvents++;
			_timeMin = std::min(_timeMin, t);
			_timeMax = std::max(_timeMax, t);
		}

		/** this method merges this run with the one provided, the number of events are summed, the timeMin and timeMax are updated to take into accont the values of the new run*/
		void add(line& l)
		{
			_nEvents += l._nEvents;
			_timeMin = std::min(_timeMin, l._timeMin);
			_timeMax = std::max(_timeMax, l._timeMax);
			_keyMax  = std::max(_keyMax,  l._keyMax);
		};

	};
	
	std::string printline(const key_t& k, const line& l)
		{
			char range[60];
			sprintf(range, "%u-%u\t%llu\t%u-%u\t%u-%u", run(k), run(l._keyMax), l._nEvents, l._timeMin, l._timeMax, lumi(k), lumi(l._keyMax));
			return std::string(range);
		};
	
	// here define the maps storing the needed infos for manipulation
	std::map<key_t, line>  _runMap;
//	std::map<key_t, line>  _runRangeMap;

	
public:
	// this is a friendly way to print the output
	friend std::ostream& operator<<(std::ostream& os, runDivide_class& r);

	void Divide(TChain *tree,
	            std::string fileName, // limits
	            unsigned int nEvents_min = 50000,
	            std::string runNumber_branchName = "runNumber",
				std::string lumiBlock_branchName = "lumiBlock",
	            std::string time_branchName      = "eventTime");

	void PrintRunRangeEvents(std::ostream& os = std::cout)
	{
		os << *this;
		return;
	}

private:

	TCut _region;
	std::set<TString> _activeBranchList;
	std::set<run_t> limits; // set of limits given by source file

	void ReadRunRangeLimits(std::string fileName); ///< read the file with the run limits: run ranges cannot go across the limits
	void LoadRunEventNumbers(TChain *tree,
	                         std::string runNumber_branchName = "runNumber",
	                         std::string lumiBlock_branchName = "lumiBlock",
	                         std::string runTime_branchName = "eventTime"); ///< read from the tree the list of run numbers, the minimum and maximum time of the events of each run and the number of events in each run. This fills all the maps

	void FillRunLimits(unsigned int nEvents_min = 100000, float nEventsFrac_min = 0.7);


	// conversions
	inline run_t run(key_t key)
	{
		return (run_t)(key >> 32);
	}
	inline lumiBlock_t lumi(key_t key)
	{
		return (lumiBlock_t)(key);
	}
	static inline key_t key(run_t run, lumiBlock_t lumi)
	{
		return ((key_t)(run) << 32) | ((key_t)(lumi));
	}

};


#endif
