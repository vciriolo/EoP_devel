#ifndef __eleIDMap__
#define __eleIDMap__

#include<iostream>
#include<map>

class eleIDMap
{

public:
	std::map<std::string, UInt_t> eleIDmap;

	eleIDMap()
	{

		eleIDmap["fiducial"]          = 0x0001;

		/* eleIDmap["loose"]             = 0x0002; */ // removed since 29/11/2017
		/* eleIDmap["medium"]            = 0x0004; */ // removed since 29/11/2017
		/* eleIDmap["tight"]             = 0x0008; */ // removed since 29/11/2017
		/* eleIDmap["WP90PU"]            = 0x0010; */ // removed since 29/11/2017
		/* eleIDmap["WP80PU"]            = 0x0020; */ // removed since 29/11/2017
		/* eleIDmap["WP70PU"]            = 0x0040; */ // removed since 29/11/2017
		/* eleIDmap["loose25nsRun2"]     = 0x0080; */ // removed since 29/11/2017
		/* eleIDmap["medium25nsRun2"]    = 0x0100; */ // removed since 29/11/2017
		/* eleIDmap["tight25nsRun2"]     = 0x0200; */ // removed since 29/11/2017
		/* eleIDmap["loose50nsRun2"]     = 0x0400; */ // removed since 29/11/2017
		/* eleIDmap["medium50nsRun2"]    = 0x0800; */ // removed since 29/11/2017
		/* eleIDmap["tight50nsRun2"]     = 0x1000; */ // removed since 29/11/2017

		// 0T ids
		eleIDmap["medium25nsRun2Boff"] = 0x2000;
		eleIDmap["diphoton25nsRun2Boff"] = 0x4000;
		eleIDmap["diphotonIso25nsRun2Boff"] = 0x8000;

		//official eleIDs
		eleIDmap["cutBasedElectronID-Spring15-25ns-V1-standalone-veto"]    =  0x10000;
		eleIDmap["cutBasedElectronID-Spring15-25ns-V1-standalone-loose"]   =  0x20000;
		eleIDmap["cutBasedElectronID-Spring15-25ns-V1-standalone-medium"]  =  0x40000;
		eleIDmap["cutBasedElectronID-Spring15-25ns-V1-standalone-tight"]   =  0x80000;
		eleIDmap["cutBasedElectronID-Spring15-50ns-V1-standalone-veto"]    =  0x100000;
		eleIDmap["cutBasedElectronID-Spring15-50ns-V1-standalone-loose"]   =  0x200000;
		eleIDmap["cutBasedElectronID-Spring15-50ns-V1-standalone-medium"]  =  0x400000;
		eleIDmap["cutBasedElectronID-Spring15-50ns-V1-standalone-tight"]   =  0x800000;

		eleIDmap["loose25nsRun2V2016"]                                     =  0x2000000;
		eleIDmap["veto25nsRun22016Moriond"]                                =  0x4000000;
		eleIDmap["loose25nsRun22016Moriond"]                               =  0x8000000;
		eleIDmap["medium5nsRun22016Moriond"]                               =  0x10000000;
		eleIDmap["tight25nsRun22016Moriond"]                               =  0x20000000;

	}

};

#endif

