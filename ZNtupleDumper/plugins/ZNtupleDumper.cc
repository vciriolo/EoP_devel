// -*- C++ -*-
//
// Package:    ZNtupleDumper
// Class:      ZNtupleDumper
//
/// Zee and E/p ntuple dumper from patElectrons
/**\class ZNtupleDumper ZNtupleDumper.cc Calibration/ZNtupleDumper/src/ZNtupleDumper.cc
 * \author Shervin Nourbakhsh
 * Description: Zee and E/p ntuple dumper from patElectrons
 *
 * \todo
 - recHitCollection is included in the PAT electron collection, take it from there
 - take the R9 from the PAT electron: electron->r9()
 - flag for Zee or Wenu dump
 - Use MET corrections
   - https://twiki.cern.ch/twiki/bin/view/CMS/MissingET
   - https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis
   - https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMETRecipe53X
 *
 */

//
// Original Author:  Shervin Nourbakhsh,40 1-B24,+41227671643,
//         Created:  Mon Jan 16 18:20:34 CET 2012
// $Id$
//
//
//
// system include files
#include <memory>
#include <iostream>
#include <string.h>
#include <regex>
#include <map>

// root include files
#include <TTree.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TString.h>
#include <TFile.h>

// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDGetToken.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtraFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//#include "Calibration/ZNtupleDumper/interface/puWeights_class.hh"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
//#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
//#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"


#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatios.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatiosRef.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAlphas.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"

#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsMCRcd.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRefRcd.h"
#include "CondFormats/DataRecord/interface/EcalLaserAlphasRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"

#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

// HLT trigger
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include <FWCore/Common/interface/TriggerNames.h>
#include <DataFormats/Common/interface/TriggerResults.h>

// alcaSkimPaths
#include "DataFormats/Provenance/interface/ParameterSetID.h"

#include "Calibration/ZNtupleDumper/interface/eleIDMap.h"
#include "DataFormats/Common/interface/ValueMap.h"

// number of electrons in each branch (max nEle)
#define NELE 3
#define initSingleFloat     -999.
#define initSingleInt          0
#define initSingleIntCharge -100
#define initFloat     { initSingleFloat, initSingleFloat, initSingleFloat }
#define initInt       { initSingleInt, initSingleInt, initSingleInt }
#define initIntCharge { initSingleIntCharge, initSingleIntCharge, initSingleIntCharge }
//#define DEBUG


class ZNtupleDumper : public edm::EDAnalyzer
{
public:
	explicit ZNtupleDumper(const edm::ParameterSet&);
	~ZNtupleDumper();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
	virtual void beginJob() ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;

	virtual void beginRun(edm::Run const&, edm::EventSetup const&);
	virtual void endRun(edm::Run const&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

	// ----------member data ---------------------------
private:
	const edm::Event *_pEvent;      // to make the event content available in private methods
	const edm::EventSetup *_pSetup; // to make the event setup available in private methods

	bool _isMC;
	bool _isPartGun;
	double _doHighEta_LowerEtaCut;

	//Handles and inputTags
private:
	//------------------------------ Handles
	edm::ESHandle<CaloTopology>  _topologyHandle;
	edm::ESHandle<EcalPedestals> _ped;

	//--------------- for main _tree
	edm::Handle<std::vector<pat::Electron> >       _electronsHandle;
	edm::Handle<std::vector<pat::Muon> >           _muonsHandle;
	edm::Handle<std::vector<pat::Photon> >         _photonsHandle;
	edm::Handle<std::vector< reco::GenParticle> >  _genParticlesHandle;
	edm::Handle<std::vector<reco::SuperCluster>>   _EESuperClustersHandle; //used only for high-eta
	edm::Handle<reco::VertexCollection>            _primaryVertexHandle; // for _nPV
	edm::Handle<double>                            _rhoHandle;
	edm::Handle<std::vector< PileupSummaryInfo > > _PupInfo;
	edm::Handle< GenEventInfoProduct >             _GenEventInfoHandle;
	edm::Handle< reco::PFMETCollection >           _metHandle;
	edm::Handle<edm::TriggerResults>               _triggerResultsHandle;
	edm::Handle<edm::TriggerResults>               _WZSkimResultsHandle;
	edm::Handle<EcalRecHitCollection>              _ESRechitsHandle;

	//--------------- for eleIDtree
	edm::Handle<reco::BeamSpot> _bsHandle;
	edm::Handle<reco::ConversionCollection> _conversionsHandle;

	//--------------- for _extraCalibTree
	edm::Handle<EcalUncalibratedRecHitCollection> _pEBUncRecHits;
	edm::Handle<EcalUncalibratedRecHitCollection> _pEEUncRecHits;

	//------------------------------ Input Tags
	// input tag for primary vertex
	edm::EDGetTokenT<GenEventInfoProduct>              _generatorInfoToken;
	edm::EDGetTokenT<reco::VertexCollection>           _vtxCollectionToken;
	edm::EDGetTokenT<reco::BeamSpot>                   _beamSpotToken;
	edm::EDGetTokenT<pat::ElectronCollection>          _electronsToken;
	edm::EDGetTokenT<pat::MuonCollection>              _muonsToken;
	edm::EDGetTokenT<pat::PhotonCollection>            _photonsToken;
	edm::EDGetTokenT<std::vector<reco::GenParticle> >  _genParticlesToken;

	edm::EDGetTokenT<EcalRecHitCollection>             _recHitCollectionEBToken;
	edm::EDGetTokenT<EcalRecHitCollection>             _recHitCollectionEEToken;
	edm::EDGetTokenT<EcalRecHitCollection>             _recHitCollectionESToken;
	edm::EDGetTokenT<EcalUncalibratedRecHitCollection> _ebURHToken;
	edm::EDGetTokenT<EcalUncalibratedRecHitCollection> _eeURHToken;

	edm::EDGetTokenT<std::vector<reco::SuperCluster> > _EESuperClustersToken;

	// input _rho
	edm::EDGetTokenT<double>                          _rhoToken;
	edm::EDGetTokenT<std::vector<PileupSummaryInfo> > _pileupInfoToken;
	edm::EDGetTokenT<reco::ConversionCollection>      _conversionsProducerToken;
	edm::EDGetTokenT<reco::PFMETCollection>           _metToken;
	edm::EDGetTokenT<edm::TriggerResults>             _triggerResultsToken;
	edm::EDGetTokenT<edm::TriggerResults>             _WZSkimResultsToken;
	edm::InputTag triggerResultsTAG,                  _WZSkimResultsTAG;
	std::vector< std::string> hltPaths,               _SelectEvents;
private:
	std::string _ntupleFileName;

	bool _doExtraCalibTree; ///< bool to activate the dump of the extra calib _tree for E/p ICs
	bool _doExtraStudyTree; ///< bool to activate the dump of the extra _tree for study with values for single recHits
	bool _doEleIDTree;      ///< bool to activate the dump of the electronID variables in a separate _tree
	bool _doTrackTree;      ///< bool to activate the dump of the track _tree for studies with electron tracks

	TTree * _tree;                   //< output file for standard ntuple

	edm::Timestamp _eventTimeStamp;

	// ntuple members, private to make them visible in doxygen
private:
	/**
		\addtogroup NTUPLESTRUCTURE
		@{
	 */
	UInt_t    _runNumber;     ///< run number
	UShort_t  _lumiBlock;     ///< lumi section
	Long64_t  _eventNumber;   ///< event number
	UInt_t    _eventTime;     ///< unix time of the event
	UShort_t  _nBX;           ///< bunch crossing

	Float_t   _mcGenWeight;   ///< weight in generator for MC

	std::vector< std::string >  _HLTNames[1];   ///< List of HLT names
	std::vector<Bool_t>         _HLTResults[1]; ///< 0 = fail, 1=fire
	std::map<std::string, bool> _HLTBits;
	Bool_t                      _HLTfire;       ///< true if pass the triggers indicated by hltPaths in cfg

	//pileup
	Float_t  _rho;    ///< _rho fast jet
	UChar_t  _nPV;    ///< nVtx
	UChar_t  _nPU;    ///< number of PU (filled only for MC)

	// primary vertex
	Float_t         _vtxX;           ///< primary vertex x coordinate
	Float_t         _vtxY;           ///< primary vertex y coordinate
	Float_t         _vtxZ;           ///< primary vertex z coordinate

	// selection
	UInt_t _eleID[NELE] = initInt;      ///< bit mask for _eleID: 1=fiducial, 2=loose, 6=medium, 14=tight, 16=WP90PU, 48=WP80PU, 112=WP70PU, 128=loose25nsRun2, 384=medium25nsRun2, 896=tight25nsRun2, 1024=loose50nsRun2, 3072=medium50nsRun2, 7168=tight50nsRun2. Selection from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification#Electron_ID_Working_Points

	Short_t  _chargeEle[NELE]    = initIntCharge; ///< -100: no electron, 0: SC or photon, -1 or +1:electron or muon //Char_t is interpreted as char and not as integer
	UChar_t  _recoFlagsEle[NELE] = initInt;       ///< 1=trackerDriven, 2=ecalDriven only, 3=tracker and ecal driven
	Float_t  _etaEle[NELE]       = initFloat;
	Float_t  _phiEle[NELE]       = initFloat;     ///< phi of the electron (electron object)
	Float_t  _R9Ele[NELE]        = initFloat;     ///< e3x3/_rawEnergySCEle
	Float_t  _fbremEle[NELE]     = initFloat;

	Bool_t _isMustacheSC[NELE] = initInt;
	Float_t _etaSCEle[NELE]    = initFloat;
	Float_t _phiSCEle[NELE]    = initFloat; ///< phi of the SC

	// seed of the SC
	Short_t _xSeedSC[NELE]                    = initInt;    ///< ieta(ix) of the SC seed in EB(EE)
	Short_t _ySeedSC[NELE]                    = initInt;    ///< iphi(iy) of the SC seed in EB(EE)
  UChar_t _gainSeedSC[NELE]                       = {9,9,9};    ///< Gain switch 0==gain12, 1==gain6, 2==gain1; gain status of the seed of the SC
	Float_t _energySeedSC[NELE]               = initFloat;  ///< energy of the rechit seeding the SC
	Float_t _timeSeedSC[NELE]                 = initFloat;  ///< time of the rechit seeding the SC
	Float_t _icSeedSC[NELE]                   = initFloat;  ///< inter-calibration coefficient of the SC seed crystal
	Float_t _laserSeedSC[NELE]                = initFloat;  ///< laser correction of the SC seed crystal
	Float_t _avgLCSC[NELE]                    = initFloat;  ///< average laser correction for the SC
	Float_t _alphaSeedSC[NELE]                = initFloat;  ///< alpha of the seed
	Float_t _pedestalSeedSC[NELE]             = initFloat;  ///< pedestal of the seed of the SC
	Float_t _noiseSeedSC[NELE]                = initFloat;  ///< noise of the seed of the SC

	Float_t _amplitudeSeedSC[NELE]            = initFloat;  ///< DetId amplitude with the second highest energy in the SC
	// second to seed, i.e. the second highest energy crystal of the SC
	Float_t _energySecondToSeedSC[NELE]       = initFloat;  ///< energy of the rechit with the second highest energy in the SC
	Float_t _timeSecondToSeedSC[NELE]         = initFloat;  ///< time of the rechit with the second highest energy in the SC
	Float_t _amplitudeSecondToSeedSC[NELE]    = initFloat;  ///< time of the rechit with the second highest energy in the SC

	Float_t _energyEle[NELE]                  = initFloat;  ///< electron.energy(), not changed by rereco
	Float_t _rawEnergySCEle[NELE]             = initFloat;  ///< SC energy without cluster corrections
	Float_t _energy_ECAL_ele[NELE]            = initFloat;  ///< ele-tuned regression energy: mustache for rereco and correctedEcalEnergy for official reco
	Float_t _energy_ECAL_pho[NELE]            = initFloat;  ///< pho-tuned regression energy: mustache for rereco and correctedEcalEnergy for official reco
	Float_t _energyUncertainty_ECAL_ele[NELE] = initFloat;  ///< ele-tuned regression energy: mustache for rereco and correctedEcalEnergy for official reco
	Float_t _energyUncertainty_ECAL_pho[NELE] = initFloat;  ///< pho-tuned regression energy: mustache for rereco and correctedEcalEnergy for official reco

	Float_t _esEnergySCEle[NELE]              = initFloat;  ///< pre-shower energy associated to the electron
	Float_t _esEnergyPlane1SCEle[NELE]        = initFloat;  ///< energy associate to the electron in the first plane of ES
	Float_t _esEnergyPlane2SCEle[NELE]        = initFloat;  ///< energy associate to the electron in the second plane of ES
	Float_t _rawESEnergyPlane1SCEle[NELE]     = initFloat;  ///< pre-shower rechit energy sum of Plane 1 associated to the electron
	Float_t _rawESEnergyPlane2SCEle[NELE]     = initFloat;  ///< pre-shower recHit energy sum of Plane 2 associated to the electron

	Float_t _energy_3x3SC[NELE]               = initFloat;  //< sum of the recHit energy in 3x3 matrix centered at the seed of the SC
	Float_t _energy_5x5SC[NELE]               = initFloat;  ///< sum of the recHit energy in 5x5 matrix centered at the seed of the SC
	Float_t _eBCseedEle[NELE]                 = initFloat;  ///< energy of the basic cluster seeding the SC
	Float_t _pModeGsfEle[NELE]                = initFloat;  ///< track momentum from Gsf Track (mode)
	Float_t _pAtVtxGsfEle[NELE]               = initFloat;  ///< momentum estimated at the vertex
	Float_t _trackMomentumErrorEle[NELE]      = initFloat;  ///< track momentum error from standard electron method
	Float_t _pNormalizedChi2Ele[NELE]         = initFloat;  ///< track normalized chi2 of the fit (GSF)
	Float_t _gsfTrackLengthFromVtxP[NELE]     = initFloat;  ///< track length computed from the vertex position and momentum to the extrapolated impact on the calorimeter
	Float_t _gsfTrackLengthFromTangents[NELE] = initFloat;  ///< track length computed from the vertex position to the extrapolated impact on the calorimeter, using at each tracker layer its mometum estimate (and position)

	Float_t _invMass;
	Float_t _invMass_rawSC;
	Float_t _invMass_rawSC_esSC;
	Float_t _invMass_ECAL_ele;   ///< invariant mass using ECAL energy, this is mustache ele-tuned regression if rereco, and correctedEcalEnergy if official reco
	Float_t _invMass_ECAL_pho;   ///< invariant mass using ECAL energy, this is mustache pho-tuned regression if rereco, and correctedEcalEnergy if official reco
	//   Float_t invMass_e3x3;
	Float_t _invMass_5x5SC;
	Float_t _invMass_highEta;   ///< invariant mass using same energy as invMass_ECAL_ele for the first electron and the e5x5 for the electron in the high eta region (|eta|>2.5)

	Float_t _invMass_mumu;

	Float_t  _etaMCEle[NELE]    = initFloat;
	Float_t  _phiMCEle[NELE]    = initFloat;
	Float_t  _energyMCEle[NELE] = initFloat;  ///< Electron MC true energy
	Float_t  _invMass_MC;
	Bool_t	  _ZEvent = false;


	//============================== ExtraCalibTree (for E/p)
	TFile *_extraCalibTreeFile;
	TTree *_extraCalibTree;
	Int_t _nRecHitsEle[NELE] = initInt;
	Int_t _nHitsSCEle[NELE] = initInt;
	std::vector<unsigned int>  _rawIdRecHitSCEle[NELE];
	std::vector<int>           _XRecHitSCEle[NELE];
	std::vector<int>           _YRecHitSCEle[NELE];
	std::vector<int>           _ZRecHitSCEle[NELE];
	std::vector<float>         _energyRecHitSCEle[NELE];
	std::vector<float>         _fracRecHitSCEle[NELE];
	std::vector<int>           _recoFlagRecHitSCEle[NELE];
	//==============================

	//============================== ExtraStudyTree
	TFile *_extraStudyTreeFile;
	TTree *_extraStudyTree;
	std::vector<float>               _LCRecHitSCEle[NELE];
	std::vector<float>               _AlphaRecHitSCEle[NELE];
	std::vector<float>               _ICRecHitSCEle[NELE];
	std::vector<std::vector<float> > _ootAmplisUncalibRecHitSCEle[NELE];
	std::vector<float>               _ampliUncalibRecHitSCEle[NELE];
	std::vector<float>               _ampliErrUncalibRecHitSCEle[NELE];
	std::vector<float>               _pedEUncalibRecHitSCEle[NELE];
	std::vector<float>               _jitterUncalibRecHitSCEle[NELE];
	std::vector<float>               _jitterErrUncalibRecHitSCEle[NELE];
	std::vector<float>               _chi2UncalibRecHitSCEle[NELE];
	std::vector<uint32_t>            _flagsUncalibRecHitSCEle[NELE];

	//============================== check ele-id and iso
	TFile *_eleIDTreeFile;
	TTree *_eleIDTree;
	Float_t _sigmaIEtaIEtaSCEle[NELE]             = initFloat;
	Float_t _sigmaIPhiIPhiSCEle[NELE]             = initFloat;
	Float_t _hOverE[NELE]                         = initFloat;
	Float_t _hOverEBC[NELE]                       = initFloat;
	Float_t _dr03TkSumPt[NELE]                    = initFloat;
	Float_t _dr03EcalRecHitSumEt[NELE]            = initFloat;
	Float_t _dr03HcalTowerSumEt[NELE]             = initFloat;
	Float_t _deltaPhiSuperClusterTrackAtVtx[NELE] = initFloat;
	Float_t _deltaEtaSuperClusterTrackAtVtx[NELE] = initFloat;
	Float_t _E1x5[NELE]                           = initFloat;
	Float_t _E1x3[NELE]                           = initFloat;
	Float_t _E2x2[NELE]                           = initFloat;
	Float_t _E2x5Max[NELE]                        = initFloat;
	Float_t _S4[NELE]                             = initFloat;
	Float_t _etaWidth[NELE]                       = initFloat;
	Float_t _phiWidth[NELE]                       = initFloat;
	Bool_t _hasMatchedConversion[NELE]            = initInt;
	Int_t _maxNumberOfExpectedMissingHits[NELE]   = initInt;
	Float_t _pfMVA[NELE]                          = initFloat;
	Float_t _eleIDloose[NELE]                     = initFloat;
	Float_t _eleIDmedium[NELE]                    = initFloat;
	Float_t _eleIDtight[NELE]                     = initFloat;
	//==============================

	//============================== TrackTree
	TFile *_trackTreeFile;
	TTree *_trackTree;
	Float_t   _gsfTrackLengthFromOuterP[NELE] = initFloat;
	Float_t   _gsfTrackLengthFromEleP[NELE] = initFloat;
	Float_t   _gsfTrackOuterPt[NELE] = initFloat;
	Float_t   _gsfTrackOuterEta[NELE] = initFloat;
	Float_t   _gsfTrackOuterPhi[NELE] = initFloat;
	Float_t   _gsfTrackVtxPt[NELE] = initFloat;
	Float_t   _gsfTrackVtxEta[NELE] = initFloat;
	Float_t   _gsfTrackVtxPhi[NELE] = initFloat;
	Float_t   _gsfTrackVtxX[NELE] = initFloat;
	Float_t   _gsfTrackVtxY[NELE] = initFloat;
	Float_t   _gsfTrackVtxZ[NELE] = initFloat;
	Float_t   _gsfTrackCaloX[NELE] = initFloat;
	Float_t   _gsfTrackCaloY[NELE] = initFloat;
	Float_t   _gsfTrackCaloZ[NELE] = initFloat;
	std::vector<float> _gsfTrackTangentsPt[NELE];
	std::vector<float> _gsfTrackTangentsEta[NELE];
	std::vector<float> _gsfTrackTangentsPhi[NELE];
	std::vector<float> _gsfTrackTangentsDeltaP[NELE];
	std::vector<float> _gsfTrackTangentsDeltaPErr[NELE];
	std::vector<float> _gsfTrackTangentsX[NELE];
	std::vector<float> _gsfTrackTangentsY[NELE];
	std::vector<float> _gsfTrackTangentsZ[NELE];

	/**
	 @}
	*/
	//==============================
private:
	TFile * _tree_file;
	void InitNewTree(void);
	void ResetMainTreeVar();

	void TreeSetSingleSCVar(const reco::SuperCluster& sc, int index);

	void TreeSetSingleElectronVar(const pat::Electron& electron1, int index);
	void TreeSetSingleElectronVar(const reco::SuperCluster& sc, int index);
	void TreeSetSinglePhotonVar(const pat::Photon& photon, int index);
	void TreeSetSingleMuonVar(const pat::Muon& muon1, int index);
	void TreeSetDiElectronVar(const pat::Electron& electron1, const pat::Electron& electron2);
	void TreeSetDiElectronVar(const pat::Electron& electron1, const reco::SuperCluster& electron2);
	void TreeSetMuMuGammaVar(const pat::Photon& photon, const pat::Muon& muon1, const pat::Muon& muon2);
	/*
	 * sort the hits belonging to the SuperCluster in descending energy and return the DetId and energy of the one ranked `rank`
	 */
	std::pair<DetId, float> findEnergySortedHit(const reco::SuperCluster& cluster, const EcalRecHitCollection * recHits, size_t rank = 0);

	void InitExtraCalibTree(void);
	void ResetExtraCalibVar();
	void TreeSetExtraCalibVar(const pat::Electron& electron1, int index);
	void TreeSetExtraCalibVar(const reco::SuperCluster& electron1, int index);
	void TreeSetExtraCalibVar(const pat::Photon& photon, int index);
	void TreeSetExtraCalibVar(const pat::Muon& muon1, int index);
	void TreeSetExtraCalibVar(const pat::Electron& electron1, const pat::Electron& electron2);
	void TreeSetExtraCalibVar(const pat::Electron& electron1, const reco::SuperCluster& electron2);
	void TreeSetExtraCalibVar(const pat::Photon& photon, const pat::Muon& muon1, const pat::Muon& muon2);
	void TreeSetExtraCalibVar(const std::vector<std::pair<DetId, float> > & hitsFracs, int index, bool isEB);

	void InitExtraStudyTree(void); // the extra study _tree uses the method of the extracalibtree
	void ResetExtraStudyVar();

	void InitEleIDTree(void);
	void ResetEleIDVar();
	void TreeSetEleIDVar(const pat::Electron& electron1, int index);
	void TreeSetEleIDVar(const pat::Electron& electron1, const pat::Electron& electron2);
	void TreeSetEleIDVar(const pat::Photon& photon, const pat::Muon& muon1, const pat::Muon& muon2);
	void TreeSetEleIDVar(const pat::Photon& photon, int index);
	void TreeSetEleIDVar(const pat::Muon& muon1, int index);

	void InitTrackTree(void);
	void ResetTrackVar();
	void TreeSetTrackVar(const pat::Electron & electron1, const pat::Electron & electron2);
	void TreeSetTrackVar(const pat::Electron & electron1, int index);

	//  void Tree_output(TString filename);
	void TreeSetEventSummaryVar(const edm::Event& iEvent);
	void TreeSetPileupVar(void);
	void TreeSetPrimaryVertex(void);
	float GetESPlaneRawEnergy(const reco::SuperCluster& sc, unsigned int planeIndex) const;

	bool elePreselection(const pat::Electron& electron) const;
	//puWeights_class puWeights;

private:

	// --------------- selection cuts
private:

	//------------------------------
	// cluster tools
	noZS::EcalClusterLazyTools *_clustertools;

	std::set<unsigned int> _alcaSkimPathIndexes;
	edm::ParameterSetID _alcaSkimPathID;
	//
	// static data member definitions
	//
	bool _presel;
	std::string _eleID_loose, _eleID_medium, _eleID_tight;

	typedef enum {
		ZEE = 0,
		WENU,
		ZSC,
		ZMMG,
		PARTGUN,
		UNKNOWN,
		SINGLEELE, //no skim, no preselection and no selection are applied
		DEBUG_T
	} eventType_t;


	eventType_t _eventType;

	bool _isMINIAOD;
	bool _hasURecHits;
};


//
// constructors and destructor
//
ZNtupleDumper::ZNtupleDumper(const edm::ParameterSet& iConfig):
	//  _isMC(iConfig.getParameter<bool>("isMC")),
	_isPartGun(iConfig.getParameter<bool>("isPartGun")),
	_doHighEta_LowerEtaCut(iConfig.getParameter<double>("doHighEta_LowerEtaCut")),
	_generatorInfoToken(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
	_vtxCollectionToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
	_beamSpotToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotCollection"))),
	_electronsToken(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronCollection"))),
	_muonsToken(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonCollection"))),
	_photonsToken(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonCollection"))),
	_genParticlesToken(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"))),
	_recHitCollectionEBToken(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>( "recHitCollectionEB" ))),
	_recHitCollectionEEToken(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>( "recHitCollectionEE" ))),
	_recHitCollectionESToken(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("recHitCollectionES"))),
	_ebURHToken(mayConsume<EcalUncalibratedRecHitCollection>(iConfig.getParameter<edm::InputTag>( "uncalibRecHitCollectionEB" ))),
	_eeURHToken(mayConsume<EcalUncalibratedRecHitCollection>(iConfig.getParameter<edm::InputTag>( "uncalibRecHitCollectionEE" ))),
	_EESuperClustersToken(consumes<reco::SuperClusterCollection>(iConfig.getParameter< edm::InputTag>("EESuperClusterCollection"))),
	_rhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastJet"))),
	_pileupInfoToken(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupInfo"))),
	_conversionsProducerToken(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversionCollection"))),
	_metToken(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metCollection"))),
	_triggerResultsToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResultsCollection"))),
	_WZSkimResultsToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("WZSkimResultsCollection"))),
	triggerResultsTAG(iConfig.getParameter<edm::InputTag>("triggerResultsCollection")),
	_WZSkimResultsTAG(iConfig.getParameter<edm::InputTag>("WZSkimResultsCollection")),
	hltPaths(iConfig.getParameter< std::vector<std::string> >("hltPaths")),
	_SelectEvents(iConfig.getParameter<std::vector<std::string> >("SelectEvents")),
	_ntupleFileName(iConfig.getParameter<std::string>("foutName")),
	_doExtraCalibTree(iConfig.getParameter<bool>("doExtraCalibTree")),
	_doExtraStudyTree(iConfig.getParameter<bool>("doExtraStudyTree")),
	_doEleIDTree(iConfig.getParameter<bool>("doEleIDTree")),
	_doTrackTree(iConfig.getParameter<bool>("doTrackTree")),
	_presel(iConfig.getParameter<bool>("useIDforPresel")),
	// used for preselection and event type determination
	_eleID_loose(iConfig.getParameter< std::string>("eleID_loose")),
	_eleID_medium(iConfig.getParameter< std::string>("eleID_medium")),
	_eleID_tight(iConfig.getParameter< std::string>("eleID_tight")),
	_eventType(ZEE),
	_isMINIAOD(!(iConfig.getParameter<edm::InputTag>("recHitCollectionEB") == edm::InputTag("alCaIsolatedElectrons", "alcaBarrelHits"))),
	_hasURecHits(!(iConfig.getParameter<edm::InputTag>("uncalibRecHitCollectionEB") == edm::InputTag("", "")))
{
	if(_isMINIAOD) std::cout << "[INFO ZntupleDumper] running on MINIAOD" << std::endl;
	//  current_dir.cd();
}


ZNtupleDumper::~ZNtupleDumper()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

	std::cout << "[STATUS] Calling the destructor" << std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void ZNtupleDumper::beginJob()
{
#ifdef DEBUG
	std::cout << "[DEBUG] Starting creation of ntuples" << std::endl;
#endif

	_tree_file = new TFile(_ntupleFileName.c_str(), "recreate");
	if(_tree_file->IsZombie()) {
		throw cms::Exception("OutputError") <<  "Output _tree not created (Zombie): " << _ntupleFileName;
		return;
	}
	_tree_file->cd();

	//now do what ever initialization is needed
	_tree = new TTree("selected", "selected");
	//_tree = fs->make<TTree>("selected","selected"); //no otherwise you have the extraCalib in the same file

	_tree->SetDirectory(_tree_file);
	// controllo su _tree==NULL

	InitNewTree();  // inizializzo il _tree dell'ntupla ridotta selezionata

	if(_doExtraCalibTree) {
		//_extraCalibTree = fs->make<TTree>("extraCalibTree","");
		// put the _extraCalibTree into the default outfile
		_extraCalibTreeFile = new TFile("extraCalibTree.root", "recreate");
		if(_extraCalibTreeFile->IsZombie()) {
			throw cms::Exception("OutputError") <<  "Output _tree for extra calib not created (Zombie): " << _ntupleFileName;
			return;
		}
		_extraCalibTreeFile->cd();

		_extraCalibTree = new TTree("extraCalibTree", "extraCalibTree");
		_extraCalibTree->SetDirectory(_extraCalibTreeFile);
		InitExtraCalibTree();
	}

	if(_doExtraStudyTree) {
		_extraStudyTreeFile = new TFile("extraStudyTree.root", "recreate");
		if(_extraStudyTreeFile->IsZombie()) {
			throw cms::Exception("OutputError") <<  "Output _tree for extra study not created (Zombie): " << _ntupleFileName;
			return;
		}
		_extraStudyTreeFile->cd();

		_extraStudyTree = new TTree("extraStudyTree", "extraStudyTree");
		_extraStudyTree->SetDirectory(_extraStudyTreeFile);
		InitExtraStudyTree();
	}


	if(_doEleIDTree) {
		_eleIDTreeFile = new TFile("eleIDTree.root", "recreate");
		if(_eleIDTreeFile->IsZombie()) {
			throw cms::Exception("OutputError") <<  "Output _tree for extra calib not created (Zombie): " << _ntupleFileName;
			return;
		}
		_eleIDTreeFile->cd();
		_eleIDTree = new TTree("eleIDTree", "eleIDTree");
		_eleIDTree->SetDirectory(_eleIDTreeFile);
		//_eleIDTree = fs->make<TTree>("eleIDTree","");
		InitEleIDTree();
	}


	if(_doTrackTree) {
		_trackTreeFile = new TFile("trackTree.root", "recreate");
		if(_trackTreeFile->IsZombie()) {
			throw cms::Exception("OutputError") <<  "Output _tree for track studies not created (Zombie): " << _ntupleFileName;
			return;
		}
		_trackTreeFile->cd();

		_trackTree = new TTree("trackTree", "trackTree");
		_trackTree->SetDirectory(_trackTreeFile);
		InitTrackTree();
	}
#ifdef DEBUG
	std::cout << "[DEBUG] End creation of ntuples" << std::endl;
#endif
}


//
// member functions
//

// ------------ method called for each event  ------------
void ZNtupleDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	ResetMainTreeVar();
	ResetExtraCalibVar();
	ResetExtraStudyVar();
	ResetEleIDVar();
	ResetTrackVar();

	iSetup.get<EcalPedestalsRcd>().get(_ped);

	//  using namespace edm;
	_eventType = _isPartGun ? PARTGUN : UNKNOWN;

	_chargeEle[0] = initSingleIntCharge;
	_chargeEle[1] = initSingleIntCharge;
	_chargeEle[2] = initSingleIntCharge;
	_invMass_mumu = 0;


	_energyEle[0] = initSingleFloat;
	_energyEle[1] = initSingleFloat;

	_pEvent = &iEvent;
	_pSetup = &iSetup;

	if( !iEvent.isRealData() ) {
		iEvent.getByToken(_pileupInfoToken, _PupInfo);
		iEvent.getByToken(_generatorInfoToken, _GenEventInfoHandle);
		_isMC = true;
	} else _isMC = false;

	if(_isMC) {
		// Gen Particles
		iEvent.getByToken(_genParticlesToken, _genParticlesHandle);
		_ZEvent = false;
		for (auto& p : *_genParticlesHandle) {
			if(p.pdgId() == 23) {
				_ZEvent = true;
				break;
			}
		}
	}

	//------------------------------ HLT
	/// \todo check why
	if(triggerResultsTAG.label() != "") iEvent.getByToken(_triggerResultsToken, _triggerResultsHandle);
	if(_WZSkimResultsTAG.label() != "") {
		iEvent.getByToken(_WZSkimResultsToken,  _WZSkimResultsHandle); //else it is not produced with ALCARECO selection
		//then the type of event has to be defined

		//Check if it is Wenu, Z or ZSC event according to triggerResults
		edm::TriggerNames alcaSkimPathNames = iEvent.triggerNames(*_WZSkimResultsHandle);

		if(!_SelectEvents.empty()) {
			// If the alca skim paths are not changing, this is done only once
			if(_alcaSkimPathID != alcaSkimPathNames.parameterSetID()) { //order of trigger results is changed
				_alcaSkimPathID = alcaSkimPathNames.parameterSetID();    //update the map of trigger index
				_alcaSkimPathIndexes.clear();
				unsigned int alcaSkimPathNameSize = alcaSkimPathNames.size(); // should have the same size of _WZSkimResultsHandle

				for(unsigned int i = 0; i < alcaSkimPathNameSize; i++) { // look over alcaSkimPaths
					std::string trgName = alcaSkimPathNames.triggerName(i);

					for(std::vector<std::string>::const_iterator selectEvents_itr = _SelectEvents.begin();
					        selectEvents_itr != _SelectEvents.end();
					        selectEvents_itr++) {
						if(std::regex_match(trgName, std::regex(*selectEvents_itr))) {
							_alcaSkimPathIndexes.insert(i);
							std::cout << " - Trigger path saved in ntuples: " << trgName << "\t" << i << std::endl;
							break;
						}

					}
					//if(_alcaSkimPathIndexes.count(i)==0){
					//std::cout << " -! Trigger path not saved in ntuples: " << trgName <<  "\t" << i << std::endl;
					//}
				}
			}

			_eventType = DEBUG_T;
			bool skipEvent = true;
			for(std::set<unsigned int>::const_iterator alcaSkimPath_itr = _alcaSkimPathIndexes.begin();
			        alcaSkimPath_itr != _alcaSkimPathIndexes.end() && skipEvent == true;
			        alcaSkimPath_itr++) {
				//std::cout << *alcaSkimPath_itr << std::endl;
				if(_WZSkimResultsHandle->accept(*alcaSkimPath_itr)) {
					skipEvent = false;
					std::string hltName_str(alcaSkimPathNames.triggerName(*alcaSkimPath_itr));
					if(hltName_str.find("WElectron") != std::string::npos)
						_eventType = WENU;
					else if(hltName_str.find("ZSCElectron") != std::string::npos)
						_eventType = ZSC;
					else if(hltName_str.find("ZElectron") != std::string::npos)
						_eventType = ZEE;
					else if(hltName_str.find("SingleElectron") != std::string::npos)
						_eventType = SINGLEELE;
					else if(hltName_str.find("Zmmg") != std::string::npos)
						_eventType = ZMMG;
					else
						_eventType = UNKNOWN;
					// this paths are exclusive, then we can skip the check of the others
					//
					//	std::cout << alcaSkimPathNames.triggerName(*alcaSkimPath_itr) << "\t" << _eventType << std::endl;
					break;
				}
			}
			//if(_alcaSkimPathIndexes.size()==0){
//			skipEvent = false;
//			_eventType = UNKNOWN;
			//}
			//std::cout << "skip event: " << skipEvent << "\t" << _eventType << std::endl;
			//assert(!skipEvent);
			if(skipEvent) return; // event not coming from any skim or paths
		}
	}
	//------------------------------ CONVERSIONS
	iEvent.getByToken(_conversionsProducerToken, _conversionsHandle);

	//------------------------------
	_clustertools = new noZS::EcalClusterLazyTools(iEvent, iSetup, _recHitCollectionEBToken,
	        _recHitCollectionEEToken);

	if(_hasURecHits) {
		iEvent.getByToken(_ebURHToken, _pEBUncRecHits);
		iEvent.getByToken(_eeURHToken, _pEEUncRecHits);
	}
	//------------------------------ electrons
	if (_eventType == ZMMG) {
		//------------------------------ muons
		iEvent.getByToken(_muonsToken, _muonsHandle);
		//------------------------------ photons
		iEvent.getByToken(_photonsToken, _photonsHandle);
	}	else {
		iEvent.getByToken(_electronsToken, _electronsHandle);
	}

	//------------------------------ SuperClusters (for high Eta studies)
	iEvent.getByToken(_EESuperClustersToken, _EESuperClustersHandle);

	// for conversions with full vertex fit
	//------------------------------  VERTEX
	iEvent.getByToken(_vtxCollectionToken, _primaryVertexHandle);
	iEvent.getByToken(_beamSpotToken, _bsHandle);
	iEvent.getByToken(_rhoToken, _rhoHandle);

	iEvent.getByToken(_metToken, _metHandle);
	iEvent.getByToken(_recHitCollectionESToken, _ESRechitsHandle);
	//if(_metHandle.isValid()==false) iEvent.getByType(_metHandle);
	reco::PFMET met = _metHandle.isValid() ? ((*_metHandle))[0] : reco::PFMET(); /// \todo use corrected phi distribution


	//Here the _HLTBits are filled. TriggerResults
	TreeSetEventSummaryVar(iEvent);
	TreeSetPileupVar(); // this can be filled once per event
	TreeSetPrimaryVertex(); // this can be filled once per event


	// at least one of the triggers
	_HLTfire = false;
	if(!hltPaths.empty()) {
		for(std::vector<std::string>::const_iterator hltPath_itr = hltPaths.begin();
		        hltPath_itr != hltPaths.end();
		        hltPath_itr++) {
			if(hltPath_itr->empty()) continue;
			std::map<std::string, bool>::const_iterator it = _HLTBits.find(*hltPath_itr);
			if(it != _HLTBits.end()) {
				_HLTfire += it->second;
				// std::cout << "Not empty:" << hltPaths[0] << "\t" << it->first << "\t" << it->second << "\t" << triggerFire << std::endl;
				//}else{
				//for(std::map<std::string, bool>::const_iterator it = _HLTBits.begin();
				//    it!=_HLTBits.end();
				//    it++){
				//    std::cout  << "\t" << it->first << "\t" << it->second << std::endl;
				//  }
				//edm::LogError("ZNtupleDumper") << "HLT path required but not find in TriggerResults" << " " << _HLTBits.size();
				//edm::LogWarning("ZNtupleDumper") << "HLT path " << *hltPath_itr << " required but not found in TriggerResults" << " " << _HLTBits.size();
				//exit(1);
			}
		}
	}
#ifdef DEBUG
	std::cout << "After trigger selection" << std::endl;
#endif

	// count electrons: needed to avoid double counting events in Wenu and Zee samples
	// in Wenu is required ONLY ONE tight electron
	// in Zee at least two loose electrons
	// in particle gun case, the matching with the gen particle is required

	int nTight   = 0; //number of electrons passing only the tight  ID for preselection
	int nMedium  = 0; //number of electrons passing only the medium ID for preselection
	int nLoose   = 0; //number of electrons passing only the loose  ID for preselection

	//if (_eventType!=ZMMG) { // count the number of electrons passing the different IDs for preselection and event type determination
	if (_eventType != UNKNOWN) { // count the number of electrons passing the different IDs for preselection and event type determination
		for( pat::ElectronCollection::const_iterator eleIter1 = _electronsHandle->begin();
		        eleIter1 != _electronsHandle->end();
		        eleIter1++) {
			if( eleIter1->electronID(_eleID_tight) )          ++nTight;
			else if( eleIter1->electronID(_eleID_medium) ) ++nMedium;
			else if( eleIter1->electronID(_eleID_loose) )  ++nLoose;
		}
	}

	bool doFill = false;

	if(_eventType == PARTGUN) {
		pat::ElectronCollection::const_iterator eleIter1 = _electronsHandle->begin();
		pat::ElectronCollection::const_iterator eleIter2 = eleIter1;
		for(eleIter1 = _electronsHandle->begin();
		        eleIter1 != _electronsHandle->end() && eleIter1->genLepton() == 0;
		        eleIter1++) {
		}
		//if no electron matching the gen particles then the event is skipped
		//if(eleIter1 == _electronsHandle->end()) return;
		if(eleIter1 == _electronsHandle->end()) eleIter1 = _electronsHandle->begin();

		//in order to not put duplicate electrons, remove the ones not matching the genparticle
		for(eleIter2 = eleIter1, eleIter2++;
		        eleIter2 != _electronsHandle->end() && eleIter2->genLepton() == 0;
		        eleIter2++) {
		}
		if(eleIter2 == _electronsHandle->end()) {
			if(eleIter1->genLepton() != 0 || _electronsHandle->size() < NELE) eleIter2 = eleIter1;
			else eleIter2 = eleIter1 + 1;
		}

		//if one electron matching the gen particles then eleIter2 = eleIter1
		//else we have two electrons
		TreeSetDiElectronVar(*eleIter1, *eleIter2);
		doFill = true;
		if(_doExtraCalibTree || _doExtraStudyTree) {
			TreeSetExtraCalibVar(*eleIter1, *eleIter2);
		}
		if(_doEleIDTree) {
			TreeSetEleIDVar(*eleIter1, *eleIter2);
		}
		// always set track vars, two of them are in the main _tree
		TreeSetTrackVar(*eleIter1, *eleIter2);
	} else if(_eventType == ZEE || _eventType == WENU || _eventType == UNKNOWN) {
#ifdef DEBUG
		std::cout << "[DEBUG] Electrons in the event: " << _electronsHandle->size() << std::endl;
#endif

		for( pat::ElectronCollection::const_iterator eleIter1 = _electronsHandle->begin();
		        eleIter1 != _electronsHandle->end();
		        eleIter1++) {

			if(! elePreselection(*eleIter1)) continue;
			if(!eleIter1->ecalDrivenSeed()) continue; // skip tracker-driven only electrons
			if(_eventType == WENU) {
				if(! (eleIter1->electronID(_eleID_tight)) ) continue;
				if( nTight != 1 || nLoose > 0 ) continue; //to be a Wenu event request only 1 ele WP70 in the event

				// MET/MT selection
				if(  met.et() < 25. ) continue;
				if( sqrt( 2.*eleIter1->et()*met.et() * (1 - cos(eleIter1->phi() - met.phi()))) < 50. ) continue;
				if( eleIter1->et() < 30) continue;

				doFill = true;
				if(_eventType == UNKNOWN) _eventType = WENU;
				TreeSetSingleElectronVar(*eleIter1, 0);  //fill first electron
				TreeSetSingleElectronVar(*eleIter1, -1); // fill fake second electron

				if(_doExtraCalibTree || _doExtraStudyTree) {
					TreeSetExtraCalibVar(*eleIter1, 0);
					TreeSetExtraCalibVar(*eleIter1, -1);
				}
				if(_doEleIDTree) {
					TreeSetEleIDVar(*eleIter1, 0);
					TreeSetEleIDVar(*eleIter1, -1);
				}
				// always set track vars, two of them are in the main _tree
				TreeSetTrackVar(*eleIter1, 0);
			} else { //ZEE or UNKNOWN
				// take only the fist di-electron pair (highest pt)
				for(pat::ElectronCollection::const_iterator eleIter2 = eleIter1 + 1;
				        eleIter2 != _electronsHandle->end() && doFill == false;
				        eleIter2++) {

					if(! elePreselection(*eleIter2)) continue;
					if(!eleIter2->ecalDrivenSeed()) continue; // skip tracker-driven only electrons
#ifdef DEBUG
					std::cout << "[DEBUG] Electron passing preselection" << std::endl;
#endif
					//	  float mass=(eleIter1->p4()+eleIter2->p4()).mass();

					//calculate the invariant mass
					double t1 = TMath::Exp(-eleIter1->eta());
					double t1q = t1 * t1;
					double t2 = TMath::Exp(-eleIter2->eta());
					double t2q = t2 * t2;
//					if(!eleIter2->parentSuperCluster().isNonnull()) continue;
					double angle = 1 -
					               ( (1 - t1q) * (1 - t2q) + 4 * t1 * t2 * cos(eleIter1->phi() - eleIter2->phi())) / (
					                   (1 + t1q) * (1 + t2q)
					               );
					double mass = sqrt(2 * eleIter1->energy() * eleIter2->energy() * angle); //use default electron energy, in order to have the same number of events between alcareco and alcarereco ntuples

					//	  std::cout<<" ele1 SC: "<<eleIter1->superCluster()->energy()<<" ele1 SC must: "<<eleIter1->parentSuperCluster()->energy()<<" eta1: "<<eleIter1->eta()<<" phi1: "<<eleIter1->phi()<<std::endl
					//		   <<" ele2 SC: "<<eleIter2->superCluster()->energy()<<" ele2 SC must: "<<eleIter2->parentSuperCluster()->energy()<<" eta2: "<<eleIter2->eta()<<" phi2: "<<eleIter2->phi()<<"mass: "<<mass<<std::endl;

					if(mass < 55 ) continue;
					doFill = true;

					if(_eventType == UNKNOWN) _eventType = ZEE;
					TreeSetDiElectronVar(*eleIter1, *eleIter2);
					if(_eleID[0] < 2 || (abs(_chargeEle[1]) == 1 && _eleID[1] < 2)) {
						// this event is not passing any _eleID, skip it
						doFill = false;
						continue;
					}
					if(_doExtraCalibTree || _doExtraStudyTree) {
						TreeSetExtraCalibVar(*eleIter1, *eleIter2);
					}
					if(_doEleIDTree) {
						TreeSetEleIDVar(*eleIter1, *eleIter2);
					}
					// always set track vars, two of them are in the main _tree
					TreeSetTrackVar(*eleIter1, *eleIter2);

					// if(_electronsHandle->size() < NELE &&  _eventType == SINGLEELE){

					// 	doFill=true;
					// 	TreeSetSingleElectronVar(*eleIter1, 0);  //fill first electron
					// 	TreeSetSingleElectronVar(*eleIter1, -1); // fill fake second electron

					// 	if(_doExtraCalibTree){
					// 		TreeSetExtraCalibVar(*eleIter1, 0);
					// 		TreeSetExtraCalibVar(*eleIter1, -1);
					// 	}
					// 	if(_doEleIDTree){
					// 		TreeSetEleIDVar(*eleIter1, 0);
					// 		TreeSetEleIDVar(*eleIter1, -1);
					// 	}
					// }
				}
			}
		}

	}	  else if (_eventType == ZMMG) {
		for( pat::MuonCollection::const_iterator muIter1 = _muonsHandle->begin();
		        muIter1 != _muonsHandle->end() && doFill == false;
		        muIter1++) {

			for(pat::MuonCollection::const_iterator muIter2 = muIter1 + 1;
			        muIter2 != _muonsHandle->end() && doFill == false;
			        muIter2++) {

				// should exit when muIter1 == end-1
				//if(! muIter2->electronID("loose") ) continue;
				for( pat::PhotonCollection::const_iterator phoIter1 = _photonsHandle->begin();
				        phoIter1 != _photonsHandle->end() && doFill == false;
				        phoIter1++) {

					float mass = (muIter1->p4() + muIter2->p4() + phoIter1->p4()).mass();

					if (phoIter1->pt() < 10) continue;
					if((mass < 55 || mass > 125)) continue;
					doFill = true;

					TreeSetMuMuGammaVar(*phoIter1, *muIter1, *muIter2);

					if(_doExtraCalibTree || _doExtraStudyTree) {
						TreeSetExtraCalibVar(*phoIter1, *muIter1, *muIter2);
					}
					if(_doEleIDTree) {
						TreeSetEleIDVar(*phoIter1, *muIter1, *muIter2);
					}
				}

			}
		}
	}

	if (_eventType == ZSC) { // || _eventType==UNKNOWN){ removed for miniAOD, to be put back in

		//leading pt electron in EB (patElectrons should be pt-ordered)
		// iterators storing pat Electons and HighEta SCs
		pat::ElectronCollection::const_iterator PatEle1 = _electronsHandle->begin();
		// iterators storing HighEta SCs
		// select the highest pt SC in the highEta region
		reco::SuperClusterCollection::const_iterator HighEtaSC1 = _EESuperClustersHandle->end();

		for( PatEle1 = _electronsHandle->begin();
		        //stop when HighEtaSC1 is a valid SC (this means that there is a pair matching the Z mass
		        PatEle1 != _electronsHandle->end();
		        PatEle1++) {
			if(!PatEle1->ecalDrivenSeed()) continue;

			// you have the first electrons candidate satifying the electrons criteria
			// now look for a SC matching the Z invariant mass. If not SC is found, let's look to another electrons candidate

			double HighEtaSCPt = 0;
			double t1 = TMath::Exp(-PatEle1->eta());
			double t1q = t1 * t1;
			for( reco::SuperClusterCollection::const_iterator iter = _EESuperClustersHandle->begin();
			        iter != _EESuperClustersHandle->end();
			        iter++) {
				// only SCs in the high eta region
				if (fabs(iter->eta()) < _doHighEta_LowerEtaCut) continue;

				//calculate the invariant mass
				double t2 = TMath::Exp(-iter->eta());
				double t2q = t2 * t2;
				double angle = 1 -
				               ( (1 - t1q) * (1 - t2q) + 4 * t1 * t2 * cos(PatEle1->phi() - iter->phi())) / (
				                   (1 + t1q) * (1 + t2q)
				               );
				float mass = sqrt(2 * PatEle1->energy() * iter->energy() * angle);
				if((mass < 55 || mass > 125)) continue;

				//take the highest pt SC matching the Z mass
				double pt = iter->energy() / cosh(iter->eta());
				if(HighEtaSCPt < pt) {
					HighEtaSCPt = pt;
					HighEtaSC1 = iter;
				}
			}

			if(HighEtaSC1 != _EESuperClustersHandle->end()) break;
		}

		// if you have found an ele-SC pair matching the high eta criteria,
		// save the event in the _tree
		if(HighEtaSC1 != _EESuperClustersHandle->end()) {
			doFill = true;
			TreeSetDiElectronVar(*PatEle1, *HighEtaSC1);
			if(_doExtraCalibTree || _doExtraStudyTree) {
				TreeSetExtraCalibVar(*PatEle1, *HighEtaSC1);
			}
			// always set track vars, two of them are in the main _tree
			TreeSetTrackVar(*PatEle1, 0);
		}
	}

	if(doFill) {
		_tree->Fill();
		if(_doExtraCalibTree) _extraCalibTree->Fill();
		if(_doEleIDTree)      _eleIDTree->Fill();
		if(_doExtraStudyTree) _extraStudyTree->Fill();
		if(_doTrackTree) _trackTree->Fill();
	}
	delete _clustertools;
	return;
}


// ------------ method called once each job just after ending the event loop  ------------
void ZNtupleDumper::endJob()
{
	if(_tree->GetEntries() > 0) {
		_tree->BuildIndex("runNumber", "eventNumber");
		if(_doEleIDTree)       _eleIDTree->BuildIndex("runNumber", "eventNumber");
		if(_doExtraCalibTree) _extraCalibTree->BuildIndex("runNumber", "eventNumber");
		if(_doExtraStudyTree) _extraStudyTree->BuildIndex("runNumber", "eventNumber");
		if(_doTrackTree) _trackTree->BuildIndex("runNumber", "eventNumber");
	}
	// save the _tree into the file
	_tree_file->cd();
	_tree->Write();
	_tree_file->Close();

	if(_doExtraCalibTree) {
		_extraCalibTreeFile->cd();
		_extraCalibTree->Write();
		_extraCalibTreeFile->Close();
	}
	if(_doExtraStudyTree) {
		_extraStudyTreeFile->cd();
		_extraStudyTree->Write();
		_extraStudyTreeFile->Close();
	}
	if(_doEleIDTree) {
		_eleIDTreeFile->cd();
		_eleIDTree->Write();
		_eleIDTreeFile->Close();
	}
	if(_doTrackTree) {
		_trackTreeFile->cd();
		_trackTree->Write();
		_trackTreeFile->Close();
	}
}


// ------------ method called when starting to processes a run  ------------
void ZNtupleDumper::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{
}


// ------------ method called when ending the processing of a run  ------------
void ZNtupleDumper::endRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when starting to processes a luminosity block  ------------
void ZNtupleDumper::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a luminosity block  ------------
void ZNtupleDumper::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ZNtupleDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


// output tree initialization
void ZNtupleDumper::InitNewTree()
{
	//  _tree = new TTree("selected",fChain->GetTitle());
	std::cout << "[STATUS] InitNewTree" << std::endl;
	if(_tree == NULL) return;
	_tree->Branch("runNumber",     &_runNumber,   "runNumber/i");
	_tree->Branch("lumiBlock",     &_lumiBlock,   "lumiBlock/s");
	_tree->Branch("eventNumber",   &_eventNumber, "eventNumber/l");
	_tree->Branch("eventTime",       &_eventTime,     "eventTime/i");
	_tree->Branch("nBX",           &_nBX,         "nBX/s");

	_tree->Branch("mcGenWeight",   &_mcGenWeight, "mcGenWeight/F");

	_tree->Branch("HLTfire", &_HLTfire, "HLTfire/B");
	//_tree->Branch("HLTNames",&(_HLTNames[0]));
	//_tree->Branch("HLTResults",&(_HLTResults[0]));

	// pu
	_tree->Branch("rho", &_rho, "rho/F");
	_tree->Branch("nPV", &_nPV, "nPV/b");
	_tree->Branch("nPU", &_nPU, "nPU/b");

	// vertex
	_tree->Branch("vtxX", &_vtxX, "vtxX/F");
	_tree->Branch("vtxY", &_vtxY, "vtxY/F");
	_tree->Branch("vtxZ", &_vtxZ, "vtxZ/F");

	// ele
	_tree->Branch("eleID", _eleID, "eleID[3]/i");
	_tree->Branch("chargeEle",   _chargeEle,    "chargeEle[3]/S");
	_tree->Branch("recoFlagsEle", _recoFlagsEle, "recoFlagsEle[3]/b");
	_tree->Branch("etaEle",      _etaEle,       "etaEle[3]/F");
	_tree->Branch("phiEle",      _phiEle,       "phiEle[3]/F");
	_tree->Branch("fbremEle",    _fbremEle,     "fbremEle[3]/F");
	_tree->Branch("R9Ele", _R9Ele, "R9Ele[3]/F");

	// ele track lenths
	_tree->Branch("gsfTrackLengthFromVtxP", _gsfTrackLengthFromVtxP, "gsfTrackLengthFromVtxP[3]/F");
	_tree->Branch("gsfTrackLengthFromTangents", _gsfTrackLengthFromTangents, "gsfTrackLengthFromTangents[3]/F");

	// SC
	_tree->Branch("etaSCEle",      _etaSCEle,       "etaSCEle[3]/F");
	_tree->Branch("phiSCEle",      _phiSCEle,       "phiSCEle[3]/F");
	_tree->Branch("nHitsSCEle", _nHitsSCEle, "nHitsSCEle[3]/I");
	_tree->Branch("avgLCSC",  _avgLCSC,       "avgLCSC[3]/F");
//	_tree->Branch("isMustacheSC", _isMustacheSC, "isMustacheSC[3]/b");
	_tree->Branch("rawEnergySCEle", _rawEnergySCEle, "rawEnergySCEle[3]/F");
	_tree->Branch("energy_ECAL_ele", _energy_ECAL_ele, "energy_ECAL_ele[3]/F"); ///< correctedEcalEnergy from MINIAOD or mustache regression if rereco
	_tree->Branch("energy_ECAL_pho", _energy_ECAL_pho, "energy_ECAL_pho[3]/F");
	_tree->Branch("energyUncertainty_ECAL_ele", _energyUncertainty_ECAL_ele, "energyUncertainty_ECAL_ele[3]/F"); ///< correctedEcalEnergy from MINIAOD or mustache regression if rereco
	_tree->Branch("energyUncertainty_ECAL_pho", _energyUncertainty_ECAL_pho, "energyUncertainty_ECAL_pho[3]/F");

	_tree->Branch("energy_5x5SC", _energy_5x5SC, "energy_5x5SC[3]/F");
	_tree->Branch("pModeGsfEle", _pModeGsfEle, "pModeGsfEle[3]/F");
	_tree->Branch("pAtVtxGsfEle", _pAtVtxGsfEle, "pAtVtxGsfEle[3]/F");
	_tree->Branch("pNormalizedChi2Ele", _pNormalizedChi2Ele, "pNormalizedChi2Ele[3]/F");
	_tree->Branch("trackMomentumErrorEle", _trackMomentumErrorEle, "trackMomentumErrorEle[3]/F");

	// seed recHit
	_tree->Branch("xSeedSC",            _xSeedSC,            "xSeedSC[3]/S");
	_tree->Branch("ySeedSC",            _ySeedSC,            "ySeedSC[3]/S");
	_tree->Branch("gainSeedSC", _gainSeedSC, "gainSeedSC[3]/b");
	_tree->Branch("energySeedSC",       _energySeedSC,       "energySeedSC[3]/F");
	_tree->Branch("energySecondToSeedSC",       _energySecondToSeedSC,       "energySecondToSeedSC[3]/F");
	_tree->Branch("amplitudeSeedSC",       _amplitudeSeedSC,       "amplitudeSeedSC[3]/F");
	_tree->Branch("amplitudeSecondToSeedSC",       _amplitudeSecondToSeedSC,       "amplitudeSecondToSeedSC[3]/F");
	_tree->Branch("timeSeedSC",       _timeSeedSC,       "timeSeedSC[3]/F");
	_tree->Branch("timeSecondToSeedSC",       _timeSecondToSeedSC,       "timeSecondToSeedSC[3]/F");
	_tree->Branch("icSeedSC",        _icSeedSC,        "icSeedSC[3]/F");
	_tree->Branch("laserSeedSC",     _laserSeedSC,     "laserSeedSC[3]/F");
	_tree->Branch("alphaSeedSC",     _alphaSeedSC,     "alphaSeedSC[3]/F");
	_tree->Branch("pedestalSeedSC",  _pedestalSeedSC,  "pedestalSeedSC[3]/F");
	_tree->Branch("noiseSeedSC",  _noiseSeedSC,  "noiseSeedSC[3]/F");

	// ES
	_tree->Branch("esEnergySCEle", _esEnergySCEle, "esEnergySCEle[3]/F");
	_tree->Branch("esEnergyPlane2SCEle", _esEnergyPlane2SCEle, "esEnergyPlane2SCEle[3]/F");
	_tree->Branch("esEnergyPlane1SCEle", _esEnergyPlane1SCEle, "esEnergyPlane1SCEle[3]/F");
	_tree->Branch("rawESEnergyPlane2SCEle", _rawESEnergyPlane2SCEle, "rawESEnergyPlane2SCEle[3]/F");
	_tree->Branch("rawESEnergyPlane1SCEle", _rawESEnergyPlane1SCEle, "rawESEnergyPlane1SCEle[3]/F");

	// MC truth
	_tree->Branch("etaMCEle",      _etaMCEle,       "etaMCEle[3]/F");	//[nEle]
	_tree->Branch("phiMCEle",      _phiMCEle,       "phiMCEle[3]/F");	//[nEle]
	_tree->Branch("invMass_MC", &_invMass_MC, "invMass_MC/F");

	_tree->Branch("ZEvent",      &_ZEvent,       "ZEvent/O");	//[nEle]

	// invariant mass
	_tree->Branch("invMass",    &_invMass,      "invMass/F");
	_tree->Branch("invMass_ECAL_ele", &_invMass_ECAL_ele, "invMass_ECAL_ele/F"); ///< using correctedEcalEnergy or using mustache SC dedicated regression
	_tree->Branch("invMass_ECAL_pho", &_invMass_ECAL_pho,   "invMass_ECAL_pho/F");
	//   _tree->Branch("invMass_e3x3",    &invMass_e3x3,      "invMass_e3x3/F");
	_tree->Branch("invMass_5x5SC",    &_invMass_5x5SC,      "invMass_5x5SC/F");
	//_tree->Branch("invMass_efull5x5",    &invMass_efull5x5,      "invMass_efull5x5/F");
	_tree->Branch("invMass_rawSC", &_invMass_rawSC,   "invMass_rawSC/F");
	_tree->Branch("invMass_rawSC_esSC", &_invMass_rawSC_esSC,   "invMass_rawSC_esSC/F");
	_tree->Branch("invMass_highEta", &_invMass_highEta,   "invMass_highEta/F");

	return;
}


void ZNtupleDumper::ResetMainTreeVar()
{
	// limit the reset to variables that may not always be filled
	// ele
	for (int i = 0; i < NELE; ++i) {
		_eleID[i] = initSingleInt;

		_chargeEle[i] = initSingleIntCharge;
		_recoFlagsEle[i] = 0;
		_etaEle[i] = initSingleFloat;
		_phiEle[i] = initSingleFloat;
		_R9Ele[i] = initSingleFloat;
		_fbremEle[i] = initSingleFloat;

		_isMustacheSC[i] = initSingleInt;
		_etaSCEle[i] = initSingleFloat;
		_phiSCEle[i] = initSingleFloat;

		_xSeedSC[i] = -999;
		_ySeedSC[i] = -999;
		_gainSeedSC[i] = 9;
		_energySeedSC[i] = initSingleFloat;
		_timeSeedSC[i] = initSingleFloat;
		_icSeedSC[i] = initSingleFloat;
		_laserSeedSC[i] = initSingleFloat;
		_avgLCSC[i] = initSingleFloat;
		_alphaSeedSC[i] = initSingleFloat;
		_pedestalSeedSC[i] = initSingleFloat;
		_noiseSeedSC[i] = initSingleFloat;
		_amplitudeSeedSC[i] = initSingleFloat;
		_energySecondToSeedSC[i] = initSingleFloat;
		_timeSecondToSeedSC[i] = initSingleFloat;
		_amplitudeSecondToSeedSC[i] = initSingleFloat;

		_energyEle[i] = initSingleFloat;
		_rawEnergySCEle[i] = initSingleFloat;
		_energy_ECAL_ele[i] = initSingleFloat;
		_energy_ECAL_pho[i] = initSingleFloat;
		_energyUncertainty_ECAL_ele[i] = initSingleFloat;
		_energyUncertainty_ECAL_pho[i] = initSingleFloat;

		_esEnergySCEle[i] = initSingleFloat;
		_esEnergyPlane1SCEle[i] = initSingleFloat;
		_esEnergyPlane2SCEle[i] = initSingleFloat;
		_rawESEnergyPlane1SCEle[i] = initSingleFloat;
		_rawESEnergyPlane2SCEle[i] = initSingleFloat;

		_energy_3x3SC[i] = initSingleFloat;
		_energy_5x5SC[i] = initSingleFloat;
		_eBCseedEle[i] = initSingleFloat;
		_pModeGsfEle[i] = initSingleFloat;
		_pAtVtxGsfEle[i] = initSingleFloat;
		_trackMomentumErrorEle[i] = initSingleFloat;
		_pNormalizedChi2Ele[i] = initSingleFloat;
		_gsfTrackLengthFromVtxP[i] = initSingleFloat;
		_gsfTrackLengthFromTangents[i] = initSingleFloat;
		_etaMCEle[i] = initSingleFloat;
		_phiMCEle[i] = initSingleFloat;
		_energyMCEle[i] = initSingleFloat;
	}

	_invMass = initSingleFloat;
	_invMass_rawSC = initSingleFloat;
	_invMass_rawSC_esSC = initSingleFloat;
	_invMass_ECAL_ele = initSingleFloat;
	_invMass_ECAL_pho = initSingleFloat;
	_invMass_5x5SC = initSingleFloat;
	_invMass_highEta = initSingleFloat;
	_invMass_mumu = initSingleFloat;
	_invMass_MC = initSingleFloat;
}


void ZNtupleDumper::TreeSetEventSummaryVar(const edm::Event& iEvent)
{
	_eventTimeStamp   =  iEvent.eventAuxiliary().time();
	_eventTime = (UInt_t) _eventTimeStamp.unixTime();
	_runNumber = (UInt_t) iEvent.run();
	_eventNumber = (Long64_t) iEvent.id().event();
//	_nBX = (UShort_t) (iEvent.isRealData()) ? iEvent.bunchCrossing() : 0;
	_nBX = (UShort_t)  iEvent.bunchCrossing();

	_lumiBlock = (UShort_t) iEvent.luminosityBlock();

	if(!hltPaths.empty()) {
		edm::TriggerNames hltNames = iEvent.triggerNames(*_triggerResultsHandle);
		int hltCount = _triggerResultsHandle->size();
		_HLTNames[0].clear();
		_HLTBits.clear();
		for (int i = 0 ; i < hltCount; ++i) {
			std::string hltName_str(hltNames.triggerName(i));
			(_HLTNames[0]).push_back(hltName_str);
			(_HLTResults[0]).push_back(_triggerResultsHandle->accept(i));
			_HLTBits.insert(std::pair<std::string, bool>( hltName_str, _triggerResultsHandle->accept(i)));
		} // for i
	}

	return;
}


void ZNtupleDumper::TreeSetPrimaryVertex(void)
{
	_vtxX = initSingleFloat;
	_vtxY = initSingleFloat;
	_vtxZ = initSingleFloat;
	for(reco::VertexCollection::const_iterator v = _primaryVertexHandle->begin(); v != _primaryVertexHandle->end(); ++v) {
		_vtxX = (*v).x();
		_vtxY = (*v).y();
		_vtxZ = (*v).z();
		break;
	}
}


void ZNtupleDumper::TreeSetPileupVar(void)
{
	_rho = *_rhoHandle;
	_nPV = 255;
	_nPU = 255;
	_mcGenWeight = -1;

	if(_primaryVertexHandle->size() > 0) {
		for(reco::VertexCollection::const_iterator v = _primaryVertexHandle->begin(); v != _primaryVertexHandle->end(); ++v) {
			//if((*v).tracksSize() > 0)
			_nPV++; // non mi ricordo perche' ho fatto cosi'....
		}
	}

	if(_isMC) {
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		for(PVI = _PupInfo->begin(); PVI != _PupInfo->end(); ++PVI) {
			int BX = PVI->getBunchCrossing();
			if(BX == 0) { // in-time pu
				_nPU = PVI->getTrueNumInteractions();
			}
		}

		if(!_GenEventInfoHandle->weights().empty()) {
			_mcGenWeight = (_GenEventInfoHandle->weights())[0];
		}
	}
	return;
}


// a negative index means that the corresponding electron does not exist, fill with 0
void ZNtupleDumper::TreeSetSingleElectronVar(const pat::Electron& electron, int index)
{
	if(index < 0) {
		_chargeEle[-index] = -100;
		_etaEle[-index]    = 0;
		_phiEle[-index]    = 0;
		_recoFlagsEle[-index] = -1;
		return;
	}
	assert(electron.ecalDrivenSeed());

	_energyEle[index] = electron.energy();
	_chargeEle[index] = (Char_t)electron.charge();
	_etaEle[index]    = electron.eta(); // degli elettroni
	_phiEle[index]    = electron.phi();
	_fbremEle[index]  = electron.fbrem(); // fatta con pIn-pOut

	if(electron.ecalDrivenSeed()) {
		if(electron.trackerDrivenSeed()) _recoFlagsEle[index] = 3;
		else _recoFlagsEle[index] = 2;
	} else _recoFlagsEle[index] = 1;

	_isMustacheSC[index] = !electron.superCluster().isNonnull();
	assert(_isMustacheSC);
	const reco::SuperClusterRef& sc = _isMustacheSC[index] ?  electron.parentSuperCluster() : electron.superCluster();
	assert(sc.isNonnull()); // at least one of the SCs has to be defined!

	TreeSetSingleSCVar(*sc, index);
	_energy_ECAL_ele[index]			  = (_isMINIAOD) ? electron.correctedEcalEnergy()     : electron.userFloat("energySCEleMust");
	//assert(!(_energy_ECAL_ele[index]<0 && _recoFlagsEle[index]>1));
	_energyUncertainty_ECAL_ele[index] = (_isMINIAOD) ? electron.correctedEcalEnergyError() : electron.userFloat("energySCEleMustVar");

	_energy_ECAL_pho[index]			  = electron.userFloat("energySCElePho");
	_energyUncertainty_ECAL_pho[index] = electron.userFloat("energySCElePhoVar");

	// change in an electron properties please, EleNewEnergyProducer
	_energy_3x3SC[index] = _clustertools->e3x3(*sc->seed());
	_energy_5x5SC[index] = electron.full5x5_e5x5();

	Float_t efull5x5SCEle = _clustertools->e5x5(*sc->seed()); // this is the official one, in rerecoes it is the one of the prompt reco, not updated by the rereco
	if(fabs(efull5x5SCEle - _energy_5x5SC[index]) > 1e-6) {
		std::cout << "[WARNING] " << efull5x5SCEle << "\t" << _energy_5x5SC[index] << "\t" << efull5x5SCEle - _energy_5x5SC[index] << "\t" << _clustertools->n5x5(*sc->seed()) << "\t" << _etaEle[index] << "\t" << _runNumber << "\t" << _eventNumber << std::endl;
	}
	//if(efull5x5SCEle>15) assert(fabs(efull5x5SCEle - _energy_5x5SC[index])< 1e-6);

	_pModeGsfEle[index] = electron.gsfTrack()->pMode();
	_trackMomentumErrorEle[index] = electron.trackMomentumError();
	_pNormalizedChi2Ele[index] = electron.gsfTrack()->normalizedChi2();
	_pAtVtxGsfEle[index] = electron.trackMomentumAtVtx().R();

	//_R9Ele[index] = e3x3SCEle[index] / sc->rawEnergy(); //already commented
	_R9Ele[index] = electron.full5x5_r9();
	if(fabs(_R9Ele[index] - _energy_3x3SC[index] / sc->rawEnergy()) > 1e-4) {
		std::cerr << "[WARNING] R9 different: " << _runNumber << "\t" << _eventNumber << "\t" << _etaEle[index] << "\t" << _R9Ele[index] << "\t" << _energy_3x3SC[index] / sc->rawEnergy() << std::endl;
	}

	eleIDMap eleID_map;

	_eleID[index] = 0;
	for (std::map<std::string, UInt_t>::iterator it = eleID_map.eleIDmap.begin(); it != eleID_map.eleIDmap.end(); ++it) {
#ifdef check_newID
		if((it->first).compare("loose25nsRun2V2016") == 0) {
			std::cout << "eleID is " << it->first << "; isAvailable is " << electron.isElectronIDAvailable(it->first) << std::endl;
		}
#endif
		if(electron.isElectronIDAvailable(it->first)) { //
			if ((bool) electron.electronID(it->first))  _eleID[index] |= it->second;//
		}//
	}
//	return;//To be fixed
	//  _sigmaIEtaIEtaSCEle[index] = sqrt(_clustertools->localCovariances(*(electron.superCluster()->seed()))[index]);
	//  _sigmaIEtaIEtaSCEle[1] = sqrt(_clustertools->localCovariances(*(electron2.superCluster()->seed()))[index]);
	//  sigmaIPhiIPhiBC = sqrt(_clustertools->localCovariances(*b)[3]);
	//  sigmaIEtaIPhiBC = _clustertools->localCovariances(*b)[1];

	if(electron.genLepton() != 0) {
		_energyMCEle[index]  = electron.genLepton()->energy();
		_etaMCEle[index]  = electron.genLepton()->eta();
		_phiMCEle[index]  = electron.genLepton()->phi();
	} else {
		_energyMCEle[index] = 0;
		_etaMCEle[index] = 0;
		_phiMCEle[index] = 0;
	}

	return;
}


void ZNtupleDumper::TreeSetTrackVar(const pat::Electron& electron1, const pat::Electron & electron2)
{
	TreeSetTrackVar(electron1, 0);
	TreeSetTrackVar(electron2, 1);
}


// can very likely find a better usage of GlobalPoint and GlobalVector
void ZNtupleDumper::TreeSetTrackVar(const pat::Electron& electron, int index)
{
	reco::GsfTrackRef track = electron.gsfTrack();
	edm::ESHandle<MagneticField> magFieldHandle;
	_pSetup->get<IdealMagneticFieldRecord>().get(magFieldHandle);
	auto magField = magFieldHandle.product();

	// momentum and position at vertex
	GlobalPoint gx_vtx(track->vx(), track->vy(), track->vz());
	GlobalVector gp_vtx(track->px(), track->py(), track->pz());
	GlobalVector gp_ele(electron.px(), electron.py(), electron.pz());

	// fill track momentum at vertex and vertex position
	_gsfTrackVtxPt[index]  = track->pt();
	_gsfTrackVtxEta[index] = track->eta();
	_gsfTrackVtxPhi[index] = track->phi();
	_gsfTrackVtxX[index]   = track->vx();
	_gsfTrackVtxY[index]   = track->vy();
	_gsfTrackVtxZ[index]   = track->vz();

	// position extrapolated at the impact point on the calorimeter
	auto tpc = electron.trackPositionAtCalo();
	GlobalPoint gx_calo(tpc.x(), tpc.y(), tpc.z());
	_gsfTrackCaloX[index] = tpc.x();
	_gsfTrackCaloY[index] = tpc.y();
	_gsfTrackCaloZ[index] = tpc.z();

	//fprintf(stderr, "--> pt: %f  eta: %f  phi: %f  %u  fbrem: %f\n", electron.pt(), electron.eta(), electron.phi(), track->gsfExtra()->tangentsSize(), electron.fbrem());

	SteppingHelixPropagator helix(magField);
	// path length in one go from momentum at vertex
	FreeTrajectoryState trajectory(gx_vtx, gp_vtx, track->charge(), magField);
	auto tpOne = helix.propagateWithPath(trajectory, gx_calo);
	//fprintf(stderr, "path length in one go: %f (in: %f %f %f, pt = %f)\n", tpOne.second, gp_vtx.x(), gp_vtx.y(), gp_vtx.z(), sqrt(gp_vtx.x()*gp_vtx.x() + gp_vtx.y()*gp_vtx.y()));
	_gsfTrackLengthFromVtxP[index] = tpOne.second;

	// path length in one go from electron momentum
	FreeTrajectoryState trajectoryEle(gx_vtx, gp_ele, track->charge(), magField);
	auto tpOneEle = helix.propagateWithPath(trajectoryEle, gx_calo);
	_gsfTrackLengthFromEleP[index] = tpOneEle.second;

	// path length in one go from outer momentum
	GlobalVector pOut(electron.trackMomentumOut().x(), electron.trackMomentumOut().y(), electron.trackMomentumOut().z());

	// fill track momentum at outermost layer
	_gsfTrackOuterPt[index]  = electron.trackMomentumOut().rho();
	_gsfTrackOuterEta[index] = electron.trackMomentumOut().eta();
	_gsfTrackOuterPhi[index] = electron.trackMomentumOut().phi();

	// propagate from vtx with outer momentum
	FreeTrajectoryState trajectoryOut(gx_vtx, pOut, track->charge(), magField);
	auto tpOut = helix.propagateWithPath(trajectoryOut, gx_calo);
	//fprintf(stderr, "path length in one go: %f (out: %f %f %f, pt = %f)\n", tpOut.second, pOut.x(), pOut.y(), pOut.z(), sqrt(pOut.x()*pOut.x() + pOut.y()*pOut.y()));
	_gsfTrackLengthFromOuterP[index] = tpOut.second;

	if(!_isMINIAOD) {
		// path length in steps
		double pathLength = 0.;
		auto & tangents = track->gsfExtra()->tangents();
		if (!tangents.size()) {
			return;
		}
		reco::GsfTrackExtra::Point start_x, stop_x;
		reco::GsfTrackExtra::Vector start_p;
		auto & tg = tangents[0];

		// from vtx to first tangent
		stop_x = tg.position();
		GlobalPoint gx_stop(stop_x.x(), stop_x.y(), stop_x.z());
		FreeTrajectoryState traj_beg(gx_vtx, gp_vtx, track->charge(), magField);
		auto tp = helix.propagateWithPath(traj_beg, gx_stop);
		double length = tp.second;
		pathLength += length;

		// for next step (the first inside the loop)
		start_x = stop_x;
		start_p = tg.momentum();

		// fill first tangent
		_gsfTrackTangentsPt[index].push_back(start_p.rho());
		_gsfTrackTangentsEta[index].push_back(start_p.eta());
		_gsfTrackTangentsPhi[index].push_back(start_p.phi());
		_gsfTrackTangentsDeltaP[index].push_back(tg.deltaP().value());
		_gsfTrackTangentsDeltaPErr[index].push_back(tg.deltaP().error());
		_gsfTrackTangentsX[index].push_back(start_x.x());
		_gsfTrackTangentsY[index].push_back(start_x.y());
		_gsfTrackTangentsZ[index].push_back(start_x.z());

		for (size_t it = 1; it < tangents.size(); ++it) {
			auto tg = tangents[it];
			stop_x = tg.position();
			GlobalPoint gx_start(start_x.x(), start_x.y(), start_x.z());
			GlobalPoint gx_stop(stop_x.x(), stop_x.y(), stop_x.z());
			GlobalVector gp_start(start_p.x(), start_p.y(), start_p.z());
			FreeTrajectoryState traj(gx_start, gp_start, track->charge(), magField);
			tp = helix.propagateWithPath(traj, gx_stop);
			length = tp.second;
			pathLength += length;
			//fprintf(stderr, "segment length: %f (deltaP = %f +/- %f)  total length: %f\n", length, tg.deltaP().value(), tg.deltaP().error(), pathLength);
			start_x = stop_x;
			start_p = tg.momentum();
			// fill tangents
			_gsfTrackTangentsPt[index].push_back(start_p.rho());
			_gsfTrackTangentsEta[index].push_back(start_p.eta());
			_gsfTrackTangentsPhi[index].push_back(start_p.phi());
			_gsfTrackTangentsDeltaP[index].push_back(tg.deltaP().value());
			_gsfTrackTangentsDeltaPErr[index].push_back(tg.deltaP().error());
			_gsfTrackTangentsX[index].push_back(start_x.x());
			_gsfTrackTangentsY[index].push_back(start_x.y());
			_gsfTrackTangentsZ[index].push_back(start_x.z());
		}
		GlobalPoint gx_start(start_x.x(), start_x.y(), start_x.z());
		GlobalVector gp_start(start_p.x(), start_p.y(), start_p.z());
		FreeTrajectoryState traj_end(gx_start, gp_start, track->charge(), magField);
		tp = helix.propagateWithPath(traj_end, gx_calo);
		length = tp.second;
		pathLength += length;

		// fill track length from steps
		_gsfTrackLengthFromTangents[index] = pathLength;
		//fprintf(stderr, "segment length: %f  total length: %f\n", length, pathLength);
	} // else it is initSingleFloat
	return;
}


void ZNtupleDumper::TreeSetSingleSCVar(const reco::SuperCluster& sc, int index)
{
	_etaSCEle[index] = sc.eta();
	_phiSCEle[index] = sc.phi();

	_nHitsSCEle[index] = sc.size();

	_rawEnergySCEle[index]	   = sc.rawEnergy();
	_eBCseedEle[index] = sc.seed()->energy();

	_esEnergySCEle[index]	   = sc.preshowerEnergy();
	_esEnergyPlane1SCEle[index] = sc. preshowerEnergyPlane1();
	_esEnergyPlane2SCEle[index] = sc. preshowerEnergyPlane2();
	_rawESEnergyPlane1SCEle[index] = GetESPlaneRawEnergy(sc, 1);
	_rawESEnergyPlane2SCEle[index] = GetESPlaneRawEnergy(sc, 2);

	// change in an electron properties please, EleNewEnergyProducer
	_energy_3x3SC[index] = _clustertools->e3x3(*sc.seed());
	_energy_5x5SC[index] = _clustertools->e5x5(*sc.seed());

	_R9Ele[index] = _energy_3x3SC[index] / _rawEnergySCEle[index]; //already commented

	DetId seedDetId = sc.seed()->seed();
	if(seedDetId.subdetId() == EcalBarrel) {
		EBDetId seedDetIdEcal(seedDetId);
		_xSeedSC[index] = seedDetIdEcal.ieta();
		_ySeedSC[index] = seedDetIdEcal.iphi();
	} else {
		EEDetId seedDetIdEcal(seedDetId);
		_xSeedSC[index] = seedDetIdEcal.ix();
		_ySeedSC[index] = seedDetIdEcal.iy();
	}

	bool isEB = (seedDetId.subdetId() == EcalBarrel);

	const EcalUncalibratedRecHitCollection * uncHits = isEB ? _pEBUncRecHits.product() : _pEEUncRecHits.product();
	_amplitudeSeedSC[index] = initSingleFloat;
	if (uncHits) {
		auto uHitSeed = uncHits->find(seedDetId) ;
		if(uHitSeed == uncHits->end()) {
			edm::LogError("ZNtupleDumper") << "No uncalibrated recHit found for xtal "  << seedDetId.rawId()
			                               << " in subdetector " << seedDetId.subdetId() << "; continuing...";
		} else {
			_amplitudeSeedSC[index] = uHitSeed->amplitude();
		}
	} else {
		fprintf(stderr, "NO UNCHITS - seed\n");
	}

	const EcalRecHitCollection *recHits = (seedDetId.subdetId() == EcalBarrel) ?  _clustertools->getEcalEBRecHitCollection() : _clustertools->getEcalEERecHitCollection();
	EcalRecHitCollection::const_iterator seedRecHit = recHits->find(seedDetId) ;
	assert(seedRecHit != recHits->end());
	_energySeedSC[index] = seedRecHit->energy();
	_timeSeedSC[index]   = seedRecHit->time();

	// second highest energy crystal
	auto snd_seed = findEnergySortedHit(sc, recHits, 1);
	_energySecondToSeedSC[index] = snd_seed.second;
	auto snd_rh = recHits->find(snd_seed.first);
	_timeSecondToSeedSC[index]   = snd_rh->time();
	_amplitudeSecondToSeedSC[index] = initSingleFloat;
	if (uncHits) {
		auto uHitSeed = uncHits->find(snd_seed.first) ;
		if(uHitSeed == uncHits->end()) {
			edm::LogError("ZNtupleDumper") << "No uncalibrated recHit found for xtal "  << seedDetId.rawId()
			                               << " in subdetector " << seedDetId.subdetId() << "; continuing...";
		} else {
			_amplitudeSecondToSeedSC[index] = uHitSeed->amplitude();
		}
	} else {
		fprintf(stderr, "NO UNCHITS - second to seed\n");
	}


	const EcalIntercalibConstantMap & icalMap = _clustertools->getEcalIntercalibConstants();
	_icSeedSC[index] = *(icalMap.find(seedDetId)); // unsafe, but never found a missing IC in years of running ;)

	const edm::ESHandle<EcalLaserDbService>& laserHandle = _clustertools->getLaserHandle();
	_laserSeedSC[index] = laserHandle->getLaserCorrection(seedDetId, _eventTimeStamp);

	edm::ESHandle<EcalLaserDbService> aSetup; //ALPHA PART
	_pSetup->get<EcalLaserDbRecord>().get( aSetup );
	const EcalLaserAlphas* myalpha =  aSetup->getAlphas();
	const EcalLaserAlphaMap& laserAlphaMap =  myalpha->getMap();
	_alphaSeedSC[index] = *(laserAlphaMap.find( seedDetId ));

	_gainSeedSC[index] = 0;
	if(seedRecHit->checkFlag(EcalRecHit::kHasSwitchToGain6)) _gainSeedSC[index] |= 0x01;
	if(seedRecHit->checkFlag(EcalRecHit::kHasSwitchToGain1)) _gainSeedSC[index] |= 0x02;


	_pedestalSeedSC[index] = _ped->find(seedDetId)->mean(1);
	_noiseSeedSC[index]    = _ped->find(seedDetId)->rms(1);

	float sumLC_E = 0.;
	float sumE = 0.;

	std::vector< std::pair<DetId, float> > hitsAndFractions = sc.hitsAndFractions();
	for (const auto&  detitr : hitsAndFractions) {
		EcalRecHitCollection::const_iterator oneHit = recHits->find( (detitr.first) ) ;
#ifdef DEBUG
		assert(oneHit != recHits->end()); // DEBUG
#endif
		sumLC_E += laserHandle->getLaserCorrection(detitr.first, _eventTimeStamp) * oneHit->energy();
		sumE    += oneHit->energy();
	}
	_avgLCSC[index] = sumLC_E / sumE;


	return;
}


// a negative index means that the corresponding electron does not exist, fill with 0
void ZNtupleDumper::TreeSetSingleElectronVar(const reco::SuperCluster& electron1, int index)
{
	if(index < 0) {
		_chargeEle[-index] = -100;
		_etaEle[-index]    = 0;
		_phiEle[-index]    = 0;
		return;
	}

	//checks

	_chargeEle[index] = 0; // dont know the charge for SC
	_etaEle[index]    = electron1.eta(); // eta = etaSC
	_phiEle[index]    = electron1.phi();

	_recoFlagsEle[index] = -1; // define -1 as a SC

	_fbremEle[index] = -1; // no bremstrahlung for SC

	_etaSCEle[index] = electron1.eta(); // itself is a SC
	_phiSCEle[index] = electron1.phi();

	// no MC matching has been considered yet for SCV
	_energyMCEle[index] = -100;
	_etaMCEle[index] = -100;
	_phiMCEle[index] = -100;


	_pModeGsfEle[index] = -1; // no track, though ..
	_trackMomentumErrorEle[index] = -1;
	_pNormalizedChi2Ele[index] = -1;
	_pAtVtxGsfEle[index] = -1;


	// temporary ignore the id and classification
	_eleID[index] = -100;

	TreeSetSingleSCVar(electron1, index);
	return;
}


// a negative index means that the corresponding muon does not exist, fill with 0
void ZNtupleDumper::TreeSetSingleMuonVar(const pat::Muon& muon1, int index)
{
	if(index < 0) {
		_chargeEle[-index] = -100;
		_etaEle[-index]    = 0;
		_phiEle[-index]    = 0;
		return;
	}

	_chargeEle[index] = muon1.charge();
	_etaEle[index]    = muon1.eta(); // degli elettroni
	_phiEle[index]    = muon1.phi();


	if(muon1.genLepton() != 0) {
		_energyMCEle[index]  = muon1.genLepton()->energy();
		_etaMCEle[index]  = muon1.genLepton()->eta();
		_phiMCEle[index]  = muon1.genLepton()->phi();
	} else {
		_energyMCEle[index] = 0;
		_etaMCEle[index] = 0;
		_phiMCEle[index] = 0;
	}


	// why the hell this does not work????
	//  _eleID[index] = ((bool) muon1.muonID("soft")) << 0;
	//  _eleID[index] += ((bool) muon1.muonID("loose")) << 1;
	//  _eleID[index] += ((bool) muon1.muonID("highPT")) << 2;
	//  _eleID[index] += ((bool) muon1.muonID("tight")) << 3;


	return;
}


// a negative index means that the corresponding electron does not exist, fill with 0
void ZNtupleDumper::TreeSetSinglePhotonVar(const pat::Photon& photon, int index)
{
	if(index < 0) {
		_chargeEle[-index] = -100;
		_etaEle[-index]    = 0;
		_phiEle[-index]    = 0;
		return;
	}

	_chargeEle[index] = (Char_t) photon.charge();
	_etaEle[index]    = photon.eta();
	_phiEle[index]    = photon.phi();

	_etaSCEle[index] = photon.superCluster()->eta();
	_phiSCEle[index] = photon.superCluster()->phi();

	TreeSetSingleSCVar(*(photon.superCluster()), index);

	if(photon.genParticle() != 0) {
		_energyMCEle[index]  = photon.genParticle()->energy();
		_etaMCEle[index]  = photon.genParticle()->eta();
		_phiMCEle[index]  = photon.genParticle()->phi();
	} else {
		_energyMCEle[index] = 0;
		_etaMCEle[index] = 0;
		_phiMCEle[index] = 0;
	}


	_rawEnergySCEle[index]  = photon.superCluster()->rawEnergy();
	_esEnergySCEle[index] = photon.superCluster()->preshowerEnergy();
	_esEnergyPlane1SCEle[index] = photon.superCluster()-> preshowerEnergyPlane1();
	_esEnergyPlane2SCEle[index] = photon.superCluster()-> preshowerEnergyPlane2();

	//  energySCEle_corr[index] = photon.scEcalEnergy(); //but, I don't think this is the correct energy..

	// change in an electron properties please, EleNewEnergyProducer
	_energy_3x3SC[index] = _clustertools->e3x3(*photon.superCluster()->seed());
	_energy_5x5SC[index] = _clustertools->e5x5(*photon.superCluster()->seed());
	_eBCseedEle[index] = photon.superCluster()->seed()->energy();
	_R9Ele[index] = _energy_3x3SC[index] / photon.superCluster()->rawEnergy(); //original

	//   if(_isMC){
	//     if(photon.isEB())
	//       _R9Ele[index] = _R9Ele[index]*1.0045+0.0010;
	//     else
	//       _R9Ele[index] = _R9Ele[index]*1.0086-0.0007;
	//   }

	// _eleID[index] = ((bool) photon.photonID("fiducial")) << 0;
	// _eleID[index] += ((bool) photon.photonID("loose")) << 1;
	// _eleID[index] += ((bool) photon.photonID("medium")) << 2;
	// _eleID[index] += ((bool) photon.photonID("tight")) << 3;
	// _eleID[index] += ((bool) photon.photonID("loose25nsRun2")) << 4;
	// _eleID[index] += ((bool) photon.photonID("medium25nsRun2")) << 5;
	// _eleID[index] += ((bool) photon.photonID("tight25nsRun2")) << 6;

	return;
}


void ZNtupleDumper:: TreeSetDiElectronVar(const pat::Electron& electron1, const pat::Electron& electron2)
{
	TreeSetSingleElectronVar(electron1, 0);
	TreeSetSingleElectronVar(electron2, 1);

	if (_doTrackTree) {
		TreeSetTrackVar(electron1, 0);
		TreeSetTrackVar(electron2, 1);
	}

	double t1 = TMath::Exp(-_etaEle[0]);
	double t2 = TMath::Exp(-_etaEle[1]);
	double t1q = t1 * t1;
	double t2q = t2 * t2;

	double angle = 1 -
	               ( (1 - t1q) * (1 - t2q) + 4 * t1 * t2 * cos(_phiEle[0] - _phiEle[1])) / (
	                   (1 + t1q) * (1 + t2q)
	               );


	_invMass		= sqrt(2 * _energyEle[0] * _energyEle[1] * angle);
	_invMass_5x5SC = sqrt(2 * _energy_5x5SC[0] * _energy_5x5SC[1] * angle);
	_invMass_highEta = sqrt(2 * _energy_ECAL_ele[0] * _energy_5x5SC[1] * angle);

	_invMass_ECAL_ele = sqrt(2 * _energy_ECAL_ele[0] * _energy_ECAL_ele[1] * angle);
	_invMass_ECAL_pho = sqrt(2 * _energy_ECAL_pho[0] * _energy_ECAL_pho[1] * angle);

	_invMass_rawSC = sqrt(2 * _rawEnergySCEle[0] * _rawEnergySCEle[1] * angle);

	_invMass_rawSC_esSC = sqrt(2 * (_rawEnergySCEle[0] + _esEnergySCEle[0]) *
	                           (_rawEnergySCEle[1] + _esEnergySCEle[1]) *
	                           angle);

	if(electron1.genLepton() != 0 && electron2.genLepton() != 0) {
		_invMass_MC     = sqrt(2 * electron1.genLepton()->energy() * electron2.genLepton()->energy() *
		                       angle);
	} else _invMass_MC = 0;
	//  invMass_genMC     = (electron1.genLepton()->p4 + electron2.genLepton()->p4()).M();


	// se non hanno fatto match con il MC?
	// qual e' la frazione di Z selezionate che non matchano il MC?

	//#ifdef shervin
	//  r9weight[0]=r9Weight(_etaEle[0], _R9Ele[0]);
	//  r9weight[1]=r9Weight(_etaEle[1], _R9Ele[1]);
	//#endif

	return;
}


void ZNtupleDumper::TreeSetDiElectronVar(const pat::Electron& electron1, const reco::SuperCluster& electron2)
{
	TreeSetSingleElectronVar(electron1, 0);
	TreeSetSingleElectronVar(electron2, 1);

	double t1 = TMath::Exp(-_etaEle[0]);
	double t2 = TMath::Exp(-_etaEle[1]);
	double t1q = t1 * t1;
	double t2q = t2 * t2;

	double angle = 1 - ( (1 - t1q) * (1 - t2q) + 4 * t1 * t2 * cos(_phiEle[0] - _phiEle[1])) / ( (1 + t1q) * (1 + t2q) );

	_invMass = sqrt(2 * electron1.energy() * electron2.energy() * angle);
	_invMass_5x5SC   = sqrt(2 * _energy_5x5SC[0] * _energy_5x5SC[1] * /// full 5x5
	                        angle);
	_invMass_highEta   = sqrt(2 * _energy_ECAL_ele[0] * _energy_5x5SC[1] * /// full 5x5
							  angle);

	_invMass_ECAL_ele = sqrt(2 * _energy_ECAL_ele[0] * _energy_ECAL_ele[1] * angle);
	_invMass_ECAL_pho = sqrt(2 * _energy_ECAL_pho[0] * _energy_ECAL_pho[1] * angle);

	_invMass_rawSC = sqrt(2 * _rawEnergySCEle[0] * _rawEnergySCEle[1] *
	                      angle);

	_invMass_rawSC_esSC = sqrt(2 * (_rawEnergySCEle[0] + _esEnergySCEle[0]) *
	                           (_rawEnergySCEle[1] + _esEnergySCEle[1]) *
	                           angle);


	_invMass_MC = -100; // temporary set it to be -100 for SC

	return;
}


void ZNtupleDumper:: TreeSetMuMuGammaVar(const pat::Photon& photon, const pat::Muon& muon1, const pat::Muon& muon2)
{
	TreeSetSinglePhotonVar(photon, 0);
	TreeSetSingleMuonVar(muon1, 1);
	TreeSetSingleMuonVar(muon2, 2);

	TLorentzVector *Z = new TLorentzVector();
	TLorentzVector *ph = new TLorentzVector();
	TLorentzVector *m1 = new TLorentzVector();
	TLorentzVector *m2 = new TLorentzVector();
//	ph->SetPtEtaPhiE(PtEle[0], _etaEle[0], _phiEle[0], photon.energy());
//	m1->SetPtEtaPhiE(PtEle[1], _etaEle[1], _phiEle[1], muon1.energy());
//	m2->SetPtEtaPhiE(PtEle[2], _etaEle[2], _phiEle[2], muon2.energy());
	*Z = *ph + *m1 + *m2;
	_invMass = Z->M();

	Z->SetE(photon.e5x5() + muon1.energy() + muon2.energy());
	_invMass_5x5SC = Z->M();

	Z->SetE(_rawEnergySCEle[0] + muon1.energy() + muon2.energy());
	_invMass_rawSC = Z->M();

	Z->SetE(_rawEnergySCEle[0] + _esEnergySCEle[0] + muon1.energy() + muon2.energy());
	_invMass_rawSC_esSC = Z->M();

	if(photon.genPhoton() != 0 && muon1.genLepton() != 0 && muon2.genLepton() != 0) {
		Z->SetE(photon.genPhoton()->energy() + muon1.genLepton()->energy() + muon2.genLepton()->energy());
		_invMass_MC = Z->M();
	} else _invMass_MC = 0;

	double t1 = TMath::Exp(-_etaEle[1]);
	double t2 = TMath::Exp(-_etaEle[2]);
	double t1q = t1 * t1;
	double t2q = t2 * t2;

	double angle = 1 -
	               ( (1 - t1q) * (1 - t2q) + 4 * t1 * t2 * cos(_phiEle[1] - _phiEle[2])) / (
	                   (1 + t1q) * (1 + t2q)
	               );

	_invMass     = sqrt(2 * muon1.energy() * muon2.energy() *
	                    angle);

	delete Z;
	delete ph;
	delete m1;
	delete m2;

	return;
}


//#============================== extra calib tree
void ZNtupleDumper::InitExtraCalibTree()
{
	//  _tree = new TTree("selected",fChain->GetTitle());
	std::cout << "[STATUS] InitExtraCalibTree" << std::endl;
	if(_extraCalibTree == NULL) return;

	_extraCalibTree->Branch("runNumber",     &_runNumber,   "runNumber/i");
	_extraCalibTree->Branch("lumiBlock",     &_lumiBlock,   "lumiBlock/s");
	_extraCalibTree->Branch("eventNumber",   &_eventNumber, "eventNumber/l");
	_extraCalibTree->Branch("eventTime",       &_eventTime,     "eventTime/i");

	_extraCalibTree->Branch("nHitsSCEle", _nHitsSCEle, "nHitsSCEle[3]/I");

	_extraCalibTree->Branch("recoFlagRecHitSCEle1", &(_recoFlagRecHitSCEle[0]));
	_extraCalibTree->Branch("recoFlagRecHitSCEle2", &(_recoFlagRecHitSCEle[1]));
	_extraCalibTree->Branch("rawIdRecHitSCEle1", &(_rawIdRecHitSCEle[0]));
	_extraCalibTree->Branch("rawIdRecHitSCEle2", &(_rawIdRecHitSCEle[1]));
	_extraCalibTree->Branch("XRecHitSCEle1", &(_XRecHitSCEle[0]));
	_extraCalibTree->Branch("XRecHitSCEle2", &(_XRecHitSCEle[1]));
	_extraCalibTree->Branch("YRecHitSCEle1", &(_YRecHitSCEle[0]));
	_extraCalibTree->Branch("YRecHitSCEle2", &(_YRecHitSCEle[1]));
	_extraCalibTree->Branch("ZRecHitSCEle1", &(_ZRecHitSCEle[0]));
	_extraCalibTree->Branch("ZRecHitSCEle2", &(_ZRecHitSCEle[1]));
	_extraCalibTree->Branch("energyRecHitSCEle1", &(_energyRecHitSCEle[0]));
	_extraCalibTree->Branch("energyRecHitSCEle2", &(_energyRecHitSCEle[1]));

	_extraCalibTree->Branch("fracRecHitSCEle1", &(_fracRecHitSCEle[0]));
	_extraCalibTree->Branch("fracRecHitSCEle2", &(_fracRecHitSCEle[1]));

	return;
}


void ZNtupleDumper::ResetExtraCalibVar()
{
	for (int i = 0; i < NELE; ++i) {
		_nRecHitsEle[i] = initSingleFloat;
		_nHitsSCEle[i] = initSingleFloat;
		_rawIdRecHitSCEle[i].clear();
		_XRecHitSCEle[i].clear();
		_YRecHitSCEle[i].clear();
		_ZRecHitSCEle[i].clear();
		_energyRecHitSCEle[i].clear();
		_fracRecHitSCEle[i].clear();
		_recoFlagRecHitSCEle[i].clear();
	}
}


void ZNtupleDumper::TreeSetExtraCalibVar(const pat::Electron& electron1, const pat::Electron& electron2)
{
	TreeSetExtraCalibVar(electron1, 0);
	TreeSetExtraCalibVar(electron2, 1);
	return;
}


void ZNtupleDumper::TreeSetExtraCalibVar(const pat::Electron& electron1, const reco::SuperCluster& electron2)
{
	TreeSetExtraCalibVar(electron1, 0);
	TreeSetExtraCalibVar(electron2, 1);
	return;
}


void ZNtupleDumper::TreeSetExtraCalibVar(const std::vector<std::pair<DetId, float> > & hitsFracs, int index, bool isEB)
{
	const EcalRecHitCollection *recHits = isEB ? _clustertools->getEcalEBRecHitCollection() : _clustertools->getEcalEERecHitCollection();
	const EcalUncalibratedRecHitCollection * uncHits = isEB ? _pEBUncRecHits.product() : _pEEUncRecHits.product();
	const EcalIntercalibConstantMap& icalMap = _clustertools->getEcalIntercalibConstants();
	const edm::ESHandle<EcalLaserDbService>& laserHandle = _clustertools->getLaserHandle();
	for (std::vector<std::pair<DetId, float> >::const_iterator detitr = hitsFracs.begin();
	        detitr != hitsFracs.end(); detitr++ ) {
		//      EcalRecHitCollection::const_iterator theSeedHit = recHits->find (id); // trash this
		EcalRecHitCollection::const_iterator oneHit = recHits->find( (detitr -> first) ) ;
		if(oneHit == recHits->end()) {
			edm::LogError("ZNtupleDumper") << "No intercalib const found for xtal "  << (detitr->first).rawId()
			                               << " in subdetector " << (detitr->first).subdetId() << " bailing out";
			//assert(0);
			continue;
		}
		_recoFlagRecHitSCEle[index].push_back(oneHit->recoFlag());
		_rawIdRecHitSCEle[index].push_back(detitr->first.rawId());
		if(isEB) {
			EBDetId recHitId(detitr->first);
			_XRecHitSCEle[index].push_back(recHitId.ieta());
			_YRecHitSCEle[index].push_back(recHitId.iphi());
			_ZRecHitSCEle[index].push_back(recHitId.zside());
		} else {
			EEDetId recHitId(detitr->first);
			_XRecHitSCEle[index].push_back(recHitId.ix());
			_YRecHitSCEle[index].push_back(recHitId.iy());
			_ZRecHitSCEle[index].push_back(recHitId.zside());
		}
		_energyRecHitSCEle[index].push_back(oneHit->energy());
		_fracRecHitSCEle[index].push_back(detitr->second);
		EcalUncalibratedRecHitCollection::const_iterator oneUHit = uncHits->find( (detitr -> first) ) ;
		if(oneUHit == uncHits->end()) {
			edm::LogError("ZNtupleDumper") << "No uncalibRecHit found for xtal "  << (detitr->first).rawId()
			                               << " in subdetector " << (detitr->first).subdetId() << " bailing out";
			//assert(0);
			continue;
		}
		// UncalibRecHit's information on OOT amplitudes
		float amplis[EcalDataFrame::MAXSAMPLES];
		for (int i = 0; i < EcalDataFrame::MAXSAMPLES; ++i) amplis[i] = oneUHit->outOfTimeAmplitude(i);
		_ootAmplisUncalibRecHitSCEle[index].push_back(std::vector<float>(amplis, amplis + EcalDataFrame::MAXSAMPLES));
		_ampliUncalibRecHitSCEle[index].push_back(oneUHit->amplitude());
		_ampliErrUncalibRecHitSCEle[index].push_back(oneUHit->amplitudeError());
		_pedEUncalibRecHitSCEle[index].push_back(oneUHit->pedestal());
		_jitterUncalibRecHitSCEle[index].push_back(oneUHit->jitter());
		_jitterErrUncalibRecHitSCEle[index].push_back(oneUHit->jitterError());
		_chi2UncalibRecHitSCEle[index].push_back(oneUHit->chi2());
		_flagsUncalibRecHitSCEle[index].push_back(oneUHit->flags());
		// in order to get back the ADC counts from the recHit energy, three ingredients are necessary:
		// 1) get laser correction coefficient
		_LCRecHitSCEle[index].push_back(laserHandle->getLaserCorrection(detitr->first, _eventTimeStamp));
		//laserHandle->
		// 2) get intercalibration
		EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detitr->first);
		EcalIntercalibConstant icalconst = 1.;
		if( icalit != icalMap.end() ) {
			icalconst = (*icalit);
			// std::cout << "icalconst set to: " << icalconst << std::endl;
		} else {
			edm::LogError("ZNtupleDumper") << "No intercalib const found for xtal "  << (detitr->first).rawId() << "bailing out";
			//assert(0);
			continue;
		}
		// 3) get adc2GeV
		//float adcToGeV = ( (detitr -> first).subdetId() == EcalBarrel ) ?
		// float(adcToGeVHandle->getEBValue()) : float(adcToGeVHandle->getEEValue());
		_ICRecHitSCEle[index].push_back(icalconst);
	}
}


void ZNtupleDumper::TreeSetExtraCalibVar(const pat::Electron& electron1, int index)
{
	if (index < 0) return;

	//  EcalIntercalibConstantMap icMap = icHandle->get()
	std::vector< std::pair<DetId, float> > hitsAndFractions_ele1 = electron1.superCluster()->hitsAndFractions();
	_nHitsSCEle[index] = hitsAndFractions_ele1.size();
	TreeSetExtraCalibVar(hitsAndFractions_ele1, index, electron1.isEB());
	return;
}


void ZNtupleDumper::TreeSetExtraCalibVar(const reco::SuperCluster& electron1, int index)
{
	if (index < 0) return;

	std::vector< std::pair<DetId, float> > hitsAndFractions_ele1 = electron1.hitsAndFractions();
	_nHitsSCEle[index] = hitsAndFractions_ele1.size();
	TreeSetExtraCalibVar(hitsAndFractions_ele1, index, electron1.seed()->seed().subdetId() == EcalBarrel);
	return;
}


void ZNtupleDumper::TreeSetExtraCalibVar(const pat::Photon& photon, const pat::Muon& muon1, const pat::Muon& muon2)
{
	TreeSetExtraCalibVar(photon, 0);
	TreeSetExtraCalibVar(muon1, -1);
	TreeSetExtraCalibVar(muon2, -2);
	return;
}


void ZNtupleDumper::TreeSetExtraCalibVar(const pat::Photon& photon, int index)
{
	if (index < 0) return;

	//  EcalIntercalibConstantMap icMap = icHandle->get()
	std::vector< std::pair<DetId, float> > hitsAndFractions_ele1 = photon.superCluster()->hitsAndFractions();
	_nHitsSCEle[index] = hitsAndFractions_ele1.size();
	TreeSetExtraCalibVar(hitsAndFractions_ele1, index, photon.superCluster()->seed()->seed().subdetId() == EcalBarrel);
	return;
}


void ZNtupleDumper::TreeSetExtraCalibVar(const pat::Muon& muon1, int index)
{
	if (index < 0) return;
	return;
}


//#============================== extra study tree
void ZNtupleDumper::InitExtraStudyTree()
{
	//  _tree = new TTree("selected",fChain->GetTitle());
	std::cout << "[STATUS] InitExtraStudyTree" << std::endl;
	if(_extraStudyTree == NULL) return;

	_extraStudyTree->Branch("runNumber",     &_runNumber,   "runNumber/i");
	_extraStudyTree->Branch("lumiBlock",     &_lumiBlock,   "lumiBlock/s");
	_extraStudyTree->Branch("eventNumber",   &_eventNumber, "eventNumber/l");
	_extraStudyTree->Branch("eventTime",       &_eventTime,     "eventTime/i");

	_extraStudyTree->Branch("LCRecHitSCEle1", &(_LCRecHitSCEle[0]));
	_extraStudyTree->Branch("LCRecHitSCEle2", &(_LCRecHitSCEle[1]));
	_extraStudyTree->Branch("AlphaRecHitSCEle1", &(_AlphaRecHitSCEle[0]));
	_extraStudyTree->Branch("AlphaRecHitSCEle2", &(_AlphaRecHitSCEle[1]));
	_extraStudyTree->Branch("ICRecHitSCEle1", &(_ICRecHitSCEle[0]));
	_extraStudyTree->Branch("ICRecHitSCEle2", &(_ICRecHitSCEle[1]));
	_extraStudyTree->Branch("ampliErrUncalibRecHitSCEle1", &(_ampliErrUncalibRecHitSCEle[0]));
	_extraStudyTree->Branch("ampliErrUncalibRecHitSCEle2", &(_ampliErrUncalibRecHitSCEle[1]));
	_extraStudyTree->Branch("ampliUncalibRecHitSCEle1", &(_ampliUncalibRecHitSCEle[0]));
	_extraStudyTree->Branch("ampliUncalibRecHitSCEle2", &(_ampliUncalibRecHitSCEle[1]));
	_extraStudyTree->Branch("chi2UncalibRecHitSCEle1", &(_chi2UncalibRecHitSCEle[0]));
	_extraStudyTree->Branch("chi2UncalibRecHitSCEle2", &(_chi2UncalibRecHitSCEle[1]));
	_extraStudyTree->Branch("flagsUncalibRecHitSCEle1", &(_flagsUncalibRecHitSCEle[0]));
	_extraStudyTree->Branch("flagsUncalibRecHitSCEle2", &(_flagsUncalibRecHitSCEle[1]));
	_extraStudyTree->Branch("jitterErrUncalibRecHitSCEle1", &(_jitterErrUncalibRecHitSCEle[0]));
	_extraStudyTree->Branch("jitterErrUncalibRecHitSCEle2", &(_jitterErrUncalibRecHitSCEle[1]));
	_extraStudyTree->Branch("jitterUncalibRecHitSCEle1", &(_jitterUncalibRecHitSCEle[0]));
	_extraStudyTree->Branch("jitterUncalibRecHitSCEle2", &(_jitterUncalibRecHitSCEle[1]));
	_extraStudyTree->Branch("ootAmplitudesUncalibRecHitSCEle1", &(_ootAmplisUncalibRecHitSCEle[0]));
	_extraStudyTree->Branch("ootAmplitudesUncalibRecHitSCEle2", &(_ootAmplisUncalibRecHitSCEle[1]));
	_extraStudyTree->Branch("pedUncalibRecHitSCEle1", &(_pedEUncalibRecHitSCEle[0]));
	_extraStudyTree->Branch("pedUncalibRecHitSCEle2", &(_pedEUncalibRecHitSCEle[1]));

	return;
}


void ZNtupleDumper::ResetExtraStudyVar()
{
	for (int i = 0; i < NELE; ++i) {
		_LCRecHitSCEle[i].clear();
		_AlphaRecHitSCEle[i].clear();
		_ICRecHitSCEle[i].clear();
		_ootAmplisUncalibRecHitSCEle[i].clear();
		_ampliUncalibRecHitSCEle[i].clear();
		_ampliErrUncalibRecHitSCEle[i].clear();
		_pedEUncalibRecHitSCEle[i].clear();
		_jitterUncalibRecHitSCEle[i].clear();
		_jitterErrUncalibRecHitSCEle[i].clear();
		_chi2UncalibRecHitSCEle[i].clear();
		_flagsUncalibRecHitSCEle[i].clear();
	}
}


//#============================== Ele ID tree
void ZNtupleDumper::InitEleIDTree()
{
	//  _tree = new TTree("selected",fChain->GetTitle());
	std::cout << "[STATUS] InitEleIDTree" << std::endl;
	if(_eleIDTree == NULL) {
		return;
	}

	_eleIDTree->Branch("runNumber",     &_runNumber,   "runNumber/i");
	_eleIDTree->Branch("lumiBlock",     &_lumiBlock,   "lumiBlock/s");
	_eleIDTree->Branch("eventNumber",   &_eventNumber, "eventNumber/l");
	_eleIDTree->Branch("eventTime",       &_eventTime,     "eventTime/i");

	_eleIDTree->Branch("dr03TkSumPt", _dr03TkSumPt, "dr03TkSumPt[3]/F");
	_eleIDTree->Branch("dr03EcalRecHitSumEt", _dr03EcalRecHitSumEt, "dr03EcalRecHitSumEt[3]/F");
	_eleIDTree->Branch("dr03HcalTowerSumEt", _dr03HcalTowerSumEt, "dr03HcalTowerSumEt[3]/F");
	_eleIDTree->Branch("sigmaIEtaIEtaSCEle", _sigmaIEtaIEtaSCEle, "sigmaIEtaIEtaSCEle[3]/F");

	_eleIDTree->Branch("E1x5",    _E1x5,    "E1x5[3]/F");
	_eleIDTree->Branch("E2x5Max", _E2x5Max, "E2x5Max[3]/F");
	_eleIDTree->Branch("E1x3",    _E1x3,    "E1x3[3]/F");
	_eleIDTree->Branch("E2x2",    _E2x2,    "E2x2[3]/F");
	_eleIDTree->Branch("S4", _S4, "S4[3]/F");
	_eleIDTree->Branch("etaWidth", _etaWidth, "etaWidth[3]/F");
	_eleIDTree->Branch("phiWidth", _phiWidth, "phiWidth[3]/F");
	//  _eleIDTree->Branch("sigmaIPhiIPhiSCEle", _sigmaIPhiIPhiSCEle, "sigmaIPhiIPhiSCEle[3]/F");
	_eleIDTree->Branch("deltaEtaSuperClusterTrackAtVtx", _deltaEtaSuperClusterTrackAtVtx, "deltaEtaSuperClusterTrackAtVtx[3]/F");
	_eleIDTree->Branch("deltaPhiSuperClusterTrackAtVtx", _deltaPhiSuperClusterTrackAtVtx, "deltaPhiSuperClusterTrackAtVtx[3]/F");
	_eleIDTree->Branch("hOverE", _hOverE, "hOverE[3]/F");
	_eleIDTree->Branch("pfMVA", _pfMVA, "pfMVA[3]/F");
	_eleIDTree->Branch("hasMatchedConversion", _hasMatchedConversion, "hasMatchedConversion[3]/b");
	_eleIDTree->Branch("maxNumberOfExpectedMissingHits", _maxNumberOfExpectedMissingHits, "maxNumberOfExpectedMissingHits[3]/I");
	return;
}


void ZNtupleDumper::ResetEleIDVar()
{
	for (int i = 0; i < NELE; ++i) {
		_sigmaIEtaIEtaSCEle[i]             = initSingleFloat;
		_sigmaIPhiIPhiSCEle[i]             = initSingleFloat;
		_hOverE[i]                         = initSingleFloat;
		_hOverEBC[i]                       = initSingleFloat;
		_dr03TkSumPt[i]                    = initSingleFloat;
		_dr03EcalRecHitSumEt[i]            = initSingleFloat;
		_dr03HcalTowerSumEt[i]             = initSingleFloat;
		_deltaPhiSuperClusterTrackAtVtx[i] = initSingleFloat;
		_deltaEtaSuperClusterTrackAtVtx[i] = initSingleFloat;
		_E1x5[i]                           = initSingleFloat;
		_E1x3[i]                           = initSingleFloat;
		_E2x2[i]                           = initSingleFloat;
		_E2x5Max[i]                        = initSingleFloat;
		_S4[i]                             = initSingleFloat;
		_etaWidth[i]                       = initSingleFloat;
		_phiWidth[i]                       = initSingleFloat;
		_maxNumberOfExpectedMissingHits[i] = initSingleFloat;
		_pfMVA[i]                          = initSingleFloat;
		_hasMatchedConversion[i]           = initSingleInt;
		_eleIDloose[i]                     = initSingleInt;
		_eleIDmedium[i]                    = initSingleInt;
		_eleIDtight[i]                     = initSingleInt;
	}
}


void ZNtupleDumper::TreeSetEleIDVar(const pat::Electron& electron1, const pat::Electron& electron2)
{
	TreeSetEleIDVar(electron1, 0);
	TreeSetEleIDVar(electron2, 1);
	return;
}


void ZNtupleDumper::TreeSetEleIDVar(const pat::Electron& ele, int index)
{
	if(index < 0) {
		_hOverE[-index] = -1;
		return;
	}

	_dr03TkSumPt[index]					  = ele.dr03TkSumPt();
	_dr03EcalRecHitSumEt[index]			  = ele.dr03EcalRecHitSumEt();
	_dr03HcalTowerSumEt[index]			  = ele.dr03HcalTowerSumEt();
	_sigmaIEtaIEtaSCEle[index]			  = ele.full5x5_sigmaIetaIeta();
	_sigmaIPhiIPhiSCEle[index]			  = ele.full5x5_sigmaIphiIphi();
	_E1x5[index]							  = ele.full5x5_e1x5();
	_E2x5Max[index]						  = ele.full5x5_e2x5Max();
	_deltaPhiSuperClusterTrackAtVtx[index] = ele.deltaPhiSuperClusterTrackAtVtx();
	_deltaEtaSuperClusterTrackAtVtx[index] = ele.deltaEtaSuperClusterTrackAtVtx();
	_hOverE[index]						  = ele.full5x5_hcalOverEcal();
	_hOverEBC[index]                       = ele.full5x5_hcalOverEcalBc();

	_pfMVA[index]   = ele.mva_e_pi();
	_hasMatchedConversion[index] = ConversionTools::hasMatchedConversion(ele, _conversionsHandle, _bsHandle->position());
	_maxNumberOfExpectedMissingHits[index] = ele.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);

	const reco::CaloClusterPtr seed_clu = ele.superCluster()->seed();
	_E1x3[index] = _clustertools->e1x3( *seed_clu );
	_E2x2[index] = _clustertools->e2x2( *seed_clu );
	_S4[index]   = _clustertools->e2x2( *seed_clu ) / _clustertools->e5x5( *seed_clu );
	_etaWidth[index] = ele.superCluster()->etaWidth();
	_phiWidth[index] = ele.superCluster()->phiWidth();

	return;
}


void ZNtupleDumper::TreeSetEleIDVar(const pat::Photon& photon, const pat::Muon& muon1, const pat::Muon& muon2)
{
	TreeSetEleIDVar(photon, 0);
	TreeSetEleIDVar(muon1, -1);
	TreeSetEleIDVar(muon2, -2);
	return;
}


void ZNtupleDumper::TreeSetEleIDVar(const pat::Photon& photon, int index)
{
	if(index < 0) {
		return;
	}

	_sigmaIEtaIEtaSCEle[index]  = photon.sigmaIetaIeta(); // alcarereco
	_hOverE[index] = photon.hadronicOverEm();

	const reco::SuperCluster photonSC = *(photon.superCluster());
	_hasMatchedConversion[index] = ConversionTools::hasMatchedConversion(photonSC, _conversionsHandle, _bsHandle->position());

	// _eleIDloose[index]  = photon.photonID("loose");
	// _eleIDmedium[index] = photon.photonID("medium");
	// _eleIDtight[index]  = photon.photonID("tight");
	return;
}


void ZNtupleDumper::TreeSetEleIDVar(const pat::Muon& muon1, int index)
{
	if(index < 0) {
		return;
	}

	//  _eleIDloose[index]  = muon1.muonID("loose");
	//  _eleIDmedium[index] = muon1.muonID("medium");
	//  _eleIDtight[index]  = muon1.muonID("tight");
	return;
}


//#============================== Track study tree
void ZNtupleDumper::InitTrackTree()
{
	std::cout << "[STATUS] InitTrackTree" << std::endl;
	if(_trackTree == NULL) return;

	_trackTree->Branch("runNumber",     &_runNumber,   "runNumber/i");
	_trackTree->Branch("lumiBlock",     &_lumiBlock,   "lumiBlock/s");
	_trackTree->Branch("eventNumber",   &_eventNumber, "eventNumber/l");
	_trackTree->Branch("eventTime",     &_eventTime,     "eventTime/i");

	_trackTree->Branch("gsfTrackLengthFromOuterP", _gsfTrackLengthFromOuterP, "gsfTrackLengthFromOuterP[3]/F");
	_trackTree->Branch("gsfTrackLengthFromEleP", _gsfTrackLengthFromEleP, "gsfTrackLengthFromEleP[3]/F");
	_trackTree->Branch("gsfTrackOuterPt",  _gsfTrackOuterPt,  "gsfTrackOuterPt[3]/F");
	_trackTree->Branch("gsfTrackOuterEta", _gsfTrackOuterEta, "gsfTrackOuterEta[3]/F");
	_trackTree->Branch("gsfTrackOuterPhi", _gsfTrackOuterPhi, "gsfTrackOuterPhi[3]/F");
	_trackTree->Branch("gsfTrackVtxPt",    _gsfTrackVtxPt,    "gsfTrackVtxPt[3]/F");
	_trackTree->Branch("gsfTrackVtxEta",   _gsfTrackVtxEta,   "gsfTrackVtxEta[3]/F");
	_trackTree->Branch("gsfTrackVtxPhi",   _gsfTrackVtxPhi,   "gsfTrackVtxPhi[3]/F");
	_trackTree->Branch("gsfTrackVtxX",     _gsfTrackVtxX,     "gsfTrackVtxX[3]/F");
	_trackTree->Branch("gsfTrackVtxY",     _gsfTrackVtxY,     "gsfTrackVtxY[3]/F");
	_trackTree->Branch("gsfTrackVtxZ",     _gsfTrackVtxZ,     "gsfTrackVtxZ[3]/F");
	_trackTree->Branch("gsfTrackCaloX",    _gsfTrackCaloX,    "gsfTrackCaloX[3]/F");
	_trackTree->Branch("gsfTrackCaloY",    _gsfTrackCaloY,    "gsfTrackCaloY[3]/F");
	_trackTree->Branch("gsfTrackCaloZ",    _gsfTrackCaloZ,    "gsfTrackCaloZ[3]/F");
	_trackTree->Branch("gsfTrackTangentsPt", &(_gsfTrackTangentsPt[0]));
	_trackTree->Branch("gsfTrackTangentsEta", &(_gsfTrackTangentsEta[0]));
	_trackTree->Branch("gsfTrackTangentsPhi", &(_gsfTrackTangentsPhi[0]));
	_trackTree->Branch("gsfTrackTangentsDeltaP", &(_gsfTrackTangentsDeltaP[0]));
	_trackTree->Branch("gsfTrackTangentsDeltaPErr", &(_gsfTrackTangentsDeltaPErr[0]));
	_trackTree->Branch("gsfTrackTangentsX", &(_gsfTrackTangentsX[0]));
	_trackTree->Branch("gsfTrackTangentsY", &(_gsfTrackTangentsY[0]));
	_trackTree->Branch("gsfTrackTangentsZ", &(_gsfTrackTangentsZ[0]));
	return;
}


void ZNtupleDumper::ResetTrackVar()
{
	for (int i = 0; i < NELE; ++i) {
		_gsfTrackLengthFromOuterP[i] = initSingleFloat;
		_gsfTrackLengthFromVtxP[i] = initSingleFloat;
		_gsfTrackLengthFromEleP[i] = initSingleFloat;
		_gsfTrackLengthFromTangents[i] = initSingleFloat;
		_gsfTrackOuterPt[i] = initSingleFloat;
		_gsfTrackOuterEta[i] = initSingleFloat;
		_gsfTrackOuterPhi[i] = initSingleFloat;
		_gsfTrackVtxPt[i] = initSingleFloat;
		_gsfTrackVtxEta[i] = initSingleFloat;
		_gsfTrackVtxPhi[i] = initSingleFloat;
		_gsfTrackVtxX[i] = initSingleFloat;
		_gsfTrackVtxY[i] = initSingleFloat;
		_gsfTrackVtxZ[i] = initSingleFloat;
		_gsfTrackCaloX[i] = initSingleFloat;
		_gsfTrackCaloY[i] = initSingleFloat;
		_gsfTrackCaloZ[i] = initSingleFloat;
		_gsfTrackTangentsPt[i].clear();
		_gsfTrackTangentsEta[i].clear();
		_gsfTrackTangentsPhi[i].clear();
		_gsfTrackTangentsDeltaP[i].clear();
		_gsfTrackTangentsDeltaPErr[i].clear();
		_gsfTrackTangentsX[i].clear();
		_gsfTrackTangentsY[i].clear();
		_gsfTrackTangentsZ[i].clear();
	}
}


// method to get the raw energy of one plane of ES summing the energy of only recHits associated to the electron SC
///\todo highly inefficient: instead of the loop over the recHits should use a ->find() method, it should return both energies of both planes
float ZNtupleDumper::GetESPlaneRawEnergy(const reco::SuperCluster& sc, unsigned int planeIndex) const
{
	double RawenergyPlane = 0.;
	double pfRawenergyPlane = 0.;
//	if(!_ESRechitsHandle.isValid())
//		return RawenergyPlane;
	if (!sc.preshowerClusters().isAvailable()) //protection for miniAOD
		return initSingleFloat;

	for(auto iES = sc.preshowerClustersBegin(); iES != sc.preshowerClustersEnd(); ++iES) {//loop over preshower clusters
		const std::vector< std::pair<DetId, float> > hits = (*iES)->hitsAndFractions();
		for(std::vector<std::pair<DetId, float> >::const_iterator rh = hits.begin(); rh != hits.end(); ++rh) { // loop over recHits of the cluster
			//      std::cout << "print = " << (*iES)->printHitAndFraction(iCount);
			//      ++iCount;
			for(ESRecHitCollection::const_iterator esItr = _ESRechitsHandle->begin(); esItr != _ESRechitsHandle->end(); ++esItr) {//loop over ES rechits to find the one in the cluster
				ESDetId rhid = ESDetId(esItr->id());
				if(rhid == (*rh).first) { // found ESrechit
					// std::cout << " ES energy = " << esItr->energy() << " pf energy = " << (*rh).second << std::endl;
					if((int) rhid.plane() == (int) planeIndex) {
						RawenergyPlane += esItr->energy();
						pfRawenergyPlane += rh->second;
					}
					break;
				}
			}
		}
	}

	if (pfRawenergyPlane) ; // avoid compilation error for unused var
	if (RawenergyPlane);
	//std::cout << "LC DEBUG RawenergyPlane "<< RawenergyPlane << ", pfRawenergyPlane " << pfRawenergyPlane << std::endl;
	return RawenergyPlane;
}


std::pair<DetId, float> ZNtupleDumper::findEnergySortedHit(const reco::SuperCluster& cluster, const EcalRecHitCollection * recHits, size_t rank)
{
	std::vector< std::pair<DetId, float> > hitsFractions = cluster.hitsAndFractions();
	std::vector< std::pair<DetId, float> > ordered;
	for(auto & it : hitsFractions) {
		auto rh = recHits->find(it.first);
		if (rh != recHits->end()) {
			ordered.push_back(std::make_pair(it.first, it.second * rh->energy()));
		} else {
			std::cerr << "[ERROR findEnergySortedHit] RecHit not found for DetID: " << it.first.rawId() << std::endl;
		}
	}
	std::sort(ordered.begin(), ordered.end(), [](auto & a, auto & b) {
		return a.second > b.second;
	});
#ifdef DEBUG
	for (auto & it : ordered) {
		std::cout << "[DEBUG findEnergySortedHit]  DetId: " << it.first.rawId() << " energy: " << it.second << "\n";
	}
#endif
	size_t s = ordered.size();
	if (rank >= s) {
		std::cerr << "[ERROR findEnergySortedHit] Requesting rank " << rank << " for a SC with only " << s << " crystal(s).\n";
		return std::make_pair(DetId(0), 0);
	}
	return ordered[rank];
}


bool ZNtupleDumper::elePreselection(const pat::Electron& electron) const
{
	if(electron.et() < 10) return false; // minimum Et cut in preselection

	//to make alcareco/alcarereco ntuples coeherent
	if(!electron.ecalDriven()) return false;

	//if(eleIter1->parentSuperCluster().isNull()) continue;
	if(_presel) {
		if(! (electron.electronID(_eleID_loose))) return false;
	}
	return true;
}


//define this as a plug-in
DEFINE_FWK_MODULE(ZNtupleDumper);

//  LocalWords:  pileupInfoTAG conversionsProducerTAG triggerResultsTAG
