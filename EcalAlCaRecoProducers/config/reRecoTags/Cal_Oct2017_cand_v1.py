import FWCore.ParameterSet.Config as cms

from CondCore.DBCommon.CondDBSetup_cfi import *

RerecoGlobalTag = cms.ESSource("PoolDBESSource",
                               CondDBSetup,
                               connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
                               globaltag = cms.string('92X_dataRun2_2017Repro_Candidate_2017_10_31_16_58_36'),
                               toGet = cms.VPSet( 
        cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
                 tag = cms.string("EcalIntercalibConstants_Cal_Oct2017_PNcorrection_eop_v1"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                 ),
        ),
                               )
