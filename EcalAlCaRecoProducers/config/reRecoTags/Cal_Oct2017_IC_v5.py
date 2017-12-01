import FWCore.ParameterSet.Config as cms

from CondCore.DBCommon.CondDBSetup_cfi import *

RerecoGlobalTag = cms.ESSource("PoolDBESSource",
                               CondDBSetup,
                               connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
                               globaltag = cms.string('92X_dataRun2_2017Repro_Candidate_2017_11_10_15_04_54'),
                               snapshotTime = cms.string("9999-12-31 23:59:59.000"),
                               toGet = cms.VPSet(
        cms.PSet(record = cms.string("EcalPulseShapesRcd"),
                 tag = cms.string("EcalPulseShapes_Run2017F"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                 ),
        cms.PSet(record = cms.string("EcalPedestalsRcd"),
                 tag = cms.string("EcalPedestals_Legacy2017_v1"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                 ),
        cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
                 tag = cms.string("EcalIntercalibConstants_eopPNEB_etaScalePNEE_v5_RunFonly"),
                 connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS"),
                 ),
        ),
                               )
