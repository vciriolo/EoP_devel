import FWCore.ParameterSet.Config as cms

from CondCore.DBCommon.CondDBSetup_cfi import *

RerecoGlobalTag = cms.ESSource("PoolDBESSource",
                               CondDBSetup,
                               connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
                               globaltag = cms.string('92X_dataRun2_Prompt_v9'),
                               toGet = cms.VPSet(
        cms.PSet(record = cms.string("EcalPedestalsRcd"),
                 tag = cms.string("EcalPedestals_Legacy2017_v1"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                 ),
        # cms.PSet(record = cms.string("EcalPulseShapesRcd"),
        #          tag = cms.string("EcalPulseShapes_October2017_rereco_v1"),
        #          connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        #          ),
        cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
                 tag = cms.string("EcalIntercalibConstants_2017_2015_at_high_eta"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                 ),
        cms.PSet(record = cms.string("EcalADCToGeVConstantRcd"),
                 tag = cms.string("EcalADCToGeVConstant_plus_2.4prct_in_EE"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                 ),
        ),
                               )
