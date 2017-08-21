import FWCore.ParameterSet.Config as cms

from CondCore.ESSources.CondDBESSource_cfi import * 
#CondDBConnection.connect = cms.string( 'frontier://FrontierProd/CMS_CONDITIONS' )
RerecoGlobalTag = GlobalTag.clone(
    globaltag = cms.string('92X_dataRun2_Prompt_v5'),
    toGet = cms.VPSet(
        cms.PSet(record = cms.string("ESIntercalibConstantsRcd"),
                 tag = cms.string("ESIntercalibConstants_LG_Run2017A"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                 ),
        # cms.PSet(record = cms.string("ESEEIntercalibConstantsRcd"),
        #          tag = cms.string("ESEEIntercalibConstants_LG_offline_data_default"),
        #          connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        #          ),
        ),
    )
