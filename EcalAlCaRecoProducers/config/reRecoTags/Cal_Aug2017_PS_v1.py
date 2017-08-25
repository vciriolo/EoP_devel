import FWCore.ParameterSet.Config as cms


from CondCore.ESSources.CondDBESSource_cfi import * 
RerecoGlobalTag = GlobalTag.clone(
    globaltag = cms.string('92X_dataRun2_Prompt_v8'),
    #toGet = cms.VPSet( ),   # hook to override or add single payloads
    toGet = cms.VPSet(
        cms.PSet(record = cms.string("EcalPulseShapesRcd"),
                 tag = cms.string("EcalPulseShapes_October2017_rereco_v1"),
                 connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
                 ),
        ),
    )
