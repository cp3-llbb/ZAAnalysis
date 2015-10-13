import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework

process = Framework.create(True, eras.Run2_50ns, '74X_dataRun2_v2', redoJEC=True)

process.framework.analyzers.dilepton = cms.PSet(
        type = cms.string('dilepton_analyzer'),
        prefix = cms.string('dilepton_'),
        enable = cms.bool(True),
        categories_parameters = cms.PSet(
            mll_cut = cms.untracked.double(20)
            )
        )

process.framework.analyzers.za = cms.PSet(
        type = cms.string('za_analyzer'),
        prefix = cms.string('za_'),
        enable = cms.bool(True),
        categories_parameters = cms.PSet(
            mll_cut = cms.untracked.double(20)
            )
        )
    
    

Framework.schedule(process, ['dilepton','za'])

process.source.fileNames = cms.untracked.vstring(
        '/store/data/Run2015B/DoubleMuon/MINIAOD/17Jul2015-v1/30000/D8ED75E7-C12E-E511-8CBF-0025905A608C.root'
        )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
