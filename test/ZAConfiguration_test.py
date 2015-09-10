import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework

process = Framework.create(None, 'MCRUN2_74_V9', cms.PSet(
    za = cms.PSet(
        type = cms.string('za_analyzer'),
        prefix = cms.string('za_'),
        enable = cms.bool(True),
        categories_parameters = cms.PSet(
            mll_cut = cms.untracked.double(20)
            )
        ),

    test = cms.PSet(
        type = cms.string('test_analyzer'),
        prefix = cms.string('test_'),
        enable = cms.bool(True)
        )
    )
    )

Framework.schedule(process, ['za', 'test'])

process.source.fileNames = cms.untracked.vstring(
        'file:///home/fynu/sbrochet/storage/MINIAODSIM/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9_reduced.root'
        )

