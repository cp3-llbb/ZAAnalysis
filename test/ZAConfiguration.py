
import FWCore.ParameterSet.Config as cms

process = cms.Process("ZAAna")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.GlobalTag.globaltag = "MCRUN2_74_V9"

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.source = cms.Source("PoolSource")
process.source.fileNames = cms.untracked.vstring(
        'file:///home/fynu/sbrochet/storage/MINIAODSIM/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9_reduced.root'
        )
##### #####
# Electron ID recipe
##### #####
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff','RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff', 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']
#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

##### #####
# Services
##### #####
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Producers
from cp3_llbb.Framework import EventProducer
from cp3_llbb.Framework import GenParticlesProducer
from cp3_llbb.Framework import HLTProducer
from cp3_llbb.Framework import JetsProducer
from cp3_llbb.Framework import METProducer
from cp3_llbb.Framework import MuonsProducer
from cp3_llbb.Framework import ElectronsProducer
from cp3_llbb.Framework import VerticesProducer

process.framework = cms.EDProducer("ExTreeMaker",
        output = cms.string('output_mc.root'),

        producers = cms.PSet(

            event = EventProducer.default_configuration,

            hlt = HLTProducer.default_configuration,

            gen_particles = GenParticlesProducer.default_configuration,

            jets = JetsProducer.default_configuration.clone(
                parameters = cms.PSet(
                    cut = cms.untracked.string("pt > 25")
                    )
                ),

            met = METProducer.default_configuration,

            muons = MuonsProducer.default_configuration.clone(
                parameters = cms.PSet(
                    cut = cms.untracked.string("pt > 20 && isMediumMuon == 1"),
                    src = cms.untracked.InputTag('slimmedMuons'),
                    ea_R03 = cms.untracked.FileInPath('cp3_llbb/Framework/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons.txt'),
                    ea_R04 = cms.untracked.FileInPath('cp3_llbb/Framework/data/effAreaMuons_cone04_pfNeuHadronsAndPhotons.txt'),
                    scale_factors = cms.untracked.PSet(
                        id_loose = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_Loose_id.json'),
                        id_soft = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_Soft_id.json'),
                        id_tight = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_Tight_id.json'),

                        iso_012_tight = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_combRelIsoPF04dBeta<012_Tight_iso.json'),
                        iso_02_loose = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_combRelIsoPF04dBeta<02_Loose_iso.json'),
                        iso_02_tight = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Muon_combRelIsoPF04dBeta<02_Tight_iso.json'),
                        )
                    )
                ),

            electrons = ElectronsProducer.default_configuration.clone(
                parameters = cms.PSet(
                    cut = cms.untracked.string("pt > 20"),
                    src = cms.untracked.InputTag('slimmedElectrons'),
                    ea_R03 = cms.untracked.FileInPath('RecoEgamma/ElectronIdentification/data/PHYS14/effAreaElectrons_cone03_pfNeuHadronsAndPhotons.txt'),
                    ea_R04 = cms.untracked.FileInPath('cp3_llbb/Framework/data/effAreaElectrons_cone04_pfNeuHadronsAndPhotons.txt'),
                    ids = cms.untracked.VInputTag(
                        'egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80',
                        'egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90',
                        'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose',
                        'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium',
                        'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight',
                        'egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto',
                        'egmGsfElectronIDs:heepElectronID-HEEPV51'
                    ),
                    scale_factors = cms.untracked.PSet(
                        id_veto = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_sfVETO_id.json'),
                        id_loose = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_sfLOOSE_id.json'),
                        id_medium = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_sfMEDIUM_id.json'),
                        id_tight = cms.untracked.FileInPath('cp3_llbb/Framework/data/ScaleFactors/Electron_sfTIGHT_id.json'),
                        )
                    )
                ),

            vertices = VerticesProducer.default_configuration,
            ),

        analyzers = cms.PSet(
            za = cms.PSet(
                type = cms.string('za_analyzer'),
                prefix = cms.string('za_'),
                enable = cms.bool(True)
                )
            )
        )

process.p = cms.Path(
        process.egmGsfElectronIDSequence *
        process.framework
        )


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Output and Log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.options   = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
process.options.allowUnscheduled = cms.untracked.bool(True)
