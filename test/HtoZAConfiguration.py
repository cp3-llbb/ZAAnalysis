
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework
from cp3_llbb.Framework import METProducer
from cp3_llbb.Framework.CmdLine import CmdLine

options = CmdLine()
runOnData = options.runOnData == 1

globalTag_ = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
processName_ = 'PAT'
if runOnData :
    globalTag_ = '80X_dataRun2_2016SeptRepro_v7'
    processName_ = 'RECO'

# WARNING: prescale factor for data should be set to -1
options = CmdLine(defaults=dict(runOnData=0, era="25ns", globalTag='80X_dataRun2_2016SeptRepro_v7'))

framework = Framework.Framework(options)

framework.addAnalyzer('hZA_analyzer', cms.PSet(
        type = cms.string('hZA_analyzer'),
        prefix = cms.string('hZA_'),
        enable = cms.bool(True),
        categories_parameters = cms.PSet(
            # Per-category lepton pt cuts
            mumu_leadingLeptonPtCut = cms.untracked.double(20), # muon
            mumu_subleadingLeptonPtCut = cms.untracked.double(10), # muon
            elel_leadingLeptonPtCut = cms.untracked.double(25), # electron
            elel_subleadingLeptonPtCut = cms.untracked.double(15), # electron
            muel_leadingLeptonPtCut = cms.untracked.double(25), # muon
            muel_subleadingLeptonPtCut = cms.untracked.double(15), # electron
            elmu_leadingLeptonPtCut = cms.untracked.double(25), # electron
            elmu_subleadingLeptonPtCut = cms.untracked.double(10), # muon
        ),
        parameters = cms.PSet(
            prescaleFactor = cms.untracked.int32(-1),
            # Producers
            electronsProducer = cms.string('electrons'),
            muonsProducer = cms.string('muons'),
            jetsProducer = cms.string('jets'),
            metProducer = cms.string('met'),
            nohfMETProducer = cms.string('nohf_met'),

            # Pre-selection pt cut, applied to all leptons
            leadingElectronPtCut = cms.untracked.double(25),
            subleadingElectronPtCut = cms.untracked.double(15),
            leadingMuonPtCut = cms.untracked.double(20),
            subleadingMuonPtCut = cms.untracked.double(10),

            electronEtaCut = cms.untracked.double(2.5),
            muonLooseIsoCut = cms.untracked.double(.25), # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO 
            muonTightIsoCut = cms.untracked.double(.15), # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO 
            muonEtaCut = cms.untracked.double(2.4),
            electrons_loose_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-loose"), # https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#HLT_safe_selection_for_2016_data
            electrons_medium_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-medium"),
            electrons_tight_wp_name = cms.untracked.string("cutBasedElectronID-Summer16-80X-V1-tight"),
            electrons_hlt_safe_wp_name = cms.untracked.string("cutBasedElectronHLTPreselection-Summer16-V1"),
            electrons_mva_wp80_name = cms.untracked.string("mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
            electrons_mva_wp90_name = cms.untracked.string("mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
            electrons_mva_HZZ_loose_wp_name = cms.untracked.string("mvaEleID-Spring16-HZZ-V1-wpLoose"),
            
            jetEtaCut = cms.untracked.double(2.4),
            jetPtCut = cms.untracked.double(20),

            # BTAG INFO
            # Working points from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
            discr_name_cmva = cms.untracked.string("pfCombinedMVAV2BJetTags"),
            discr_cut_cMVAv2_loose =  cms.untracked.double(-0.5884),
            discr_cut_cMVAv2_medium =  cms.untracked.double(0.4432),
            discr_cut_cMVAv2_tight =  cms.untracked.double(0.9432),

            discr_name_deepcsv_probb = cms.untracked.string("pfDeepCSVJetTags:probb"),
            discr_name_deepcsv_probbb = cms.untracked.string("pfDeepCSVJetTags:probbb"),
            discr_cut_deepCSV_loose = cms.untracked.double(0.2219),
            discr_cut_deepCSV_medium = cms.untracked.double(0.6324),
            discr_cut_deepCSV_tight = cms.untracked.double(0.8958),


            minDR_l_j_Cut = cms.untracked.double(0.3),
            hltDRCut = cms.untracked.double(0.1),
            hltDPtCut = cms.untracked.double(0.5),  # cut will be DPt/Pt < hltDPtCut
            applyBJetRegression = cms.untracked.bool(False), # BE SURE TO ACTIVATE computeRegression FLAG BELOW

            hlt_efficiencies = cms.untracked.PSet(

                    IsoMu17leg = cms.untracked.FileInPath('cp3_llbb/ZAAnalysis/data/Efficiencies/Muon_DoubleIsoMu17Mu8_IsoMu17leg.json'),
                    IsoMu8orIsoTkMu8leg = cms.untracked.FileInPath('cp3_llbb/ZAAnalysis/data/Efficiencies/Muon_DoubleIsoMu17TkMu8_IsoMu8legORTkMu8leg.json'),

                    DoubleEleHighPtleg = cms.untracked.FileInPath('cp3_llbb/ZAAnalysis/data/Efficiencies/Electron_IsoEle23Leg.json'),
                    DoubleEleLowPtleg = cms.untracked.FileInPath('cp3_llbb/ZAAnalysis/data/Efficiencies/Electron_IsoEle12Leg.json'),

                    EleMuHighPtleg = cms.untracked.FileInPath('cp3_llbb/ZAAnalysis/data/Efficiencies/Electron_IsoEle23Leg.json'),
                    MuEleLowPtleg = cms.untracked.FileInPath('cp3_llbb/ZAAnalysis/data/Efficiencies/Electron_IsoEle12Leg.json'),

                    IsoMu8leg = cms.untracked.FileInPath('cp3_llbb/ZAAnalysis/data/Efficiencies/Muon_XPathIsoMu8leg.json'),
                    IsoMu23leg = cms.untracked.FileInPath('cp3_llbb/ZAAnalysis/data/Efficiencies/Muon_XPathIsoMu23leg.json'),
            )
        )
    )
)

# Remove fat jets
framework.removeProducer('fat_jets')

framework.getProducer('hlt').parameters.triggers = cms.untracked.FileInPath('cp3_llbb/ZAAnalysis/data/triggers.xml')
# framework.getProducer('jets').parameters.cut = cms.untracked.string("pt > 20")
#framework.getProducer('jets').parameters.computeRegression = cms.untracked.bool(True)

framework.getProducer('electrons').parameters.scale_factors.id_mediumplushltsafe_hh = cms.untracked.FileInPath('cp3_llbb/ZAAnalysis/data/ScaleFactors/Electron_MediumPlusHLTSafeID_moriond17.json')

from cp3_llbb.Framework.JetsProducer import discriminators_deepFlavour
framework.redoJEC(addBtagDiscriminators=discriminators_deepFlavour)

framework.applyMuonCorrection('rochester')

framework.applyElectronRegression()
framework.applyElectronSmearing()

if not runOnData:
    framework.smearJets(resolutionFile='cp3_llbb/Framework/data/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt', scaleFactorFile='cp3_llbb/Framework/data/Spring16_25nsV10_MC_SF_AK4PFchs.txt')
    framework.doSystematics(['jec', 'jer'], jec={'uncertaintiesFile': 'cp3_llbb/ZAAnalysis/data/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt', 'splitBySources': True})

process = framework.create()

if runOnData: 
    process.source.fileNames = cms.untracked.vstring(
            '/store/data/Run2016F/DoubleMuon/MINIAOD/23Sep2016-v1/50000/040EDEBA-0490-E611-A424-008CFA110C68.root'
        )
else: 
    process.framework.treeFlushSize = cms.untracked.uint64(5 * 1024 * 1024)

    process.source.fileNames = cms.untracked.vstring(
	    
            # Signal: H->ZA - mH=1000, mA=200
            # '/store/mc/RunIISummer16MiniAODv2/HToZATo2L2B_MH-1000_MA-200_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/A2BC660D-17DD-E611-9B31-001517F7F950.root'

            # Signal: H->ZA - mH=500, mA=100
            '/store/mc/RunIISummer16MiniAODv2/HToZATo2L2B_MH-500_MA-100_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/9475A8FB-25DF-E611-AFA8-20CF3019DF02.root'

            # Background: TT
            # '/store/mc/RunIISummer16MiniAODv2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/00ED79D3-CFC1-E611-B748-3417EBE64483.root'

            # Background: DY
            # '/store/mc/RunIISummer16MiniAODv2/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/100000/00CEFB4F-C1D2-E611-BBF4-7845C4FC3C11.root'
        )

#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.source.skipEvents = cms.untracked.uint32(10)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
