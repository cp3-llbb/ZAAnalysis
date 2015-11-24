#pragma once

#include <string>
#include <utility>
#include <vector>
#include <limits>


#include <cp3_llbb/Framework/interface/Analyzer.h>
#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/FatJetsProducer.h>


#include <cp3_llbb/ZAAnalysis/interface/ZATypes.h>
#include <cp3_llbb/ZAAnalysis/interface/Tools.h>

class ZAAnalyzer: public Framework::Analyzer {
    public:
        ZAAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config),

            m_electronPtCut( config.getUntrackedParameter<double>("electronPtCut", 20) ),
            m_electronEtaCut( config.getUntrackedParameter<double>("electronEtaCut", 2.5) ),
            m_electronVetoIDName( config.getUntrackedParameter<std::string>("electronVetoIDName") ),
            m_electronLooseIDName( config.getUntrackedParameter<std::string>("electronLooseIDName") ),
            m_electronMediumIDName( config.getUntrackedParameter<std::string>("electronMediumIDName") ),
            m_electronTightIDName( config.getUntrackedParameter<std::string>("electronTightIDName") ),
//            m_electronSelectionID( config.getUntrackedParameter<std::string>("electronSelectionID") ),
            
            m_muonPtCut( config.getUntrackedParameter<double>("muonPtCut", 20) ),
            m_muonEtaCut( config.getUntrackedParameter<double>("muonEtaCut", 2.4) ),
//            m_muonSelectionID( config.getUntrackedParameter<std::string>("muonSelectionID","tight") ),
            m_muonLooseIsoCut( config.getUntrackedParameter<double>("muonLooseIsoCut", 0.25) ),
            m_muonTightIsoCut( config.getUntrackedParameter<double>("muonTightIsoCut", 0.15) ),
//            m_muonSelectionIsoCut( config.getUntrackedParameter<double>("muonSelectionIsoCut", 0.12) ),
            
            m_jetPtCut( config.getUntrackedParameter<double>("jetPtCut", 30) ),
            m_jetEtaCut( config.getUntrackedParameter<double>("jetEtaCut", 2.5) ),
            m_jetPUID( config.getUntrackedParameter<double>("jetPUID", std::numeric_limits<float>::min()) ),
            m_jetDRleptonCut( config.getUntrackedParameter<double>("jetDRleptonCut", 0.3) ),
            m_jetID( config.getUntrackedParameter<std::string>("jetID", "tight") ),
            m_jetCSVv2Name( config.getUntrackedParameter<std::string>("jetCSVv2Name", "pfCombinedInclusiveSecondaryVertexV2BJetTags") ),
            m_jetCSVv2L( config.getUntrackedParameter<double>("jetCSVv2L", 0.605) ),
            m_jetCSVv2M( config.getUntrackedParameter<double>("jetCSVv2M", 0.89) ),
            m_jetCSVv2T( config.getUntrackedParameter<double>("jetCSVv2T", 0.97) ),
            
            m_hltDRCut( config.getUntrackedParameter<double>("hltDRCut", std::numeric_limits<float>::max()) ),
            m_hltDPtCut( config.getUntrackedParameter<double>("hltDPtCut", std::numeric_limits<float>::max()) )
        {
        }
        

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const AnalyzersManager&, const CategoryManager&) override;
        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet&) override;


        BRANCH(leptons, std::vector<ZAAnalysis::Lepton>);
        BRANCH(isolatedElectrons, std::vector<ZAAnalysis::Lepton>);
        BRANCH(isolatedMuons, std::vector<ZAAnalysis::Lepton>);
        BRANCH(vetoLeptons, std::vector<ZAAnalysis::Lepton>);
        BRANCH(selJets, std::vector<ZAAnalysis::Jet>);
        BRANCH(selFatJets, std::vector<ZAAnalysis::Jet>);
        BRANCH(selBjetsM, std::vector<ZAAnalysis::Jet>);

        BRANCH(dilep_ptOrdered, std::vector<ZAAnalysis::Lepton>);  //FIXME remove the vectors and define "clearing mechanism", or use pair?
        BRANCH(dijet_ptOrdered, std::vector<ZAAnalysis::Jet>);
        BRANCH(dijet_CSVv2Ordered, std::vector<ZAAnalysis::Jet>);

        BRANCH(diLeptons, std::vector<ZAAnalysis::DiLepton>);
        BRANCH(diJets, std::vector<ZAAnalysis::DiJet>);
        BRANCH(diLepDiJets, std::vector<ZAAnalysis::DiLepDiJet>);
        //Branch(diLepFatJets, std::vector<ZAAnalysis::DiLepFatJet>);

    private:

        const float m_electronPtCut, m_electronEtaCut;
        const std::string m_electronVetoIDName;
        const std::string m_electronLooseIDName;
        const std::string m_electronMediumIDName;
        const std::string m_electronTightIDName;
        const std::string m_electronSelectionID;

//        const std::string m_muonSelectionID;
        const float m_muonPtCut, m_muonEtaCut, m_muonLooseIsoCut, m_muonTightIsoCut;

        const float m_jetPtCut, m_jetEtaCut, m_jetPUID, m_jetDRleptonCut;
        const std::string m_jetID, m_jetCSVv2Name;
        const float m_jetCSVv2L, m_jetCSVv2M, m_jetCSVv2T;

        const float m_hltDRCut, m_hltDPtCut;

        static inline bool muonIDAccessor(const MuonsProducer& muons, const uint16_t index, const std::string& muonID){
            if(index >= muons.p4.size())
              throw edm::Exception(edm::errors::StdException, "Invalid muon index passed to ID accessor");

            if(muonID == "loose")
                return muons.isLoose[index];
            
            if(muonID == "medium")
                return muons.isMedium[index];
            
            if(muonID == "tight")
                return muons.isTight[index];
            
            throw edm::Exception(edm::errors::NotFound, "Unknown muonID passed to analyzer");
        }
        
        static inline bool jetIDAccessor(const JetsProducer& jets, const uint16_t index, const std::string& jetID){
            if(index >= jets.p4.size())
              throw edm::Exception(edm::errors::StdException, "Invalid jet index passed to ID accessor");
            
            if(jetID == "loose")
                return jets.passLooseID[index];
            
            if(jetID == "tight")
                return jets.passTightID[index];
            
            if(jetID == "tightLeptonVeto")
                return jets.passTightLeptonVetoID[index];
            
            throw edm::Exception(edm::errors::NotFound, "Unknown jetID passed to analyzer");
        }

};
