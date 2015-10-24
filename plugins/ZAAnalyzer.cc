#include <cp3_llbb/ZAAnalysis/interface/Types.h>
#include <cp3_llbb/ZAAnalysis/interface/Tools.h>
//#include <cp3_llbb/ZAAnalysis/interface/GenStatusFlags.h>
#include <cp3_llbb/ZAAnalysis/interface/ZAAnalyzer.h>
#include <cp3_llbb/ZAAnalysis/interface/ZADileptonCategories.h>

#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/METProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>
#include <cp3_llbb/Framework/interface/GenParticlesProducer.h>

#include <Math/PtEtaPhiE4D.h>
#include <Math/LorentzVector.h>
#include <Math/VectorUtil.h>

using namespace ROOT::Math;

using namespace ZAAnalysis;

float ZAAnalysis::DeltaEta(const myLorentzVector& v1, const myLorentzVector& v2){
  return abs(v1.Eta() - v2.Eta());
}


void ZAAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup, const ProducersManager& producers, const AnalyzersManager& analyzers, const CategoryManager&) {

  #ifdef _ZA_DEBUG_
    std::cout << "Begin event." << std::endl;
  #endif

  electrons_IDIso.resize( LepID::Count * LepIso::Count );
  muons_IDIso.resize( LepID::Count * LepIso::Count );
  leptons_IDIso.resize( LepID::Count * LepIso::Count );

  diLeptons_IDIso.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count );
  
  selJets_selID_DRCut.resize( LepID::Count * LepIso::Count );
  selBJets_DRCut_BWP_PtOrdered.resize( LepID::Count * LepIso::Count * BWP::Count );
  selBJets_DRCut_BWP_CSVv2Ordered.resize( LepID::Count * LepIso::Count * BWP::Count );

/*
  diJets_DRCut.resize( LepID::Count * LepIso::Count );
  diBJets_DRCut_BWP_PtOrdered.resize( LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  diBJets_DRCut_BWP_CSVv2Ordered.resize( LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  
  diLepDiJets_DRCut.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count );
  diLepDiBJets_DRCut_BWP_PtOrdered.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  diLepDiBJets_DRCut_BWP_CSVv2Ordered.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
*/  
  ///////////////////////////
  //       ELECTRONS       //
  ///////////////////////////

  #ifdef _ZA_DEBUG_
    std::cout << "Electrons" << std::endl;
  #endif

  const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");

  for(uint16_t ielectron = 0; ielectron < electrons.p4.size(); ielectron++){
    if( electrons.p4[ielectron].Pt() > m_electronPtCut && abs(electrons.p4[ielectron].Eta()) < m_electronEtaCut ){
      
      Lepton m_lepton(
          electrons.p4[ielectron], 
          ielectron, 
          electrons.charge[ielectron], 
          true, false,
          electrons.ids[ielectron][m_electronVetoIDName],
          electrons.ids[ielectron][m_electronLooseIDName],
          electrons.ids[ielectron][m_electronMediumIDName],
          electrons.ids[ielectron][m_electronTightIDName],
          electrons.relativeIsoR03_withEA[ielectron]
      );
      
      for(const LepID::LepID& id: LepID::it){
        for(const LepIso::LepIso& iso: LepIso::it){ // loop over Iso not really needed since not considered for electrons (for the moment)
          uint16_t idx = LepIDIso(id, iso);
          if( m_lepton.ID[id] && m_lepton.iso[iso] )
            electrons_IDIso[idx].push_back(ielectron);
        }
      }
      
      leptons.push_back(m_lepton);
    }
  }

  ///////////////////////////
  //       MUONS           //
  ///////////////////////////
  
  #ifdef _ZA_DEBUG_
    std::cout << "Muons" << std::endl;
  #endif

  const MuonsProducer& muons = producers.get<MuonsProducer>("muons");

  for(uint16_t imuon = 0; imuon < muons.p4.size(); imuon++){
    if(muons.p4[imuon].Pt() > m_muonPtCut && abs(muons.p4[imuon].Eta()) < m_muonEtaCut ){
      
      Lepton m_lepton(
          muons.p4[imuon], 
          imuon,
          muons.charge[imuon], 
          false, true,
          muons.isLoose[imuon], // isVeto => for muons, re-use isLoose
          muons.isLoose[imuon],
          muons.isMedium[imuon],
          muons.isTight[imuon],
          muons.relativeIsoR04_withEA[imuon],
          muons.relativeIsoR04_withEA[imuon] > m_muonLooseIsoCut,
          muons.relativeIsoR04_withEA[imuon] > m_muonTightIsoCut
      );

      for(const LepID::LepID& id: LepID::it){
        for(const LepIso::LepIso& iso: LepIso::it){
          uint16_t idx = LepIDIso(id, iso);
          if( m_lepton.ID[id] && m_lepton.iso[iso] )
            muons_IDIso[idx].push_back(imuon);
        }
      }

      leptons.push_back(m_lepton);
    }
  }


  ///////////////////////////
  //       DILEPTONS       //
  ///////////////////////////

  #ifdef _ZA_DEBUG_
    std::cout << "Dileptons" << std::endl;
  #endif

  for(uint16_t i1 = 0; i1 < leptons.size(); i1++){
    for(uint16_t i2 = i1 + 1; i2 < leptons.size(); i2++){
      const Lepton& l1 = leptons[i1];
      const Lepton& l2 = leptons[i2];

      DiLepton m_diLepton;

      m_diLepton.p4 = l1.p4 + l2.p4; 
      m_diLepton.idxs = std::make_pair(l1.idx, l2.idx); 
      m_diLepton.lidxs = std::make_pair(i1, i2); 
      m_diLepton.isElEl = l1.isEl && l2.isEl;
      m_diLepton.isElMu = l1.isEl && l2.isMu;
      m_diLepton.isMuEl = l1.isMu && l2.isEl;
      m_diLepton.isMuMu = l1.isMu && l2.isMu;
      m_diLepton.isOS = l1.charge != l2.charge;
      m_diLepton.isSF = m_diLepton.isElEl || m_diLepton.isMuMu;
 
      // Save the combination of IDs
      for(const LepID::LepID& id1: LepID::it){
        for(const LepID::LepID& id2: LepID::it){
          uint16_t idx = LepLepID(id1, id2);
          m_diLepton.ID[idx] = l1.ID[id1] && l2.ID[id2];
        }
      }

      // Save the combination of isolations
      for(const LepIso::LepIso& iso1: LepIso::it){
        for(const LepIso::LepIso& iso2: LepIso::it){
          uint16_t idx = LepLepIso(iso1, iso2);
          m_diLepton.iso[idx] = l1.iso[iso1] && l2.iso[iso2];
        }
      }
      
      m_diLepton.DR = VectorUtil::DeltaR(l1.p4, l2.p4);
      m_diLepton.DEta = ZAAnalysis::DeltaEta(l1.p4, l2.p4);
      m_diLepton.DPhi = VectorUtil::DeltaPhi(l1.p4, l2.p4);

      diLeptons.push_back(m_diLepton);
    }
  }

  // Save indices to DiLeptons for the combinations of IDs & Isolationss
  for(uint16_t i = 0; i < diLeptons.size(); i++){
    const DiLepton& m_diLepton = diLeptons[i];
    
    for(const LepID::LepID& id1: LepID::it){
      for(const LepID::LepID& id2: LepID::it){
        for(const LepIso::LepIso& iso1: LepIso::it){
          for(const LepIso::LepIso& iso2: LepIso::it){
            
            uint16_t idx_ids = LepLepID(id1, id2);
            uint16_t idx_isos = LepLepIso(iso1, iso2);
            uint16_t idx_comb = LepLepIDIso(id1, iso1, id2, iso2);

            if(m_diLepton.ID[idx_ids] && m_diLepton.iso[idx_isos])
              diLeptons_IDIso[idx_comb].push_back(i);
          
          }
        }
      }
    }
    
  }

}


void ZAAnalyzer::registerCategories(CategoryManager& manager, const edm::ParameterSet& config) {
  manager.new_category<ZAAnalysis::ElElCategory>("elel", "Category with leading leptons as two electrons", config);
  // manager.new_category<TTAnalysis::ElMuCategory>("elmu", "Category with leading leptons as electron, muon", config);
  //manager.new_category<TTAnalysis::MuElCategory>("muel", "Category with leading leptons as muon, electron", config);
  manager.new_category<ZAAnalysis::MuMuCategory>("mumu", "Category with leading leptons as two muons", config);
}

#include <FWCore/PluginManager/interface/PluginFactory.h>
DEFINE_EDM_PLUGIN(ExTreeMakerAnalyzerFactory, ZAAnalyzer, "za_analyzer");
