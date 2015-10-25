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


  diJets_DRCut.resize( LepID::Count * LepIso::Count );
  diBJets_DRCut_BWP_PtOrdered.resize( LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  diBJets_DRCut_BWP_CSVv2Ordered.resize( LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  
  diLepDiJets_DRCut.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count );
  diLepDiBJets_DRCut_BWP_PtOrdered.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  diLepDiBJets_DRCut_BWP_CSVv2Ordered.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  
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
          muons.relativeIsoR04_withEA[imuon] < m_muonLooseIsoCut,
          muons.relativeIsoR04_withEA[imuon] < m_muonTightIsoCut
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

  // Save indices to DiLeptons for the combinations of IDs & Isolations
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

  ///////////////////////////
  //       JETS            //
  ///////////////////////////

  #ifdef _ZA_DEBUG_
    std::cout << "Jets" << std::endl;
  #endif

  const JetsProducer& jets = producers.get<JetsProducer>("jets");

  // First find the jets passing kinematic cuts and save them as Jet objects

  uint16_t jetCounter(0);
  for(uint16_t ijet = 0; ijet < jets.p4.size(); ijet++){
    // Save the jets that pass the kinematic cuts
    if( abs(jets.p4[ijet].Eta()) < m_jetEtaCut && jets.p4[ijet].Pt() > m_jetPtCut){
      Jet m_jet;
      
      m_jet.p4 = jets.p4[ijet];
      m_jet.idx = ijet;
      m_jet.ID[JetID::L] = jets.passLooseID[ijet]; 
      m_jet.ID[JetID::T] = jets.passTightID[ijet]; 
      m_jet.ID[JetID::TLV] = jets.passTightLeptonVetoID[ijet];
      m_jet.CSVv2 = jets.getBTagDiscriminant(ijet, m_jetCSVv2Name);
      m_jet.BWP[BWP::L] = m_jet.CSVv2 > m_jetCSVv2L;
      m_jet.BWP[BWP::M] = m_jet.CSVv2 > m_jetCSVv2M;
      m_jet.BWP[BWP::T] = m_jet.CSVv2 > m_jetCSVv2T;
      
      // Save minimal DR(l,j) using selected leptons, for each Lepton ID/Iso
      for(const LepID::LepID& id: LepID::it){
        for(const LepIso::LepIso& iso: LepIso::it){
              
          uint16_t idx_comb = LepIDIso(id, iso);
          
          for(const uint16_t& lepIdx: leptons_IDIso[idx_comb]){
            const Lepton& m_lepton = leptons[lepIdx];
            float DR = (float) VectorUtil::DeltaR(jets.p4[ijet], m_lepton.p4);
            if( DR < m_jet.minDRjl_lepIDIso[idx_comb] )
              m_jet.minDRjl_lepIDIso[idx_comb] = DR;
          }
          
          // Save the indices to Jets passing the selected jetID and minDRjl > cut for this lepton ID/Iso
          if( m_jet.minDRjl_lepIDIso[idx_comb] > m_jetDRleptonCut && jetIDAccessor(jets, ijet, m_jetID) ){
            selJets_selID_DRCut[idx_comb].push_back(jetCounter);

            // Out of these, save the indices for different b-tagging working points
            for(const BWP::BWP& wp: BWP::it){
              uint16_t idx_comb_b = LepIDIsoJetBWP(id, iso, wp);
              if(m_jet.BWP[wp])
                selBJets_DRCut_BWP_PtOrdered[idx_comb_b].push_back(jetCounter);
            }
          }            
        
        }
      }
      
      if(jetIDAccessor(jets, ijet, m_jetID)) // Save the indices to Jets passing the selected jet ID
        selJets_selID.push_back(jetCounter);
      
      selJets.push_back(m_jet);
      
      jetCounter++;
    }
  }

  // Sort the b-jets according to decreasing CSVv2 value
  selBJets_DRCut_BWP_CSVv2Ordered = selBJets_DRCut_BWP_PtOrdered;
  for(const LepID::LepID& id: LepID::it){
    for(const LepIso::LepIso& iso: LepIso::it){
      for(const BWP::BWP& wp: BWP::it){ 
        uint16_t idx_comb_b = LepIDIsoJetBWP(id, iso, wp);
        std::sort(selBJets_DRCut_BWP_CSVv2Ordered[idx_comb_b].begin(), selBJets_DRCut_BWP_CSVv2Ordered[idx_comb_b].end(), jetBTagDiscriminantSorter(jets, m_jetCSVv2Name, selJets));
      }
    }
  }
        
  ///////////////////////////
  //       DIJETS          //
  ///////////////////////////

  #ifdef _ZA_DEBUG_
    std::cout << "Dijets" << std::endl;
  #endif

  // Next, construct DiJets out of selected jets with selected ID (not accounting for minDRjl here)

  uint16_t diJetCounter(0);

  for(uint16_t j1 = 0; j1 < selJets_selID.size(); j1++){
    for(uint16_t j2 = j1 + 1; j2 < selJets_selID.size(); j2++){
      const uint16_t jidx1 = selJets_selID[j1];
      const Jet& jet1 = selJets[jidx1];
      const uint16_t jidx2 = selJets_selID[j2];
      const Jet& jet2 = selJets[jidx2];

      DiJet m_diJet; 
      m_diJet.p4 = jet1.p4 + jet2.p4;
      m_diJet.idxs = std::make_pair(jet1.idx, jet2.idx);
      m_diJet.jidxs = std::make_pair(jidx1, jidx2);
      
      m_diJet.DR = VectorUtil::DeltaR(jet1.p4, jet2.p4);
      m_diJet.DEta = DeltaEta(jet1.p4, jet2.p4);
      m_diJet.DPhi = VectorUtil::DeltaPhi(jet1.p4, jet2.p4);
     
      for(const BWP::BWP& wp1: BWP::it){
        for(const BWP::BWP& wp2: BWP::it){
          uint16_t comb = JetJetBWP(wp1, wp2);
          m_diJet.BWP[comb] = jet1.BWP[wp1] && jet2.BWP[wp2];
        }
      }
      
      for(const LepID::LepID& id: LepID::it){
        for(const LepIso::LepIso& iso: LepIso::it){
          uint16_t combIDIso = LepIDIso(id, iso);

          m_diJet.minDRjl_lepIDIso[combIDIso] = std::min(jet1.minDRjl_lepIDIso[combIDIso], jet2.minDRjl_lepIDIso[combIDIso]);
          
          // Save the DiJets which have minDRjl>cut, for each leptonIDIso
          if(m_diJet.minDRjl_lepIDIso[combIDIso] > m_jetDRleptonCut){
            diJets_DRCut[combIDIso].push_back(diJetCounter);

            // Out of these, save di-b-jets for each combination of b-tagging working points
            for(const BWP::BWP& wp1: BWP::it){
              for(const BWP::BWP& wp2: BWP::it){
                uint16_t combB = JetJetBWP(wp1, wp2);
                uint16_t combAll = LepIDIsoJetJetBWP(id, iso, wp1, wp2);
                if(m_diJet.BWP[combB])
                  diBJets_DRCut_BWP_PtOrdered[combAll].push_back(diJetCounter);
              }
            }
          
          }
        
        }
      }
      
      diJets.push_back(m_diJet); 
      diJetCounter++;
    }
  }

  // Order selected di-b-jets according to decreasing CSVv2 discriminant
  diBJets_DRCut_BWP_CSVv2Ordered = diBJets_DRCut_BWP_PtOrdered;
  for(const LepID::LepID& id: LepID::it){
    for(const LepIso::LepIso& iso: LepIso::it){
      for(const BWP::BWP& wp1: BWP::it){ 
        for(const BWP::BWP& wp2: BWP::it){ 
          uint16_t idx_comb_b = LepIDIsoJetJetBWP(id, iso, wp1, wp2);
          std::sort(diBJets_DRCut_BWP_CSVv2Ordered[idx_comb_b].begin(), diBJets_DRCut_BWP_CSVv2Ordered[idx_comb_b].end(), diJetBTagDiscriminantSorter(jets, m_jetCSVv2Name, diJets));
        }
      }
    }
  }


  ///////////////////////////
  //    EVENT VARIABLES    //
  ///////////////////////////
  
  #ifdef _ZA_DEBUG_
    std::cout << "Dileptons-dijets" << std::endl;
  #endif

  // leptons-(b-)jets

  uint16_t diLepDiJetCounter(0);

  for(uint16_t dilep = 0; dilep < diLeptons.size(); dilep++){
    const DiLepton& m_diLepton = diLeptons[dilep];
    
    for(uint16_t dijet = 0; dijet < diJets.size(); dijet++){
      const DiJet& m_diJet =  diJets[dijet];
      
      DiLepDiJet m_diLepDiJet(m_diLepton, dilep, m_diJet, dijet);

      m_diLepDiJet.minDRjl = std::min( {
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );
      m_diLepDiJet.maxDRjl = std::max( {
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );
      m_diLepDiJet.minDEtajl = std::min( {
          DeltaEta(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          DeltaEta(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          DeltaEta(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          DeltaEta(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );
      m_diLepDiJet.maxDEtajl = std::max( {
          DeltaEta(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          DeltaEta(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          DeltaEta(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          DeltaEta(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );
      m_diLepDiJet.minDPhijl = std::min( {
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );
      m_diLepDiJet.maxDPhijl = std::max( {
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );

      diLepDiJets.push_back(m_diLepDiJet);

      for(const LepID::LepID& id1: LepID::it){
        for(const LepID::LepID& id2: LepID::it){
          for(const LepIso::LepIso& iso1: LepIso::it){
            for(const LepIso::LepIso& iso2: LepIso::it){
              
              uint16_t combID = LepLepID(id1, id2);
              LepID::LepID minID = std::min(id1, id2);
              
              uint16_t combIso = LepLepIso(iso1, iso2);
              LepIso::LepIso minIso = std::min(iso1, iso2);

              uint16_t minCombIDIso = LepIDIso(minID, minIso);
              uint16_t diLepCombIDIso = LepLepIDIso(id1, iso1, id2, iso2);
             
              // Store objects for each combined lepton ID/Iso, with jets having minDRjl>cut for leptons corresponding to the loosest combination of the aforementioned ID/Iso
              if(m_diLepton.ID[combID] && m_diLepton.iso[combIso] && m_diJet.minDRjl_lepIDIso[minCombIDIso] > m_jetDRleptonCut){
                diLepDiJets_DRCut[diLepCombIDIso].push_back(diLepDiJetCounter);
                
                // Out of these, store combinations of b-tagging working points
                for(const BWP::BWP& wp1: BWP::it){
                  for(const BWP::BWP& wp2: BWP::it){
                    uint16_t combB = JetJetBWP(wp1, wp2);
                    uint16_t combAll = LepLepIDIsoJetJetBWP(id1, iso1, id2, iso2, wp1, wp2);
                    if(m_diJet.BWP[combB])
                      diLepDiBJets_DRCut_BWP_PtOrdered[combAll].push_back(diLepDiJetCounter);
                  }
                } // end b-jet loops

              } // end minDRjl>cut

            }
          } // end lepton iso loops
        }
      } // end lepton ID loops

      diLepDiJetCounter++;
    } // end dijet loop
  } // end dilepton loop

  // Order selected di-lepton-di-b-jets according to decreasing CSVv2 discriminant
  diLepDiBJets_DRCut_BWP_CSVv2Ordered = diLepDiBJets_DRCut_BWP_PtOrdered;
  
  for(const LepID::LepID& id1: LepID::it){
    for(const LepID::LepID& id2: LepID::it){
      
      for(const LepIso::LepIso& iso1: LepIso::it){
        for(const LepIso::LepIso& iso2: LepIso::it){
          
          for(const BWP::BWP& wp1: BWP::it){ 
            for(const BWP::BWP& wp2: BWP::it){ 
              
              uint16_t idx_comb_all = LepLepIDIsoJetJetBWP(id1, iso1, id2, iso2, wp1, wp2);
              std::sort(diLepDiBJets_DRCut_BWP_CSVv2Ordered[idx_comb_all].begin(), diLepDiBJets_DRCut_BWP_CSVv2Ordered[idx_comb_all].end(), diJetBTagDiscriminantSorter(jets, m_jetCSVv2Name, diLepDiJets));
            
            }
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
