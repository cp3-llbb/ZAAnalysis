#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cp3_llbb/ZAAnalysis/interface/ZADileptonCategories.h>
#include <cp3_llbb/ZAAnalysis/interface/ZAAnalyzer.h>

#include <cp3_llbb/ZAAnalysis/interface/Types.h>
#include <cp3_llbb/ZAAnalysis/interface/Indices.h>

using namespace ZAAnalysis;

// ***** ***** *****
// Dilepton El-El category
// ***** ***** *****
bool ElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
  const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
  return electrons.p4.size() >= 2;
}

bool ElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const ZAAnalyzer& za = analyzers.get<ZAAnalyzer>("za");

  // If at least one DiLepton of highest Pt and of type ElEl among all ID pairs is found, keep event in this category

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
          if(za.diLeptons_IDIso[comb].size() >= 1) {
            if( za.diLeptons[ za.diLeptons_IDIso[comb][0] ].isElEl )
              return true;
          }
        
        }
      }
    }
  }

  return false;
}

void ElElCategory::register_cuts(CutManager& manager) {
  
  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          
          manager.new_cut(baseStrCategory + postFix, baseStrCategory + postFix);
          manager.new_cut(baseStrExtraDiLeptonVeto + postFix, baseStrExtraDiLeptonVeto + postFix);
          manager.new_cut(baseStrDiLeptonTriggerMatch + postFix, baseStrDiLeptonTriggerMatch + postFix);
          manager.new_cut(baseStrMllCut + postFix, baseStrMllCut + postFix);
          manager.new_cut(baseStrDiLeptonIsOS + postFix, baseStrDiLeptonIsOS + postFix);

        }
      }
    }
  }

}

void ElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const ZAAnalyzer& za = analyzers.get<ZAAnalyzer>("za");
  const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
    
          if(za.diLeptons_IDIso[comb].size() >= 1) {
            const DiLepton& m_diLepton = za.diLeptons[ za.diLeptons_IDIso[comb][0] ];
            
            if(m_diLepton.isElEl) {
              manager.pass_cut(baseStrCategory + postFix);

              if(m_diLepton.hlt_idxs.first >= 0 && m_diLepton.hlt_idxs.second >= 0){
                // We have fired a trigger. Now, check that it is actually a DoubleEG trigger
                if( checkHLT(hlt, m_diLepton.hlt_idxs.first, m_diLepton.hlt_idxs.second, HLT::DoubleEG) )
                  manager.pass_cut(baseStrDiLeptonTriggerMatch + postFix);
              }
              
              if(m_diLepton.p4.M() > m_MllCutSF)
                manager.pass_cut(baseStrMllCut + postFix);
              
              if(m_diLepton.isOS)
              {
                manager.pass_cut(baseStrDiLeptonIsOS + postFix);
                if (m_diLepton.p4.M() > m_lowLooseZcut  && m_diLepton.p4.M() < m_highLooseZcut )
                {
                  manager.pass_cut(baseStrLooseZCandidate + postFix);
                  if (m_diLepton.p4.M() > m_lowTightZcut  && m_diLepton.p4.M() < m_highTightZcut )
                    manager.pass_cut(baseStrTightZCandidate + postFix);
                }
              }
            }
          }
          
          // For electrons, in principe only veto using VetoID.
          // But since the user can access any cut he wants, he can take the IDVV_IsoWhatever cut.
          if(za.diLeptons_IDIso[comb].size() >= 2) { 
            manager.pass_cut(baseStrExtraDiLeptonVeto + postFix);
          }

        }
      }
    }
  }

}


// ***** ***** *****
// Dilepton Mu-Mu category
// ***** ***** *****
bool MuMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
  const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
  return muons.p4.size() >= 2;
}

bool MuMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const ZAAnalyzer& za = analyzers.get<ZAAnalyzer>("za");

  // It at least one DiLepton of highest Pt and of type MuMu among all ID pairs is found, keep event in this category

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
          if(za.diLeptons_IDIso[comb].size() >= 1) {
            if( za.diLeptons[ za.diLeptons_IDIso[comb][0] ].isMuMu )
              return true;
          }
        
        }
      }
    }
  }

  return false;
}

void MuMuCategory::register_cuts(CutManager& manager) {
  
  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          
          manager.new_cut(baseStrCategory + postFix, baseStrCategory + postFix);
          manager.new_cut(baseStrExtraDiLeptonVeto + postFix, baseStrExtraDiLeptonVeto + postFix);
          manager.new_cut(baseStrDiLeptonTriggerMatch + postFix, baseStrDiLeptonTriggerMatch + postFix);
          manager.new_cut(baseStrMllCut + postFix, baseStrMllCut + postFix);
          manager.new_cut(baseStrDiLeptonIsOS + postFix, baseStrDiLeptonIsOS + postFix);

        }
      }
    }
  }

}

void MuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const ZAAnalyzer& za = analyzers.get<ZAAnalyzer>("za");
  const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
    
          if(za.diLeptons_IDIso[comb].size() >= 1) {
            const DiLepton& m_diLepton = za.diLeptons[ za.diLeptons_IDIso[comb][0] ];
            
            if(m_diLepton.isMuMu) {
              manager.pass_cut(baseStrCategory + postFix);

              if(m_diLepton.hlt_idxs.first >= 0 && m_diLepton.hlt_idxs.second >= 0){
                // We have fired a trigger. Now, check that it is actually a DoubleMuon trigger
                if( checkHLT(hlt, m_diLepton.hlt_idxs.first, m_diLepton.hlt_idxs.second, HLT::DoubleMuon) )
                  manager.pass_cut(baseStrDiLeptonTriggerMatch + postFix);
              }
              
              if(m_diLepton.p4.M() > m_MllCutSF)
                manager.pass_cut(baseStrMllCut + postFix);
              
              if(m_diLepton.isOS)
              {
                manager.pass_cut(baseStrDiLeptonIsOS + postFix);
                if (m_diLepton.p4.M() > m_lowLooseZcut  && m_diLepton.p4.M() < m_highLooseZcut )
                {
                  manager.pass_cut(baseStrLooseZCandidate + postFix);
                  if (m_diLepton.p4.M() > m_lowTightZcut  && m_diLepton.p4.M() < m_highTightZcut )
                    manager.pass_cut(baseStrTightZCandidate + postFix);
                }
              } 
            }
          }
          
          if(za.diLeptons_IDIso[comb].size() >= 2) { 
            manager.pass_cut(baseStrExtraDiLeptonVeto + postFix);
          }

        }
      }
    }
  }

}
