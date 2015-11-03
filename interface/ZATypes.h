#pragma once

#include <utility>
#include <vector>
#include <limits>

#include <Math/PtEtaPhiE4D.h>
#include <Math/LorentzVector.h>
#include <Math/VectorUtil.h>

//#include <cp3_llbb/ZAAnalysis/interface/Indices.h>

// Needed because of gcc bug when using typedef and std::map
#define myLorentzVector ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>

namespace ZAAnalysis {

  // Forward declaration to use this here
  float DeltaEta(const myLorentzVector &v1, const myLorentzVector &v2);

  struct BaseObject {
    BaseObject(myLorentzVector p4): p4(p4) {}
    BaseObject() {}
    virtual ~BaseObject() {}

    myLorentzVector p4;
  };

  struct Lepton: BaseObject {
    Lepton()
    {}
    Lepton(myLorentzVector p4, uint16_t idx, uint16_t charge, bool isEl, bool isMu, bool isVeto = false, bool isLoose = false, bool isMedium = false, bool isTight = false, float isoValue = 0, bool isoLoose = false, bool isoTight = false):
      BaseObject(p4), 
      idx(idx), 
      charge(charge),
      isoValue(isoValue),
      isEl(isEl), 
      isMu(isMu),
      isIDVeto(false),
      isIDLoose(false),
      isIDMedium(false),
      isIDTight(false),
      isIsoLoose(false),
      isIsoTight(false) 
      {
        if(isEl)
          isIDVeto = isVeto;
        else
          isIDVeto = isLoose; // for muons, re-use Loose as Veto ID
        isIDLoose = isLoose;
        isIDMedium = isMedium;
        isIDTight = isTight;

        if(isMu){
          isIsoLoose = isoLoose;
          isIsoTight = isoTight;
        }else{
          isIsoLoose = true;
        }
      }
    
    uint16_t idx; // stores index to electron/muon arrays
    uint16_t charge;
    float isoValue;
    int16_t hlt_idx = -1; // Index to the matched HLT object. -1 if no match
    bool isEl;
    bool isMu;
    bool isIDVeto, isIDLoose, isIDMedium, isIDTight; // lepton ID: veto-loose-medium-tight
    bool isIsoLoose, isIsoTight; // lepton Iso: loose-tight (only for muons -> electrons only have loose)

    bool hlt_already_matched = false; // Internal flag; if true, it means this lepton has already been matched to an online object, even if no match has been found.

    int8_t pdg_id() const {
        int8_t id = (isEl) ? 11 : 13;
        return charge * id;
    }
  };
  
  struct DiLepton: BaseObject {
    DiLepton()
    {}
    DiLepton(const Lepton& _l1, const Lepton& _l2):
      BaseObject(_l1.p4 + _l2.p4),
      DR( ROOT::Math::VectorUtil::DeltaR(_l1.p4, _l2.p4) ),
      DEta( DeltaEta(_l1.p4, _l2.p4) ),
      DPhi( ROOT::Math::VectorUtil::DeltaPhi(_l1.p4, _l2.p4)),
      isLL( _l1.isIDLoose &&  _l2.isIDLoose ),
      isMM( _l1.isIDMedium &&  _l2.isIDMedium ),
      isTT( _l1.isIDTight &&  _l2.isIDTight ),
      isIsoLL( _l1.isIsoLoose && _l2.isIsoLoose),
      isIsoTT( _l1.isIsoTight && _l2.isIsoTight)
      {}
    bool isElEl, isElMu, isMuEl, isMuMu;
    bool isOS; // opposite sign
    bool isSF; // same flavour
    float DR;
    float DEta;
    float DPhi;
    bool isLL, isMM, isTT;
    bool isIsoLL, isIsoTT;

  };
 
  struct Jet: BaseObject {
    Jet()
      {}
    uint16_t idx; // index to jet array
    bool isIDLoose;
    bool isIDTight;
    bool isTLV;
    float minDRjl;  
    float CSVv2;
    bool isBWPL;
    bool isBWPM;
    bool isBWPT;
  };
  
  struct DiJet: BaseObject {
    DiJet()
      {}
    DiJet(const Jet& _jet1, const Jet& _jet2):
      BaseObject(_jet1.p4 + _jet2.p4),
      DR( ROOT::Math::VectorUtil::DeltaR(_jet1.p4, _jet2.p4) ),
      DEta( DeltaEta(_jet1.p4, _jet2.p4) ),
      DPhi( ROOT::Math::VectorUtil::DeltaPhi(_jet1.p4, _jet2.p4)),
      isML( _jet1.isBWPM && _jet2.isBWPL),
      isMM( _jet1.isBWPM && _jet2.isBWPM),
      isTM( _jet1.isBWPT && _jet2.isBWPM)
      {}
    float DR;
    float DEta;
    float DPhi;
    bool isML; // b-tagging working points
    bool isMM; // b-tagging working points
    bool isTM; // b-tagging working points
  };

  struct DiLepDiJet: BaseObject {
    DiLepDiJet()
      {}
    DiLepDiJet(const DiLepton& _diLepton, const DiJet& _diJet):
      BaseObject(_diLepton.p4 + _diJet.p4),
      diLepton(&_diLepton),
      diJet(&_diJet),
      DR_ll_jj( ROOT::Math::VectorUtil::DeltaR(_diLepton.p4, _diJet.p4) ),
      DEta_ll_jj( DeltaEta(_diLepton.p4, _diJet.p4) ),
      DPhi_ll_jj( ROOT::Math::VectorUtil::DeltaPhi(_diLepton.p4, _diJet.p4) )
      {}

    const DiLepton* diLepton;
    uint16_t diLepIdx;
    const DiJet* diJet;
    uint16_t diJetIdx;

    float DR_ll_jj, DEta_ll_jj, DPhi_ll_jj;
    
    float minDRjl, maxDRjl;
    float minDEtajl, maxDEtajl;
    float minDPhijl, maxDPhijl;
  };


}

