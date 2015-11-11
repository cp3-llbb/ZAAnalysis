#pragma once

#include <cp3_llbb/ZAAnalysis/interface/ZATypes.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>

namespace ZAAnalysis {
  
  float DeltaEta(const myLorentzVector &v1, const myLorentzVector &v2);

  bool sortByBtag(const ZAAnalysis::Jet& _jet1, const ZAAnalysis::Jet& _jet2);
  bool sortByPt(const ZAAnalysis::Jet& _jet1, const ZAAnalysis::Jet& _jet2); 
}

