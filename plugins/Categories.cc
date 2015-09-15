#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/Category.h>
#include <cp3_llbb/Framework/interface/DileptonCategories.h>
#include <cp3_llbb/Framework/interface/DileptonAnalyzer.h>
#include <cp3_llbb/ZAAnalysis/interface/ZAAnalyzer.h>
#include <cp3_llbb/ZAAnalysis/interface/Categories.h>

/*
class MuMuJetJetCategory: public Category {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override {
        const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
        const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
        const JetsProducer& jets = producers.get<JetsProducer>("jets");
        if( muons.p4.size() >= 2 && jets.p4.size() >=2)
        {
            if( electrons.p4.size() >= 1 ) // if there is electrons at all, check the muons are the leading leptons
            {
                if( muons.p4[1].Pt() > electrons.p4[0].Pt() )
                    return true;
            } else {
                return true;
            }
        }
        return false;
    }
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override { return true; }
};
*/

void ZAElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {

        ElElCategory::evaluate_cuts_post_analyzers(manager,producers,analyzers);

        const ZAAnalyzer& za_analyzer = analyzers.get<ZAAnalyzer>("za");
        if (za_analyzer.selectedjets.size() >= 2 ) {manager.pass_cut("two_jets");}
        if (za_analyzer.selectedbjets.size() >= 2 ) {manager.pass_cut("two_bjets");}
        if (za_analyzer.isolatedElectrons.size() >= 2) {manager.pass_cut("two_isolated_el");}
        if (za_analyzer.isolatedMuons.size() >= 2) {manager.pass_cut("two_isolated_mu");}
};

void ZAMuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {

        MuMuCategory::evaluate_cuts_post_analyzers(manager,producers,analyzers);

        const ZAAnalyzer& za_analyzer = analyzers.get<ZAAnalyzer>("za");
        if (za_analyzer.selectedjets.size() >= 2 ) {manager.pass_cut("two_jets");}
        if (za_analyzer.selectedbjets.size() >= 2 ) {manager.pass_cut("two_bjets");}
        if (za_analyzer.isolatedElectrons.size() >= 2) {manager.pass_cut("two_isolated_el");}
        if (za_analyzer.isolatedMuons.size() >= 2) {manager.pass_cut("two_isolated_mu");}
};



