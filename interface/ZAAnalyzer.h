#ifndef ZAANALYZER_H
#define ZAANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>
#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>

#include <cp3_llbb/ZAAnalysis/interface/Categories.h>
#include <cp3_llbb/Framework/interface/DileptonAnalyzer.h>

class ZAAnalyzer: public Framework::Analyzer {
    public:
        ZAAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config) {

        }

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const CategoryManager&) override;

        virtual void registerCategories(CategoryManager& manager) {
            //manager.new_category<MuMuCategory>("mumujj", "Category with two tight muons and two jets");
            manager.new_category<ZAElElCategory>("elel", "Category with two electrons");
        }

    private:
};

#endif
