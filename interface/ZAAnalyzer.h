#ifndef TTANALYZER_H
#define TTANALYZER_H

#include <cp3_llbb/Framework/interface/Analyzer.h>
#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>

#include <cp3_llbb/ZAAnalysis/interface/Categories.h>

class ZAAnalyzer: public Framework::Analyzer {
    public:
        ZAAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config) {

        }

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&) override;

        virtual void registerCategories(CategoryManager& manager) {
            manager.new_category<MuMuCategory>("mumu", "Category with leading leptons as two muons");
            manager.new_category<ElElCategory>("elel", "Category with leading leptons as two electrons");
        }

    private:
};

#endif
