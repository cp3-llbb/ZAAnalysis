#pragma once

#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/Category.h>



class MuMuCategory: public Category {
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

class ElElCategory: public Category {
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override {
        const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
        const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
        const JetsProducer& jets = producers.get<JetsProducer>("jets");
        if( electrons.p4.size() >= 2 )
        {
            if( muons.p4.size() >= 1 ) // if there is muons at all, check the electrons are the leading leptons
            {
                if( electrons.p4[1].Pt() > muons.p4[0].Pt() )
                    return true;
            } else {
                return true;
            }
        }
        return false;
    }
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override { return true; }
};


