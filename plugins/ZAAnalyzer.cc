#include <cp3_llbb/ZAAnalysis/interface/ZAAnalyzer.h>


void ZAAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup, const ProducersManager& producers, const CategoryManager&) {

/*
    BRANCH(selectedjets,std::vector<LorentzVector>);
    BRANCH(selectedbjets,std::vector<LorentzVector>);
    BRANCH(dijets,LorentzVector);
    BRANCH(dibjets,LorentzVector);
    BRANCH(njets, int);
    BRANCH(nbjets, int);
*/

    std::cout << "analyse " << std::endl;
    const JetsProducer& jets = producers.get<JetsProducer>("jets");


    for( unsigned int ijet = 0 ; ijet < jets.p4.size() ; ijet++ ){

        if (jets.p4[ijet].pt() > 30 && abs(jets.p4[ijet].eta()) < 2.4) {
            selectedjets.push_back(jets.p4[ijet]);
            if (jets.getBTagDiscriminant(ijet,"combinedSecondaryVertexBTags") > 0.679 ) {
                selectedbjets.push_back(jets.p4[ijet]);                
            }
        } 
    }

    njets = selectedjets.size();
    nbjets = selectedbjets.size();

    // temporary definition only, one should keep each combination of b, and not only the highest pt one, I guess.

    if (selectedjets.size() > 2) {
        dijets = selectedjets[0]+selectedjets[1];
        if (selectedbjets.size() > 2) {
            dibjets = selectedbjets[0]+selectedbjets[1];
        }
    }



/*
    const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
    const MuonsProducer& muons = producers.get<MuonsProducer>("muons");

    BRANCH(dileptons_mue,std::vector<LorentzVector>);
    BRANCH(dileptons_mue_indices, std::vector<std::pair<int,int>>);

    for( unsigned int imuon = 0 ; imuon < muons.p4.size() ; imuon++ ){
        muons.EA_R03[imuon];
    }



    //Dilepton muon-electron

    std::cout << " electron size : " << electrons.p4.size() << std::endl;
    std::cout << " muon size     : " << muons.p4.size() << std::endl;

    for( unsigned int imuon = 0 ; imuon < muons.p4.size() ; imuon++ ){
        std::cout << "muon is loose? : " << std::endl;
        if(muons.isLoose[imuon]){
            std::cout << "   yes " << std::endl;
            for( unsigned int jmuon = 0 ; jmuon < muons.p4.size() ; jmuon++ ){

                if( muons.p4[jmuon].pt() < muons.p4[imuon].pt() &&
                    muons.isLoose[jmuon] &&
                    muons.charge[imuon] * muons.charge[jmuon] < 0 ){

                    LorentzVector dilepton = muons.p4[imuon] + muons.p4[jmuon];

                    if( dilepton.mass() > 20. )
                    {

                        std::cout << " di lepton candidate " << std::endl;

                        dileptons_mue.push_back(dilepton);
                        dileptons_mue_indices.push_back(std::make_pair(imuon, jmuon));

                    }

                }

            }

        }

    }
*/
}

#include <FWCore/PluginManager/interface/PluginFactory.h>
DEFINE_EDM_PLUGIN(ExTreeMakerAnalyzerFactory, ZAAnalyzer, "za_analyzer");
