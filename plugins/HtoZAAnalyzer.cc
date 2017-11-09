#include <cp3_llbb/ZAAnalysis/interface/HtoZAAnalyzer.h>
#include <cp3_llbb/Framework/interface/BTagsAnalyzer.h>
#include <cp3_llbb/ZAAnalysis/interface/Categories.h>
#include <cp3_llbb/ZAAnalysis/interface/GenStatusFlags.h>

#include <cp3_llbb/Framework/interface/EventProducer.h>
#include <cp3_llbb/Framework/interface/GenParticlesProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/LeptonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/METProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cmath>

#include <TH1F.h>
#include <TCanvas.h>

#define HtoZA_GEN_DEBUG (false)
#define TT_GEN_DEBUG (false)

void HtoZAAnalyzer::registerCategories(CategoryManager& manager, const edm::ParameterSet& config) {
    edm::ParameterSet newconfig = edm::ParameterSet(config);
    newconfig.addUntrackedParameter("m_analyzer_name", this->m_name);
    manager.new_category<MuMuCategory>("mumu", "Category with leading leptons as two muons", newconfig);
    manager.new_category<ElElCategory>("elel", "Category with leading leptons as two electrons", newconfig);
    manager.new_category<ElMuCategory>("elmu", "Category with leading leptons as electron, subleading as muon", newconfig);
    manager.new_category<MuElCategory>("muel", "Category with leading leptons as muon, subleading as electron", newconfig);
}



void HtoZAAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&, const ProducersManager& producers, const AnalyzersManager&, const CategoryManager&) {

    //Reset event
    leptons.clear();
    ll.clear();
    jj.clear();
    lljj.clear();
    met.clear();
    

    const JetsProducer& alljets = producers.get<JetsProducer>(m_jets_producer);
    const ElectronsProducer& allelectrons = producers.get<ElectronsProducer>(m_electrons_producer);
    const MuonsProducer& allmuons = producers.get<MuonsProducer>(m_muons_producer);
    const EventProducer& fwevent = producers.get<EventProducer>("event");
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    const METProducer& pf_met = producers.get<METProducer>(m_met_producer);


    if(!event.isRealData()){
    
    size_t n_taus = 0;
    const GenParticlesProducer& gp = producers.get<GenParticlesProducer>("gen_particles");

#if HtoZA_GEN_DEBUG
    std::function<void(size_t)> print_mother_chain = [&gp, &print_mother_chain](size_t p) {
    
        if (gp.pruned_mothers_index[p].empty()) {
            std::cout << std::endl;
            return;
        }
 
        // Get the index in gp.pruned_p4 of the particle that is a mother
        size_t index = gp.pruned_mothers_index[p][0];
            std::cout << " <- #" << index << "(" << gp.pruned_pdg_id[index] << ")" << std::endl;
            print_mother_chain(index);
    };
#endif


    std::function<bool(size_t, size_t)> pruned_decays_from = [&pruned_decays_from, &gp](size_t particle_index, size_t mother_index) -> bool {
      // Iterator over all pruned particles to find if the particle `particle_index` has `mother_index` in its decay history
      if (gp.pruned_mothers_index[particle_index].empty())
          return false;

      size_t index = gp.pruned_mothers_index[particle_index][0];

      if (index == mother_index) {
          return true;
      }

      if(pruned_decays_from(index, mother_index))
          return true;

      return false;
    };

    std::function<bool(size_t, size_t, bool)> pruned_decays_from_pdg_id = [&pruned_decays_from_pdg_id, &gp](size_t particle_index, uint64_t pdg_id, bool direct) -> bool {
      // Iterator over all pruned particles to find if the particle `particle_index` decays from a particle with pdg id == pdg_id
      if (gp.pruned_mothers_index[particle_index].empty())
          return false;

      size_t index = gp.pruned_mothers_index[particle_index][0];

      if (std::abs(gp.pruned_pdg_id[index]) == pdg_id) {
          return true;
      }

      if (!direct && pruned_decays_from_pdg_id(index, pdg_id, direct))
          return true;

      return false;
    };


    // Construct signal gen info
    
    gen_iH = -1;
    gen_iA = -1;
    gen_iA_afterFSR = -1;
    gen_iZ = -1;
    gen_iZ_afterFSR = -1;
    gen_ih = -1;
    gen_ih_afterFSR = -1;
    gen_iB = gen_iBbar = -1;
    gen_iB_afterFSR = gen_iBbar_afterFSR = -1;
    gen_iW = -1;
    gen_iW_afterFSR = -1;
    gen_iLminus = gen_iLplus = -1;
    gen_iLminus_afterFSR = gen_iLplus_afterFSR = -1;
    gen_iNu1 = gen_iNu2 = -1;

    for (unsigned int ip = 0; ip < gp.pruned_p4.size(); ip++) {
        std::bitset<15> flags (gp.pruned_status_flags[ip]);

        if (!flags.test(8))
            continue;

        int64_t pdg_id = gp.pruned_pdg_id[ip];
        //int8_t status = gp.pruned_status[ip];

#if HtoZA_GEN_DEBUG
        std::cout << "[" << ip << "] pdg id: " << pdg_id << " flags: " << flags << " p = " << gp.pruned_p4[ip] << " status: " << (int) status << std::endl;
        std::cout << "gp.pruned_pdg_id.size(): " << gp.pruned_pdg_id.size() << std::endl;
        std::cout << "gp.pruned_p4.size(): " << gp.pruned_p4.size() << std::endl;
        print_mother_chain(ip);
#endif

        auto p4 = gp.pruned_p4[ip];
        if (std::abs(pdg_id) == 35) {
            ASSIGN_HtoZA_GEN_INFO_NO_FSR(H, "H");
#if HtoZA_GEN_DEBUG
            std::cout << "ONE H FOUND WITH gen_iH: " << (int) gen_iH << std::endl;
#endif
        } else if (pdg_id == 36) {
            ASSIGN_HtoZA_GEN_INFO(A, "A");
#if HtoZA_GEN_DEBUG
            std::cout << "ONE A FOUND WITH gen_iA: " << (int) gen_iA << std::endl;
#endif
        } else if (pdg_id == 23) {
            ASSIGN_HtoZA_GEN_INFO(Z, "Z"); 
#if HtoZA_GEN_DEBUG
            std::cout << "ONE Z FOUND WITH gen_iZ: " << (int) gen_iZ << std::endl;
#endif
        } else if (pdg_id == 25) {
            ASSIGN_HtoZA_GEN_INFO(h, "SM Higgs");
        }

        // Only look for A and Z decays (not for the H decay because it's often off-shell) -- check it
        if ((gen_iA == -1) || (gen_iZ == -1))
            continue;

        // And if the particles come directly from Z and A decays (and H???)
        bool from_Z_decay = pruned_decays_from(ip, gen_iZ);
        bool from_A_decay = pruned_decays_from(ip, gen_iA);

        // Only keep particles coming from Z and A decays
        if (!from_Z_decay && !from_A_decay)
            continue;

        if (pdg_id == 5) {
            ASSIGN_HtoZA_GEN_INFO(B, "B");
        } else if (pdg_id == -5) {
            ASSIGN_HtoZA_GEN_INFO(Bbar, "Bbar");        
        }

        if (pruned_decays_from_pdg_id(ip, 5, false))
            continue;

        // Count the number of tau coming directly from a W or a Z
        if ((std::abs(pdg_id) == 15) && (pruned_decays_from_pdg_id(ip, 24, true) || pruned_decays_from_pdg_id(ip, 23, true))) {
            n_taus++;
        }

        // Electron, muon, tauon
        if ((pdg_id == 11) || (pdg_id == 13) || (pdg_id == 15)) {
            ASSIGN_HtoZA_GEN_INFO(Lminus, "L-");
            //std::cout << "AFTER ASSIGNEMENT: gen_iLminus = " << (int) gen_iLminus << std::endl;
        } else if ((pdg_id == -11) || (pdg_id == -13) || (pdg_id == -15)) {
            ASSIGN_HtoZA_GEN_INFO(Lplus, "L+");
            //std::cout << "AFTER ASSIGNEMENT: gen_iLplus = " << (int) gen_iLplus << std::endl;
        } else if (std::abs(pdg_id) == 24) {
            ASSIGN_HtoZA_GEN_INFO(W, "W");
        } else if ((std::abs(pdg_id) == 12) || (std::abs(pdg_id) == 14) || (std::abs(pdg_id) == 16)) {
            ASSIGN_HtoZA_GEN_INFO_2_NO_FSR(Nu1, Nu2, "neutrinos");
          }
    } // End of loop

#if HtoZA_GEN_DEBUG
    if (gen_iH == -1 ) {
        offshellH_counter++;
    }
    std::cout << "offshellH_counter: " << offshellH_counter << std::endl;
    events_counter++;
    std::cout << "events_counter: " << events_counter << std::endl;
#endif

      // Get the invariant mass of Z+A
      if ((gen_iZ != -1) && (gen_iA != -1)) {
          gen_mZA = (gen_Z + gen_A).M();
          gen_costhetastar = getCosThetaStar_CS(gen_Z, gen_A);
      }

 #if HtoZA_GEN_DEBUG
      std::cout << "*************************" << std::endl;
      PRINT_PARTICULE(H);
      PRINT_RESONANCE(Z, A);
      PRINT_RESONANCE(B, Bbar);
      PRINT_RESONANCE(Lminus, Lplus);
      PRINT_RESONANCE_NO_FSR(Nu1, Nu2);
      std::cout << "*************************" << std::endl;

      //Rebuild resonances for consistency checks
      auto LplusLminus = gen_Lplus + gen_Lminus;
      std::cout << "    gen_(L+ L-).M() = " << LplusLminus.M() << std::endl;

      auto LplusLminus_afterFSR = gen_Lplus_afterFSR + gen_Lminus_afterFSR;
      std::cout << "    gen_(L+ L-)_afterFSR.M() = " << LplusLminus_afterFSR.M() << std::endl;

      auto BBbar = gen_B + gen_Bbar;
      std::cout << "    gen_(B Bbar).M() = " << BBbar.M() << std::endl;

      auto BBbar_afterFSR = gen_B_afterFSR + gen_Bbar_afterFSR;
      std::cout << "    gen_(B Bbar)_afterFSR.M() = " << BBbar_afterFSR.M() << std::endl;

      auto LLBB = gen_Lplus + gen_Lminus + gen_B + gen_Bbar;
      std::cout << "    gen_(L+ L- B Bbar).M() = " << LLBB.M() << std::endl;

      auto LLBB_afterFSR = gen_Lplus_afterFSR + gen_Lminus_afterFSR + gen_B_afterFSR + gen_Bbar_afterFSR;
      std::cout << "    gen_(L+ L- B Bbar)_afterFSR.M() = " << LLBB_afterFSR.M() << std::endl;
      std::cout << "+++++++++++++++++++++++++++" << std::endl;

#endif


        // ***** ***** *****
        // Matching
        // ***** ***** *****
        for (auto p4: alljets.gen_p4) {
            gen_deltaR_jet_B.push_back(deltaR(p4, gen_B));
            gen_deltaR_jet_Bbar.push_back(deltaR(p4, gen_Bbar));
            gen_deltaR_jet_B_afterFSR.push_back(deltaR(p4, gen_B_afterFSR));
            gen_deltaR_jet_Bbar_afterFSR.push_back(deltaR(p4, gen_Bbar_afterFSR));
        }
        for (auto p4: allelectrons.gen_p4) {
            gen_deltaR_electron_L1.push_back(deltaR(p4, gen_Lminus));
            gen_deltaR_electron_L2.push_back(deltaR(p4, gen_Lplus));
            gen_deltaR_electron_L1_afterFSR.push_back(deltaR(p4, gen_Lminus_afterFSR));
            gen_deltaR_electron_L2_afterFSR.push_back(deltaR(p4, gen_Lplus_afterFSR));
        }
        for (auto p4: allmuons.gen_p4) {
            gen_deltaR_muon_L1.push_back(deltaR(p4, gen_Lminus));
            gen_deltaR_muon_L2.push_back(deltaR(p4, gen_Lplus));
            gen_deltaR_muon_L1_afterFSR.push_back(deltaR(p4, gen_Lminus_afterFSR));
            gen_deltaR_muon_L2_afterFSR.push_back(deltaR(p4, gen_Lplus_afterFSR));
        }
    }


    LorentzVector null_p4(0., 0., 0., 0.);
    float event_weight = fwevent.weight;
    float tmp_count_has2leptons = 0.;
    float tmp_count_has2leptons_elel = 0.;
    float tmp_count_has2leptons_elmu = 0.;
    float tmp_count_has2leptons_muel = 0.;
    float tmp_count_has2leptons_mumu = 0.;
    float tmp_count_has2leptons_1lljj = 0.;
    float tmp_count_has2leptons_elel_1lljj = 0.;
    float tmp_count_has2leptons_elmu_1lljj = 0.;
    float tmp_count_has2leptons_muel_1lljj = 0.;
    float tmp_count_has2leptons_mumu_1lljj = 0.;
    float tmp_count_has2leptons_1lljj_2btagM = 0.;
    float tmp_count_has2leptons_elel_1lljj_2btagM = 0.;
    float tmp_count_has2leptons_elmu_1lljj_2btagM = 0.;
    float tmp_count_has2leptons_muel_1lljj_2btagM = 0.;
    float tmp_count_has2leptons_mumu_1lljj_2btagM = 0.;

    // ***** ***** *****
    // Trigger Matching
    // ***** ***** *****

    // the actual trigger matching to dilepton HLT paths happens only once we have a dilepton candidate to consider in the event
    
    // ********** 
    // Leptons and dileptons
    // ********** 



    // Define a function that is true whenever the electron passes the trigger
    static auto electron_pass_HLT_ID = [&allelectrons, this](size_t index) {
        auto electron = allelectrons.products[index];

        // Use POG HLT-safe id
        
        bool result = allelectrons.ids[index][m_electron_hlt_safe_wp_name];

        // Add dxy and dz cuts described at https://github.com/latinos/LatinoAnalysis/blob/de026c531cec33d6f1cb999ec720e52692596116/Gardener/python/variables/LeptonSel_cfg.py
        if (electron->isEB()) {
            result &= std::abs(allelectrons.dz[index]) < 0.1;
            result &= std::abs(allelectrons.dxy[index]) < 0.05;
            result &= std::abs(allelectrons.relativeIsoR04_withEA[index]) < 0.05880;
        } else {
            result &= std::abs(allelectrons.dz[index]) < 0.2;
            result &= std::abs(allelectrons.dxy[index]) < 0.1;
            result &= std::abs(allelectrons.relativeIsoR04_withEA[index]) < 0.0571;
        }

        return result;
    };


    // Fill lepton structures
    // allelectrons is a vector with the reco electrons but it also contains gen info (in allelectrons.gen_p4)
    for (unsigned int ielectron = 0; ielectron < allelectrons.p4.size(); ielectron++)
    {
        if (allelectrons.p4[ielectron].Pt() > m_subleadingElectronPtCut && fabs(allelectrons.p4[ielectron].Eta()) < m_electronEtaCut)
        {
            // some selection
            // Ask for loose ID
            if (!allelectrons.ids[ielectron][m_electron_mva_wp90_name])
                continue;

            HtoZA::Lepton ele;
            ele.p4 = allelectrons.p4[ielectron];
            ele.charge = allelectrons.charge[ielectron];
            ele.idx = ielectron;
            ele.isMu = false;
            ele.isEl = true;
            ele.ele_hlt_id = electron_pass_HLT_ID(ielectron);

            ele.gen_matched = allelectrons.matched[ielectron];
            ele.gen_p4 = ele.gen_matched ? allelectrons.gen_p4[ielectron] : null_p4;
            ele.gen_DR = ele.gen_matched ? ROOT::Math::VectorUtil::DeltaR(ele.p4, ele.gen_p4): -1.;
            ele.gen_DPtOverPt = ele.gen_matched ? (ele.p4.Pt() - ele.gen_p4.Pt()) / ele.p4.Pt() : -10.;
            ele.hlt_leg1 = false;
            ele.hlt_leg2 = false;

            ele.sc_eta = allelectrons.products[ielectron]->superCluster()->eta();

            leptons.push_back(ele);
        }
    } // end of loop over electrons

    for (unsigned int imuon = 0; imuon < allmuons.p4.size(); imuon++)
    {
        if (allmuons.p4[imuon].Pt() > m_subleadingMuonPtCut && fabs(allmuons.p4[imuon].Eta()) < m_muonEtaCut)
        {
            // Ask for tight ID & tight ISO
            if (!allmuons.isTight[imuon] || allmuons.relativeIsoR04_deltaBeta[imuon] >= m_muonTightIsoCut)
                continue;

            HtoZA::Lepton mu;
            mu.p4 = allmuons.p4[imuon];
            mu.charge = allmuons.charge[imuon];
            mu.idx = imuon;
            mu.isMu = true;
            mu.isEl = false;
            mu.gen_matched = allmuons.matched[imuon];
            mu.gen_p4 = mu.gen_matched ? allmuons.gen_p4[imuon] : null_p4;
            mu.gen_DR = mu.gen_matched ? ROOT::Math::VectorUtil::DeltaR(mu.p4, mu.gen_p4) : -1.;
            mu.gen_DPtOverPt = mu.gen_matched ? (mu.p4.Pt() - mu.gen_p4.Pt()) / mu.p4.Pt() : -10.;
            mu.hlt_leg1 = false;
            mu.hlt_leg2 = false;
        
            leptons.push_back(mu);
        }
    } // end of loop over muons



    // sort leptons by pt (ignoring flavour, id and iso)
    std::sort(leptons.begin(), leptons.end(), [](const HtoZA::Lepton& lep1, const HtoZA::Lepton& lep2) { return lep1.p4.Pt() > lep2.p4.Pt();});

    for (unsigned int ilep1 = 0; ilep1 < leptons.size(); ilep1++)
    {
        if ((leptons[ilep1].isMu && leptons[ilep1].p4.Pt() < m_leadingMuonPtCut) || (leptons[ilep1].isEl && leptons[ilep1].p4.Pt() < m_leadingElectronPtCut))
            continue;

        for (unsigned int ilep2 = ilep1+1; ilep2 < leptons.size(); ilep2++)
        {
            HtoZA::Dilepton dilep;
            dilep.p4 = leptons[ilep1].p4 + leptons[ilep2].p4;
            // indices for dilepton: simply the pair (ilep1, ilep2)
            dilep.idxs = std::make_pair(leptons[ilep1].idx, leptons[ilep2].idx);
            dilep.ilep1 = ilep1;
            dilep.ilep2 = ilep2;
            // if the two leptons have opposite sign
            dilep.isOS = leptons[ilep1].charge * leptons[ilep2].charge < 0;

            dilep.isPlusMinus = leptons[ilep1].charge > 0 && leptons[ilep2].charge < 0;
            dilep.isMinusPlus = leptons[ilep1].charge < 0 && leptons[ilep2].charge > 0;
            dilep.isMuMu = leptons[ilep1].isMu && leptons[ilep2].isMu;
            dilep.isElEl = leptons[ilep1].isEl && leptons[ilep2].isEl;
            dilep.isElMu = leptons[ilep1].isEl && leptons[ilep2].isMu;
            dilep.isMuEl = leptons[ilep1].isMu && leptons[ilep2].isEl;
            // if the two leptons have same flavour
            dilep.isSF = dilep.isMuMu || dilep.isElEl;

            dilep.DR_l_l = ROOT::Math::VectorUtil::DeltaR(leptons[ilep1].p4, leptons[ilep2].p4);
            dilep.DPhi_l_l = fabs(ROOT::Math::VectorUtil::DeltaPhi(leptons[ilep1].p4, leptons[ilep2].p4));
            dilep.ht_l_l = leptons[ilep1].p4.Pt() + leptons[ilep2].p4.Pt();
            // We are selecting all the possible combinations of dileptons in an event. But which one is the
            // one that triggered the event? You set this by calling the function "matchOfflineLepton()".
            if (!hlt.paths.empty()) {
                matchOfflineLepton(hlt, dilep);
                dilep.hlt_idxs = std::make_pair(leptons[dilep.ilep1].hlt_idx, leptons[dilep.ilep2].hlt_idx);
            }
            dilep.gen_matched = leptons[ilep1].gen_matched && leptons[ilep2].gen_matched;
            dilep.gen_p4 = dilep.gen_matched ? leptons[ilep1].gen_p4 + leptons[ilep2].gen_p4 : null_p4;
            dilep.gen_DR = dilep.gen_matched ? ROOT::Math::VectorUtil::DeltaR(dilep.p4, dilep.gen_p4) : -1.;
            dilep.gen_DPtOverPt = dilep.gen_matched ? (dilep.p4.Pt() - dilep.gen_p4.Pt()) / dilep.p4.Pt() : -10.;


            if (event.isRealData()) {
                dilep.trigger_efficiency = 1.;
                dilep.trigger_efficiency_downVariated = 1.;
                dilep.trigger_efficiency_upVariated = 1.;
            } else {
                fillTriggerEfficiencies(leptons[ilep1], leptons[ilep2], dilep);
            }
            
            // Some selection
            // Note that ID and isolation criteria are in both electron and muon loops
            // Require that the leptons have opposite sign and same flavor
            //if (!dilep.isOS || !dilep.isSF)
            if (!dilep.isOS)
			    continue;

            // FIXME L1 EMTF bug mitigation -- cut the overlap on data if it's a run affected by the bug
            // On MC, apply the fraction of lumi the bug was not present
            if (dilep.isMuMu && isCSCWithOverlap(leptons[ilep1], leptons[ilep2])) {
                if (event.isRealData() && fwevent.run < 278167) {
                    continue;
                } else if (!event.isRealData()) {
                    dilep.trigger_efficiency *= 0.5265;
                    dilep.trigger_efficiency_downVariated *= 0.5265;
                    dilep.trigger_efficiency_upVariated *= 0.5265;
                }
            }


            // Throw event if there is no matched dilepton trigger path (only on data)
            if (event.isRealData() && !((leptons[dilep.ilep1].hlt_leg1 && leptons[dilep.ilep2].hlt_leg2)
                                  || (leptons[dilep.ilep1].hlt_leg2 && leptons[dilep.ilep2].hlt_leg1))) {
                continue;
            }

            // Counters
            tmp_count_has2leptons = event_weight;
            if (dilep.isElEl)
                tmp_count_has2leptons_elel = event_weight;
            if (dilep.isElMu)
                tmp_count_has2leptons_elmu = event_weight;
            if (dilep.isMuEl)
                tmp_count_has2leptons_muel = event_weight;
            if (dilep.isMuMu)
                tmp_count_has2leptons_mumu = event_weight;

            // Fill
            ll.push_back(dilep);
        } //end loop with lep2
    } //end loop with lep1

    // have the ll collection sorted by ht:
    // "sort" outputs the sorted vector
    std::sort(ll.begin(), ll.end(), [&](HtoZA::Dilepton& a, HtoZA::Dilepton& b){return a.ht_l_l > b.ht_l_l;});

    // Keep only the first ll candidate
    if (ll.size() > 1) {
        ll.resize(1); // discard all the elements beyond the first one
    }


    // ***** 
    // Adding MET(s) 
    // *****
    HtoZA::Met mymet;
    mymet.p4 = pf_met.p4;
    mymet.gen_matched = false;
    mymet.gen_p4 = null_p4;
    mymet.gen_DR = -1.;
    mymet.gen_DPhi = -1.;
    mymet.gen_DPtOverPt = -10.;
    if (!event.isRealData())
    { // genMet is not constructed in the framework, so construct it manually out of the neutrinos hanging around the mc particles
        const GenParticlesProducer& gp = producers.get<GenParticlesProducer>("gen_particles");
        for (unsigned int ip = 0; ip < gp.pruned_p4.size(); ip++) {
            std::bitset<15> flags (gp.pruned_status_flags[ip]);
            if (!flags.test(13)) continue; // take the last copies
            if (abs(gp.pruned_pdg_id[ip]) == 12 || abs(gp.pruned_pdg_id[ip]) == 14 || abs(gp.pruned_pdg_id[ip]) == 16)
            {
                mymet.gen_matched = true;
                mymet.gen_p4 += gp.pruned_p4[ip];
            }
        }
        mymet.gen_DR = mymet.gen_matched ? ROOT::Math::VectorUtil::DeltaR(mymet.p4, mymet.gen_p4) : -1.;
        mymet.gen_DPhi = mymet.gen_matched ? fabs(ROOT::Math::VectorUtil::DeltaPhi(mymet.p4, mymet.gen_p4)) : -1.;
        mymet.gen_DPtOverPt = mymet.gen_matched ? (mymet.p4.Pt() - mymet.gen_p4.Pt()) / mymet.p4.Pt() : -10.;
    }
    met.push_back(mymet);


    // ***** 
    // Jets and dijets 
    // *****

    for (unsigned int ijet = 0; ijet < alljets.p4.size(); ijet++)
    {
        float correctionFactor = m_applyBJetRegression ? alljets.regPt[ijet] / alljets.p4[ijet].Pt() : 1.;

        if ((alljets.p4[ijet].Pt() * correctionFactor > m_jetPtCut) && (fabs(alljets.p4[ijet].Eta()) < m_jetEtaCut))
        {
        
            if (!alljets.passLooseID[ijet])
                continue;

            HtoZA::Jet myjet;
            myjet.p4 = alljets.p4[ijet] * correctionFactor;
            myjet.idx = ijet;

            myjet.CSV = alljets.getBTagDiscriminant(ijet, "pfCombinedInclusiveSecondaryVertexV2BJetTags");
            myjet.CMVAv2 = alljets.getBTagDiscriminant(ijet, m_jet_bDiscrName_cMVAv2);
            myjet.deepCSV = alljets.getBTagDiscriminant(ijet, m_jet_bDiscrName_deepCSV_probb) + alljets.getBTagDiscriminant(ijet, m_jet_bDiscrName_deepCSV_probbb);
            // Ask for medium WP discr. cut for both cMVAv2 and deepCSV
            myjet.btag_cMVAv2_M = myjet.CMVAv2 > m_jet_bDiscrCut_cMVAv2_medium;
            myjet.btag_deepCSV_M = myjet.deepCSV > m_jet_bDiscrCut_deepCSV_medium;
            myjet.gen_matched_bParton = (std::abs(alljets.partonFlavor[ijet]) == 5);
            myjet.gen_matched_bHadron = (alljets.hadronFlavor[ijet]) == 5;
            myjet.gen_matched = alljets.matched[ijet];
            myjet.gen_p4 = myjet.gen_matched ? alljets.gen_p4[ijet] : null_p4;
            myjet.gen_DR = myjet.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myjet.p4, myjet.gen_p4) : -1.;
            myjet.gen_DPtOverPt = myjet.gen_matched ? (myjet.p4.Pt() - myjet.gen_p4.Pt()) / myjet.p4.Pt() : -10.;
            myjet.gen_b = (alljets.hadronFlavor[ijet]) == 5; // redundant with gen_matched_bHadron defined above
            myjet.gen_c = (alljets.hadronFlavor[ijet]) == 4;
            //light jets:
            myjet.gen_l = (alljets.hadronFlavor[ijet]) < 4;

            // We don't want a lepton close to the jet
            bool isThereACloseSelectedLepton = false;
            for (auto& mylepton: leptons) {
                if (ROOT::Math::VectorUtil::DeltaR(myjet.p4, mylepton.p4) < m_minDR_l_j_Cut) {
                    isThereACloseSelectedLepton = true;
                    break;
                }
            }

            if (isThereACloseSelectedLepton)
                continue;


            jets.push_back(myjet);
        }
    } // end of loop over jets

    // Do NOT change the loop logic here: we expect [0] to be made out of the leading jets
    for (unsigned int ijet1 = 0; ijet1 < jets.size(); ijet1++)
    {
        for (unsigned int ijet2 = ijet1 + 1; ijet2 < jets.size(); ijet2++)
        {
            HtoZA::Dijet myjj;
            myjj.p4 = jets[ijet1].p4 + jets[ijet2].p4;
            myjj.idxs = std::make_pair(jets[ijet1].idx, jets[ijet2].idx);
            myjj.ijet1 = ijet1;
            myjj.ijet2 = ijet2;
            // Commented things
            myjj.btag_cMVAv2_MM = jets[ijet1].btag_cMVAv2_M && jets[ijet2].btag_cMVAv2_M;
            myjj.btag_deepCSV_MM = jets[ijet1].btag_deepCSV_M && jets[ijet2].btag_deepCSV_M;
            // Commented things
            myjj.sumCSV = jets[ijet1].CSV + jets[ijet2].CSV;
            myjj.sumCMVAv2 = jets[ijet1].CMVAv2 + jets[ijet2].CMVAv2;
            myjj.sumDeepCSV = jets[ijet1].deepCSV + jets[ijet2].deepCSV;
            myjj.DR_j_j = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, jets[ijet2].p4);
            myjj.DPhi_j_j = fabs(ROOT::Math::VectorUtil::DeltaPhi(jets[ijet1].p4, jets[ijet2].p4));
            myjj.ht_j_j = jets[ijet1].p4.Pt() + jets[ijet2].p4.Pt();
            myjj.gen_matched_bbPartons = jets[ijet1].gen_matched_bParton && jets[ijet2].gen_matched_bParton;
            myjj.gen_matched_bbHadrons = jets[ijet1].gen_matched_bHadron && jets[ijet2].gen_matched_bHadron;
            myjj.gen_matched = jets[ijet1].gen_matched && jets[ijet2].gen_matched;
            myjj.gen_p4 = myjj.gen_matched ? jets[ijet1].gen_p4 + jets[ijet2].gen_p4 : null_p4;
            myjj.gen_DR = myjj.gen_matched ? ROOT::Math::VectorUtil::DeltaR(myjj.p4, myjj.gen_p4) : -1.;
            myjj.gen_DPtOverPt = myjj.gen_matched ? (myjj.p4.Pt() - myjj.gen_p4.Pt()) / myjj.p4.Pt() : -10.;
            myjj.gen_bb = (jets[ijet1].gen_b && jets[ijet2].gen_b);
            myjj.gen_bc = (jets[ijet1].gen_b && jets[ijet2].gen_c) || (jets[ijet1].gen_c && jets[ijet2].gen_b);
            myjj.gen_bl = (jets[ijet1].gen_b && jets[ijet2].gen_l) || (jets[ijet1].gen_l && jets[ijet2].gen_b);
            myjj.gen_cc = (jets[ijet1].gen_c && jets[ijet2].gen_c);
            myjj.gen_cl = (jets[ijet1].gen_c && jets[ijet2].gen_l) || (jets[ijet1].gen_l && jets[ijet2].gen_c);
            myjj.gen_ll = (jets[ijet1].gen_l && jets[ijet2].gen_l);
            jj.push_back(myjj);
        } //end of loop wtih jet2
    } //end of loop with jet1

    // have the jj collection sorted by ht
    std::sort(jj.begin(), jj.end(), [&](HtoZA::Dijet& a, HtoZA::Dijet& b){return a.p4.Pt() > b.p4.Pt();});


    // ********** 
    // lljj, llbb
    // **********
    for (unsigned int ill = 0; ill < ll.size(); ill++)
    {
        for (unsigned int ijj = 0; ijj < jj.size(); ijj++)
        {
            unsigned int ijet1 = jj[ijj].ijet1;
            unsigned int ijet2 = jj[ijj].ijet2;
            unsigned int ilep1 = ll[ill].ilep1;
            unsigned int ilep2 = ll[ill].ilep2;
            HtoZA::DileptonDijet mylljj;
            mylljj.p4 = ll[ill].p4 + jj[ijj].p4;
            mylljj.lep1_p4 = leptons[ilep1].p4;
            mylljj.lep2_p4 = leptons[ilep2].p4;
            mylljj.jet1_p4 = jets[ijet1].p4;
            mylljj.jet2_p4 = jets[ijet2].p4;
            mylljj.ll_p4 = ll[ill].p4;
            mylljj.jj_p4 = jj[ijj].p4;
            // gen info
            mylljj.gen_matched = ll[ill].gen_matched && jj[ijj].gen_matched;
            mylljj.gen_p4 = mylljj.gen_matched ? ll[ill].gen_p4 + jj[ijj].gen_p4 : null_p4;
            mylljj.gen_DR = mylljj.gen_matched ? ROOT::Math::VectorUtil::DeltaR(mylljj.p4, mylljj.gen_p4) : -1.;
            mylljj.gen_DPhi = mylljj.gen_matched ? fabs(ROOT::Math::VectorUtil::DeltaPhi(mylljj.p4, mylljj.gen_p4)) : -1.;
            mylljj.gen_DPtOverPt = mylljj.gen_matched ? (mylljj.p4.Pt() - mylljj.gen_p4.Pt()) / mylljj.p4.Pt() : -10.;
            mylljj.gen_lep1_p4 = leptons[ilep1].gen_p4;
            mylljj.gen_lep2_p4 = leptons[ilep2].gen_p4;
            mylljj.gen_jet1_p4 = jets[ijet1].gen_p4;
            mylljj.gen_jet2_p4 = jets[ijet2].gen_p4;
            mylljj.gen_ll_p4 = ll[ill].gen_p4;
            mylljj.gen_jj_p4 = jj[ijj].gen_p4;
            // blind copy of the jj content
            mylljj.ijet1 = jj[ijj].ijet1;
            mylljj.ijet2 = jj[ijj].ijet2;
            mylljj.btag_cMVAv2_MM = jj[ijj].btag_cMVAv2_MM;
            mylljj.btag_deepCSV_MM = jj[ijj].btag_deepCSV_MM;
            mylljj.sumCSV = jj[ijj].sumCSV;
            mylljj.sumCMVAv2 = jj[ijj].sumCMVAv2;
            mylljj.sumDeepCSV = jj[ijj].sumDeepCSV;
            mylljj.DR_j_j = jj[ijj].DR_j_j;
            mylljj.DPhi_j_j = jj[ijj].DPhi_j_j;
            mylljj.ht_j_j = jj[ijj].ht_j_j;
            mylljj.gen_matched_bbPartons = jj[ijj].gen_matched_bbPartons;
            mylljj.gen_matched_bbHadrons = jj[ijj].gen_matched_bbHadrons;
            mylljj.gen_bb = jj[ijj].gen_bb;
            mylljj.gen_bc = jj[ijj].gen_bc;
            mylljj.gen_bl = jj[ijj].gen_bl;
            mylljj.gen_cc = jj[ijj].gen_cc;
            mylljj.gen_cl = jj[ijj].gen_cl;
            mylljj.gen_ll = jj[ijj].gen_ll;
            // blind copy of the ll content
            mylljj.ilep1 = ll[ill].ilep1;
            mylljj.ilep2 = ll[ill].ilep2;
            mylljj.isOS = ll[ill].isOS;
            mylljj.isPlusMinus = ll[ill].isPlusMinus;
            mylljj.isMinusPlus = ll[ill].isMinusPlus;
            mylljj.isMuMu = ll[ill].isMuMu;
            mylljj.isElEl = ll[ill].isElEl;
            mylljj.isElMu = ll[ill].isElMu;
            mylljj.isMuEl = ll[ill].isMuEl;
            mylljj.isSF = ll[ill].isSF;
            mylljj.DR_l_l = ll[ill].DR_l_l;
            mylljj.DPhi_l_l = ll[ill].DPhi_l_l;
            mylljj.ht_l_l = ll[ill].ht_l_l;
            mylljj.trigger_efficiency = ll[ill].trigger_efficiency;
            mylljj.trigger_efficiency_downVariated = ll[ill].trigger_efficiency_downVariated;
            mylljj.trigger_efficiency_upVariated = ll[ill].trigger_efficiency_upVariated;
        
            float DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2;
            DR_j1l1 = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, leptons[ilep1].p4);
            DR_j1l2 = ROOT::Math::VectorUtil::DeltaR(jets[ijet1].p4, leptons[ilep2].p4);
            DR_j2l1 = ROOT::Math::VectorUtil::DeltaR(jets[ijet2].p4, leptons[ilep1].p4);
            DR_j2l2 = ROOT::Math::VectorUtil::DeltaR(jets[ijet2].p4, leptons[ilep2].p4);
            mylljj.maxDR_l_j = std::max({DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2});
            mylljj.minDR_l_j = std::min({DR_j1l1, DR_j1l2, DR_j2l1, DR_j2l2});
            mylljj.DR_ll_jj = ROOT::Math::VectorUtil::DeltaR(ll[ill].p4, jj[ijj].p4);
            mylljj.DPhi_ll_jj = fabs(ROOT::Math::VectorUtil::DeltaPhi(ll[ill].p4, jj[ijj].p4));
            mylljj.MT_fullsystem = mylljj.p4.Mt();
            mylljj.melaAngles = getMELAAngles(ll[ill].p4, jj[ijj].p4, leptons[ilep1].p4, leptons[ilep2].p4, jets[ijet1].p4, jets[ijet2].p4);
            mylljj.visMelaAngles = getMELAAngles(ll[ill].p4, jj[ijj].p4, leptons[ilep1].p4, leptons[ilep2].p4, jets[ijet1].p4, jets[ijet2].p4);
            
            // Compute MT2. See https://arxiv.org/pdf/1309.6318v1.pdf and https://arxiv.org/pdf/1411.4312v5.pdf
            // No need to compute MT in this case, right? It uses invisible momenta that in this case we don't have


            //Counters
            tmp_count_has2leptons_1lljj = event_weight;
            if (mylljj.isElEl)
                tmp_count_has2leptons_elel_1lljj = event_weight;
            if (mylljj.isElMu)
                tmp_count_has2leptons_elmu_1lljj = event_weight;
            if (mylljj.isMuEl)
                tmp_count_has2leptons_muel_1lljj = event_weight;
            if (mylljj.isMuMu)
                tmp_count_has2leptons_mumu_1lljj = event_weight;
            if (mylljj.btag_deepCSV_MM)
            {
                tmp_count_has2leptons_1lljj_2btagM = event_weight;
                if (mylljj.isElEl)
                    tmp_count_has2leptons_elel_1lljj_2btagM = event_weight;
                if (mylljj.isElMu)
                    tmp_count_has2leptons_elmu_1lljj_2btagM = event_weight;
                if (mylljj.isMuEl)
                    tmp_count_has2leptons_muel_1lljj_2btagM = event_weight;
                if (mylljj.isMuMu)
                    tmp_count_has2leptons_mumu_1lljj_2btagM = event_weight;
            }
            // Fill
            lljj.push_back(mylljj);
        } //end of loop over jj
    } //end of loop over ll

    std::sort(lljj.begin(), lljj.end(), [&](HtoZA::DileptonDijet& a, const HtoZA::DileptonDijet& b){ return a.sumCMVAv2 > b.sumCMVAv2; });

    // Keep only the first candidate
    if (lljj.size() > 1){
        lljj.resize(1);
    }


    // ***** ***** *****
    // Event variables
    // ***** ***** *****

    // HT: the two selected leptons - if present - plus all selected jets
    HT = 0;
    if (lljj.size() > 0)
        // take the first lljj since it's been resized: it has only one entry
        HT += lljj[0].lep1_p4.Pt() + lljj[0].lep2_p4.Pt();
    for (unsigned int ijet=0; ijet < jets.size(); ijet++) {
        HT += jets[ijet].p4.Pt();
    }

    nJetsL = jets.size();  //We selected the jets to be all loose
	// Count the number of b-tagged jets
    if (! doingSystematics()) {
        nBJetsM = 0;
        for (unsigned int ijet = 0; ijet < jets.size(); ijet++) {
            if (jets[ijet].btag_deepCSV_M)
                nBJetsM++;
        }


        nMuonsT = 0; // We selected the muons as tight
        nElectronsM = 0;  // We selected the electrons as medium
        for (unsigned int ilepton = 0; ilepton < leptons.size(); ilepton++)
        {
            if (leptons[ilepton].isMu) {
                nMuonsT++;
            }

            if (leptons[ilepton].isEl) {
                nElectronsM++;
            }
        }

        count_has2leptons += tmp_count_has2leptons;
        count_has2leptons_elel += tmp_count_has2leptons_elel;
        count_has2leptons_elmu += tmp_count_has2leptons_elmu;
        count_has2leptons_muel += tmp_count_has2leptons_muel;
        count_has2leptons_mumu += tmp_count_has2leptons_mumu;
        count_has2leptons_1lljj += tmp_count_has2leptons_1lljj;
        count_has2leptons_elel_1lljj += tmp_count_has2leptons_elel_1lljj;
        count_has2leptons_elmu_1lljj += tmp_count_has2leptons_elmu_1lljj;
        count_has2leptons_muel_1lljj += tmp_count_has2leptons_muel_1lljj;
        count_has2leptons_mumu_1lljj += tmp_count_has2leptons_mumu_1lljj;
        count_has2leptons_1lljj_2btagM += tmp_count_has2leptons_1lljj_2btagM;
        count_has2leptons_elel_1lljj_2btagM += tmp_count_has2leptons_elel_1lljj_2btagM;
        count_has2leptons_elmu_1lljj_2btagM += tmp_count_has2leptons_elmu_1lljj_2btagM;
        count_has2leptons_muel_1lljj_2btagM += tmp_count_has2leptons_muel_1lljj_2btagM;
        count_has2leptons_mumu_1lljj_2btagM += tmp_count_has2leptons_mumu_1lljj_2btagM;
    } //end if (! doingSystematics())




    if (!event.isRealData() && !doingSystematics())
	{
// ***** ***** *****
// Get the MC truth information on the hard process
// ***** ***** *****
// from https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/HepMCCandidate/interface/GenStatusFlags.h
//    enum StatusBits {
//0      kIsPrompt = 0,
//1      kIsDecayedLeptonHadron,
//2      kIsTauDecayProduct,
//3      kIsPromptTauDecayProduct,
//4      kIsDirectTauDecayProduct,
//5      kIsDirectPromptTauDecayProduct,
//6      kIsDirectHadronDecayProduct,
//7      kIsHardProcess,
//8      kFromHardProcess,
//9      kIsHardProcessTauDecayProduct,
//10      kIsDirectHardProcessTauDecayProduct,
//11      kFromHardProcessBeforeFSR,
//12      kIsFirstCopy,
//13      kIsLastCopy,
//14      kIsLastCopyBeforeFSR
//    };


    // TTBAR MC TRUTH
    const GenParticlesProducer& gen_particles = producers.get<GenParticlesProducer>("gen_particles");
		
    // 'Pruned' particles are from the hard process
    // 'Packed' particles are stable particles


#if TT_GEN_DEBUG
    std::function<void(size_t)> print_mother_chain = [&gen_particles, &print_mother_chain](size_t p) {

	if (gen_particles.pruned_mothers_index[p].empty()) {
	    std::cout << std::endl;
		return;
	}

	size_t index = gen_particles.pruned_mothers_index[p][0];
	    std::cout << " <- #" << index << "(" << gen_particles.pruned_pdg_id[index] << ")";
		print_mother_chain(index);
    };
#endif

    std::function<bool(size_t, size_t)> pruned_decays_from = [&pruned_decays_from, &gen_particles](size_t particle_index, size_t mother_index) -> bool {
        // Iterator over all pruned particles to find if the particle `particle_index` has `mother_index` in its decay history
        if (gen_particles.pruned_mothers_index[particle_index].empty())
            return false;

        size_t index = gen_particles.pruned_mothers_index[particle_index][0];

        if (index == mother_index) {
            return true;
        }

        if (pruned_decays_from(index, mother_index))
            return true;

        return false;
    };

#define ASSIGN_INDEX( X ) \
    if (flags.isLastCopy()) { \
	    gen_##X = i; \
	}\
	if (flags.isFirstCopy()) { \
	    gen_##X##_beforeFSR = i; \
	}

// Assign index to X if it's empty, or Y if not
#define ASSIGN_INDEX2(X, Y, ERROR) \
    if (flags.isLastCopy()) { \
	    if (gen_##X == 0) \
		    gen_##X = i; \
		else if (gen_##Y == 0)\
		    gen_##Y = i; \
		else \
		    std::cout << ERROR << std::endl; \
	} \
	if (flags.isFirstCopy()) { \
	    if (gen_##X##_beforeFSR == 0) \
		    gen_##X##_beforeFSR = i; \
		else if (gen_##Y##_beforeFSR == 0)\
		    gen_##Y##_beforeFSR = i; \
		else \
		    std::cout << ERROR << std::endl; \
	}

    gen_t = 0; // Index of the top quark
	gen_t_beforeFSR = 0; // Index of the top quark, before any FSR
	gen_tbar = 0; // Index of the anti-top quark
	gen_tbar_beforeFSR = 0; // Index of the anti-top quark, before any FSR

	gen_b = 0; // Index of the b quark coming from the top decay
	gen_b_beforeFSR = 0; // Index of the b quark coming from the top decay, before any FSR
	gen_bbar = 0; // Index of the anti-b quark coming from the anti-top decay
	gen_bbar_beforeFSR = 0; // Index of the anti-b quark coming from the anti-top decay, before any FSR

	gen_jet1_t = 0; // Index of the first jet from the top decay chain
	gen_jet1_t_beforeFSR = 0; // Index of the first jet from the top decay chain, before any FSR
	gen_jet2_t = 0; // Index of the second jet from the top decay chain
	gen_jet2_t_beforeFSR = 0; // Index of the second jet from the top decay chain, before any FSR

	gen_jet1_tbar = 0; // Index of the first jet from the anti-top decay chain
	gen_jet1_tbar_beforeFSR = 0; // Index of the first jet from the anti-top decay chain, before any FSR
	gen_jet2_tbar = 0; // Index of the second jet from the anti-top decay chain
	gen_jet2_tbar_beforeFSR = 0; // Index of the second jet from the anti-top decay chain, before any FSR

	gen_lepton_t = 0; // Index of the lepton from the top decay chain
	gen_lepton_t_beforeFSR = 0; // Index of the lepton from the top decay chain, before any FSR
	gen_neutrino_t = 0; // Index of the neutrino from the top decay chain
	gen_neutrino_t_beforeFSR = 0; // Index of the neutrino from the top decay chain, before any FSR

	gen_lepton_tbar = 0; // Index of the lepton from the anti-top decay chain
	gen_lepton_tbar_beforeFSR = 0; // Index of the lepton from the anti-top decay chain, before any FSR
	gen_neutrino_tbar = 0; // Index of the neutrino from the anti-top decay chain
	gen_neutrino_tbar_beforeFSR = 0; // Index of the neutrino from the anti-top decay chain, before any FSR
	for (size_t i = 0; i < gen_particles.pruned_pdg_id.size(); i++) {

	    int16_t pdg_id = gen_particles.pruned_pdg_id[i];
		uint16_t a_pdg_id = std::abs(pdg_id);

		// We only care of particles with PDG id <= 16 (16 is neutrino tau)
		if (a_pdg_id > 16)
		    continue;

		GenStatusFlags flags(gen_particles.pruned_status_flags[i]);

		if (! flags.isLastCopy() && ! flags.isFirstCopy())
		    continue;

		if (! flags.fromHardProcess())
		    continue;

#if TT_GEN_DEBUG
        std::cout << "---" << std::endl;
	    std::cout << "Gen particle #" << i << ": PDG id: " << gen_particles.pruned_pdg_id[i];
	    print_mother_chain(i);
	    flags.dump();
#endif

        if (pdg_id == 6) {
		    ASSIGN_INDEX(t);
			continue;
		} else if (pdg_id == -6) {
		    ASSIGN_INDEX(tbar);
			continue;
		}

		if (gen_t == 0 || gen_tbar == 0) {
		    // Don't bother if we don't have found the tops
			continue;
		}

		bool from_t_decay = pruned_decays_from(i, gen_t);
		bool from_tbar_decay = pruned_decays_from(i, gen_tbar);

		// Only keep particles coming from the tops decay
		if (! from_t_decay && ! from_tbar_decay)
		    continue;

		if (pdg_id == 5) {
		    // Maybe it's a b coming from the W decay
			if (!flags.isFirstCopy() && flags.isLastCopy() && gen_b == 0) {

			    // This can be a B decaying from a W
				// However, we can't rely on the presence of the W in the decay chain, as it may be generator specific
				// Since it's the last copy (ie, after FSR), we can check if this B comes from the B assigned to the W decay (ie, gen_jet1_t_beforeFSR, gen_jet2_t_beforeFSR)
				// If yes, then it's not the B coming directly from the top decay
				if ((gen_jet1_t_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet1_t_beforeFSR]) == 5) ||
				    (gen_jet2_t_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet2_t_beforeFSR]) == 5) ||
					(gen_jet1_tbar_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet1_tbar_beforeFSR]) == 5) ||
					(gen_jet2_tbar_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet2_tbar_beforeFSR]) == 5)) {


#if TT_GEN_DEBUG
                    std::cout << "A quark coming from W decay is a b" << std::endl;
#endif

                    if (! (gen_jet1_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_tbar_beforeFSR)) &&
					    ! (gen_jet2_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_tbar_beforeFSR)) &&
						! (gen_jet1_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_t_beforeFSR)) &&
						! (gen_jet2_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_t_beforeFSR))) {
#if TT_GEN_DEBUG
                        std::cout << "This after-FSR b quark is not coming from a W decay" << std::endl;
#endif
                        gen_b = i;
					    continue;
				    }
#if TT_GEN_DEBUG
                    else {
					    std::cout << "This after-FSR b quark comes from a W decay" << std::endl;
					}
#endif
                } else {
#if TT_GEN_DEBUG
                    std::cout << "Assigning gen_b" << std::endl;
#endif
                    gen_b = i;
					continue;
				}
			} else if (flags.isFirstCopy() && gen_b_beforeFSR == 0) {
			    gen_b_beforeFSR = i;
				continue;
			} else {
#if TT_GEN_DEBUG
                std::cout << "This should not happen!" << std::endl;
#endif
            }
		} else if (pdg_id == -5) {
		    if (!flags.isFirstCopy() && flags.isLastCopy() && gen_bbar == 0) {

			    // This can be a B decaying from a W
				// However, we can't rely on the presence of the W in the decay chain, as it may be generator specific
				// Since it's the last copy (ie, after FSR), we can check if this B comes from the B assigned to the W decay (ie, gen_jet1_t_beforeFSR, gen_jet2_t_beforeFSR)
				// If yes, then it's not the B coming directly from the top decay
				if ((gen_jet1_t_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet1_t_beforeFSR]) == 5) ||
				    (gen_jet2_t_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet2_t_beforeFSR]) == 5) ||
					(gen_jet1_tbar_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet1_tbar_beforeFSR]) == 5) ||
					(gen_jet2_tbar_beforeFSR != 0 && std::abs(gen_particles.pruned_pdg_id[gen_jet2_tbar_beforeFSR]) == 5)) {

#if TT_GEN_DEBUG
                    std::cout << "A quark coming from W decay is a bbar" << std::endl;
#endif

                    if (! (gen_jet1_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_tbar_beforeFSR)) &&
					    ! (gen_jet2_tbar_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_tbar_beforeFSR)) &&
						! (gen_jet1_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet1_t_beforeFSR)) &&
						! (gen_jet2_t_beforeFSR != 0 && pruned_decays_from(i, gen_jet2_t_beforeFSR))) {
#if TT_GEN_DEBUG
                        std::cout << "This after-fsr b anti-quark is not coming from a W decay" << std::endl;
#endif
                        gen_bbar = i;
						continue;
					}
#if TT_GEN_DEBUG
                    else {
					    std::cout << "This after-fsr b anti-quark comes from a W decay" << std::endl;
					}
#endif
                } else {
#if TT_GEN_DEBUG
                    std::cout << "Assigning gen_bbar" << std::endl;
#endif
                    gen_bbar = i;
					continue;
				}
			} else if (flags.isFirstCopy() && gen_bbar_beforeFSR == 0) {
			    gen_bbar_beforeFSR = i;
				continue;
			}
		}

		if ((gen_tbar == 0) || (gen_t == 0))
		    continue;

		if (gen_t != 0 && from_t_decay) {
#if TT_GEN_DEBUG
        std::cout << "Coming from the top chain decay" << std::endl;
#endif
            if (a_pdg_id >= 1 && a_pdg_id <= 5) {
			    ASSIGN_INDEX2(jet1_t, jet2_t, "Error: more than two quarks coming from top decay");
			} else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
			    ASSIGN_INDEX(lepton_t);
		    } else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
			    ASSIGN_INDEX(neutrino_t);
			} else {
			    std::cout << "Error: unknown particle coming from top decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
			}
		} else if (gen_tbar != 0 && from_tbar_decay) {
#if TT_GEN_DEBUG
        std::cout << "Coming from the anti-top chain decay" << std::endl;
#endif
            if (a_pdg_id >= 1 && a_pdg_id <= 5) {
			    ASSIGN_INDEX2(jet1_tbar, jet2_tbar, "Error: more than two quarks coming from anti-top decay");
		    } else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
			    ASSIGN_INDEX(lepton_tbar);
			} else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
			    ASSIGN_INDEX(neutrino_tbar);
		    } else {
			    std::cout << "Error: unknown particle coming from anti-top decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
			}
		}
	}

	if (!gen_t || !gen_tbar) {
#if TT_GEN_DEBUG
        std::cout << "This is not a ttbar event" << std::endl;
#endif
        gen_ttbar_decay_type = NotTT;
		return;
	}

	if ((gen_jet1_t != 0) && (gen_jet2_t != 0) && (gen_jet1_tbar != 0) && (gen_jet2_tbar != 0)) {
#if TT_GEN_DEBUG
        std::cout << "Hadronic ttbar decay" << std::endl;
#endif
        gen_ttbar_decay_type = Hadronic;
	} else if (
	        ((gen_lepton_t != 0) && (gen_lepton_tbar == 0)) ||
			((gen_lepton_t == 0) && (gen_lepton_tbar != 0))
			) {

#if TT_GEN_DEBUG
        std::cout << "Semileptonic ttbar decay" << std::endl;
#endif

        uint16_t lepton_pdg_id;
		if (gen_lepton_t != 0)
		    lepton_pdg_id = std::abs(gen_particles.pruned_pdg_id[gen_lepton_t]);
		else
		    lepton_pdg_id = std::abs(gen_particles.pruned_pdg_id[gen_lepton_tbar]);

		if (lepton_pdg_id == 11)
		    gen_ttbar_decay_type = Semileptonic_e;
		else if (lepton_pdg_id == 13)
		    gen_ttbar_decay_type = Semileptonic_mu;
		else
		    gen_ttbar_decay_type = Semileptonic_tau;
	} else if (gen_lepton_t != 0 && gen_lepton_tbar != 0) {
	    uint16_t lepton_t_pdg_id = std::abs(gen_particles.pruned_pdg_id[gen_lepton_t]);
		uint16_t lepton_tbar_pdg_id = std::abs(gen_particles.pruned_pdg_id[gen_lepton_tbar]);

#if TT_GEN_DEBUG
        std::cout << "Dileptonic ttbar decay" << std::endl;
#endif

        if (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 11)
		    gen_ttbar_decay_type = Dileptonic_ee;
	    else if (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 13)
		    gen_ttbar_decay_type = Dileptonic_mumu;
	    else if (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 15)
		    gen_ttbar_decay_type = Dileptonic_tautau;
		else if (
		        (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 13) ||
				(lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 11)
				) {
			gen_ttbar_decay_type = Dileptonic_mue;
		}
		else if (
		        (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 15) ||
				(lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 11)
				) {
			gen_ttbar_decay_type = Dileptonic_etau;
		}
		else if (
		        (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 15) ||
				(lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 13)
				) {
			gen_ttbar_decay_type = Dileptonic_mutau;
		} else {
		    std::cout << "Error: unknown dileptonic ttbar decay." << std::endl;
			gen_ttbar_decay_type = NotTT;
			return;
		}
	} else {
	    std::cout << "Error: unknown ttbar decay." << std::endl;
		gen_ttbar_decay_type = UnknownTT;
	}
	} // end of if !event.isRealData()

}

void HtoZAAnalyzer::endJob(MetadataManager& metadata) {
    if (! doingSystematics()) {
        metadata.add(this->m_name + "_count_has2leptons", count_has2leptons);
        metadata.add(this->m_name + "_count_has2leptons_elel", count_has2leptons_elel);
        metadata.add(this->m_name + "_count_has2leptons_elmu", count_has2leptons_elmu);
        metadata.add(this->m_name + "_count_has2leptons_muel", count_has2leptons_muel);
        metadata.add(this->m_name + "_count_has2leptons_mumu", count_has2leptons_mumu);
        metadata.add(this->m_name + "_count_has2leptons_1lljj", count_has2leptons_1lljj);
        metadata.add(this->m_name + "_count_has2leptons_elel_1lljj", count_has2leptons_elel_1lljj);
        metadata.add(this->m_name + "_count_has2leptons_elmu_1lljj", count_has2leptons_elmu_1lljj);
        metadata.add(this->m_name + "_count_has2leptons_muel_1lljj", count_has2leptons_muel_1lljj);
        metadata.add(this->m_name + "_count_has2leptons_mumu_1lljj", count_has2leptons_mumu_1lljj);
        metadata.add(this->m_name + "_count_has2leptons_1lljj_2btagM", count_has2leptons_1lljj_2btagM);
        metadata.add(this->m_name + "_count_has2leptons_elel_1lljj_2btagM", count_has2leptons_elel_1lljj_2btagM);
        metadata.add(this->m_name + "_count_has2leptons_elmu_1lljj_2btagM", count_has2leptons_elmu_1lljj_2btagM);
        metadata.add(this->m_name + "_count_has2leptons_muel_1lljj_2btagM", count_has2leptons_muel_1lljj_2btagM);
        metadata.add(this->m_name + "_count_has2leptons_mumu_1lljj_2btagM", count_has2leptons_mumu_1lljj_2btagM);
    }
}
