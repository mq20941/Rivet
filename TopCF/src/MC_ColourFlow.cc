// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Colour flow in hadronic top decay at 13 TeV
  class MC_ColourFlow : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_ColourFlow);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections of particles/jets before the run
    void init() { 

      const FinalState fs;

      PromptFinalState promptFs(fs);
      promptFs.acceptTauDecays(true); //Accept leptons from decays of prompt taus as themselves being prompt
      promptFs.acceptMuonDecays(false); //Accept leptons from decays of prompt muons as themselves being prompt

      IdentifiedFinalState neutrino_fs(promptFs); //accept all neutrinos in final state
      neutrino_fs.acceptNeutrinos();
      declare(neutrino_fs, "NEUTRINO_FS");

      IdentifiedFinalState Photon(fs); 
      Photon.acceptIdPair(PID::PHOTON); //accept photons in final state

      IdentifiedFinalState bare_muons_fs(promptFs); 
      bare_muons_fs.acceptIdPair(PID::MUON); //accept muons and antimuons in final state

      IdentifiedFinalState bare_elecs_fs(promptFs);
      bare_elecs_fs.acceptIdPair(PID::ELECTRON); //accept electrons and positrons in final state

      Cut lep_cuts = (Cuts::abseta < 2.5) & (Cuts::pT > 1*MeV);
      DressedLeptons muons(Photon, bare_muons_fs, 0.1, lep_cuts); //select muons
      declare(muons, "MUONS");

      DressedLeptons elecs(Photon, bare_elecs_fs, 0.1, lep_cuts); //select electrons
      declare(elecs, "ELECS");

      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(muons); //veto all muons which don't satisfy lep_cuts from final state
      vfs.addVetoOnThisFinalState(elecs); //veto all electrons which don't satisfy lep_cuts from final state
      vfs.addVetoOnThisFinalState(neutrino_fs);

      FastJets fjets(vfs, FastJets::ANTIKT, 0.4);
      fjets.useInvisibles(); //include invinisible particles in jet construction
      declare(fjets, "jets"); 

      //book histograms
      book(h_pull_12_all, "pull_12_all", 6, 0, 1);
      book(h_pull_21_all, "pull_21_all", 6, 0, 1);
      book(h_pull_both_all, "pull_both_all", 6, 0, 1);
      book(h_bpull_12_all, "bpull_12_all", 6, 0, 1);
      book(h_bpull_21_all, "bpull_21_all", 6, 0, 1);
      book(h_pull_angle_12_jetpTbins, "pull_angle_12_jetpTbins", 10, 20, 120);
      book(h_pull_angle_21_jetpTbins, "pull_angle_21_jetpTbins", 10, 20, 120);
      book(h_pull_12_all_lowpT, "pull_12_all_lowpT", 6, 0, 1); 
      book(h_pull_12_all_highpT, "pull_12_all_highpT", 6, 0, 1);
      book(h_pull_12_all_veryhighpT, "pull_12_all_veryhighpT", 6, 0, 1);
      book(h_wJet1_pT, "wJet1_pT", 6,  20, 100); //20 GeV cut on all jets
      book(h_wJet2_pT, "wJet2_pT", 6,  20, 100);
      book(h_bJet1_pT, "bJet1_pT", 6,  20, 140);
      book(h_bJet2_pT, "bJet2_pT", 6,  20, 140);
      book(h_pullmagnitudej1_all, "pullmagnitudej1_all", 5, 0, 0.01);
      book(h_bpull_both_all, "bpull_both_all", 6, 0, 1);
      book(h_pullvectorj1_perpproj_all, "pullvectorj1_perpproj_all", 5, 0, 0.01);
      book(h_pullvectorj1_paraproj_all, "pullvectorj1_paraproj_all", 5, 0, 0.01);
      book(h_pull_both_all_lowpT, "pull_both_all_lowpT", 6, 0, 1);
      book(h_pull_both_all_highpT, "pull_both_all_highpT", 6, 0, 1);
      book(h_pull_both_all_veryhighpT, "pull_both_all_veryhighpT", 6, 0, 1);
      book(h_bpullmagnitudej1_all, "bpullmagnitudej1_all", 5, 0, 0.01);
      book(h_leading_w_b_pull_all, "leading_w_b_pull_all", 5, 0, 1);
      book(h_pull_12_w, "pull_12_w", 5, 0, 1);
      book(h_pull_12_lowpT_w, "pull_12_lowpT_w", 5, 0, 1);
      book(h_pull_12_highpT_w, "pull_12_highpT_w", 5, 0, 1);
      
      //book profiles
      book(p_paraproj_asymm, "paraproj_asymm", 3, 0, 0.01);
      book(p_paraproj_asymm_tail, "paraproj_asymm_tail", 6, 0, 0.02);
      book(p_paraproj_asymm_onebin, "paraproj_asymm_onebin", 1, 0, 0.01);
      book(p_perpproj_asymm, "perpproj_asymm", 3, 0, 0.01);
      book(p_perpproj_asymm_onebin, "perpproj_asymm_onebin", 1, 0, 0.01);
      book(p_paraproj_asymm_M, "paraproj_asymm_M", 3, 0, 0.01);
      book(p_perpproj_asymm_M, "perpproj_asymm_M", 3, 0, 0.01);
      book(p_pull_angle_paraproj_asymm_cor, "pull_angle_paraproj_asymm_cor", 3, 0, 0.5);
      book(p_pull_angle_paraproj_asymm_cor_binone, "pull_angle_paraproj_asymm_cor_binone", 1, 0, 0.16666667);
      book(p_pullmagnitude_paraproj_asymm_cor, "pullmagnitude_paraproj_asymm_cor", 3, 0, 0.01);
      book(p_paraproj_pull_angle_cor, "paraproj_pull_angle_cor", 3, 0, 0.01);
      book(p_paraproj_pull_angle_cor_signed, "paraproj_pull_angle_cor_signed", 6, -0.01, 0.01);
      book(p_paraproj_pullmagnitude_cor, "paraproj_pullmagnitude_cor", 3, 0, 0.01);
      book(p_pull_angle_asymm, "pull_angle_asymm", 3, 0, 0.5);
      book(p_bparaproj_asymm, "bparaproj_asymm", 3, 0, 0.01);
      book(p_bperpproj_asymm, "bperpproj_asymm", 3, 0, 0.01);
      book(p_bparaproj_asymm_onebin, "bparaproj_asymm_onebin", 1, 0, 0.01);
      book(p_bperpproj_asymm_onebin, "bperpproj_asymm_onebin", 1, 0, 0.01);
      book(p_paraproj_asymm_j2, "paraproj_asymm_j2", 3, 0, 0.01);
      book(p_paraproj_asymm_j2_onebin, "paraproj_asymm_j2_onebin", 1, 0, 0.01);
      book(p_paraproj_asymm_jetpTbins, "paraproj_asymm_jetpTbins", 3, 20, 100);
      book(p_paraproj_asymm_j2_jetpTbins, "paraproj_asymm_j2_jetpTbins", 3, 20, 100);
      book(p_jetpT_pullmagnitude_cor, "jetpT_pullmagnitude_cor", 3, 20, 100);
      book(p_jetpT_paraproj_cor, "jetpT_paraproj_cor", 3, 20, 100);
      book(p_paraproj_asymm_lowpT, "paraproj_asymm_lowpT", 3, 0, 0.01);
      book(p_paraproj_asymm_highpT, "paraproj_asymm_highpT", 3, 0, 0.01);
      book(p_paraproj_asymm_veryhighpT, "paraproj_asymm_veryhighpT", 3, 0, 0.01);
      book(p_paraproj_asymm_onebin_lowpT, "paraproj_asymm_onebin_lowpT", 1, 0, 0.01);
      book(p_paraproj_asymm_onebin_highpT, "paraproj_asymm_onebin_highpT", 1, 0, 0.01);
      book(p_paraproj_asymm_onebin_veryhighpT, "paraproj_asymm_onebin_veryhighpT", 1, 0, 0.01);
      book(p_paraproj_asymm_onebin_jetpTbins, "paraproj_asymm_onebin_jetpTbins", 5, 20, 140);
      book(p_paraproj_asymm_j2_onebin_jetpTbins, "paraproj_asymm_j2_onebin_jetpTbins", 5, 20, 140);
      book(p_paraproj_asymm_onebin_jetpTbins_twobins, "paraproj_asymm_onebin_jetpTbins_twobins", 2, 20, 100);
      book(p_paraproj_asymm_j2_onebin_jetpTbins_twobins, "paraproj_asymm_j2_onebin_jetpTbins_twobins", 2, 20, 100);
      book(p_paraproj_asymm_both, "paraproj_asymm_both", 3, 0, 0.01);
      book(p_paraproj_asymm_both_onebin, "paraproj_asymm_both_onebin", 1, 0, 0.01);
      book(p_bpull_angle_paraproj_asymm_cor, "bpull_angle_paraproj_asymm_cor", 3, 0, 0.5);
      book(p_bparaproj_asymm_onebin_jetpTbins, "bparaproj_asymm_onebin_jetpTbins", 5, 20, 140);
      book(p_bparaproj_asymm_onebin_M, "bparaproj_asymm_onebin_M", 1, 0, 0.01);
      book(p_bpull_12_21_cor, "bpull_12_21_cor", 6, 0, 1);
      book(p_wjetpurity, "wjetpurity", 1, 0, 2);
      book(p_wjetpurity_pTbins, "wjetpurity_pTbins", 2, 20, 100);

    
    } 

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /**************
       *    JETS    *
       **************/
      const Jets& allJets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20.0*GeV && Cuts::absrap < 2.5); //jet pT,rap cut (high energy jets, within detector)
      const vector<DressedLepton>& all_elecs = apply<DressedLeptons>(event, "ELECS").dressedLeptons();
      const vector<DressedLepton>& all_muons = apply<DressedLeptons>(event, "MUONS").dressedLeptons();
      Jets goodJets;
      for (const Jet j : allJets) {
        bool keep = true;
        for (const DressedLepton el : all_elecs)  keep &= deltaR(j, el) >= 0.2;
        if (keep)  goodJets += j;
      }
      if ( goodJets.size() < 4 )  vetoEvent; //need at least 4 jets (so we have semileptonic decay channel)
      
      /****************
       *    LEPTONS   *
       ****************/
      vector<DressedLepton> muons, vetoMuons;
      for (const DressedLepton mu : all_muons) {
        bool keep = true;
        for (const Jet j : goodJets)  keep &= deltaR(j, mu) >= 0.4;
        if (keep && mu.pt() > 15*GeV) {
          vetoMuons.push_back(mu);
          if (mu.pt() > 20*GeV)  muons.push_back(mu); //muon pT cut 
        } //vetoMuons looser cut than muons, muon subsection of vetomuons
      }

      vector<DressedLepton> elecs, vetoElecs;
      for (const DressedLepton el : all_elecs) {
        bool keep = true;
        for (const Jet j : goodJets)  keep &= deltaR(j, el) >= 0.4;
        if (keep && el.pt() > 15*GeV) {
          vetoElecs.push_back(el);
          if (el.pt() > 20*GeV)  elecs.push_back(el); //electron pT cut
        }
      }
      //&& is and operator, || is or operator
      if (muons.empty() && elecs.empty())  vetoEvent; //no muons or leptons -> veto event 

      bool muCandidate = !( muons.size() < 1 || vetoMuons.size() > 1 || vetoElecs.size() > 0 ); //exactly 1 high pt muon, no electrons
      bool elCandidate = !( elecs.size() < 1 || vetoElecs.size() > 1 || vetoMuons.size() > 0 ); //exactly 1 high pt electron, no muons

      if (!elCandidate && !muCandidate)  vetoEvent;

      /******************************
       *    ELECTRON-MUON OVERLAP   *
       ******************************/
      for (const DressedLepton electron : elecs) {
        for (const DressedLepton muon : muons) {
          double d_theta = fabs(muon.theta() - electron.theta());
          double d_phi = deltaPhi(muon.phi(), electron.phi());
          if (d_theta < 0.005 && d_phi < 0.005)  vetoEvent;
        }
      }

      /****************
       *  NEUTRINOS   *
       ****************/
      const Particles& neutrinos = apply<IdentifiedFinalState>(event, "NEUTRINO_FS").particlesByPt();
      FourMomentum metVector = FourMomentum(0.,0.,0.,0.);
      for (const Particle& n : neutrinos) {
        metVector += n.momentum();
      }
      double met = metVector.pt();
      if (met <= 20*GeV)  vetoEvent;
      //working out missing transverse energy, at least 20GeV 
      
      //if ( (_mT(muCandidate? muons[0] : elecs[0], metVector) + met) <= 60. )  vetoEvent; 
      
      /****************
       *    B-JETS    *
       ****************/
      // Sort the jets into b-jets and light jets. We expect one hard b-jet from
      // each top decay, so our 4 hardest jets should include two b-jets. The
      // Jet::bTagged() method is equivalent to perfect experimental
      // b-tagging, in a generator-independent way.
      Jets bJets, wJets;
      for(Jet j : goodJets) {
        bool b_tagged = false;
        Particles bTags = j.bTags();
        for ( Particle b : bTags ) {
          b_tagged |= b.pT() > 5*GeV; //5 GeV pT cut on b particles
        }
        if (b_tagged)  bJets += j; //appending into bJet
        if (!b_tagged) wJets += j; //appending into wJet (assumption that if not b jet, is w jet)
      }

      if ( bJets.size() < 2 || wJets.size() < 2 )  vetoEvent; // need 2 of each jet in an event for semileptonic ttbar 

      Jets truewJets;
      for(Jet j : wJets) {
        bool w_tagged = false;
        Particles& constituents = j.particles();
        for (Particle p : constituents) {
          w_tagged |= p.hasAncestor(24, false);
          w_tagged |= p.hasAncestor(-24, false);
        }
        if (w_tagged) truewJets += j;
      }
      
      if (truewJets.size() > 0) {
        for (size_t i = 0; i < truewJets.size(); ++i) {
          double purity = CalculateJetPurity(truewJets[i]);
          p_wjetpurity->fill(1, purity);
          p_wjetpurity_pTbins->fill(truewJets[i].pT(), purity);
        }
      }

      if(truewJets.size() > 1) {
        for(size_t i = 0; i < truewJets.size()-1; ++i) {
          for(size_t j = i+1; j < truewJets.size(); ++j) {
            double pull_12_w = CalculatePullAngle(truewJets[i], truewJets[j], 0);
            h_pull_12_w->fill(pull_12_w/Rivet::PI);
            if(truewJets[i].pT() < 60*GeV) {
              double pull_12_lowpT_w = CalculatePullAngle(truewJets[i], truewJets[j], 0);
              h_pull_12_lowpT_w->fill(pull_12_lowpT_w/Rivet::PI);
            }
            if(truewJets[i].pT() > 60*GeV) {
              double pull_12_highpT_w = CalculatePullAngle(truewJets[i], truewJets[j], 0);
              h_pull_12_highpT_w->fill(pull_12_highpT_w/Rivet::PI);
            }
          }
        }
      }

     
      // Construct the hadronically decaying W momentum 4-vector from pairs of
      // non-b-tagged jets. The pair(s) within 20GeV of W mass is used. 
      for (size_t i = 0; i < wJets.size()-1; ++i) { 
       for (size_t j = i + 1; j < wJets.size(); ++j) {
         double Wjetsinvariantmass = fabs(CalculateInvariantMass(wJets[i], wJets[j]));
         
         if (Wjetsinvariantmass > 60.4*GeV && Wjetsinvariantmass < 100.4*GeV) {
           
           double pull_angle_12_all = fabs(CalculatePullAngle(wJets[i], wJets[j], 0)); 
           h_pull_12_all->fill(pull_angle_12_all / Rivet::PI);
           double pull_angle_21_all = fabs(CalculatePullAngle(wJets[j], wJets[i], 0));
           h_pull_21_all->fill(pull_angle_21_all / Rivet::PI); 
       
           h_pull_both_all->fill(pull_angle_12_all / Rivet::PI); 
           h_pull_both_all->fill(pull_angle_21_all / Rivet::PI); 
           
           double pullvectorj1_all = fabs(CalculatePullMagnitude(wJets[i], 0)); 
           h_pullmagnitudej1_all->fill(pullvectorj1_all);
           double pullvectorj1_perpproj_all = fabs(CalculatePullPerpProjection(wJets[i], wJets[j], 0));
           double pullvectorj1_perpproj_all_signed = CalculatePullPerpProjectionSigned(wJets[i], wJets[j], 0);
           double pullvectorj1_paraproj_all = fabs(CalculatePullParaProjection(wJets[i], wJets[j], 0));
           double pullvectorj1_paraproj_all_signed = CalculatePullParaProjectionSigned(wJets[i], wJets[j], 0);

           h_pullvectorj1_perpproj_all->fill(pullvectorj1_perpproj_all);
           h_pullvectorj1_paraproj_all->fill(pullvectorj1_paraproj_all);

           p_perpproj_asymm->fill(fabs(pullvectorj1_perpproj_all_signed), double(pullvectorj1_perpproj_all_signed > 0) - double(pullvectorj1_perpproj_all_signed < 0));
           p_perpproj_asymm_M->fill(fabs(pullvectorj1_perpproj_all_signed), 1000*pullvectorj1_perpproj_all_signed);
           p_perpproj_asymm_onebin->fill(fabs(pullvectorj1_perpproj_all_signed), double(pullvectorj1_perpproj_all_signed > 0) - double(pullvectorj1_perpproj_all_signed < 0));
           p_paraproj_asymm->fill(fabs(pullvectorj1_paraproj_all_signed), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
           p_paraproj_asymm_tail->fill(fabs(pullvectorj1_paraproj_all_signed), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
           p_paraproj_asymm_onebin->fill(fabs(pullvectorj1_paraproj_all_signed), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
           p_paraproj_asymm_M->fill(fabs(pullvectorj1_paraproj_all_signed), 1000*pullvectorj1_paraproj_all_signed); 

           double pull_angle_12_all_modhalfpi = CalculatePullAngle_ModHalfPI(wJets[i], wJets[j], 0); 
           p_pull_angle_paraproj_asymm_cor->fill((pull_angle_12_all_modhalfpi / Rivet::PI), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
           p_pullmagnitude_paraproj_asymm_cor->fill(pullvectorj1_all, double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
           
           p_paraproj_pull_angle_cor->fill(pullvectorj1_paraproj_all, (pull_angle_12_all_modhalfpi / Rivet::PI));
           p_paraproj_pull_angle_cor_signed->fill(pullvectorj1_paraproj_all_signed, (pull_angle_12_all/ Rivet::PI));
           p_paraproj_pullmagnitude_cor->fill(pullvectorj1_paraproj_all, pullvectorj1_all);
           
           p_pull_angle_asymm->fill((pull_angle_12_all_modhalfpi / Rivet::PI), double((pull_angle_12_all / Rivet::PI) < 0.5) - ((pull_angle_12_all / Rivet::PI) > 0.5)); 
           
           if ((pull_angle_12_all / Rivet::PI) < 0.3333333) {
             h_pull_angle_12_jetpTbins->fill(wJets[i].pT());
           }

           if ((pull_angle_21_all / Rivet::PI) < 0.3333333) {
             h_pull_angle_21_jetpTbins->fill(wJets[j].pT());
           }

           if ((pull_angle_12_all_modhalfpi / Rivet::PI) < 0.166666666667){
             p_pull_angle_paraproj_asymm_cor_binone->fill((pull_angle_12_all_modhalfpi / Rivet::PI), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
           }
           
           double pullvectorj2_paraproj_all_signed = CalculatePullParaProjectionSigned(wJets[j], wJets[i], 0);
           p_paraproj_asymm_j2->fill(fabs(pullvectorj2_paraproj_all_signed), double(pullvectorj2_paraproj_all_signed > 0) - double(pullvectorj2_paraproj_all_signed < 0));
           p_paraproj_asymm_j2_onebin->fill(fabs(pullvectorj2_paraproj_all_signed), double(pullvectorj2_paraproj_all_signed > 0) - double(pullvectorj2_paraproj_all_signed < 0));
           
           p_paraproj_asymm_both->fill(fabs(pullvectorj1_paraproj_all_signed), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
           p_paraproj_asymm_both->fill(fabs(pullvectorj2_paraproj_all_signed), double(pullvectorj2_paraproj_all_signed > 0) - double(pullvectorj2_paraproj_all_signed < 0));
           p_paraproj_asymm_both_onebin->fill(fabs(pullvectorj1_paraproj_all_signed), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
           p_paraproj_asymm_both_onebin->fill(fabs(pullvectorj2_paraproj_all_signed), double(pullvectorj2_paraproj_all_signed > 0) - double(pullvectorj2_paraproj_all_signed < 0));

           if (fabs(pullvectorj1_paraproj_all_signed) < 0.003333333) {
             p_paraproj_asymm_jetpTbins->fill(wJets[i].pT(), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
           }

           if (fabs(pullvectorj2_paraproj_all_signed) < 0.003333333) {
             p_paraproj_asymm_j2_jetpTbins->fill(wJets[j].pT(), double(pullvectorj2_paraproj_all_signed > 0) - double(pullvectorj2_paraproj_all_signed < 0));
           }

           p_paraproj_asymm_onebin_jetpTbins->fill(wJets[i].pT(), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
           p_paraproj_asymm_j2_onebin_jetpTbins->fill(wJets[j].pT(), double(pullvectorj2_paraproj_all_signed > 0) - double(pullvectorj2_paraproj_all_signed < 0));
           p_paraproj_asymm_onebin_jetpTbins_twobins->fill(wJets[i].pT(), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
           p_paraproj_asymm_j2_onebin_jetpTbins_twobins->fill(wJets[j].pT(), double(pullvectorj2_paraproj_all_signed > 0) - double(pullvectorj2_paraproj_all_signed < 0));

           p_jetpT_pullmagnitude_cor->fill(wJets[i].pT(), pullvectorj1_all);
           p_jetpT_paraproj_cor->fill(wJets[i].pT(), pullvectorj1_paraproj_all);
           
           if (wJets[i].pT() < 60*GeV){
             p_paraproj_asymm_lowpT->fill(fabs(pullvectorj1_paraproj_all_signed), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
             p_paraproj_asymm_onebin_lowpT->fill(fabs(pullvectorj1_paraproj_all_signed), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
             h_pull_12_all_lowpT->fill(pull_angle_12_all / Rivet::PI);
             
             h_pull_both_all_lowpT->fill(pull_angle_12_all / Rivet::PI);
             h_pull_both_all_lowpT->fill(pull_angle_21_all / Rivet::PI);
           }

           if (wJets[i].pT() > 60*GeV){
             p_paraproj_asymm_highpT->fill(fabs(pullvectorj1_paraproj_all_signed), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
             p_paraproj_asymm_onebin_highpT->fill(fabs(pullvectorj1_paraproj_all_signed), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
             h_pull_12_all_highpT->fill(pull_angle_12_all / Rivet::PI);

             h_pull_both_all_highpT->fill(pull_angle_12_all / Rivet::PI);
             h_pull_both_all_highpT->fill(pull_angle_21_all / Rivet::PI);
           }

           if (wJets[i].pT() > 80*GeV){
             p_paraproj_asymm_veryhighpT->fill(fabs(pullvectorj1_paraproj_all_signed), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
             p_paraproj_asymm_onebin_veryhighpT->fill(fabs(pullvectorj1_paraproj_all_signed), double(pullvectorj1_paraproj_all_signed > 0) - double(pullvectorj1_paraproj_all_signed < 0));
             h_pull_12_all_veryhighpT->fill(pull_angle_12_all / Rivet::PI);

             h_pull_both_all_veryhighpT->fill(pull_angle_12_all / Rivet::PI);
             h_pull_both_all_veryhighpT->fill(pull_angle_21_all / Rivet::PI);
           }

           h_wJet1_pT->fill(wJets[i].pT());
           h_wJet2_pT->fill(wJets[j].pT());
         }


       
        }
      }

     //non W-jets plots

     double bpull_angle_12_all = fabs(CalculatePullAngle(bJets[0], bJets[1], 0)); //b jet pull angle (1,2)
     h_bpull_12_all->fill(bpull_angle_12_all / Rivet::PI);
     double bpull_angle_21_all = fabs(CalculatePullAngle(bJets[1], bJets[0], 0)); //b jet pull angle (2, 1)
     h_bpull_21_all->fill(bpull_angle_21_all / Rivet::PI);

     h_bpull_both_all->fill(bpull_angle_12_all / Rivet::PI);
     h_bpull_both_all->fill(bpull_angle_21_all / Rivet::PI);

     p_bpull_12_21_cor->fill((bpull_angle_12_all / Rivet::PI), (bpull_angle_21_all / Rivet::PI));

     double pullvectorb1_all = CalculatePullMagnitude(bJets[0], 0);
     h_bpullmagnitudej1_all->fill(fabs(pullvectorb1_all));

     double pullvectorb1_paraproj_all_signed = CalculatePullParaProjectionSigned(bJets[0], bJets[1], 0);
     double pullvectorb1_perpproj_all_signed = CalculatePullPerpProjectionSigned(bJets[0], bJets[1], 0);
     p_bparaproj_asymm->fill(fabs(pullvectorb1_paraproj_all_signed), double(pullvectorb1_paraproj_all_signed > 0) - double(pullvectorb1_paraproj_all_signed < 0));
     p_bparaproj_asymm_onebin->fill(fabs(pullvectorb1_paraproj_all_signed), double(pullvectorb1_paraproj_all_signed > 0) - double(pullvectorb1_paraproj_all_signed < 0));
     p_bperpproj_asymm->fill(fabs(pullvectorb1_perpproj_all_signed), double(pullvectorb1_perpproj_all_signed > 0) - double(pullvectorb1_perpproj_all_signed < 0));
     p_bperpproj_asymm_onebin->fill(fabs(pullvectorb1_perpproj_all_signed), double(pullvectorb1_perpproj_all_signed > 0) - double(pullvectorb1_perpproj_all_signed < 0));
     p_bparaproj_asymm_onebin_M->fill(fabs(pullvectorb1_paraproj_all_signed), 1000*pullvectorb1_paraproj_all_signed);

     double bpull_angle_12_all_modhalfpi = CalculatePullAngle_ModHalfPI(bJets[0], bJets[1], 0); 
     p_bpull_angle_paraproj_asymm_cor->fill((bpull_angle_12_all_modhalfpi / Rivet::PI), double(pullvectorb1_paraproj_all_signed > 0) - double(pullvectorb1_paraproj_all_signed < 0));
     
     p_bparaproj_asymm_onebin_jetpTbins->fill(bJets[0].pT(), double(pullvectorb1_paraproj_all_signed > 0) - double(pullvectorb1_paraproj_all_signed < 0));

     h_bJet1_pT->fill(bJets[0].pT());
     h_bJet2_pT->fill(bJets[1].pT());

     double leading_w_b_pull_all = CalculatePullAngle(wJets[0], bJets[0], 0);
     h_leading_w_b_pull_all->fill(leading_w_b_pull_all / Rivet::PI);


    }//end of per-event analysis 

    Vector3 CalculatePull(Jet& jet, bool &isCharged) { 
      Vector3 pull(0.0, 0.0, 0.0);
      double PT = jet.pT();
      Particles& constituents = jet.particles(); //getting particles in this jet
      Particles charged_constituents;
      if (isCharged) { //makes subset of jet constituents that are charged
        for (Particle p : constituents) {
          if (p.charge3() != 0)  charged_constituents += p;
        }
        constituents = charged_constituents;
      }
      // calculate axis
      FourMomentum axis;
      for (Particle p : constituents)  axis += p.momentum(); // Assigns axis as 4-momentum of individual particles
      Vector3 J(axis.rap(), axis.phi(MINUSPI_PLUSPI), 0.0); // 'MINUSPI_PLUSPI' --> want phi between -pi and pi, not 0 and 2pi. Calculates J 3-vector which defines jet axis
      // calculate pull
      for (Particle p : constituents) {
        Vector3 ri = Vector3(p.rap(), p.phi(MINUSPI_PLUSPI), 0.0) - J; //position of jet constituents relative to jet axis
        while (ri.y() >  Rivet::PI) ri.setY(ri.y() - Rivet::TWOPI); //making sure rapidity is between -pi and pi
        while (ri.y() < -Rivet::PI) ri.setY(ri.y() + Rivet::TWOPI); // 
        pull.setX(pull.x() + (ri.mod() * ri.x() * p.pT()) / PT); //calculate pull vector x component summing over all particles in jet
        pull.setY(pull.y() + (ri.mod() * ri.y() * p.pT()) / PT); //calculate pull vector y component summing over all particles in jet 
      }
      return pull;
    }

    double CalculatePullMagnitude(Jet& jet, bool isCharged) { 
      Vector3 pull_vector = CalculatePull(jet, isCharged);
      return pull_vector.mod(); 
    }

    double CalculatePullPerpProjection(Jet& jet1, Jet& jet2, bool isCharged) {
      Vector3 pull_vector = CalculatePull(jet1, isCharged);
      double drap = jet2.rap() - jet1.rap(); //calculating components of jet connection vector
      double dphi = jet2.phi(MINUSPI_PLUSPI) - jet1.phi(MINUSPI_PLUSPI);
      Vector3 jetconnection_vector(drap, dphi, 0.0); 
      double modC = jetconnection_vector.mod();
      Vector3 cross_vector = cross(pull_vector, jetconnection_vector);
      return (cross_vector.mod() / modC);

    } 

    double CalculatePullPerpProjectionSigned(Jet& jet1, Jet& jet2, bool isCharged) {
      Vector3 pull_vector = CalculatePull(jet1, isCharged);
      double drap = jet2.rap() - jet1.rap(); //calculating components of jet connection vector
      double dphi = jet2.phi(MINUSPI_PLUSPI) - jet1.phi(MINUSPI_PLUSPI);
      Vector3 jetconnection_vector(drap, dphi, 0.0);
      Vector3 planarnormal_vector(0.0, 0.0, 1.0);
      Vector3 axis_para = (jetconnection_vector/jetconnection_vector.mod());
      Vector3 axis_perp = cross(axis_para, planarnormal_vector);
      return (dot(pull_vector, axis_perp));

    }

    double CalculatePullParaProjection(Jet& jet1, Jet& jet2, bool isCharged) {
      Vector3 pull_vector = CalculatePull(jet1, isCharged);
      double drap = jet2.rap() - jet1.rap(); //calculating components of jet connection vector
      double dphi = jet2.phi(MINUSPI_PLUSPI) - jet1.phi(MINUSPI_PLUSPI);
      Vector3 jetconnection_vector(drap, dphi, 0.0); 
      double modC = jetconnection_vector.mod();
      return (dot(pull_vector, jetconnection_vector) / modC);

    }

    double CalculatePullParaProjectionSigned(Jet& jet1, Jet& jet2, bool isCharged) {
      Vector3 pull_vector = CalculatePull(jet1, isCharged);
      double drap = jet2.rap() - jet1.rap(); //calculating components of jet connection vector
      double dphi = jet2.phi(MINUSPI_PLUSPI) - jet1.phi(MINUSPI_PLUSPI);
      Vector3 jetconnection_vector(drap, dphi, 0.0);
      Vector3 axis_para = (jetconnection_vector/jetconnection_vector.mod());
      return (dot(pull_vector, axis_para));

    }

    double CalculatePullAngle(Jet& jet1, Jet& axisjet, bool isCharged) { //inputs: jet1 is jet whose pull vector are calculating, axisjet is other jet
      Vector3 pull_vector = CalculatePull(jet1, isCharged);
      pull_vector = Vector3(1000.*pull_vector.x(), 1000.*pull_vector.y(), 0.); //why 1000 x each component?
      double drap = axisjet.rap() - jet1.rap(); //calculating components of jet connection vector
      double dphi = axisjet.phi(MINUSPI_PLUSPI) - jet1.phi(MINUSPI_PLUSPI);
      Vector3 j2_vector(drap, dphi, 0.0); 
      return mapAngleMPiToPi(deltaPhi(pull_vector, j2_vector)); 
    }

    double CalculatePullAngle_ModHalfPI(Jet& jet1, Jet& axisjet, bool isCharged) { //inputs: jet1 is jet whose pull vector are calculating, axisjet is other jet
      Vector3 pull_vector = CalculatePull(jet1, isCharged);
      pull_vector = Vector3(1000.*pull_vector.x(), 1000.*pull_vector.y(), 0.); //why 1000 x each component?
      double drap = axisjet.rap() - jet1.rap(); //calculating components of jet connection vector
      double dphi = axisjet.phi(MINUSPI_PLUSPI) - jet1.phi(MINUSPI_PLUSPI);
      Vector3 j2_vector(drap, dphi, 0.0); 
      double pull_angle = fabs(mapAngleMPiToPi(deltaPhi(pull_vector, j2_vector))); 
      if ((pull_angle / Rivet::PI) < 0.5) {
        return pull_angle;
      }
      else {
        return (Rivet::PI - pull_angle);
      }
    }

   double CalculateInvariantMass(Jet& jet1, Jet& jet2) { //takes jets as non-composite objects, then add them together to get invariant mass
     FourMomentum totalmom = jet1.momentum();
     totalmom += jet2.momentum();
     
     return totalmom.mass();
   }

   double CalculateJetPurity(Jet& jet) {
     bool w_tagged;
     double Npure = 0;
     double Nimpure = 0;
     Particles& constituents = jet.particles();

     for (Particle p : constituents) {
       w_tagged = false;
       w_tagged |= p.hasAncestor(24, false);
       w_tagged |= p.hasAncestor(-24, false);
       if (w_tagged) Npure += 1;
       if (!w_tagged) Nimpure += 1;
      }

     return Npure/(Npure+Nimpure);
   }
 
    // Void finalize() is called after a run is finished. 
    // Here the analysis class should do whatever manipulations are necessary on the histograms.
    void finalize() {
      normalize(h_pull_12_all); 
      normalize(h_pull_21_all);
      normalize(h_pull_both_all);
      normalize(h_bpull_12_all);
      normalize(h_bpull_21_all);
      normalize(h_pull_angle_12_jetpTbins);
      normalize(h_pull_angle_21_jetpTbins);
      normalize(h_pull_12_all_lowpT);
      normalize(h_pull_12_all_highpT);
      normalize(h_pull_12_all_veryhighpT);
      normalize(h_wJet1_pT);
      normalize(h_wJet2_pT);
      normalize(h_bJet1_pT);
      normalize(h_bJet2_pT);
      normalize(h_pullmagnitudej1_all);
      normalize(h_bpull_both_all);
      normalize(h_pullvectorj1_perpproj_all);
      normalize(h_pullvectorj1_paraproj_all);
      normalize(h_pull_both_all_lowpT);
      normalize(h_pull_both_all_highpT);
      normalize(h_pull_both_all_veryhighpT);
      normalize(h_bpullmagnitudej1_all);
      normalize(h_leading_w_b_pull_all);
      normalize(h_pull_12_w);
      normalize(h_pull_12_lowpT_w);
      normalize(h_pull_12_highpT_w);

    }

    //@}


  private:

    Profile1DPtr p_pull_angle_paraproj_asymm_cor;
    Profile1DPtr p_pullmagnitude_paraproj_asymm_cor;
    Profile1DPtr p_paraproj_pull_angle_cor;
    Profile1DPtr p_paraproj_pull_angle_cor_signed;
    Profile1DPtr p_paraproj_pullmagnitude_cor;
    Profile1DPtr p_paraproj_asymm;
    Profile1DPtr p_paraproj_asymm_onebin;
    Profile1DPtr p_paraproj_asymm_M;
    Profile1DPtr p_paraproj_asymm_tail;
    Profile1DPtr p_perpproj_asymm;
    Profile1DPtr p_perpproj_asymm_onebin;
    Profile1DPtr p_perpproj_asymm_M;
    Profile1DPtr p_pull_angle_asymm; 
    Profile1DPtr p_bparaproj_asymm;
    Profile1DPtr p_bparaproj_asymm_onebin;
    Profile1DPtr p_bperpproj_asymm;
    Profile1DPtr p_bperpproj_asymm_onebin;
    Profile1DPtr p_paraproj_asymm_j2;
    Profile1DPtr p_paraproj_asymm_j2_onebin;
    Profile1DPtr p_paraproj_asymm_jetpTbins;
    Profile1DPtr p_paraproj_asymm_j2_jetpTbins;
    Profile1DPtr p_jetpT_pullmagnitude_cor;
    Profile1DPtr p_jetpT_paraproj_cor;
    Profile1DPtr p_paraproj_asymm_lowpT;
    Profile1DPtr p_paraproj_asymm_onebin_jetpTbins;
    Profile1DPtr p_paraproj_asymm_j2_onebin_jetpTbins;
    Profile1DPtr p_paraproj_asymm_onebin_jetpTbins_twobins;
    Profile1DPtr p_paraproj_asymm_j2_onebin_jetpTbins_twobins;
    Profile1DPtr p_paraproj_asymm_both;
    Profile1DPtr p_paraproj_asymm_highpT;
    Profile1DPtr p_paraproj_asymm_veryhighpT;
    Profile1DPtr p_bpull_angle_paraproj_asymm_cor;
    Profile1DPtr p_bparaproj_asymm_onebin_jetpTbins;
    Profile1DPtr p_paraproj_asymm_onebin_lowpT;
    Profile1DPtr p_paraproj_asymm_onebin_highpT;
    Profile1DPtr p_paraproj_asymm_onebin_veryhighpT;
    Profile1DPtr p_bparaproj_asymm_onebin_M;
    Profile1DPtr p_pull_angle_paraproj_asymm_cor_binone;
    Profile1DPtr p_paraproj_asymm_both_onebin;
    Profile1DPtr p_bpull_12_21_cor;
    Profile1DPtr p_wjetpurity;
    Profile1DPtr p_wjetpurity_pTbins;

    Histo1DPtr h_pull_12_all;
    Histo1DPtr h_pull_21_all;
    Histo1DPtr h_pull_both_all;
    Histo1DPtr h_bpull_12_all;
    Histo1DPtr h_bpull_21_all;
    Histo1DPtr h_pull_angle_12_jetpTbins;
    Histo1DPtr h_pull_angle_21_jetpTbins;
    Histo1DPtr h_pull_12_all_lowpT;
    Histo1DPtr h_pull_12_all_highpT;
    Histo1DPtr h_pull_12_all_veryhighpT;
    Histo1DPtr h_wJet1_pT;
    Histo1DPtr h_wJet2_pT;
    Histo1DPtr h_bJet1_pT;
    Histo1DPtr h_bJet2_pT;
    Histo1DPtr h_pullmagnitudej1_all;
    Histo1DPtr h_bpull_both_all;
    Histo1DPtr h_pullvectorj1_perpproj_all;
    Histo1DPtr h_pullvectorj1_paraproj_all;
    Histo1DPtr h_pull_both_all_lowpT;
    Histo1DPtr h_pull_both_all_highpT;
    Histo1DPtr h_pull_both_all_veryhighpT;
    Histo1DPtr h_bpullmagnitudej1_all;
    Histo1DPtr h_leading_w_b_pull_all;
    Histo1DPtr h_pull_12_w;
    Histo1DPtr h_pull_12_lowpT_w;
    Histo1DPtr h_pull_12_highpT_w;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ColourFlow);

}
