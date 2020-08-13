// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Using colour information to identify jets in semileptonic top decay at 13 TeV
  class MC_ColourReconstruction_tttt : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_ColourReconstruction_tttt);


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
      book(h_numberwjets, "numberwjets", 8, -0.5, 7.5);
      book(h_numberbjets, "numberbjets", 5, 2.5, 8.5);
      book(h_numberlightjets, "numberlightjets", 6, 5.5, 11.5);
      book(h_truewjetsinvariantmass, "truewjetsinvariantmass", 150, 0, 300);
      book(h_numberojets, "numberojets", 8, -0.5, 7.5);

      //book profiles
      book(p_PT_s1, "PT_s1", 20, 0, 1);
      book(p_TP_s1, "TP_s1", 20, 0, 1);
      book(p_FP_s1, "FP_s1", 20, 0, 1);
      book(p_PT_s2, "PT_s2", 20, 0, 1);
      book(p_TP_s2, "TP_s2", 20, 0, 1);
      book(p_FP_s2, "FP_s2", 20, 0, 1);

    } 

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /**************
       *    JETS    *
       **************/
      const Jets& allJets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30.0*GeV && Cuts::absrap < 2.5); //jet pT,rap cut - high jet pT cut makes pull angle smaller
      const vector<DressedLepton>& all_elecs = apply<DressedLeptons>(event, "ELECS").dressedLeptons();
      const vector<DressedLepton>& all_muons = apply<DressedLeptons>(event, "MUONS").dressedLeptons();
      Jets goodJets;
      for (const Jet j : allJets) {
        bool keep = true;
        for (const DressedLepton el : all_elecs)  keep &= deltaR(j, el) >= 0.2; //ensure no overlap between jets and electrons
        if (keep)  goodJets += j;
      }
      if ( goodJets.size() < 9 )  vetoEvent; //jet multiplicity minimum 
      
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
        } //vetoMuons looser cut than muons, muons subsection of vetomuons
      }

      vector<DressedLepton> elecs, vetoElecs;
      for (const DressedLepton el : all_elecs) {
        bool keep = true;
        for (const Jet j : goodJets)  keep &= deltaR(j, el) >= 0.4;
        if (keep && el.pt() > 15*GeV) {
          vetoElecs.push_back(el);;
          if (el.pt() > 20*GeV)  elecs.push_back(el); //electron pT cut 
        }
      }

      if (muons.empty() && elecs.empty())  vetoEvent; //mono-leptonic 

      bool muCandidate = !( muons.size() < 1 || vetoMuons.size() > 1 || vetoElecs.size() > 0 ); //1 muon, 0 electrons
      bool elCandidate = !( elecs.size() < 1 || vetoElecs.size() > 1 || vetoMuons.size() > 0 ); //1 electron, 0 muons 

      if (!elCandidate && !muCandidate)  vetoEvent; //mono-leptonic 

      /******************************
       *    ELECTRON-MUON OVERLAP   *
       ******************************/
      for (const DressedLepton electron : elecs) {
        for (const DressedLepton muon : muons) {
          double d_theta = fabs(muon.theta() - electron.theta());
          double d_phi = deltaPhi(muon.phi(), electron.phi());
          if (d_theta < 0.005 && d_phi < 0.005)  vetoEvent; //
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
      if (met <= 50*GeV)  vetoEvent; //supress residual QCD background, mono-leptonic events 

      /****************
       *    B-JETS    *
       ****************/
      // Sort the jets into b-jets and light jets. We expect one hard b-jet from
      // each top decay, so our 4 hardest jets should include two b-jets. The
      // Jet::bTagged() method is equivalent to perfect experimental
      // b-tagging, in a generator-independent way.
      Jets bJets, lJets;
      double Ht = 0;
      for(Jet j : goodJets) {
        bool b_tagged = false;
        Particles bTags = j.bTags();
        for ( Particle b : bTags ) {
          b_tagged |= b.pT() > 5*GeV; //5 GeV pT cut on b particles
        }
        if (b_tagged)  bJets += j;
        if (!b_tagged) lJets += j;
        Ht += j.pT();
      }

      if (Ht < 500*GeV) vetoEvent;

      if ( bJets.size() < 3 || lJets.size() < 6 )  vetoEvent; //b jets multiplicity cut
      h_numberlightjets->fill(lJets.size());
      h_numberbjets->fill(bJets.size());
      

      //Find W jets
      Jets wJets, oJets;
      for(Jet j : lJets) {
        bool w_tagged = false;
        Particles& constituents = j.particles();
        for (Particle p : constituents) {
          w_tagged |= p.hasAncestor(24, false);
          w_tagged |= p.hasAncestor(-24, false);
        }
        if (w_tagged) {
          double purity = CalculateWJetPurity(j);
          if (purity > 0.1) wJets += j;
        }
        else {
          oJets += j;
        }
      }
      
      h_numberwjets->fill(wJets.size());
      h_numberojets->fill(oJets.size());

      if (wJets.size() > 1) {
        for (size_t i = 0; i < wJets.size()-1; ++i) {
          for (size_t j = i+1; j < wJets.size(); ++j) {
            double wjetsinvariantmass = CalculateInvariantMass(wJets[i], wJets[j]);
            h_truewjetsinvariantmass->fill(wjetsinvariantmass);
          }
        }
      }

      double Ncuts = 60;
      //pT cuts

      //PULL ANGLE SELECTION 1 (loose):
      //Reconstructing W from light jets with colour connection information via pull angle 
      vector<double> s1_Nsignal(Ncuts);
      vector<double> s1_Nbackground(Ncuts);
      for (size_t i = 0; i < lJets.size()-1; ++i) { 
        for (size_t j = i + 1; j < lJets.size(); ++j) {
          if (lJets[i].pT() > 80*GeV && lJets[j].pT() > 80*GeV) {
            double pull_12 = CalculatePullAngle(lJets[i], lJets[j], 0);
            double pull_21 = CalculatePullAngle(lJets[j], lJets[i], 0);
            for (double k = 0; k < Ncuts; ++k) {
              double cut = (k+1)/Ncuts;
              if ((pull_12/Rivet::PI) < cut || (pull_21/Rivet::PI) < cut) {
                double wpurity_j1 = CalculateWJetPurity(lJets[i]);
                double wpurity_j2 = CalculateWJetPurity(lJets[j]);
                if (wpurity_j1 > 0.1 && wpurity_j2 > 0.1) {
                  s1_Nsignal[k] += 1;
                }
                if (wpurity_j1 < 0.1 || wpurity_j2 < 0.1) {
                  s1_Nbackground[k] += 1;
                }
              }
            }
          }
        }
      }

      vector<double> s1_signal_efficiency(Ncuts);
      for (double i = 0; i < Ncuts; ++i) {
        double cut = (i+1)/Ncuts;
        if (s1_Nsignal[i] > 0 || s1_Nbackground[i] > 0) {
          s1_signal_efficiency[i] = s1_Nsignal[i]/(s1_Nsignal[i]+s1_Nbackground[i]);
          p_PT_s1->fill(cut, s1_signal_efficiency[i]);
        }
      }

      vector<double> s1_TP(Ncuts);
      vector<double> s1_TN(Ncuts);
      if (wJets.size() > 1) {
        for (size_t i = 0; i < wJets.size()-1; ++i) { 
          for (size_t j = i + 1; j < wJets.size(); ++j) {
            if (wJets[i].pT() > 80*GeV && wJets[j].pT() > 80*GeV) {
              bool same = FromSameW(wJets[i], wJets[j]);
              if (same) {
                double pull_12 = CalculatePullAngle(wJets[i], wJets[j], 0);
                double pull_21 = CalculatePullAngle(wJets[j], wJets[i], 0);
                for (double k = 0; k < Ncuts; ++k) {
                  double cut = (k+1)/Ncuts;
                  if ((pull_12/Rivet::PI) < cut || (pull_21/Rivet::PI) < cut) {
                    s1_TP[k] += 1;
                  }
                  if ((pull_12/Rivet::PI) > cut && (pull_21/Rivet::PI) > cut) {
                    s1_TN[k] += 1;
                  }
                }
              }
            }
          }
        }
      }

      vector<double> s1_TP_efficiency(Ncuts);
      for (double i = 0; i < Ncuts; ++i) {
        double cut = (i+1)/Ncuts;
        if (s1_TP[i] > 0 || s1_TN[i] > 0) {
          s1_TP_efficiency[i] = s1_TP[i]/(s1_TP[i]+s1_TN[i]);
          p_TP_s1->fill(cut, s1_TP_efficiency[i]);
        }
      }

      vector<double> s1_FP(Ncuts);
      vector<double> s1_FN(Ncuts);
      if (oJets.size() > 1) {
        for (size_t i = 0; i < oJets.size()-1; ++i) { 
          for (size_t j = i + 1; j < oJets.size(); ++j) {
            if (oJets[i].pT() > 80*GeV && oJets[j].pT() > 80*GeV) {
              double pull_12 = CalculatePullAngle(oJets[i], oJets[j], 0);
              double pull_21 = CalculatePullAngle(oJets[j], oJets[i], 0);
              for (double k = 0; k < Ncuts; ++k) {
                double cut = (k+1)/Ncuts;
                if ((pull_12/Rivet::PI) < cut || (pull_21/Rivet::PI) < cut) {
                  s1_FP[k] += 1;
                }
                if ((pull_12/Rivet::PI) > cut && (pull_21/Rivet::PI) > cut) {
                  s1_FN[k] += 1;
                }
              }
            }
          }
        }
      }

      vector<double> s1_FP_efficiency(Ncuts);
      for (double i = 0; i < Ncuts; ++i) {
        double cut = (i+1)/Ncuts;
        if (s1_FP[i] > 0 || s1_FN[i] > 0) {
          s1_FP_efficiency[i] = s1_FP[i]/(s1_FP[i]+s1_FN[i]);
          p_FP_s1->fill(cut, s1_FP_efficiency[i]);
        }
      }

      //PULL ANGLE SELECTION 2 (tight):
      //Reconstructing W from light jets with colour connection information via pull angle 
      vector<double> s2_Nsignal(Ncuts);
      vector<double> s2_Nbackground(Ncuts);
      for (size_t i = 0; i < lJets.size()-1; ++i) { 
        for (size_t j = i + 1; j < lJets.size(); ++j) {
          if (lJets[i].pT() > 80*GeV && lJets[j].pT() > 80*GeV) {
            double pull_12 = CalculatePullAngle(lJets[i], lJets[j], 0);
            double pull_21 = CalculatePullAngle(lJets[j], lJets[i], 0);
            for (double k = 0; k < Ncuts; ++k) {
              double cut = (k+1)/Ncuts;
              if ((pull_12/Rivet::PI) < cut && (pull_21/Rivet::PI) < cut) {
                double wpurity_j1 = CalculateWJetPurity(lJets[i]);
                double wpurity_j2 = CalculateWJetPurity(lJets[j]);
                if (wpurity_j1 > 0.1 && wpurity_j2 > 0.1) {
                  s2_Nsignal[k] += 1;
                }
                if (wpurity_j1 < 0.1 || wpurity_j2 < 0.1) {
                  s2_Nbackground[k] += 1;
                }
              }
            }
          }
        }
      }

      vector<double> s2_signal_efficiency(Ncuts);
      for (double i = 0; i < Ncuts; ++i) {
        double cut = (i+1)/Ncuts;
        if (s2_Nsignal[i] > 0 || s2_Nbackground[i] > 0) {
          s2_signal_efficiency[i] = s2_Nsignal[i]/(s2_Nsignal[i]+s2_Nbackground[i]);
          p_PT_s2->fill(cut, s2_signal_efficiency[i]);
        }
      }

      vector<double> s2_TP(Ncuts);
      vector<double> s2_TN(Ncuts);
      if (wJets.size() > 1) {
        for (size_t i = 0; i < wJets.size()-1; ++i) { 
          for (size_t j = i + 1; j < wJets.size(); ++j) {
            if (wJets[i].pT() > 80*GeV && wJets[j].pT() > 80*GeV) {
              bool same = FromSameW(wJets[i], wJets[j]);
              if (same) {
                double pull_12 = CalculatePullAngle(wJets[i], wJets[j], 0);
                double pull_21 = CalculatePullAngle(wJets[j], wJets[i], 0);
                for (double k = 0; k < Ncuts; ++k) {
                  double cut = (k+1)/Ncuts;
                  if ((pull_12/Rivet::PI) < cut || (pull_21/Rivet::PI) < cut) {
                    s2_TP[k] += 1;
                  }
                  if ((pull_12/Rivet::PI) > cut && (pull_21/Rivet::PI) > cut) {
                    s2_TN[k] += 1;
                  }
                }
              }
            }
          }
        }
      }

      vector<double> s2_TP_efficiency(Ncuts);
      for (double i = 0; i < Ncuts; ++i) {
        double cut = (i+1)/Ncuts;
        if (s2_TP[i] > 0 || s2_TN[i] > 0) {
          s2_TP_efficiency[i] = s2_TP[i]/(s2_TP[i]+s2_TN[i]);
          p_TP_s2->fill(cut, s2_TP_efficiency[i]);
        }
      }

      vector<double> s2_FP(Ncuts);
      vector<double> s2_FN(Ncuts);
      if (oJets.size() > 1) {
        for (size_t i = 0; i < oJets.size()-1; ++i) { 
          for (size_t j = i + 1; j < oJets.size(); ++j) {
            if (oJets[i].pT() > 80*GeV && oJets[j].pT() > 80*GeV) {
              double pull_12 = CalculatePullAngle(oJets[i], oJets[j], 0);
              double pull_21 = CalculatePullAngle(oJets[j], oJets[i], 0);
              for (double k = 0; k < Ncuts; ++k) {
                double cut = (k+1)/Ncuts;
                if ((pull_12/Rivet::PI) < cut && (pull_21/Rivet::PI) < cut) {
                  s2_FP[k] += 1;
                }
                if ((pull_12/Rivet::PI) > cut || (pull_21/Rivet::PI) > cut) {
                  s2_FN[k] += 1;
                }
              }
            }
          }
        }
      }

      vector<double> s2_FP_efficiency(Ncuts);
      for (double i = 0; i < Ncuts; ++i) {
        double cut = (i+1)/Ncuts;
        if (s2_FP[i] > 0 || s2_FN[i] > 0) {
          s2_FP_efficiency[i] = s2_FP[i]/(s2_FP[i]+s2_FN[i]);
          p_FP_s2->fill(cut, s2_FP_efficiency[i]);
        }
      }
  


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

    double CalculatePullAngle(Jet& jet1, Jet& axisjet, bool isCharged) { //inputs: jet1 is jet whose pull vector are calculating, axisjet is other jet
      Vector3 pull_vector = CalculatePull(jet1, isCharged);
      pull_vector = Vector3(1000.*pull_vector.x(), 1000.*pull_vector.y(), 0.); //why 1000 x each component?
      double drap = axisjet.rap() - jet1.rap(); //calculating components of jet connection vector
      double dphi = axisjet.phi(MINUSPI_PLUSPI) - jet1.phi(MINUSPI_PLUSPI);
      Vector3 j2_vector(drap, dphi, 0.0); 
      return mapAngleMPiToPi(deltaPhi(pull_vector, j2_vector)); 
    }

   double CalculateInvariantMass(Jet& jet1, Jet& jet2) { //takes jets as non-composite objects, then add them together to get invariant mass
     FourMomentum totalmom = jet1.momentum();
     totalmom += jet2.momentum();
     
     return totalmom.mass();
   }

   double CalculateWJetPurity(Jet& jet) {
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

   bool FromSameW(Jet& jet1, Jet& jet2) {
     bool same = false;

     Particles W_j1;
     Particles& constituents_j1 = jet1.particles();
     for (Particle p_j1 : constituents_j1) {
       Particles ancestors_j1 = p_j1.ancestors(Cuts::OPEN, false);
       for (Particle a_j1 : ancestors_j1) {
         if (a_j1.abspid() == 24) W_j1 += a_j1;
       }
     }
     Particles W_j2;
     Particles& constituents_j2 = jet2.particles();
     for (Particle p_j2 : constituents_j2) {
       Particles ancestors_j2 = p_j2.ancestors(Cuts::OPEN, false);
       for (Particle a_j2 : ancestors_j2) {
         if (a_j2.abspid() == 24) W_j2 += a_j2;
       }
     }

     if (W_j1[0].isSame(W_j2[0])) same = true;
     return same;
    }
   
    // Void finalize() is called after a run is finished. 
    // Here the analysis class should do whatever manipulations are necessary on the histograms.
    void finalize() {
      normalize(h_numberwjets);
      normalize(h_numberlightjets);
      normalize(h_numberbjets);
      normalize(h_truewjetsinvariantmass);
      normalize(h_numberojets);
    }

    //@}

  private:

    Histo1DPtr h_numberwjets;
    Histo1DPtr h_numberlightjets;
    Histo1DPtr h_numberbjets;
    Histo1DPtr h_truewjetsinvariantmass;
    Histo1DPtr h_numberojets;

    Profile1DPtr p_PT_s1;
    Profile1DPtr p_TP_s1;
    Profile1DPtr p_FP_s1;
    Profile1DPtr p_PT_s2;
    Profile1DPtr p_TP_s2;
    Profile1DPtr p_FP_s2;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ColourReconstruction_tttt);

}
