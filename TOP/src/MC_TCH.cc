// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/AnalysisLoader.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_TCH : public Analysis {
  public:

    /// Constructor
    MC_TCH() : Analysis("MC_TCH")
    {
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // generic final state
      FinalState fs(-5.0,5.0,0*GeV);
      
      // leptons
      ChargedLeptons lfs(FinalState((Cuts::abseta < 2.1 && Cuts::pT > 35*GeV && Cuts::abspid == PID::ELECTRON) || (Cuts::abseta < 2.4 && Cuts::pT > 26*GeV && Cuts::abspid == PID::MUON)));
      declare(lfs,"LFS");

      // neutrinos
      IdentifiedFinalState neutrinos(fs);
      neutrinos.acceptNeutrinos();
      declare(neutrinos, "Neutrinos");

      // jets
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(lfs);
      declare(FastJets(jet_fs, FastJets::ANTIKT, 0.4), "Jets");

      // MET
      declare(MissingMomentum(fs), "MissingET");

      // Book histograms
      //
      _h_njets = bookHisto1D("jet_mult", 11, -0.5, 10.5);
      _h_jet_HT = bookHisto1D("jet_HT", logspace(50, 100.0, 2000.0));
      //
      _h_lep_pT = bookHisto1D("lep_pT", logspace(20, 20.0, 200.0));
      _h_lep_eta = bookHisto1D("lep_eta", 25, -5.0, 5.0);
      //  
      _h_jet_1_pT = bookHisto1D("jet_1_pT", logspace(25, 20.0, 500.0));
      _h_jet_2_pT = bookHisto1D("jet_2_pT", logspace(25, 20.0, 400.0));
      _h_jet_1_eta = bookHisto1D("jet_1_eta", 25, -5.0, 5.0);
      _h_jet_2_eta = bookHisto1D("jet_2_eta", 25, -5.0, 5.0);
      //
      _h_bjet_pT = bookHisto1D("jetb_pT", logspace(25, 20.0, 400.0));
      _h_ljet_pT = bookHisto1D("jetl_pT", logspace(25, 20.0, 400.0)); 
      _h_bjet_eta = bookHisto1D("jetb_eta", 25, -5.0, 5.0);
      _h_ljet_eta = bookHisto1D("jetl_eta", 25, -5.0, 5.0);
      //
      _h_W_mass = bookHisto1D("W_mass", 25, 30, 180);
      _h_t_mass = bookHisto1D("t_mass", 25, 130, 430);
      //
      _h_t_pT = bookHisto1D("t_pT", logspace(25, 20.0, 500.0));
      //
      _h_jetb_jetl_mass = bookHisto1D("jetb_jetl_mass", 25, 0.0, 500.0);
      //
      _h_jetb_jetl_dR   = bookHisto1D("jetb_jetl_dR", 20, 0.0, 7.0);
      _h_jetb_jetl_deta = bookHisto1D("jetb_jetl_deta", 20, 0.0, 7.0);
      _h_jetb_jetl_dphi = bookHisto1D("jetb_jetl_dphi", 20, 0.0, M_PI);
      _h_jetb_W_dR        = bookHisto1D("jetb_W_dR", 20, 0.0, 7.0);
      _h_jetb_W_deta      = bookHisto1D("jetb_W_deta", 20, 0.0, 7.0);
      _h_jetb_W_dphi      = bookHisto1D("jetb_W_dphi", 20, 0.0, M_PI);
      _h_jetb_l_dR        = bookHisto1D("jetb_l_dR", 20, 0.0, 7.0);
      _h_jetb_l_deta      = bookHisto1D("jetb_l_deta", 20, 0.0, 7.0);
      _h_jetb_l_dphi      = bookHisto1D("jetb_l_dphi", 20, 0.0, M_PI);
      _h_jetb_l_mass      = bookHisto1D("jetb_l_mass", 25, 0.0, 500.0); 
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here
      const double weight = event.weight();

      // require one hard lepton and one matching neutrino, construct W boson
      const double wmass = 80.4;

      const ChargedLeptons& lfs = apply<ChargedLeptons>(event, "LFS");
      MSG_DEBUG("Charged lepton multiplicity = " << lfs.chargedLeptons().size());
      if (lfs.chargedLeptons().empty()) vetoEvent;
      if (lfs.chargedLeptons().size() != 1) vetoEvent;     
      Particle lepton = lfs.chargedLeptons()[0];

      const Particles neutrinos = apply<IdentifiedFinalState>(event, "Neutrinos").particlesByPt();
      Particle neutrino;
      FourMomentum W;
      for (const Particle& nu : neutrinos) {
        if (nu.abspid() == lepton.abspid()+1 and nu.pid()*lepton.pid()<0) {
          W = lepton.momentum() + nu.momentum();
          _h_W_mass->fill(W.mass(), weight);
        }
        else
          continue;
        if (inRange(W.mass(), wmass-5., wmass+5.)) {
          neutrino = nu;
          break;
        }
      }
      if (neutrino.pid() == 0)
        vetoEvent;

      // fill lepton histograms
      _h_lep_pT->fill(lepton.pT()/GeV, weight);
      _h_lep_eta->fill(lepton.eta(), weight);

      // missing ET
      const MissingMomentum& met = apply<MissingMomentum>(event, "MissingET");
      MSG_DEBUG("Vector ET = " << met.vectorEt().mod() << " GeV");

      // jets projection: consider all jets with pT > 40 GeV and abs(eta) < 4.7
      Cut cut_jets =  Cuts::abseta < 4.7 && Cuts::pT > 40*GeV;
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(cut_jets);
      _h_njets->fill(jets.size(), weight);
      double ht = 0.0;
      foreach (const Jet& j, jets) { ht += j.pT(); }
      _h_jet_HT->fill(ht/GeV, weight);
      Jets bjets, ljets;
      // only events with jets that do not overlap with hard lepton selected
      foreach (const Jet& jet, jets) {
        bool isolated = true;
        if (deltaR(jet.momentum(), lepton.momentum()) < 0.4) {
          isolated = false;
          break;
        }
        if (!isolated) {
          MSG_DEBUG("Jet failed lepton isolation cut");
          break;
        }
      // find b-tagged and light jet
      if (jet.bTagged(Cuts::abseta < 2.4)) bjets.push_back(jet);
      else ljets.push_back(jet);
      }
      // only events with exact 1 b-tagged jet and 1 light jet selected
      if (bjets.size() != 1) {
        MSG_DEBUG("Event failed b jet multiplicity cut");
        vetoEvent;
      }
      if (ljets.size() != 1) {
        MSG_DEBUG("Event failed light jet multiplicity cut");
        vetoEvent;
      }

      // reconstruct top quark
      FourMomentum t = W + bjets[0].momentum();
      _h_t_mass->fill(t.mass(), weight);
      _h_t_pT->fill(t.pT()/GeV, weight);
           
      // fill jet histograms
      _h_jet_1_pT->fill(jets[0].pT()/GeV, weight);
      _h_jet_2_pT->fill(jets[1].pT()/GeV, weight);
      _h_jet_1_eta->fill(jets[0].eta(), weight);
      _h_jet_2_eta->fill(jets[1].eta(), weight);
      _h_bjet_pT->fill(bjets[0].pT()/GeV, weight);
      _h_ljet_pT->fill(ljets[0].pT()/GeV, weight);
      _h_bjet_eta->fill(bjets[0].eta(), weight);
      _h_ljet_eta->fill(ljets[0].eta(), weight);

      // dijet mass histogram
      _h_jetb_jetl_mass->fill(FourMomentum(bjets[0].momentum()+ljets[0].momentum()).mass(), weight);

      // more histograms 
      _h_jetb_jetl_dR->fill(deltaR(bjets[0].momentum(), ljets[0].momentum()),weight);
      _h_jetb_jetl_deta->fill(fabs(bjets[0].eta()-ljets[0].eta()),weight);
      _h_jetb_jetl_dphi->fill(deltaPhi(bjets[0].momentum(),ljets[0].momentum()),weight);
      _h_jetb_W_dR->fill(deltaR(bjets[0].momentum(), W),weight);
      _h_jetb_W_deta->fill(fabs(bjets[0].eta()-W.eta()),weight);
      _h_jetb_W_dphi->fill(deltaPhi(bjets[0].momentum(),W),weight);
      FourMomentum l = lepton.momentum();
      _h_jetb_l_dR->fill(deltaR(bjets[0].momentum(), l),weight);
      _h_jetb_l_deta->fill(fabs(bjets[0].eta()-l.eta()),weight);
      _h_jetb_l_dphi->fill(deltaPhi(bjets[0].momentum(),l),weight);
      _h_jetb_l_mass->fill(FourMomentum(bjets[0].momentum()+l).mass(), weight);

    }


    /// Normalise histograms etc., after the run

    void finalize() {

      normalize(_h_njets);
      normalize(_h_jet_HT);
      normalize(_h_lep_pT);
      normalize(_h_lep_eta);
      normalize(_h_jet_1_pT);
      normalize(_h_jet_2_pT);
      normalize(_h_jet_1_eta);
      normalize(_h_jet_2_eta);
      normalize(_h_bjet_pT);
      normalize(_h_ljet_pT);
      normalize(_h_bjet_eta);
      normalize(_h_ljet_eta);
      normalize(_h_W_mass);
      normalize(_h_t_mass);
      normalize(_h_t_pT);
      normalize(_h_jetb_jetl_mass);
      normalize(_h_jetb_jetl_dR);
      normalize(_h_jetb_jetl_deta);
      normalize(_h_jetb_jetl_dphi);
      normalize(_h_jetb_W_dR);
      normalize(_h_jetb_W_deta);
      normalize(_h_jetb_W_dphi);
      normalize(_h_jetb_l_dR);
      normalize(_h_jetb_l_deta);
      normalize(_h_jetb_l_dphi);
      normalize(_h_jetb_l_mass);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_njets;
    Histo1DPtr _h_jet_HT;
    Histo1DPtr _h_lep_pT, _h_lep_eta;
    Histo1DPtr _h_jet_1_pT, _h_jet_2_pT;
    Histo1DPtr _h_jet_1_eta, _h_jet_2_eta;
    Histo1DPtr _h_bjet_pT, _h_ljet_pT;
    Histo1DPtr _h_bjet_eta, _h_ljet_eta;
    Histo1DPtr _h_W_mass;
    Histo1DPtr _h_t_mass;
    Histo1DPtr _h_t_pT;
    Histo1DPtr _h_jetb_jetl_mass;
    Histo1DPtr _h_jetb_jetl_dR, _h_jetb_jetl_deta, _h_jetb_jetl_dphi;
    Histo1DPtr _h_jetb_W_dR, _h_jetb_W_deta, _h_jetb_W_dphi;
    Histo1DPtr _h_jetb_l_dR, _h_jetb_l_deta, _h_jetb_l_dphi,_h_jetb_l_mass;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TCH);


}
