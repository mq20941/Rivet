#include "Rivet/Analysis.hh"
//#include "Rivet/AnalysisLoader.hh"
//#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {

  namespace { //< only visible in this compilation unit
  
    // @brief Pseudo top finder
    // 
    // Find top quark in the particle level.
    // The definition is based on the agreement at the LHC working group.
    class PseudoTop : public FinalState {
      public:
        /// @name Standard constructors and destructors.
        //@{

        /// The default constructor. May specify the minimum and maximum
        /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
        PseudoTop(double lepR = 0.1, double lepMinPt = 15, double lepMaxEta = 2.5,
            double jetR = 0.4, double jetMinPt = 30, double jetMaxEta = 2.4)
          : FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV),
          _lepR(lepR), _lepMinPt(lepMinPt), _lepMaxEta(lepMaxEta),
          _jetR(jetR), _jetMinPt(jetMinPt), _jetMaxEta(jetMaxEta)
        {
          setName("PseudoTop");
        }

        enum TTbarMode {CH_NONE=-1, CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON};
        enum DecayMode {CH_HADRON = 0, CH_MUON, CH_ELECTRON};

        TTbarMode mode() const {
          if (!_isValid) return CH_NONE;
          if (_mode1 == CH_HADRON && _mode2 == CH_HADRON) return CH_FULLHADRON;
          else if ( _mode1 != CH_HADRON && _mode2 != CH_HADRON) return CH_FULLLEPTON;
          else return CH_SEMILEPTON;
        }
        DecayMode mode1() const {return _mode1;}
        DecayMode mode2() const {return _mode2;}
        
        /// Clone on the heap.
        virtual unique_ptr<Projection> clone() const {
          return unique_ptr<Projection>(new PseudoTop(*this));
        }

        //@}

      public:
        Particle t1() const {return _t1;}
        Particle t2() const {return _t2;}
        Particle b1() const {return _b1;}
        Particle b2() const {return _b2;}
        Particle w1() const {return _w1;}
        Particle w2() const {return _w2;}
        ParticleVector wDecays1() const {return _wDecays1;}
        ParticleVector wDecays2() const {return _wDecays2;}
        Jets jets() const {return _jets;}
        Jets bjets() const {return _bjets;}
        Jets ljets() const {return _ljets;}

      protected:
        // Apply the projection to the event
        void project(const Event& e) override;
        void cleanup(std::map<double, std::pair<size_t, size_t> >& v, const bool doCrossCleanup=false) const;

      private:
        const double _lepR, _lepMinPt, _lepMaxEta;
        const double _jetR, _jetMinPt, _jetMaxEta;

        constexpr static double _tMass = 172.5;
        constexpr static double _wMass = 80.4;

      private:
        bool _isValid;
        DecayMode _mode1, _mode2;

        Particle _t1, _t2;
        Particle _b1, _b2;
        Particle _w1, _w2;
        ParticleVector _wDecays1, _wDecays2;
        Jets _jets, _bjets, _ljets;

    };
  
  }

  class CMS_LesHouches2015 : public Analysis {
    public:
      CMS_LesHouches2015() : Analysis("CMS_LesHouches2015") {
      }

      void init() {
        addProjection(PseudoTop(), "ttbar");

        // Booking of histograms
        _h_njets  = bookHisto1D("jet_mult", 11, -0.5, 10.5);
        _h_nbjets = bookHisto1D("bjet_mult", 6, -0.5, 5.5);
        _h_nljets = bookHisto1D("ljet_mult", 6, -0.5, 5.5);

        _h_jet_1_pT = bookHisto1D("jet_1_pT", logspace(50, 20., 1000.));
        _h_jet_2_pT = bookHisto1D("jet_2_pT", logspace(50, 20., 700.));
        _h_jet_3_pT = bookHisto1D("jet_3_pT", logspace(50, 20., 500.));
        _h_jet_4_pT = bookHisto1D("jet_4_pT", logspace(50, 20., 400.));
        _h_bjet_1_pT = bookHisto1D("jetb_1_pT", logspace(50, 20., 1000.));
        _h_bjet_2_pT = bookHisto1D("jetb_2_pT", logspace(50, 20., 700.));
        _h_bjet_3_pT = bookHisto1D("jetb_3_pT", logspace(50, 20., 500.));
        _h_bjet_4_pT = bookHisto1D("jetb_4_pT", logspace(50, 20., 400.));
        _h_ljet_1_pT = bookHisto1D("jetl_1_pT", logspace(50, 20., 1000.));
        _h_ljet_2_pT = bookHisto1D("jetl_2_pT", logspace(50, 20., 700.));
        _h_ljet_3_pT = bookHisto1D("jetl_3_pT", logspace(50, 20., 500.));
        _h_ljet_4_pT = bookHisto1D("jetl_4_pT", logspace(50, 20., 400.));

        _h_top_pT = bookHisto1D("top_pT", logspace(50, 10., 500));
        _h_top_mass = bookHisto1D("top_mass", 100, 0., 500.);
        _h_ttbar_pT = bookHisto1D("ttbar_pT", logspace(50, 10., 1000.));
        _h_ttbar_mass = bookHisto1D("ttbar_mass", 200, 0., 2000.);
        _h_ttjet_pT = bookHisto1D("ttjet_pT", logspace(50, 10., 1000.));
        _h_b1b2_mass = bookHisto1D("b1b2_mass", 100, 0, 1000);

        _h_lep_1_lep_2_dphi = bookHisto1D("lep_1_lep_2_dphi", 32, 0.0, 3.2);
        _h_jet_1_jet_2_dphi = bookHisto1D("jet_1_jet_2_dphi", 32, 0.0, 3.2);

        _h_nLepton_deta = bookHisto1D("n_lepton_deta", 2, -1., 1.);
        _h_nTtbar_dy = bookHisto1D("n_ttbar_dy", 2, -1., 1.);
        _h_nNoAddJet = bookHisto1D("n_noaddjet", 2, -0.5, 1.5);

        _h_eta_lplus = bookHisto1D("eta_lplus",50,-2.5,2.5);
        _h_eta_lminus = bookHisto1D("eta_lminus",50,-2.5,2.5);
        _h_y_tplus = bookHisto1D("eta_tplus",50,-5.,5.);
        _h_y_tminus = bookHisto1D("eta_tminus",50,-5.,5.);

        _h_ttbar_mass_deta_pos = bookHisto1D("ttbar_mass_deta_pos", 20, 0., 2000.);
        _h_ttbar_mass_deta_neg = bookHisto1D("ttbar_mass_deta_neg", 20, 0., 2000.);
        _h_ttbar_mass_dy_pos = bookHisto1D("ttbar_mass_dy_pos", 20, 0., 2000.);
        _h_ttbar_mass_dy_neg = bookHisto1D("ttbar_mass_dy_neg", 20, 0., 2000.);
        _h_ttbar_mass_noaddjet = bookHisto1D("ttbar_mass_noaddjet", 20, 0., 2000.);
        _h_ttbar_mass_addjet = bookHisto1D("ttbar_mass_addjet", 20, 0., 2000.);

        _s_lepton_asym_ttbar_mass = bookScatter2D("lepton_asym_ttbar_mass",20,0.,2000.);
        _s_ttbar_asym_ttbar_mass = bookScatter2D("ttbar_asym_ttbar_mass",20,0.,2000.);
        _s_gap_fraction_ttbar_mass = bookScatter2D("gap_fraction_ttbar_mass",20,0.,2000.);
      };

      void analyze(const Event& event) {
        const double weight = event.weight();

        // Get the parton level ttbar candidate
        const PseudoTop& ttbar = applyProjection<PseudoTop>(event, "ttbar");
        if (ttbar.mode() != PseudoTop::CH_FULLLEPTON) {
          MSG_DEBUG("Event fail channel topology cuts");
          vetoEvent;
        }

        _h_nNoAddJet->fill(1., weight);
        _h_ttbar_mass_addjet->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
        if (ttbar.ljets().size() < 1) {
          _h_nNoAddJet->fill(0., weight);
          _h_ttbar_mass_noaddjet->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
        }

        _h_njets->fill(ttbar.jets().size(), weight);
        _h_nbjets->fill(ttbar.bjets().size(), weight);
        _h_nljets->fill(ttbar.ljets().size(), weight);

        if (ttbar.jets().size() > 0) {
          _h_jet_1_pT->fill(ttbar.jets()[0].momentum().pT(), weight);
          if (ttbar.jets().size() > 1) {
            _h_jet_2_pT->fill(ttbar.jets()[1].momentum().pT(), weight);
            if (ttbar.jets().size() > 2) {
              _h_jet_3_pT->fill(ttbar.jets()[2].momentum().pT(), weight);
              if (ttbar.jets().size() > 3) {
                _h_jet_4_pT->fill(ttbar.jets()[3].momentum().pT(), weight);
              }
            }
          }
        }
        if (ttbar.bjets().size() > 0) {
          _h_bjet_1_pT->fill(ttbar.bjets()[0].momentum().pT(), weight);
          if (ttbar.bjets().size() > 1) {
            _h_bjet_2_pT->fill(ttbar.bjets()[1].momentum().pT(), weight);
            if (ttbar.bjets().size() > 2) {
              _h_bjet_3_pT->fill(ttbar.bjets()[2].momentum().pT(), weight);
              if (ttbar.bjets().size() > 3) {
                _h_bjet_4_pT->fill(ttbar.bjets()[3].momentum().pT(), weight);
              }
            }
          }
        }
        if (ttbar.ljets().size() > 0) {
          _h_ljet_1_pT->fill(ttbar.ljets()[0].momentum().pT(), weight);
          if (ttbar.ljets().size() > 1) {
            _h_ljet_2_pT->fill(ttbar.ljets()[1].momentum().pT(), weight);
            if (ttbar.ljets().size() > 2) {
              _h_ljet_3_pT->fill(ttbar.ljets()[2].momentum().pT(), weight);
              if (ttbar.ljets().size() > 3) {
                _h_ljet_4_pT->fill(ttbar.ljets()[3].momentum().pT(), weight);
              }
            }
          }
        }

        _h_top_pT->fill(ttbar.t1().momentum().pT(), weight);
        _h_top_pT->fill(ttbar.t2().momentum().pT(), weight);
        _h_top_mass->fill(ttbar.t1().momentum().mass(), weight);
        _h_top_mass->fill(ttbar.t2().momentum().mass(), weight);
        _h_ttbar_pT->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).pT(), weight);
        _h_ttbar_mass->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
        if (ttbar.ljets().size() > 0) 
          _h_ttjet_pT->fill((ttbar.t1().momentum()+ttbar.t2().momentum()+ttbar.ljets()[0].momentum()).pT(), weight);
        _h_b1b2_mass->fill((ttbar.b1().momentum()+ttbar.b2().momentum()).mass(), weight);

        _h_lep_1_lep_2_dphi->fill(deltaPhi(ttbar.wDecays1()[0].momentum(),ttbar.wDecays2()[0].momentum()), weight);
        if (ttbar.jets().size() > 1)
          _h_jet_1_jet_2_dphi->fill(deltaPhi(ttbar.jets()[0].momentum(),ttbar.jets()[1].momentum()), weight);

        double lepton_deta = 0.;
        if (ttbar.wDecays1()[0].charge() > 0) { 
          lepton_deta = ttbar.wDecays1()[0].momentum().abseta() - ttbar.wDecays2()[0].momentum().abseta();
          _h_eta_lplus->fill(ttbar.wDecays1()[0].momentum().eta());
          _h_eta_lminus->fill(ttbar.wDecays2()[0].momentum().eta());
        }
        else if (ttbar.wDecays2()[0].charge() > 0) { 
          lepton_deta = ttbar.wDecays2()[0].momentum().abseta() - ttbar.wDecays1()[0].momentum().abseta();
          _h_eta_lplus->fill(ttbar.wDecays2()[0].momentum().eta());
          _h_eta_lminus->fill(ttbar.wDecays1()[0].momentum().eta());
        }
        if (lepton_deta > 1e-6) {
          _h_nLepton_deta->fill(0.5, weight);
          _h_ttbar_mass_deta_pos->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
        }
        if (lepton_deta < -1e-6) {
          _h_nLepton_deta->fill(-0.5, weight);
          _h_ttbar_mass_deta_neg->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
        }

        double ttbar_dy = 0.;
        if (ttbar.t1().charge() > 0) {
          ttbar_dy = ttbar.t1().momentum().absrapidity() - ttbar.t2().momentum().absrapidity();
          _h_y_tplus->fill(ttbar.t1().momentum().rapidity());
          _h_y_tminus->fill(ttbar.t2().momentum().rapidity());
        }
        else if (ttbar.t2().charge() > 0) {
          ttbar_dy = ttbar.t2().momentum().absrapidity() - ttbar.t1().momentum().absrapidity();
          _h_y_tplus->fill(ttbar.t2().momentum().rapidity());
          _h_y_tminus->fill(ttbar.t1().momentum().rapidity());
        }
        if (ttbar_dy > 1e-6) {
          _h_nTtbar_dy->fill(0.5, weight);
          _h_ttbar_mass_dy_pos->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
        }
        if (ttbar_dy < -1e-6) {
          _h_nTtbar_dy->fill(-0.5, weight);
          _h_ttbar_mass_dy_neg->fill((ttbar.t1().momentum()+ttbar.t2().momentum()).mass(), weight);
        }
      };

      void finalize() {
        asymm(_h_ttbar_mass_deta_pos, _h_ttbar_mass_deta_neg, _s_lepton_asym_ttbar_mass);
        asymm(_h_ttbar_mass_dy_pos, _h_ttbar_mass_dy_neg, _s_ttbar_asym_ttbar_mass);
        divide(_h_ttbar_mass_noaddjet, _h_ttbar_mass_addjet, _s_gap_fraction_ttbar_mass);

        normalize(_h_njets);
        normalize(_h_nbjets);
        normalize(_h_nljets);
        normalize(_h_jet_1_pT);
        normalize(_h_jet_2_pT);
        normalize(_h_jet_3_pT);
        normalize(_h_jet_4_pT);
        normalize(_h_bjet_1_pT);
        normalize(_h_bjet_2_pT);
        normalize(_h_bjet_3_pT);
        normalize(_h_bjet_4_pT);
        normalize(_h_ljet_1_pT);
        normalize(_h_ljet_2_pT);
        normalize(_h_ljet_3_pT);
        normalize(_h_ljet_4_pT);
        normalize(_h_eta_lplus);
        normalize(_h_eta_lminus);
        normalize(_h_y_tplus);
        normalize(_h_y_tminus);
        normalize(_h_top_pT);
        normalize(_h_top_mass);
        normalize(_h_ttbar_pT);
        normalize(_h_ttbar_mass);
        normalize(_h_ttjet_pT);
        normalize(_h_b1b2_mass);
        normalize(_h_lep_1_lep_2_dphi);
        normalize(_h_jet_1_jet_2_dphi);
      };

    private:

      Histo1DPtr _h_njets, _h_nbjets, _h_nljets;

      Histo1DPtr _h_jet_1_pT, _h_jet_2_pT, _h_jet_3_pT, _h_jet_4_pT;
      Histo1DPtr _h_ljet_1_pT, _h_ljet_2_pT, _h_ljet_3_pT, _h_ljet_4_pT;
      Histo1DPtr _h_bjet_1_pT, _h_bjet_2_pT, _h_bjet_3_pT, _h_bjet_4_pT;

      Histo1DPtr _h_top_pT, _h_top_mass, _h_ttbar_pT, _h_ttbar_mass, _h_ttjet_pT, _h_b1b2_mass;

      Histo1DPtr _h_lep_1_lep_2_dphi, _h_jet_1_jet_2_dphi;

      Histo1DPtr _h_nLepton_deta, _h_eta_lplus, _h_eta_lminus, _h_ttbar_mass_deta_pos, _h_ttbar_mass_deta_neg;
      Histo1DPtr _h_nTtbar_dy, _h_y_tplus, _h_y_tminus, _h_ttbar_mass_dy_pos, _h_ttbar_mass_dy_neg;
      Histo1DPtr _h_nNoAddJet, _h_ttbar_mass_noaddjet, _h_ttbar_mass_addjet;

      Scatter2DPtr _s_lepton_asym_ttbar_mass, _s_ttbar_asym_ttbar_mass, _s_gap_fraction_ttbar_mass; 

  };

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CMS_LesHouches2015> plugin_CMS_LesHouches2015;
  
  ///////////////

  // PseudoTop implementation
  // TODO: This code should be simplified using standard Rivet projections
  // TODO: Drop lepton+jets, not used in this plugin
  namespace {
  
    struct GreaterByPt
    {
      bool operator()(const Particle& a, const Particle& b) {
        return a.pt() > b.pt();
      }
    };

    void PseudoTop::cleanup(std::map<double, std::pair<size_t, size_t> >& v, const bool doCrossCleanup) const
    {
      std::vector<std::map<double, std::pair<size_t, size_t> >::const_iterator> toErase;
      std::set<size_t> usedLeg1, usedLeg2;
      if ( !doCrossCleanup ) {
        for (auto key = v.begin(); key != v.end(); ++key) {
          const size_t leg1 = key->second.first;
          const size_t leg2 = key->second.second;
          if (usedLeg1.find(leg1) == usedLeg1.end() and
              usedLeg2.find(leg2) == usedLeg2.end()) {
            usedLeg1.insert(leg1);
            usedLeg2.insert(leg2);
          } else {
            toErase.push_back(key);
          }
        }
      }
      else {
        for (auto key = v.begin(); key != v.end(); ++key) {
          const size_t leg1 = key->second.first;
          const size_t leg2 = key->second.second;
          if (usedLeg1.find(leg1) == usedLeg1.end() and
              usedLeg1.find(leg2) == usedLeg1.end()) {
            usedLeg1.insert(leg1);
            usedLeg1.insert(leg2);
          } else {
            toErase.push_back(key);
          }
        }
      }
      for (auto& key : toErase) v.erase(key);
    }

    void PseudoTop::project(const Event& e) {
      // Leptons : do the lepton clustering anti-kt R=0.1 using stable photons and leptons not from hadron decay
      // Neutrinos : neutrinos not from hadron decay
      // MET : vector sum of all invisible particles in x-y plane 
      // Jets : anti-kt R=0.4 using all particles excluding neutrinos and particles used in lepton clustering
      //        add ghost B hadrons during the jet clustering to identify B jets.

      // W->lv : dressed lepton and neutrino pairs
      // W->jj : light flavored dijet
      // W candidate : select lv or jj pairs which minimise |mW1-80.4|+|mW2-80.4|
      //               lepton-neutrino pair will be selected with higher priority

      // t->Wb : W candidate + b jet
      // t candidate : select Wb pairs which minimise |mtop1-172.5|+|mtop2-172.5|

      _isValid = false;
      _theParticles.clear();
      _wDecays1.clear();
      _wDecays2.clear();
      _jets.clear();
      _bjets.clear();
      _ljets.clear();
      _mode1 = _mode2 = CH_HADRON;

      // Collect final state particles
      Particles pForLep, pForJet;
      Particles neutrinos; // Prompt neutrinos
      foreach (const GenParticle* p, Rivet::particles(e.genEvent())) {
        const int status = p->status();
        const int pdgId = p->pdg_id();
        if (status == 1) {
          Particle rp(*p);
          if (!PID::isHadron(pdgId) && !rp.fromHadron()) {
            // Collect particles not from hadron decay
            if (rp.isNeutrino()) {
              // Prompt neutrinos are kept in separate collection
              neutrinos.push_back(rp);
            } else if (pdgId == 22 || rp.isLepton()) {
              // Leptons and photons for the dressing
              pForLep.push_back(rp);
            }
          } else if (!rp.isNeutrino()) {
            // Use all particles from hadron decay
            pForJet.push_back(rp);
          }
        } else if (PID::isHadron(pdgId) && PID::hasBottom(pdgId)) {
          // NOTE: Consider B hadrons with pT > 5GeV - not in CMS proposal
          //if ( p->momentum().perp() < 5 ) continue; 

          // Do unstable particles, to be used in the ghost B clustering
          // Use last B hadrons only
          bool isLast = true;
          foreach (GenParticle* pp, Rivet::particles(p->end_vertex(), HepMC::children)) {
            if (PID::hasBottom(pp->pdg_id())) {
              isLast = false;
              break;
            }
          }
          if (!isLast) continue;

          // Rescale momentum by 10^-20
          Particle ghost(pdgId, FourMomentum(p->momentum())*1e-20/p->momentum().rho());
          pForJet.push_back(ghost);
        }
      }

      // Start object building from trivial thing - prompt neutrinos
      std::sort(neutrinos.begin(), neutrinos.end(), GreaterByPt());

      // Proceed to lepton dressing
      FastJets fjLep(FinalState(), FastJets::ANTIKT, _lepR);
      fjLep.calc(pForLep);
      Jets leptons;
      std::vector<int> leptonsId;
      std::set<int> dressedIdxs;
      foreach (const Jet& lep, fjLep.jetsByPt(_lepMinPt)) {
        if (std::abs(lep.eta()) > _lepMaxEta) continue;

        double leadingPt = -1;
        int leptonId = 0;
        foreach (const Particle& p, lep.particles()) {
          dressedIdxs.insert(p.genParticle()->barcode());
          if (p.isLepton() && p.momentum().pt() > leadingPt) {
            leadingPt = p.momentum().pt();
            leptonId = p.pdgId();
          }
        }
        if (leptonId == 0) continue;
        leptons.push_back(lep);
        leptonsId.push_back(leptonId);
      }

      // Re-use particles not used in lepton dressing
      foreach (const Particle& rp, pForLep) {
        const int barcode = rp.genParticle()->barcode();
        // Skip if the particle is used in dressing
        if (dressedIdxs.find(barcode) != dressedIdxs.end()) continue;

        // Put back to be used in jet clustering
        pForJet.push_back(rp);
      }

      // Then do the jet clustering
      FastJets fjJet(FinalState(), FastJets::ANTIKT, _jetR);
      //fjJet.useInvisibles(); // NOTE: CMS proposal to remove neutrinos
      fjJet.calc(pForJet);
      // Jets _bjets, _ljets; FIME
      foreach (const Jet& jet, fjJet.jetsByPt(_jetMinPt)) {
        if (std::abs(jet.eta()) > _jetMaxEta) continue;
        _jets.push_back(jet);

        bool isBJet = false;
        foreach (const Particle& rp, jet.particles()) {
          if (PID::hasBottom(rp.pdgId())) {
            isBJet = true;
            break;
          }
        }

        if ( isBJet ) _bjets.push_back(jet);
        else _ljets.push_back(jet);
      }

      // Every building blocks are ready. Continue to pseudo-W and pseudo-top combination

      if (_bjets.size() < 2) return; // Ignore single top for now
      std::map<double, std::pair<size_t, size_t> > wLepCandIdxs;
      std::map<double, std::pair<size_t, size_t> > wHadCandIdxs;

      // Collect leptonic-decaying W's
      for (size_t iLep = 0, nLep = leptons.size(); iLep < nLep; ++iLep) {
        const Jet& lep = leptons.at(iLep);
        for (size_t iNu = 0, nNu = neutrinos.size(); iNu < nNu; ++iNu) {
          const Particle& nu = neutrinos.at(iNu);
          const double m = (lep.momentum()+nu.momentum()).mass();
          const double dm = std::abs(m-_wMass);
          wLepCandIdxs[dm] = make_pair(iLep, iNu);
        }
      }

      // Continue to hadronic decaying W's
      for (size_t i = 0, nLjet = _ljets.size(); i < nLjet; ++i) {
        const Jet& ljet1 = _ljets[i];
        for (size_t j = i+1; j < nLjet; ++j) {
          const Jet& ljet2 = _ljets[j];
          const double m = (ljet1.momentum()+ljet2.momentum()).mass();
          const double dm = std::abs(m-_wMass);
          wHadCandIdxs[dm] = make_pair(i, j);
        }
      }

      // Cleanup W candidate, choose pairs with minimum dm if they share decay products
      cleanup(wLepCandIdxs);
      cleanup(wHadCandIdxs, true);
      const size_t nWLepCand = wLepCandIdxs.size();
      const size_t nWHadCand = wHadCandIdxs.size();

      if (nWLepCand + nWHadCand < 2) return; // We skip single top

      int w1Q = 1, w2Q = -1;
      int w1dau1Id = 1, w2dau1Id = -1;
      FourMomentum w1dau1LVec, w1dau2LVec;
      FourMomentum w2dau1LVec, w2dau2LVec;
      if (nWLepCand == 0) { // Full hadronic case
        const auto& idPair1 = wHadCandIdxs.begin()->second;
        const auto& idPair2 = std::next(wHadCandIdxs.begin())->second;
        const auto& w1dau1 = _ljets[idPair1.first];
        const auto& w1dau2 = _ljets[idPair1.second];
        const auto& w2dau1 = _ljets[idPair2.first];
        const auto& w2dau2 = _ljets[idPair2.second];

        w1dau1LVec = w1dau1.momentum();
        w1dau2LVec = w1dau2.momentum();
        w2dau1LVec = w2dau1.momentum();
        w2dau2LVec = w2dau2.momentum();
      } else if (nWLepCand == 1) { // Semi-leptonic case
        const auto& idPair1 = wLepCandIdxs.begin()->second;
        const auto& idPair2 = wHadCandIdxs.begin()->second;
        const auto& w1dau1 = leptons[idPair1.first];
        const auto& w1dau2 = neutrinos[idPair1.second];
        const auto& w2dau1 = _ljets[idPair2.first];
        const auto& w2dau2 = _ljets[idPair2.second];

        w1dau1LVec = w1dau1.momentum();
        w1dau2LVec = w1dau2.momentum();
        w2dau1LVec = w2dau1.momentum();
        w2dau2LVec = w2dau2.momentum();
        w1dau1Id = leptonsId[idPair1.first];
        w1Q = w1dau1Id > 0 ? -1 : 1;
        w2Q = -w1Q;

        switch (w1dau1Id) {
          case 13: case -13: _mode1 = CH_MUON; break;
          case 11: case -11: _mode1 = CH_ELECTRON; break;
        }
      } else { // Full leptonic case
        const auto& idPair1 = wLepCandIdxs.begin()->second;
        const auto& idPair2 = std::next(wLepCandIdxs.begin())->second;
        const auto& w1dau1 = leptons[idPair1.first];
        const auto& w1dau2 = neutrinos[idPair1.second];
        const auto& w2dau1 = leptons[idPair2.first];
        const auto& w2dau2 = neutrinos[idPair2.second];

        w1dau1LVec = w1dau1.momentum();
        w1dau2LVec = w1dau2.momentum();
        w2dau1LVec = w2dau1.momentum();
        w2dau2LVec = w2dau2.momentum();
        w1dau1Id = leptonsId[idPair1.first];
        w2dau1Id = leptonsId[idPair2.first];
        w1Q = w1dau1Id > 0 ? -1 : 1;
        w2Q = w2dau1Id > 0 ? -1 : 1;

        switch (w1dau1Id) {
          case 13: case -13: _mode1 = CH_MUON; break;
          case 11: case -11: _mode1 = CH_ELECTRON; break;
        }
        switch (w2dau1Id) {
          case 13: case -13: _mode2 = CH_MUON; break;
          case 11: case -11: _mode2 = CH_ELECTRON; break;
        }
      }
      const auto w1LVec = w1dau1LVec+w1dau2LVec;
      const auto w2LVec = w2dau1LVec+w2dau2LVec;

      // Combine b jets
      double sumDm = 1e9;
      FourMomentum b1LVec, b2LVec;
      for (size_t i = 0, n = _bjets.size(); i < n; ++i) {
        const Jet& bjet1 = _bjets[i];
        const double mtop1 = (w1LVec+bjet1.momentum()).mass();
        const double dmtop1 = std::abs(mtop1-_tMass);
        for (size_t j=0; j<n; ++j) {
          if (i == j) continue;
          const Jet& bjet2 = _bjets[j];
          const double mtop2 = (w2LVec+bjet2.momentum()).mass();
          const double dmtop2 = std::abs(mtop2-_tMass);

          if (sumDm <= dmtop1+dmtop2) continue;

          sumDm = dmtop1+dmtop2;
          b1LVec = bjet1.momentum();
          b2LVec = bjet2.momentum();
        }
      }
      if (sumDm >= 1e9) return; // Failed to make top, but this should not happen.

      const auto t1LVec = w1LVec + b1LVec;
      const auto t2LVec = w2LVec + b2LVec;

      // Put all of them into candidate collection
      _t1 = Particle(w1Q*6, t1LVec);
      _b1 = Particle(w1Q*5, b1LVec);
      _w1 = Particle(w1Q*24, w1LVec);
      _wDecays1.push_back(Particle(w1dau1Id, w1dau1LVec));
      _wDecays1.push_back(Particle(-w1dau1Id+w1Q, w1dau2LVec));

      _t2 = Particle(w2Q*6, t2LVec);
      _b2 = Particle(w2Q*5, b2LVec);
      _w2 = Particle(w2Q*24, w2LVec);
      _wDecays2.push_back(Particle(w2dau1Id, w2dau1LVec));
      _wDecays2.push_back(Particle(-w2dau1Id+w2Q, w2dau2LVec));

      _isValid = true;
    }
  
  }

}
