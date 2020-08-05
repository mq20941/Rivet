#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/AnalysisLoader.hh"

namespace Rivet {

  class MC_TOPMASS_LJETS : public Analysis {
    public:

      /// Minimal constructor
      MC_TOPMASS_LJETS() : Analysis("MC_TOPMASS_LJETS")
    {
    }


      /// @name Analysis methods
      //@{

      
      void bookHistos(std::string postfix = "") {
        _h["weight"+postfix] = bookHisto1D("weight"+postfix, 200, -2., 2.);
        _h["scale"+postfix]  = bookHisto1D("scale"+postfix, logspace(50, 100.0, 1000.0));
        _h["stage"+postfix]  = bookHisto1D("stage"+postfix, 10, 0, 10);
        //
        _h["nbjet"+postfix]    = bookHisto1D("nbjet"+postfix, 10, 0, 10);
        _h["nljet"+postfix]    = bookHisto1D("nljet"+postfix, 10, 0, 10);
        _h["delta_mt"+postfix] = bookHisto1D("delta_mt"+postfix, 100, -100., 100.);
        //
        _h["wlep_mass"+postfix]       = bookHisto1D("wlep_mass"+postfix, 40, 60.4, 100.4);
        _h["wlep_pt"+postfix]         = bookHisto1D("wlep_pt"+postfix, 50, 0., 400.);
        _h["wlep_pt_scaled"+postfix]  = bookHisto1D("wlep_pt_scaled"+postfix, 50, 0., 400.);
        _h["whad_mass"+postfix]       = bookHisto1D("whad_mass"+postfix, 75, 30., 180.);
        _h["whad_pt"+postfix]         = bookHisto1D("whad_pt"+postfix, 50, 0., 400.);
        _h["whad_pt_scaled"+postfix]  = bookHisto1D("whad_pt_scaled"+postfix, 50, 0., 400.);
        _h["tlep_mass"+postfix]       = bookHisto1D("tlep_mass"+postfix, 100, 100., 300.);
        _h["tlep_pt"+postfix]         = bookHisto1D("tlep_pt"+postfix, 50, 0., 400.);
        _h["thad_mass"+postfix]       = bookHisto1D("thad_mass"+postfix, 100, 100., 300.);
        _h["thad_pt"+postfix]         = bookHisto1D("thad_pt"+postfix, 50, 0., 400.);
        _h["ttbar_mass"+postfix]      = bookHisto1D("ttbar_mass"+postfix, 50, 0., 1000.);
        _h["ttbar_pt"+postfix]        = bookHisto1D("ttbar_pt"+postfix, 50, 0., 400.);
      }
      
      /// Set up projections
      void init() {

        // Cuts
        Cut particle_cut = (Cuts::abseta < 5.0) and (Cuts::pT > 0.0*MeV);
        Cut lepton_cut   = (Cuts::abseta < 2.1) and (Cuts::pT > 30.*GeV);
        // Jet cut needs to be applied in analyze() function
        // Cut jet_cut = (Cuts::abseta < 2.4) and (Cuts::pT > 30.*GeV);
        
        // Generic final state
        FinalState fs(particle_cut);
        
        // Dressed leptons
        ChargedLeptons charged_leptons(fs);
        PromptFinalState prompt_leptons(charged_leptons);
        IdentifiedFinalState photons(fs);
        photons.acceptIdPair(PID::PHOTON);
        PromptFinalState prompt_photons(photons); //photons not from hadrons
        DressedLeptons dressed_leptons(prompt_photons, prompt_leptons, 0.1, lepton_cut, true);
        declare(dressed_leptons, "DressedLeptons");
        
        // Neutrinos
        IdentifiedFinalState neutrinos(fs);
        neutrinos.acceptNeutrinos();
        PromptFinalState prompt_neutrinos(neutrinos);
        declare(prompt_neutrinos, "Neutrinos");
        
        // Jets
        VetoedFinalState jet_fs(fs);
        jet_fs.addVetoOnThisFinalState(dressed_leptons);
        declare(FastJets(jet_fs, FastJets::ANTIKT, 0.4), "Jets");
        declare(FastJets(jet_fs, FastJets::ANTIKT, 0.4, JetAlg::ALL_MUONS, JetAlg::DECAY_INVISIBLES), "Jets_WithNu");
        
        // MET
        declare(MissingMomentum(fs), "MET");

        // Booking of histograms
        bookHistos();
        bookHistos("_WithNu");
      }


      void analyzeWithPostfix(const Event& event, std::string postfix = "") {
        const double wmass = 80.4;
        
        const double weight = event.weight();
        _h["weight"+postfix]->fill(weight);
        _h["scale"+postfix]->fill(event.genEvent()->event_scale(), weight);
        int stage = 0;
        _h["stage"+postfix]->fill(stage++, weight);

        // Get analysis objects from projections
        
        // Exactly 1 lepton and 1 matching neutrino for lepton+jets channel
        const std::vector<DressedLepton>& leptons = applyProjection<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
        if (leptons.size() != 1)
          vetoEvent;
        _h["stage"+postfix]->fill(stage++, weight);
        const DressedLepton lepton = leptons[0];

        const Particles neutrinos = apply<PromptFinalState>(event, "Neutrinos").particlesByPt();
        Particle neutrino;
        FourMomentum wlep;
        for (const Particle& nu : neutrinos) {
          if (nu.abspid() == lepton.abspid()+1 and nu.pid()*lepton.pid()<0) {
            wlep = lepton.momentum() + nu.momentum();
            _h["wlep_mass"+postfix]->fill(wlep.mass(), weight);
          }
          else
            continue;
          if (inRange(wlep.mass(), wmass-5., wmass+5.)) {
            neutrino = nu;
            break;
          }
        }
        if (neutrino.pid() == 0)
          vetoEvent;
        _h["stage"+postfix]->fill(stage++, weight);
        
        // Select events with at least 4 jets (2 b-tagged, 2 untagged)
        Cut jet_cut = (Cuts::abseta < 2.4) and (Cuts::pT > 30.*GeV);
        const Jets& jets = apply<FastJets>(event, "Jets"+postfix).jetsByPt(jet_cut);
        Jets bjets, ljets;
        for (const Jet& jet : jets) {
          if (jet.bTagged()) bjets.push_back(jet);
          else               ljets.push_back(jet);
        }
        _h["nbjet"+postfix]->fill(bjets.size(), weight);
        _h["nljet"+postfix]->fill(ljets.size(), weight);
        if (bjets.size() < 2 or ljets.size() < 2)
          vetoEvent;
        _h["stage"+postfix]->fill(stage++, weight);
        
        // Reconstruct top quarks
        // Constrain W candidate masses to 80.4 GeV
        // Top candidates should have similar masses
        FourMomentum whad, thad, tlep;
        double kmin = numeric_limits<double>::max();
        for (const Jet& ljet1 : ljets) {
          for (const Jet& ljet2 : ljets) {
            if (&ljet1 == &ljet2 or ljet1.pt() < ljet2.pt())
              continue;
            FourMomentum whad_cand(ljet1.momentum() + ljet2.momentum());
            for (const Jet& bjethad : bjets) {
              FourMomentum thad_cand(whad_cand * wmass/whad_cand.mass() + bjethad.momentum());
              for (const Jet& bjetlep : bjets) {
                if (&bjethad == &bjetlep)
                  continue;
                FourMomentum tlep_cand(wlep * wmass/wlep.mass() + bjetlep.momentum());

                double K = pow(whad_cand.mass() - wmass, 2) + pow(thad_cand.mass() - tlep_cand.mass(), 2);
                if(K < kmin) {
                  kmin = K;
                  whad = whad_cand;
                  thad = thad_cand;
                  tlep = tlep_cand;
                }
              }
            }
          }
        }
        FourMomentum ttbar = thad + tlep;
        
        _h["whad_mass"+postfix]->fill(whad.mass(), weight);
        double delta_mt = thad.mass() - tlep.mass();
        _h["delta_mt"+postfix]->fill(delta_mt, weight);
        if (not(inRange(whad.mass(), wmass-10., wmass+10.) and abs(delta_mt)<20.))
          vetoEvent;
        _h["stage"+postfix]->fill(stage++, weight);
        
        // Fill plots
        _h["wlep_pt"+postfix]->fill(wlep.pt(), weight);
        _h["wlep_pt_scaled"+postfix]->fill(wlep.pt() * wmass/wlep.mass(), weight);
        _h["whad_pt"+postfix]->fill(whad.pt(), weight);
        _h["whad_pt_scaled"+postfix]->fill(whad.pt() * wmass/whad.mass(), weight);
        _h["tlep_mass"+postfix]->fill(tlep.mass(), weight);
        _h["tlep_pt"+postfix]->fill(tlep.pt(), weight);
        _h["thad_mass"+postfix]->fill(thad.mass(), weight);
        _h["thad_pt"+postfix]->fill(thad.pt(), weight);
        _h["ttbar_mass"+postfix]->fill(ttbar.mass(), weight);
        _h["ttbar_pt"+postfix]->fill(ttbar.pt(), weight);

      }
      
      void analyze(const Event& event) {
        analyzeWithPostfix(event);
        analyzeWithPostfix(event, "_WithNu");
      }


      void finalize() {
      }

      //@}


    private:

      // @name Histogram data members
      //@{

      map<string, Histo1DPtr> _h;

      //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TOPMASS_LJETS);

}
