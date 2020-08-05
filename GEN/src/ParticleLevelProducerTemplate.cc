#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {

  // This analysis template can be used to obtain the same physics objects
  // as the ParticleLevelProducer added to GeneratorInterface/RivetInterface.
  // See CMSSW pull requests #18402 (80X) and #18404 (master).
  class ParticleLevelProducerTemplate : public Analysis {

  public:
    ParticleLevelProducerTemplate() : Analysis("ParticleLevelProducerTemplate")
    {
    }

    // Initialize Rivet projections
    void init() {
      // TODO: Please shorten the code below to match your analysis settings
      //       * replace the config variables by their actual values
      //       * remove unused parts of if-else statements
      //       * remove unused projections
      //       * remove unused includes on top of this file
      
      bool _usePromptFinalStates = true;
      bool _excludePromptLeptonsFromJetClustering = false;
      bool _excludeNeutrinosFromJetClustering = true;
      
      double _particleMinPt = 0.;
      double _particleMaxEta = 5.;
      
      double _lepConeSize = 0.1;
      double _lepMinPt = 15.;
      double _lepMaxEta = 2.5;
      
      double _jetConeSize = 0.4;
      double _fatJetConeSize = 0.8;
      
      // Cuts
      Cut particle_cut = (Cuts::abseta < _particleMaxEta) and (Cuts::pT > _particleMinPt*GeV);
      Cut lepton_cut   = (Cuts::abseta < _lepMaxEta)      and (Cuts::pT > _lepMinPt*GeV);
      
      // Generic final state
      FinalState fs(particle_cut);
      
      // Dressed leptons
      ChargedLeptons charged_leptons(fs);
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      
      PromptFinalState prompt_leptons(charged_leptons);
      prompt_leptons.acceptMuonDecays(true);
      prompt_leptons.acceptTauDecays(true);
      
      PromptFinalState prompt_photons(photons);
      prompt_photons.acceptMuonDecays(true);
      prompt_photons.acceptTauDecays(true);
      
      // useDecayPhotons=true allows for photons with tau ancestor,
      // photons from hadrons are vetoed by the PromptFinalState;
      // will be default DressedLeptons behaviour for Rivet >= 2.5.4
      DressedLeptons dressed_leptons(prompt_photons, prompt_leptons, _lepConeSize, 
                     lepton_cut, /*useDecayPhotons*/ true);
      if (not _usePromptFinalStates)
        dressed_leptons = DressedLeptons(photons, charged_leptons, _lepConeSize, 
                          lepton_cut, /*useDecayPhotons*/ true);
      declare(dressed_leptons, "DressedLeptons");
      
      // Photons
      if (_usePromptFinalStates) {
        // We remove the photons used up for lepton dressing in this case
        VetoedFinalState vetoed_prompt_photons(prompt_photons);
        vetoed_prompt_photons.addVetoOnThisFinalState(dressed_leptons);
        declare(vetoed_prompt_photons, "Photons");
      }
      else
        declare(photons, "Photons");
      
      // Jets
      VetoedFinalState fsForJets(fs);
      if (_usePromptFinalStates and _excludePromptLeptonsFromJetClustering)
        fsForJets.addVetoOnThisFinalState(dressed_leptons);
      JetAlg::Invisibles invisiblesStrategy = JetAlg::Invisibles::DECAY;
      if (_excludeNeutrinosFromJetClustering)
        invisiblesStrategy = JetAlg::Invisibles::NONE;
      declare(FastJets(fsForJets, FastJets::ANTIKT, _jetConeSize,
                             JetAlg::Muons::ALL, invisiblesStrategy), "Jets");
      
      // FatJets
      declare(FastJets(fsForJets, FastJets::ANTIKT, _fatJetConeSize), "FatJets");
      
      // Neutrinos
      IdentifiedFinalState neutrinos(fs);
      neutrinos.acceptNeutrinos();
      if (_usePromptFinalStates) {
        PromptFinalState prompt_neutrinos(neutrinos);
        prompt_neutrinos.acceptMuonDecays(true);
        prompt_neutrinos.acceptTauDecays(true);
        declare(prompt_neutrinos, "Neutrinos");
      }
      else
        declare(neutrinos, "Neutrinos");
      
      // MET
      declare(MissingMomentum(fs), "MET");
      
      // Booking of histograms
      book(_h["n_lep"] ,  10, -0.5,   9.5);
      book(_h["lep_pt"], 100,  0.0, 200.0);
      book(_h["n_jet"] ,  10, -0.5,   9.5);
      book(_h["jet_pt"], 100,  0.0, 500.0);
      book(_h["n_bjet"],  10, -0.5,   9.5);
    };

    // Apply Rivet projections
    void analyze(const Event& event) {
      // TODO: Please shorten the code below to match your analysis settings
      //       * replace the config variables by their actual values
      //       * remove unused projections
      
      double _jetMinPt = 30.;
      double _jetMaxEta = 2.4;
      
      double _fatJetMinPt = 170.;
      double _fatJetMaxEta = 2.4;
      
      // Actual analyze code begins here
      
      // Leptons
      const std::vector<DressedLepton>& leptons = applyProjection<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      _h["n_lep"]->fill(leptons.size());
      for (const DressedLepton& lepton : leptons) {
        _h["lep_pt"]->fill(lepton.pt());
      }
      
      // Jets
      Cut jet_cut = (Cuts::abseta < _jetMaxEta) and (Cuts::pT > _jetMinPt*GeV);
      const Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(jet_cut);
      _h["n_jet"]->fill(jets.size());
      Jets bjets;
      for (const Jet& jet : jets) {
        _h["jet_pt"]->fill(jet.pt());
        if (jet.bTagged()) bjets.push_back(jet);
      }
      _h["n_bjet"]->fill(bjets.size());
      
      Cut fatjet_cut = (Cuts::abseta < _fatJetMaxEta) and (Cuts::pT > _fatJetMinPt*GeV);
      const Jets fatjets   = applyProjection<FastJets>(event, "FatJets").jetsByPt(fatjet_cut);
      
      // Other objects
      const Particles photons   = applyProjection<FinalState>(event, "Photons").particlesByPt();
      const Particles neutrinos = applyProjection<FinalState>(event, "Neutrinos").particlesByPt();
      const FourMomentum met    = applyProjection<MissingMomentum>(event, "MET").missingMomentum();
    };

    // Finalize histograms here
    void finalize() {
      // Just normalize
      // normalize({_h["n_lep"], _h["lep_pt"], _h["n_jet"], _h["jet_pt"], _h["n_bjet"]});
      
      // Scale according to cross section
      scale({_h["n_lep"], _h["lep_pt"], _h["n_jet"], _h["jet_pt"], _h["n_bjet"]}, crossSection()/picobarn / sumOfWeights());
    };

    private:

      // @name Histogram data members
      //@{

      map<string, Histo1DPtr> _h;

      //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ParticleLevelProducerTemplate);
}
