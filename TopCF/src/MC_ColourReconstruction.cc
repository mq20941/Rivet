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
  class MC_ColourReconstruction : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_ColourReconstruction);


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
      //book(h_numberwjets, "numberwjets", 5, -0.5, 4.5);
      book(h_numberbjets, "numberbjets", 5, -0.5, 4.5);
      book(h_numberlightjets, "numberlightjets", 10, -0.5, 9.5);
      //book(h_numberbjets_sanity, "numberbjets_sanity", 5, -0.5, 4.5);
      book(h_leadingjetsinvariantmass, "leadingjetsinvariantmass", 150, 0, 300);
      book(h_leadingjets_pull_12, "leadingjets_pull_12", 5, 0, 1);
      book(h_jets13invariantmass, "jets13invariantmass", 30, 0, 300);
      book(h_jets13_pull_12, "jets13_pull_12", 5, 0, 1);
      book(h_jets23invariantmass, "jets23invariantmass", 30, 0, 300);
      book(h_jets23_pull_12, "jets23_pull_12", 5, 0, 1);
      book(h_wjetsinvariantmass, "wjetsinvariantmass", 150, 0, 300);
      book(h_wnonwjets_pull_12, "wnonwjets_pull_12", 5, 0, 1);
      book(h_nonwjets_pull_12, "nonwjets_pull_12", 5, 0, 1);
      book(h_wjet1index, "wjet1index", 10, -0.5, 9.5);
      book(h_wjet2index, "wjet2index", 10, -0.5, 9.5);
      book(h_wjets_pull_12, "wjets_pull_12", 5, 0, 1);
      book(h_lightjetsinvariantmass, "lightjetsinvariantmass", 300, 0, 300);
      book(h_lightjets_pull_12, "lightjets_pull_12", 5, 0, 1);

    } 

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /**************
       *    JETS    *
       **************/
      const Jets& allJets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20.0*GeV && Cuts::absrap < 2.5); //jet pT,rap cut //80pT jet pT cut makes pull angle smaller
      const vector<DressedLepton>& all_elecs = apply<DressedLeptons>(event, "ELECS").dressedLeptons();
      const vector<DressedLepton>& all_muons = apply<DressedLeptons>(event, "MUONS").dressedLeptons();
      Jets goodJets;
      for (const Jet j : allJets) {
        bool keep = true;
        for (const DressedLepton el : all_elecs)  keep &= deltaR(j, el) >= 0.2; //ensure no overlap between jets and electrons
        if (keep)  goodJets += j;
      }
      if ( goodJets.size() < 4 )  vetoEvent; //jet multiplicity minimum 
      
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
      if (met <= 20*GeV)  vetoEvent; //supress residual QCD background, mono-leptonic events 

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

      if ( bJets.size() < 2 || lJets.size() < 2 )  vetoEvent; //b jets multiplicity cut
      h_numberlightjets->fill(lJets.size());
      h_numberbjets->fill(bJets.size());
      

      //Find W jets
      //Jets wJets;
      //for(Jet j : lJets) {
      //  if (fromW(j)) wJets +=j;
      //}
      //h_numberwjets->fill(wJets.size());

      //Jets bJets_sanity;
      //for (Jet j : goodJets) {
      //  if (fromBottom(j)) bJets_sanity +=j;
      //}
      //h_numberbjets_sanity->fill(bJets_sanity.size());

      double leadingjetsinvariantmass = CalculateInvariantMass(lJets[0], lJets[1]);
      h_leadingjetsinvariantmass->fill(leadingjetsinvariantmass);
      double leadingjets_pull_12 = CalculatePullAngle(lJets[0], lJets[1], 0);
      h_leadingjets_pull_12->fill(leadingjets_pull_12/Rivet::PI);

      if (lJets.size() >= 3) {
        double jets13invariantmass = CalculateInvariantMass(lJets[0], lJets[2]);
        h_jets13invariantmass->fill(jets13invariantmass);
        double jets13_pull_12 = CalculatePullAngle(lJets[0], lJets[2], 0);
        h_jets13_pull_12->fill(jets13_pull_12/Rivet::PI);

        double jets23invariantmass = CalculateInvariantMass(lJets[1], lJets[2]);
        h_jets23invariantmass->fill(jets23invariantmass);
        double jets23_pull_12 = CalculatePullAngle(lJets[1], lJets[2], 0);
        h_jets23_pull_12->fill(jets23_pull_12/Rivet::PI);
      }

      Jets wJets;
      double Wcandmass = 1000*GeV;
      size_t wjet1index = 100;
      size_t wjet2index = 100;
      for (size_t i = 0; i < lJets.size()-1; ++i) {
        for (size_t j = i + 1; j < lJets.size(); ++j) {
          double dijetinvariantmass = CalculateInvariantMass(lJets[i], lJets[j]);
          h_lightjetsinvariantmass->fill(dijetinvariantmass);
          double lightjets_pull_12 = CalculatePullAngle(lJets[i], lJets[j], 0);
          h_lightjets_pull_12->fill(lightjets_pull_12/Rivet::PI);
          if (fabs(dijetinvariantmass - 80.4*GeV) < fabs(Wcandmass - 80.4*GeV)) {
            Wcandmass = dijetinvariantmass;
            wjet1index = i;
            wjet2index = j;
          }
        }
      }
      h_wjetsinvariantmass->fill(fabs(Wcandmass));
      h_wjet1index->fill(wjet1index);
      h_wjet2index->fill(wjet2index);

      if (wjet1index != 100 && wjet2index != 100) {
        wJets += lJets[wjet1index];
        wJets += lJets[wjet2index];

        double wjets_pull_12 = CalculatePullAngle(wJets[0], wJets[1], 0);
        h_wjets_pull_12->fill(wjets_pull_12/Rivet::PI);

        for (size_t i = 0; i < lJets.size()-1; ++i) {
          if (i != wjet1index && i != wjet2index) {
            double wnonwjets_pull_12 = CalculatePullAngle(wJets[0], lJets[i], 0);
            h_wnonwjets_pull_12->fill(wnonwjets_pull_12/Rivet::PI);
            for (size_t j = i + 1; j < lJets.size(); ++j) {
              if (j != wjet1index && j != wjet2index) {
                double nonwjets_pull_12 = CalculatePullAngle(lJets[i], lJets[j], 0);
                h_nonwjets_pull_12->fill(nonwjets_pull_12/Rivet::PI);
              }
            }
          }
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

   bool fromW(Jet& jet) {
     bool w_tagged = false;
     Particles& constituents = jet.particles();
     for (Particle p : constituents) {
       Particles ancestors = p.ancestors();
       for (Particle m : ancestors) {
         bool ContainsW = false;
         const PdgId pid = m.pid();
         if (pid == PID::WMINUS || pid == PID::WPLUS) ContainsW = true;
         w_tagged|= ContainsW;
       }
     }
     return w_tagged;
    }

    bool fromBottom(Jet& jet) {
     bool b_tagged = false;
     Particles& constituents = jet.particles();
     for (Particle p : constituents) {
       Particles ancestors = p.ancestors();
       for (Particle m : ancestors) {
         bool Containsb = false;
         const PdgId pid = m.pid();
         if (PID::isHadron(pid) && PID::hasBottom(pid)) Containsb = true;
         b_tagged |= Containsb;
       }
     }
     return b_tagged;
    }
   

    // Void finalize() is called after a run is finished. 
    // Here the analysis class should do whatever manipulations are necessary on the histograms.
    void finalize() {
      //normalize(h_numberwjets);
      normalize(h_numberlightjets);
      normalize(h_numberbjets);
      //normalize(h_numberbjets_sanity);
      normalize(h_leadingjetsinvariantmass);
      normalize(h_leadingjets_pull_12);
      normalize(h_jets13invariantmass);
      normalize(h_jets13_pull_12);
      normalize(h_jets23invariantmass);
      normalize(h_jets23_pull_12);
      normalize(h_wjetsinvariantmass);
      normalize(h_wnonwjets_pull_12);
      normalize(h_nonwjets_pull_12);
      normalize(h_wjet1index);
      normalize(h_wjet2index);
      normalize(h_wjets_pull_12);
      normalize(h_lightjetsinvariantmass);
      normalize(h_lightjets_pull_12);
    }

    //@}

  private:

    //Histo1DPtr h_numberwjets;
    Histo1DPtr h_numberlightjets;
    Histo1DPtr h_numberbjets;
    //Histo1DPtr h_numberbjets_sanity;
    Histo1DPtr h_leadingjetsinvariantmass;
    Histo1DPtr h_leadingjets_pull_12;
    Histo1DPtr h_jets13invariantmass;
    Histo1DPtr h_jets13_pull_12;
    Histo1DPtr h_jets23invariantmass;
    Histo1DPtr h_jets23_pull_12;
    Histo1DPtr h_wjetsinvariantmass;
    Histo1DPtr h_wnonwjets_pull_12;
    Histo1DPtr h_nonwjets_pull_12;
    Histo1DPtr h_wjet1index;
    Histo1DPtr h_wjet2index;
    Histo1DPtr h_wjets_pull_12;
    Histo1DPtr h_lightjetsinvariantmass;
    Histo1DPtr h_lightjets_pull_12;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ColourReconstruction);

}
