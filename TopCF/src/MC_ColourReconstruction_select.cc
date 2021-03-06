// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Using colour information to identify jets in all-hadronic top decay at 13 TeV
  class MC_ColourReconstruction_select : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_ColourReconstruction_select);


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
      book(h_numberlightjets, "numberlightjets", 5, 4.5, 10.5);
      book(h_numberbjets, "numberbjets", 4, 3.5, 7.5);
      book(h_numberwjets, "numberwjets", 11, -0.5, 10.5);
      book(h_numberojets, "numberojets", 11, -0.5, 10.5);
      book(h_wjetspT, "wjetspT", 7, 20, 160);

      book(h_wcandinvariantmass, "wcandinvariantmass", 30, 40.5, 120);
      book(h_numberwjetpairs, "numberwjetpairs", 6, -0.5, 5.5);

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
          if (mu.pt() > 26*GeV)  muons.push_back(mu); //muon pT cut 
        } //vetoMuons looser cut than muons, muons subsection of vetomuons
      }

      vector<DressedLepton> elecs, vetoElecs;
      for (const DressedLepton el : all_elecs) {
        bool keep = true;
        for (const Jet j : goodJets)  keep &= deltaR(j, el) >= 0.4;
        if (keep && el.pt() > 15*GeV) {
          vetoElecs.push_back(el);;
          if (el.pt() > 35*GeV)  elecs.push_back(el); //electron pT cut 
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

      if (Ht < 500*GeV) vetoEvent; //sum of all jets pT cut, supress tt+jets background 

      if ( bJets.size() < 4 || lJets.size() < 5 )  vetoEvent; //b jets multiplicity cut
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
          double wpurity = CalculateWJetPurity(j);
          if (wpurity > 0.1) wJets += j;
        }
        else {
          oJets += j;
        }
      }
      
      h_numberwjets->fill(wJets.size());
      h_numberojets->fill(oJets.size());

      for (size_t i = 0; i < wJets.size(); ++i) {
        h_wjetspT->fill(wJets[i].pT());
      }

      Jets W1cands, W2cands, W3cands, W4cands, W5cands, W6cands;
      set<size_t> doublecount;

      double W1candmass = 40.4;
      size_t W1index1 = 100;
      size_t W1index2 = 100;
      for (size_t i = 0; i < lJets.size()-1; ++i) {
        for (size_t j = 0; j < lJets.size(); ++j) {
          double dijetinvariantmass = CalculateInvariantMass(lJets[i], lJets[j]);
          if (fabs(dijetinvariantmass - 80.4) < fabs(W1candmass - 80.4)) {
            W1candmass = dijetinvariantmass;
            W1index1 = i;
            W1index2 = j;
          }
        }
      }

      h_wcandinvariantmass->fill(W1candmass);
      
      if (W1index1 != 100 && W1index2 != 100) {
        W1cands += lJets[W1index1];
        W1cands += lJets[W1index2];
        doublecount.insert(W1index1);
        doublecount.insert(W1index2);
      }

      double W2candmass = 40.4;
      size_t W2index1 = 100;
      size_t W2index2 = 100;
      for (size_t i = 0; i < lJets.size()-1; ++i) {
        for (size_t j = 0; j < lJets.size(); ++j) {
          double dijetinvariantmass = CalculateInvariantMass(lJets[i], lJets[j]);
          if (fabs(dijetinvariantmass - 80.4) < fabs(W2candmass - 80.4) && !doublecount.count(i) && !doublecount.count(j)) {
            W2candmass = dijetinvariantmass;
            W2index1 = i;
            W2index2 = j;
          }
        }
      }

      h_wcandinvariantmass->fill(W2candmass);
      
      if (W2index1 != 100 && W2index2 != 100) {
        W2cands += lJets[W2index1];
        W2cands += lJets[W2index2];
        doublecount.insert(W2index1);
        doublecount.insert(W2index2);
      }

      double W3candmass = 40.4;
      size_t W3index1 = 100;
      size_t W3index2 = 100;
      for (size_t i = 0; i < lJets.size()-1; ++i) {
        for (size_t j = 0; j < lJets.size(); ++j) {
          double dijetinvariantmass = CalculateInvariantMass(lJets[i], lJets[j]);
          if (fabs(dijetinvariantmass - 80.4) < fabs(W3candmass - 80.4) && !doublecount.count(i) && !doublecount.count(j)) {
            W3candmass = dijetinvariantmass;
            W3index1 = i;
            W3index2 = j;
          }
        }
      }

      h_wcandinvariantmass->fill(W3candmass);
      
      if (W3index1 != 100 && W3index2 != 100) {
        W3cands += lJets[W3index1];
        W3cands += lJets[W3index2];
        doublecount.insert(W3index1);
        doublecount.insert(W3index2);
      }

      double W4candmass = 40.4;
      size_t W4index1 = 100;
      size_t W4index2 = 100;
      for (size_t i = 0; i < lJets.size()-1; ++i) {
        for (size_t j = 0; j < lJets.size(); ++j) {
          double dijetinvariantmass = CalculateInvariantMass(lJets[i], lJets[j]);
          if (fabs(dijetinvariantmass - 80.4) < fabs(W4candmass - 80.4) && !doublecount.count(i) && !doublecount.count(j)) {
            W4candmass = dijetinvariantmass;
            W4index1 = i;
            W4index2 = j;
          }
        }
      }

      h_wcandinvariantmass->fill(W4candmass);
      
      if (W4index1 != 100 && W4index2 != 100) {
        W4cands += lJets[W4index1];
        W4cands += lJets[W4index2];
        doublecount.insert(W4index1);
        doublecount.insert(W4index2);
      }

      double W5candmass = 40.4;
      size_t W5index1 = 100;
      size_t W5index2 = 100;
      for (size_t i = 0; i < lJets.size()-1; ++i) {
        for (size_t j = 0; j < lJets.size(); ++j) {
          double dijetinvariantmass = CalculateInvariantMass(lJets[i], lJets[j]);
          if (fabs(dijetinvariantmass - 80.4) < fabs(W5candmass - 80.4) && !doublecount.count(i) && !doublecount.count(j)) {
            W5candmass = dijetinvariantmass;
            W5index1 = i;
            W5index2 = j;
          }
        }
      }

      h_wcandinvariantmass->fill(W5candmass);
      
      if (W5index1 != 100 && W5index2 != 100) {
        W5cands += lJets[W5index1];
        W5cands += lJets[W5index2];
        doublecount.insert(W5index1);
        doublecount.insert(W5index2);
      }

      double W6candmass = 40.4;
      size_t W6index1 = 100;
      size_t W6index2 = 100;
      for (size_t i = 0; i < lJets.size()-1; ++i) {
        for (size_t j = 0; j < lJets.size(); ++j) {
          double dijetinvariantmass = CalculateInvariantMass(lJets[i], lJets[j]);
          if (fabs(dijetinvariantmass - 80.4) < fabs(W6candmass - 80.4) && !doublecount.count(i) && !doublecount.count(j)) {
            W5candmass = dijetinvariantmass;
            W6index1 = i;
            W6index2 = j;
          }
        }
      }
      
      h_wcandinvariantmass->fill(W6candmass);

      if (W6index1 != 100 && W6index2 != 100) {
        W6cands += lJets[W6index1];
        W6cands += lJets[W6index2];
        doublecount.insert(W6index1);
        doublecount.insert(W6index2);
      }
      

      int numberwjetpairs = 0;

      if (W1cands.size() > 1) {
        double pull_12 = CalculatePullAngle(W1cands[0], W1cands[1], 0);
        double pull_21 = CalculatePullAngle(W1cands[1], W1cands[0], 0);
        if (W1cands[0].pT() > 80*GeV && W1cands[1].pT() > 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.725 || (pull_21/Rivet::PI) < 0.725) {
            numberwjetpairs += 1;
          }
        }
        if (W1cands[0].pT() < 80*GeV || W1cands[1].pT() < 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.94 || (pull_21/Rivet::PI) < 0.94) {
            numberwjetpairs += 1;
          }
        }
      }

      if (W2cands.size() > 1) {
        double pull_12 = CalculatePullAngle(W2cands[0], W2cands[1], 0);
        double pull_21 = CalculatePullAngle(W2cands[1], W2cands[0], 0);
        if (W2cands[0].pT() > 80*GeV && W2cands[1].pT() > 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.725 || (pull_21/Rivet::PI) < 0.725) {
            numberwjetpairs += 1;
          }
        }
        if (W2cands[0].pT() < 80*GeV || W2cands[1].pT() < 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.94 || (pull_21/Rivet::PI) < 0.94) {
            numberwjetpairs += 1;
          }
        }
      }

      if (W3cands.size() > 1) {
        double pull_12 = CalculatePullAngle(W3cands[0], W3cands[1], 0);
        double pull_21 = CalculatePullAngle(W3cands[1], W3cands[0], 0);
        if (W3cands[0].pT() > 80*GeV && W3cands[1].pT() > 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.725 || (pull_21/Rivet::PI) < 0.725) {
            numberwjetpairs += 1;
          }
        }
        if (W3cands[0].pT() < 80*GeV || W3cands[1].pT() < 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.94 || (pull_21/Rivet::PI) < 0.94) {
            numberwjetpairs += 1;
          }
        }
      }

      if (W4cands.size() > 1) {
        double pull_12 = CalculatePullAngle(W4cands[0], W4cands[1], 0);
        double pull_21 = CalculatePullAngle(W4cands[1], W4cands[0], 0);
        if (W4cands[0].pT() > 80*GeV && W4cands[1].pT() > 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.725 || (pull_21/Rivet::PI) < 0.725) {
            numberwjetpairs += 1;
          }
        }
        if (W4cands[0].pT() < 80*GeV || W4cands[1].pT() < 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.94 || (pull_21/Rivet::PI) < 0.94) {
            numberwjetpairs += 1;
          }
        }
      }

      if (W5cands.size() > 1) {
        double pull_12 = CalculatePullAngle(W5cands[0], W5cands[1], 0);
        double pull_21 = CalculatePullAngle(W5cands[1], W5cands[0], 0);
        if (W5cands[0].pT() > 80*GeV && W5cands[1].pT() > 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.725 || (pull_21/Rivet::PI) < 0.725) { //signal efficiency 90%, 22% background removed
            numberwjetpairs += 1;
          }
        }
        if (W5cands[0].pT() < 80*GeV || W5cands[1].pT() < 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.94 || (pull_21/Rivet::PI) < 0.94) { //signal efficiency 90%, 15% background removed
            numberwjetpairs += 1;
          }
        }
      }

      if (W6cands.size() > 1) {
        double pull_12 = CalculatePullAngle(W6cands[0], W6cands[1], 0);
        double pull_21 = CalculatePullAngle(W6cands[1], W6cands[0], 0);
        if (W6cands[0].pT() > 80*GeV && W6cands[1].pT() > 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.725 || (pull_21/Rivet::PI) < 0.725) {
            numberwjetpairs += 1;
          }
        }
        if (W6cands[0].pT() < 80*GeV || W6cands[1].pT() < 80*GeV) {
          if ((pull_12/Rivet::PI) < 0.94 || (pull_21/Rivet::PI) < 0.94) {
            numberwjetpairs += 1;
          }
        }
      }

      h_numberwjetpairs->fill(numberwjetpairs);
      int numberwjetcands = 2*numberwjetpairs;
      h_numberwjetcands->fill(numberwjetcands);
      

      

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
 
    // Void finalize() is called after a run is finished. 
    // Here the analysis class should do whatever manipulations are necessary on the histograms.
    void finalize() {
      normalize(h_numberwjets);
      normalize(h_numberlightjets);
      normalize(h_numberbjets);
      normalize(h_numberojets);
      normalize(h_wjetspT);
      normalize(h_wcandinvariantmass);
      normalize(h_numberwjetpairs);

    }

    //@}


  private: 

    Histo1DPtr h_numberwjets;
    Histo1DPtr h_numberlightjets;
    Histo1DPtr h_numberbjets;
    Histo1DPtr h_numberojets;
    Histo1DPtr h_wjetspT;
    
    Histo1DPtr h_wcandinvariantmass;
    Histo1DPtr h_numberwjetpairs;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ColourReconstruction_select);

}
