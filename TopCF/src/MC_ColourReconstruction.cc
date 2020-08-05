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
      book(h_lightjetsinvariantmass, "lightjetsinvariantmass", 15, 0, 300);
      book(h_leading_lightjetsinvariantmass, "leading_lightjetsinvariantmass", 15, 0, 300);
      book(h_ccjetsinvariantmass, "ccjetsinvariantmass", 15, 0, 300);
      book(h_numberwjetpairs, "numberwjetpairs", 11, -0.5, 10.5);
      book(h_numberccjetpairs, "numberccjetpairs", 11, -0.5, 10.5);
      book(h_numberlightjets, "numberlightjets", 10, 4.5, 14.5);
      book(h_numberbjets, "numberbjets", 6, 1.5, 7.5);
      book(h_pull_12_min, "pull_12_min", 10, 0, 1);

      book(h_pull_12_j1j2, "pull_12_j1j2", 5, 0, 1);
      book(h_pull_12_j1j3, "pull_12_j1j3", 5, 0, 1);
      book(h_pull_12_j1j4, "pull_12_j1j4", 5, 0, 1);
      book(h_pull_12_j1j5, "pull_12_j1j5", 5, 0, 1);
      book(h_pull_12_j2j3, "pull_12_j2j3", 5, 0, 1);
      book(h_pull_12_j2j4, "pull_12_j2j4", 5, 0, 1);
      book(h_pull_12_j2j5, "pull_12_j2j5", 5, 0, 1);
      book(h_pull_12_j3j4, "pull_12_j3j4", 5, 0, 1);
      book(h_pull_12_j3j5, "pull_12_j3j5", 5, 0, 1);
      book(h_pull_12_j4j5, "pull_12_j4j5", 5, 0, 1);
    
    } 

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /**************
       *    JETS    *
       **************/
      const Jets& allJets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20.0*GeV && Cuts::absrap < 2.5); //jet pT,rap cut 
      const vector<DressedLepton>& all_elecs = apply<DressedLeptons>(event, "ELECS").dressedLeptons();
      const vector<DressedLepton>& all_muons = apply<DressedLeptons>(event, "MUONS").dressedLeptons();
      Jets goodJets;
      for (const Jet j : allJets) {
        bool keep = true;
        for (const DressedLepton el : all_elecs)  keep &= deltaR(j, el) >= 0.2; //ensure no overlap between jets and electrons
        if (keep)  goodJets += j;
      }
      if ( goodJets.size() < 7 )  vetoEvent; //jet minimum 
      
      /****************
       *    LEPTONS   *
       ****************/
      vector<DressedLepton> muons, vetoMuons;
      for (const DressedLepton mu : all_muons) {
        bool keep = true;
        for (const Jet j : goodJets)  keep &= deltaR(j, mu) >= 0.4;
        if (keep && mu.pt() > 15*GeV) {
          vetoMuons.push_back(mu);
          if (mu.pt() > 30*GeV)  muons.push_back(mu); //muon pT cut 
        } //vetoMuons looser cut than muons, muons subsection of vetomuons
      }

      vector<DressedLepton> elecs, vetoElecs;
      for (const DressedLepton el : all_elecs) {
        bool keep = true;
        for (const Jet j : goodJets)  keep &= deltaR(j, el) >= 0.4;
        if (keep && el.pt() > 15*GeV) {
          vetoElecs.push_back(el);;
          if (el.pt() > 30*GeV)  elecs.push_back(el); //electron pT cut
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
      if (met <= 20*GeV)  vetoEvent; //mono-leptonic 

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
        if (b_tagged)  bJets += j;
        if (!b_tagged) wJets += j; 
      }

      if ( bJets.size() < 2 || wJets.size() < 5 )  vetoEvent; //at least 2 b jets, at least 5 light jets 
      h_numberlightjets->fill(wJets.size());
      h_numberbjets->fill(bJets.size());

      double numberwjetpairs = 0;
      for (size_t i = 0; i < wJets.size()-1; ++i) { 
       for (size_t j = i + 1; j < wJets.size(); ++j) {
         double lightjetsinvariantmass = fabs(CalculateInvariantMass(wJets[i], wJets[j]));
         h_lightjetsinvariantmass->fill(lightjetsinvariantmass);
         if (lightjetsinvariantmass > 60*GeV && lightjetsinvariantmass < 100*GeV) { 
           numberwjetpairs += 1;
         }
       }
      }

      h_numberwjetpairs->fill(numberwjetpairs);

      double leading_lightjetsinvariantmass = fabs(CalculateInvariantMass(wJets[0], wJets[1]));
      h_leading_lightjetsinvariantmass->fill(leading_lightjetsinvariantmass);

      //Finding colour connected jets from light jets (non-b)
      double numberccjetpairs = 0;
      //vector<size_t> ccjet2index;
      set<size_t> doublecount; //to prevent double counting of jet pairs in loop: if jet a is cc with jet b, then jet b cannot be cc with jet c 
      for (size_t i = 0; i < wJets.size()-1; ++i) { 
        double pull_12_min = Rivet::PI; 
        size_t min_pull_12_jet2index = 0;
        for (size_t j = 0; j < wJets.size()-1; ++j) {
          if (i != j) {
            double pull_12 = CalculatePullAngle(wJets[i], wJets[j], 0);
            double pull_21 = CalculatePullAngle(wJets[j], wJets[i], 0);
            if (pull_12 < pull_12_min && (pull_21 / Rivet::PI) < 0.5) {
              pull_12_min = pull_12;
              min_pull_12_jet2index = j;
            }
          }
       }
       h_pull_12_min->fill(pull_12_min / Rivet::PI);
        //ccjet2index.push_back(min_pull_12_jet2index);
        if ((pull_12_min / Rivet::PI) < 0.5 && !doublecount.count(i) && !doublecount.count(min_pull_12_jet2index)) { // if cc jet candidate found for jet i 
          numberccjetpairs += 1; 
          double ccjetsinvariantmass = fabs(CalculateInvariantMass(wJets[i], wJets[min_pull_12_jet2index]));
          h_ccjetsinvariantmass->fill(ccjetsinvariantmass);
          doublecount.insert(i);
          doublecount.insert(min_pull_12_jet2index);
        }
      }


      double pull_12_j1j2 = CalculatePullAngle(wJets[0], wJets[1], 0);
      h_pull_12_j1j2->fill(pull_12_j1j2 / Rivet::PI);
      double pull_12_j1j3 = CalculatePullAngle(wJets[0], wJets[2], 0);
      h_pull_12_j1j3->fill(pull_12_j1j3 / Rivet::PI);
      double pull_12_j1j4 = CalculatePullAngle(wJets[0], wJets[3], 0);
      h_pull_12_j1j4->fill(pull_12_j1j4 / Rivet::PI);
      double pull_12_j1j5 = CalculatePullAngle(wJets[0], wJets[4], 0);
      h_pull_12_j1j5->fill(pull_12_j1j5 / Rivet::PI);
      double pull_12_j2j3 = CalculatePullAngle(wJets[1], wJets[2], 0);
      h_pull_12_j2j3->fill(pull_12_j2j3 / Rivet::PI);
      double pull_12_j2j4 = CalculatePullAngle(wJets[1], wJets[3], 0);
      h_pull_12_j2j4->fill(pull_12_j2j4 / Rivet::PI);
      double pull_12_j2j5 = CalculatePullAngle(wJets[1], wJets[4], 0);
      h_pull_12_j2j5->fill(pull_12_j2j5 / Rivet::PI);
      double pull_12_j3j4 = CalculatePullAngle(wJets[2], wJets[3], 0);
      h_pull_12_j3j4->fill(pull_12_j3j4 / Rivet::PI);
      double pull_12_j3j5 = CalculatePullAngle(wJets[2], wJets[4], 0);
      h_pull_12_j3j5->fill(pull_12_j3j5 / Rivet::PI);
      double pull_12_j4j5 = CalculatePullAngle(wJets[3], wJets[4], 0);
      h_pull_12_j4j5->fill(pull_12_j4j5 / Rivet::PI);

      //Make sure both jets are identified as colour connected with each other both ways 
      //Then calculate dijet invariant mass
      //set<size_t> doublecount; //to prevent double counting of jet pairs in loop: if jet a is cc with jet b, then jet b cannot be cc with jet c 
      //for (size_t i = 0; i < wJets.size()-1; ++i) {
        //size_t j = ccjet2index[i];
        //if (j != 0) { //colour connected jet candidate found 
          //if (ccjet2index[j] == i && !doublecount.count(i) && !doublecount.count(j)) {
            //numberccjetpairs += 1; 
            //double ccjetsinvariantmass = fabs(CalculateInvariantMass(wJets[i], wJets[j]));
            //h_ccjetsinvariantmass->fill(ccjetsinvariantmass);
            //doublecount.insert(i);
            //doublecount.insert(j);
          //}
        //}
      //}
      
      h_numberccjetpairs->fill(numberccjetpairs);

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
 
    // Void finalize() is called after a run is finished. 
    // Here the analysis class should do whatever manipulations are necessary on the histograms.
    void finalize() {
      normalize(h_lightjetsinvariantmass);
      normalize(h_leading_lightjetsinvariantmass);
      normalize(h_ccjetsinvariantmass);
      normalize(h_numberwjetpairs);
      normalize(h_numberccjetpairs);
      normalize(h_numberlightjets);
      normalize(h_numberbjets);
      normalize(h_pull_12_min);

      normalize(h_pull_12_j1j2);
      normalize(h_pull_12_j1j3);
      normalize(h_pull_12_j1j4);
      normalize(h_pull_12_j1j5);
      normalize(h_pull_12_j2j3);
      normalize(h_pull_12_j2j4);
      normalize(h_pull_12_j2j5);
      normalize(h_pull_12_j3j4);
      normalize(h_pull_12_j3j5);
      normalize(h_pull_12_j4j5);

    }

    //@}


  private:

    Histo1DPtr h_lightjetsinvariantmass;
    Histo1DPtr h_leading_lightjetsinvariantmass;
    Histo1DPtr h_ccjetsinvariantmass;
    Histo1DPtr h_numberwjetpairs;
    Histo1DPtr h_numberccjetpairs;
    Histo1DPtr h_numberlightjets;
    Histo1DPtr h_numberbjets;
    Histo1DPtr h_pull_12_min;

    Histo1DPtr h_pull_12_j1j2;
    Histo1DPtr h_pull_12_j1j3;
    Histo1DPtr h_pull_12_j1j4;
    Histo1DPtr h_pull_12_j1j5;
    Histo1DPtr h_pull_12_j2j3;
    Histo1DPtr h_pull_12_j2j4;
    Histo1DPtr h_pull_12_j2j5;
    Histo1DPtr h_pull_12_j3j4;
    Histo1DPtr h_pull_12_j3j5;
    Histo1DPtr h_pull_12_j4j5;

    //const vector<double> pull_bins = {0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1};

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ColourReconstruction);

}
