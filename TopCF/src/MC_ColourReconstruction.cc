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
      book(h_lightjetsinvariantmass, "lightjetsinvariantmass", 20, 0, 300);
      book(h_leading_lightjetsinvariantmass, "leading_lightjetsinvariantmass", 20, 0, 300);
      book(h_numberwjetpairs, "numberwjetpairs", 11, -0.5, 10.5);
      book(h_numberccjetpairs_s1, "numberccjetpairs_s1", 6, -0.5, 5.5);
      book(h_numberccjetpairs_s2, "numberccjetpairs_s2", 6, -0.5, 5.5);
      book(h_numberccjetpairs_s3, "numberccjetpairs_s3", 6, -0.5, 5.5);
      book(h_numberlightjets, "numberlightjets", 6, 2.5, 8.5);
      book(h_numberbjets, "numberbjets", 4, 1.5, 5.5);
      book(h_numberwjets, "numberwjets", 7, -0.5, 6.5);
      book(h_pull_wj1wj2, "pull_wj1wj2", 6, 0, 1);
      book(h_Ht, "Ht", 20, 500, 1500);
      book(h_ccjetsinvariantmass_s1, "ccjetsinvariantmass_s1", 20, 0, 300);
      book(h_ccjetsinvariantmass_s2, "ccjetsinvariantmass_s2", 20, 0, 300);
      book(h_ccjetsinvariantmass_s3, "ccjetsinvariantmass_s3", 20, 0, 300);
      book(h_pull_12_s1, "pull_12_s1", 6, 0, 1);
      book(h_pull_21_s1, "pull_21_s1", 6, 0, 1);
      book(h_pull_12_s2, "pull_12_s2", 6, 0, 1);
      book(h_pull_21_s2, "pull_21_s2", 6, 0, 1);
      book(h_pull_12_s3, "pull_12_s3", 6, 0, 1);
      book(h_pull_21_s3, "pull_21_s3", 6, 0, 1);

    
    } 

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /**************
       *    JETS    *
       **************/
      const Jets& allJets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 80.0*GeV && Cuts::absrap < 2.5); //jet pT,rap cut //80pT jet pT cut makes pull angle smaller
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

      if (Ht < 500*GeV) vetoEvent; //sum of all jets pT cut, supress tt+jets background 
      h_Ht->fill(Ht); 

      if ( bJets.size() < 2 || lJets.size() < 2 )  vetoEvent; //number b jets cut
      h_numberlightjets->fill(lJets.size());
      h_numberbjets->fill(bJets.size());
      

      //Find W jets
      Jets wJets;
      for(Jet j : lJets) {
        bool w_tagged = false;
        Particles& constituents = j.particles();
        for ( Particle p : constituents) {
          w_tagged |= p.hasAncestorWith(Cuts::abspid == PID::TQUARK);
        }
        if (w_tagged) wJets += j;
      }

      if (wJets.size() >= 2) {
        double pull_wj1wj2 = CalculatePullAngle(wJets[0], wJets[1], 0);
        h_pull_wj1wj2->fill(pull_wj1wj2 / Rivet::PI);
      }
      h_numberwjets->fill(wJets.size());



      double numberwjetpairs = 0;
      for (size_t i = 0; i < lJets.size()-1; ++i) { 
       for (size_t j = i + 1; j < lJets.size(); ++j) {
         double lightjetsinvariantmass = fabs(CalculateInvariantMass(lJets[i], lJets[j]));
         h_lightjetsinvariantmass->fill(lightjetsinvariantmass);
         if (lightjetsinvariantmass > 60*GeV && lightjetsinvariantmass < 100*GeV) { 
           numberwjetpairs += 1;
         }
       }
      }

      h_numberwjetpairs->fill(numberwjetpairs);

      double leading_lightjetsinvariantmass = fabs(CalculateInvariantMass(lJets[0], lJets[1]));
      h_leading_lightjetsinvariantmass->fill(leading_lightjetsinvariantmass);
      
      //PULL ANGLE SELECTION 1:
      //Finding colour connected jets from light jets (non-b)
      //Must be colour connected in both directions: if jet a cc with jet b, then jet b must be cc with jet a
      double numberccjetpairs_s1 = 0;
      vector<size_t> ccjet1index;
      vector<size_t> ccjet2index;
      set<size_t> doublecount_s1; //to prevent double counting of jet pairs in loop 
      for (size_t i = 0; i < lJets.size(); ++i) {
        double pull_12_min = Rivet::PI;
        size_t jet1index = i;
        size_t jet2index = 100;
        for(size_t j = 0; j < lJets.size(); ++j) {
          if (j != i) {
            double pull_12 = CalculatePullAngle(lJets[i], lJets[j], 0);
            if (pull_12 < pull_12_min) {
              pull_12_min = pull_12;
              jet2index = j;
            }
          }

        }
        if (jet2index != 100) { //if cc jets found
          ccjet1index.push_back(jet1index);
          ccjet2index.push_back(jet2index);
        }
      }

      for (size_t i = 0; i < ccjet1index.size(); ++i) {
        size_t j = ccjet2index[i];
        size_t k = ccjet2index[j];
        double ccjetsinvariantmass = fabs(CalculateInvariantMass(lJets[i], lJets[j]));
        if (!doublecount_s1.count(i) && !doublecount_s1.count(j)) {
          numberccjetpairs_s1 += 1;
          h_ccjetsinvariantmass_s1->fill(ccjetsinvariantmass);
          double pull_12 = CalculatePullAngle(lJets[i], lJets[j], 0);
          h_pull_12_s1->fill(pull_12 / Rivet::PI);
          double pull_21 = CalculatePullAngle(lJets[j], lJets[i], 0);
          h_pull_21_s1->fill(pull_21 / Rivet::PI);
          doublecount_s1.insert(i);
          doublecount_s1.insert(j);
        }
      }
      h_numberccjetpairs_s1->fill(numberccjetpairs_s1);

      //PULL ANGLE SELECTION 2 (tight):
      //Reconstructing W from light jets with colour connection information via pull angle 
      double numberccjetpairs_s2 = 0;
      for (size_t i = 0; i < lJets.size()-1; ++i) { 
       for (size_t j = i + 1; j < lJets.size(); ++j) {
         double pull_12 = CalculatePullAngle(lJets[i], lJets[j], 0);
         double pull_21 = CalculatePullAngle(lJets[j], lJets[i], 0);
         double ccjetsinvariantmass = fabs(CalculateInvariantMass(lJets[i], lJets[j]));
         if ((pull_12 / Rivet::PI) < 0.5 && (pull_21 / Rivet::PI) < 0.5) {
           numberccjetpairs_s2 += 1; 
           h_ccjetsinvariantmass_s2->fill(ccjetsinvariantmass);
           h_pull_12_s2->fill(pull_12 / Rivet::PI);
           h_pull_21_s2->fill(pull_21 / Rivet::PI);
         }
       }
      }
      h_numberccjetpairs_s2->fill(numberccjetpairs_s2);
      
      
      //PULL ANGLE SELECTION 3 (loose):
      //Reconstructing W from light jets with colour connection information via pull angle 
      double numberccjetpairs_s3 = 0;
      set<size_t> doublecount_s3;
      for (size_t i = 0; i < lJets.size()-1; ++i) { 
       for (size_t j = i + 1; j < lJets.size(); ++j) {
         double pull_12 = CalculatePullAngle(lJets[i], lJets[j], 0);
         double pull_21 = CalculatePullAngle(lJets[j], lJets[i], 0);
         double ccjetsinvariantmass = fabs(CalculateInvariantMass(lJets[i], lJets[j]));
         if ((pull_12 / Rivet::PI) < 0.5 || (pull_21 / Rivet::PI) < 0.5) {
           numberccjetpairs_s3 += 1; 
           h_ccjetsinvariantmass_s3->fill(ccjetsinvariantmass);
           h_pull_12_s3->fill(pull_12 / Rivet::PI);
           h_pull_21_s3->fill(pull_21 / Rivet::PI);
         }
       }
      }
      h_numberccjetpairs_s3->fill(numberccjetpairs_s3);
      

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
      normalize(h_numberwjetpairs);
      normalize(h_numberccjetpairs_s1);
      normalize(h_numberccjetpairs_s2);
      normalize(h_numberccjetpairs_s3);
      normalize(h_numberlightjets);
      normalize(h_numberbjets);
      normalize(h_numberwjets);
      normalize(h_pull_wj1wj2);
      normalize(h_Ht);
      normalize(h_ccjetsinvariantmass_s1);
      normalize(h_ccjetsinvariantmass_s2);
      normalize(h_ccjetsinvariantmass_s3);
      normalize(h_pull_12_s1);
      normalize(h_pull_21_s1);
      normalize(h_pull_12_s2);
      normalize(h_pull_21_s2);
      normalize(h_pull_12_s3);
      normalize(h_pull_21_s3);

    }

    //@}


  private:

    Histo1DPtr h_lightjetsinvariantmass;
    Histo1DPtr h_leading_lightjetsinvariantmass;
    Histo1DPtr h_numberwjetpairs;
    Histo1DPtr h_numberccjetpairs_s1;
    Histo1DPtr h_numberccjetpairs_s2;
    Histo1DPtr h_numberccjetpairs_s3;
    Histo1DPtr h_numberlightjets;
    Histo1DPtr h_numberbjets;
    Histo1DPtr h_Ht;
    Histo1DPtr h_numberwjets;
    Histo1DPtr h_pull_wj1wj2;
    Histo1DPtr h_ccjetsinvariantmass_s1;
    Histo1DPtr h_ccjetsinvariantmass_s2;
    Histo1DPtr h_ccjetsinvariantmass_s3;
    Histo1DPtr h_pull_12_s1;
    Histo1DPtr h_pull_21_s1;
    Histo1DPtr h_pull_12_s2;
    Histo1DPtr h_pull_21_s2;
    Histo1DPtr h_pull_12_s3;
    Histo1DPtr h_pull_21_s3;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ColourReconstruction);

}
