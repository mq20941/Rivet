// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// @brief Using colour information to identify jets in all-hadronic top decay at 13 TeV
  class MC_ColourReconstruction_allhad : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_ColourReconstruction_allhad);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections of particles/jets before the run
    void init() { 

      // A FinalState is used to select particles within |eta| < 2.5 and with pT
      // > 20 GeV, out of which the ChargedLeptons projection picks only the
      // electrons and muons, to be accessed later as "LFS".
      ChargedLeptons lfs(FinalState(Cuts::abseta < 2.5 && Cuts::pT > 20*GeV));
      declare(lfs, "LFS");

      // A second FinalState is used to select all particles in |eta| < 2.5,
      // with no pT cut. This is used to construct jets and measure missing
      // transverse energy.
      VetoedFinalState fs(FinalState(Cuts::abseta < 2.5));
      fs.addVetoOnThisFinalState(lfs);
      fs.vetoNeutrinos();
      FastJets fjets(fs, FastJets::ANTIKT, 0.4);
      declare(fjets, "Jets");
      declare(MissingMomentum(fs), "MissingET");


      //book histograms
      book(h_lightjetsinvariantmass, "lightjetsinvariantmass", 20, 0, 300);
      book(h_leading_lightjetsinvariantmass, "leading_lightjetsinvariantmass", 20, 0, 300);
      book(h_numberwjetpairs, "numberwjetpairs", 8, -0.5, 7.5);
      book(h_numberlightjets, "numberlightjets", 6, 3.5, 9.5);
      book(h_numberccjetpairs_s1, "numberccjetpairs_s1", 6, -0.5, 5.5);
      book(h_numberccjetpairs_s2, "numberccjetpairs_s2", 6, -0.5, 5.5);
      book(h_numberccjetpairs_s3, "numberccjetpairs_s3", 6, -0.5, 5.5);
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


      //event baseline cuts


      const ChargedLeptons& lfs = apply<ChargedLeptons>(event, "LFS");
      size_t nLeps = lfs.chargedLeptons().size();
      if (nLeps != 0) vetoEvent; // all-hadronic

      // Use a missing ET cut to bias against events with a hard neutrino from
      // the leptonically decaying W. This helps to reduce leptonic backgrounds.
      const Vector3& met = apply<MissingMomentum>(event, "MissingET").vectorMissingPt();
      if (met.mod() > 60*GeV) vetoEvent;

      // "Jets" projection to get jets with pT > 25Gev
      // cut on jet multiplicity 
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25.0*GeV && Cuts::absrap < 2.5);
      if (jets.size() < 6)  vetoEvent; // all-hadronic
      for (int i = 0; i<5; ++i) {
        if (jets[i].pT() < 80*GeV) vetoEvent;
      }
      
      //for (size_t i = 0; i < jets.size()-1; ++i) { 
       //for (size_t j = i + 1; j < jets.size(); ++j) {
         //if (deltaR(jets[i], jets[j]) < 0.6) vetoEvent;
       //}
      //}

      // Sort the jets into b jets and light jets
      Jets bjets, ljets;
      for (const Jet& jet : jets) {
        if (jet.bTagged())  bjets += jet;
        else                ljets += jet;
      }

      if (bjets.size() != 2) vetoEvent;
      if (ljets.size() < 4)  vetoEvent;
      h_numberlightjets->fill(ljets.size());

      if (deltaPhi(bjets[0].pTvec(), bjets[1].pTvec()) < 1.5) vetoEvent;




      //analysis on jets 

      double numberwjetpairs = 0;
      for (size_t i = 0; i < ljets.size()-1; ++i) { 
       for (size_t j = i + 1; j < ljets.size(); ++j) {
         double lightjetsinvariantmass = fabs(CalculateInvariantMass(ljets[i], ljets[j]));
         h_lightjetsinvariantmass->fill(lightjetsinvariantmass);
         if (lightjetsinvariantmass > 60*GeV && lightjetsinvariantmass < 100*GeV) { 
           numberwjetpairs += 1;
         }
       }
      }

      h_numberwjetpairs->fill(numberwjetpairs);

      double leading_lightjetsinvariantmass = fabs(CalculateInvariantMass(ljets[0], ljets[1]));
      h_leading_lightjetsinvariantmass->fill(leading_lightjetsinvariantmass);


      //PULL ANGLE SELECTION 1:
      //Finding colour connected jets from light jets (non-b)
      //Must be colour connected in both directions: if jet a cc with jet b, then jet b must be cc with jet a
      double numberccjetpairs_s1 = 0;
      vector<size_t> ccjet1index;
      vector<size_t> ccjet2index;
      set<size_t> doublecount_s1; //to prevent double counting of jet pairs in loop 
      for (size_t i = 0; i < ljets.size(); ++i) {
        vector<double> pull_values = {};
        double pull_12_min = Rivet::PI;
        size_t jet1index = i;
        size_t jet2index = 100;
        for(size_t j = 0; j < ljets.size(); ++j) {
          if (j != i) {
            double pull_12 = CalculatePullAngle(ljets[i], ljets[j], 0);
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
        double ccjetsinvariantmass = fabs(CalculateInvariantMass(ljets[i], ljets[j]));
        if (!doublecount_s1.count(i) && !doublecount_s1.count(j)) {
          numberccjetpairs_s1 += 1;
          h_ccjetsinvariantmass_s1->fill(ccjetsinvariantmass);
          double pull_12 = CalculatePullAngle(ljets[i], ljets[j], 0);
          h_pull_12_s1->fill(pull_12 / Rivet::PI);
          double pull_21 = CalculatePullAngle(ljets[j], ljets[i], 0);
          h_pull_21_s1->fill(pull_21 / Rivet::PI);
          doublecount_s1.insert(i);
          doublecount_s1.insert(j);
        }
      }
      h_numberccjetpairs_s1->fill(numberccjetpairs_s1);


      //PULL ANGLE SELECTION 2 (tight):
      //Reconstructing W from light jets with colour connection information via pull angle 
      double numberccjetpairs_s2 = 0;
      set<size_t> doublecount_s2;
      for (size_t i = 0; i < ljets.size()-1; ++i) { 
       for (size_t j = i + 1; j < ljets.size(); ++j) {
         double pull_12 = CalculatePullAngle(ljets[i], ljets[j], 0);
         double pull_21 = CalculatePullAngle(ljets[j], ljets[i], 0);
         double ccjetsinvariantmass = fabs(CalculateInvariantMass(ljets[i], ljets[j]));
         if ((pull_12 / Rivet::PI) < 0.5 && (pull_21 / Rivet::PI) < 0.5 && !doublecount_s2.count(i) && !doublecount_s2.count(j)) {
           numberccjetpairs_s2 += 1; 
           h_ccjetsinvariantmass_s2->fill(ccjetsinvariantmass);
           h_pull_12_s2->fill(pull_12 / Rivet::PI);
           h_pull_21_s2->fill(pull_21 / Rivet::PI);
           doublecount_s2.insert(i);
           doublecount_s2.insert(j);
         }
       }
      }
      h_numberccjetpairs_s2->fill(numberccjetpairs_s2);
      
      
      //PULL ANGLE SELECTION 3 (loose):
      //Reconstructing W from light jets with colour connection information via pull angle 
      double numberccjetpairs_s3 = 0;
      set<size_t> doublecount_s3;
      for (size_t i = 0; i < ljets.size()-1; ++i) { 
       for (size_t j = i + 1; j < ljets.size(); ++j) {
         double pull_12 = CalculatePullAngle(ljets[i], ljets[j], 0);
         double pull_21 = CalculatePullAngle(ljets[j], ljets[i], 0);
         double ccjetsinvariantmass = fabs(CalculateInvariantMass(ljets[i], ljets[j]));
         if (((pull_12 / Rivet::PI) < 0.1 || (pull_21 / Rivet::PI) < 0.1) && !doublecount_s3.count(i) && !doublecount_s3.count(j)) {
           numberccjetpairs_s3 += 1; 
           h_ccjetsinvariantmass_s3->fill(ccjetsinvariantmass);
           h_pull_12_s3->fill(pull_12 / Rivet::PI);
           h_pull_21_s3->fill(pull_21 / Rivet::PI);
           doublecount_s3.insert(i);
           doublecount_s3.insert(j);
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
      normalize(h_numberlightjets);
      normalize(h_numberwjetpairs);
      normalize(h_numberccjetpairs_s1);
      normalize(h_numberccjetpairs_s2);
      normalize(h_numberccjetpairs_s3);
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
    Histo1DPtr h_numberlightjets;
    Histo1DPtr h_numberwjetpairs;
    Histo1DPtr h_numberccjetpairs_s1;
    Histo1DPtr h_numberccjetpairs_s2;
    Histo1DPtr h_numberccjetpairs_s3;
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
  DECLARE_RIVET_PLUGIN(MC_ColourReconstruction_allhad);

}
