// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <iostream>
namespace Rivet {


   /// CMS Azimuthal corellations delta_phi_min2j at 13TeV
   class Azimuthal_Threejets_ak4_delta_phi_min2j : public Analysis {
  public:

    Azimuthal_Threejets_ak4_delta_phi_min2j() : Analysis("Azimuthal_Threejets_ak4_delta_phi_min2j") {}


    void init() {
      FinalState fs;
      FastJets akt(fs, FastJets::ANTIKT, 0.4);
      addProjection(akt, "antikT");

      _h_deltaPhi.addHistogram( 200.,  300., bookHisto1D(1, 1, 1));
      _h_deltaPhi.addHistogram( 300.,  400., bookHisto1D(2, 1, 1));
      _h_deltaPhi.addHistogram( 400.,  500., bookHisto1D(3, 1, 1));
      _h_deltaPhi.addHistogram( 500.,  600., bookHisto1D(4, 1, 1));
      _h_deltaPhi.addHistogram( 600.,  700., bookHisto1D(5, 1, 1));
      _h_deltaPhi.addHistogram( 700.,  800., bookHisto1D(6, 1, 1));
      _h_deltaPhi.addHistogram( 800.,  1000., bookHisto1D(7, 1, 1));
      _h_deltaPhi.addHistogram( 1000., 7000., bookHisto1D(8, 1, 1));

    }


    void analyze(const Event & event) {

      const double weight = event.weight();

      const Jets& jets = applyProjection<JetAlg>(event, "antikT").jetsByPt();
      
      if ( jets.size() < 3 ) vetoEvent;

      if ( fabs(jets[0].rap()) > 2.5 || jets[0].pT() < 200.*GeV ) vetoEvent;
      
      if ( fabs(jets[1].rap()) > 2.5 || jets[1].pT() < 100.*GeV ) vetoEvent;
      
      if ( fabs(jets[2].rap()) > 2.5 || jets[2].pT() < 100.*GeV ) vetoEvent;

      double dphi01 = deltaPhi(jets[0].phi(), jets[1].phi());

      if (dphi01 < PI/2. ) vetoEvent;

      double dphi02 = deltaPhi(jets[0].phi(), jets[2].phi());
     
      double dphi12 = deltaPhi(jets[1].phi(), jets[2].phi());
         
       /// evaluate DPhi2Jmin

       double Dphis2J[3]={dphi01,dphi02,dphi12};

       double DPhi2Jmin=Dphis2J[0];

       for(int gg=1; gg<3; ++gg){ if(DPhi2Jmin>Dphis2J[gg]){DPhi2Jmin=Dphis2J[gg];} }
      
      _h_deltaPhi.fill(jets[0].pT(), DPhi2Jmin, weight);
    }


    void finalize() {
      
      foreach (Histo1DPtr histo, _h_deltaPhi.getHistograms()) {
        normalize(histo); 
      }
    }

  private:

    BinnedHistogram<double> _h_deltaPhi;

  };

  // A hook for the plugin system
  DECLARE_RIVET_PLUGIN(Azimuthal_Threejets_ak4_delta_phi_min2j);

}
