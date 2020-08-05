#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/PartonicTops.hh"


namespace Rivet {

  class CMS_2015_I1370682_parton : public Analysis {
    public:
      CMS_2015_I1370682_parton() : Analysis("CMS_2015_I1370682_parton") {}

      void init() {
        // Parton level top quarks
        declare(PartonicTops(PartonicTops::E_MU, false), "LeptonicPartonTops");
        declare(PartonicTops(PartonicTops::HADRONIC),    "HadronicPartonTops");

        _hSL_topPt         = bookHisto1D("d15-x01-y01");
        _hSL_topPtTtbarSys = bookHisto1D("d16-x01-y01");
        _hSL_topY          = bookHisto1D("d17-x01-y01");
        _hSL_ttbarDelPhi   = bookHisto1D("d18-x01-y01");
        _hSL_topPtLead     = bookHisto1D("d19-x01-y01");
        _hSL_topPtSubLead  = bookHisto1D("d20-x01-y01");
        _hSL_ttbarPt       = bookHisto1D("d21-x01-y01");
        _hSL_ttbarY        = bookHisto1D("d22-x01-y01");
        _hSL_ttbarMass     = bookHisto1D("d23-x01-y01");

        _hDL_topPt           = bookHisto1D("d24-x01-y01");
        _hDL_topPtTtbarSys   = bookHisto1D("d25-x01-y01");
        _hDL_topY            = bookHisto1D("d26-x01-y01");
        _hDL_ttbarDelPhi     = bookHisto1D("d27-x01-y01");
        _hDL_topPtLead       = bookHisto1D("d28-x01-y01");
        _hDL_topPtSubLead    = bookHisto1D("d29-x01-y01");
        _hDL_ttbarPt         = bookHisto1D("d30-x01-y01");
        _hDL_ttbarY          = bookHisto1D("d31-x01-y01");
        _hDL_ttbarMass       = bookHisto1D("d32-x01-y01");

      };

      void analyze(const Event& event) {
        const double weight = event.weight();
        
        bool isDilepton   = false;
        bool isSemilepton = false;
        
        // Do the analysis only for the ttbar full leptonic or semileptonic channel, without tau decay
        const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
        const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt();
        if (leptonicpartontops.size() == 1 and hadronicpartontops.size() == 1)
          isSemilepton = true;
        else if (leptonicpartontops.size() == 2)
          isDilepton = true;
        else
          vetoEvent;

        // Parton level at full phase space
        // Fill top quarks defined in the parton level, full phase space
        FourMomentum t1P4;
        FourMomentum t2P4;
        if (isSemilepton) {
          t1P4 = leptonicpartontops[0].momentum();
          t2P4 = hadronicpartontops[0].momentum();
        }
        else if (isDilepton) {
          t1P4 = leptonicpartontops[0].momentum();
          t2P4 = leptonicpartontops[1].momentum();
        }
        double t1Pt = t1P4.pT(), t2Pt = t2P4.pT();
        FourMomentum ttbarP4 = t1P4+t2P4;
        FourMomentum t1P4AtCM = LorentzTransform::mkFrameTransformFromBeta(ttbarP4.betaVec()).transform(t1P4);
        double dPhi = deltaPhi(t1P4.phi(), t2P4.phi());

        if (isSemilepton) {
          _hSL_topPt->fill(t1Pt, weight);
          _hSL_topPt->fill(t2Pt, weight);
          _hSL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
          _hSL_topY->fill(t1P4.rapidity(), weight);
          _hSL_topY->fill(t2P4.rapidity(), weight);
          _hSL_ttbarDelPhi->fill(dPhi, weight);
          _hSL_topPtLead->fill(std::max(t1Pt, t2Pt), weight);
          _hSL_topPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
          _hSL_ttbarPt->fill(ttbarP4.pT(), weight);
          _hSL_ttbarY->fill(ttbarP4.rapidity(), weight);
          _hSL_ttbarMass->fill(ttbarP4.mass(), weight);
        }
        else if (isDilepton) {
          _hDL_topPt->fill(t1Pt, weight);
          _hDL_topPt->fill(t2Pt, weight);
          _hDL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
          _hDL_topY->fill(t1P4.rapidity(), weight);
          _hDL_topY->fill(t2P4.rapidity(), weight);
          _hDL_ttbarDelPhi->fill(dPhi, weight);
          _hDL_topPtLead->fill(std::max(t1Pt, t2Pt), weight);
          _hDL_topPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
          _hDL_ttbarPt->fill(ttbarP4.pT(), weight);
          _hDL_ttbarY->fill(ttbarP4.rapidity(), weight);
          _hDL_ttbarMass->fill(ttbarP4.mass(), weight);
        }
      };

      void finalize() {
        normalize(_hSL_topPt        );
        normalize(_hSL_topPtTtbarSys);
        normalize(_hSL_topY         );
        normalize(_hSL_ttbarDelPhi  );
        normalize(_hSL_topPtLead    );
        normalize(_hSL_topPtSubLead );
        normalize(_hSL_ttbarPt      );
        normalize(_hSL_ttbarY       );
        normalize(_hSL_ttbarMass    );

        normalize(_hDL_topPt        );
        normalize(_hDL_topPtTtbarSys);
        normalize(_hDL_topY         );
        normalize(_hDL_ttbarDelPhi  );
        normalize(_hDL_topPtLead    );
        normalize(_hDL_topPtSubLead );
        normalize(_hDL_ttbarPt      );
        normalize(_hDL_ttbarY       );
        normalize(_hDL_ttbarMass    );

      };

    private:
      Histo1DPtr _hSL_topPt        ;
      Histo1DPtr _hSL_topPtTtbarSys;
      Histo1DPtr _hSL_topY         ;
      Histo1DPtr _hSL_ttbarDelPhi  ;
      Histo1DPtr _hSL_topPtLead    ;
      Histo1DPtr _hSL_topPtSubLead ;
      Histo1DPtr _hSL_ttbarPt      ;
      Histo1DPtr _hSL_ttbarY       ;
      Histo1DPtr _hSL_ttbarMass    ;

      Histo1DPtr _hDL_topPt        ;
      Histo1DPtr _hDL_topPtTtbarSys;
      Histo1DPtr _hDL_topY         ;
      Histo1DPtr _hDL_ttbarDelPhi  ;
      Histo1DPtr _hDL_topPtLead    ;
      Histo1DPtr _hDL_topPtSubLead ;
      Histo1DPtr _hDL_ttbarPt      ;
      Histo1DPtr _hDL_ttbarY       ;
      Histo1DPtr _hDL_ttbarMass    ;

  };

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CMS_2015_I1370682_parton> plugin_CMS_2015_I1370682_parton;

}
