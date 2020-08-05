// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {

  class CMS_FSQ_15_008 : public Analysis {
  public:

    /// Constructor
    CMS_FSQ_15_008()
      : Analysis("CMS_FSQ_15_008")
    {    }


    /// @name Analysis methods

    /// Book histograms and initialise projections before the run
    void init() {

        /// @todo Initialise and register projections here
        const FinalState fsa(-5.0, 5.0, 0.0000*GeV);
        addProjection(fsa, "FSA");
        const ChargedFinalState cfs(-2.4, 2.4, 0.5*GeV); //2.5 -> 2.4
        addProjection(cfs, "CFS");

        _Nevt_noCuts  = 0;
        _Nevt_incl    = 0;
        _Nevt_after_cuts_or   = 0;
        _Nevt_after_cuts_and  = 0;
        _Nevt_after_cuts_xor  = 0;
        _Nevt_after_cuts_xorm = 0; 
        _Nevt_after_cuts_xorp = 0; 
        
        /// @todo Book histograms here
        _h_dNch_dEta_incl       = bookHisto1D("dn_deta_incl",30,-3,3);  //40 -> 30 bins
        _h_dNch_all_dEta_incl   = bookHisto1D(1,1,1);
        _h_dNch_all_dEta_OR     = bookHisto1D(1,2,1);
        _h_dNch_all_dEta_AND    = bookHisto1D(1,3,1);
        _h_dNch_all_dEta_XOR    = bookHisto1D(1,4,1);
        _h_dNch_all_dEta_XORm   = bookHisto1D(1,5,1);
        _h_dNch_all_dEta_XORp   = bookHisto1D(1,6,1);
        _h_dNch_all_dEta_XORpm  = bookHisto1D(1,7,1);

        _h_dNch_all_dpt_incl   = bookHisto1D(2,1,1);
        _h_dNch_all_dpt_OR     = bookHisto1D(2,2,1);
        _h_dNch_all_dpt_AND    = bookHisto1D(2,3,1);
        _h_dNch_all_dpt_XOR    = bookHisto1D(2,4,1);
        _h_dNch_all_dpt_XORm   = bookHisto1D(2,5,1);
        _h_dNch_all_dpt_XORp   = bookHisto1D(2,6,1);

        _h_dNch_leading_dpt_incl   = bookHisto1D(3,1,1);
        _h_dNch_leading_dpt_OR     = bookHisto1D(3,2,1);
        _h_dNch_leading_dpt_AND    = bookHisto1D(3,3,1);
        _h_dNch_leading_dpt_XOR    = bookHisto1D(3,4,1);
        _h_dNch_leading_dpt_XORm   = bookHisto1D(3,5,1);
        _h_dNch_leading_dpt_XORp   = bookHisto1D(3,6,1);

        _h_integrated_leading_pt_incl   = bookHisto1D(4,1,1);
        _h_integrated_leading_pt_OR     = bookHisto1D(4,2,1);
        _h_integrated_leading_pt_AND    = bookHisto1D(4,3,1);
        _h_integrated_leading_pt_XOR    = bookHisto1D(4,4,1);
        _h_integrated_leading_pt_XORm   = bookHisto1D(4,5,1);
        _h_integrated_leading_pt_XORp   = bookHisto1D(4,6,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const ChargedFinalState& charged = applyProjection<ChargedFinalState>(event, "CFS");
      const FinalState& fsa = applyProjection<FinalState>(event, "FSA");
      int count_plus = 0, count_minus = 0;
      
      double Ecut = 5.;
      foreach (const Particle& p, fsa.particles()) {
        if(p.energy() >= Ecut ){ 
          if (inRange(p.eta(),  3.0,  5.0)) count_plus++;
          if (inRange(p.eta(), -5.0, -3.0)) count_minus++;
        }
      }

      // Cut combinations
      const bool cutsor   = (  count_plus  > 0 || count_minus  > 0);
      const bool cutsand  = (  count_plus  > 0 && count_minus  > 0);
      const bool cutsxor  = ( (count_plus  > 0 && count_minus == 0) || (count_plus == 0 && count_minus > 0) );
      const bool cutsxorm = ( (count_plus == 0 && count_minus  > 0) ); 
      const bool cutsxorp = ( (count_plus  > 0 && count_minus == 0) ); 

      // Increment counters and fill histos
      _Nevt_noCuts += weight;
      //find pttrackmax
      double track_ptmax = -99.9;
      
      if ( charged.size() >= 1 ) { 
          _Nevt_incl += weight;
          if (cutsor)   _Nevt_after_cuts_or   += weight;
          if (cutsand)  _Nevt_after_cuts_and  += weight;
          if (cutsxor)  _Nevt_after_cuts_xor  += weight; 
          if (cutsxorm) _Nevt_after_cuts_xorm += weight; 
          if (cutsxorp) _Nevt_after_cuts_xorp += weight; 
          // cout << " norm factors " << _Nevt_incl << " " << _Nevt_after_cuts_or << " " << _Nevt_after_cuts_and << " " << _Nevt_after_cuts_xor << endl;
      }
     
      //Loop over partcles
      foreach (const Particle& p, charged.particles()) {
          int id = p.pid();
          
          // continue if particle is a proton, a kaon or a pion
          if (( (abs(id) == 211) || (abs(id) == 321) || (abs(id) == 2212)))  {
              _h_dNch_dEta_incl ->fill(p.eta(), weight); 
          }

          _h_dNch_all_dEta_incl -> fill(p.eta(), weight); 
          if (cutsor)   _h_dNch_all_dEta_OR   -> fill(p.eta(), weight);
          if (cutsand)  _h_dNch_all_dEta_AND  -> fill(p.eta(), weight);
          if (cutsxor)  _h_dNch_all_dEta_XOR  -> fill(p.eta(), weight);
          if (cutsxorm) _h_dNch_all_dEta_XORm -> fill(p.eta(), weight); 
          if (cutsxorp) _h_dNch_all_dEta_XORp -> fill(p.eta(), weight); 
          //Average xorm & xorp
          if (cutsxorm) _h_dNch_all_dEta_XORpm -> fill(p.eta(), weight);   //It has to be normilized to _Nevt_after_cuts_xorm + _Nevt_after_cuts_xorp
          if (cutsxorp) _h_dNch_all_dEta_XORpm -> fill(-1*p.eta(), weight);

          _h_dNch_all_dpt_incl -> fill(p.pT(), weight); 
          if (cutsor)   _h_dNch_all_dpt_OR   -> fill(p.pT(), weight);
          if (cutsand)  _h_dNch_all_dpt_AND  -> fill(p.pT(), weight);
          if (cutsxor)  _h_dNch_all_dpt_XOR  -> fill(p.pT(), weight);
          if (cutsxorm) _h_dNch_all_dpt_XORm -> fill(p.pT(), weight); 
          if (cutsxorp) _h_dNch_all_dpt_XORp -> fill(p.pT(), weight); 

          // Find the leading track of the event
          if (inRange(p.eta(),  -2.4,  2.4)) { 
              if (p.pT() > track_ptmax) track_ptmax = p.pT();
          } 
      }//End loop over particles
    

      //Filling pt leading distributions
                    _h_dNch_leading_dpt_incl -> fill(track_ptmax, weight); 
      if (cutsor)   _h_dNch_leading_dpt_OR   -> fill(track_ptmax, weight);
      if (cutsand)  _h_dNch_leading_dpt_AND  -> fill(track_ptmax, weight);
      if (cutsxor)  _h_dNch_leading_dpt_XOR  -> fill(track_ptmax, weight);
      if (cutsxorm) _h_dNch_leading_dpt_XORm -> fill(track_ptmax, weight); 
      if (cutsxorp) _h_dNch_leading_dpt_XORp -> fill(track_ptmax, weight); 


      //Integration of pt leading
      for (size_t i = 0; i < _h_integrated_leading_pt_incl->numBins(); ++i) {

        double binlimitlow_t = _h_integrated_leading_pt_incl->bin(i).xMin(),
               weightbw_t    = _h_integrated_leading_pt_incl->bin(i).xWidth()*weight,
               xbin_t        = _h_integrated_leading_pt_incl->bin(i).xMid();
        // cout << " track bin " << i << " " << binlimitlow_t << " " <<_h_tracks->bin(i).xWidth() << " " << xbin_t << endl;
        
        if (track_ptmax > binlimitlow_t) { 
                          _h_integrated_leading_pt_incl -> fill(xbin_t, weightbw_t);
            if (cutsor)   _h_integrated_leading_pt_OR   -> fill(xbin_t, weightbw_t);
            if (cutsand)  _h_integrated_leading_pt_AND  -> fill(xbin_t, weightbw_t); 
            if (cutsxor)  _h_integrated_leading_pt_XOR  -> fill(xbin_t, weightbw_t); 
            if (cutsxorm) _h_integrated_leading_pt_XORm -> fill(xbin_t, weightbw_t);  
            if (cutsxorp) _h_integrated_leading_pt_XORp -> fill(xbin_t, weightbw_t);  
        }   
      }
    }//End Analyze


    /// Normalise histograms etc., after the run
    void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here

      // scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
      // normalize(_h_YYYY); // normalize to unity
      
      cout << " norm factors: All =  " << _Nevt_noCuts <<", Inclusive = "<<_Nevt_incl << ", Inelastic = " << _Nevt_after_cuts_or << ", NSD = " << _Nevt_after_cuts_and << ", Xor = " << _Nevt_after_cuts_xor << ", Xorm = " << _Nevt_after_cuts_xorm<< ", Xorp = " << _Nevt_after_cuts_xorp << endl; 
      
      scale(_h_dNch_dEta_incl,      1./_Nevt_noCuts);
      scale(_h_dNch_all_dEta_incl,  1./_Nevt_incl);
      scale(_h_dNch_all_dEta_OR,    1./_Nevt_after_cuts_or);
      scale(_h_dNch_all_dEta_AND,   1./_Nevt_after_cuts_and);
      scale(_h_dNch_all_dEta_XOR,   1./_Nevt_after_cuts_xor);
      scale(_h_dNch_all_dEta_XORm,  1./_Nevt_after_cuts_xorm); 
      scale(_h_dNch_all_dEta_XORp,  1./_Nevt_after_cuts_xorp); 
      scale(_h_dNch_all_dEta_XORpm, 1./(_Nevt_after_cuts_xorm + _Nevt_after_cuts_xorp)); 

      scale(_h_dNch_all_dpt_incl, 1./_Nevt_incl);
      scale(_h_dNch_all_dpt_OR,   1./_Nevt_after_cuts_or);
      scale(_h_dNch_all_dpt_AND,  1./_Nevt_after_cuts_and);
      scale(_h_dNch_all_dpt_XOR,  1./_Nevt_after_cuts_xor);
      scale(_h_dNch_all_dpt_XORm, 1./_Nevt_after_cuts_xorm); 
      scale(_h_dNch_all_dpt_XORp, 1./_Nevt_after_cuts_xorp); 
      
      scale(_h_dNch_leading_dpt_incl, 1./_Nevt_incl);
      scale(_h_dNch_leading_dpt_OR,   1./_Nevt_after_cuts_or);
      scale(_h_dNch_leading_dpt_AND,  1./_Nevt_after_cuts_and);
      scale(_h_dNch_leading_dpt_XOR,  1./_Nevt_after_cuts_xor);
      scale(_h_dNch_leading_dpt_XORm, 1./_Nevt_after_cuts_xorm); 
      scale(_h_dNch_leading_dpt_XORp, 1./_Nevt_after_cuts_xorp); 
          
      scale(_h_integrated_leading_pt_incl, 1./_Nevt_incl);
      scale(_h_integrated_leading_pt_OR,   1./_Nevt_after_cuts_or);
      scale(_h_integrated_leading_pt_AND,  1./_Nevt_after_cuts_and);
      scale(_h_integrated_leading_pt_XOR,  1./_Nevt_after_cuts_xor);
      scale(_h_integrated_leading_pt_XORm, 1./_Nevt_after_cuts_xorm);
      scale(_h_integrated_leading_pt_XORp, 1./_Nevt_after_cuts_xorp);
     
      //double norm_t0 =  _h_integrated_leading_pt_incl->bin(7).height()/2.056170e-03; // ???
      ////double norm_t1 =  _h_tracks_incl->bin(7).sumW()/2.056170e-03;
      //if (norm_t0 > 0 ) {
      //    scale(_h_integrated_leading_pt_incl, 1./ norm_t0);
      //    scale(_h_integrated_leading_pt_OR,   1./ norm_t0);
      //    scale(_h_integrated_leading_pt_AND,  1./ norm_t0);
      //    scale(_h_integrated_leading_pt_XOR,  1./ norm_t0);
      //    scale(_h_integrated_leading_pt_XORm, 1./ norm_t0); 
      //    scale(_h_integrated_leading_pt_XORp, 1./ norm_t0); 
      //}
    }//End Finalize


  private:
    // Data members like post-cuts event weight counters go here
    Histo1DPtr
        _h_dNch_dEta_incl,
        _h_dNch_all_dEta_incl,  
        _h_dNch_all_dEta_AND,   
        _h_dNch_all_dEta_OR,    
        _h_dNch_all_dEta_XOR,   
        _h_dNch_all_dEta_XORm,  
        _h_dNch_all_dEta_XORp, 
        _h_dNch_all_dEta_XORpm; 

    Histo1DPtr                                                              
        _h_dNch_all_dpt_incl,
        _h_dNch_all_dpt_AND, 
        _h_dNch_all_dpt_OR,  
        _h_dNch_all_dpt_XOR,   
        _h_dNch_all_dpt_XORm,  
        _h_dNch_all_dpt_XORp;

    Histo1DPtr                                     
        _h_dNch_leading_dpt_incl,
        _h_dNch_leading_dpt_AND, 
        _h_dNch_leading_dpt_OR, 
        _h_dNch_leading_dpt_XOR,  
        _h_dNch_leading_dpt_XORm, 
        _h_dNch_leading_dpt_XORp;

    Histo1DPtr                                     
        _h_integrated_leading_pt_incl,
        _h_integrated_leading_pt_AND,
        _h_integrated_leading_pt_OR, 
        _h_integrated_leading_pt_XOR,  
        _h_integrated_leading_pt_XORm, 
        _h_integrated_leading_pt_XORp;


    double _Nevt_noCuts,
           _Nevt_incl,
           _Nevt_after_cuts_and,
           _Nevt_after_cuts_or,
           _Nevt_after_cuts_xor,
           _Nevt_after_cuts_xorp,
           _Nevt_after_cuts_xorm;


    /// @name Histograms
    //@{
    Profile1DPtr _h_XXXX;
    Histo1DPtr _h_YYYY;
    //@}


  };
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_FSQ_15_008);
}
