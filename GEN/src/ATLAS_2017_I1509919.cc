// -*- C++ -*-
#include <sstream> 

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class ATLAS_2017_I1509919 : public Analysis {
  public:

   enum regionID {
      kToward  = 0,
      kAway,
      kTrans,
      kTransMax,
      kTransMin,
      kTransDiff,
      kNregions
    };

    const static int    kNptleadCuts = 3;
    const static double ptCut      [kNptleadCuts];


    ATLAS_2017_I1509919() : Analysis("ATLAS_2017_I1509919") {    
      for (int iC=0; iC < kNptleadCuts; ++iC) {
        _sumw[iC] = 0.0 ;
      }
    }


    void init() {

      const ChargedFinalState cfs500(-2.5, 2.5, 0.5*GeV);
      addProjection(cfs500, "CFS500");
      //addProjection(ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT > 500*MeV), "CFS500");

      for (int iR=0; iR < kNregions; ++iR) {
        
        // Nch profiles vs. pT_lead 
        _hist_nch[iR] = bookProfile1D(1, 1 + iR, 1); 

        // pTsum profiles vs. pT_lead 
        _hist_ptsum[iR] = bookProfile1D(2, 1 + iR, 1); 

        // <pT> profiles vs pT_lead (not measured for trans diff)
        if (iR != kTransDiff)  _hist_ptavg[iR] = bookProfile1D(3, 1 + iR, 1);
        
        // <pT> profiles vs. Nch (not measured for trans diff)
        if (iR != kTransDiff)  _hist_dn_dpt[iR] = bookProfile1D(4, 1 + iR, 1);
        //only measured for trans max/min
        if ( (iR == kTransMax) || (iR == kTransMin) )  _hist_dn_dpt2[iR] = bookProfile1D(5, 1 + iR, 1); 
        
      } //end loop over regions
      
      for (int iC=0; iC < kNptleadCuts; ++iC) {
        // Nch vs. Delta(phi) profiles 
        _hist_N_vs_dPhi[iC] = bookProfile1D(6, 1 + iC, 1);
        
        // pT vs. Delta(phi) profiles 
        _hist_pT_vs_dPhi[iC] = bookProfile1D(7, 1 + iC, 1);

        //ptLead histos only for 1 and 5 GeV cuts
        if ( (iC == 0) || (iC == 1) )  _hist_ptLead[iC] = bookHisto1D(8, 1 + iC, 1); 
      }
      
    } 

    void analyze(const Event& event) {

      const double weight = event.weight();

      // These are the charged particles (tracks) with pT > 500 MeV
      const ChargedFinalState& charged500 = applyProjection<ChargedFinalState>(event, "CFS500");
      
      MSG_DEBUG("charged500 size : " << charged500.size() );

      const Particles& particlesAll = charged500.particlesByPt();
      MSG_DEBUG("particles size : " << particlesAll.size() );
      
      //const Cut& pcut = (
      //     (Cuts::abspid!=PID::SIGMAMINUS) && (Cuts::abspid!=PID::SIGMAPLUS) &&
      //     (Cuts::abspid!=PID::XIMINUS)    && (Cuts::abspid!=PID::OMEGAMINUS));
      const Particles& particles = charged500.particlesByPt();
      MSG_DEBUG("particles size (without strange baryons): " <<  particles.size() );

      //fill pt lead histograms first
      if (particles.size()) {
        PdgId pdg = particles[0].abspid();
        if ( pdg == 3112 || // Sigma-
             pdg == 3222 || // Sigma+
             pdg == 3312 || // Xi-
             pdg == 3334 ){  // Omega-
          for (int iC = 0; iC < kNptleadCuts; ++iC) {
            if ( (iC == 0) || (iC == 1) )  {
              if (particles[0].perp()/GeV >= ptCut[iC]) {
                _sumw[iC] += weight;
                _hist_ptLead[iC]->fill( particles[0].perp()/GeV, weight);
              }
            }
          }
        }
      }

      // Require at least one track in the event with pT >= 1 GeV
      if (particles.empty())             vetoEvent;
      if (particles[0].perp()/GeV < 1.0) vetoEvent;

      // Identify leading track and its phi and pT 
      Particle p_lead      = particles[0];
      const double philead = p_lead.phi();
      const double etalead = p_lead.eta();
      const double pTlead  = p_lead.perp();
      MSG_DEBUG("Leading track: pT = " << pTlead << ", eta = " << etalead << ", phi = " << philead);

      // Iterate over all particles and count particles and scalar pTsum in three basic regions
      vector<double> num(kNregions, 0), ptSum(kNregions, 0.0), avgpt(kNregions, 0.0);

      // Temporary histos that bin Nch and pT in dPhi.
      Histo1D hist_num_dphi (*_hist_N_vs_dPhi[0].get(), "/hist_num_dphi");
      Histo1D hist_pt_dphi  (*_hist_pT_vs_dPhi[0].get(), "/hist_pt_dphi");

      hist_num_dphi.reset();
      hist_pt_dphi .reset();

      int    tmpnch[2]   = {0,0};
      double tmpptsum[2] = {0,0};

      foreach (const Particle& p, particles) {

        PdgId pdg = particles[0].abspid();
        if ( pdg == 3112 || // Sigma-
             pdg == 3222 || // Sigma+
             pdg == 3312 || // Xi-
             pdg == 3334 )  // Omega-
          continue;

        const double pT   = p.pT();
        const double dPhi = deltaPhi(philead, p.phi());  //inrange (0,pi)
        const int    ir   = region_index(dPhi); //gives just toward/away/trans
        
        //toward/away/trans region: just count 
        num  [ir] += 1;
        ptSum[ir] += pT;
        
        //determine which transverse side
        if (ir == kTrans) {
          const size_t iside = (mapAngleMPiToPi(p.phi() - philead) > 0) ? 0 : 1;
          tmpnch  [iside] += 1;
          tmpptsum[iside] += p.pT();
        }
        
        // Fill temp histos to bin Nch and pT in dPhi
        if (p.genParticle() != p_lead.genParticle()) { // We don't want to fill all those zeros from the leading track...
          hist_num_dphi.fill(dPhi * 180. / Rivet::pi, 1);
          hist_pt_dphi .fill(dPhi * 180. / Rivet::pi, pT);
        }
      }//end loop over particles

      //construct max/min/diff regions
      num[kTransMax ]   = std::max(tmpnch[0], tmpnch[1]);
      num[kTransMin ]   = std::min(tmpnch[0], tmpnch[1]);
      num[kTransDiff]   = num[kTransMax ] - num[kTransMin ];
      ptSum[kTransMax ] = std::max(tmpptsum[0], tmpptsum[1]);
      ptSum[kTransMin ] = std::min(tmpptsum[0], tmpptsum[1]);
      ptSum[kTransDiff] = ptSum[kTransMax ] - ptSum[kTransMin ];
      avgpt[kToward]    = (num[kToward] > 0 ) ? ptSum[kToward] / num[kToward] : 0. ;
      avgpt[kAway]      = (num[kAway  ] > 0 ) ? ptSum[kAway]   / num[kAway]   : 0. ;
      avgpt[kTrans]     = (num[kTrans ] > 0 ) ? ptSum[kTrans]  / num[kTrans]  : 0. ;
      //avg pt max/min regions determined according sumpt max/min
      int sumptMaxRegID = (tmpptsum[0] >  tmpptsum[1]) ? 0 : 1 ; 
      int sumptMinRegID = (sumptMaxRegID == 0) ? 1 : 0; 
      avgpt[kTransMax ] = (tmpnch[sumptMaxRegID] > 0) ? tmpptsum[sumptMaxRegID] / tmpnch[sumptMaxRegID] : 0.;
      avgpt[kTransMin ] = (tmpnch[sumptMinRegID] > 0) ? tmpptsum[sumptMinRegID] / tmpnch[sumptMinRegID] : 0.;
      avgpt[kTransDiff] = ((tmpnch[sumptMaxRegID] > 0) && (tmpnch[sumptMinRegID] > 0)) ? avgpt[kTransMax ] - avgpt[kTransMin ] : 0.;
        
      // Now fill underlying event histograms

      // The densities are calculated by dividing the UE properties by dEta*dPhi
      // -- each basic region has a dPhi of 2*PI/3 and dEta is two times 2.5
      // min/max/diff regions are only half of that
      const double dEtadPhi[kNregions] = { 2*2.5 * 2*PI/3.0, 2*2.5 * 2*PI/3.0, 2*2.5 * 2*PI/3.0, 
                                           2*2.5 *   PI/3.0, 2*2.5 *   PI/3.0, 2*2.5 *   PI/3.0 };
      
      for (int iR=0; iR < kNregions; ++iR) {
        
        _hist_nch  [iR]->fill(pTlead/GeV, num[iR]   /dEtadPhi[iR]     , weight);
        _hist_ptsum[iR]->fill(pTlead/GeV, ptSum[iR] /GeV/dEtadPhi[iR] , weight);
        
        // <pT> profiles vs. pT_lead
        switch (iR) {
          //first 3 are the same!!!
          case kToward    : 
          case kAway      : 
          case kTrans     : 
            if (num[iR] > 0) {
              _hist_ptavg[iR]->fill(pTlead/GeV, avgpt[iR] /GeV, weight); 
            }
                  break;
          case kTransMax  : 
            if (tmpnch[sumptMaxRegID] > 0) {
              _hist_ptavg[iR]->fill(pTlead/GeV, avgpt[iR] /GeV, weight); 
            }
            break;
          case kTransMin  : 
            if (tmpnch[sumptMinRegID] > 0) {
              _hist_ptavg[iR]->fill(pTlead/GeV, avgpt[iR] /GeV, weight);
            }
            break;
          case kTransDiff : 
            break;
          default : //should not get here!!!
            MSG_INFO("unknown region in <pT> profiles vs.pt lead switch!!! : " <<  iR);
        } //end switch over regions
        
        
        // <pT> profiles vs. Nch
        switch (iR) {
          //first 3 are the same!!!
          case kToward    : 
          case kAway      : 
          case kTrans     : 
            if (num[iR] > 0) {
              _hist_dn_dpt[iR]->fill(num[iR] , avgpt[iR] /GeV, weight); 
            }
            break;
          case kTransMax  : 
            if (tmpnch[sumptMaxRegID] > 0) {
              _hist_dn_dpt [iR]->fill(num[kTrans]          , avgpt[iR] /GeV, weight); 
              _hist_dn_dpt2[iR]->fill(tmpnch[sumptMaxRegID], avgpt[iR] /GeV, weight); 
            }
            break;
          case kTransMin  : 
            if (tmpnch[sumptMinRegID] > 0) {
              _hist_dn_dpt [iR]->fill(num[kTrans]          , avgpt[iR] /GeV, weight);
              _hist_dn_dpt2[iR]->fill(tmpnch[sumptMinRegID], avgpt[iR] /GeV, weight);
            }
            break;
          case kTransDiff : 
            break;
          default : //should not get here!!!
            MSG_INFO("unknown region in <pT> profiles vs. nch switch!!! : " <<  iR);
        } //end switch over regions
        
      } //end loop over regions
      
      
      // Update the "proper" dphi profile histograms
      // Note that we fill dN/dEtadPhi: dEta = 2*2.5, dPhi = 2*PI/nBins
      //Rivet::pi / 180. --> to convert normalization to units of radians instead of degrees
      const double dEtadPhi2 = (2*2.5 * 2) * Rivet::pi / 180.;
      
      for (unsigned i = 0; i < hist_num_dphi.numBins() ; ++i) {  
        
        // First Nch
        double mean = hist_num_dphi.bin(i).xMid() ;
        double value = 0.;
        if (hist_num_dphi.bin(i).numEntries() > 0) {
          mean  = hist_num_dphi.bin(i).xMean() ;
          value = hist_num_dphi.bin(i).area()/hist_num_dphi.bin(i).xWidth()/dEtadPhi2;
        }
        for (int iC=0; iC < kNptleadCuts; ++iC) {
          if (pTlead/GeV >= ptCut[iC]) _hist_N_vs_dPhi[iC] ->fill(mean, value, weight);
        }
        
        // Then pT
        mean = hist_pt_dphi.bin(i).xMid() ;
        value = 0.;
        if (hist_pt_dphi.bin(i).numEntries() > 0) {
          mean  = hist_pt_dphi.bin(i).xMean() ;
          value = hist_pt_dphi.bin(i).area()/hist_pt_dphi.bin(i).xWidth()/dEtadPhi2;
        }
        for (int iC=0; iC < kNptleadCuts; ++iC) {
          if (pTlead/GeV >= ptCut[iC]) _hist_pT_vs_dPhi[iC] ->fill(mean, value, weight);
        }
      } //end loop over bins
            
    } //end analyze()
          
          
    void finalize() {

      for (int iC=0; iC < kNptleadCuts; ++iC) {
        if ( (iC == 0) || (iC == 1) )  scale(_hist_ptLead[iC], 1.0/_sumw[iC]);
      }
    }


  private:

    // Little helper function to identify basic Delta(phi) regions: toward/away/trans
    inline int region_index(double dphi) {
      assert(inRange(dphi, 0.0, PI, CLOSED, CLOSED));
      if (dphi < PI/3.0)    return kToward;
      if (dphi < 2*PI/3.0)  return kTrans;
      return kAway;
    }

  private:

    //nch, sumpt, avgpt profiles vs. pt lead
    Profile1DPtr _hist_nch   [kNregions]; //for regions: all 6 regions  
    Profile1DPtr _hist_ptsum [kNregions]; //for regions: all 6 regions  
    Profile1DPtr _hist_ptavg [kNregions]; //for regions: trans towards/away/all/min/max 

    //nch, sumpt, avgpt profiles vs. nch
    Profile1DPtr _hist_dn_dpt       [kNregions]; //regions: towards/away/ vs nch(region) & trans all/min/max vs nch(trans) 
    Profile1DPtr _hist_dn_dpt2      [kNregions]; //regions: trans min/max vs. nch(region) 

    Profile1DPtr _hist_N_vs_dPhi [kNptleadCuts]; 
    Profile1DPtr _hist_pT_vs_dPhi[kNptleadCuts]; 

    Histo1DPtr   _hist_ptLead[kNptleadCuts]; //for 1,5 GeV cuts only 

    double _sumw[kNptleadCuts];

  };

  //constants
  const double ATLAS_2017_I1509919::ptCut[] = { 1., 5., 10.};
 
  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1509919);
}
