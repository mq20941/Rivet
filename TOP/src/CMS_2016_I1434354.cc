// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include "HepMC/GenParticle.h"

#include <vector>

namespace Rivet
{


  class CMS_2016_I1434354 : public Analysis
  {
    private:
      Particles _muelphs;
      Particles _neutrinos;
      Particles _jetparticles;
      Jets _dressLeps;
      vector<const Jet*> _bJets;
      const Jet*_bl;
      const Jet* _bh;
      const Jet* _wja;
      const Jet* _wjb;
      FourMomentum _tl;
      FourMomentum _th;
      FourMomentum _wl;
      FourMomentum _wh;
      FourMomentum _nusum;

      Histo1DPtr _hist_thadpt;
      Histo1DPtr _hist_thady;
      Histo1DPtr _hist_tleppt;
      Histo1DPtr _hist_tlepy;
      Histo1DPtr _hist_ttpt;
      Histo1DPtr _hist_tty;
      Histo1DPtr _hist_ttm;
      Histo1DPtr _hist_njet;

    public:

      double _jetetamax;
      double _jetptmin;
      double _lepetamax;
      double _lepptmin;


      /// Constructor
      CMS_2016_I1434354()
        : Analysis("CMS_2016_I1434354"), _jetetamax(2.5), _jetptmin(25.), _lepetamax(2.5), _lepptmin(30.)
      {    }

      /// Book histograms and initialise projections before the run
      void init()
      {
        const FinalState fs(Cuts::pT > 0. && Cuts::abseta < 6.);
        addProjection(fs, "FS");

        //book hists
        _hist_thadpt = bookHisto1D("d01-x01-y01");
        _hist_thady = bookHisto1D("d02-x01-y01");
        _hist_tleppt = bookHisto1D("d03-x01-y01");
        _hist_tlepy = bookHisto1D("d04-x01-y01");
        _hist_ttpt = bookHisto1D("d05-x01-y01");
        _hist_tty = bookHisto1D("d06-x01-y01");
        _hist_ttm = bookHisto1D("d07-x01-y01");
        _hist_njet = bookHisto1D("d08-x01-y01");

      }


      /// Perform the per-event analysis
      void analyze(const Event& event)
      {
        const double weight = event.weight();
        _muelphs.clear();
        _neutrinos.clear();
        _dressLeps.clear();
        _bJets.clear();
        _jetparticles.clear();

        const ParticleVector& FSpars = applyProjection<FinalState>(event, "FS").particles();
        for(ParticleVector::const_iterator itpar = FSpars.begin() ; itpar != FSpars.end() ; ++itpar)
        {
          const Particle& par = *itpar;
          if(isFromHadron(par.genParticle())) continue;
          int pdgid = abs(par.pdgId());
          if(pdgid == 12 || pdgid == 14 || pdgid == 16)
          {
            _neutrinos.push_back(par);
          }
          else if(pdgid == 11 || pdgid == 13 || pdgid == 22)
          {
            _muelphs.push_back(par);
          }
        }

        FastJets _fjDressLep(FinalState(), FastJets::ANTIKT, 0.1);
        _fjDressLep.calc(_muelphs);
        const Jets dressLeps = _fjDressLep.jets(Cuts::abseta < _lepetamax && Cuts::pT > _lepptmin*GeV);
        for(Jets::const_iterator itdl = dressLeps.begin() ; itdl != dressLeps.end() ; ++itdl)
        {
          Particles jetconst = itdl->particles();
          for(Particles::const_iterator itpar = jetconst.begin() ; itpar != jetconst.end(); ++itpar)
          {
            if((abs(itpar->pdgId()) == 11 || abs(itpar->pdgId()) == 13) && itpar->pt()/itdl->pt() > 0.5)
            {
              _dressLeps.push_back(*itdl);
            }
          }
        }

        if(_dressLeps.size() != 1) return;

        const vector<const HepMC::GenParticle*>& allgenpars = Rivet::particles(event.genEvent());
        for(vector<const HepMC::GenParticle*>::const_iterator itgenpar = allgenpars.begin() ; itgenpar != allgenpars.end() ; ++itgenpar)
        {
          const HepMC::GenParticle* genpar = *itgenpar;
          const int pdgid = genpar->pdg_id();
          if(isFinalBHadron(genpar))
          {
            _jetparticles.push_back(Particle(555555555, 1E-20*FourMomentum(genpar->momentum())));
            continue;
          }
          if(genpar->status() != 1) continue;

          //remove selected neutrinos
          bool isneutrino = false;
          for(Particles::const_iterator itnu = _neutrinos.begin() ; itnu != _neutrinos.end(); ++itnu)
          {
            if(itnu->genParticle() == genpar)
            {
              isneutrino = true;
              break;
            }
          }
          if(isneutrino) continue;

          //remove dressed lepton constituents
          bool islepton = false;
          for(Jets::const_iterator itdl = _dressLeps.begin() ; itdl != _dressLeps.end() ; ++itdl)
          {
            Particles jetconst = itdl->particles();
            for(Particles::const_iterator itpar = jetconst.begin() ; itpar != jetconst.end(); ++itpar)
            {
              if(itpar->genParticle() == genpar)
              {
                islepton = true;
                break;
              }
            }
          }
          if(islepton) continue;

          _jetparticles.push_back(Particle(pdgid, FourMomentum(genpar->momentum())));
        }

        FastJets _fjJets(FinalState(), FastJets::ANTIKT, 0.4);
        _fjJets.calc(_jetparticles);
        const Jets allJets = _fjJets.jets(Cuts::abseta < _jetetamax && Cuts::pT > _jetptmin*GeV);
        for(Jets::const_iterator itjet = allJets.begin() ; itjet != allJets.end() ; ++itjet)
        {
          const Particles jetconst = itjet->particles();
          for(Particles::const_iterator itpar = jetconst.begin() ; itpar != jetconst.end(); ++itpar)
          {
            if(itpar->pdgId() == 555555555)
            {
              _bJets.push_back(&(*itjet));
              break;
            }
          }
        }

        if(_bJets.size() < 2 || allJets.size() < 4) return;


        _nusum = FourMomentum(0., 0., 0., 0.);
        for(Particles::const_iterator itnu = _neutrinos.begin() ; itnu != _neutrinos.end(); ++itnu)
        {
          _nusum += itnu->momentum();
        }
        _wl = _nusum + _dressLeps[0].momentum();

        //construct top quark proxies
        double Kmin = numeric_limits<double>::max();
        for(Jets::const_iterator itaj = allJets.begin() ; itaj != allJets.end() ; ++itaj)
        {
          for(Jets::const_iterator itbj = allJets.begin() ; itbj != itaj ; ++itbj)
          {
            FourMomentum wh(itaj->momentum() + itbj->momentum());
            for(vector<const Jet*>::const_iterator ithbj = _bJets.begin() ; ithbj != _bJets.end() ; ++ithbj)
            {
              FourMomentum th(wh + (*ithbj)->momentum());
              if(&(*itaj) == *ithbj || &(*itbj) == *ithbj) continue;
              for(vector<const Jet*>::const_iterator itlbj = _bJets.begin() ; itlbj != _bJets.end() ; ++itlbj)
              {
                if(&(*itaj) == *itlbj || &(*itbj) == *itlbj || ithbj == itlbj) continue;
                FourMomentum tl(_wl + (*itlbj)->momentum());

                double K = pow(wh.mass() - 80.4, 2) + pow(th.mass() - 172.5, 2) + pow(tl.mass() - 172.5, 2);
                if(K < Kmin)
                {
                  Kmin = K;
                  _bl = *itlbj;
                  _bh = *ithbj;
                  _wja = &(*itaj); 
                  _wjb = &(*itbj);
                  _tl = tl;
                  _th = th;
                  _wh = wh;
                }
              }
            }
          }
        }

        _hist_thadpt->fill(_th.pt(), weight);
        _hist_thady->fill(abs(_th.rapidity()) , weight);
        _hist_tleppt->fill(_tl.pt() , weight);
        _hist_tlepy->fill(abs(_tl.rapidity()) , weight);
        FourMomentum tt(_tl+_th);
        _hist_ttpt->fill(tt.pt() , weight);
        _hist_tty->fill(abs(tt.rapidity()) , weight);
        _hist_ttm->fill(tt.mass() , weight);
        _hist_njet->fill(min(allJets.size()-4., 4.), weight);

      }


      /// Normalise histograms etc., after the run
      void finalize()
      {
        scale(_hist_thadpt, crossSection()/sumOfWeights());
        scale(_hist_thady, crossSection()/sumOfWeights());
        scale(_hist_tleppt, crossSection()/sumOfWeights());
        scale(_hist_tlepy, crossSection()/sumOfWeights());
        scale(_hist_ttpt, crossSection()/sumOfWeights());
        scale(_hist_tty, crossSection()/sumOfWeights());
        scale(_hist_ttm, crossSection()/sumOfWeights());
        scale(_hist_njet, crossSection()/sumOfWeights());

      }

      bool isFinalBHadron(const HepMC::GenParticle* genp)
      {
        if(isBHadron(genp->pdg_id()) == false) return false;

        const GenVertex* prodVtx = genp->production_vertex();
        if (prodVtx == nullptr) return false;
        const vector<const GenParticle*> ancestors(particles(prodVtx, HepMC::ancestors));
        for(vector<const GenParticle*>::const_iterator itan = ancestors.begin() ; itan != ancestors.end() ; ++itan)
        {
          if(isBHadron((*itan)->pdg_id()) == true) return false;
        }
        return true;
      }

      bool isBHadron(int pdgid)
      {
        int abspdgid = abs(pdgid);
        if(abspdgid <= 100) return false;

        int nq3 = (abspdgid / 10) % 10;
        int nq2 = (abspdgid / 100) % 10;
        int nq1 = (abspdgid / 1000) % 10;

        if ( nq3 == 0 ) return false;
        if ( nq1 == 5 or nq2 == 5 ) return true;

        return false;
      }

      bool isFromHadron(const HepMC::GenParticle* genp)
      {
        const GenVertex* prodVtx = genp->production_vertex();
        const vector<const GenParticle*> ancestors(particles(prodVtx, HepMC::ancestors));
        for(vector<const GenParticle*>::const_iterator itan = ancestors.begin() ; itan != ancestors.end() ; ++itan)
        {
          if((*itan)->status() == 2 && abs((*itan)->pdg_id()) > 100) return true;
        }
        return false;
      }


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1434354);


}

