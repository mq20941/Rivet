#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {

  namespace { //< only visible in this compilation unit

    class CMSGenParticle : public FinalState {
      public:
        /// @name Standard constructors and destructors.
        //@{

        /// The default constructor.
        CMSGenParticle()
          : FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV),
          _vetoIds({1000022,1000012, 1000014, 1000016,
              2000012, 2000014, 2000016, 1000039, 5100039,
              4000012, 4000014, 4000016, 9900012, 9900014, 9900016, 39}),
          _vetoIdsFromResonances({12, 13, 14, 16})
        {
          setName("CMSGenParticle");

          std::sort(_vetoIds.begin(), _vetoIds.end());
          std::sort(_vetoIdsFromResonances.begin(), _vetoIdsFromResonances.end());
        }

        //@}

        // partonicFinalState = false
        // exclude Resonances  = true
        // tauAsjets = false
      protected:
        // Apply the projection to the event
        void project(const Event& e) override;

        bool isParton(int pdgId) const {
          const int iid = abs(pdgId) % 10000;
          return (iid > 0 && iid < 6) || iid == 7 || iid == 9 || iid == 21;
        }
        bool isResonance(const int pdgId) const {
          const int iid = abs(pdgId) % 10000;
          return (iid > 21 && iid <= 42) || iid == 6 || iid == 8;
        }
        bool isIgnored(const unsigned int absId) const {
          auto pos = std::lower_bound(_vetoIds.begin(), _vetoIds.end(), absId);
          return pos != _vetoIds.end() && *pos == absId;
        }
        bool isExcludedFromResonance(const unsigned int absId) const {
          auto pos = std::lower_bound(_vetoIdsFromResonances.begin(), _vetoIdsFromResonances.end(), absId);
          return pos != _vetoIdsFromResonances.end() && *pos == absId;
        }

        int fromResonance(std::set<const GenParticle*>& invalid, std::vector<const GenParticle*>& pv, const GenParticle* p) const;

      protected:
        std::vector<unsigned int> _vetoIds;
        std::vector<unsigned int> _vetoIdsFromResonances;

    };
    
    // @brief Parton level top quark finder
    // 
    // Find top quark in the parton level directly tracking particle history.
    // This does not fit with the Rivet philosophy and can be generator dependent,
    // so please use this with your own risks.
    class PartonTop : public FinalState {
      public:
        enum TTbarMode { CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON };
        enum DecayMode { CH_HADRON = 0, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

        /// @name Standard constructors and destructors.
        //@{

        /// The default constructor.
        PartonTop() : FinalState(-MAXDOUBLE, MAXDOUBLE, 0.0*GeV)
        {
          setName("PartonTop");
        }
        
        /// Clone on the heap.
        virtual unique_ptr<Projection> clone() const {
          return unique_ptr<Projection>(new PartonTop(*this));
        }

        //@}

      public:
        TTbarMode mode() const {
          const bool isLepton1 = _mode1%3 != 0;
          const bool isLepton2 = _mode2%3 != 0;
          if      (  isLepton1 &&  isLepton2 ) return CH_FULLLEPTON;
          else if ( !isLepton1 && !isLepton2 ) return CH_FULLHADRON;
          return CH_SEMILEPTON;
        }
        DecayMode mode1() const { return _mode1; }
        DecayMode mode2() const { return _mode2; }

        Particle t1() const { return _t1; }
        Particle t2() const { return _t2; }
        Particle b1() const { return _b1; }
        Particle b2() const { return _b2; }
        ParticleVector wDecays1() const { return _wDecays1; }
        ParticleVector wDecays2() const { return _wDecays2; }
        Particle lepton1() const { return findLepton(_wDecays1); };
        Particle lepton2() const { return findLepton(_wDecays2); };

      protected:
        // Apply the projection to the event
        void project(const Event& e);

        Particle findLepton(const ParticleVector& v) const;

      private:
        DecayMode _mode1, _mode2;
        Particle _t1, _t2;
        Particle _b1, _b2;
        ParticleVector _wDecays1, _wDecays2;
    };

    // More implementation below the analysis code

  }

  class CMS_2015_I1370682_internal : public Analysis {
    public:
      CMS_2015_I1370682_internal() : Analysis("CMS_2015_I1370682_internal"),
      _doShapeOnly(false) {
      }

      void init() {
        // Parton level top quarks
        declare(PartonTop(), "ttbar");

        CMSGenParticle fsForJets;
        FastJets fj(fsForJets, FastJets::ANTIKT, 0.5);
        fj.useInvisibles();
        addProjection(fj, "Jets");

        _hSL_lepPt    = bookHisto1D("d01-x01-y01");
        _hSL_lepEta   = bookHisto1D("d02-x01-y01");
        _hSL_bqPt     = bookHisto1D("d03-x01-y01");
        _hSL_bqEta    = bookHisto1D("d04-x01-y01");
        _hSL_bbbarPt  = bookHisto1D("d05-x01-y01");
        _hSL_bbbarMass= bookHisto1D("d06-x01-y01");

        _hDL_lepPt      = bookHisto1D("d07-x01-y01");
        _hDL_lepEta     = bookHisto1D("d08-x01-y01");
        _hDL_dilepPt    = bookHisto1D("d09-x01-y01");
        _hDL_dilepMass  = bookHisto1D("d10-x01-y01");
        _hDL_bqPt       = bookHisto1D("d11-x01-y01");
        _hDL_bqEta      = bookHisto1D("d12-x01-y01");
        _hDL_bbbarPt    = bookHisto1D("d13-x01-y01");
        _hDL_bbbarMass  = bookHisto1D("d14-x01-y01");

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

      void findBAncestors(const ParticleVector& ps, std::set<int>& bAncestors) {
        bAncestors.clear();

        foreach (const Particle& p, ps) {
          GenVertex* v = p.genParticle()->production_vertex();
          if ( !v ) continue;

          foreach (const GenParticle* ancestor, Rivet::particles(v, HepMC::ancestors)) {
            if ( ancestor->status() != 2 ) continue;
            const PdgId pid = ancestor->pdg_id();
            //if ( !PID::isHadron(pid) or !PID::hasBottom(pid) ) continue;
            if ( !PID::hasBottom(pid) ) continue;

            GenVertex* av = ancestor->production_vertex();
            if ( !av ) continue;

            bool isDuplicated = false;
            foreach (const GenParticle* ap, Rivet::particles(av, HepMC::parents)) {
              if ( p.genParticle() != ap && pid == ap->pdg_id() ) {
                isDuplicated = true;
                break;
              }
            }

            if ( !isDuplicated ) bAncestors.insert(ancestor->barcode());
          }
        }
      }

      void analyze(const Event& event) {
        const double weight = event.weight();

        // Get the parton level ttbar candidate
        const PartonTop& ttbarState = applyProjection<PartonTop>(event, "ttbar");
        // Do the analysis only for the ttbar full leptonic or semileptonic channel, without tau decay
        if ( ttbarState.mode() != PartonTop::CH_SEMILEPTON and
            ttbarState.mode() != PartonTop::CH_FULLLEPTON ) vetoEvent;
        if ( ttbarState.mode1() >= PartonTop::CH_TAU_HADRON ||
            ttbarState.mode2() >= PartonTop::CH_TAU_HADRON ) vetoEvent;

        // Parton level at full phase space
        // Fill top quarks they are defined in the parton level, full phase space
        const FourMomentum& t1P4 = ttbarState.t1().momentum();
        const FourMomentum& t2P4 = ttbarState.t2().momentum();
        const double t1Pt = t1P4.pT(), t2Pt = t2P4.pT();
        const FourMomentum ttbarP4 = t1P4+t2P4;
        const FourMomentum t1P4AtCM = LorentzTransform::mkFrameTransformFromBeta(ttbarP4.betaVec()).transform(t1P4);
        const double dPhi = deltaPhi(t1P4.phi(), t2P4.phi());

        if ( ttbarState.mode() == PartonTop::CH_SEMILEPTON ) {
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
        else if ( ttbarState.mode() == PartonTop::CH_FULLLEPTON ) {
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

        // Find leptons
        Particles lCands;
        if ( ttbarState.mode() == PartonTop::CH_SEMILEPTON ) {
          const Particle l1 = ttbarState.lepton1();
          const Particle l2 = ttbarState.lepton2();
          if      ( l1.pdgId() ) lCands.push_back(l1);
          else if ( l2.pdgId() ) lCands.push_back(l2);
          else vetoEvent; // This never happens if PartonTop runs properly

          // Apply the particle level phase space cut
          if ( lCands[0].pT() <= 33 or std::abs(lCands[0].eta()) >= 2.1 ) vetoEvent;
        }
        else if ( ttbarState.mode() == PartonTop::CH_FULLLEPTON ) {
          lCands.push_back(ttbarState.lepton1());
          lCands.push_back(ttbarState.lepton2());
          if ( !lCands[0].pdgId() or !lCands[1].pdgId() ) vetoEvent; // This never happens if PartonTop runs properly

          if ( lCands[0].pT() < lCands[1].pT() ) std::swap(lCands[0], lCands[1]);
          const double l1Pt = lCands[0].pT(), l1Abseta = std::abs(lCands[0].eta());
          const double l2Pt = lCands[1].pT(), l2Abseta = std::abs(lCands[1].eta());

          // Apply the particle level phase space cut
          if ( l1Pt <= 20 or l1Abseta >= 2.4 or l2Pt <= 20 or l2Abseta >= 2.4 ) vetoEvent;
          if ( (lCands[0].momentum()+lCands[1].momentum()).mass() < 20 ) vetoEvent;
        }

        // Build genJets
        const Jets& jetsIn = applyProjection<JetAlg>(event, "Jets").jetsByPt(30*GeV);
        Jets jets;
        foreach ( const Jet& jet, jetsIn ) {
          if ( std::abs(jet.eta()) > 2.4 ) continue;
          bool isOverlapped = false;
          foreach ( const Particle& lCand, lCands ) {
            if ( deltaR(lCand.momentum(), jet.momentum()) < 0.3 ) {
              isOverlapped = true;
              break;
            }
          }
          if ( isOverlapped ) continue;
          jets.push_back(jet);
        }
        if ( ttbarState.mode() == PartonTop::CH_SEMILEPTON and jets.size() < 4 ) vetoEvent;
        else if ( ttbarState.mode() == PartonTop::CH_FULLLEPTON and jets.size() < 2 ) vetoEvent;

        // Get Leading two jets
        std::set<int> bAncestors;
        const Jet* bjet1 = 0, * bjet2 = 0;
        foreach (const Jet& jet, jets) {
          if ( !bjet1 ) {
            findBAncestors(jet.particles(), bAncestors);
            if ( !bAncestors.empty() ) bjet1 = &jet;
          } else if ( !bjet2 ) {
            std::set<int> bAncestors2;
            findBAncestors(jet.particles(), bAncestors2);
            bool hasUniqueB = false;
            foreach (int barCode, bAncestors2) {
              if ( bAncestors.find(barCode) == bAncestors.end() ) {
                hasUniqueB = true;
                break;
              }
            }
            if ( hasUniqueB ) {
              bjet2 = &jet;
              break;
            }
          } else break;
        }
        // Require at least 2 b jets.
        // It is safe with checking bjet2 only. But let's check bjet1 for clearity.
        if ( !bjet2 or !bjet1 ) vetoEvent;

        const FourMomentum& b1P4 = bjet1->momentum();
        const FourMomentum& b2P4 = bjet2->momentum();
        const FourMomentum bbP4 = b1P4+b2P4;
        const double b1Pt = b1P4.pT(), b1Eta = b1P4.eta();
        const double b2Pt = b2P4.pT(), b2Eta = b2P4.eta();
        const double bbPt = bbP4.pT(), bbMass = bbP4.mass();

        // Find leptons, apply channel dependent phase space cuts, fill histograms
        if ( ttbarState.mode() == PartonTop::CH_SEMILEPTON ) {
          // Do the semileptonic channel
          const FourMomentum& lP4 = lCands[0].momentum();
          const double lPt = lP4.pT(), lEta = lP4.eta();

          _hSL_lepPt->fill(lPt, weight);
          _hSL_lepEta->fill(lEta, weight);
          _hSL_bqPt->fill(b1Pt, weight);
          _hSL_bqPt->fill(b2Pt, weight);
          _hSL_bqEta->fill(b1Eta, weight);
          _hSL_bqEta->fill(b2Eta, weight);
          _hSL_bbbarPt->fill(bbPt, weight);
          _hSL_bbbarMass->fill(bbMass, weight);
        }
        else if ( ttbarState.mode() == PartonTop::CH_FULLLEPTON ) {
          // Do the dileptonic channel
          const FourMomentum& l1P4 = lCands[0].momentum();
          const FourMomentum& l2P4 = lCands[1].momentum();
          const FourMomentum dilP4 = l1P4+l2P4;
          const double l1Pt = l1P4.pT(), l1Eta = l1P4.eta();
          const double l2Pt = l2P4.pT(), l2Eta = l2P4.eta();

          _hDL_lepPt->fill(l1Pt, weight);
          _hDL_lepPt->fill(l2Pt, weight);
          _hDL_lepEta->fill(l1Eta, weight);
          _hDL_lepEta->fill(l2Eta, weight);
          _hDL_dilepPt->fill(dilP4.pT(), weight);
          _hDL_dilepMass->fill(dilP4.mass(), weight);

          _hDL_bqPt->fill(b1Pt, weight);
          _hDL_bqPt->fill(b2Pt, weight);
          _hDL_bqEta->fill(b1Eta, weight);
          _hDL_bqEta->fill(b2Eta, weight);
          _hDL_bbbarPt->fill(bbPt, weight);
          _hDL_bbbarMass->fill(bbMass, weight);
        }
      };

      void finalize() {
        if ( _doShapeOnly ) {
          normalize(_hSL_lepPt    );
          normalize(_hSL_lepEta   );
          normalize(_hSL_bqPt     );
          normalize(_hSL_bqEta    );
          normalize(_hSL_bbbarPt  );
          normalize(_hSL_bbbarMass);

          normalize(_hDL_lepPt      );
          normalize(_hDL_lepEta     );
          normalize(_hDL_dilepPt    );
          normalize(_hDL_dilepMass  );
          normalize(_hDL_bqPt       );
          normalize(_hDL_bqEta      );
          normalize(_hDL_bbbarPt    );
          normalize(_hDL_bbbarMass  );

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
        }
        else {
          const double s = 1./sumOfWeights();
          scale(_hSL_lepPt    ,s);
          scale(_hSL_lepEta   ,s);
          scale(_hSL_bqPt     ,s);
          scale(_hSL_bqEta    ,s);
          scale(_hSL_bbbarPt  ,s);
          scale(_hSL_bbbarMass,s);

          scale(_hDL_lepPt      ,s);
          scale(_hDL_lepEta     ,s);
          scale(_hDL_dilepPt    ,s);
          scale(_hDL_dilepMass  ,s);
          scale(_hDL_bqPt       ,s);
          scale(_hDL_bqEta      ,s);
          scale(_hDL_bbbarPt    ,s);
          scale(_hDL_bbbarMass  ,s);

          scale(_hSL_topPt        ,s);
          scale(_hSL_topPtTtbarSys,s);
          scale(_hSL_topY         ,s);
          scale(_hSL_ttbarDelPhi  ,s);
          scale(_hSL_topPtLead    ,s);
          scale(_hSL_topPtSubLead ,s);
          scale(_hSL_ttbarPt      ,s);
          scale(_hSL_ttbarY       ,s);
          scale(_hSL_ttbarMass    ,s);

          scale(_hDL_topPt        ,s);
          scale(_hDL_topPtTtbarSys,s);
          scale(_hDL_topY         ,s);
          scale(_hDL_ttbarDelPhi  ,s);
          scale(_hDL_topPtLead    ,s);
          scale(_hDL_topPtSubLead ,s);
          scale(_hDL_ttbarPt      ,s);
          scale(_hDL_ttbarY       ,s);
          scale(_hDL_ttbarMass    ,s);
        }

      };

    private:
      const bool _doShapeOnly;

      Histo1DPtr _hSL_lepPt;
      Histo1DPtr _hSL_lepEta;
      Histo1DPtr _hSL_bqPt;
      Histo1DPtr _hSL_bqEta;
      Histo1DPtr _hSL_bbbarPt;
      Histo1DPtr _hSL_bbbarMass;

      Histo1DPtr _hDL_lepPt;
      Histo1DPtr _hDL_lepEta;
      Histo1DPtr _hDL_dilepPt;
      Histo1DPtr _hDL_dilepMass;
      Histo1DPtr _hDL_bqPt;
      Histo1DPtr _hDL_bqEta;
      Histo1DPtr _hDL_bbbarPt;
      Histo1DPtr _hDL_bbbarMass;

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
  AnalysisBuilder<CMS_2015_I1370682_internal> plugin_CMS_2015_I1370682_internal;
  
  ///////////////

  // Additional classes implementation
  namespace {
  
    void CMSGenParticle::project(const Event& e) {
      _theParticles.clear();

      std::vector<const GenParticle*> particles = Rivet::particles(e.genEvent());
      std::set<const GenParticle*> selected, invalid;

      foreach (const GenParticle* p, particles ) {
        if ( invalid.find(p) != invalid.end() ) continue;
        if ( p->status() == 1 ) selected.insert(p);
      }

      foreach (const GenParticle* p, selected ) {
        if ( invalid.find(p) != invalid.end() ) continue;

        if ( fromResonance(invalid, particles, p) ) {
          invalid.insert(p);
          continue;
        }

        if ( isIgnored(p->pdg_id()) ) continue;

        _theParticles.push_back(p);
      }
    }

    int CMSGenParticle::fromResonance(std::set<const GenParticle*>& invalid, std::vector<const GenParticle*>& pv, const GenParticle* p) const {
      const int id = p->pdg_id();
      const unsigned int aid = std::abs(id);

      if ( invalid.find(p) != invalid.end() ) return 2;
      if ( isResonance(aid) && p->status() == 3 ) return 1;
      if ( !isIgnored(aid) && isParton(aid) ) return 0;

      GenVertex* vtx = p->production_vertex();
      if ( !vtx ) return 0;

      std::vector<GenParticle*> mothers = Rivet::particles(vtx, HepMC::parents);
      if ( mothers.empty() ) return 0;

      foreach (const GenParticle* mother, mothers) {
        const int result = fromResonance(invalid, pv, mother);
        switch ( result ) {
          case 0: break;
          case 1:
                  if ( mother->pdg_id() == id or isResonance(id) ) return 1;
                  if ( !isExcludedFromResonance(aid) ) break;
          case 2: return 2;
        }
      }
      return 0;
    }
    
    void PartonTop::project(const Event& e) {
      _theParticles.clear();
      _wDecays1.clear();
      _wDecays2.clear();
      _mode1 = _mode2 = CH_HADRON; // Set default decay mode to full-hadronic
      _t1 = _t2 = _b1 = _b2 = Particle();

      const double ptmin = 0;
      const double etamin = -MAXDOUBLE, etamax = MAXDOUBLE;

      int nTop = 0;
      bool isTau1 = false, isTau2 = false;
      foreach (const GenParticle* p, Rivet::particles(e.genEvent())) {
        const int pdgId = p->pdg_id();
        const int absId = abs(pdgId);
        if ( absId > 20 ) continue; // We are only interested in quarks and leptons
        //if ( PID::isHadron(pdgId) ) continue; // skip hadrons
        //if ( pdgId == 22 ) continue; // skip photons
        //if ( pdgId == 91 or pdgId == 92 ) continue; // Skip cluster, strings

        if ( isZero(p->momentum().perp()) || p->momentum().perp() < ptmin ) continue;
        if ( !inRange(p->momentum().eta(), etamin, etamax) ) continue;

        // Avoid double counting by skipping if particle ID == parent ID
        std::vector<GenParticle*> pps;
        if ( absId == 6 and p->end_vertex() != 0 ) {
          pps = Rivet::particles(p->end_vertex(), HepMC::children);
        }
        else if ( absId != 6 and p->production_vertex() != 0 )
        {
          pps = Rivet::particles(p->production_vertex(), HepMC::parents);
        }
        else continue;

        bool isDuplicated = false;
        foreach (GenParticle* pp, pps) {
          if ( p != pp && p->pdg_id() == pp->pdg_id() ) {
            isDuplicated = true;
            break;
          }
        }
        if ( isDuplicated ) continue;

        // Build Rivet::Particle
        Particle rp(*p);
        // Skip particles from hadronization (and keep tau decay)
        if ( rp.fromDecay() and !rp.hasAncestor(15) and !rp.hasAncestor(-15) ) continue;

        if      ( pdgId ==  6 ) { nTop++; _t1 = rp; }
        else if ( pdgId == -6 ) { nTop++; _t2 = rp; }
        else if ( pdgId ==  5 and rp.pT() > _b1.pT() ) _b1 = rp;
        else if ( pdgId == -5 and rp.pT() > _b2.pT() ) _b2 = rp;
        else if ( absId <= 16 && rp.hasAncestor( 24) ) {
          if ( pdgId == -15 ) isTau1 = true;
          else if ( pdgId == -11 ) _mode1 = CH_ELECTRON;
          else if ( pdgId == -13 ) _mode1 = CH_MUON;
          _wDecays1.push_back(rp);
        }
        else if ( absId <= 16 && rp.hasAncestor(-24) ) {
          if ( pdgId == 15 ) isTau2 = true;
          else if ( pdgId == 11 ) _mode2 = CH_ELECTRON;
          else if ( pdgId == 13 ) _mode2 = CH_MUON;
          _wDecays2.push_back(rp);
        }
      }
      if ( isTau1 ) _mode1 = static_cast<DecayMode>(_mode1+3);
      if ( isTau2 ) _mode2 = static_cast<DecayMode>(_mode2+3);

    }

    Particle PartonTop::findLepton(const ParticleVector& v) const {
      Particle pp = Particle();
      foreach (const Particle& p, v) {
        const int aid = std::abs(p.pdgId());
        if ( (aid == 11 or aid == 13) and  p.pT() > pp.pT() ) pp = p;
      }
      return pp;
    }
  
  }

}
