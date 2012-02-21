// -*- C++ -*-
//
// Package:    MyAnalyzer
// Class:      MyAnalyzer
//
/**\class MyAnalyzer MyAnalyzer.cc MyAnalysis/MyAnalyzer/src/MyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dinko Ferencek
//         Created:  Mon Sep 12 15:06:41 CDT 2011
// $Id: MyAnalyzer_MainAnalysis_DijetBBTag_2011.cc,v 1.12 2012/02/09 18:50:30 ferencek Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/MergeableCounter.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

// BaseClass
#include "MyAnalysis/MyAnalyzer/interface/BaseClass.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT
#include <TF1.h>
#include <TLorentzVector.h>


using namespace std;

//
// class declaration
//

class BTagScaleFactorCalculator
{
 public:
   BTagScaleFactorCalculator();
   void init(const double TCHEL_SFb, const double TCHEL_SFl, const double TCHPT_SFb, const double TCHPT_SFl, const double SSVHPT_SFb, const double SSVHPT_SFl);
   double scaleFactor(const int partonFlavor, const int btagger);
   double scaleFactor(const int partonFlavor, const double jetPt, const double jetEta, const int btagger);
   double scaleFactorBC_TCHEL(const double jetPt, const double jetEta);
   double scaleFactorUDSG_TCHEL(const double jetPt, const double jetEta);

 private:
   double TCHEL_SFb_;
   double TCHEL_SFl_;
   double TCHPT_SFb_;
   double TCHPT_SFl_;
   double SSVHPT_SFb_;
   double SSVHPT_SFl_;
   TF1 *TCHEL_SFb_0to2p4;
   TF1 *TCHEL_SFl_0to2p4;
   TF1 *TCHEL_SFl_0to0p5;
   TF1 *TCHEL_SFl_0p5to1p0;
   TF1 *TCHEL_SFl_1p0to1p5;
   TF1 *TCHEL_SFl_1p5to2p4;
};

class MyAnalyzer : public BaseClass, public edm::EDFilter {
   public:
      explicit MyAnalyzer(const edm::ParameterSet&);
      ~MyAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      double bTagEventWeight(const vector<double>& SFsForBTaggedJets, const unsigned int nBTags);

      // ----------member data ---------------------------
      HLTConfigProvider hltConfig;
      edm::InputTag     hltInputTag;
      edm::LumiReWeighting LumiWeights;
      BTagScaleFactorCalculator sfCalculator;
};

//
// constructors and destructor
//

MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig) :
  BaseClass(iConfig)
{
   //now do whatever initialization is needed
   hltInputTag = iConfig.getParameter<edm::InputTag>("HLTInputTag");

   sfCalculator.init(getPreCutValue1("TCHEL_SFb"),getPreCutValue1("TCHEL_SFl"),getPreCutValue1("TCHPT_SFb"),getPreCutValue1("TCHPT_SFl"),getPreCutValue1("SSVHPT_SFb"),getPreCutValue1("SSVHPT_SFl"));
}

MyAnalyzer::~MyAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void
MyAnalyzer::beginJob()
{
   //#####################################################################################
   //########################### User's code starts here #################################/
   
   // book your histograms here
   CreateUserTH1D("h1_J1J2PartonFlavor;Parton Flavor (PDG ID);Entries", 51, -0.5, 50.5);
   CreateUserTH1D("h1_J1J2HeavyFlavor;Heavy Flavor;Entries", 2, -0.5, 1.5);
   CreateUserTH1D("h1_nMuons_vs_DijetMass_pretag;Dijet Mass [GeV];nMuons", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_nMuons_vs_DijetMass;Dijet Mass [GeV];nMuons", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   
   CreateUserTH2D("h2_EtaJ2_vs_EtaJ1;#eta_{1};#eta_{2}", getHistoNBins("EtaJ1"), getHistoMin("EtaJ1"), getHistoMax("EtaJ1"), getHistoNBins("EtaJ1"), getHistoMin("EtaJ1"), getHistoMax("EtaJ1"));

   int doEventBins = int(getPreCutValue1("doEventBins"));
   int doSFReweighting = int(getPreCutValue1("doSFReweighting"));
   if( doEventBins )
   {
     // histograms for in-situ b-tag SF measurement
     // 0 muon, maxEta<1.2 case:
     CreateUserTH2D("h2_n0_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     CreateUserTH2D("h2_n1_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     CreateUserTH2D("h2_n2_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     if( doSFReweighting )
     {
       CreateUserTH2D("h2_N22_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N21_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N20_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N12_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N11_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N10_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N02_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N01_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N00_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     }
     else
     {
       CreateUserTH2D("h2_N220_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N210_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N200_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N111_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N110_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N101_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N100_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N002_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N001_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N000_0mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     }
     // 0 muon, maxEta>=1.2 case:
     CreateUserTH2D("h2_n0_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     CreateUserTH2D("h2_n1_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     CreateUserTH2D("h2_n2_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     if( doSFReweighting )
     {
       CreateUserTH2D("h2_N22_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N21_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N20_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N12_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N11_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N10_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N02_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N01_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N00_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     }
     else
     {
       CreateUserTH2D("h2_N220_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N210_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N200_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N111_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N110_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N101_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N100_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N002_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N001_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N000_0mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     }
     // >=1 muon, maxEta<1.2 case:
     CreateUserTH2D("h2_n0_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     CreateUserTH2D("h2_n1_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     CreateUserTH2D("h2_n2_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     if( doSFReweighting )
     {
       CreateUserTH2D("h2_N22_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N21_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N20_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N12_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N11_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N10_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N02_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N01_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N00_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     }
     else
     {
       CreateUserTH2D("h2_N220_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N210_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N200_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N111_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N110_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N101_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N100_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N002_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N001_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N000_ge1mu_maxEta_lt_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     }
     // >=1 muon, maxEta<1.2 case:
     CreateUserTH2D("h2_n0_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     CreateUserTH2D("h2_n1_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     CreateUserTH2D("h2_n2_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     if( doSFReweighting )
     {
       CreateUserTH2D("h2_N22_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N21_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N20_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N12_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N11_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N10_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N02_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N01_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N00_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     }
     else
     {
       CreateUserTH2D("h2_N220_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N210_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N200_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N111_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N110_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N101_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N100_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N002_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N001_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
       CreateUserTH2D("h2_N000_ge1mu_maxEta_ge_1p2;Dijet Mass [GeV];nPV",getHistoNBins("DijetMass"),getHistoMin("DijetMass"),getHistoMax("DijetMass"),getHistoNBins("nGoodVertices"),getHistoMin("nGoodVertices"),getHistoMax("nGoodVertices"));
     }
   }
   
   // initialize your variables here

   // Summer11 PU_S4 distribution obtained by only looking at the in-time crossing.  This is the "spike+smear" distribution, RECOMMENDED FOR REWEIGHTING.
   double PileUpDistMC_ObservedBX0_d[35] = {1.45346E-01,6.42802E-02,6.95255E-02,6.96747E-02,6.92955E-02,6.84997E-02,6.69528E-02,6.45515E-02,6.09865E-02,5.63323E-02,5.07322E-02,4.44681E-02,3.79205E-02,3.15131E-02,2.54220E-02,2.00184E-02,1.53776E-02,1.15387E-02,8.47608E-03,6.08715E-03,4.28255E-03,2.97185E-03,2.01918E-03,1.34490E-03,8.81587E-04,5.69954E-04,3.61493E-04,2.28692E-04,1.40791E-04,8.44606E-05,5.10204E-05,3.07802E-05,1.81401E-05,1.00201E-05,5.80004E-06};
   vector<float> PileUpDistMC_ObservedBX0(PileUpDistMC_ObservedBX0_d, PileUpDistMC_ObservedBX0_d + sizeof(PileUpDistMC_ObservedBX0_d) / sizeof(double) );
   // Run2011A+Run2011B pile-up distribution
   double PileUpDistData_Observed_d[35] = {13446512.0, 59065300.0, 140902672.0, 241301168.0, 333744896.0, 398710976.0, 430106432.0, 432283008.0, 413820192.0, 382845984.0, 345163680.0, 304343808.0, 262555024.0, 221330752.0, 181982560.0, 145689760.0, 113413440.0, 85778864.0, 63012392.0, 44959588.0, 31169040.0, 21007856.0, 13775880.0, 8796407.0, 5474417.5, 3323775.5, 1970637.75, 1142040.5, 647538.625, 359547.1875, 195673.15625, 104459.9453125, 54745.15234375, 28185.56640625, 28005.548828125};
   vector<float> PileUpDistData_Observed(PileUpDistData_Observed_d, PileUpDistData_Observed_d + sizeof(PileUpDistData_Observed_d) / sizeof(double) );
   // Run2011A pile-up distribution
   double PileUpDistDataRun2011A_Observed_d[35] = {12965370.0, 55851368.0, 129329360.0, 212133600.0, 276137728.0, 303603552.0, 293257504.0, 255632864.0, 204970016.0, 153263664.0, 107935616.0, 72100608.0, 45912988.0, 27970044.0, 16342576.0, 9175983.0, 4958610.0, 2582392.75, 1297695.75, 629975.0625, 295784.25, 134469.671875, 59260.0703125, 25343.8671875, 10530.08984375, 4255.04833984375, 1673.949462890625, 641.7764892578125, 240.02249145507812, 87.650428771972656, 31.280984878540039, 10.919528007507324, 3.7314565181732178, 1.2492282390594482, 0.60236752033233643};
   vector<float> PileUpDistDataRun2011A_Observed(PileUpDistDataRun2011A_Observed_d, PileUpDistDataRun2011A_Observed_d + sizeof(PileUpDistDataRun2011A_Observed_d) / sizeof(double) );
   // Run2011B pile-up distribution
   double PileUpDistDataRun2011B_Observed_d[35] = {481141.72472439631, 3213932.6435495433, 11573306.879849812, 29167566.222834267, 57607165.266923025, 95107416.620279759, 136848927.47142315, 176650126.03691837, 208850159.2230134, 229582323.64430469, 237228075.21969861, 232243229.29014349, 216642041.8648839, 193360704.19149143, 165639997.39323151, 136513771.36073488, 108454824.41062045, 83196475.076050699, 61714695.179401994, 44329612.988056183, 30873256.083125703, 20873387.048466656, 13716620.52313938, 8771062.8496078029, 5463887.6202946259, 3319520.4784049694, 1968963.8927261904, 1141398.7709201651, 647298.63969274936, 359459.56447747251, 195641.88301312848, 104449.02504469089, 54741.419995713957, 28184.317687585703, 28004.947183609398};
   vector<float> PileUpDistDataRun2011B_Observed(PileUpDistDataRun2011B_Observed_d, PileUpDistDataRun2011B_Observed_d + sizeof(PileUpDistDataRun2011B_Observed_d) / sizeof(double) );
   
   int puReweightingEra = int(getPreCutValue1("puReweightingEra"));
   if( puReweightingEra==1 ) LumiWeights = edm::LumiReWeighting(PileUpDistMC_ObservedBX0, PileUpDistDataRun2011A_Observed);
   else if( puReweightingEra==2 ) LumiWeights = edm::LumiReWeighting(PileUpDistMC_ObservedBX0, PileUpDistDataRun2011B_Observed);
   else LumiWeights = edm::LumiReWeighting(PileUpDistMC_ObservedBX0, PileUpDistData_Observed);
   
   //############################# User's code ends here #################################
   //#####################################################################################
}

// ------------ method called when starting to processes a run  ------------
bool
MyAnalyzer::beginRun(edm::Run& iRun, const edm::EventSetup& iSetup)
{
  bool changed = true;
  if (hltConfig.init(iRun, iSetup, hltInputTag.process(), changed)) {
    // if init returns TRUE, initialisation has succeeded!
    edm::LogInfo("MyAnalyzer::beginRun") << "HLT config with process name " << hltInputTag.process() << " successfully extracted";
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("MyAnalyzer::beginRun") << "Error! HLT config extraction with process name " << hltInputTag.process() << " failed";
    // In this case, all access methods will return empty values!
    return false;
  }
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool
MyAnalyzer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called on each new Event  ------------
bool
MyAnalyzer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   bool ret = false;
   // event weight (by default set to 1)
   double eventWeight = 1;

   //#####################################################################################
   //########################### User's code starts here #################################

   int doPUReweighting = int(getPreCutValue1("doPUReweighting"));
   int doSFReweighting = int(getPreCutValue1("doSFReweighting"));
   int useFixedSFs = int(getPreCutValue1("useFixedSFs"));
   int btagger = int(getPreCutValue1("btagger"));
   int useHFkFactor = int(getPreCutValue1("useHFkFactor"));
   int matchingType = int(getPreCutValue1("matchingType"));
   double matchingRadius = getPreCutValue1("matchingRadius");
   int doEventBins = int(getPreCutValue1("doEventBins"));
   int doEventPrintout = int(getPreCutValue1("doEventPrintout"));
   double METoSumET_cut = getPreCutValue1("METoSumET_cut");
   double DeltaPhiJ1J2_cut = getPreCutValue1("DeltaPhiJ1J2_cut");
   
   // grab necessary objects from the event
//    edm::Handle<edm::TriggerResults> triggerResults;
//    iEvent.getByLabel(hltInputTag, triggerResults);

   edm::Handle<vector<bool> > PVIsFake;
   iEvent.getByLabel(edm::InputTag("Vertices:IsFake"), PVIsFake);
   edm::Handle<vector<double> > PVX;
   iEvent.getByLabel(edm::InputTag("Vertices:X"), PVX);
   edm::Handle<vector<double> > PVY;
   iEvent.getByLabel(edm::InputTag("Vertices:Y"), PVY);
   edm::Handle<vector<double> > PVZ;
   iEvent.getByLabel(edm::InputTag("Vertices:Z"), PVZ);
   edm::Handle<vector<double> > PVNDF;
   iEvent.getByLabel(edm::InputTag("Vertices:NDF"), PVNDF);
   
   edm::Handle<vector<unsigned int> > NPU;
   iEvent.getByLabel(edm::InputTag("GenEventInfo:PileUpNumberOfInteractions"), NPU);
   edm::Handle<vector<int> > BX;
   iEvent.getByLabel(edm::InputTag("GenEventInfo:PileUpBunchCrossing"), BX);

   edm::Handle<vector<double> > GenParticlePt;
   iEvent.getByLabel(edm::InputTag("GenParticles:Pt"), GenParticlePt);
   edm::Handle<vector<double> > GenParticleEta;
   iEvent.getByLabel(edm::InputTag("GenParticles:Eta"), GenParticleEta);
   edm::Handle<vector<double> > GenParticlePhi;
   iEvent.getByLabel(edm::InputTag("GenParticles:Phi"), GenParticlePhi);
   edm::Handle<vector<double> > GenParticleE;
   iEvent.getByLabel(edm::InputTag("GenParticles:Energy"), GenParticleE);
   edm::Handle<vector<int> > GenParticlePdgId;
   iEvent.getByLabel(edm::InputTag("GenParticles:PdgId"), GenParticlePdgId);
   edm::Handle<vector<int> > GenParticleStatus;
   iEvent.getByLabel(edm::InputTag("GenParticles:Status"), GenParticleStatus);
   
   edm::Handle<bool> passHBHENoiseFilter;
   iEvent.getByLabel(edm::InputTag("EventSelection:PassHBHENoiseFilter"), passHBHENoiseFilter);
   edm::Handle<bool> passBeamHaloFilterTight;
   iEvent.getByLabel(edm::InputTag("EventSelection:PassBeamHaloFilterTight"), passBeamHaloFilterTight);
   edm::Handle<bool> passTrackingFailure;
   iEvent.getByLabel(edm::InputTag("EventSelection:PassTrackingFailure"), passTrackingFailure);
   edm::Handle<bool> passEcalMaskedCellDRFilter;
   iEvent.getByLabel(edm::InputTag("EventSelection:PassEcalMaskedCellDRFilter"), passEcalMaskedCellDRFilter);
   edm::Handle<bool> passCaloBoundaryDRFilter;
   iEvent.getByLabel(edm::InputTag("EventSelection:PassCaloBoundaryDRFilter"), passCaloBoundaryDRFilter);

   edm::Handle<vector<double> > MET;
   iEvent.getByLabel(edm::InputTag("PFMET:Mag"), MET);
   edm::Handle<vector<double> > SumET;
   iEvent.getByLabel(edm::InputTag("PFMET:SumET"), SumET);
   
   edm::Handle<vector<double> > PFJetPt_;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Pt"), PFJetPt_);
   edm::Handle<vector<double> > PFJetPtRaw;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:PtRaw"), PFJetPtRaw);
   edm::Handle<vector<double> > PFJetEta;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Eta"), PFJetEta);
   edm::Handle<vector<double> > PFJetPhi;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Phi"), PFJetPhi);
   edm::Handle<vector<double> > PFJetE_;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Energy"), PFJetE_);
   edm::Handle<vector<double> > PFJetUnc;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:JECUnc"), PFJetUnc);
   edm::Handle<vector<int> > PFJetPassJetID;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:PassTightID"), PFJetPassJetID);
   edm::Handle<vector<double> > PFJetSSVHE;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:SimpleSecondaryVertexHighEffBTag"), PFJetSSVHE);
   edm::Handle<vector<double> > PFJetSSVHP;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:SimpleSecondaryVertexHighPurBTag"), PFJetSSVHP);
   edm::Handle<vector<double> > PFJetTCHE;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:TrackCountingHighEffBTag"), PFJetTCHE);
   edm::Handle<vector<double> > PFJetTCHP;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:TrackCountingHighPurBTag"), PFJetTCHP);
   edm::Handle<vector<int> > PFJetPartonFlavor;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:PartonFlavor"), PFJetPartonFlavor);
   
   edm::Handle<vector<double> > MuonPt;
   iEvent.getByLabel(edm::InputTag("Muons:Pt"), MuonPt);
   edm::Handle<vector<double> > MuonEta;
   iEvent.getByLabel(edm::InputTag("Muons:Eta"), MuonEta);
   edm::Handle<vector<double> > MuonPhi;
   iEvent.getByLabel(edm::InputTag("Muons:Phi"), MuonPhi);
   edm::Handle<vector<bool> > MuonIsGlobal;
   iEvent.getByLabel(edm::InputTag("Muons:IsGlobal"), MuonIsGlobal);
   edm::Handle<vector<int> > MuonNValidTrackerHits;
   iEvent.getByLabel(edm::InputTag("Muons:NHitsTracker"), MuonNValidTrackerHits);
   edm::Handle<vector<int> > MuonNValidPixelHits;
   iEvent.getByLabel(edm::InputTag("Muons:NHitsPixel"), MuonNValidPixelHits);
   edm::Handle<vector<int> > MuonNLostTrackerHitsOut;
   iEvent.getByLabel(edm::InputTag("Muons:NLostHitsTrackerOut"), MuonNLostTrackerHitsOut);
   edm::Handle<vector<int> > MuonNValidMuonHits;
   iEvent.getByLabel(edm::InputTag("Muons:NHitsMuon"), MuonNValidMuonHits);
   edm::Handle<vector<int> > MuonNMatches;
   iEvent.getByLabel(edm::InputTag("Muons:NMatchedChambers"), MuonNMatches);
   edm::Handle<vector<double> > MuonChi2;
   iEvent.getByLabel(edm::InputTag("Muons:Chi2"), MuonChi2);
   edm::Handle<vector<double> > MuonNdof;
   iEvent.getByLabel(edm::InputTag("Muons:Ndof"), MuonNdof);
   edm::Handle<vector<double> > MuonTrkChi2;
   iEvent.getByLabel(edm::InputTag("Muons:TrkChi2"), MuonTrkChi2);
   edm::Handle<vector<double> > MuonTrkNdof;
   iEvent.getByLabel(edm::InputTag("Muons:TrkNdof"), MuonTrkNdof);
   edm::Handle<vector<double> > MuonRefPtZ;
   iEvent.getByLabel(edm::InputTag("Muons:RefPtZ"), MuonRefPtZ);

   // apply pile-up reweighting
   if( !iEvent.isRealData() && doPUReweighting )
   {
     int npu = -1;
     for( vector<int>::const_iterator bxIt = BX->begin(); bxIt != BX->end(); ++bxIt )
     {
        if( *bxIt == 0 ) {
          npu = NPU->at( distance(BX->begin(),bxIt) );
          break;
        }
     }
     // set the event weight
     eventWeight = LumiWeights.weight( npu );
   }

   int nStatus3_bQuarks = 0;
   int nStatus2_bQuarks = 0;

   for(size_t i=0; i<GenParticlePt->size(); i++)
   {
     if( abs(GenParticlePdgId->at(i))==5 && GenParticleStatus->at(i)==3 ) nStatus3_bQuarks++;
     if( abs(GenParticlePdgId->at(i))==5 && GenParticleStatus->at(i)==2 ) nStatus2_bQuarks++;
   }

   auto_ptr<std::vector<double> >  PFJetPt ( new std::vector<double>() );
   auto_ptr<std::vector<double> >  PFJetE  ( new std::vector<double>() );

   for(size_t i=0; i<PFJetPt_->size(); i++)
   {
     double JES_ScaleFactor = 1.;
     if( !iEvent.isRealData() ) JES_ScaleFactor = 1. + getPreCutValue1("JES_Shift")*PFJetUnc->at(i);

     PFJetPt->push_back( PFJetPt_->at(i)*JES_ScaleFactor );
     PFJetE ->push_back( PFJetE_ ->at(i)*JES_ScaleFactor );
   }
   
   // loop over PFJets and select PFJets that pass JetID requirements
   vector<int> v_idx_pfjet_JetID;
   for(size_t i=0; i<PFJetPt->size(); i++)
   {
     if( !PFJetPassJetID->at(i) ) continue;
     v_idx_pfjet_JetID.push_back(i);
   }

   // loop over primary vertices and select good ones
   vector<int> v_idx_goodPV;
   for(size_t i=0; i<PVX->size(); i++)
   {
     double rho = sqrt(PVX->at(i)*PVX->at(i) + PVY->at(i)*PVY->at(i));

     if( PVIsFake->at(i) ) continue;
     if( !(PVNDF->at(i) > 4) ) continue;
     if( !(fabs(PVZ->at(i)) <= 24) ) continue;
     if( !(rho <= 2) ) continue;
     v_idx_goodPV.push_back(i);
   }
   
   // loop over muons and select muons passing tight muon ID
   vector<int> v_idx_muon_tight;
   for(size_t i=0; i<MuonPt->size(); i++)
   {
     double normChi2 = (MuonNdof->at(i) != 0 ? MuonChi2->at(i) / MuonNdof->at(i) : MuonChi2->at(i) * 1e6);
     double normTrkChi2 = (MuonTrkNdof->at(i) != 0 ? MuonTrkChi2->at(i) / MuonTrkNdof->at(i) : MuonTrkChi2->at(i) * 1e6);
     double deltaZ = (v_idx_goodPV.size()>0 ? fabs( MuonRefPtZ->at(i) - PVZ->at(v_idx_goodPV[0]) ) : fabs( MuonRefPtZ->at(i) - PVZ->at(0) ));

     if( !MuonIsGlobal->at(i) ) continue;
     if( !(MuonPt->at(i) > 5) ) continue;
     if( !(fabs(MuonEta->at(i)) < 2.4) ) continue;
     if( !(MuonNValidTrackerHits->at(i) > 10) ) continue;
     if( !(MuonNValidPixelHits->at(i) > 1) ) continue;
     if( !(MuonNLostTrackerHitsOut->at(i) < 3) ) continue;
     if( !(MuonNValidMuonHits->at(i) > 0) ) continue;
     if( !(MuonNMatches->at(i) > 1) ) continue;
     if( !(normChi2 < 10) ) continue;
     if( !(normTrkChi2 < 10) ) continue;
     if( !(deltaZ < 1) ) continue;
     v_idx_muon_tight.push_back(i);
   }

   int passEEAnomJetFilter = 1;
   if( v_idx_pfjet_JetID.size() > 0 )
   {
     if( PFJetPt->at(v_idx_pfjet_JetID[0]) > 15000 ) passEEAnomJetFilter = 0;
   }

   int nBTaggedJets = 0;
   vector<double> scaleFactors;
   int nHeavyFlavorJets = 0;
   int nBTaggedHeavyFlavorJets = 0;
   int nMuons = 0;

   if( v_idx_pfjet_JetID.size() >= 2 )
   {
     // jet, GenParticle, and muon 4-vectors
     TLorentzVector v_j, v_gp, v_m;

     // loop over two leading jets
     for(size_t i=0; i<2; ++i)
     {
       int partonFlavor = 0;

       // set jet 4-vector
       v_j.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_JetID[i]),PFJetEta->at(v_idx_pfjet_JetID[i]),PFJetPhi->at(v_idx_pfjet_JetID[i]),PFJetE->at(v_idx_pfjet_JetID[i]));

       if( !iEvent.isRealData() )
       {
         if( matchingType==0 ) // parton-based matching
         {
           partonFlavor = abs(PFJetPartonFlavor->at(v_idx_pfjet_JetID[i]));
          
           if( abs(PFJetPartonFlavor->at(v_idx_pfjet_JetID[i]))==5 ) ++nHeavyFlavorJets;
         }
         else if( matchingType!=0 ) // hadron-based matching
         {
           double minDeltaR = 999.;

           // loop over GenParticles
           for(size_t j=0; j<GenParticlePt->size(); ++j)
           {
             int pdgID = abs(GenParticlePdgId->at(j));

             if( pdgID==511 || pdgID==521 || pdgID==531 || pdgID==541 || pdgID==5122 || pdgID==5132 || pdgID==5232 || pdgID==5332
                 || pdgID==411 || pdgID==421 || pdgID==431 || pdgID==4122 || pdgID==4132 || pdgID==4232 || pdgID==4332 )
             {
               // set GenParticle 4-vector
               v_gp.SetPtEtaPhiE(GenParticlePt->at(j),GenParticleEta->at(j),GenParticlePhi->at(j),GenParticleE->at(j));
               double deltaR = v_j.DeltaR(v_gp);

               if( deltaR < minDeltaR ) minDeltaR = deltaR;
             }
           }

           if( minDeltaR < matchingRadius )
           {
             ++nHeavyFlavorJets;
             partonFlavor = 5; // This is not necessarily true since hadron-based matching cannot distinguish b- and c-jets. However, since the same scale factors are applied to b- and c-jets, this is a reasonable default.
           }
         }
       }

       if( (btagger==0 && PFJetTCHE->at(v_idx_pfjet_JetID[i])>getPreCutValue1("TCHEL_WP")) ||
           (btagger==1 && PFJetTCHE->at(v_idx_pfjet_JetID[i])>getPreCutValue1("TCHEM_WP")) ||
           (btagger==2 && PFJetSSVHE->at(v_idx_pfjet_JetID[i])>getPreCutValue1("SSVHEM_WP")) ||
           (btagger==3 && PFJetTCHP->at(v_idx_pfjet_JetID[i])>getPreCutValue1("TCHPT_WP")) ||
           (btagger==4 && PFJetSSVHP->at(v_idx_pfjet_JetID[i])>getPreCutValue1("SSVHPT_WP")) )
       {
         ++nBTaggedJets;
         if( partonFlavor==5 ) ++nBTaggedHeavyFlavorJets;
         // if MC, get b-tag scale factor
         if( !iEvent.isRealData() )
         {
           if( useFixedSFs ) scaleFactors.push_back(sfCalculator.scaleFactor(partonFlavor,btagger));
           else scaleFactors.push_back(sfCalculator.scaleFactor(partonFlavor,PFJetPt->at(v_idx_pfjet_JetID[i]),PFJetEta->at(v_idx_pfjet_JetID[i]),btagger));
         }
       }

       // loop over all tight muons and find those that are inside the jet (DeltaR<0.4)
       for(size_t j=0; j<v_idx_muon_tight.size(); ++j)
       {
         // set muon 4-vector
         v_m.SetPtEtaPhiM(MuonPt->at(v_idx_muon_tight[j]),MuonEta->at(v_idx_muon_tight[j]),MuonPhi->at(v_idx_muon_tight[j]),0);
         if( v_j.DeltaR(v_m) < 0.4 ) ++nMuons;
       }
     }
   }

   // apply heavy flavor k-factor
   if( !iEvent.isRealData() && useHFkFactor )
   {
     if( nHeavyFlavorJets==0 ) eventWeight *= (1+(1-getPreCutValue1("HFkFactor"))*(getPreCutValue1("HFFraction")/(1-getPreCutValue1("HFFraction"))));
     else                      eventWeight *= getPreCutValue1("HFkFactor");
   }

   double pretagWeight = eventWeight;
   double tagWeight = pretagWeight;
   
   // Set the evaluation of the cuts to false and clear the variable values and filled status
   resetCuts();
   
   fillVariableWithValue("PassHBHENoiseFilter", ( *passHBHENoiseFilter ? 1 : 0 ), pretagWeight );
   fillVariableWithValue("PassBeamHaloFltTight", ( *passBeamHaloFilterTight ? 1 : 0 ), pretagWeight );
   fillVariableWithValue("PassTrackingFailure", ( *passTrackingFailure ? 1 : 0 ), pretagWeight );
   fillVariableWithValue("PassEcalMskCellDRFlt", ( *passEcalMaskedCellDRFilter ? 1 : 0 ), pretagWeight );
   fillVariableWithValue("PassCaloBndDRFlt", ( *passCaloBoundaryDRFilter ? 1 : 0 ), pretagWeight );
   fillVariableWithValue("PassEEAnomJetFilter", passEEAnomJetFilter, pretagWeight );

   fillVariableWithValue("METoSumET_pretag", MET->front()/SumET->front(), pretagWeight );

   fillVariableWithValue("nJets_all", PFJetPt->size(), pretagWeight);
   fillVariableWithValue("nJets_JetID", v_idx_pfjet_JetID.size(), pretagWeight);

   fillVariableWithValue("nGoodVertices_pretag", v_idx_goodPV.size(), pretagWeight );

   if( v_idx_pfjet_JetID.size() >= 1 )
   {
     fillVariableWithValue( "PhiJ1_pretag", PFJetPhi->at(v_idx_pfjet_JetID[0]), pretagWeight );
     fillVariableWithValue( "absEtaJ1", fabs( PFJetEta->at(v_idx_pfjet_JetID[0]) ), pretagWeight );
     fillVariableWithValue( "EtaJ1_pretag", PFJetEta->at(v_idx_pfjet_JetID[0]), pretagWeight );
     fillVariableWithValue( "PtJ1_cut", PFJetPt->at(v_idx_pfjet_JetID[0]), pretagWeight );
     fillVariableWithValue( "PtJ1_pretag", getVariableValue("PtJ1_cut"), pretagWeight );
   }
   if( v_idx_pfjet_JetID.size() >= 2 )
   {
     fillVariableWithValue( "PhiJ2_pretag", PFJetPhi->at(v_idx_pfjet_JetID[1]), pretagWeight );
     fillVariableWithValue( "absEtaJ2", fabs( PFJetEta->at(v_idx_pfjet_JetID[1]) ), pretagWeight );
     fillVariableWithValue( "EtaJ2_pretag", PFJetEta->at(v_idx_pfjet_JetID[1]), pretagWeight );
     fillVariableWithValue( "PtJ2_cut", PFJetPt->at(v_idx_pfjet_JetID[1]), pretagWeight );
     fillVariableWithValue( "PtJ2_pretag", getVariableValue("PtJ2_cut"), pretagWeight );
    
     // calculate |DeltaEta(j1,j2)|
     fillVariableWithValue( "absDeltaEtaJ1J2", fabs( PFJetEta->at(v_idx_pfjet_JetID[0]) - PFJetEta->at(v_idx_pfjet_JetID[1]) ), pretagWeight );
     fillVariableWithValue( "DeltaEtaJ1J2_pretag", getVariableValue("absDeltaEtaJ1J2"), pretagWeight );
     
     TLorentzVector v_j1j2, v_j1, v_j2;
     v_j1.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_JetID[0]),PFJetEta->at(v_idx_pfjet_JetID[0]),PFJetPhi->at(v_idx_pfjet_JetID[0]),PFJetE->at(v_idx_pfjet_JetID[0]));
     v_j2.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_JetID[1]),PFJetEta->at(v_idx_pfjet_JetID[1]),PFJetPhi->at(v_idx_pfjet_JetID[1]),PFJetE->at(v_idx_pfjet_JetID[1]));
     // calculate M_j1j2
     v_j1j2 = v_j1 + v_j2;
     
     fillVariableWithValue( "DijetMassThreshold", v_j1j2.M(), pretagWeight );
     
     fillVariableWithValue( "absDeltaPhiJ1J2", fabs( v_j1.DeltaPhi(v_j2) ), pretagWeight );
     fillVariableWithValue( "DeltaPhiJ1J2_pretag", getVariableValue("absDeltaPhiJ1J2"), pretagWeight );
     
     fillVariableWithValue( "DijetMass_pretag", getVariableValue("DijetMassThreshold"), pretagWeight );

     FillUserTH2D("h2_EtaJ2_vs_EtaJ1", getVariableValue("EtaJ1_pretag"), getVariableValue("EtaJ2_pretag"), pretagWeight);

     
     fillVariableWithValue( "nMuons_pretag", nMuons, pretagWeight );

     // in MC, if the b-tag scale factor reweighting is enabled, apply the b-tag event weight
     if( !iEvent.isRealData() && doSFReweighting )
     {
       tagWeight = eventWeight*bTagEventWeight(scaleFactors,0);
       fillVariableWithValue( "nJets_btag", 0, tagWeight );
     }
     else
     {
       fillVariableWithValue( "nJets_btag", nBTaggedJets, tagWeight );
     }
     fillVariableWithValue( "nGoodVertices", getVariableValue("nGoodVertices_pretag"), tagWeight );
     fillVariableWithValue( "PhiJ1", getVariableValue("PhiJ1_pretag"), tagWeight );
     fillVariableWithValue( "PhiJ2", getVariableValue("PhiJ2_pretag"), tagWeight );
     fillVariableWithValue( "DeltaPhiJ1J2", getVariableValue("DeltaPhiJ1J2_pretag"), tagWeight );
     fillVariableWithValue( "EtaJ1", getVariableValue("EtaJ1_pretag"), tagWeight );
     fillVariableWithValue( "EtaJ2", getVariableValue("EtaJ2_pretag"), tagWeight );
     fillVariableWithValue( "DeltaEtaJ1J2", getVariableValue("DeltaEtaJ1J2_pretag"), tagWeight );
     fillVariableWithValue( "PtJ1", getVariableValue("PtJ1_pretag"), tagWeight );
     fillVariableWithValue( "PtJ2", getVariableValue("PtJ2_pretag"), tagWeight );
     fillVariableWithValue( "DijetMass", getVariableValue("DijetMass_pretag"), tagWeight );
     fillVariableWithValue( "nMuons", getVariableValue("nMuons_pretag"), tagWeight );
     fillVariableWithValue( "METoSumET", getVariableValue("METoSumET_pretag"), tagWeight );
   }
   
   // Evaluate cuts (but do not apply them)
   evaluateCuts();

   
   if(passedAllPreviousCuts("DijetMass_pretag"))
   {
     FillUserTH1D("h1_J1J2PartonFlavor", abs( PFJetPartonFlavor->at(v_idx_pfjet_JetID[0]) ), pretagWeight );
     FillUserTH1D("h1_J1J2PartonFlavor", abs( PFJetPartonFlavor->at(v_idx_pfjet_JetID[1]) ), pretagWeight );
     
     if( nHeavyFlavorJets==2 )
       FillUserTH1D("h1_J1J2HeavyFlavor", 1, 2.*pretagWeight );
     else if( nHeavyFlavorJets==1 )
     {
       FillUserTH1D("h1_J1J2HeavyFlavor", 1, pretagWeight );
       FillUserTH1D("h1_J1J2HeavyFlavor", 0, pretagWeight );
     }
     else
       FillUserTH1D("h1_J1J2HeavyFlavor", 0, 2.*pretagWeight );
   }
   if(passedAllPreviousCuts("nMuons_pretag")) FillUserTH1D("h1_nMuons_vs_DijetMass_pretag", getVariableValue("DijetMass_pretag"), double(nMuons)*pretagWeight );
   if(passedAllPreviousCuts("nMuons")) FillUserTH1D("h1_nMuons_vs_DijetMass", getVariableValue("DijetMass_pretag"), double(nMuons)*tagWeight );

   // select only those events that pass the full selection
   if( passedCut("all") ) ret = true;

   // cuts to be reset in MC when b-tag event reweighting is enabled
   vector<string> cutNames;
   cutNames.push_back("nJets_btag"); cutNames.push_back("nGoodVertices"); cutNames.push_back("PhiJ1"); cutNames.push_back("PhiJ2"); cutNames.push_back("DeltaPhiJ1J2");
   cutNames.push_back("EtaJ1"); cutNames.push_back("EtaJ2"); cutNames.push_back("DeltaEtaJ1J2"); cutNames.push_back("PtJ1"); cutNames.push_back("PtJ2");
   cutNames.push_back("DijetMass"); cutNames.push_back("nMuons"); cutNames.push_back("METoSumET");

   if( v_idx_pfjet_JetID.size() >= 2 && !iEvent.isRealData() && doSFReweighting )
   {
     for( int nbtags=1; nbtags<=2; ++nbtags )
     {
       // Set the evaluation of the cuts to false and clear the variable values and filled status for a subset of cuts defined by the cutNames vector
       resetCuts(cutNames);

       tagWeight = eventWeight*bTagEventWeight(scaleFactors,nbtags);

       fillVariableWithValue( "nJets_btag", nbtags, tagWeight );
       fillVariableWithValue( "nGoodVertices", getVariableValue("nGoodVertices_pretag"), tagWeight );
       fillVariableWithValue( "PhiJ1", getVariableValue("PhiJ1_pretag"), tagWeight );
       fillVariableWithValue( "PhiJ2", getVariableValue("PhiJ2_pretag"), tagWeight );
       fillVariableWithValue( "DeltaPhiJ1J2", getVariableValue("DeltaPhiJ1J2_pretag"), tagWeight );
       fillVariableWithValue( "EtaJ1", getVariableValue("EtaJ1_pretag"), tagWeight );
       fillVariableWithValue( "EtaJ2", getVariableValue("EtaJ2_pretag"), tagWeight );
       fillVariableWithValue( "DeltaEtaJ1J2", getVariableValue("DeltaEtaJ1J2_pretag"), tagWeight );
       fillVariableWithValue( "PtJ1", getVariableValue("PtJ1_pretag"), tagWeight );
       fillVariableWithValue( "PtJ2", getVariableValue("PtJ2_pretag"), tagWeight );
       fillVariableWithValue( "DijetMass", getVariableValue("DijetMass_pretag"), tagWeight );
       fillVariableWithValue( "nMuons", getVariableValue("nMuons_pretag"), tagWeight );
       fillVariableWithValue( "METoSumET", getVariableValue("METoSumET_pretag"), tagWeight );

       // Evaluate cuts (but do not apply them)
       evaluateCuts();

       if(passedAllPreviousCuts("nMuons")) FillUserTH1D("h1_nMuons_vs_DijetMass", getVariableValue("DijetMass_pretag"), double(nMuons)*tagWeight );

       // select only those events that pass the full selection
       if( passedCut("all") ) ret = true;
     }
   }
   
   // fill EventBin histograms
   if( doEventBins && passedAllPreviousCuts("DijetMass_pretag") )
   {
     if( nMuons==0 && max(getVariableValue("absEtaJ1"),getVariableValue("absEtaJ2"))<1.2 )
     {
       if( iEvent.isRealData() )
       {
         if( nBTaggedJets==0 ) FillUserTH2D("h2_n0_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
         if( nBTaggedJets==1 ) FillUserTH2D("h2_n1_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
         if( nBTaggedJets==2 ) FillUserTH2D("h2_n2_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
       }
       else if( !iEvent.isRealData() )
       {
         if( doSFReweighting )
         {
           if( nHeavyFlavorJets==2 )
           {
             FillUserTH2D("h2_N22_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N21_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N20_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
           else if( nHeavyFlavorJets==1 )
           {
             FillUserTH2D("h2_N12_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N11_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N10_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
           else if( nHeavyFlavorJets==0 )
           {
             FillUserTH2D("h2_N02_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N01_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N00_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
         }
         else
         {
           if( nHeavyFlavorJets==2 )
           {
             if( nBTaggedJets==2 )      FillUserTH2D("h2_N220_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 ) FillUserTH2D("h2_N210_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else                       FillUserTH2D("h2_N200_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
           else if( nHeavyFlavorJets==1 )
           {
             if( nBTaggedJets==2 )                                    FillUserTH2D("h2_N111_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 && nBTaggedHeavyFlavorJets==1 ) FillUserTH2D("h2_N110_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 && nBTaggedHeavyFlavorJets==0 ) FillUserTH2D("h2_N101_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==0 )                               FillUserTH2D("h2_N100_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
           else if( nHeavyFlavorJets==0 )
           {
             if( nBTaggedJets==2 )      FillUserTH2D("h2_N002_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 ) FillUserTH2D("h2_N001_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else                       FillUserTH2D("h2_N000_0mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
         }
       }
     }
     else if( nMuons==0 && max(getVariableValue("absEtaJ1"),getVariableValue("absEtaJ2"))>=1.2 )
     {
       if( iEvent.isRealData() )
       {
         if( nBTaggedJets==0 ) FillUserTH2D("h2_n0_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
         if( nBTaggedJets==1 ) FillUserTH2D("h2_n1_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
         if( nBTaggedJets==2 ) FillUserTH2D("h2_n2_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
       }
       else if( !iEvent.isRealData() )
       {
         if( doSFReweighting )
         {
           if( nHeavyFlavorJets==2 )
           {
             FillUserTH2D("h2_N22_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N21_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N20_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
           else if( nHeavyFlavorJets==1 )
           {
             FillUserTH2D("h2_N12_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N11_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N10_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
           else if( nHeavyFlavorJets==0 )
           {
             FillUserTH2D("h2_N02_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N01_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N00_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
         }
         else
         {
           if( nHeavyFlavorJets==2 )
           {
             if( nBTaggedJets==2 )      FillUserTH2D("h2_N220_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 ) FillUserTH2D("h2_N210_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else                       FillUserTH2D("h2_N200_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
           else if( nHeavyFlavorJets==1 )
           {
             if( nBTaggedJets==2 )                                    FillUserTH2D("h2_N111_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 && nBTaggedHeavyFlavorJets==1 ) FillUserTH2D("h2_N110_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 && nBTaggedHeavyFlavorJets==0 ) FillUserTH2D("h2_N101_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==0 )                               FillUserTH2D("h2_N100_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
           else if( nHeavyFlavorJets==0 )
           {
             if( nBTaggedJets==2 )      FillUserTH2D("h2_N002_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 ) FillUserTH2D("h2_N001_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else                       FillUserTH2D("h2_N000_0mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
         }
       }
     }
     else if( nMuons>=1 && max(getVariableValue("absEtaJ1"),getVariableValue("absEtaJ2"))<1.2 )
     {
       if( iEvent.isRealData() )
       {
         if( nBTaggedJets==0 ) FillUserTH2D("h2_n0_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
         if( nBTaggedJets==1 ) FillUserTH2D("h2_n1_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
         if( nBTaggedJets==2 ) FillUserTH2D("h2_n2_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
       }
       else if( !iEvent.isRealData() )
       {
         if( doSFReweighting )
         {
           if( nHeavyFlavorJets==2 )
           {
             FillUserTH2D("h2_N22_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N21_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N20_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
           else if( nHeavyFlavorJets==1 )
           {
             FillUserTH2D("h2_N12_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N11_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N10_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
           else if( nHeavyFlavorJets==0 )
           {
             FillUserTH2D("h2_N02_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N01_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N00_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
         }
         else
         {
           if( nHeavyFlavorJets==2 )
           {
             if( nBTaggedJets==2 )      FillUserTH2D("h2_N220_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 ) FillUserTH2D("h2_N210_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else                       FillUserTH2D("h2_N200_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
           else if( nHeavyFlavorJets==1 )
           {
             if( nBTaggedJets==2 )                                    FillUserTH2D("h2_N111_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 && nBTaggedHeavyFlavorJets==1 ) FillUserTH2D("h2_N110_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 && nBTaggedHeavyFlavorJets==0 ) FillUserTH2D("h2_N101_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==0 )                               FillUserTH2D("h2_N100_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
           else if( nHeavyFlavorJets==0 )
           {
             if( nBTaggedJets==2 )      FillUserTH2D("h2_N002_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 ) FillUserTH2D("h2_N001_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else                       FillUserTH2D("h2_N000_ge1mu_maxEta_lt_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
         }
       }
     }
     else if( nMuons>=1 && max(getVariableValue("absEtaJ1"),getVariableValue("absEtaJ2"))>=1.2 )
     {
       if( iEvent.isRealData() )
       {
         if( nBTaggedJets==0 ) FillUserTH2D("h2_n0_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
         if( nBTaggedJets==1 ) FillUserTH2D("h2_n1_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
         if( nBTaggedJets==2 ) FillUserTH2D("h2_n2_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
       }
       else if( !iEvent.isRealData() )
       {
         if( doSFReweighting )
         {
           if( nHeavyFlavorJets==2 )
           {
             FillUserTH2D("h2_N22_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N21_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N20_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
           else if( nHeavyFlavorJets==1 )
           {
             FillUserTH2D("h2_N12_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N11_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N10_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
           else if( nHeavyFlavorJets==0 )
           {
             FillUserTH2D("h2_N02_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,2));
             FillUserTH2D("h2_N01_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,1));
             FillUserTH2D("h2_N00_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),eventWeight*bTagEventWeight(scaleFactors,0));
           }
         }
         else
         {
           if( nHeavyFlavorJets==2 )
           {
             if( nBTaggedJets==2 )      FillUserTH2D("h2_N220_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 ) FillUserTH2D("h2_N210_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else                       FillUserTH2D("h2_N200_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
           else if( nHeavyFlavorJets==1 )
           {
             if( nBTaggedJets==2 )                                    FillUserTH2D("h2_N111_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 && nBTaggedHeavyFlavorJets==1 ) FillUserTH2D("h2_N110_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 && nBTaggedHeavyFlavorJets==0 ) FillUserTH2D("h2_N101_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==0 )                               FillUserTH2D("h2_N100_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
           else if( nHeavyFlavorJets==0 )
           {
             if( nBTaggedJets==2 )      FillUserTH2D("h2_N002_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else if( nBTaggedJets==1 ) FillUserTH2D("h2_N001_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
             else                       FillUserTH2D("h2_N000_ge1mu_maxEta_ge_1p2",getVariableValue("DijetMass_pretag"),getVariableValue("nGoodVertices_pretag"),pretagWeight);
           }
         }
       }
     }
   }

   // ##################### Event Printout - START #####################

   // METoverSumET_gt_0.45 event printout
   if( doEventPrintout && iEvent.isRealData() && passedAllPreviousCuts("DijetMass_pretag") && getVariableValue("METoSumET_pretag")>METoSumET_cut )
   {
     string category = "METoverSumET_gt_0.45";
     
     cout << category << ": ----------- START ------------" << endl;
     cout << category << ": Run, lumi, event: "<< iEvent.id().run() << ", "
                                               << iEvent.luminosityBlock() << ", "
                                               << iEvent.id().event() << endl;

     // loop over the two leading ak7 PFJets that pass jet ID requirements
     for (size_t i=0; i<2; ++i)
     {
       cout << category << ": PassJetID PFJet "<< i << " Pt, PtRaw, eta, phi: " << PFJetPt->at(v_idx_pfjet_JetID[i]) << ", "
                                                                                  << PFJetPtRaw->at(v_idx_pfjet_JetID[i]) << ", "
                                                                                  << PFJetEta->at(v_idx_pfjet_JetID[i]) << ", "
                                                                                  << PFJetPhi->at(v_idx_pfjet_JetID[i]) << endl;
     }
     cout << category << ": |DeltaEtaJ1J2|: "<< getVariableValue("absDeltaEtaJ1J2") << endl;
     cout << category << ": |DeltaPhiJ1J2|: "<< getVariableValue("absDeltaPhiJ1J2") << endl;
     cout << category << ": MET: "<< MET->front() << endl;
     cout << category << ": MET/SumET: "<< getVariableValue("METoSumET_pretag") << endl;
     cout << category << ": Dijet Mass: "<< getVariableValue("DijetMass_pretag") << endl;
     cout << category << ": ------------ END -------------" << endl;
   }

   // DeltaPhiJ1J2_lt_1.5 event printout
   if( doEventPrintout && iEvent.isRealData() && passedAllPreviousCuts("DijetMass_pretag") && getVariableValue("absDeltaPhiJ1J2")<DeltaPhiJ1J2_cut )
   {
     string category = "DeltaPhiJ1J2_lt_1.5";
    
     cout << category << ": ----------- START ------------" << endl;
     cout << category << ": Run, lumi, event: "<< iEvent.id().run() << ", "
                                               << iEvent.luminosityBlock() << ", "
                                               << iEvent.id().event() << endl;

     // loop over the two leading ak7 PFJets that pass jet ID requirements
     for (size_t i=0; i<2; ++i)
     {
       cout << category << ": PassJetID PFJet "<< i << " Pt, PtRaw, eta, phi: " << PFJetPt->at(v_idx_pfjet_JetID[i]) << ", "
                                                                                  << PFJetPtRaw->at(v_idx_pfjet_JetID[i]) << ", "
                                                                                  << PFJetEta->at(v_idx_pfjet_JetID[i]) << ", "
                                                                                  << PFJetPhi->at(v_idx_pfjet_JetID[i]) << endl;
     }
     cout << category << ": |DeltaEtaJ1J2|: "<< getVariableValue("absDeltaEtaJ1J2") << endl;
     cout << category << ": |DeltaPhiJ1J2|: "<< getVariableValue("absDeltaPhiJ1J2") << endl;
     cout << category << ": MET: "<< MET->front() << endl;
     cout << category << ": MET/SumET: "<< getVariableValue("METoSumET_pretag") << endl;
     cout << category << ": Dijet Mass: "<< getVariableValue("DijetMass_pretag") << endl;
     cout << category << ": ------------ END -------------" << endl;
   }

   // ##################### Event Printout - END #######################
   
   //############################# User's code ends here #################################
   //#####################################################################################
   
   // increment event counters
   eventCountBeforeWeight++;
   eventCount += eventWeight;
   
   return ret;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool
MyAnalyzer::endLuminosityBlock(edm::LuminosityBlock& iLumi, edm::EventSetup const& iSetup)
{
  if ( skimWasMade_ )
  {
      edm::Handle<edm::MergeableCounter> eventCounter;

      if (iLumi.getByLabel(eventCounterInputTag_, eventCounter) && eventCounter.isValid())
      {
          NEvtTotBeforeWeight_ += (double) eventCounter->value;
      }
      else
      {
          edm::LogError("MyAnalyzer::endLuminosityBlock") << "Can't get the product " << eventCounterInputTag_ <<". Please make sure the skimWasMade and eventCounterInputTag parameters were set correctly.";
          exit(1);
      }
  }

  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool
MyAnalyzer::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called once each job just after ending the event loop  ------------
void
MyAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ------------ method that calculates the event weight based on the number of b-tagged jets in MC and the expected number of b-tags among the two leading jets  ------------
double
MyAnalyzer::bTagEventWeight(const vector<double>& SFsForBTaggedJets, const unsigned int nBTags)
{
  if( SFsForBTaggedJets.size() > 2 )
  {
    edm::LogError("MyAnalyzer::bTagEventWeight") << "Only two leading jets are considered. Hence, the number of b-tagged jets cannot exceed 2.";
    exit(1);
  }
  if( nBTags > 2 )
  {
    edm::LogError("MyAnalyzer::bTagEventWeight") << "Only two leading jets are considered. Hence, the number of b-tags cannot exceed 2.";
    exit(1);
  }
  /*
    ##################################################################
    Event weight matrix:
    ------------------------------------------------------------------
    nBTags\b-tagged jets  |    0        1             2
    ------------------------------------------------------------------
      0                   |    1      1-SF      (1-SF1)(1-SF2)
                          |
      1                   |    0       SF    SF1(1-SF2)+(1-SF1)SF2
                          |
      2                   |    0        0           SF1SF2
    ##################################################################
  */
  
  if( nBTags > SFsForBTaggedJets.size() ) return 0;

  if( nBTags==0 && SFsForBTaggedJets.size()==0 ) return 1;

  double weight = 0;

  if( SFsForBTaggedJets.size()==1 )
  {
    double SF = SFsForBTaggedJets.at(0);

    for( unsigned int i=0; i<=1; ++i )
    {
      if( i != nBTags ) continue;

      weight += pow(SF,i)*pow(1-SF,1-i);
    }
  }
  else if( SFsForBTaggedJets.size()==2 )
  {
    double SF1 = SFsForBTaggedJets.at(0);
    double SF2 = SFsForBTaggedJets.at(1);
    
    for( unsigned int i=0; i<=1; ++i )
    {
      for( unsigned int j=0; j<=1; ++j )
      {
        if( (i+j) != nBTags ) continue;

        weight += pow(SF1,i)*pow(1-SF1,1-i)*pow(SF2,j)*pow(1-SF2,1-j);
      }
    }
  }
  return weight;
}

// BTagScaleFactorCalculator constructor
BTagScaleFactorCalculator::BTagScaleFactorCalculator()
{
  TCHEL_SFb_ = 1.;
  TCHEL_SFl_ = 1.;
  TCHPT_SFb_ = 1.;
  TCHPT_SFl_ = 1.;
  SSVHPT_SFb_ = 1.;
  SSVHPT_SFl_ = 1.;
  TCHEL_SFb_0to2p4 = new TF1("TCHEL_SFb_0to2p4","0.603913*((1.+(0.286361*x))/(1.+(0.170474*x)))", 30.,670.);
  TCHEL_SFl_0to2p4 = new TF1("TCHEL_SFl_0to2p4","(1.10649*((1+(-9.00297e-05*x))+(2.32185e-07*(x*x))))+(-4.04925e-10*(x*(x*(x/(1+(-0.00051036*x))))))", 20.,670.);
  TCHEL_SFl_0to0p5 = new TF1("TCHEL_SFl_0to0p5","(1.13615*((1+(-0.00119852*x))+(1.17888e-05*(x*x))))+(-9.8581e-08*(x*(x*(x/(1+(0.00689317*x))))))", 20.,670.);
  TCHEL_SFl_0p5to1p0 = new TF1("TCHEL_SFl_0p5to1p0","(1.13277*((1+(-0.00084146*x))+(3.80313e-06*(x*x))))+(-8.75061e-09*(x*(x*(x/(1+(0.00118695*x))))))", 20.,670.);
  TCHEL_SFl_1p0to1p5 = new TF1("TCHEL_SFl_1p0to1p5","(1.17163*((1+(-0.000828475*x))+(3.0769e-06*(x*x))))+(-4.692e-09*(x*(x*(x/(1+(0.000337759*x))))))", 20.,670.);
  TCHEL_SFl_1p5to2p4 = new TF1("TCHEL_SFl_1p5to2p4","(1.14554*((1+(-0.000128043*x))+(4.10899e-07*(x*x))))+(-2.07565e-10*(x*(x*(x/(1+(-0.00118618*x))))))", 20.,670.);
}

// ------------ method that initializes the BTagScaleFactorCalculator class  ------------
void
BTagScaleFactorCalculator::init(const double TCHEL_SFb, const double TCHEL_SFl, const double TCHPT_SFb, const double TCHPT_SFl, const double SSVHPT_SFb, const double SSVHPT_SFl)
{
  TCHEL_SFb_ = TCHEL_SFb;
  TCHEL_SFl_ = TCHEL_SFl;
  TCHPT_SFb_ = TCHPT_SFb;
  TCHPT_SFl_ = TCHPT_SFl;
  SSVHPT_SFb_ = SSVHPT_SFb;
  SSVHPT_SFl_ = SSVHPT_SFl;
}

// ------------ method that returns the b-tag efficiency scale factor  ------------
double
BTagScaleFactorCalculator::scaleFactor(const int partonFlavor, const int btagger)
{
  if( partonFlavor==5 || partonFlavor==4 )
  {
    if(btagger==0)      return TCHEL_SFb_;
    else if(btagger==3) return TCHPT_SFb_;
    else if(btagger==4) return SSVHPT_SFb_;
    else                return 1.;
  }
  else
  {
    if(btagger==0)      return TCHEL_SFl_;
    else if(btagger==3) return TCHPT_SFl_;
    else if(btagger==4) return SSVHPT_SFl_;
    else                return 1.;
  }
}


// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor  ------------
double
BTagScaleFactorCalculator::scaleFactor(const int partonFlavor, const double jetPt, const double jetEta, const int btagger)
{
  if( partonFlavor==5 || partonFlavor==4 )
  {
    if(btagger==0)  return scaleFactorBC_TCHEL(jetPt,jetEta);
    else            return 1.;
  }
  else
  {
    if(btagger==0)  return scaleFactorUDSG_TCHEL(jetPt,jetEta);
    else            return 1.;
  }
}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for b- and c-jets and TCHEL tagger  ------------
double
BTagScaleFactorCalculator::scaleFactorBC_TCHEL(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  // for scale factor extrapolation
  if(Pt<30) Pt = 30;
  if(Pt>670) Pt = 670;

  return TCHEL_SFb_0to2p4->Eval(Pt);

}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for light flavor jets and TCHEL tagger ------------
double
BTagScaleFactorCalculator::scaleFactorUDSG_TCHEL(const double jetPt, const double jetEta)
{
  double SF = 1.;
  double Pt = jetPt;
  double eta = fabs(jetEta);
  // for scale factor extrapolation
  if(Pt<20) Pt = 20;

  if(eta<0.5)
  {
    if( Pt>670 ) SF = TCHEL_SFl_0to2p4->Eval(670);
    else         SF = TCHEL_SFl_0to0p5->Eval(Pt);
  }
  else if(eta>=0.5 && eta<1.)
  {
    if( Pt>670 ) SF = TCHEL_SFl_0to2p4->Eval(670);
    else         SF = TCHEL_SFl_0p5to1p0->Eval(Pt);
  }
  else if(eta>=1. && eta<1.5)
  {
    if( Pt>670 ) SF = TCHEL_SFl_0to2p4->Eval(670);
    else         SF = TCHEL_SFl_1p0to1p5->Eval(Pt);
  }
  else
  {
    if( Pt>670 ) SF = TCHEL_SFl_0to2p4->Eval(670);
    else         SF = TCHEL_SFl_1p5to2p4->Eval(Pt);
  }

  return SF;
}


DEFINE_FWK_MODULE(MyAnalyzer);
