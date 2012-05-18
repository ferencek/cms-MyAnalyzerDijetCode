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
// $Id: MyAnalyzer_MainAnalysis_DijetBBTag_2011.cc,v 1.18 2012/04/13 19:54:15 ferencek Exp $
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
#include <TH1D.h>
#include <TLorentzVector.h>


using namespace std;

//
// class declaration
//

class BTagScaleFactorCalculator
{
 public:
   BTagScaleFactorCalculator();
   void init(const double SFb_shift, const double SFl_shift, const double TCHEL_SFb, const double TCHEL_SFl, const double TCHPT_SFb, const double TCHPT_SFl, const double SSVHPT_SFb, const double SSVHPT_SFl);
   double scaleFactor(const int partonFlavor, const int btagger);
   double scaleFactor(const int partonFlavor, const double jetPt, const double jetEta, const int btagger);

   double scaleFactorB_TCHEL(const double jetPt, const double jetEta);
   double scaleFactorC_TCHEL(const double jetPt, const double jetEta);
   double scaleFactorUDSG_TCHEL(const double jetPt, const double jetEta);

   double scaleFactorB_CSVL(const double jetPt, const double jetEta);
   double scaleFactorC_CSVL(const double jetPt, const double jetEta);
   double scaleFactorUDSG_CSVL(const double jetPt, const double jetEta);

   double scaleFactorB_CSVM(const double jetPt, const double jetEta);
   double scaleFactorC_CSVM(const double jetPt, const double jetEta);
   double scaleFactorUDSG_CSVM(const double jetPt, const double jetEta);

 private:
   double SFb_shift_;
   double SFl_shift_;
   double TCHEL_SFb_;
   double TCHEL_SFl_;
   double TCHPT_SFb_;
   double TCHPT_SFl_;
   double SSVHPT_SFb_;
   double SSVHPT_SFl_;

   // TCHEL
   TF1  *TCHEL_SFb_0to2p4;
   TH1D *TCHEL_SFb_errors;

   TF1 *TCHEL_SFl_0to2p4;
   TF1 *TCHEL_SFl_0to0p5;
   TF1 *TCHEL_SFl_0p5to1p0;
   TF1 *TCHEL_SFl_1p0to1p5;
   TF1 *TCHEL_SFl_1p5to2p4;

   TF1 *TCHEL_SFl_0to2p4_min;
   TF1 *TCHEL_SFl_0to0p5_min;
   TF1 *TCHEL_SFl_0p5to1p0_min;
   TF1 *TCHEL_SFl_1p0to1p5_min;
   TF1 *TCHEL_SFl_1p5to2p4_min;

   TF1 *TCHEL_SFl_0to2p4_max;
   TF1 *TCHEL_SFl_0to0p5_max;
   TF1 *TCHEL_SFl_0p5to1p0_max;
   TF1 *TCHEL_SFl_1p0to1p5_max;
   TF1 *TCHEL_SFl_1p5to2p4_max;
   // CSVL
   TF1  *CSVL_SFb_0to2p4;
   TH1D *CSVL_SFb_errors;

   TF1 *CSVL_SFl_0to2p4;
   TF1 *CSVL_SFl_0to0p5;
   TF1 *CSVL_SFl_0p5to1p0;
   TF1 *CSVL_SFl_1p0to1p5;
   TF1 *CSVL_SFl_1p5to2p4;

   TF1 *CSVL_SFl_0to2p4_min;
   TF1 *CSVL_SFl_0to0p5_min;
   TF1 *CSVL_SFl_0p5to1p0_min;
   TF1 *CSVL_SFl_1p0to1p5_min;
   TF1 *CSVL_SFl_1p5to2p4_min;

   TF1 *CSVL_SFl_0to2p4_max;
   TF1 *CSVL_SFl_0to0p5_max;
   TF1 *CSVL_SFl_0p5to1p0_max;
   TF1 *CSVL_SFl_1p0to1p5_max;
   TF1 *CSVL_SFl_1p5to2p4_max;
   // CSVM
   TF1  *CSVM_SFb_0to2p4;
   TH1D *CSVM_SFb_errors;

   TF1 *CSVM_SFl_0to2p4;
   TF1 *CSVM_SFl_0to0p8;
   TF1 *CSVM_SFl_0p8to1p6;
   TF1 *CSVM_SFl_1p6to2p4;

   TF1 *CSVM_SFl_0to2p4_min;
   TF1 *CSVM_SFl_0to0p8_min;
   TF1 *CSVM_SFl_0p8to1p6_min;
   TF1 *CSVM_SFl_1p6to2p4_min;

   TF1 *CSVM_SFl_0to2p4_max;
   TF1 *CSVM_SFl_0to0p8_max;
   TF1 *CSVM_SFl_0p8to1p6_max;
   TF1 *CSVM_SFl_1p6to2p4_max;
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

   sfCalculator.init(getPreCutValue1("SFb_Shift"),getPreCutValue1("SFl_Shift"),getPreCutValue1("TCHEL_SFb"),getPreCutValue1("TCHEL_SFl"),getPreCutValue1("TCHPT_SFb"),getPreCutValue1("TCHPT_SFl"),getPreCutValue1("SSVHPT_SFb"),getPreCutValue1("SSVHPT_SFl"));
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
   CreateUserTH1D("h1_DijetMass_0tag;Dijet Mass [GeV]", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_1tag;Dijet Mass [GeV]", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_2tag;Dijet Mass [GeV]", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   
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
   int useWideJets = int(getPreCutValue1("useWideJets"));
   
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
   edm::Handle<vector<int> > GenParticleMotherIndex;
   iEvent.getByLabel(edm::InputTag("GenParticles:MotherIndex"), GenParticleMotherIndex);
   
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
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:Pt"), PFJetPt_);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:Pt"), PFJetPt_);
   edm::Handle<vector<double> > PFJetPtRaw;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:PtRaw"), PFJetPtRaw);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:PtRaw"), PFJetPtRaw);
   edm::Handle<vector<double> > PFJetEta;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:Eta"), PFJetEta);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:Eta"), PFJetEta);
   edm::Handle<vector<double> > PFJetPhi;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:Phi"), PFJetPhi);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:Phi"), PFJetPhi);
   edm::Handle<vector<double> > PFJetE_;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:Energy"), PFJetE_);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:Energy"), PFJetE_);
   edm::Handle<vector<double> > PFJetUnc;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:JECUnc"), PFJetUnc);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:JECUnc"), PFJetUnc);
   edm::Handle<vector<int> > PFJetPassLooseID;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:PassLooseID"), PFJetPassLooseID);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:PassLooseID"), PFJetPassLooseID);
   edm::Handle<vector<int> > PFJetPassTightID;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:PassTightID"), PFJetPassTightID);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:PassTightID"), PFJetPassTightID);
   edm::Handle<vector<double> > PFJetSSVHE;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:SimpleSecondaryVertexHighEffBTag"), PFJetSSVHE);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:SimpleSecondaryVertexHighEffBTag"), PFJetSSVHE);
   edm::Handle<vector<double> > PFJetSSVHP;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:SimpleSecondaryVertexHighPurBTag"), PFJetSSVHP);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:SimpleSecondaryVertexHighPurBTag"), PFJetSSVHP);
   edm::Handle<vector<double> > PFJetTCHE;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:TrackCountingHighEffBTag"), PFJetTCHE);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:TrackCountingHighEffBTag"), PFJetTCHE);
   edm::Handle<vector<double> > PFJetTCHP;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:TrackCountingHighPurBTag"), PFJetTCHP);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:TrackCountingHighPurBTag"), PFJetTCHP);
   edm::Handle<vector<double> > PFJetJP;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:JetProbabilityBTag"), PFJetJP);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:JetProbabilityBTag"), PFJetJP);
   edm::Handle<vector<double> > PFJetCSV;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:CombinedSecondaryVertexBJetTag"), PFJetCSV);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:CombinedSecondaryVertexBJetTag"), PFJetCSV);
   edm::Handle<vector<int> > PFJetPartonFlavor;
   if(useWideJets) iEvent.getByLabel(edm::InputTag("AK5PFJets:PartonFlavor"), PFJetPartonFlavor);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:PartonFlavor"), PFJetPartonFlavor);
   
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

   // for RS graviton samples
   int nSt3_b_fromRSG = 0, nSt3_c_fromRSG = 0, nSt3_q_fromRSG = 0;

   for(size_t i=0; i<GenParticlePt->size(); i++)
   {
     if( abs(GenParticlePdgId->at(i))==5 && GenParticleStatus->at(i)==3 && GenParticleMotherIndex->at(i)>=0 )
     {
       if( abs(GenParticlePdgId->at(GenParticleMotherIndex->at(i)))==5000039 ) ++nSt3_b_fromRSG;
     }
     if( abs(GenParticlePdgId->at(i))==4 && GenParticleStatus->at(i)==3 && GenParticleMotherIndex->at(i)>=0 )
     {
       if( abs(GenParticlePdgId->at(GenParticleMotherIndex->at(i)))==5000039 ) ++nSt3_c_fromRSG;
     }
     if( abs(GenParticlePdgId->at(i))!=21 && GenParticleStatus->at(i)==3 && GenParticleMotherIndex->at(i)>=0 )
     {
       if( abs(GenParticlePdgId->at(GenParticleMotherIndex->at(i)))==5000039 ) ++nSt3_q_fromRSG;
     }
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
   if( PFJetPt->size() > 0 )
   {
     if( PFJetPt->at(0) > 15000 ) passEEAnomJetFilter = 0;
   }

   int nBTaggedJets = 0;
   vector<double> scaleFactors;
   int nHeavyFlavorJets = 0;
   int nBTaggedHeavyFlavorJets = 0;
   int nMuons = 0;

   if( PFJetPt->size() >= 2 )
   {
     // jet, GenParticle, and muon 4-vectors
     TLorentzVector v_j, v_gp, v_m;

     // loop over two leading jets
     for(size_t i=0; i<2; ++i)
     {
       int partonFlavor = 0;

       // set jet 4-vector
       v_j.SetPtEtaPhiE(PFJetPt->at(i),PFJetEta->at(i),PFJetPhi->at(i),PFJetE->at(i));

       if( !iEvent.isRealData() )
       {
         if( matchingType==0 ) // parton-based matching
         {
           partonFlavor = abs(PFJetPartonFlavor->at(i));
          
           if( abs(PFJetPartonFlavor->at(i))==5 ) ++nHeavyFlavorJets;
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

       if( (btagger==0 && PFJetTCHE->at(i)>getPreCutValue1("TCHEL_WP")) ||
           (btagger==1 && PFJetTCHE->at(i)>getPreCutValue1("TCHEM_WP")) ||
           (btagger==2 && PFJetTCHP->at(i)>getPreCutValue1("TCHPT_WP")) ||
           (btagger==3 && PFJetSSVHE->at(i)>getPreCutValue1("SSVHEM_WP")) ||
           (btagger==4 && PFJetSSVHP->at(i)>getPreCutValue1("SSVHPT_WP")) ||
           (btagger==5 && PFJetJP->at(i)>getPreCutValue1("JPL_WP")) ||
           (btagger==6 && PFJetJP->at(i)>getPreCutValue1("JPM_WP")) ||
           (btagger==7 && PFJetJP->at(i)>getPreCutValue1("JPT_WP")) ||
           (btagger==8 && PFJetCSV->at(i)>getPreCutValue1("CSVL_WP")) ||
           (btagger==9 && PFJetCSV->at(i)>getPreCutValue1("CSVM_WP")) ||
           (btagger==10 && PFJetCSV->at(i)>getPreCutValue1("CSVT_WP"))
         )
       {
         ++nBTaggedJets;
         if( partonFlavor==5 ) ++nBTaggedHeavyFlavorJets;
         // if MC, get b-tag scale factor
         if( !iEvent.isRealData() )
         {
           if( useFixedSFs ) scaleFactors.push_back(sfCalculator.scaleFactor(partonFlavor,btagger));
           else scaleFactors.push_back(sfCalculator.scaleFactor(partonFlavor,PFJetPt->at(i),PFJetEta->at(i),btagger));
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

   fillVariableWithValue( "nSt3_q_fromRSG", nSt3_q_fromRSG, pretagWeight );
   fillVariableWithValue( "nSt3_c_fromRSG", nSt3_c_fromRSG, pretagWeight );
   fillVariableWithValue( "nSt3_b_fromRSG", nSt3_b_fromRSG, pretagWeight );
   
   fillVariableWithValue( "PassHBHENoiseFilter", ( *passHBHENoiseFilter ? 1 : 0 ), pretagWeight );
   fillVariableWithValue( "PassBeamHaloFltTight", ( *passBeamHaloFilterTight ? 1 : 0 ), pretagWeight );
   fillVariableWithValue( "PassTrackingFailure", ( *passTrackingFailure ? 1 : 0 ), pretagWeight );
   fillVariableWithValue( "PassEcalMskCellDRFlt", ( *passEcalMaskedCellDRFilter ? 1 : 0 ), pretagWeight );
   fillVariableWithValue( "PassCaloBndDRFlt", ( *passCaloBoundaryDRFilter ? 1 : 0 ), pretagWeight );
   fillVariableWithValue( "PassEEAnomJetFilter", passEEAnomJetFilter, pretagWeight );

   fillVariableWithValue( "MET_pretag", MET->front(), eventWeight );
   fillVariableWithValue( "SumET_pretag", SumET->front(), eventWeight );
   fillVariableWithValue( "METoSumET_pretag", MET->front()/SumET->front(), pretagWeight );

   fillVariableWithValue( "nJets", PFJetPt->size(), pretagWeight);

   fillVariableWithValue( "nGoodVertices_pretag", v_idx_goodPV.size(), pretagWeight );

   if( PFJetPt->size() >= 2 )
   {
     TLorentzVector dijet, jet1, jet2;
     jet1.SetPtEtaPhiE(PFJetPt->at(0),PFJetEta->at(0),PFJetPhi->at(0),PFJetE->at(0));
     jet2.SetPtEtaPhiE(PFJetPt->at(1),PFJetEta->at(1),PFJetPhi->at(1),PFJetE->at(1));

     // wide jets
     if( useWideJets )
     {
       TLorentzVector jet1_ = jet1, jet2_ = jet2, subjet;

       for(unsigned j=2; j<PFJetPt->size(); ++j)
       {
         if( fabs( PFJetEta->at(j) ) > getPreCutValue1("subleadingEtaCut") || !PFJetPassLooseID->at(j) ||
             PFJetPt->at(j) < getPreCutValue1("subleadingPtCut") ) continue;

         subjet.SetPtEtaPhiE(PFJetPt->at(j),PFJetEta->at(j),PFJetPhi->at(j),PFJetE->at(j));

         double dR1 = subjet.DeltaR(jet1_);
         double dR2 = subjet.DeltaR(jet2_);

         if (dR1 < getPreCutValue1("wideJetDeltaR") && dR1 < dR2)  jet1 += subjet;
         if (dR2 < getPreCutValue1("wideJetDeltaR") && dR2 <= dR1) jet2 += subjet;
       }
     }
         
     fillVariableWithValue( "passJetIdJ1", ( PFJetPassTightID->at(0) ? 1 : 0 ), pretagWeight );
     fillVariableWithValue( "PhiJ1_pretag", jet1.Phi(), pretagWeight );
     fillVariableWithValue( "absEtaJ1", fabs( PFJetEta->at(0) ), pretagWeight ); // even with wide jets, |eta| cut is still applied to AK5 PF jets
     fillVariableWithValue( "EtaJ1_pretag", jet1.Eta(), pretagWeight );
     fillVariableWithValue( "PtJ1_cut", jet1.Pt(), pretagWeight );
     fillVariableWithValue( "PtJ1_pretag", getVariableValue("PtJ1_cut"), pretagWeight );

     fillVariableWithValue( "passJetIdJ2", ( PFJetPassTightID->at(1) ? 1 : 0 ), pretagWeight );
     fillVariableWithValue( "PhiJ2_pretag", jet2.Phi(), pretagWeight );
     fillVariableWithValue( "absEtaJ2", fabs( PFJetEta->at(1) ), pretagWeight ); // even with wide jets, |eta| cut is still applied to AK5 PF jets
     fillVariableWithValue( "EtaJ2_pretag", jet2.Eta(), pretagWeight );
     fillVariableWithValue( "PtJ2_cut", jet2.Pt(), pretagWeight );
     fillVariableWithValue( "PtJ2_pretag", getVariableValue("PtJ2_cut"), pretagWeight );
    
     // calculate |DeltaEta(j1,j2)|
     fillVariableWithValue( "absDeltaEtaJ1J2", fabs( jet1.Eta() - jet2.Eta() ), pretagWeight );
     fillVariableWithValue( "DeltaEtaJ1J2_pretag", getVariableValue("absDeltaEtaJ1J2"), pretagWeight );
     
     // calculate M_jj
     dijet = jet1 + jet2;
     
     fillVariableWithValue( "DijetMassThreshold", dijet.M(), pretagWeight );
     
     fillVariableWithValue( "absDeltaPhiJ1J2", fabs( jet1.DeltaPhi(jet2) ), pretagWeight );
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
     fillVariableWithValue( "MET", getVariableValue("MET_pretag"), tagWeight );
     fillVariableWithValue( "SumET", getVariableValue("SumET_pretag"), tagWeight );
     fillVariableWithValue( "METoSumET", getVariableValue("METoSumET_pretag"), tagWeight );
   }
   
   // Evaluate cuts (but do not apply them)
   evaluateCuts();

   
   if(passedAllPreviousCuts("DijetMass_pretag"))
   {
     FillUserTH1D("h1_J1J2PartonFlavor", abs( PFJetPartonFlavor->at(0) ), pretagWeight );
     FillUserTH1D("h1_J1J2PartonFlavor", abs( PFJetPartonFlavor->at(1) ), pretagWeight );
     
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
   cutNames.push_back("DijetMass"); cutNames.push_back("nMuons"); cutNames.push_back("MET"); cutNames.push_back("SumET"); cutNames.push_back("METoSumET");

   if( PFJetPt->size() >= 2 && !iEvent.isRealData() && doSFReweighting )
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
       fillVariableWithValue( "MET", getVariableValue("MET_pretag"), tagWeight );
       fillVariableWithValue( "SumET", getVariableValue("SumET_pretag"), tagWeight );
       fillVariableWithValue( "METoSumET", getVariableValue("METoSumET_pretag"), tagWeight );

       // Evaluate cuts (but do not apply them)
       evaluateCuts();

       if(passedAllPreviousCuts("nMuons")) FillUserTH1D("h1_nMuons_vs_DijetMass", getVariableValue("DijetMass_pretag"), double(nMuons)*tagWeight );

       // select only those events that pass the full selection
       if( passedCut("all") ) ret = true;
     }
   }

   if( passedAllPreviousCuts("DijetMass_pretag") )
   {
     if( doSFReweighting && !iEvent.isRealData() )
     {
       FillUserTH1D("h1_DijetMass_0tag", getVariableValue("DijetMass_pretag"), eventWeight*bTagEventWeight(scaleFactors,0) );
       FillUserTH1D("h1_DijetMass_1tag", getVariableValue("DijetMass_pretag"), eventWeight*bTagEventWeight(scaleFactors,1) );
       FillUserTH1D("h1_DijetMass_2tag", getVariableValue("DijetMass_pretag"), eventWeight*bTagEventWeight(scaleFactors,2) );
     }
     else
     {
       if( nBTaggedJets==0 ) FillUserTH1D("h1_DijetMass_0tag", getVariableValue("DijetMass_pretag"), pretagWeight );
       if( nBTaggedJets==1 ) FillUserTH1D("h1_DijetMass_1tag", getVariableValue("DijetMass_pretag"), pretagWeight );
       if( nBTaggedJets==2 ) FillUserTH1D("h1_DijetMass_2tag", getVariableValue("DijetMass_pretag"), pretagWeight );
     }
   }
   
   // fill EventBin histograms
   if( doEventBins && passedAllPreviousCuts("DijetMass_pretag") )
   {
     if( nMuons==0 && max(fabs(getVariableValue("EtaJ1_pretag")),fabs(getVariableValue("EtaJ2_pretag")))<1.2 )
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
     else if( nMuons==0 && max(fabs(getVariableValue("EtaJ1_pretag")),fabs(getVariableValue("EtaJ2_pretag")))>=1.2 )
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
     else if( nMuons>=1 && max(fabs(getVariableValue("EtaJ1_pretag")),fabs(getVariableValue("EtaJ2_pretag")))<1.2 )
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
     else if( nMuons>=1 && max(fabs(getVariableValue("EtaJ1_pretag")),fabs(getVariableValue("EtaJ2_pretag")))>=1.2 )
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
       cout << category << ": PassJetID PFJet "<< i << " Pt, PtRaw, eta, phi: " << PFJetPt->at(i) << ", "
                                                                                  << PFJetPtRaw->at(i) << ", "
                                                                                  << PFJetEta->at(i) << ", "
                                                                                  << PFJetPhi->at(i) << endl;
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
       cout << category << ": PassJetID PFJet "<< i << " Pt, PtRaw, eta, phi: " << PFJetPt->at(i) << ", "
                                                                                  << PFJetPtRaw->at(i) << ", "
                                                                                  << PFJetEta->at(i) << ", "
                                                                                  << PFJetPhi->at(i) << endl;
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
  SFb_shift_ = 0.;
  SFl_shift_ = 0.;
  TCHEL_SFb_ = 1.;
  TCHEL_SFl_ = 1.;
  TCHPT_SFb_ = 1.;
  TCHPT_SFl_ = 1.;
  SSVHPT_SFb_ = 1.;
  SSVHPT_SFl_ = 1.;

  double PtBins_b[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670};
  // TCHEL
  TCHEL_SFb_0to2p4 = new TF1("TCHEL_SFb_0to2p4","0.603913*((1.+(0.286361*x))/(1.+(0.170474*x)))", 30.,670.);

  TCHEL_SFb_errors = new TH1D("TCHEL_SFb_errors", "TCHEL_SFb_errors", 14, PtBins_b);
  TCHEL_SFb_errors->SetBinContent( 0,0.12);
  TCHEL_SFb_errors->SetBinContent( 1,0.0244956);
  TCHEL_SFb_errors->SetBinContent( 2,0.0237293);
  TCHEL_SFb_errors->SetBinContent( 3,0.0180131);
  TCHEL_SFb_errors->SetBinContent( 4,0.0182411);
  TCHEL_SFb_errors->SetBinContent( 5,0.0184592);
  TCHEL_SFb_errors->SetBinContent( 6,0.0106444);
  TCHEL_SFb_errors->SetBinContent( 7,0.011073);
  TCHEL_SFb_errors->SetBinContent( 8,0.0106296);
  TCHEL_SFb_errors->SetBinContent( 9,0.0175259);
  TCHEL_SFb_errors->SetBinContent(10,0.0161566);
  TCHEL_SFb_errors->SetBinContent(11,0.0158973);
  TCHEL_SFb_errors->SetBinContent(12,0.0186782);
  TCHEL_SFb_errors->SetBinContent(13,0.0371113);
  TCHEL_SFb_errors->SetBinContent(14,0.0289788);
  TCHEL_SFb_errors->SetBinContent(15,(2*0.0289788));

  TCHEL_SFl_0to2p4 =   new TF1("TCHEL_SFl_0to2p4","(1.10649*((1+(-9.00297e-05*x))+(2.32185e-07*(x*x))))+(-4.04925e-10*(x*(x*(x/(1+(-0.00051036*x))))))", 20.,670.);
  TCHEL_SFl_0to0p5 =   new TF1("TCHEL_SFl_0to0p5","(1.13615*((1+(-0.00119852*x))+(1.17888e-05*(x*x))))+(-9.8581e-08*(x*(x*(x/(1+(0.00689317*x))))))", 20.,670.);
  TCHEL_SFl_0p5to1p0 = new TF1("TCHEL_SFl_0p5to1p0","(1.13277*((1+(-0.00084146*x))+(3.80313e-06*(x*x))))+(-8.75061e-09*(x*(x*(x/(1+(0.00118695*x))))))", 20.,670.);
  TCHEL_SFl_1p0to1p5 = new TF1("TCHEL_SFl_1p0to1p5","(1.17163*((1+(-0.000828475*x))+(3.0769e-06*(x*x))))+(-4.692e-09*(x*(x*(x/(1+(0.000337759*x))))))", 20.,670.);
  TCHEL_SFl_1p5to2p4 = new TF1("TCHEL_SFl_1p5to2p4","(1.14554*((1+(-0.000128043*x))+(4.10899e-07*(x*x))))+(-2.07565e-10*(x*(x*(x/(1+(-0.00118618*x))))))", 20.,670.);

  TCHEL_SFl_0to2p4_min =   new TF1("TCHEL_SFl_0to2p4_min","(1.01541*((1+(-6.04627e-05*x))+(1.38195e-07*(x*x))))+(-2.83043e-10*(x*(x*(x/(1+(-0.000633609*x))))))", 20.,670.);
  TCHEL_SFl_0to0p5_min =   new TF1("TCHEL_SFl_0to0p5_min","(1.0369*((1+(-0.000945578*x))+(7.73273e-06*(x*x))))+(-4.47791e-08*(x*(x*(x/(1+(0.00499343*x))))))", 20.,670.);
  TCHEL_SFl_0p5to1p0_min = new TF1("TCHEL_SFl_0p5to1p0_min","(0.983748*((1+(7.13613e-05*x))+(-1.08648e-05*(x*x))))+(2.96162e-06*(x*(x*(x/(1+(0.282104*x))))))", 20.,670.);
  TCHEL_SFl_1p0to1p5_min = new TF1("TCHEL_SFl_1p0to1p5_min","(1.0698*((1+(-0.000731877*x))+(2.56922e-06*(x*x))))+(-3.0318e-09*(x*(x*(x/(1+(5.04118e-05*x))))))", 20.,670.);
  TCHEL_SFl_1p5to2p4_min = new TF1("TCHEL_SFl_1p5to2p4_min","(1.04766*((1+(-6.87499e-05*x))+(2.2454e-07*(x*x))))+(-1.18395e-10*(x*(x*(x/(1+(-0.00128734*x))))))", 20.,670.);

  TCHEL_SFl_0to2p4_max =   new TF1("TCHEL_SFl_0to2p4_max","(1.19751*((1+(-0.000114197*x))+(3.08558e-07*(x*x))))+(-5.27598e-10*(x*(x*(x/(1+(-0.000422372*x))))))", 20.,670.);
  TCHEL_SFl_0to0p5_max =   new TF1("TCHEL_SFl_0to0p5_max","(1.22179*((1+(-0.000946228*x))+(7.37821e-06*(x*x))))+(-4.8451e-08*(x*(x*(x/(1+(0.0047976*x))))))", 20.,670.);
  TCHEL_SFl_0p5to1p0_max = new TF1("TCHEL_SFl_0p5to1p0_max","(1.22714*((1+(-0.00085562*x))+(3.74425e-06*(x*x))))+(-8.91028e-09*(x*(x*(x/(1+(0.00109346*x))))))", 20.,670.);
  TCHEL_SFl_1p0to1p5_max = new TF1("TCHEL_SFl_1p0to1p5_max","(1.27351*((1+(-0.000911891*x))+(3.5465e-06*(x*x))))+(-6.69625e-09*(x*(x*(x/(1+(0.000590847*x))))))", 20.,670.);
  TCHEL_SFl_1p5to2p4_max = new TF1("TCHEL_SFl_1p5to2p4_max","(1.24367*((1+(-0.000182494*x))+(5.92637e-07*(x*x))))+(-3.3745e-10*(x*(x*(x/(1+(-0.00107694*x))))))", 20.,670.);

  // CSVL
  CSVL_SFb_0to2p4 = new TF1("CSVL_SFb_0to2p4","1.02658*((1.+(0.0195388*x))/(1.+(0.0209145*x)))", 30.,670.);

  CSVL_SFb_errors = new TH1D("CSVL_SFb_errors", "CSVL_SFb_errors", 14, PtBins_b);
  CSVL_SFb_errors->SetBinContent( 0,0.12);
  CSVL_SFb_errors->SetBinContent( 1,0.0188743);
  CSVL_SFb_errors->SetBinContent( 2,0.0161816);
  CSVL_SFb_errors->SetBinContent( 3,0.0139824);
  CSVL_SFb_errors->SetBinContent( 4,0.0152644);
  CSVL_SFb_errors->SetBinContent( 5,0.0161226);
  CSVL_SFb_errors->SetBinContent( 6,0.0157396);
  CSVL_SFb_errors->SetBinContent( 7,0.0161619);
  CSVL_SFb_errors->SetBinContent( 8,0.0168747);
  CSVL_SFb_errors->SetBinContent( 9,0.0257175);
  CSVL_SFb_errors->SetBinContent(10,0.026424);
  CSVL_SFb_errors->SetBinContent(11,0.0264928);
  CSVL_SFb_errors->SetBinContent(12,0.0315127);
  CSVL_SFb_errors->SetBinContent(13,0.030734);
  CSVL_SFb_errors->SetBinContent(14,0.0438259);
  CSVL_SFb_errors->SetBinContent(15,(2*0.0438259));

  CSVL_SFl_0to2p4 =   new TF1("CSVL_SFl_0to2p4","((1.0344+(0.000962994*x))+(-3.65392e-06*(x*x)))+(3.23525e-09*(x*(x*x)))", 20.,670.);
  CSVL_SFl_0to0p5 =   new TF1("CSVL_SFl_0to0p5","((1.07536+(0.000175506*x))+(-8.63317e-07*(x*x)))+(3.27516e-10*(x*(x*x)))", 20.,670.);
  CSVL_SFl_0p5to1p0 = new TF1("CSVL_SFl_0p5to1p0","((1.07846+(0.00032458*x))+(-1.30258e-06*(x*x)))+(8.50608e-10*(x*(x*x)))", 20.,670.);
  CSVL_SFl_1p0to1p5 = new TF1("CSVL_SFl_1p0to1p5","((1.08294+(0.000474818*x))+(-1.43857e-06*(x*x)))+(1.13308e-09*(x*(x*x)))", 20.,670.);
  CSVL_SFl_1p5to2p4 = new TF1("CSVL_SFl_1p5to2p4","((1.0617+(0.000173654*x))+(-5.29009e-07*(x*x)))+(5.55931e-10*(x*(x*x)))", 20.,670.);

  CSVL_SFl_0to2p4_min =   new TF1("CSVL_SFl_0to2p4_min","((0.956023+(0.000825106*x))+(-3.18828e-06*(x*x)))+(2.81787e-09*(x*(x*x)))", 20.,670.);
  CSVL_SFl_0to0p5_min =   new TF1("CSVL_SFl_0to0p5_min","((0.994425+(-8.66392e-05*x))+(-3.03813e-08*(x*x)))+(-3.52151e-10*(x*(x*x)))", 20.,670.);
  CSVL_SFl_0p5to1p0_min = new TF1("CSVL_SFl_0p5to1p0_min","((0.998088+(6.94916e-05*x))+(-4.82731e-07*(x*x)))+(1.63506e-10*(x*(x*x)))", 20.,670.);
  CSVL_SFl_1p0to1p5_min = new TF1("CSVL_SFl_1p0to1p5_min","((1.00294+(0.000289844*x))+(-7.9845e-07*(x*x)))+(5.38525e-10*(x*(x*x)))", 20.,670.);
  CSVL_SFl_1p5to2p4_min = new TF1("CSVL_SFl_1p5to2p4_min","((0.979816+(0.000138797*x))+(-3.14503e-07*(x*x)))+(2.38124e-10*(x*(x*x)))", 20.,670.);

  CSVL_SFl_0to2p4_max =   new TF1("CSVL_SFl_0to2p4_max","((1.11272+(0.00110104*x))+(-4.11956e-06*(x*x)))+(3.65263e-09*(x*(x*x)))", 20.,670.);
  CSVL_SFl_0to0p5_max =   new TF1("CSVL_SFl_0to0p5_max","((1.15628+(0.000437668*x))+(-1.69625e-06*(x*x)))+(1.00718e-09*(x*(x*x)))", 20.,670.);
  CSVL_SFl_0p5to1p0_max = new TF1("CSVL_SFl_0p5to1p0_max","((1.15882+(0.000579711*x))+(-2.12243e-06*(x*x)))+(1.53771e-09*(x*(x*x)))", 20.,670.);
  CSVL_SFl_1p0to1p5_max = new TF1("CSVL_SFl_1p0to1p5_max","((1.16292+(0.000659848*x))+(-2.07868e-06*(x*x)))+(1.72763e-09*(x*(x*x)))", 20.,670.);
  CSVL_SFl_1p5to2p4_max = new TF1("CSVL_SFl_1p5to2p4_max","((1.14357+(0.00020854*x))+(-7.43519e-07*(x*x)))+(8.73742e-10*(x*(x*x)))", 20.,670.);

  // CSVM
  CSVM_SFb_0to2p4 = new TF1("CSVM_SFb_0to2p4","0.6981*((1.+(0.414063*x))/(1.+(0.300155*x)))", 30.,670.);

  CSVM_SFb_errors = new TH1D("CSVM_SFb_errors", "CSVM_SFb_errors", 14, PtBins_b);
  CSVM_SFb_errors->SetBinContent( 0,0.12);
  CSVM_SFb_errors->SetBinContent( 1,0.0295675);
  CSVM_SFb_errors->SetBinContent( 2,0.0295095);
  CSVM_SFb_errors->SetBinContent( 3,0.0210867);
  CSVM_SFb_errors->SetBinContent( 4,0.0219349);
  CSVM_SFb_errors->SetBinContent( 5,0.0227033);
  CSVM_SFb_errors->SetBinContent( 6,0.0204062);
  CSVM_SFb_errors->SetBinContent( 7,0.0185857);
  CSVM_SFb_errors->SetBinContent( 8,0.0256242);
  CSVM_SFb_errors->SetBinContent( 9,0.0383341);
  CSVM_SFb_errors->SetBinContent(10,0.0409675);
  CSVM_SFb_errors->SetBinContent(11,0.0420284);
  CSVM_SFb_errors->SetBinContent(12,0.0541299);
  CSVM_SFb_errors->SetBinContent(13,0.0578761);
  CSVM_SFb_errors->SetBinContent(14,0.0655432);
  CSVM_SFb_errors->SetBinContent(15,(2*0.0655432));

  CSVM_SFl_0to2p4 =   new TF1("CSVM_SFl_0to2p4","((1.04318+(0.000848162*x))+(-2.5795e-06*(x*x)))+(1.64156e-09*(x*(x*x)))", 20.,670.);
  CSVM_SFl_0to0p8 =   new TF1("CSVM_SFl_0to0p8","((1.06182+(0.000617034*x))+(-1.5732e-06*(x*x)))+(3.02909e-10*(x*(x*x)))", 20.,670.);
  CSVM_SFl_0p8to1p6 = new TF1("CSVM_SFl_0p8to1p6","((1.111+(-9.64191e-06*x))+(1.80811e-07*(x*x)))+(-5.44868e-10*(x*(x*x)))", 20.,670.);
  CSVM_SFl_1p6to2p4 = new TF1("CSVM_SFl_1p6to2p4","((1.08498+(-0.000701422*x))+(3.43612e-06*(x*x)))+(-4.11794e-09*(x*(x*x)))", 20.,670.);

  CSVM_SFl_0to2p4_min =   new TF1("CSVM_SFl_0to2p4_min","((0.962627+(0.000448344*x))+(-1.25579e-06*(x*x)))+(4.82283e-10*(x*(x*x)))", 20.,670.);
  CSVM_SFl_0to0p8_min =   new TF1("CSVM_SFl_0to0p8_min","((0.972455+(7.51396e-06*x))+(4.91857e-07*(x*x)))+(-1.47661e-09*(x*(x*x)))", 20.,670.);
  CSVM_SFl_0p8to1p6_min = new TF1("CSVM_SFl_0p8to1p6_min","((1.02055+(-0.000378856*x))+(1.49029e-06*(x*x)))+(-1.74966e-09*(x*(x*x)))", 20.,670.);
  CSVM_SFl_1p6to2p4_min = new TF1("CSVM_SFl_1p6to2p4_min","((0.983476+(-0.000607242*x))+(3.17997e-06*(x*x)))+(-4.01242e-09*(x*(x*x)))", 20.,670.);

  CSVM_SFl_0to2p4_max =   new TF1("CSVM_SFl_0to2p4_max","((1.12368+(0.00124806*x))+(-3.9032e-06*(x*x)))+(2.80083e-09*(x*(x*x)))", 20.,670.);
  CSVM_SFl_0to0p8_max =   new TF1("CSVM_SFl_0to0p8_max","((1.15116+(0.00122657*x))+(-3.63826e-06*(x*x)))+(2.08242e-09*(x*(x*x)))", 20.,670.);
  CSVM_SFl_0p8to1p6_max = new TF1("CSVM_SFl_0p8to1p6_max","((1.20146+(0.000359543*x))+(-1.12866e-06*(x*x)))+(6.59918e-10*(x*(x*x)))", 20.,670.);
  CSVM_SFl_1p6to2p4_max = new TF1("CSVM_SFl_1p6to2p4_max","((1.18654+(-0.000795808*x))+(3.69226e-06*(x*x)))+(-4.22347e-09*(x*(x*x)))", 20.,670.);
}

// ------------ method that initializes the BTagScaleFactorCalculator class  ------------
void
BTagScaleFactorCalculator::init(const double SFb_shift, const double SFl_shift, const double TCHEL_SFb, const double TCHEL_SFl, const double TCHPT_SFb, const double TCHPT_SFl, const double SSVHPT_SFb, const double SSVHPT_SFl)
{
  SFb_shift_ = SFb_shift;
  SFl_shift_ = SFl_shift;
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
    else if(btagger==2) return TCHPT_SFb_;
    else if(btagger==4) return SSVHPT_SFb_;
    else                return 1.;
  }
  else
  {
    if(btagger==0)      return TCHEL_SFl_;
    else if(btagger==2) return TCHPT_SFl_;
    else if(btagger==4) return SSVHPT_SFl_;
    else                return 1.;
  }
}


// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor  ------------
double
BTagScaleFactorCalculator::scaleFactor(const int partonFlavor, const double jetPt, const double jetEta, const int btagger)
{
  if( partonFlavor==5 )
  {
    if(btagger==0)      return scaleFactorB_TCHEL(jetPt,jetEta);
    else if(btagger==8) return scaleFactorB_CSVL(jetPt,jetEta);
    else if(btagger==9) return scaleFactorB_CSVM(jetPt,jetEta);
    else                return 1.;
  }
  else if( partonFlavor==4 )
  {
    if(btagger==0)      return scaleFactorC_TCHEL(jetPt,jetEta);
    else if(btagger==8) return scaleFactorC_CSVL(jetPt,jetEta);
    else if(btagger==9) return scaleFactorC_CSVM(jetPt,jetEta);
    else                return 1.;
  }
  else
  {
    if(btagger==0)      return scaleFactorUDSG_TCHEL(jetPt,jetEta);
    else if(btagger==8) return scaleFactorUDSG_CSVL(jetPt,jetEta);
    else if(btagger==9) return scaleFactorUDSG_CSVM(jetPt,jetEta);
    else                return 1.;
  }
}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for b-jets and TCHEL tagger  ------------
double
BTagScaleFactorCalculator::scaleFactorB_TCHEL(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  // for scale factor extrapolation
  if(Pt<30) Pt = 30;
  if(Pt>670) Pt = 670;

  return TCHEL_SFb_0to2p4->Eval(Pt) + SFb_shift_*TCHEL_SFb_errors->GetBinContent(TCHEL_SFb_errors->GetXaxis()->FindBin(jetPt));
}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for c-jets and TCHEL tagger  ------------
double
BTagScaleFactorCalculator::scaleFactorC_TCHEL(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  // for scale factor extrapolation
  if(Pt<30) Pt = 30;
  if(Pt>670) Pt = 670;

  return TCHEL_SFb_0to2p4->Eval(Pt) + 2*SFb_shift_*TCHEL_SFb_errors->GetBinContent(TCHEL_SFb_errors->GetXaxis()->FindBin(jetPt));
}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for light flavor jets and TCHEL tagger ------------
double
BTagScaleFactorCalculator::scaleFactorUDSG_TCHEL(const double jetPt, const double jetEta)
{
  double SF = 1.;
  double Pt = jetPt;
  double absEta = fabs(jetEta);
  // for scale factor extrapolation
  if(Pt<20) Pt = 20;


  if( Pt>670 )
    SF = TCHEL_SFl_0to2p4->Eval(670) + 2*fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (TCHEL_SFl_0to2p4_max->Eval(670) - TCHEL_SFl_0to2p4->Eval(670)) : (TCHEL_SFl_0to2p4_min->Eval(670) - TCHEL_SFl_0to2p4->Eval(670)) );
  else
  {
    if(absEta<0.5)
      SF = TCHEL_SFl_0to0p5->Eval(Pt) + fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (TCHEL_SFl_0to0p5_max->Eval(Pt) - TCHEL_SFl_0to0p5->Eval(Pt)) : (TCHEL_SFl_0to0p5_min->Eval(Pt) - TCHEL_SFl_0to0p5->Eval(Pt)) );
    else if(absEta>=0.5 && absEta<1.)
      SF = TCHEL_SFl_0p5to1p0->Eval(Pt) + fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (TCHEL_SFl_0p5to1p0_max->Eval(Pt) - TCHEL_SFl_0p5to1p0->Eval(Pt)) : (TCHEL_SFl_0p5to1p0_min->Eval(Pt) - TCHEL_SFl_0p5to1p0->Eval(Pt)) );
    else if(absEta>=1. && absEta<1.5)
      SF = TCHEL_SFl_1p0to1p5->Eval(Pt) + fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (TCHEL_SFl_1p0to1p5_max->Eval(Pt) - TCHEL_SFl_1p0to1p5->Eval(Pt)) : (TCHEL_SFl_1p0to1p5_min->Eval(Pt) - TCHEL_SFl_1p0to1p5->Eval(Pt)) );
    else
      SF = TCHEL_SFl_1p5to2p4->Eval(Pt) + fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (TCHEL_SFl_1p5to2p4_max->Eval(Pt) - TCHEL_SFl_1p5to2p4->Eval(Pt)) : (TCHEL_SFl_1p5to2p4_min->Eval(Pt) - TCHEL_SFl_1p5to2p4->Eval(Pt)) );
  }

  return SF;
}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for b-jets and CSVL tagger  ------------
double
BTagScaleFactorCalculator::scaleFactorB_CSVL(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  // for scale factor extrapolation
  if(Pt<30) Pt = 30;
  if(Pt>670) Pt = 670;

  return CSVL_SFb_0to2p4->Eval(Pt) + SFb_shift_*CSVL_SFb_errors->GetBinContent(CSVL_SFb_errors->GetXaxis()->FindBin(jetPt));
}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for c-jets and CSVL tagger  ------------
double
BTagScaleFactorCalculator::scaleFactorC_CSVL(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  // for scale factor extrapolation
  if(Pt<30) Pt = 30;
  if(Pt>670) Pt = 670;

  return CSVL_SFb_0to2p4->Eval(Pt) + 2*SFb_shift_*CSVL_SFb_errors->GetBinContent(CSVL_SFb_errors->GetXaxis()->FindBin(jetPt));
}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for light flavor jets and CSVL tagger ------------
double
BTagScaleFactorCalculator::scaleFactorUDSG_CSVL(const double jetPt, const double jetEta)
{
  double SF = 1.;
  double Pt = jetPt;
  double absEta = fabs(jetEta);
  // for scale factor extrapolation
  if(Pt<20) Pt = 20;


  if( Pt>670 )
    SF = CSVL_SFl_0to2p4->Eval(670) + 2*fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (CSVL_SFl_0to2p4_max->Eval(670) - CSVL_SFl_0to2p4->Eval(670)) : (CSVL_SFl_0to2p4_min->Eval(670) - CSVL_SFl_0to2p4->Eval(670)) );
  else
  {
    if(absEta<0.5)
      SF = CSVL_SFl_0to0p5->Eval(Pt) + fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (CSVL_SFl_0to0p5_max->Eval(Pt) - CSVL_SFl_0to0p5->Eval(Pt)) : (CSVL_SFl_0to0p5_min->Eval(Pt) - CSVL_SFl_0to0p5->Eval(Pt)) );
    else if(absEta>=0.5 && absEta<1.)
      SF = CSVL_SFl_0p5to1p0->Eval(Pt) + fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (CSVL_SFl_0p5to1p0_max->Eval(Pt) - CSVL_SFl_0p5to1p0->Eval(Pt)) : (CSVL_SFl_0p5to1p0_min->Eval(Pt) - CSVL_SFl_0p5to1p0->Eval(Pt)) );
    else if(absEta>=1. && absEta<1.5)
      SF = CSVL_SFl_1p0to1p5->Eval(Pt) + fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (CSVL_SFl_1p0to1p5_max->Eval(Pt) - CSVL_SFl_1p0to1p5->Eval(Pt)) : (CSVL_SFl_1p0to1p5_min->Eval(Pt) - CSVL_SFl_1p0to1p5->Eval(Pt)) );
    else
      SF = CSVL_SFl_1p5to2p4->Eval(Pt) + fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (CSVL_SFl_1p5to2p4_max->Eval(Pt) - CSVL_SFl_1p5to2p4->Eval(Pt)) : (CSVL_SFl_1p5to2p4_min->Eval(Pt) - CSVL_SFl_1p5to2p4->Eval(Pt)) );
  }

  return SF;
}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for b-jets and CSVM tagger  ------------
double
BTagScaleFactorCalculator::scaleFactorB_CSVM(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  // for scale factor extrapolation
  if(Pt<30) Pt = 30;
  if(Pt>670) Pt = 670;

  return CSVM_SFb_0to2p4->Eval(Pt) + SFb_shift_*CSVM_SFb_errors->GetBinContent(CSVM_SFb_errors->GetXaxis()->FindBin(jetPt));
}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for c-jets and CSVM tagger  ------------
double
BTagScaleFactorCalculator::scaleFactorC_CSVM(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  // for scale factor extrapolation
  if(Pt<30) Pt = 30;
  if(Pt>670) Pt = 670;

  return CSVM_SFb_0to2p4->Eval(Pt) + 2*SFb_shift_*CSVM_SFb_errors->GetBinContent(CSVM_SFb_errors->GetXaxis()->FindBin(jetPt));
}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for light flavor jets and CSVM tagger ------------
double
BTagScaleFactorCalculator::scaleFactorUDSG_CSVM(const double jetPt, const double jetEta)
{
  double SF = 1.;
  double Pt = jetPt;
  double absEta = fabs(jetEta);
  // for scale factor extrapolation
  if(Pt<20) Pt = 20;


  if( Pt>670 )
    SF = CSVM_SFl_0to2p4->Eval(670) + 2*fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (CSVM_SFl_0to2p4_max->Eval(670) - CSVM_SFl_0to2p4->Eval(670)) : (CSVM_SFl_0to2p4_min->Eval(670) - CSVM_SFl_0to2p4->Eval(670)) );
  else
  {
    if(absEta<0.8)
      SF = CSVM_SFl_0to0p8->Eval(Pt) + fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (CSVM_SFl_0to0p8_max->Eval(Pt) - CSVM_SFl_0to0p8->Eval(Pt)) : (CSVM_SFl_0to0p8_min->Eval(Pt) - CSVM_SFl_0to0p8->Eval(Pt)) );
    else if(absEta>=0.8 && absEta<1.6)
      SF = CSVM_SFl_0p8to1p6->Eval(Pt) + fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (CSVM_SFl_0p8to1p6_max->Eval(Pt) - CSVM_SFl_0p8to1p6->Eval(Pt)) : (CSVM_SFl_0p8to1p6_min->Eval(Pt) - CSVM_SFl_0p8to1p6->Eval(Pt)) );
    else
      SF = CSVM_SFl_1p6to2p4->Eval(Pt) + fabs(SFl_shift_)*( SFl_shift_ >= 0. ? (CSVM_SFl_1p6to2p4_max->Eval(Pt) - CSVM_SFl_1p6to2p4->Eval(Pt)) : (CSVM_SFl_1p6to2p4_min->Eval(Pt) - CSVM_SFl_1p6to2p4->Eval(Pt)) );
  }

  return SF;
}


DEFINE_FWK_MODULE(MyAnalyzer);
