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
// $Id: MyAnalyzer_bTaggingEfficiency_DijetBBTag_2011.cc,v 1.7 2012/05/10 05:34:41 ferencek Exp $
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

// BaseClass
#include "MyAnalysis/MyAnalyzer/interface/BaseClass.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT
#include <TLorentzVector.h>


using namespace std;

//
// class declaration
//

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

      // ----------member data ---------------------------

};

//
// constructors and destructor
//

MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig) :
  BaseClass(iConfig)
{
   //now do whatever initialization is needed

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
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_All;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_TCHEL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_TCHEM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_TCHPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_TCHPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_SSVHEM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_SSVHPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_JPL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_JPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_JPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_CSVL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_CSVM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_B_CSVT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);

   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_All;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_TCHEL;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_TCHEM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_TCHPM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_TCHPT;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_SSVHEM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_SSVHPT;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_JPL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_JPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_JPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_CSVL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_CSVM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_B_CSVT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);

   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_All;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_TCHEL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_TCHEM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_TCHPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_TCHPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_SSVHEM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_SSVHPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_JPL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_JPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_JPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_CSVL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_CSVM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_C_CSVT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);

   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_All;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_TCHEL;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_TCHEM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_TCHPM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_TCHPT;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_SSVHEM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_SSVHPT;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_JPL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_JPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_JPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_CSVL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_CSVM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_C_CSVT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);

   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_All;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_TCHEL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_TCHEM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_TCHPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_TCHPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_SSVHEM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_SSVHPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_JPL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_JPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_JPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_CSVL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_CSVM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_CSVT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);

   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_All;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_TCHEL;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_TCHEM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_TCHPM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_TCHPT;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_SSVHEM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_SSVHPT;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_JPL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_JPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_JPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_CSVL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_CSVM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_CSVT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);

   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_All;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_TCHEL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_TCHEM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_TCHPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_TCHPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_SSVHEM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_SSVHPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_JPL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_JPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_JPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_CSVL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_CSVM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_G_CSVT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);

   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_All;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_TCHEL;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_TCHEM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_TCHPM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_TCHPT;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_SSVHEM;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_SSVHPT;p_{T,2} [GeV];#eta_{2}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_JPL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_JPM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_JPT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_CSVL;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_CSVM;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_G_CSVT;p_{T,1} [GeV];#eta_{1}", 600, 0, 6000, 100, -5, 5);

   // initialize your variables here

   //############################# User's code ends here #################################
   //#####################################################################################
}

// ------------ method called when starting to processes a run  ------------
bool
MyAnalyzer::beginRun(edm::Run& iRun, const edm::EventSetup& iSetup)
{
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

   int matchingType = int(getPreCutValue1("matchingType"));
   double matchingRadius = getPreCutValue1("matchingRadius");
   int useAK5Jets = int(getPreCutValue1("useAK5Jets"));
   
   // grab necessary objects from the event
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
   
   edm::Handle<vector<double> > PFJetPt_;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:Pt"), PFJetPt_);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:Pt"), PFJetPt_);
   edm::Handle<vector<double> > PFJetPtRaw;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:PtRaw"), PFJetPtRaw);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:PtRaw"), PFJetPtRaw);
   edm::Handle<vector<double> > PFJetEta;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:Eta"), PFJetEta);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:Eta"), PFJetEta);
   edm::Handle<vector<double> > PFJetPhi;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:Phi"), PFJetPhi);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:Phi"), PFJetPhi);
   edm::Handle<vector<double> > PFJetE_;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:Energy"), PFJetE_);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:Energy"), PFJetE_);
   edm::Handle<vector<double> > PFJetUnc;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:JECUnc"), PFJetUnc);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:JECUnc"), PFJetUnc);
   edm::Handle<vector<int> > PFJetPassLooseID;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:PassLooseID"), PFJetPassLooseID);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:PassLooseID"), PFJetPassLooseID);
   edm::Handle<vector<int> > PFJetPassTightID;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:PassTightID"), PFJetPassTightID);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:PassTightID"), PFJetPassTightID);
   edm::Handle<vector<double> > PFJetSSVHE;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:SimpleSecondaryVertexHighEffBTag"), PFJetSSVHE);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:SimpleSecondaryVertexHighEffBTag"), PFJetSSVHE);
   edm::Handle<vector<double> > PFJetSSVHP;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:SimpleSecondaryVertexHighPurBTag"), PFJetSSVHP);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:SimpleSecondaryVertexHighPurBTag"), PFJetSSVHP);
   edm::Handle<vector<double> > PFJetTCHE;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:TrackCountingHighEffBTag"), PFJetTCHE);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:TrackCountingHighEffBTag"), PFJetTCHE);
   edm::Handle<vector<double> > PFJetTCHP;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:TrackCountingHighPurBTag"), PFJetTCHP);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:TrackCountingHighPurBTag"), PFJetTCHP);
   edm::Handle<vector<double> > PFJetJP;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:JetProbabilityBTag"), PFJetJP);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:JetProbabilityBTag"), PFJetJP);
   edm::Handle<vector<double> > PFJetCSV;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:CombinedSecondaryVertexBJetTag"), PFJetCSV);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:CombinedSecondaryVertexBJetTag"), PFJetCSV);
   edm::Handle<vector<int> > PFJetPartonFlavor;
   if(useAK5Jets) iEvent.getByLabel(edm::InputTag("AK5PFJets:PartonFlavor"), PFJetPartonFlavor);
   else            iEvent.getByLabel(edm::InputTag("AK7PFJets:PartonFlavor"), PFJetPartonFlavor);

   edm::Handle<unsigned int> ProcessID;
   iEvent.getByLabel(edm::InputTag("GenEventInfo:ProcessID"), ProcessID);

   
   auto_ptr<std::vector<double> >  PFJetPt ( new std::vector<double>() );
   auto_ptr<std::vector<double> >  PFJetE  ( new std::vector<double>() );

   for(size_t i=0; i<PFJetPt_->size(); i++)
   {
       double JES_ScaleFactor = 1.;
       if( !iEvent.isRealData() ) JES_ScaleFactor = 1. + getPreCutValue1("JES_Shift")*PFJetUnc->at(i);

       PFJetPt->push_back( PFJetPt_->at(i)*JES_ScaleFactor );
       PFJetE ->push_back( PFJetE_ ->at(i)*JES_ScaleFactor );
   }

   // Set the evaluation of the cuts to false and clear the variable values and filled status
   resetCuts();

   int gg_FinalState = 0;

   if( *ProcessID==13 || *ProcessID==68 ) gg_FinalState = 1;

   fillVariableWithValue( "gg_FinalState", gg_FinalState );
   
   int nSt3_b = 0, nSt2_b = 0, nSt3_q_fromRSG = 0, nSt3_c_fromRSG = 0, nSt3_b_fromRSG = 0;

   for(size_t i=0; i<GenParticlePt->size(); i++)
   {
     //cout << "Index= " << i << " PdgId=" << GenParticlePdgId->at(i) << " Status=" << GenParticleStatus->at(i) << " Mother index=" << GenParticleMotherIndex->at(i) << endl;
     if( abs(GenParticlePdgId->at(i))==5 && GenParticleStatus->at(i)==3 ) ++nSt3_b;
     if( abs(GenParticlePdgId->at(i))==5 && GenParticleStatus->at(i)==2 ) ++nSt2_b;
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

   fillVariableWithValue( "nSt3_b", nSt3_b );
   fillVariableWithValue( "nSt2_b", nSt2_b );
   fillVariableWithValue( "nSt3_q_fromRSG", nSt3_q_fromRSG );
   fillVariableWithValue( "nSt3_c_fromRSG", nSt3_c_fromRSG );
   fillVariableWithValue( "nSt3_b_fromRSG", nSt3_b_fromRSG );

   fillVariableWithValue( "nJets", PFJetPt->size() );

   if( PFJetPt->size() >= 1 )
   {
     fillVariableWithValue( "passJetIdJ1", ( PFJetPassTightID->at(0) ? 1 : 0 ) );
     fillVariableWithValue( "absEtaJ1", fabs( PFJetEta->at(0) ));
   }
   if( PFJetPt->size() >= 2 )
   {
     fillVariableWithValue( "passJetIdJ2", ( PFJetPassTightID->at(1) ? 1 : 0 ) );
     fillVariableWithValue( "absEtaJ2", fabs( PFJetEta->at(1) ) );
     
     TLorentzVector v_j1j2, v_j1, v_j2;
     v_j1.SetPtEtaPhiE(PFJetPt->at(0),PFJetEta->at(0),PFJetPhi->at(0),PFJetE->at(0));
     v_j2.SetPtEtaPhiE(PFJetPt->at(1),PFJetEta->at(1),PFJetPhi->at(1),PFJetE->at(1));
     // calculate |DeltaEta(j1,j2)|
     fillVariableWithValue( "absDeltaEtaJ1J2", fabs( PFJetEta->at(0) - PFJetEta->at(1) ) );
     // calculate M_j1j2
     v_j1j2 = v_j1 + v_j2;
     fillVariableWithValue( "DijetMass", v_j1j2.M() );
   }
   
   // Evaluate cuts (but do not apply them)
   evaluateCuts();
   

   int nHeavyFlavorJets = 0;
   
   if( PFJetPt->size() >= 2 )
   {
     // jet and GenParticle 4-vectors
     TLorentzVector v_j, v_gp;

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

       if( passedAllPreviousCuts("nJets") )
       {
         if( partonFlavor==5 )
         {
           if(i==0)
           {
             FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_All", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEL_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_TCHEL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEM_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_TCHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPM_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_TCHPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPT_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_TCHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHE->at(i)>getPreCutValue1("SSVHEM_WP") ) FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_SSVHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHP->at(i)>getPreCutValue1("SSVHPT_WP") ) FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_SSVHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPL_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_JPL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPM_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_JPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPT_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_JPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVL_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_CSVL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVM_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_CSVM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVT_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_B_CSVT", PFJetPt->at(i), PFJetEta->at(i));
           }
           if(i==1)
           {
             FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_All", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEL_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_TCHEL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEM_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_TCHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPM_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_TCHPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPT_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_TCHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHE->at(i)>getPreCutValue1("SSVHEM_WP") ) FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_SSVHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHP->at(i)>getPreCutValue1("SSVHPT_WP") ) FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_SSVHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPL_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_JPL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPM_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_JPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPT_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_JPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVL_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_CSVL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVM_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_CSVM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVT_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_B_CSVT", PFJetPt->at(i), PFJetEta->at(i));
           }
         }
         else if( partonFlavor==4 )
         {
           if(i==0)
           {
             FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_All", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEL_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_TCHEL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEM_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_TCHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPM_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_TCHPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPT_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_TCHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHE->at(i)>getPreCutValue1("SSVHEM_WP") ) FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_SSVHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHP->at(i)>getPreCutValue1("SSVHPT_WP") ) FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_SSVHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPL_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_JPL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPM_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_JPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPT_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_JPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVL_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_CSVL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVM_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_CSVM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVT_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_C_CSVT", PFJetPt->at(i), PFJetEta->at(i));
           }
           if(i==1)
           {
             FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_All", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEL_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_TCHEL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEM_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_TCHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPM_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_TCHPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPT_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_TCHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHE->at(i)>getPreCutValue1("SSVHEM_WP") ) FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_SSVHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHP->at(i)>getPreCutValue1("SSVHPT_WP") ) FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_SSVHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPL_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_JPL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPM_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_JPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPT_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_JPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVL_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_CSVL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVM_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_CSVM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVT_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_C_CSVT", PFJetPt->at(i), PFJetEta->at(i));
           }
         }
         else if( partonFlavor==21 )
         {
           if(i==0)
           {
             FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_All", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEL_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_TCHEL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEM_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_TCHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPM_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_TCHPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPT_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_TCHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHE->at(i)>getPreCutValue1("SSVHEM_WP") ) FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_SSVHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHP->at(i)>getPreCutValue1("SSVHPT_WP") ) FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_SSVHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPL_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_JPL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPM_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_JPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPT_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_JPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVL_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_CSVL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVM_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_CSVM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVT_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_G_CSVT", PFJetPt->at(i), PFJetEta->at(i));
           }
           if(i==1)
           {
             FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_All", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEL_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_TCHEL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEM_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_TCHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPM_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_TCHPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPT_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_TCHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHE->at(i)>getPreCutValue1("SSVHEM_WP") ) FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_SSVHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHP->at(i)>getPreCutValue1("SSVHPT_WP") ) FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_SSVHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPL_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_JPL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPM_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_JPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPT_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_JPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVL_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_CSVL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVM_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_CSVM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVT_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_G_CSVT", PFJetPt->at(i), PFJetEta->at(i));
           }
         }
         else
         {
           if(i==0)
           {
             FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_All", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEL_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_TCHEL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEM_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_TCHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPM_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_TCHPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPT_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_TCHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHE->at(i)>getPreCutValue1("SSVHEM_WP") ) FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_SSVHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHP->at(i)>getPreCutValue1("SSVHPT_WP") ) FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_SSVHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPL_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_JPL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPM_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_JPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPT_WP") )       FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_JPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVL_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_CSVL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVM_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_CSVM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVT_WP") )      FillUserTH2D("h2_EtaJ1_vs_PtJ1_UDS_CSVT", PFJetPt->at(i), PFJetEta->at(i));
           }
           if(i==1)
           {
             FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_All", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEL_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_TCHEL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHE->at(i)>getPreCutValue1("TCHEM_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_TCHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPM_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_TCHPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetTCHP->at(i)>getPreCutValue1("TCHPT_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_TCHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHE->at(i)>getPreCutValue1("SSVHEM_WP") ) FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_SSVHEM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetSSVHP->at(i)>getPreCutValue1("SSVHPT_WP") ) FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_SSVHPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPL_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_JPL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPM_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_JPM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("JPT_WP") )       FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_JPT", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVL_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_CSVL", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVM_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_CSVM", PFJetPt->at(i), PFJetEta->at(i));
             if( PFJetJP->at(i)>getPreCutValue1("CSVT_WP") )      FillUserTH2D("h2_EtaJ2_vs_PtJ2_UDS_CSVT", PFJetPt->at(i), PFJetEta->at(i));
           }
         }
       }
     }
   }
   
   // select only those events that pass the full selection
   if( passedCut("all") ) ret = true;

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


DEFINE_FWK_MODULE(MyAnalyzer);
