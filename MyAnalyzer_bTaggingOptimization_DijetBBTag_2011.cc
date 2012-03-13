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
// $Id: MyAnalyzer_bTaggingOptimization_DijetBBTag_2011.cc,v 1.1 2012/02/23 20:38:10 ferencek Exp $
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
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass0to500GeV_eff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass0to500GeV_eff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass0to500GeV_miseff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass0to500GeV_miseff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass4000to1000000GeV_eff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass4000to1000000GeV_eff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass4000to1000000GeV_miseff_denom", 9, 0.5, 9.5, 9, 0.5, 9.5);
   CreateUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass4000to1000000GeV_miseff_num", 9, 0.5, 9.5, 9, 0.5, 9.5);
   
   
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass0to500GeV_eff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass0to500GeV_eff_num", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass0to500GeV_miseff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass0to500GeV_miseff_num", 9, 0.5, 9.5);
   
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_eff_num", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass500to1000GeV_miseff_num", 9, 0.5, 9.5);
   
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_eff_num", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass1000to1500GeV_miseff_num", 9, 0.5, 9.5);
   
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_eff_num", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass1500to2000GeV_miseff_num", 9, 0.5, 9.5);
   
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_eff_num", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass2000to2500GeV_miseff_num", 9, 0.5, 9.5);
   
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_eff_num", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass2500to3000GeV_miseff_num", 9, 0.5, 9.5);
   
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_eff_num", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass3000to3500GeV_miseff_num", 9, 0.5, 9.5);
   
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_eff_num", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass3500to4000GeV_miseff_num", 9, 0.5, 9.5);
   
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass4000to1000000GeV_eff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass4000to1000000GeV_eff_num", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass4000to1000000GeV_miseff_denom", 9, 0.5, 9.5);
   CreateUserTH1D("h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass4000to1000000GeV_miseff_num", 9, 0.5, 9.5);
   
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

//    int btagger = int(getPreCutValue1("btagger"));
   int matchingType = int(getPreCutValue1("matchingType"));
   double matchingRadius = getPreCutValue1("matchingRadius");
   int doInclusiveTagging = int(getPreCutValue1("doInclusiveTagging"));
   
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
   
   auto_ptr<std::vector<double> >  PFJetPt ( new std::vector<double>() );
   auto_ptr<std::vector<double> >  PFJetE  ( new std::vector<double>() );

   for(size_t i=0; i<PFJetPt_->size(); i++)
   {
       double JES_ScaleFactor = 1.;
       if( !iEvent.isRealData() ) JES_ScaleFactor = 1. + getPreCutValue1("JES_Shift")*PFJetUnc->at(i);

       PFJetPt->push_back( PFJetPt_->at(i)*JES_ScaleFactor );
       PFJetE ->push_back( PFJetE_ ->at(i)*JES_ScaleFactor );
   }

//    int nBTaggedJets = 0;
   int nHeavyFlavorJets = 0;
//    int nBTaggedHeavyFlavorJets = 0;
   
   // Set the evaluation of the cuts to false and clear the variable values and filled status
   resetCuts();

   int nSt3_b = 0, nSt3_b_fromRSG = 0;

   for(size_t i=0; i<GenParticlePt->size(); i++)
   {
     if( abs(GenParticlePdgId->at(i))==5 && GenParticleStatus->at(i)==3 && GenParticleMotherIndex->at(i)>=0 )
     {
       ++nSt3_b;
       if( abs(GenParticlePdgId->at(GenParticleMotherIndex->at(i)))==5000039 ) ++nSt3_b_fromRSG;
     }
   }

   fillVariableWithValue( "nSt3_b", nSt3_b );
   
   fillVariableWithValue( "nSt3_b_fromRSG", nSt3_b_fromRSG );

   fillVariableWithValue( "nJets", PFJetPt->size() );

   if( PFJetPt->size() >= 1 )
   {
     fillVariableWithValue( "passJetIdJ1", ( PFJetPassJetID->at(0) ? 1 : 0 ) );
     fillVariableWithValue( "absEtaJ1", fabs( PFJetEta->at(0) ));
   }
   if( PFJetPt->size() >= 2 )
   {
     fillVariableWithValue( "passJetIdJ2", ( PFJetPassJetID->at(1) ? 1 : 0 ) );
     fillVariableWithValue( "absEtaJ2", fabs( PFJetEta->at(1) ) );
     
     TLorentzVector v_j1j2, v_j1, v_j2;
     v_j1.SetPtEtaPhiE(PFJetPt->at(0),PFJetEta->at(0),PFJetPhi->at(0),PFJetE->at(0));
     v_j2.SetPtEtaPhiE(PFJetPt->at(1),PFJetEta->at(1),PFJetPhi->at(1),PFJetE->at(1));
     // calculate |DeltaEta(j1,j2)|
     fillVariableWithValue( "absDeltaEtaJ1J2", fabs( PFJetEta->at(0) - PFJetEta->at(1) ) );
     // calculate M_j1j2
     v_j1j2 = v_j1 + v_j2;
     fillVariableWithValue( "DijetMass", v_j1j2.M() );

     // jet and GenParticle 4-vectors
     TLorentzVector v_j, v_gp;

     // loop over two leading jets
     for(size_t i=0; i<2; ++i)
     {
       bool isHeavyFlavor = false;

       // set jet 4-vector
       v_j.SetPtEtaPhiE(PFJetPt->at(i),PFJetEta->at(i),PFJetPhi->at(i),PFJetE->at(i));

       if( !iEvent.isRealData() )
       {
         if( matchingType==0 && abs(PFJetPartonFlavor->at(i))==5 )
         {
           ++nHeavyFlavorJets;
           isHeavyFlavor = true;
         }
         else if( matchingType!=0 )
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
             isHeavyFlavor = true;
           }
         }
       }
        
//        if( (btagger==0 && PFJetTCHE->at(i)>getPreCutValue1("TCHEM_WP")) ||
//            (btagger==1 && PFJetSSVHE->at(i)>getPreCutValue1("SSVHEM_WP")) ||
//            (btagger==2 && PFJetTCHP->at(i)>getPreCutValue1("TCHPT_WP")) ||
//            (btagger==3 && PFJetSSVHP->at(i)>getPreCutValue1("SSVHPT_WP")) )
//        {
//          ++nBTaggedJets;
//          if( isHeavyFlavor ) ++nBTaggedHeavyFlavorJets;
//        }
     }
   }
   
   // Evaluate cuts (but do not apply them)
   evaluateCuts();
   
   for (int min=0; min<=4000; min+=500){
   int max = min+500;
   //int s = 0;
   	//if(max==4500){max="Inf";s=1;} //max is an int...
   	if(max==4500){max=1000000;}
	std::stringstream junk;
	junk << "h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass" << min << "to" << max << "GeV_miseff_num";
	string nummiseff = junk.str();
	junk.str("");
	junk << "h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass" << min << "to" << max << "GeV_eff_num";
	string numeff = junk.str();
	junk.str("");
	junk << "h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass" << min << "to" << max << "GeV_miseff_denom";
	string denommiseff = junk.str();
	junk.str("");
	junk << "h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass" << min << "to" << max << "GeV_eff_denom";
	string denomeff = junk.str();
	junk.str("");
	
	junk << "h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass" << min << "to" << max << "GeV_miseff_num";
	string nummiseff1 = junk.str();
	junk.str("");
	junk << "h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass" << min << "to" << max << "GeV_eff_num";
	string numeff1 = junk.str();
	junk.str("");
	junk << "h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass" << min << "to" << max << "GeV_miseff_denom";
	string denommiseff1 = junk.str();
	junk.str("");
	junk << "h1_TCHE_TCHP_SSVHE_SSVHP_DijetMass" << min << "to" << max << "GeV_eff_denom";
	string denomeff1 = junk.str();
	junk.str("");
	
	string titlenum[]={ nummiseff, "", numeff };
	string titledenom[]={ denommiseff, "", denomeff };
	
	string titlenum1[]={ nummiseff1, "", numeff1 };
	string titledenom1[]={ denommiseff1, "", denomeff1 };
   for (int n=0; n<=2; n+=2){	   
	   if(passedAllPreviousCuts("DijetMass") && getVariableValue("DijetMass")>min && getVariableValue("DijetMass")<max && ( doInclusiveTagging ? true : nHeavyFlavorJets==n ))
	   {
		for (int j=1; j<=9; j++){
			for (int i=1; i<=j; i++){
				FillUserTH2D(titledenom[n], j, i);
			}
		}
	//       FillUserTH2D("h2_TCHE_TCHP_SSVHE_SSVHP_DijetMass1to1p5TeV_denom", 1, 1); // (Original broad fill left as example to me.)

		//FILL COL 1 OF NUM
	       if( PFJetTCHE->at(0)>getPreCutValue1("TCHEL_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEL_WP") ){
			FillUserTH2D(titlenum[n], 1, 1);
		}
		 
		//FILL COL 2 OF NUM (BOTTOM UP)
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEL_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEM_WP") ) ||
		   ( PFJetTCHE->at(0)>getPreCutValue1("TCHEM_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEL_WP") ) ){
			FillUserTH2D(titlenum[n], 2, 1);
	       }
	       if( PFJetTCHE->at(0)>getPreCutValue1("TCHEM_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEM_WP") ){
			FillUserTH2D(titlenum[n], 2, 2);
		}
		
		//FILL COL 3 OF NUM (BOTTOM UP)
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEL_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHET_WP") ) ||
		   ( PFJetTCHE->at(0)>getPreCutValue1("TCHET_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEL_WP") ) ){
			FillUserTH2D(titlenum[n], 3, 1);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEM_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHET_WP") ) ||
		   ( PFJetTCHE->at(0)>getPreCutValue1("TCHET_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEM_WP") ) ){
			FillUserTH2D(titlenum[n], 3, 2);
	       }
	       if( PFJetTCHE->at(0)>getPreCutValue1("TCHET_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHET_WP") ){
			FillUserTH2D(titlenum[n], 3, 3);
		}
		
		//FILL COL 4 OF NUM (BOTTOM UP)
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEL_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPL_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPL_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEL_WP") ) ){
			FillUserTH2D(titlenum[n], 4, 1);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEM_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPL_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPL_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEM_WP") ) ){
			FillUserTH2D(titlenum[n], 4, 2);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHET_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPL_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPL_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHET_WP") ) ){
			FillUserTH2D(titlenum[n], 4, 3);
	       }
	       if( PFJetTCHP->at(0)>getPreCutValue1("TCHPL_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPL_WP") ){
			FillUserTH2D(titlenum[n], 4, 4);
	       } 
	       
	       	//FILL COL 5 OF NUM (BOTTOM UP)
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEL_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPM_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPM_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEL_WP") ) ){
			FillUserTH2D(titlenum[n], 5, 1);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEM_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPM_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPM_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEM_WP") ) ){
			FillUserTH2D(titlenum[n], 5, 2);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHET_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPM_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPM_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHET_WP") ) ){
			FillUserTH2D(titlenum[n], 5, 3);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPL_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPM_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPM_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPL_WP") ) ){
			FillUserTH2D(titlenum[n], 5, 4);
	       }
	       if( PFJetTCHP->at(0)>getPreCutValue1("TCHPM_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPM_WP") ){
			FillUserTH2D(titlenum[n], 5, 5);
	       } 
	       
	       	//FILL COL 6 OF NUM (BOTTOM UP)
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEL_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPT_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPT_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEL_WP") ) ){
			FillUserTH2D(titlenum[n], 6, 1);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEM_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPT_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPT_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEM_WP") ) ){
			FillUserTH2D(titlenum[n], 6, 2);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHET_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPT_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPT_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHET_WP") ) ){
			FillUserTH2D(titlenum[n], 6, 3);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPL_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPT_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPT_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPL_WP") ) ){
			FillUserTH2D(titlenum[n], 6, 4);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPM_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPT_WP") ) ||
		   ( PFJetTCHP->at(0)>getPreCutValue1("TCHPT_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPM_WP") ) ){
			FillUserTH2D(titlenum[n], 6, 5);
	       }
	       if( PFJetTCHP->at(0)>getPreCutValue1("TCHPT_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPT_WP") ){
			FillUserTH2D(titlenum[n], 6, 6);
	       } 
	       
	       	//FILL COL 7 OF NUM (BOTTOM UP)
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEL_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHEM_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHEM_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEL_WP") ) ){
			FillUserTH2D(titlenum[n], 7, 1);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEM_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHEM_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHEM_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEM_WP") ) ){
			FillUserTH2D(titlenum[n], 7, 2);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHET_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHEM_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHEM_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHET_WP") ) ){
			FillUserTH2D(titlenum[n], 7, 3);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPL_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHEM_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHEM_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPL_WP") ) ){
			FillUserTH2D(titlenum[n], 7, 4);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPM_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHEM_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHEM_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPM_WP") ) ){
			FillUserTH2D(titlenum[n], 7, 5);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPT_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHEM_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHEM_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPT_WP") ) ){
			FillUserTH2D(titlenum[n], 7, 6);
	       }
	       if( PFJetSSVHE->at(0)>getPreCutValue1("SSVHEM_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHEM_WP") ){
			FillUserTH2D(titlenum[n], 7, 7);
	       } 
	       
	       	//FILL COL 8 OF NUM (BOTTOM UP)
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEL_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHET_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHET_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEL_WP") ) ){
			FillUserTH2D(titlenum[n], 8, 1);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEM_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHET_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHET_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEM_WP") ) ){
			FillUserTH2D(titlenum[n], 8, 2);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHET_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHET_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHET_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHET_WP") ) ){
			FillUserTH2D(titlenum[n], 8, 3);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPL_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHET_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHET_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPL_WP") ) ){
			FillUserTH2D(titlenum[n], 8, 4);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPM_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHET_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHET_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPM_WP") ) ){
			FillUserTH2D(titlenum[n], 8, 5);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPT_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHET_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHET_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPT_WP") ) ){
			FillUserTH2D(titlenum[n], 8, 6);
	       }
	       if( ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHEM_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHET_WP") ) ||
		   ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHET_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHEM_WP") ) ){
			FillUserTH2D(titlenum[n], 8, 7);
	       }
	       if( PFJetSSVHE->at(0)>getPreCutValue1("SSVHET_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHET_WP") ){
			FillUserTH2D(titlenum[n], 8, 8);
	       } 
	       
	       	//FILL COL 9 OF NUM (BOTTOM UP)
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEL_WP") && PFJetSSVHP->at(1)>getPreCutValue1("SSVHPT_WP") ) ||
		   ( PFJetSSVHP->at(0)>getPreCutValue1("SSVHPT_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEL_WP") ) ){
			FillUserTH2D(titlenum[n], 9, 1);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHEM_WP") && PFJetSSVHP->at(1)>getPreCutValue1("SSVHPT_WP") ) ||
		   ( PFJetSSVHP->at(0)>getPreCutValue1("SSVHPT_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEM_WP") ) ){
			FillUserTH2D(titlenum[n], 9, 2);
	       }
	       if( ( PFJetTCHE->at(0)>getPreCutValue1("TCHET_WP") && PFJetSSVHP->at(1)>getPreCutValue1("SSVHPT_WP") ) ||
		   ( PFJetSSVHP->at(0)>getPreCutValue1("SSVHPT_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHET_WP") ) ){
			FillUserTH2D(titlenum[n], 9, 3);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPL_WP") && PFJetSSVHP->at(1)>getPreCutValue1("SSVHPT_WP") ) ||
		   ( PFJetSSVHP->at(0)>getPreCutValue1("SSVHPT_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPL_WP") ) ){
			FillUserTH2D(titlenum[n], 9, 4);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPM_WP") && PFJetSSVHP->at(1)>getPreCutValue1("SSVHPT_WP") ) ||
		   ( PFJetSSVHP->at(0)>getPreCutValue1("SSVHPT_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPM_WP") ) ){
			FillUserTH2D(titlenum[n], 9, 5);
	       }
	       if( ( PFJetTCHP->at(0)>getPreCutValue1("TCHPT_WP") && PFJetSSVHP->at(1)>getPreCutValue1("SSVHPT_WP") ) ||
		   ( PFJetSSVHP->at(0)>getPreCutValue1("SSVHPT_WP") && PFJetTCHP->at(1)>getPreCutValue1("TCHPT_WP") ) ){
			FillUserTH2D(titlenum[n], 9, 6);
	       }
	       if( ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHEM_WP") && PFJetSSVHP->at(1)>getPreCutValue1("SSVHPT_WP") ) ||
		   ( PFJetSSVHP->at(0)>getPreCutValue1("SSVHPT_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHEM_WP") ) ){
			FillUserTH2D(titlenum[n], 9, 7);
	       }
	       if( ( PFJetSSVHE->at(0)>getPreCutValue1("SSVHET_WP") && PFJetSSVHP->at(1)>getPreCutValue1("SSVHPT_WP") ) ||
		   ( PFJetSSVHP->at(0)>getPreCutValue1("SSVHPT_WP") && PFJetSSVHE->at(1)>getPreCutValue1("SSVHET_WP") ) ){
			FillUserTH2D(titlenum[n], 9, 8);
	       }
	       if( PFJetSSVHP->at(0)>getPreCutValue1("SSVHPT_WP") && PFJetSSVHP->at(1)>getPreCutValue1("SSVHPT_WP") ){
			FillUserTH2D(titlenum[n], 9, 9);
	       }   
	       
	       //FILL 1D HISTOGRAMS
	       
	       for(int m=1; m<=9; m+=1){
	       		FillUserTH1D(titledenom1[n], m);
	       }	
       		if( (PFJetTCHE->at(0)>getPreCutValue1("TCHEL_WP") || PFJetTCHE->at(1)>getPreCutValue1("TCHEL_WP")) &&
       		!(PFJetTCHE->at(0)>getPreCutValue1("TCHEL_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEL_WP"))){
			FillUserTH1D(titlenum1[n], 1);
		}
		if( (PFJetTCHE->at(0)>getPreCutValue1("TCHEM_WP") || PFJetTCHE->at(1)>getPreCutValue1("TCHEM_WP")) &&
       		!(PFJetTCHE->at(0)>getPreCutValue1("TCHEM_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHEM_WP"))){
			FillUserTH1D(titlenum1[n], 2);
		}
		if( (PFJetTCHE->at(0)>getPreCutValue1("TCHET_WP") || PFJetTCHE->at(1)>getPreCutValue1("TCHET_WP")) &&
       		!(PFJetTCHE->at(0)>getPreCutValue1("TCHET_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHET_WP"))){
			FillUserTH1D(titlenum1[n], 3);
		}
		if( (PFJetTCHE->at(0)>getPreCutValue1("TCHPL_WP") || PFJetTCHE->at(1)>getPreCutValue1("TCHPL_WP")) &&
       		!(PFJetTCHE->at(0)>getPreCutValue1("TCHPL_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHPL_WP"))){
			FillUserTH1D(titlenum1[n], 4);
		}
		if( (PFJetTCHE->at(0)>getPreCutValue1("TCHPM_WP") || PFJetTCHE->at(1)>getPreCutValue1("TCHPM_WP")) &&
       		!(PFJetTCHE->at(0)>getPreCutValue1("TCHPM_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHPM_WP"))){
			FillUserTH1D(titlenum1[n], 5);
		}
		if( (PFJetTCHE->at(0)>getPreCutValue1("TCHPT_WP") || PFJetTCHE->at(1)>getPreCutValue1("TCHPT_WP")) &&
       		!(PFJetTCHE->at(0)>getPreCutValue1("TCHPT_WP") && PFJetTCHE->at(1)>getPreCutValue1("TCHPT_WP"))){
			FillUserTH1D(titlenum1[n], 6);
		}
		if( (PFJetTCHE->at(0)>getPreCutValue1("SSVHEM_WP") || PFJetTCHE->at(1)>getPreCutValue1("SSVHEM_WP")) &&
       		!(PFJetTCHE->at(0)>getPreCutValue1("SSVHEM_WP") && PFJetTCHE->at(1)>getPreCutValue1("SSVHEM_WP"))){
			FillUserTH1D(titlenum1[n], 7);
		}
		if( (PFJetTCHE->at(0)>getPreCutValue1("SSVHET_WP") || PFJetTCHE->at(1)>getPreCutValue1("SSVHET_WP")) &&
       		!(PFJetTCHE->at(0)>getPreCutValue1("SSVHET_WP") && PFJetTCHE->at(1)>getPreCutValue1("SSVHET_WP"))){
			FillUserTH1D(titlenum1[n], 8);
		}
		if( (PFJetTCHE->at(0)>getPreCutValue1("SSVHPT_WP") || PFJetTCHE->at(1)>getPreCutValue1("SSVHPT_WP")) &&
       		!(PFJetTCHE->at(0)>getPreCutValue1("SSVHPT_WP") && PFJetTCHE->at(1)>getPreCutValue1("SSVHPT_WP"))){
			FillUserTH1D(titlenum1[n], 9);
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
