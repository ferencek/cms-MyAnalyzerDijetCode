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
// $Id: MyAnalyzer_ZprimeAnalysis_DijetBBTag_2011.cc,v 1.2 2011/12/15 02:16:55 ferencek Exp $
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
      int bQuark_0_counter, bQuark_st2_1_counter, bQuark_st2_0_bcHadron_2_counter, bQuark_st2_2_bcHadron_0_counter, bQuark_st2_0_st3_2_counter;
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
   CreateUserTH1D("h1_J1J2PartonFlavor;Parton Flavor (PDG ID);Entries", 51, -0.5, 50.5);
   CreateUserTH1D("h1_minDeltaR_bQuark_status2;min#DeltaR;Entries", 100, 0., 10.);
   CreateUserTH1D("h1_minDeltaR_bQuark_status3;min#DeltaR;Entries", 100, 0., 10.);
   CreateUserTH1D("h1_minDeltaR_bQuark;min#DeltaR;Entries", 100, 0., 10.);
   CreateUserTH1D("h1_minDeltaR_bHadron;min#DeltaR;Entries", 100, 0., 10.);
   CreateUserTH1D("h1_minDeltaR_cHadron;min#DeltaR;Entries", 100, 0., 10.);
   CreateUserTH1D("h1_minDeltaR_bcHadron;min#DeltaR;Entries", 100, 0., 10.);
   
   // initialize your variables here
   bQuark_0_counter = 0;
   bQuark_st2_1_counter = 0;
   bQuark_st2_0_bcHadron_2_counter = 0;
   bQuark_st2_2_bcHadron_0_counter = 0;
   bQuark_st2_0_st3_2_counter = 0;
   
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
   
   double matchingRadius = getPreCutValue1("matchingRadius");
   int matchingType = int(getPreCutValue1("matchingType"));
   int eventPrintout = int(getPreCutValue1("eventPrintout"));
   
   // grab necessary objects from the event
   edm::Handle<vector<unsigned int> > NPU;
   iEvent.getByLabel(edm::InputTag("GenEventInfo:PileUpNumberOfInteractions"), NPU);
   edm::Handle<vector<int> > BX;
   iEvent.getByLabel(edm::InputTag("GenEventInfo:PileUpBunchCrossing"), BX);

   edm::Handle<bool> passPrimaryVertex;
   iEvent.getByLabel(edm::InputTag("EventSelection:PassPrimaryVertex"), passPrimaryVertex);
   edm::Handle<bool> passBeamScraping;
   iEvent.getByLabel(edm::InputTag("EventSelection:PassBeamScraping"), passBeamScraping);

   edm::Handle<vector<double> > PFJetPt;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Pt"), PFJetPt);
   edm::Handle<vector<double> > PFJetEta;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Eta"), PFJetEta);
   edm::Handle<vector<double> > PFJetPhi;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Phi"), PFJetPhi);
   edm::Handle<vector<double> > PFJetE;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Energy"), PFJetE);
   edm::Handle<vector<int> > PFJetPassJetID;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:PassTightID"), PFJetPassJetID);
   edm::Handle<vector<int> > PFJetPartonFlavor;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:PartonFlavor"), PFJetPartonFlavor);

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

   
   int case1_counter = 0, case2_counter = 0, case3_counter = 0, case4_counter = 0, case5_counter = 0;
   
   // Set the evaluation of the cuts to false and clear the variable values and filled status
   resetCuts();
   
   fillVariableWithValue("PassPrimaryVertex", ( *passPrimaryVertex ? 1 : 0 ), eventWeight );
   fillVariableWithValue("PassBeamScraping", ( *passBeamScraping ? 1 : 0 ), eventWeight );

   fillVariableWithValue( "nJets", PFJetPt->size(), eventWeight );

   if( PFJetPt->size() >= 1 )
   {
     fillVariableWithValue( "passJetIdJ1", ( PFJetPassJetID->at(0) ? 1 : 0 ), eventWeight );
     fillVariableWithValue( "absEtaJ1", fabs( PFJetEta->at(0) ), eventWeight );
   }
   if( PFJetPt->size() >= 2 )
   {
     fillVariableWithValue( "passJetIdJ2", ( PFJetPassJetID->at(1) ? 1 : 0 ), eventWeight );
     fillVariableWithValue( "absEtaJ2", fabs( PFJetEta->at(1) ), eventWeight );
     // calculate |DeltaEta(j1,j2)|
     fillVariableWithValue( "absDeltaEtaJ1J2", fabs( PFJetEta->at(0) - PFJetEta->at(1) ), eventWeight );
    
     int nJets_HeavyFlavor = 0;
     TLorentzVector v_j, v_gp;
     
     // loop over two leading jets
     for(size_t i=0; i<2; ++i)
     {
       FillUserTH1D("h1_J1J2PartonFlavor", abs( PFJetPartonFlavor->at(i) ), eventWeight );
      
       if( matchingType==0 && abs(PFJetPartonFlavor->at(i))==5)
         ++nJets_HeavyFlavor;
       else if( matchingType!=0 )
       {
         double minDeltaR_bQuark_status2 = 999.;
         double minDeltaR_bQuark_status3 = 999.;
         double minDeltaR_bHadron = 999.;
         double minDeltaR_cHadron = 999.;

         // initialize jet 4-vector
         v_j.SetPtEtaPhiE(PFJetPt->at(i),PFJetEta->at(i),PFJetPhi->at(i),PFJetE->at(i));

         // loop over GenParticles
         for(size_t j=0; j<GenParticlePt->size(); ++j)
         {
           if( GenParticlePt->at(j)==0 ) continue;

           int pdgID = abs(GenParticlePdgId->at(j));
           int status = abs(GenParticleStatus->at(j));
           bool shouldBeMatched_bQuark_status2 = false;
           bool shouldBeMatched_bQuark_status3 = false;
           bool shouldBeMatched_bHadron = false;
           bool shouldBeMatched_cHadron = false;

           if( pdgID==5
               && status==2
              ) shouldBeMatched_bQuark_status2 = true;
           if( pdgID==5
               && status==3
              ) shouldBeMatched_bQuark_status3 = true;
           if( pdgID==511 || pdgID==521 || pdgID==531 || pdgID==541 || pdgID==5122 || pdgID==5132 || pdgID==5232 || pdgID==5332
              ) shouldBeMatched_bHadron = true;
           if( pdgID==411 || pdgID==421 || pdgID==431 || pdgID==4122 || pdgID==4132 || pdgID==4232 || pdgID==4332
              ) shouldBeMatched_cHadron = true;

           
           if( shouldBeMatched_bQuark_status2 )
           {
             v_gp.SetPtEtaPhiE(GenParticlePt->at(j),GenParticleEta->at(j),GenParticlePhi->at(j),GenParticleE->at(j));
             double deltaR = v_j.DeltaR(v_gp);

             if( deltaR < minDeltaR_bQuark_status2 ) minDeltaR_bQuark_status2 = deltaR;
           }
           if( shouldBeMatched_bQuark_status3 )
           {
             v_gp.SetPtEtaPhiE(GenParticlePt->at(j),GenParticleEta->at(j),GenParticlePhi->at(j),GenParticleE->at(j));
             double deltaR = v_j.DeltaR(v_gp);

             if( deltaR < minDeltaR_bQuark_status3 ) minDeltaR_bQuark_status3 = deltaR;
           }
           if( shouldBeMatched_bHadron )
           {
             v_gp.SetPtEtaPhiE(GenParticlePt->at(j),GenParticleEta->at(j),GenParticlePhi->at(j),GenParticleE->at(j));
             double deltaR = v_j.DeltaR(v_gp);

             if( deltaR < minDeltaR_bHadron ) minDeltaR_bHadron = deltaR;
           }
           if( shouldBeMatched_cHadron )
           {
             v_gp.SetPtEtaPhiE(GenParticlePt->at(j),GenParticleEta->at(j),GenParticlePhi->at(j),GenParticleE->at(j));
             double deltaR = v_j.DeltaR(v_gp);

             if( deltaR < minDeltaR_cHadron ) minDeltaR_cHadron = deltaR;
           }
         }
         
         FillUserTH1D("h1_minDeltaR_bQuark_status2", minDeltaR_bQuark_status2, eventWeight);
         FillUserTH1D("h1_minDeltaR_bQuark_status3", minDeltaR_bQuark_status3, eventWeight);
         FillUserTH1D("h1_minDeltaR_bQuark", min(minDeltaR_bQuark_status2,minDeltaR_bQuark_status3), eventWeight);
         FillUserTH1D("h1_minDeltaR_bHadron", minDeltaR_bHadron, eventWeight);
         FillUserTH1D("h1_minDeltaR_cHadron", minDeltaR_cHadron, eventWeight);
         FillUserTH1D("h1_minDeltaR_bcHadron", min(minDeltaR_bHadron,minDeltaR_cHadron), eventWeight);

         if( minDeltaR_bQuark_status2 < matchingRadius ) ++case1_counter;
         if( minDeltaR_bQuark_status3 < matchingRadius ) ++case2_counter;
         if( minDeltaR_bQuark_status2 < matchingRadius || minDeltaR_bQuark_status3 < matchingRadius ) ++case3_counter;
         if( minDeltaR_bHadron < matchingRadius ) ++case4_counter;
         if( minDeltaR_bHadron < matchingRadius || minDeltaR_cHadron < matchingRadius ) ++case5_counter;
         
         switch( matchingType )
         {
           case 1:
             if( minDeltaR_bQuark_status2 < matchingRadius ) ++nJets_HeavyFlavor;
             break;
           case 2:
             if( minDeltaR_bQuark_status3 < matchingRadius ) ++nJets_HeavyFlavor;
             break;
           case 3:
             if( minDeltaR_bQuark_status2 < matchingRadius || minDeltaR_bQuark_status3 < matchingRadius ) ++nJets_HeavyFlavor;
             break;
           case 4:
             if( minDeltaR_bHadron < matchingRadius ) ++nJets_HeavyFlavor;
             break;
           case 5:
             if( minDeltaR_bHadron < matchingRadius || minDeltaR_cHadron < matchingRadius ) ++nJets_HeavyFlavor;
             break;
         }
       }
     }
 
     fillVariableWithValue( "nJets_HeavyFlavor", nJets_HeavyFlavor, eventWeight );
   }
   
   // Evaluate cuts (but do not apply them)
   evaluateCuts();

   // select only those events that pass the full selection
   if( passedCut("all") ) ret = true;

   if( eventPrintout )
   {
     if( case3_counter==0 )
     {
       ++bQuark_0_counter;
       if( bQuark_0_counter<11 ) cout << ">> bQuark_0 --> Run:LS:Event = " << iEvent.id().run() << ":" << iEvent.luminosityBlock() << ":" << iEvent.id().event() << endl;
     }

     if( case1_counter==1 )
     {
       ++bQuark_st2_1_counter;
       if( bQuark_st2_1_counter<11 ) cout << ">> bQuark_st2_1 --> Run:LS:Event = " << iEvent.id().run() << ":" << iEvent.luminosityBlock() << ":" << iEvent.id().event() << endl;
     }

     if( case1_counter==0 && case5_counter==2 )
     {
       ++bQuark_st2_0_bcHadron_2_counter;
       if( bQuark_st2_0_bcHadron_2_counter<11 ) cout << ">> bQuark_st2_0_bcHadron_2 --> Run:LS:Event = " << iEvent.id().run() << ":" << iEvent.luminosityBlock() << ":" << iEvent.id().event() << endl;
     }

     if( case1_counter==2 && case5_counter==0 )
     {
       ++bQuark_st2_2_bcHadron_0_counter;
       if( bQuark_st2_2_bcHadron_0_counter<11 ) cout << ">> bQuark_st2_2_bcHadron_0 --> Run:LS:Event = " << iEvent.id().run() << ":" << iEvent.luminosityBlock() << ":" << iEvent.id().event() << endl;
     }
     
     if( case1_counter==0 && case2_counter==2 )
     {
       ++bQuark_st2_0_st3_2_counter;
       if( bQuark_st2_0_st3_2_counter<11 ) cout << ">> bQuark_st2_0_st3_2 --> Run:LS:Event = " << iEvent.id().run() << ":" << iEvent.luminosityBlock() << ":" << iEvent.id().event() << endl;
     }
   }
   
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
   int eventPrintout = int(getPreCutValue1("eventPrintout"));
   
   if( eventPrintout )
   {
     cout << ">> bQuark_0: " << bQuark_0_counter << " events in total" << endl
          << ">> bQuark_st2_1: " << bQuark_st2_1_counter << " events in total" << endl
          << ">> bQuark_st2_0_bcHadron_2: " << bQuark_st2_0_bcHadron_2_counter << " events in total" << endl
          << ">> bQuark_st2_2_bcHadron_0: " << bQuark_st2_2_bcHadron_0_counter << " events in total" << endl
          << ">> bQuark_st2_0_st3_2: " << bQuark_st2_0_st3_2_counter << " events in total" << endl;
   }
 
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
