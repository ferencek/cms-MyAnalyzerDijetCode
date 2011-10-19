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
// $Id: MyAnalyzer_TriggerStudy.cc,v 1.4 2011/10/13 21:52:28 ferencek Exp $
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

// BaseClass
#include "MyAnalysis/MyAnalyzer/interface/BaseClass.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT
#include <TH1D.h>
#include <TH2D.h>
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
      HLTConfigProvider hltConfig;
      edm::InputTag     hltInputTag;
};

//
// constructors and destructor
//

MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig) :
  BaseClass(iConfig)
{
   //now do whatever initialization is needed
   hltInputTag = iConfig.getParameter<edm::InputTag>("HLTInputTag");
}

MyAnalyzer::~MyAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// constants, enums and typedefs
//
typedef vector<double>::const_iterator myiter;

struct ordering {
    bool operator ()(pair<size_t, myiter> const& a, pair<size_t, myiter> const& b) {
        return *a.second > *b.second;
    }
};

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

   // initialize your variables here
   
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
   bool ret = true;
   // event weight (by default set to 1)
   double eventWeight = 1;

   //#####################################################################################
   //########################### User's code starts here #################################

   // grab necessary objects from the event
   edm::Handle<edm::TriggerResults> triggerResults;
   iEvent.getByLabel(hltInputTag, triggerResults);
   
   edm::Handle<vector<double> > PFJetPt;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Pt"), PFJetPt);
   edm::Handle<vector<double> > PFJetPtRaw;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:PtRaw"), PFJetPtRaw);
   edm::Handle<vector<double> > PFJetEta;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Eta"), PFJetEta);
   edm::Handle<vector<double> > PFJetPhi;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Phi"), PFJetPhi);
   edm::Handle<vector<double> > PFJetE;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Energy"), PFJetE);
   edm::Handle<vector<int> > PFJetPassLooseID;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:PassLooseID"), PFJetPassLooseID);

   edm::Handle<vector<double> > CaloJetPt;
   iEvent.getByLabel(edm::InputTag("AK7CaloJets:Pt"), CaloJetPt);
   edm::Handle<vector<double> > CaloJetPtRaw;
   iEvent.getByLabel(edm::InputTag("AK7CaloJets:PtRaw"), CaloJetPtRaw);
   edm::Handle<vector<double> > CaloJetEta;
   iEvent.getByLabel(edm::InputTag("AK7CaloJets:Eta"), CaloJetEta);
   edm::Handle<vector<double> > CaloJetPhi;
   iEvent.getByLabel(edm::InputTag("AK7CaloJets:Phi"), CaloJetPhi);
   edm::Handle<vector<double> > CaloJetE;
   iEvent.getByLabel(edm::InputTag("AK7CaloJets:Energy"), CaloJetE);
   edm::Handle<vector<int> > CaloJetPassLooseID;
   iEvent.getByLabel(edm::InputTag("AK7CaloJets:PassLooseID"), CaloJetPassLooseID);

   string jet300_name = "HLT_Jet240_v1"; // HLT_Jet300 was not used before run 163269
   if( iEvent.id().run() >= 163269 ) jet300_name = "HLT_Jet300_v1";
   if( iEvent.id().run() >= 165088 ) jet300_name = "HLT_Jet300_v2";
   if( iEvent.id().run() >= 165970 ) jet300_name = "HLT_Jet300_v3";
   if( iEvent.id().run() == 166346 ) jet300_name = "HLT_Jet300_v4";
   if( iEvent.id().run() >= 167078 ) jet300_name = "HLT_Jet300_v5";
   if( iEvent.id().run() >= 176545 ) jet300_name = "HLT_Jet300_v6";
   
   string jet370_name = "HLT_Jet370_v1";
   if( iEvent.id().run() >= 163269 ) jet370_name = "HLT_Jet370_v2";
   if( iEvent.id().run() >= 165088 ) jet370_name = "HLT_Jet370_v3";
   if( iEvent.id().run() >= 165970 ) jet370_name = "HLT_Jet370_v4";
   if( iEvent.id().run() == 166346 ) jet370_name = "HLT_Jet370_v5";
   if( iEvent.id().run() >= 167078 ) jet370_name = "HLT_Jet370_v6";

   // check trigger
   int jetTrigger_fired = 0;
   int useHLTJet370 = (int)getPreCutValue1("useHLTJet370");
   
   string jet_name = jet300_name;
   if(useHLTJet370) jet_name = jet370_name;
   
   size_t index = hltConfig.triggerIndex(jet_name);
   if( index < triggerResults->size() )
   {
     if( triggerResults->accept( index ) ) jetTrigger_fired = 1;
   }
   else
   {
     edm::LogWarning("MyAnalyzer::filter") << "Requested HLT path \"" << jet_name << "\" does not exist";
   }

   // select PFJets that pass loose JetID
   vector<int> v_idx_pfjet_looseID;
   for(size_t i=0; i<PFJetPt->size(); i++)
   {
       if( PFJetPassLooseID->at(i) ) v_idx_pfjet_looseID.push_back(i);
   }

   // select CaloJets that pass loose JetID
   vector<int> v_idx_jet_looseID;
   for(size_t i=0; i<CaloJetPt->size(); i++)
   {
       if( CaloJetPassLooseID->at(i) ) v_idx_jet_looseID.push_back(i);
   }

   // select CaloJets with |eta|<getPreCutValue1("jetFidRegion")
   vector<int> v_idx_jet_fidRegion;
   for(size_t i=0; i<CaloJetPt->size(); i++)
   {
       if( CaloJetEta->at(i) > getPreCutValue1("jetFidRegion") ) continue;
       v_idx_jet_fidRegion.push_back(i);
   }

//    // create a vector of CaloJet indices (for CaloJets with |eta|<getPreCutValue1("jetFidRegion")) ordered by uncorrected jet pT
//    vector<int> v_idx_jet_uncorr;
// 
//    vector<pair<size_t, myiter> > order;
//    size_t n = 0;
//    for (myiter it = CaloJetPtRaw->begin(); it != CaloJetPtRaw->end(); ++it, ++n)
//    {
//      if( fabs(CaloJetEta->at(n)) > getPreCutValue1("jetFidRegion") ) continue;
//      order.push_back(make_pair(n, it));
//    }
//    sort(order.begin(), order.end(), ordering());
// 
//    for(size_t i=0; i<order.size(); ++i)
//    {
//      v_idx_jet_uncorr.push_back(order[i].first); // CaloJets with |eta|<getPreCutValue1("jetFidRegion") ordered by ucorrected jet pT
//    }
//    
//    // check if there are any jets present within |DeltaEta|<getPreCutValue1("jetTriggerDeltaEtaCut") from the leading uncorrected jet
//    int deltaEtaJetFound = 0;
//    if(v_idx_jet_uncorr.size() >= 2)
//    {
//      for(size_t i=1; i<v_idx_jet_uncorr.size(); ++i)
//      {
//        if( CaloJetPtRaw->at(v_idx_jet_uncorr[i]) < getPreCutValue1("jetTriggerPtThresh") ) continue;
//        if( fabs(CaloJetEta->at(v_idx_jet_uncorr[0]) - CaloJetEta->at(v_idx_jet_uncorr[i])) < getPreCutValue1("jetTriggerDeltaEtaCut"))
//        {
//          deltaEtaJetFound = 1;
//          break;
//        }
//      }
//    }

   // check if there is a jet present within |DeltaEta|<getPreCutValue1("jetTriggerDeltaEtaCut") from the leading jet
   // (both inside the fiducial region) that forms the invariant mass with the leading jet greater than getPreCutValue1("dijetMassCut")
   int triggerCuts = 0;
   if(v_idx_jet_fidRegion.size() >= 2)
   {
     TLorentzVector v_j1;
     v_j1.SetPtEtaPhiE(CaloJetPt->at(v_idx_jet_fidRegion[0]),CaloJetEta->at(v_idx_jet_fidRegion[0]),CaloJetPhi->at(v_idx_jet_fidRegion[0]),CaloJetE->at(v_idx_jet_fidRegion[0]));
     
     for(size_t i=1; i<v_idx_jet_fidRegion.size(); ++i)
     {
       if( fabs(CaloJetEta->at(v_idx_jet_fidRegion[0]) - CaloJetEta->at(v_idx_jet_fidRegion[i])) > getPreCutValue1("jetTriggerDeltaEtaCut")) continue;
       
       TLorentzVector v_ji, v_j1ji;
       v_ji.SetPtEtaPhiE(CaloJetPt->at(v_idx_jet_fidRegion[i]),CaloJetEta->at(v_idx_jet_fidRegion[i]),CaloJetPhi->at(v_idx_jet_fidRegion[i]),CaloJetE->at(v_idx_jet_fidRegion[i]));
       // calculate M_j1ji
       v_j1ji = v_j1 + v_ji;
       if( v_j1ji.M() > getPreCutValue1("dijetMassCut") )
       {
         triggerCuts = 1;
         break;
       }
     }
   }

   int passEEAnomJetFilter = 1;
   if( v_idx_jet_looseID.size() > 0 )
   {
     if( CaloJetPt->at(v_idx_jet_looseID[0]) > 15000 ) passEEAnomJetFilter = 0;
   }

   // Set the evaluation of the cuts to false and clear the variable values and filled status
   resetCuts();

   fillVariableWithValue("HLT_Jet", jetTrigger_fired );
   fillVariableWithValue("PassEEAnomJetFilter", passEEAnomJetFilter );
   fillVariableWithValue("TriggerCuts", triggerCuts );

   fillVariableWithValue("nJet_all", CaloJetPt->size());
   fillVariableWithValue("nJet_looseID", v_idx_jet_looseID.size());
   
   if( v_idx_jet_looseID.size() >= 1 )
   {
       fillVariableWithValue( "absEtaJet1", fabs( CaloJetEta->at(v_idx_jet_looseID[0]) ) );
   }
   if( v_idx_jet_looseID.size() >= 2 )
   {
       fillVariableWithValue( "absEtaJet2", fabs( CaloJetEta->at(v_idx_jet_looseID[1]) ) );
       
       TLorentzVector v_j1j2, v_j1, v_j2;
       v_j1.SetPtEtaPhiE(CaloJetPt->at(v_idx_jet_looseID[0]),CaloJetEta->at(v_idx_jet_looseID[0]),CaloJetPhi->at(v_idx_jet_looseID[0]),CaloJetE->at(v_idx_jet_looseID[0]));
       v_j2.SetPtEtaPhiE(CaloJetPt->at(v_idx_jet_looseID[1]),CaloJetEta->at(v_idx_jet_looseID[1]),CaloJetPhi->at(v_idx_jet_looseID[1]),CaloJetE->at(v_idx_jet_looseID[1]));
       // calculate |DeltaEta(j1,j2)|
       fillVariableWithValue( "absDeltaEtaJ1J2", fabs( CaloJetEta->at(v_idx_jet_looseID[0]) - CaloJetEta->at(v_idx_jet_looseID[1]) ) );
       // calculate M_j1j2
       v_j1j2 = v_j1 + v_j2;
       
       fillVariableWithValue( "J1J2Mass", v_j1j2.M() );
       fillVariableWithValue( "J1J2Mass500", getVariableValue("J1J2Mass") );
       fillVariableWithValue( "J1J2Mass600", getVariableValue("J1J2Mass") );
       fillVariableWithValue( "J1J2Mass700", getVariableValue("J1J2Mass") );
       fillVariableWithValue( "J1J2Mass800", getVariableValue("J1J2Mass") );
       fillVariableWithValue( "J1J2Mass900", getVariableValue("J1J2Mass") );
       fillVariableWithValue( "J1J2Mass1000", getVariableValue("J1J2Mass") );
       fillVariableWithValue( "J1J2Mass1100", getVariableValue("J1J2Mass") );
       fillVariableWithValue( "J1J2Mass1200", getVariableValue("J1J2Mass") );
       fillVariableWithValue( "J1J2Mass1300", getVariableValue("J1J2Mass") );
       fillVariableWithValue( "J1J2Mass1400", getVariableValue("J1J2Mass") );
       fillVariableWithValue( "J1J2Mass1500", getVariableValue("J1J2Mass") );
   }

   // Evaluate cuts (but do not apply them)
   evaluateCuts();

   // select only those events that pass the full selection
   ret = passedCut("all");

   // ##################### Event Printout - START #####################

//    // PassJ1J2Mass1000FailTrigger event printout
//    string category = "PassJ1J2Mass1000FailTrigger";
//    if( passedAllPreviousCuts("J1J2Mass1000") && passedCut("J1J2Mass1000") && deltaEtaJetFound==0 )
//    {
//      cout << category << ": ----------- START ------------" << endl;
//      cout << category << ": Run, lumi, event: "<< iEvent.id().run() << ", "
//                                                << iEvent.luminosityBlock() << ", "
//                                                << iEvent.id().event() << endl;
//      // loop over ak7 CaloJets
//      for (size_t i=0; i<CaloJetPt->size(); ++i)
//      {
//        cout << category << ": CaloJet " << i << " Pt, PtRaw, eta, phi: " << CaloJetPt->at(i) << ", "
//                                                                          << CaloJetPtRaw->at(i) << ", "
//                                                                          << CaloJetEta->at(i) << ", "
//                                                                          << CaloJetPhi->at(i) << endl;
//      }
//      // loop over ucorrected ak7 CaloJets with |eta|<2.7
//      for (size_t i=0; i<v_idx_jet_uncorr.size(); ++i)
//      {
//        cout << category << ": Uncorrected CaloJet " << i << " Pt, PtRaw, eta, phi: " << CaloJetPt->at(v_idx_jet_uncorr[i]) << ", "
//                                                                                      << CaloJetPtRaw->at(v_idx_jet_uncorr[i]) << ", "
//                                                                                      << CaloJetEta->at(v_idx_jet_uncorr[i]) << ", "
//                                                                                      << CaloJetPhi->at(v_idx_jet_uncorr[i]) << endl;
//      }
//      // loop over ak7 CaloJets that pass loose jet ID requirements
//      for (size_t i=0; i<v_idx_jet_looseID.size(); ++i)
//      {
//        cout << category << ": PassLooseID CaloJet "<< i << " Pt, PtRaw, eta, phi: " << CaloJetPt->at(v_idx_jet_looseID[i]) << ", "
//                                                                                     << CaloJetPtRaw->at(v_idx_jet_looseID[i]) << ", "
//                                                                                     << CaloJetEta->at(v_idx_jet_looseID[i]) << ", "
//                                                                                     << CaloJetPhi->at(v_idx_jet_looseID[i]) << endl;
//      }
//      // loop over ak7 PFJets
//      for (size_t i=0; i<PFJetPt->size(); ++i)
//      {
//        cout << category << ": PFJet " << i << " Pt, PtRaw, eta, phi: " << PFJetPt->at(i) << ", "
//                                                                        << PFJetPtRaw->at(i) << ", "
//                                                                        << PFJetEta->at(i) << ", "
//                                                                        << PFJetPhi->at(i) << endl;
//      }
//      // loop over ak7 PFJets that pass loose jet ID requirements
//      for (size_t i=0; i<v_idx_pfjet_looseID.size(); ++i)
//      {
//        cout << category << ": PassLooseID PFJet "<< i << " Pt, PtRaw, eta, phi: " << PFJetPt->at(v_idx_jet_looseID[i]) << ", "
//                                                                                   << PFJetPtRaw->at(v_idx_jet_looseID[i]) << ", "
//                                                                                   << PFJetEta->at(v_idx_jet_looseID[i]) << ", "
//                                                                                   << PFJetPhi->at(v_idx_jet_looseID[i]) << endl;
//      }
//      cout << category << ": |DeltaEtaJ1J2|: "<< getVariableValue("absDeltaEtaJ1J2") << endl;
//      cout << category << ": J1J2Mass: "<< getVariableValue("J1J2Mass") << endl;
//      cout << category << ": ------------ END -------------" << endl;
//    }
//    
//    // PassJ1J2Mass1500FailTrigger event printout
//    category = "PassJ1J2Mass1500FailTrigger";
//    if( passedAllPreviousCuts("J1J2Mass1500") && passedCut("J1J2Mass1500") && deltaEtaJetFound==0 )
//    {
//      cout << category << ": ----------- START ------------" << endl;
//      cout << category << ": Run, lumi, event: "<< iEvent.id().run() << ", "
//                                                << iEvent.luminosityBlock() << ", "
//                                                << iEvent.id().event() << endl;
//      // loop over ak7 CaloJets
//      for (size_t i=0; i<CaloJetPt->size(); ++i)
//      {
//        cout << category << ": CaloJet " << i << " Pt, PtRaw, eta, phi: " << CaloJetPt->at(i) << ", "
//                                                                          << CaloJetPtRaw->at(i) << ", "
//                                                                          << CaloJetEta->at(i) << ", "
//                                                                          << CaloJetPhi->at(i) << endl;
//      }
//      // loop over ucorrected ak7 CaloJets with |eta|<2.7
//      for (size_t i=0; i<v_idx_jet_uncorr.size(); ++i)
//      {
//        cout << category << ": Uncorrected CaloJet " << i << " Pt, PtRaw, eta, phi: " << CaloJetPt->at(v_idx_jet_uncorr[i]) << ", "
//                                                                                      << CaloJetPtRaw->at(v_idx_jet_uncorr[i]) << ", "
//                                                                                      << CaloJetEta->at(v_idx_jet_uncorr[i]) << ", "
//                                                                                      << CaloJetPhi->at(v_idx_jet_uncorr[i]) << endl;
//      }
//      // loop over ak7 CaloJets that pass loose jet ID requirements
//      for (size_t i=0; i<v_idx_jet_looseID.size(); ++i)
//      {
//        cout << category << ": PassLooseID CaloJet "<< i << " Pt, PtRaw, eta, phi: " << CaloJetPt->at(v_idx_jet_looseID[i]) << ", "
//                                                                                     << CaloJetPtRaw->at(v_idx_jet_looseID[i]) << ", "
//                                                                                     << CaloJetEta->at(v_idx_jet_looseID[i]) << ", "
//                                                                                     << CaloJetPhi->at(v_idx_jet_looseID[i]) << endl;
//      }
//      // loop over ak7 PFJets
//      for (size_t i=0; i<PFJetPt->size(); ++i)
//      {
//        cout << category << ": PFJet " << i << " Pt, PtRaw, eta, phi: " << PFJetPt->at(i) << ", "
//                                                                        << PFJetPtRaw->at(i) << ", "
//                                                                        << PFJetEta->at(i) << ", "
//                                                                        << PFJetPhi->at(i) << endl;
//      }
//      // loop over ak7 PFJets that pass loose jet ID requirements
//      for (size_t i=0; i<v_idx_pfjet_looseID.size(); ++i)
//      {
//        cout << category << ": PassLooseID PFJet "<< i << " Pt, PtRaw, eta, phi: " << PFJetPt->at(v_idx_pfjet_looseID[i]) << ", "
//                                                                                   << PFJetPtRaw->at(v_idx_pfjet_looseID[i]) << ", "
//                                                                                   << PFJetEta->at(v_idx_pfjet_looseID[i]) << ", "
//                                                                                   << PFJetPhi->at(v_idx_pfjet_looseID[i]) << endl;
//      }
//      cout << category << ": |DeltaEtaJ1J2|: "<< getVariableValue("absDeltaEtaJ1J2") << endl;
//      cout << category << ": J1J2Mass: "<< getVariableValue("J1J2Mass") << endl;
//      cout << category << ": ------------ END -------------" << endl;
//    }

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

DEFINE_FWK_MODULE(MyAnalyzer);
