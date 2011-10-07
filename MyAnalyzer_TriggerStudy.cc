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
// $Id: MyAnalyzer.cc,v 1.1 2011/09/16 06:45:01 ferencek Exp $
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
   edm::Handle<vector<double> > PFJetEta;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Eta"), PFJetEta);
   edm::Handle<vector<double> > PFJetPhi;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Phi"), PFJetPhi);
   edm::Handle<vector<double> > PFJetE;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:Energy"), PFJetE);
   edm::Handle<vector<int> > PFJetPassLooseID;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:PassLooseID"), PFJetPassLooseID);

   edm::Handle<vector<double> > CaloJetPtRaw;
   iEvent.getByLabel(edm::InputTag("AK7CaloJets:PtRaw"), CaloJetPtRaw);
   edm::Handle<vector<double> > CaloJetEta;
   iEvent.getByLabel(edm::InputTag("AK7CaloJets:Eta"), CaloJetEta);

   // check trigger
   int jet370_fired = 0;
   string jet370_name = "HLT_Jet370_v6";
   unsigned int index = hltConfig.triggerIndex(jet370_name);
   if( index < triggerResults->size() ) {
     if( triggerResults->accept( index ) ) jet370_fired = 1;
   } else {
     edm::LogWarning("MyAnalyzer::filter") << "Requested HLT path \"" << jet370_name << "\" does not exist";
   }

   // create a vector of CaloJet indices (for CaloJets with |eta|<getPreCutValue1("jetFidRegion")) ordered by uncorrected jet pT
   vector<int> v_idx_calojet_uncorr;

   vector<pair<size_t, myiter> > order;
   size_t n = 0;
   for (myiter it = CaloJetPtRaw->begin(); it != CaloJetPtRaw->end(); ++it, ++n)
   {
     if( fabs(CaloJetEta->at(n)) > getPreCutValue1("jetFidRegion") ) continue;
     order.push_back(make_pair(n, it));
   }
   sort(order.begin(), order.end(), ordering());

   for(unsigned int i=0; i<order.size(); i++)
   {
     v_idx_calojet_uncorr.push_back(order[i].first); // all CaloJets with |eta|<getPreCutValue1("jetFidRegion") ordered by ucorrected jet pT
   }
   // check if there are any jets present within |DeltaEta|<getPreCutValue1("jetTriggerDeltaEtaCut") from the leading uncorrected jet
   int deltaEtaJetFound = 0;
   if(v_idx_calojet_uncorr.size() >= 2)
   {
     for(unsigned int i=1; i<v_idx_calojet_uncorr.size(); i++)
     {
       if( CaloJetPtRaw->at(v_idx_calojet_uncorr[i]) < getPreCutValue1("jetTriggerPtThresh") || CaloJetEta->at(v_idx_calojet_uncorr[i]) > getPreCutValue1("jetFidRegion") ) continue;
       if( fabs(CaloJetEta->at(v_idx_calojet_uncorr[0]) - CaloJetEta->at(v_idx_calojet_uncorr[i])) < getPreCutValue1("jetTriggerDeltaEtaCut"))
       {
         deltaEtaJetFound = 1;
         break;
       }
     }
   }
      
   // select jets that pass loose JetID
   vector<int> v_idx_jet_looseID;
   for(unsigned int i=0; i<PFJetPt->size(); i++) {
       // select jets inside a fiducial region
       if( PFJetPassLooseID->at(i) ) v_idx_jet_looseID.push_back(i);
   }

   // Set the evaluation of the cuts to false and clear the variable values and filled status
   resetCuts();

   fillVariableWithValue("HLT_Jet370_v6", jet370_fired );
   fillVariableWithValue("TriggerDeltaEtaCut", deltaEtaJetFound );

   fillVariableWithValue("nJet_all", PFJetPt->size());
   fillVariableWithValue("nJet_looseID", v_idx_jet_looseID.size());
   
   if( v_idx_jet_looseID.size() >= 1 ) {
       fillVariableWithValue( "absEtaJet1", fabs( PFJetEta->at(v_idx_jet_looseID[0]) ) );
   }
   if( v_idx_jet_looseID.size() >= 2 ) {
       fillVariableWithValue( "absEtaJet2", fabs( PFJetEta->at(v_idx_jet_looseID[1]) ) );
       
       TLorentzVector v_j1j2, v_j1, v_j2;
       v_j1.SetPtEtaPhiE(PFJetPt->at(v_idx_jet_looseID[0]),PFJetEta->at(v_idx_jet_looseID[0]),PFJetPhi->at(v_idx_jet_looseID[0]),PFJetE->at(v_idx_jet_looseID[0]));
       v_j2.SetPtEtaPhiE(PFJetPt->at(v_idx_jet_looseID[1]),PFJetEta->at(v_idx_jet_looseID[1]),PFJetPhi->at(v_idx_jet_looseID[1]),PFJetE->at(v_idx_jet_looseID[1]));
       // calculate |DeltaEta(j1,j2)|
       fillVariableWithValue( "absDeltaEtaJ1J2", fabs( PFJetEta->at(v_idx_jet_looseID[0]) - PFJetEta->at(v_idx_jet_looseID[1]) ) );
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
  if ( skimWasMade_ ) {
      edm::Handle<edm::MergeableCounter> eventCounter;

      if (iLumi.getByLabel(eventCounterInputTag_, eventCounter) && eventCounter.isValid()) {
          NEvtTotBeforeWeight_ += (double) eventCounter->value;
      } else {
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
