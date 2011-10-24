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
// $Id: MyAnalyzer_TriggerStudy.cc,v 1.5 2011/10/19 18:31:18 ferencek Exp $
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

   int useHLTJet370 = (int)getPreCutValue1("useHLTJet370");
   int useSSVHE = (int)getPreCutValue1("useSSVHE");
   
   // grab necessary objects from the event
   edm::Handle<edm::TriggerResults> triggerResults;
   iEvent.getByLabel(hltInputTag, triggerResults);

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
   edm::Handle<vector<double> > PFJetSSVHE;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:SimpleSecondaryVertexHighEffBTag"), PFJetSSVHE);
   edm::Handle<vector<double> > PFJetTCHE;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:TrackCountingHighEffBTag"), PFJetTCHE);

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

   // check trigger (only in data)
   int jetTrigger_fired = 0;
   
   string jet_name = jet300_name;
   if(useHLTJet370) jet_name = jet370_name;

   if( iEvent.isRealData() )
   {
     size_t index = hltConfig.triggerIndex(jet_name);
     if( index < triggerResults->size() )
     {
       if( triggerResults->accept( index ) ) jetTrigger_fired = 1;
     }
     else
     {
       edm::LogWarning("MyAnalyzer::filter") << "Requested HLT path \"" << jet_name << "\" does not exist";
     }
   }
   else jetTrigger_fired = 1;

   // loop over PFJets
   vector<int> v_idx_pfjet_looseID;
   for(size_t i=0; i<PFJetPt->size(); i++)
   {
       // select PFJets that pass loose JetID
       if( !PFJetPassLooseID->at(i) ) continue;
       v_idx_pfjet_looseID.push_back(i);
   }

   int passEEAnomJetFilter = 1;
   if( v_idx_pfjet_looseID.size() > 0 )
   {
     if( PFJetPt->at(v_idx_pfjet_looseID[0]) > 15000 ) passEEAnomJetFilter = 0;
   }

   // Set the evaluation of the cuts to false and clear the variable values and filled status
   resetCuts();

   fillVariableWithValue("PassHLT", jetTrigger_fired );
   fillVariableWithValue("PassHBHENoiseFilter", ( *passHBHENoiseFilter ? 1 : 0 ) );
   fillVariableWithValue("PassBeamHaloFltTight", ( !(*passBeamHaloFilterTight) ? 1 : 0 ) ); // there is a bug in the ntuple maker so have to take the negative of the stored flag
   fillVariableWithValue("PassTrackingFailure", ( *passTrackingFailure ? 1 : 0 ) );
   fillVariableWithValue("PassEcalMskCellDRFlt", ( *passEcalMaskedCellDRFilter ? 1 : 0 ) );
   fillVariableWithValue("PassCaloBndDRFlt", ( *passCaloBoundaryDRFilter ? 1 : 0 ) );
   fillVariableWithValue("PassEEAnomJetFilter", passEEAnomJetFilter );

   fillVariableWithValue("nJet_all", PFJetPt->size());
   fillVariableWithValue("nJet_looseID", v_idx_pfjet_looseID.size());
   
   if( v_idx_pfjet_looseID.size() >= 1 )
   {
       fillVariableWithValue( "absEtaJ1", fabs( PFJetEta->at(v_idx_pfjet_looseID[0]) ) );
   }
   if( v_idx_pfjet_looseID.size() >= 2 )
   {
       fillVariableWithValue( "absEtaJ2", fabs( PFJetEta->at(v_idx_pfjet_looseID[1]) ) );
       
       TLorentzVector v_j1j2, v_j1, v_j2;
       v_j1.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_looseID[0]),PFJetEta->at(v_idx_pfjet_looseID[0]),PFJetPhi->at(v_idx_pfjet_looseID[0]),PFJetE->at(v_idx_pfjet_looseID[0]));
       v_j2.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_looseID[1]),PFJetEta->at(v_idx_pfjet_looseID[1]),PFJetPhi->at(v_idx_pfjet_looseID[1]),PFJetE->at(v_idx_pfjet_looseID[1]));
       // calculate |DeltaEta(j1,j2)|
       fillVariableWithValue( "absDeltaEtaJ1J2", fabs( PFJetEta->at(v_idx_pfjet_looseID[0]) - PFJetEta->at(v_idx_pfjet_looseID[1]) ) );
       // calculate M_j1j2
       v_j1j2 = v_j1 + v_j2;

       TVector2 v2_j1;
       TVector2 v2_j2;
       v2_j1.SetMagPhi( 1., PFJetPhi->at(v_idx_pfjet_looseID[0]) );
       v2_j2.SetMagPhi( 1., PFJetPhi->at(v_idx_pfjet_looseID[1]) );
       fillVariableWithValue( "absDeltaPhiJ1J2", fabs( v2_j1.DeltaPhi(v2_j2) ) );

       fillVariableWithValue( "DijetMass", v_j1j2.M() );
       fillVariableWithValue( "DijetMass_HLTID", getVariableValue("DijetMass") );
       fillVariableWithValue( "DijetMass_pretag", getVariableValue("DijetMass") );
       fillVariableWithValue( "DijetMass500", getVariableValue("DijetMass") );
       fillVariableWithValue( "DijetMass750", getVariableValue("DijetMass") );
       fillVariableWithValue( "DijetMass1000", getVariableValue("DijetMass") );
       fillVariableWithValue( "DijetMass1250", getVariableValue("DijetMass") );
       fillVariableWithValue( "DijetMass1500", getVariableValue("DijetMass") );
       fillVariableWithValue( "DijetMass1750", getVariableValue("DijetMass") );
       fillVariableWithValue( "DijetMass2000", getVariableValue("DijetMass") );
       fillVariableWithValue( "DijetMass2250", getVariableValue("DijetMass") );
       fillVariableWithValue( "DijetMass2500", getVariableValue("DijetMass") );

       int nBtags = 0;
       for(int i=0; i<2; ++i)
       {
         if(useSSVHE)
         {
           if( PFJetSSVHE->at(v_idx_pfjet_looseID[i]) > getPreCutValue1("SSVHE_WP") ) ++nBtags;
         }
         else
         {
           if( PFJetTCHE->at(v_idx_pfjet_looseID[i]) > getPreCutValue1("TCHE_WP") ) ++nBtags;
         }
       }

       fillVariableWithValue("nJet_btag", nBtags);
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
