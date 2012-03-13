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
// $Id: MyAnalyzer_TriggerEfficiency_DijetBBTag_2011.cc,v 1.2 2011/11/17 02:12:50 ferencek Exp $
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
   CreateUserTH1D("h1_DijetMass_denom_HT600", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_num_HT600", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_denom_HT650", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_num_HT650", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_denom_HT750", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_num_HT750", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_denom_FatJet850", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_num_FatJet850", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   
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

   int btagger = int(getPreCutValue1("btagger"));
   
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
   edm::Handle<vector<int> > PFJetPassJetID;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:PassTightID"), PFJetPassJetID);
   edm::Handle<vector<double> > PFJetSSVHE;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:SimpleSecondaryVertexHighEffBTag"), PFJetSSVHE);
   edm::Handle<vector<double> > PFJetSSVHP;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:SimpleSecondaryVertexHighPurBTag"), PFJetSSVHP);
   edm::Handle<vector<double> > PFJetTCHE;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:TrackCountingHighEffBTag"), PFJetTCHE);

   int passEEAnomJetFilter = 1;
   if( PFJetPt->size() > 0 )
   {
     if( PFJetPt->at(0) > 15000 ) passEEAnomJetFilter = 0;
   }

   // Set the evaluation of the cuts to false and clear the variable values and filled status
   resetCuts();

   fillVariableWithValue("PassHBHENoiseFilter", ( *passHBHENoiseFilter ? 1 : 0 ) );
   fillVariableWithValue("PassBeamHaloFltTight", ( !(*passBeamHaloFilterTight) ? 1 : 0 ) ); // there is a bug in the ntuple maker so have to take the negative of the stored flag
   fillVariableWithValue("PassTrackingFailure", ( *passTrackingFailure ? 1 : 0 ) );
   fillVariableWithValue("PassEcalMskCellDRFlt", ( *passEcalMaskedCellDRFilter ? 1 : 0 ) );
   fillVariableWithValue("PassCaloBndDRFlt", ( *passCaloBoundaryDRFilter ? 1 : 0 ) );
   fillVariableWithValue("PassEEAnomJetFilter", passEEAnomJetFilter );

   fillVariableWithValue( "nJets", PFJetPt->size() );
   
   if( PFJetPt->size() >= 1 )
   {
       fillVariableWithValue( "passJetIdJ1", ( PFJetPassJetID->at(0) ? 1 : 0 ) );
       fillVariableWithValue( "absEtaJ1", fabs( PFJetEta->at(0) ) );
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

       TVector2 v2_j1;
       TVector2 v2_j2;
       v2_j1.SetMagPhi( 1., PFJetPhi->at(0) );
       v2_j2.SetMagPhi( 1., PFJetPhi->at(1) );
       fillVariableWithValue( "absDeltaPhiJ1J2", fabs( v2_j1.DeltaPhi(v2_j2) ) );

       fillVariableWithValue( "DijetMass", v_j1j2.M() );

       int nBTaggedJets = 0;
       for(int i=0; i<2; ++i)
       {
         if(btagger==1)
         {
           if( PFJetSSVHE->at(i) > getPreCutValue1("SSVHEM_WP") ) ++nBTaggedJets;
         }
         else if(btagger==2)
         {
           if( PFJetSSVHP->at(i) > getPreCutValue1("SSVHPT_WP") ) ++nBTaggedJets;
         }
         else
         {
           if( PFJetTCHE->at(i) > getPreCutValue1("TCHEM_WP") ) ++nBTaggedJets;
         }
       }

       fillVariableWithValue("nJets_btag", nBTaggedJets);
   }

   // Evaluate cuts (but do not apply them)
   evaluateCuts();

   
   // HLT_HT600 (wrt HLT_HT400)
   string trigger_pattern_num = "HLT_HT600_v*";
   string trigger_pattern_denom = "HLT_HT400_v*";
   
   vector<string> matched_num = hltConfig.matched(hltConfig.triggerNames(),trigger_pattern_num);
   vector<string> matched_denom = hltConfig.matched(hltConfig.triggerNames(),trigger_pattern_denom);

   if( matched_num.size()==1 && matched_denom.size()==1 )
   {
     int prescale = -1;
     if(hltConfig.prescaleSet(iEvent, iSetup)<0) edm::LogError("RootTupleMakerV2_TriggerError") << "Error! The prescale set index number could not be obtained";
     else prescale = hltConfig.prescaleValue(iEvent, iSetup, matched_num.front());
     
     if( prescale == 1 )
     {
//        cout<<matched_num.front()<<endl;
//        cout<<matched_denom.front()<<endl;
       int idx_num = hltConfig.triggerIndex(matched_num.front());
       int idx_denom = hltConfig.triggerIndex(matched_denom.front());

       if( passedAllPreviousCuts("DijetMass") && triggerResults->accept( idx_denom ) ) FillUserTH1D("h1_DijetMass_denom_HT600", getVariableValue("DijetMass") );
       if( passedAllPreviousCuts("DijetMass") && triggerResults->accept( idx_denom ) && triggerResults->accept( idx_num ) ) FillUserTH1D("h1_DijetMass_num_HT600", getVariableValue("DijetMass") );
     }
   }
   else if( matched_num.size()>1 || matched_denom.size()>1 )
   {
     edm::LogWarning("MyAnalyzer::filter") << "More than one matching trigger found for "<<trigger_pattern_num<<" or "<<trigger_pattern_denom;
     cout<<"Matched num:";
     for(vector<string>::const_iterator it = matched_num.begin(); it != matched_num.end(); ++it ) cout<<" "<<(*it);
     cout<<endl;
     cout<<"Matched denom:";
     for(vector<string>::const_iterator it = matched_denom.begin(); it != matched_denom.end(); ++it ) cout<<" "<<(*it);
     cout<<endl;
   }

   // HLT_HT650 (wrt HLT_HT400)
   trigger_pattern_num = "HLT_HT650_v*";
   trigger_pattern_denom = "HLT_HT400_v*";

   matched_num = hltConfig.matched(hltConfig.triggerNames(),trigger_pattern_num);
   matched_denom = hltConfig.matched(hltConfig.triggerNames(),trigger_pattern_denom);

   if( matched_num.size()==1 && matched_denom.size()==1 )
   {
     int prescale = -1;
     if(hltConfig.prescaleSet(iEvent, iSetup)<0) edm::LogError("RootTupleMakerV2_TriggerError") << "Error! The prescale set index number could not be obtained";
     else prescale = hltConfig.prescaleValue(iEvent, iSetup, matched_num.front());

     if( prescale == 1 )
     {
//        cout<<matched_num.front()<<endl;
//        cout<<matched_denom.front()<<endl;
       int idx_num = hltConfig.triggerIndex(matched_num.front());
       int idx_denom = hltConfig.triggerIndex(matched_denom.front());

       if( passedAllPreviousCuts("DijetMass") && triggerResults->accept( idx_denom ) ) FillUserTH1D("h1_DijetMass_denom_HT650", getVariableValue("DijetMass") );
       if( passedAllPreviousCuts("DijetMass") && triggerResults->accept( idx_denom ) && triggerResults->accept( idx_num ) ) FillUserTH1D("h1_DijetMass_num_HT650", getVariableValue("DijetMass") );
     }
   }
   else if( matched_num.size()>1 || matched_denom.size()>1 )
   {
     edm::LogWarning("MyAnalyzer::filter") << "More than one matching trigger found for "<<trigger_pattern_num<<" or "<<trigger_pattern_denom;
     cout<<"Matched num:";
     for(vector<string>::const_iterator it = matched_num.begin(); it != matched_num.end(); ++it ) cout<<" "<<(*it);
     cout<<endl;
     cout<<"Matched denom:";
     for(vector<string>::const_iterator it = matched_denom.begin(); it != matched_denom.end(); ++it ) cout<<" "<<(*it);
     cout<<endl;
   }

   // HLT_HT750 (wrt HLT_HT450)
   trigger_pattern_num = "HLT_HT750_v*";
   trigger_pattern_denom = "HLT_HT450_v*";

   matched_num = hltConfig.matched(hltConfig.triggerNames(),trigger_pattern_num);
   matched_denom = hltConfig.matched(hltConfig.triggerNames(),trigger_pattern_denom);

   if( matched_num.size()==1 && matched_denom.size()==1 )
   {
     int prescale = -1;
     if(hltConfig.prescaleSet(iEvent, iSetup)<0) edm::LogError("RootTupleMakerV2_TriggerError") << "Error! The prescale set index number could not be obtained";
     else prescale = hltConfig.prescaleValue(iEvent, iSetup, matched_num.front());

     if( prescale == 1 )
     {
//        cout<<matched_num.front()<<endl;
//        cout<<matched_denom.front()<<endl;
       int idx_num = hltConfig.triggerIndex(matched_num.front());
       int idx_denom = hltConfig.triggerIndex(matched_denom.front());

       if( passedAllPreviousCuts("DijetMass") && triggerResults->accept( idx_denom ) ) FillUserTH1D("h1_DijetMass_denom_HT750", getVariableValue("DijetMass") );
       if( passedAllPreviousCuts("DijetMass") && triggerResults->accept( idx_denom ) && triggerResults->accept( idx_num ) ) FillUserTH1D("h1_DijetMass_num_HT750", getVariableValue("DijetMass") );
     }
   }
   else if( matched_num.size()>1 || matched_denom.size()>1 )
   {
     edm::LogWarning("MyAnalyzer::filter") << "More than one matching trigger found for "<<trigger_pattern_num<<" or "<<trigger_pattern_denom;
     cout<<"Matched num:";
     for(vector<string>::const_iterator it = matched_num.begin(); it != matched_num.end(); ++it ) cout<<" "<<(*it);
     cout<<endl;
     cout<<"Matched denom:";
     for(vector<string>::const_iterator it = matched_denom.begin(); it != matched_denom.end(); ++it ) cout<<" "<<(*it);
     cout<<endl;
   }

   // HLT_FatJetMass850_DR1p1_Deta2p0 (wrt HLT_HT450)
   trigger_pattern_num = "HLT_FatJetMass850_DR1p1_Deta2p0_v*";
   trigger_pattern_denom = "HLT_HT450_v*";

   matched_num = hltConfig.matched(hltConfig.triggerNames(),trigger_pattern_num);
   matched_denom = hltConfig.matched(hltConfig.triggerNames(),trigger_pattern_denom);

   if( matched_num.size()==1 && matched_denom.size()==1 )
   {
     int prescale = -1;
     if(hltConfig.prescaleSet(iEvent, iSetup)<0) edm::LogError("RootTupleMakerV2_TriggerError") << "Error! The prescale set index number could not be obtained";
     else prescale = hltConfig.prescaleValue(iEvent, iSetup, matched_num.front());

     if( prescale == 1 )
     {
//        cout<<matched_num.front()<<endl;
//        cout<<matched_denom.front()<<endl;
       int idx_num = hltConfig.triggerIndex(matched_num.front());
       int idx_denom = hltConfig.triggerIndex(matched_denom.front());

       if( passedAllPreviousCuts("DijetMass") && triggerResults->accept( idx_denom ) ) FillUserTH1D("h1_DijetMass_denom_FatJet850", getVariableValue("DijetMass") );
       if( passedAllPreviousCuts("DijetMass") && triggerResults->accept( idx_denom ) && triggerResults->accept( idx_num ) ) FillUserTH1D("h1_DijetMass_num_FatJet850", getVariableValue("DijetMass") );
     }
   }
   else if( matched_num.size()>1 || matched_denom.size()>1 )
   {
     edm::LogWarning("MyAnalyzer::filter") << "More than one matching trigger found for "<<trigger_pattern_num<<" or "<<trigger_pattern_denom;
     cout<<"Matched num:";
     for(vector<string>::const_iterator it = matched_num.begin(); it != matched_num.end(); ++it ) cout<<" "<<(*it);
     cout<<endl;
     cout<<"Matched denom:";
     for(vector<string>::const_iterator it = matched_denom.begin(); it != matched_denom.end(); ++it ) cout<<" "<<(*it);
     cout<<endl;
   }
   
   // select only those events that pass the full selection
//    ret = passedCut("all");

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
