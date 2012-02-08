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
// $Id: MyAnalyzer_bTaggingEfficiency_DijetBBTag_2011.cc,v 1.3 2012/01/31 20:59:17 ferencek Exp $
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
//    int matchingType = int(getPreCutValue1("matchingType"));
//    double matchingRadius = getPreCutValue1("matchingRadius");
   
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

   // loop over PFJets and select PFJets that pass JetID requirements
   vector<int> v_idx_pfjet_JetID;
   for(size_t i=0; i<PFJetPt->size(); i++)
   {
       if( !PFJetPassJetID->at(i) ) continue;
       v_idx_pfjet_JetID.push_back(i);
   }

   
   // Set the evaluation of the cuts to false and clear the variable values and filled status
   resetCuts();

   int nSt3_q_fromRSG = 0, nSt3_b_fromRSG = 0;

   for(size_t i=0; i<GenParticlePt->size(); i++)
   {
     if( abs(GenParticlePdgId->at(i))==5 && GenParticleStatus->at(i)==3 && GenParticleMotherIndex->at(i)>=0 )
     {
       if( abs(GenParticlePdgId->at(GenParticleMotherIndex->at(i)))==5000039 ) ++nSt3_b_fromRSG;
     }
     if( abs(GenParticlePdgId->at(i))!=21 && GenParticleStatus->at(i)==3 && GenParticleMotherIndex->at(i)>=0 )
     {
       if( abs(GenParticlePdgId->at(GenParticleMotherIndex->at(i)))==5000039 ) ++nSt3_q_fromRSG;
     }
   }

   fillVariableWithValue( "nSt3_q_fromRSG", nSt3_q_fromRSG );
   fillVariableWithValue( "nSt3_b_fromRSG", nSt3_b_fromRSG );

   
   fillVariableWithValue("nJets_all", PFJetPt->size());
   fillVariableWithValue("nJets_JetID", v_idx_pfjet_JetID.size());

   if( v_idx_pfjet_JetID.size() >= 1 )
   {
     fillVariableWithValue( "absEtaJ1", fabs( PFJetEta->at(v_idx_pfjet_JetID[0]) ));
   }
   if( v_idx_pfjet_JetID.size() >= 2 )
   {
     fillVariableWithValue( "absEtaJ2", fabs( PFJetEta->at(v_idx_pfjet_JetID[1]) ) );
     
     TLorentzVector v_j1j2, v_j1, v_j2;
     v_j1.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_JetID[0]),PFJetEta->at(v_idx_pfjet_JetID[0]),PFJetPhi->at(v_idx_pfjet_JetID[0]),PFJetE->at(v_idx_pfjet_JetID[0]));
     v_j2.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_JetID[1]),PFJetEta->at(v_idx_pfjet_JetID[1]),PFJetPhi->at(v_idx_pfjet_JetID[1]),PFJetE->at(v_idx_pfjet_JetID[1]));
     // calculate |DeltaEta(j1,j2)|
     fillVariableWithValue( "absDeltaEtaJ1J2", fabs( PFJetEta->at(v_idx_pfjet_JetID[0]) - PFJetEta->at(v_idx_pfjet_JetID[1]) ) );
     // calculate M_j1j2
     v_j1j2 = v_j1 + v_j2;
     fillVariableWithValue( "DijetMass", v_j1j2.M() );

     fillVariableWithValue( "DijetMass_qqbar", getVariableValue("DijetMass") );
     fillVariableWithValue( "DijetMass_bbbar", getVariableValue("DijetMass") );
   }
   
   // Evaluate cuts (but do not apply them)
   evaluateCuts();
   

//    int nBTaggedJets = 0;
//    int nHeavyFlavorJets = 0;
//    int nBTaggedHeavyFlavorJets = 0;
//    
//    if( v_idx_pfjet_JetID.size() >= 2 )
//    {
//      // jet and GenParticle 4-vectors
//      TLorentzVector v_j, v_gp;
// 
//      // loop over two leading jets
//      for(size_t i=0; i<2; ++i)
//      {
//        bool isHeavyFlavor = false;
// 
//        // set jet 4-vector
//        v_j.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_JetID[i]),PFJetEta->at(v_idx_pfjet_JetID[i]),PFJetPhi->at(v_idx_pfjet_JetID[i]),PFJetE->at(v_idx_pfjet_JetID[i]));
// 
//        if( !iEvent.isRealData() )
//        {
//          if( matchingType==0 && abs(PFJetPartonFlavor->at(v_idx_pfjet_JetID[i]))==5 )
//          {
//            ++nHeavyFlavorJets;
//            isHeavyFlavor = true;
//          }
//          else if( matchingType!=0 )
//          {
//            double minDeltaR = 999.;
// 
//            // loop over GenParticles
//            for(size_t j=0; j<GenParticlePt->size(); ++j)
//            {
//              int pdgID = abs(GenParticlePdgId->at(j));
// 
//              if( pdgID==511 || pdgID==521 || pdgID==531 || pdgID==541 || pdgID==5122 || pdgID==5132 || pdgID==5232 || pdgID==5332
//                  || pdgID==411 || pdgID==421 || pdgID==431 || pdgID==4122 || pdgID==4132 || pdgID==4232 || pdgID==4332 )
//              {
//                // set GenParticle 4-vector
//                v_gp.SetPtEtaPhiE(GenParticlePt->at(j),GenParticleEta->at(j),GenParticlePhi->at(j),GenParticleE->at(j));
//                double deltaR = v_j.DeltaR(v_gp);
// 
//                if( deltaR < minDeltaR ) minDeltaR = deltaR;
//              }
//            }
// 
//            if( minDeltaR < matchingRadius )
//            {
//              ++nHeavyFlavorJets;
//              isHeavyFlavor = true;
//            }
//          }
//        }
// 
//        if( (btagger==0 && PFJetTCHE->at(v_idx_pfjet_JetID[i])>getPreCutValue1("TCHEM_WP")) ||
//            (btagger==1 && PFJetSSVHE->at(v_idx_pfjet_JetID[i])>getPreCutValue1("SSVHEM_WP")) ||
//            (btagger==2 && PFJetTCHP->at(v_idx_pfjet_JetID[i])>getPreCutValue1("TCHPT_WP")) ||
//            (btagger==3 && PFJetSSVHP->at(v_idx_pfjet_JetID[i])>getPreCutValue1("SSVHPT_WP")) )
//        {
//          ++nBTaggedJets;
//          if( isHeavyFlavor ) ++nBTaggedHeavyFlavorJets;
//        }
//      }
//    }
   
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
