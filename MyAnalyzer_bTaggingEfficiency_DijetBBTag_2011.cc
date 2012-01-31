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
// $Id: MyAnalyzer_bTaggingEfficiency_DijetBBTag_2011.cc,v 1.2 2012/01/26 19:30:06 ferencek Exp $
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
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_HF;p_{T,1} [GeV];#eta_{1}", 6000, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_TCHEM;p_{T,1} [GeV];#eta_{1}", 6000, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_SSVHEM;p_{T,1} [GeV];#eta_{1}", 6000, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_TCHPT;p_{T,1} [GeV];#eta_{1}", 6000, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ1_vs_PtJ1_SSVHPT;p_{T,1} [GeV];#eta_{1}", 6000, 0, 6000, 100, -5, 5);

   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_HF;p_{T,2} [GeV];#eta_{2}", 6000, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_TCHEM;p_{T,2} [GeV];#eta_{2}", 6000, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_SSVHEM;p_{T,2} [GeV];#eta_{2}", 6000, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_TCHPT;p_{T,2} [GeV];#eta_{2}", 6000, 0, 6000, 100, -5, 5);
   CreateUserTH2D("h2_EtaJ2_vs_PtJ2_SSVHPT;p_{T,2} [GeV];#eta_{2}", 6000, 0, 6000, 100, -5, 5);
   
   CreateUserTH2D("h2_TCHEL_TCHPL_DijetMass1to1p5TeV_denom", 2, 0.5, 2.5, 2, 0.5, 2.5);
   CreateUserTH2D("h2_TCHEL_TCHPL_DijetMass1to1p5TeV_num", 2, 0.5, 2.5, 2, 0.5, 2.5);
   
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

//    edm::Handle<vector<double> > GenJetPt;
//    iEvent.getByLabel(edm::InputTag("AK7GenJets:Pt"), GenJetPt);
//    edm::Handle<vector<double> > GenJetEta;
//    iEvent.getByLabel(edm::InputTag("AK7GenJets:Eta"), GenJetEta);
//    edm::Handle<vector<double> > GenJetPhi;
//    iEvent.getByLabel(edm::InputTag("AK7GenJets:Phi"), GenJetPhi);
//    edm::Handle<vector<double> > GenJetE;
//    iEvent.getByLabel(edm::InputTag("AK7GenJets:Energy"), GenJetE);
   
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

   //    if( GenJetPt->size() >= 2 )
//    {
//      TLorentzVector v_j1j2, v_j1, v_j2;
//      v_j1.SetPtEtaPhiE(GenJetPt->at(0),GenJetEta->at(0),GenJetPhi->at(0),GenJetE->at(0));
//      v_j2.SetPtEtaPhiE(GenJetPt->at(1),GenJetEta->at(1),GenJetPhi->at(1),GenJetE->at(1));
//
//      // calculate M_j1j2
//      v_j1j2 = v_j1 + v_j2;
//      fillVariableWithValue( "DijetMassGen", v_j1j2.M() );
//    }

   vector<TLorentzVector> bSt3_vectors;
   vector<TLorentzVector> bSt2_vectors;

   for(size_t i=0; i<GenParticlePt->size(); i++)
   {
     TLorentzVector v_gen;
     if( abs(GenParticlePdgId->at(i))==5 && GenParticleStatus->at(i)==3 )
     {
       v_gen.SetPtEtaPhiE(GenParticlePt->at(i),GenParticleEta->at(i),GenParticlePhi->at(i),GenParticleE->at(i));
       bSt3_vectors.push_back(v_gen);
       //cout << "PdgId=" << GenParticlePdgId->at(i) << " Status=" << GenParticleStatus->at(i) << endl;
     }
     if( abs(GenParticlePdgId->at(i))==5 && GenParticleStatus->at(i)==2 )
     {
       v_gen.SetPtEtaPhiE(GenParticlePt->at(i),GenParticleEta->at(i),GenParticlePhi->at(i),GenParticleE->at(i));
       bSt2_vectors.push_back(v_gen);
       //cout << "PdgId=" << GenParticlePdgId->at(i) << " Status=" << GenParticleStatus->at(i) << endl;
     }
   }

   fillVariableWithValue( "nStatus3_bQuarks", bSt3_vectors.size() );
   fillVariableWithValue( "nStatus2_bQuarks", bSt2_vectors.size() );

//    if( bSt3_vectors.size()==2 )
//    {
//      TLorentzVector res_vector = bSt3_vectors.at(0) + bSt3_vectors.at(1);
//      fillVariableWithValue( "DijetMassGen", res_vector.M() );
//    }
   

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
   }
   
   // Evaluate cuts (but do not apply them)
   evaluateCuts();
   

//    int nBTaggedJets = 0;
   int nHeavyFlavorJets = 0;
//    int nBTaggedHeavyFlavorJets = 0;
   
   if( v_idx_pfjet_JetID.size() >= 2 )
   {
     // jet and GenParticle 4-vectors
     TLorentzVector v_j, v_gp;

     // loop over two leading jets
     for(size_t i=0; i<2; ++i)
     {
       bool isHeavyFlavor = false;

       // set jet 4-vector
       v_j.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_JetID[i]),PFJetEta->at(v_idx_pfjet_JetID[i]),PFJetPhi->at(v_idx_pfjet_JetID[i]),PFJetE->at(v_idx_pfjet_JetID[i]));

       if( !iEvent.isRealData() )
       {
         if( matchingType==0 && abs(PFJetPartonFlavor->at(v_idx_pfjet_JetID[i]))==5 )
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

//        if( (btagger==0 && PFJetTCHE->at(v_idx_pfjet_JetID[i])>getPreCutValue1("TCHEM_WP")) ||
//            (btagger==1 && PFJetSSVHE->at(v_idx_pfjet_JetID[i])>getPreCutValue1("SSVHEM_WP")) ||
//            (btagger==2 && PFJetTCHP->at(v_idx_pfjet_JetID[i])>getPreCutValue1("TCHPT_WP")) ||
//            (btagger==3 && PFJetSSVHP->at(v_idx_pfjet_JetID[i])>getPreCutValue1("SSVHPT_WP")) )
//        {
//          ++nBTaggedJets;
//          if( isHeavyFlavor ) ++nBTaggedHeavyFlavorJets;
//        }

       if( isHeavyFlavor && passedAllPreviousCuts("nJets_all") )
       {
         if(i==0)
         {
           FillUserTH2D("h2_EtaJ1_vs_PtJ1_HF", PFJetPt->at(v_idx_pfjet_JetID[i]), PFJetEta->at(v_idx_pfjet_JetID[i]));
           if( PFJetTCHE->at(v_idx_pfjet_JetID[i])>getPreCutValue1("TCHEM_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_TCHEM", PFJetPt->at(v_idx_pfjet_JetID[i]), PFJetEta->at(v_idx_pfjet_JetID[i]));
           if( PFJetSSVHE->at(v_idx_pfjet_JetID[i])>getPreCutValue1("SSVHEM_WP") ) FillUserTH2D("h2_EtaJ1_vs_PtJ1_SSVHEM", PFJetPt->at(v_idx_pfjet_JetID[i]), PFJetEta->at(v_idx_pfjet_JetID[i]));
           if( PFJetTCHP->at(v_idx_pfjet_JetID[i])>getPreCutValue1("TCHPT_WP") )   FillUserTH2D("h2_EtaJ1_vs_PtJ1_TCHPT", PFJetPt->at(v_idx_pfjet_JetID[i]), PFJetEta->at(v_idx_pfjet_JetID[i]));
           if( PFJetSSVHP->at(v_idx_pfjet_JetID[i])>getPreCutValue1("SSVHPT_WP") ) FillUserTH2D("h2_EtaJ1_vs_PtJ1_SSVHPT", PFJetPt->at(v_idx_pfjet_JetID[i]), PFJetEta->at(v_idx_pfjet_JetID[i]));
         }
         if(i==1)
         {
           FillUserTH2D("h2_EtaJ2_vs_PtJ2_HF", PFJetPt->at(v_idx_pfjet_JetID[i]), PFJetEta->at(v_idx_pfjet_JetID[i]));
           if( PFJetTCHE->at(v_idx_pfjet_JetID[i])>getPreCutValue1("TCHEM_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_TCHEM", PFJetPt->at(v_idx_pfjet_JetID[i]), PFJetEta->at(v_idx_pfjet_JetID[i]));
           if( PFJetSSVHE->at(v_idx_pfjet_JetID[i])>getPreCutValue1("SSVHEM_WP") ) FillUserTH2D("h2_EtaJ2_vs_PtJ2_SSVHEM", PFJetPt->at(v_idx_pfjet_JetID[i]), PFJetEta->at(v_idx_pfjet_JetID[i]));
           if( PFJetTCHP->at(v_idx_pfjet_JetID[i])>getPreCutValue1("TCHPT_WP") )   FillUserTH2D("h2_EtaJ2_vs_PtJ2_TCHPT", PFJetPt->at(v_idx_pfjet_JetID[i]), PFJetEta->at(v_idx_pfjet_JetID[i]));
           if( PFJetSSVHP->at(v_idx_pfjet_JetID[i])>getPreCutValue1("SSVHPT_WP") ) FillUserTH2D("h2_EtaJ2_vs_PtJ2_SSVHPT", PFJetPt->at(v_idx_pfjet_JetID[i]), PFJetEta->at(v_idx_pfjet_JetID[i]));
         }
       }
     }
   }

   if(passedAllPreviousCuts("DijetMass") && getVariableValue("DijetMass")>1000 && getVariableValue("DijetMass")<1500 && nHeavyFlavorJets==2)
   {
       FillUserTH2D("h2_TCHEL_TCHPL_DijetMass1to1p5TeV_denom", 1, 1);
       FillUserTH2D("h2_TCHEL_TCHPL_DijetMass1to1p5TeV_denom", 1, 2);
       FillUserTH2D("h2_TCHEL_TCHPL_DijetMass1to1p5TeV_denom", 2, 1);
       FillUserTH2D("h2_TCHEL_TCHPL_DijetMass1to1p5TeV_denom", 2, 2);

       if( PFJetTCHE->at(v_idx_pfjet_JetID[0])>getPreCutValue1("TCHEL_WP") && PFJetTCHE->at(v_idx_pfjet_JetID[1])>getPreCutValue1("TCHEL_WP") )
         FillUserTH2D("h2_TCHEL_TCHPL_DijetMass1to1p5TeV_num", 1, 1);
       if( ( PFJetTCHE->at(v_idx_pfjet_JetID[0])>getPreCutValue1("TCHEL_WP") && PFJetTCHP->at(v_idx_pfjet_JetID[1])>getPreCutValue1("TCHPL_WP") ) ||
           ( PFJetTCHP->at(v_idx_pfjet_JetID[0])>getPreCutValue1("TCHPL_WP") && PFJetTCHE->at(v_idx_pfjet_JetID[1])>getPreCutValue1("TCHEL_WP") ) )
       {
         FillUserTH2D("h2_TCHEL_TCHPL_DijetMass1to1p5TeV_num", 1, 2);
         FillUserTH2D("h2_TCHEL_TCHPL_DijetMass1to1p5TeV_num", 2, 1);
       }
       if( PFJetTCHP->at(v_idx_pfjet_JetID[0])>getPreCutValue1("TCHPL_WP") && PFJetTCHP->at(v_idx_pfjet_JetID[1])>getPreCutValue1("TCHPL_WP") )
         FillUserTH2D("h2_TCHEL_TCHPL_DijetMass1to1p5TeV_num", 2, 2);
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
