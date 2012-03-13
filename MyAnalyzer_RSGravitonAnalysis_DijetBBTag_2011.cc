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
// $Id: MyAnalyzer_RSGravitonAnalysis_DijetBBTag_2011.cc,v 1.5 2012/03/03 00:47:48 ferencek Exp $
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
   void init(const double SF_shift, const double TCHEL_SFb, const double TCHEL_SFl, const double TCHPT_SFb, const double TCHPT_SFl, const double SSVHPT_SFb, const double SSVHPT_SFl);
   double scaleFactor(const int partonFlavor, const int btagger);
   double scaleFactor(const int partonFlavor, const double jetPt, const double jetEta, const int btagger);
   double scaleFactorB_TCHEL(const double jetPt, const double jetEta);
   double scaleFactorC_TCHEL(const double jetPt, const double jetEta);
   double scaleFactorUDSG_TCHEL(const double jetPt, const double jetEta);

 private:
   double SF_shift_;
   double TCHEL_SFb_;
   double TCHEL_SFl_;
   double TCHPT_SFb_;
   double TCHPT_SFl_;
   double SSVHPT_SFb_;
   double SSVHPT_SFl_;
   TF1 *TCHEL_SFb_0to2p4;
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
   sfCalculator.init(getPreCutValue1("SF_Shift"),getPreCutValue1("TCHEL_SFb"),getPreCutValue1("TCHEL_SFl"),getPreCutValue1("TCHPT_SFb"),getPreCutValue1("TCHPT_SFl"),getPreCutValue1("SSVHPT_SFb"),getPreCutValue1("SSVHPT_SFl"));
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
   CreateUserTH1D("h1_DijetMass_bbbar_0tag;Dijet Mass [GeV]", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_bbbar_1tag;Dijet Mass [GeV]", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_bbbar_2tag;Dijet Mass [GeV]", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_nonbbbar_0tag;Dijet Mass [GeV]", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_nonbbbar_1tag;Dijet Mass [GeV]", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_DijetMass_nonbbbar_2tag;Dijet Mass [GeV]", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   
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
   double resonanceMass = getPreCutValue1("resonanceMass");
   int btagger = int(getPreCutValue1("btagger"));
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
   edm::Handle<vector<int> > GenParticleMotherIndex;
   iEvent.getByLabel(edm::InputTag("GenParticles:MotherIndex"), GenParticleMotherIndex);

   edm::Handle<bool> passHBHENoiseFilter;
   iEvent.getByLabel(edm::InputTag("EventSelection:PassHBHENoiseFilter"), passHBHENoiseFilter);
   
   edm::Handle<vector<unsigned int> > NPU;
   iEvent.getByLabel(edm::InputTag("GenEventInfo:PileUpNumberOfInteractions"), NPU);
   edm::Handle<vector<int> > BX;
   iEvent.getByLabel(edm::InputTag("GenEventInfo:PileUpBunchCrossing"), BX);

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

   edm::Handle<vector<double> > MET;
   iEvent.getByLabel(edm::InputTag("PFMET:Mag"), MET);
   edm::Handle<vector<double> > SumET;
   iEvent.getByLabel(edm::InputTag("PFMET:SumET"), SumET);
   
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
   edm::Handle<vector<double> > PFJetJP;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:JetProbabilityBTag"), PFJetJP);
   edm::Handle<vector<double> > PFJetCSV;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:CombinedSecondaryVertexBJetTag"), PFJetCSV);
   edm::Handle<vector<int> > PFJetPartonFlavor;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:PartonFlavor"), PFJetPartonFlavor);

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
   
   int nBTaggedJets = 0;
   vector<double> scaleFactors;
   int nHeavyFlavorJets = 0;
   int nBTaggedHeavyFlavorJets = 0;

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
     }
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

   fillVariableWithValue( "nSt3_q_fromRSG", nSt3_q_fromRSG, eventWeight );
   fillVariableWithValue( "nSt3_b_fromRSG", nSt3_b_fromRSG, eventWeight );

   fillVariableWithValue("PassHBHENoiseFilter", ( *passHBHENoiseFilter ? 1 : 0 ), eventWeight );
   
   fillVariableWithValue("nGoodVertices", v_idx_goodPV.size(), eventWeight );

   fillVariableWithValue("METoSumET", MET->front()/SumET->front(), eventWeight );
   fillVariableWithValue("METoSumET_bbbar", getVariableValue("METoSumET"), eventWeight );
   
   fillVariableWithValue( "nJets", PFJetPt->size(), eventWeight );

   if( PFJetPt->size() >= 1 )
   {
     fillVariableWithValue( "passJetIdJ1", ( PFJetPassJetID->at(0) ? 1 : 0 ), eventWeight );
     fillVariableWithValue( "absEtaJ1", fabs( PFJetEta->at(0) ), eventWeight );
     fillVariableWithValue( "PhiJ1", PFJetPhi->at(0), eventWeight );
     fillVariableWithValue( "EtaJ1", PFJetEta->at(0), eventWeight );
     fillVariableWithValue( "PtJ1", PFJetPt->at(0), eventWeight );
     fillVariableWithValue( "PhiJ1_bbbar", getVariableValue("PhiJ1"), eventWeight );
     fillVariableWithValue( "EtaJ1_bbbar", getVariableValue("EtaJ1"), eventWeight );
     fillVariableWithValue( "PtJ1_bbbar", getVariableValue("PtJ1"), eventWeight );
   }
   if( PFJetPt->size() >= 2 )
   {
     fillVariableWithValue( "passJetIdJ2", ( PFJetPassJetID->at(1) ? 1 : 0 ), eventWeight );
     fillVariableWithValue( "absEtaJ2", fabs( PFJetEta->at(1) ), eventWeight );
     fillVariableWithValue( "PhiJ2", PFJetPhi->at(1), eventWeight );
     fillVariableWithValue( "EtaJ2", PFJetEta->at(1), eventWeight );
     fillVariableWithValue( "PtJ2", PFJetPt->at(1), eventWeight );
     fillVariableWithValue( "PhiJ2_bbbar", getVariableValue("PhiJ2"), eventWeight );
     fillVariableWithValue( "EtaJ2_bbbar", getVariableValue("EtaJ2"), eventWeight );
     fillVariableWithValue( "PtJ2_bbbar", getVariableValue("PtJ2"), eventWeight );
     
     TLorentzVector v_j1j2, v_j1, v_j2;
     v_j1.SetPtEtaPhiE(PFJetPt->at(0),PFJetEta->at(0),PFJetPhi->at(0),PFJetE->at(0));
     v_j2.SetPtEtaPhiE(PFJetPt->at(1),PFJetEta->at(1),PFJetPhi->at(1),PFJetE->at(1));
     // calculate |DeltaEta(j1,j2)|
     fillVariableWithValue( "absDeltaEtaJ1J2", fabs( PFJetEta->at(0) - PFJetEta->at(1) ), eventWeight );
     fillVariableWithValue( "DeltaEtaJ1J2", getVariableValue("absDeltaEtaJ1J2"), eventWeight );
     fillVariableWithValue( "DeltaEtaJ1J2_bbbar", getVariableValue("absDeltaEtaJ1J2"), eventWeight );

     fillVariableWithValue( "DeltaPhiJ1J2", fabs( v_j1.DeltaPhi(v_j2) ), eventWeight );
     fillVariableWithValue( "DeltaPhiJ1J2_bbbar", getVariableValue("DeltaPhiJ1J2"), eventWeight );
     
     // calculate M_j1j2
     v_j1j2 = v_j1 + v_j2;
     fillVariableWithValue( "DijetMass", v_j1j2.M(), eventWeight );

     fillVariableWithValue( "DijetMass_qqbar", getVariableValue("DijetMass"), eventWeight );
     fillVariableWithValue( "DijetMass_bbbar", getVariableValue("DijetMass"), eventWeight );
     fillVariableWithValue( "x_bbbar", (getVariableValue("DijetMass")/resonanceMass), eventWeight );
   }
   
   // Evaluate cuts (but do not apply them)
   evaluateCuts();

   if( passedAllPreviousCuts("DijetMass") )
   {
     if( nSt3_b_fromRSG==2 )
     {
       if( doSFReweighting )
       {
         FillUserTH1D("h1_DijetMass_bbbar_0tag", getVariableValue("DijetMass"), eventWeight*bTagEventWeight(scaleFactors,0) );
         FillUserTH1D("h1_DijetMass_bbbar_1tag", getVariableValue("DijetMass"), eventWeight*bTagEventWeight(scaleFactors,1) );
         FillUserTH1D("h1_DijetMass_bbbar_2tag", getVariableValue("DijetMass"), eventWeight*bTagEventWeight(scaleFactors,2) );
       }
       else
       {
         if( nBTaggedJets==0 ) FillUserTH1D("h1_DijetMass_bbbar_0tag", getVariableValue("DijetMass"), eventWeight );
         if( nBTaggedJets==1 ) FillUserTH1D("h1_DijetMass_bbbar_1tag", getVariableValue("DijetMass"), eventWeight );
         if( nBTaggedJets==2 ) FillUserTH1D("h1_DijetMass_bbbar_2tag", getVariableValue("DijetMass"), eventWeight );
       }
     }
     else
     {
       if( doSFReweighting )
       {
         FillUserTH1D("h1_DijetMass_nonbbbar_0tag", getVariableValue("DijetMass"), eventWeight*bTagEventWeight(scaleFactors,0) );
         FillUserTH1D("h1_DijetMass_nonbbbar_1tag", getVariableValue("DijetMass"), eventWeight*bTagEventWeight(scaleFactors,1) );
         FillUserTH1D("h1_DijetMass_nonbbbar_2tag", getVariableValue("DijetMass"), eventWeight*bTagEventWeight(scaleFactors,2) );
       }
       else
       {
         if( nBTaggedJets==0 ) FillUserTH1D("h1_DijetMass_nonbbbar_0tag", getVariableValue("DijetMass"), eventWeight );
         if( nBTaggedJets==1 ) FillUserTH1D("h1_DijetMass_nonbbbar_1tag", getVariableValue("DijetMass"), eventWeight );
         if( nBTaggedJets==2 ) FillUserTH1D("h1_DijetMass_nonbbbar_2tag", getVariableValue("DijetMass"), eventWeight );
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
  SF_shift_ = 0.;
  TCHEL_SFb_ = 1.;
  TCHEL_SFl_ = 1.;
  TCHPT_SFb_ = 1.;
  TCHPT_SFl_ = 1.;
  SSVHPT_SFb_ = 1.;
  SSVHPT_SFl_ = 1.;
  
  TCHEL_SFb_0to2p4 = new TF1("TCHEL_SFb_0to2p4","0.603913*((1.+(0.286361*x))/(1.+(0.170474*x)))", 30.,670.);

  double ptBins[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670};
  TCHEL_SFb_errors = new TH1D("TCHEL_SFb_errors", "TCHEL_SFb_errors", 14, ptBins);
  TCHEL_SFb_errors->SetBinContent(0,0.12);
  TCHEL_SFb_errors->SetBinContent(1,0.0244956);
  TCHEL_SFb_errors->SetBinContent(2,0.0237293);
  TCHEL_SFb_errors->SetBinContent(3,0.0180131);
  TCHEL_SFb_errors->SetBinContent(4,0.0182411);
  TCHEL_SFb_errors->SetBinContent(5,0.0184592);
  TCHEL_SFb_errors->SetBinContent(6,0.0106444);
  TCHEL_SFb_errors->SetBinContent(7,0.011073);
  TCHEL_SFb_errors->SetBinContent(8,0.0106296);
  TCHEL_SFb_errors->SetBinContent(9,0.0175259);
  TCHEL_SFb_errors->SetBinContent(10,0.0161566);
  TCHEL_SFb_errors->SetBinContent(11,0.0158973);
  TCHEL_SFb_errors->SetBinContent(12,0.0186782);
  TCHEL_SFb_errors->SetBinContent(13,0.0371113);
  TCHEL_SFb_errors->SetBinContent(14,0.0289788);
  TCHEL_SFb_errors->SetBinContent(15,2*0.0289788);
  
  TCHEL_SFl_0to2p4 = new TF1("TCHEL_SFl_0to2p4","(1.10649*((1+(-9.00297e-05*x))+(2.32185e-07*(x*x))))+(-4.04925e-10*(x*(x*(x/(1+(-0.00051036*x))))))", 20.,670.);
  TCHEL_SFl_0to0p5 = new TF1("TCHEL_SFl_0to0p5","(1.13615*((1+(-0.00119852*x))+(1.17888e-05*(x*x))))+(-9.8581e-08*(x*(x*(x/(1+(0.00689317*x))))))", 20.,670.);
  TCHEL_SFl_0p5to1p0 = new TF1("TCHEL_SFl_0p5to1p0","(1.13277*((1+(-0.00084146*x))+(3.80313e-06*(x*x))))+(-8.75061e-09*(x*(x*(x/(1+(0.00118695*x))))))", 20.,670.);
  TCHEL_SFl_1p0to1p5 = new TF1("TCHEL_SFl_1p0to1p5","(1.17163*((1+(-0.000828475*x))+(3.0769e-06*(x*x))))+(-4.692e-09*(x*(x*(x/(1+(0.000337759*x))))))", 20.,670.);
  TCHEL_SFl_1p5to2p4 = new TF1("TCHEL_SFl_1p5to2p4","(1.14554*((1+(-0.000128043*x))+(4.10899e-07*(x*x))))+(-2.07565e-10*(x*(x*(x/(1+(-0.00118618*x))))))", 20.,670.);

  TCHEL_SFl_0to2p4_min = new TF1("TCHEL_SFl_0to2p4_min","(1.01541*((1+(-6.04627e-05*x))+(1.38195e-07*(x*x))))+(-2.83043e-10*(x*(x*(x/(1+(-0.000633609*x))))))", 20.,670.);
  TCHEL_SFl_0to0p5_min = new TF1("TCHEL_SFl_0to0p5_min","(1.0369*((1+(-0.000945578*x))+(7.73273e-06*(x*x))))+(-4.47791e-08*(x*(x*(x/(1+(0.00499343*x))))))", 20.,670.);
  TCHEL_SFl_0p5to1p0_min = new TF1("TCHEL_SFl_0p5to1p0_min","(0.983748*((1+(7.13613e-05*x))+(-1.08648e-05*(x*x))))+(2.96162e-06*(x*(x*(x/(1+(0.282104*x))))))", 20.,670.);
  TCHEL_SFl_1p0to1p5_min = new TF1("TCHEL_SFl_1p0to1p5_min","(1.0698*((1+(-0.000731877*x))+(2.56922e-06*(x*x))))+(-3.0318e-09*(x*(x*(x/(1+(5.04118e-05*x))))))", 20.,670.);
  TCHEL_SFl_1p5to2p4_min = new TF1("TCHEL_SFl_1p5to2p4_min","(1.04766*((1+(-6.87499e-05*x))+(2.2454e-07*(x*x))))+(-1.18395e-10*(x*(x*(x/(1+(-0.00128734*x))))))", 20.,670.);

  TCHEL_SFl_0to2p4_max = new TF1("TCHEL_SFl_0to2p4_max","(1.19751*((1+(-0.000114197*x))+(3.08558e-07*(x*x))))+(-5.27598e-10*(x*(x*(x/(1+(-0.000422372*x))))))", 20.,670.);
  TCHEL_SFl_0to0p5_max = new TF1("TCHEL_SFl_0to0p5_max","(1.22179*((1+(-0.000946228*x))+(7.37821e-06*(x*x))))+(-4.8451e-08*(x*(x*(x/(1+(0.0047976*x))))))", 20.,670.);
  TCHEL_SFl_0p5to1p0_max = new TF1("TCHEL_SFl_0p5to1p0_max","(1.22714*((1+(-0.00085562*x))+(3.74425e-06*(x*x))))+(-8.91028e-09*(x*(x*(x/(1+(0.00109346*x))))))", 20.,670.);
  TCHEL_SFl_1p0to1p5_max = new TF1("TCHEL_SFl_1p0to1p5_max","(1.27351*((1+(-0.000911891*x))+(3.5465e-06*(x*x))))+(-6.69625e-09*(x*(x*(x/(1+(0.000590847*x))))))", 20.,670.);
  TCHEL_SFl_1p5to2p4_max = new TF1("TCHEL_SFl_1p5to2p4_max","(1.24367*((1+(-0.000182494*x))+(5.92637e-07*(x*x))))+(-3.3745e-10*(x*(x*(x/(1+(-0.00107694*x))))))", 20.,670.);
}

// ------------ method that initializes the BTagScaleFactorCalculator class  ------------
void
BTagScaleFactorCalculator::init(const double SF_shift, const double TCHEL_SFb, const double TCHEL_SFl, const double TCHPT_SFb, const double TCHPT_SFl, const double SSVHPT_SFb, const double SSVHPT_SFl)
{
  SF_shift_ = SF_shift; 
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
    else if(btagger==3) return TCHPT_SFb_;
    else if(btagger==4) return SSVHPT_SFb_;
    else                return 1.;
  }
  else
  {
    if(btagger==0)      return TCHEL_SFl_;
    else if(btagger==3) return TCHPT_SFl_;
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
    if(btagger==0)  return scaleFactorB_TCHEL(jetPt,jetEta);
    else            return 1.;
  }
  else if( partonFlavor==4 )
  {
    if(btagger==0)  return scaleFactorC_TCHEL(jetPt,jetEta);
    else            return 1.;
  }
  else
  {
    if(btagger==0)  return scaleFactorUDSG_TCHEL(jetPt,jetEta);
    else            return 1.;
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

  if( SF_shift_ > 0. )
    return TCHEL_SFb_0to2p4->Eval(Pt) + TCHEL_SFb_errors->GetBinContent(TCHEL_SFb_errors->GetXaxis()->FindBin(jetPt));
  else if( SF_shift_ < 0. )
    return TCHEL_SFb_0to2p4->Eval(Pt) - TCHEL_SFb_errors->GetBinContent(TCHEL_SFb_errors->GetXaxis()->FindBin(jetPt));
  else
    return TCHEL_SFb_0to2p4->Eval(Pt);

}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for c-jets and TCHEL tagger  ------------
double
BTagScaleFactorCalculator::scaleFactorC_TCHEL(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  // for scale factor extrapolation
  if(Pt<30) Pt = 30;
  if(Pt>670) Pt = 670;

  if( SF_shift_ > 0. )
    return TCHEL_SFb_0to2p4->Eval(Pt) + 2*TCHEL_SFb_errors->GetBinContent(TCHEL_SFb_errors->GetXaxis()->FindBin(jetPt));
  else if( SF_shift_ < 0. )
    return TCHEL_SFb_0to2p4->Eval(Pt) - 2*TCHEL_SFb_errors->GetBinContent(TCHEL_SFb_errors->GetXaxis()->FindBin(jetPt));
  else
    return TCHEL_SFb_0to2p4->Eval(Pt);

}

// ------------ method that returns pT- and eta-dependent b-tag efficiency scale factor for light flavor jets and TCHEL tagger ------------
double
BTagScaleFactorCalculator::scaleFactorUDSG_TCHEL(const double jetPt, const double jetEta)
{
  double SF = 1.;
  double Pt = jetPt;
  double eta = fabs(jetEta);
  // for scale factor extrapolation
  if(Pt<20) Pt = 20;

  if( SF_shift_ > 0. )
  {
    if( Pt>670 )
      SF = (TCHEL_SFl_0to2p4->Eval(670) + 2*(TCHEL_SFl_0to2p4_max->Eval(670) - TCHEL_SFl_0to2p4->Eval(670)));
    else
    {
      if(eta<0.5)
        SF = TCHEL_SFl_0to0p5_max->Eval(Pt);
      else if(eta>=0.5 && eta<1.)
        SF = TCHEL_SFl_0p5to1p0_max->Eval(Pt);
      else if(eta>=1. && eta<1.5)
        SF = TCHEL_SFl_1p0to1p5_max->Eval(Pt);
      else
        SF = TCHEL_SFl_1p5to2p4_max->Eval(Pt);
    }
  }
  else if( SF_shift_ < 0. )
  {
    if( Pt>670 )
      SF = (TCHEL_SFl_0to2p4->Eval(670) + 2*(TCHEL_SFl_0to2p4_min->Eval(670) - TCHEL_SFl_0to2p4->Eval(670)));
    else
    {
      if(eta<0.5)
        SF = TCHEL_SFl_0to0p5_min->Eval(Pt);
      else if(eta>=0.5 && eta<1.)
        SF = TCHEL_SFl_0p5to1p0_min->Eval(Pt);
      else if(eta>=1. && eta<1.5)
        SF = TCHEL_SFl_1p0to1p5_min->Eval(Pt);
      else
        SF = TCHEL_SFl_1p5to2p4_min->Eval(Pt);
    }
  }
  else
  {
    if( Pt>670 )
      SF = TCHEL_SFl_0to2p4->Eval(670);
    else
    {
      if(eta<0.5)
        SF = TCHEL_SFl_0to0p5->Eval(Pt);
      else if(eta>=0.5 && eta<1.)
        SF = TCHEL_SFl_0p5to1p0->Eval(Pt);
      else if(eta>=1. && eta<1.5)
        SF = TCHEL_SFl_1p0to1p5->Eval(Pt);
      else
        SF = TCHEL_SFl_1p5to2p4->Eval(Pt);
    }
  }
  
  return SF;
}


DEFINE_FWK_MODULE(MyAnalyzer);
