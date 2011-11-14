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
// $Id: MyAnalyzer_MainAnalysis_DijetBBTag_2011.cc,v 1.2 2011/11/10 02:22:32 ferencek Exp $
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

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

// BaseClass
#include "MyAnalysis/MyAnalyzer/interface/BaseClass.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT
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
   double scaleFactor(const int partonFlavor, const double jetPt, const double jetEta, const int btagger);
   double scaleFactorBC_SSVHEM(const double jetPt, const double jetEta);
   double scaleFactorBC_SSVHPT(const double jetPt, const double jetEta);
   double scaleFactorBC_TCHEM(const double jetPt, const double jetEta);
   double scaleFactorUDSG_SSVHEM(const double jetPt, const double jetEta);
   double scaleFactorUDSG_SSVHPT(const double jetPt, const double jetEta);
   double scaleFactorUDSG_TCHEM(const double jetPt, const double jetEta);

 private:
   TH1D *h1_SF_BC_TCHEM_lt1p2;
   TH1D *h1_SF_BC_TCHEM_1p2to2p4;
   TH1D *h1_SF_UDSG_TCHEM_lt0p8;
   TH1D *h1_SF_UDSG_TCHEM_0p8to1p6;
   TH1D *h1_SF_UDSG_TCHEM_1p6to2p4;
   TH1D *h1_SF_BC_SSVHEM_lt1p2;
   TH1D *h1_SF_BC_SSVHEM_1p2to2p4;
   TH1D *h1_SF_UDSG_SSVHEM_lt0p8;
   TH1D *h1_SF_UDSG_SSVHEM_0p8to1p6;
   TH1D *h1_SF_UDSG_SSVHEM_1p6to2p4;
   TH1D *h1_SF_BC_SSVHPT_lt1p2;
   TH1D *h1_SF_BC_SSVHPT_1p2to2p4;
   TH1D *h1_SF_UDSG_SSVHPT_lt2p4;
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
      HLTConfigProvider hltConfig;
      edm::InputTag     hltInputTag;
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
// Summer11 PU_S4 distribution obtained by only looking at the in-time crossing.  This is the "spike+smear" distribution, RECOMMENDED FOR REWEIGHTING.
double PileUpDistMC_ObservedBX0_d[35] = {1.45346E-01,6.42802E-02,6.95255E-02,6.96747E-02,6.92955E-02,6.84997E-02,6.69528E-02,6.45515E-02,6.09865E-02,5.63323E-02,5.07322E-02,4.44681E-02,3.79205E-02,3.15131E-02,2.54220E-02,2.00184E-02,1.53776E-02,1.15387E-02,8.47608E-03,6.08715E-03,4.28255E-03,2.97185E-03,2.01918E-03,1.34490E-03,8.81587E-04,5.69954E-04,3.61493E-04,2.28692E-04,1.40791E-04,8.44606E-05,5.10204E-05,3.07802E-05,1.81401E-05,1.00201E-05,5.80004E-06};
vector<float> PileUpDistMC_ObservedBX0(PileUpDistMC_ObservedBX0_d, PileUpDistMC_ObservedBX0_d + sizeof(PileUpDistMC_ObservedBX0_d) / sizeof(double) );
// Run2011A pile-up distribution
double PileUpDistData_Observed_d[35] = {12965370.0, 55851368.0, 129329360.0, 212133600.0, 276137728.0, 303603552.0, 293257504.0, 255632864.0, 204970016.0, 153263664.0, 107935616.0, 72100608.0, 45912988.0, 27970044.0, 16342576.0, 9175983.0, 4958610.0, 2582392.75, 1297695.75, 629975.0625, 295784.25, 134469.671875, 59260.0703125, 25343.8671875, 10530.08984375, 4255.04833984375, 1673.949462890625, 641.7764892578125, 240.02249145507812, 87.650428771972656, 31.280984878540039, 10.919528007507324, 3.7314565181732178, 1.2492282390594482, 0.60236752033233643};
vector<float> PileUpDistData_Observed(PileUpDistData_Observed_d, PileUpDistData_Observed_d + sizeof(PileUpDistData_Observed_d) / sizeof(double) );

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
   CreateUserTH1D("h1_J1J2PartonFlavor", 51, -0.5, 50.5);
   CreateUserTH1D("h1_nMuons_vs_DijetMass_pretag", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   CreateUserTH1D("h1_nMuons_vs_DijetMass", getHistoNBins("DijetMass"), getHistoMin("DijetMass"), getHistoMax("DijetMass"));
   
   // initialize your variables here
   LumiWeights = edm::LumiReWeighting(PileUpDistMC_ObservedBX0, PileUpDistData_Observed);
   
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
   bool ret = false;
   // event weight (by default set to 1)
   double eventWeight = 1;

   //#####################################################################################
   //########################### User's code starts here #################################

   int doPUReweighting = int(getPreCutValue1("doPUReweighting"));
   int doSFReweighting = int(getPreCutValue1("doSFReweighting"));
   int btagger = int(getPreCutValue1("btagger"));
   
   // grab necessary objects from the event
   edm::Handle<edm::TriggerResults> triggerResults;
   iEvent.getByLabel(hltInputTag, triggerResults);

   edm::Handle<vector<double> > PVZ;
   iEvent.getByLabel(edm::InputTag("Vertices:Z"), PVZ);
   
   edm::Handle<vector<unsigned int> > NPU;
   iEvent.getByLabel(edm::InputTag("GenEventInfo:PileUpNumberOfInteractions"), NPU);
   edm::Handle<vector<int> > BX;
   iEvent.getByLabel(edm::InputTag("GenEventInfo:PileUpBunchCrossing"), BX);
   
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
   edm::Handle<vector<int> > PFJetPartonFlavor;
   iEvent.getByLabel(edm::InputTag("AK7PFJets:PartonFlavor"), PFJetPartonFlavor);
   
   edm::Handle<vector<double> > MuonPt;
   iEvent.getByLabel(edm::InputTag("Muons:Pt"), MuonPt);
   edm::Handle<vector<double> > MuonEta;
   iEvent.getByLabel(edm::InputTag("Muons:Eta"), MuonEta);
   edm::Handle<vector<double> > MuonPhi;
   iEvent.getByLabel(edm::InputTag("Muons:Phi"), MuonPhi);
   edm::Handle<vector<bool> > MuonIsGlobal;
   iEvent.getByLabel(edm::InputTag("Muons:IsGlobal"), MuonIsGlobal);
   edm::Handle<vector<int> > MuonNValidTrackerHits;
   iEvent.getByLabel(edm::InputTag("Muons:NHitsTracker"), MuonNValidTrackerHits);
   edm::Handle<vector<int> > MuonNValidPixelHits;
   iEvent.getByLabel(edm::InputTag("Muons:NHitsPixel"), MuonNValidPixelHits);
   edm::Handle<vector<int> > MuonNValidMuonHits;
   iEvent.getByLabel(edm::InputTag("Muons:NHitsMuon"), MuonNValidMuonHits);
   edm::Handle<vector<int> > MuonNMatches;
   iEvent.getByLabel(edm::InputTag("Muons:NMatchedChambers"), MuonNMatches);
   edm::Handle<vector<double> > MuonChi2;
   iEvent.getByLabel(edm::InputTag("Muons:Chi2"), MuonChi2);
   edm::Handle<vector<double> > MuonNdof;
   iEvent.getByLabel(edm::InputTag("Muons:Ndof"), MuonNdof);
   edm::Handle<vector<double> > MuonTrkChi2;
   iEvent.getByLabel(edm::InputTag("Muons:TrkChi2"), MuonTrkChi2);
   edm::Handle<vector<double> > MuonTrkNdof;
   iEvent.getByLabel(edm::InputTag("Muons:TrkNdof"), MuonTrkNdof);
   edm::Handle<vector<double> > MuonPVXYDistance;
   iEvent.getByLabel(edm::InputTag("Muons:PVXYDistance"), MuonPVXYDistance);
   edm::Handle<vector<double> > MuonPV3DDistance;
   iEvent.getByLabel(edm::InputTag("Muons:PV3DDistance"), MuonPV3DDistance);

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
   
   double pretagWeight = eventWeight;
   // in MC, there are 3 iterations so each iteration gets an equal weight
   if( !iEvent.isRealData() ) pretagWeight = eventWeight/3.;
   double tagWeight = pretagWeight;
   
   string jetTrigger_name = "HLT_Jet370_v1";
   if( iEvent.id().run() >= 163269 ) jetTrigger_name = "HLT_Jet370_v2";
   if( iEvent.id().run() >= 165088 ) jetTrigger_name = "HLT_Jet370_v3";
   if( iEvent.id().run() >= 165970 ) jetTrigger_name = "HLT_Jet370_v4";
   if( iEvent.id().run() == 166346 ) jetTrigger_name = "HLT_Jet370_v5";
   if( iEvent.id().run() >= 167078 ) jetTrigger_name = "HLT_Jet370_v6";
   if( iEvent.id().run() >= 176545 ) jetTrigger_name = "HLT_Jet370_v7";
   if( iEvent.id().run() >= 178420 ) jetTrigger_name = "HLT_Jet370_v10";

   // check trigger (only in data)
   int jetTrigger_fired = 0;

   if( iEvent.isRealData() )
   {
     size_t index = hltConfig.triggerIndex(jetTrigger_name);
     if( index < triggerResults->size() )
     {
       if( triggerResults->accept( index ) ) jetTrigger_fired = 1;
     }
     else
     {
       edm::LogWarning("MyAnalyzer::filter") << "Requested HLT path \"" << jetTrigger_name << "\" does not exist";
     }
   }
   else jetTrigger_fired = 1;

   // loop over PFJets
   vector<int> v_idx_pfjet_JetID;
   for(size_t i=0; i<PFJetPt->size(); i++)
   {
       // select PFJets that pass JetID requirments
       if( !PFJetPassJetID->at(i) ) continue;
       v_idx_pfjet_JetID.push_back(i);
   }

   // loop over muons and select muons passing tight muon ID
   vector<int> v_idx_muon_tight;
   for(size_t i=0; i<MuonPt->size(); i++)
   {
       double normChi2 = (MuonNdof->at(i) != 0 ? MuonChi2->at(i) / MuonNdof->at(i) : MuonChi2->at(i) * 1e6);
       double normTrkChi2 = (MuonTrkNdof->at(i) != 0 ? MuonTrkChi2->at(i) / MuonTrkNdof->at(i) : MuonTrkChi2->at(i) * 1e6);
       double deltaZ = sqrt( pow(MuonPV3DDistance->at(i),2) - pow(MuonPVXYDistance->at(i),2) );

       if( !MuonIsGlobal->at(i) ) continue;
       if( !(MuonPt->at(i) > 5) ) continue;
       if( !(fabs(MuonEta->at(i)) < 2.4) ) continue;
       if( !(MuonNValidTrackerHits->at(i) > 10) ) continue;
       if( !(MuonNValidPixelHits->at(i) > 1) ) continue;
       if( !(MuonNValidMuonHits->at(i) > 0) ) continue;
       if( !(MuonNMatches->at(i) > 1) ) continue;
       if( !(normChi2 < 10) ) continue;
       if( !(normTrkChi2 < 10) ) continue;
       if( !(deltaZ < 1) ) continue;
       v_idx_muon_tight.push_back(i);
   }

   int passEEAnomJetFilter = 1;
   if( v_idx_pfjet_JetID.size() > 0 )
   {
     if( PFJetPt->at(v_idx_pfjet_JetID[0]) > 15000 ) passEEAnomJetFilter = 0;
   }

   // variables to be filled
   double absDeltaEtaJ1J2 = -99.;
   double absDeltaPhiJ1J2 = -99.;
   double DijetMass = -99.;
   int nBTaggedJets = 0;
   vector<double> scaleFactors; 
   int nMuons = 0;

   if( v_idx_pfjet_JetID.size() >= 2 )
   {
     TLorentzVector v_j1j2, v_j1, v_j2;
     v_j1.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_JetID[0]),PFJetEta->at(v_idx_pfjet_JetID[0]),PFJetPhi->at(v_idx_pfjet_JetID[0]),PFJetE->at(v_idx_pfjet_JetID[0]));
     v_j2.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_JetID[1]),PFJetEta->at(v_idx_pfjet_JetID[1]),PFJetPhi->at(v_idx_pfjet_JetID[1]),PFJetE->at(v_idx_pfjet_JetID[1]));
     // calculate |DeltaEta(j1,j2)|
     absDeltaEtaJ1J2 = fabs( PFJetEta->at(v_idx_pfjet_JetID[0]) - PFJetEta->at(v_idx_pfjet_JetID[1]) );
     // calculate M_j1j2
     v_j1j2 = v_j1 + v_j2;
     DijetMass = v_j1j2.M();

     TVector2 v2_j1, v2_j2;
     v2_j1.SetMagPhi( 1., PFJetPhi->at(v_idx_pfjet_JetID[0]) );
     v2_j2.SetMagPhi( 1., PFJetPhi->at(v_idx_pfjet_JetID[1]) );
     absDeltaPhiJ1J2 = fabs( v2_j1.DeltaPhi(v2_j2) );

     // loop over two leading jets
     for(size_t i=0; i<2; ++i)
     {
       if(btagger==1)
       {
         if( PFJetSSVHE->at(v_idx_pfjet_JetID[i]) > getPreCutValue1("SSVHEM_WP") )
         {
           ++nBTaggedJets;
           // if MC, get b-tag scale factor
           if( !iEvent.isRealData() ) scaleFactors.push_back(sfCalculator.scaleFactor(PFJetPartonFlavor->at(v_idx_pfjet_JetID[i]),PFJetPt->at(v_idx_pfjet_JetID[i]),PFJetEta->at(v_idx_pfjet_JetID[i]),btagger));
         }
       }
       else if(btagger==2)
       {
         if( PFJetSSVHP->at(v_idx_pfjet_JetID[i]) > getPreCutValue1("SSVHPT_WP") )
         {
           ++nBTaggedJets;
           // if MC, get b-tag scale factor
           if( !iEvent.isRealData() ) scaleFactors.push_back(sfCalculator.scaleFactor(PFJetPartonFlavor->at(v_idx_pfjet_JetID[i]),PFJetPt->at(v_idx_pfjet_JetID[i]),PFJetEta->at(v_idx_pfjet_JetID[i]),btagger));
         }
       }
       else
       {
         if( PFJetTCHE->at(v_idx_pfjet_JetID[i]) > getPreCutValue1("TCHEM_WP") )
         {
           ++nBTaggedJets;
           // if MC, get b-tag scale factor
           if( !iEvent.isRealData() ) scaleFactors.push_back(sfCalculator.scaleFactor(PFJetPartonFlavor->at(v_idx_pfjet_JetID[i]),PFJetPt->at(v_idx_pfjet_JetID[i]),PFJetEta->at(v_idx_pfjet_JetID[i]),btagger));
         }
       }

       // jet and muon Lorentz vectors
       TLorentzVector v_j, v_m;
       v_j.SetPtEtaPhiE(PFJetPt->at(v_idx_pfjet_JetID[i]),PFJetEta->at(v_idx_pfjet_JetID[i]),PFJetPhi->at(v_idx_pfjet_JetID[i]),PFJetE->at(v_idx_pfjet_JetID[i]));

       // loop over all tight muons and find those that are inside the jet (DeltaR<0.4)
       for(size_t j=0; j<v_idx_muon_tight.size(); ++j)
       {
         v_m.SetPtEtaPhiM(MuonPt->at(v_idx_muon_tight[j]),MuonEta->at(v_idx_muon_tight[j]),MuonPhi->at(v_idx_muon_tight[j]),0);
         if( v_j.DeltaR(v_m) < 0.4 ) ++nMuons;
       }
     }
   }

   // Set the evaluation of the cuts to false and clear the variable values and filled status
   resetCuts();
   
   for( int nbtags=0; nbtags<=2; ++nbtags )
   {
     // multiple iterations done only in MC
     if( iEvent.isRealData() && nbtags>0 ) break;
     
     fillVariableWithValue("PassHLT", jetTrigger_fired, pretagWeight );
     fillVariableWithValue("PassHBHENoiseFilter", ( *passHBHENoiseFilter ? 1 : 0 ), pretagWeight );
     fillVariableWithValue("PassBeamHaloFltTight", ( !(*passBeamHaloFilterTight) ? 1 : 0 ), pretagWeight ); // there is a bug in the ntuple maker (V00-00-01 and V00-00-02) so have to take the negative of the stored flag
     fillVariableWithValue("PassTrackingFailure", ( *passTrackingFailure ? 1 : 0 ), pretagWeight );
     fillVariableWithValue("PassEcalMskCellDRFlt", ( *passEcalMaskedCellDRFilter ? 1 : 0 ), pretagWeight );
     fillVariableWithValue("PassCaloBndDRFlt", ( *passCaloBoundaryDRFilter ? 1 : 0 ), pretagWeight );
     fillVariableWithValue("PassEEAnomJetFilter", passEEAnomJetFilter, pretagWeight );

     fillVariableWithValue("nJets_all", PFJetPt->size(), pretagWeight);
     fillVariableWithValue("nJets_JetID", v_idx_pfjet_JetID.size(), pretagWeight);

     fillVariableWithValue("nGoodVertices", PVZ->size(), pretagWeight );

     if( v_idx_pfjet_JetID.size() >= 1 )
     {
       fillVariableWithValue( "absEtaJ1", fabs( PFJetEta->at(v_idx_pfjet_JetID[0]) ), pretagWeight );
       fillVariableWithValue( "PtJ1", PFJetPt->at(v_idx_pfjet_JetID[0]), pretagWeight );
     }
     if( v_idx_pfjet_JetID.size() >= 2 )
     {
       fillVariableWithValue( "absEtaJ2", fabs( PFJetEta->at(v_idx_pfjet_JetID[1]) ), pretagWeight );
       fillVariableWithValue( "PtJ2", PFJetPt->at(v_idx_pfjet_JetID[1]), pretagWeight );
       fillVariableWithValue( "absDeltaEtaJ1J2", absDeltaEtaJ1J2, pretagWeight );
       fillVariableWithValue( "DijetMass1050", DijetMass, pretagWeight );
       fillVariableWithValue( "absDeltaPhiJ1J2", absDeltaPhiJ1J2, pretagWeight );
       fillVariableWithValue( "DijetMass_pretag", DijetMass, pretagWeight );
       fillVariableWithValue( "nMuons_pretag", nMuons, pretagWeight );
       // in MC, apply the b-tag event weight
       if( !iEvent.isRealData() && doSFReweighting )
       {
         tagWeight = eventWeight*bTagEventWeight(scaleFactors,nbtags);
         fillVariableWithValue( "nJets_btag", nbtags, tagWeight );
       }
       else
       {
         fillVariableWithValue( "nJets_btag", nBTaggedJets, tagWeight );
       }
       fillVariableWithValue( "DijetMass", DijetMass, tagWeight );
       fillVariableWithValue( "nMuons", nMuons, tagWeight );
     }

     // Evaluate cuts (but do not apply them)
     evaluateCuts();

     if(passedAllPreviousCuts("DijetMass_pretag"))
     {
       FillUserTH1D("h1_J1J2PartonFlavor", abs( PFJetPartonFlavor->at(v_idx_pfjet_JetID[0]) ), pretagWeight );
       FillUserTH1D("h1_J1J2PartonFlavor", abs( PFJetPartonFlavor->at(v_idx_pfjet_JetID[1]) ), pretagWeight );
     }
     if(passedAllPreviousCuts("nMuons_pretag")) FillUserTH1D("h1_nMuons_vs_DijetMass_pretag", DijetMass, double(nMuons)*pretagWeight );
     if(passedAllPreviousCuts("nMuons")) FillUserTH1D("h1_nMuons_vs_DijetMass", DijetMass, double(nMuons)*tagWeight );
  
     // select only those events that pass the full selection
     if( passedCut("all") ) ret = true;
     
     if( nbtags<2 ) resetCuts("sameEvent");
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

// ------------ method that calculates b-tag scale factors  ------------
double
BTagScaleFactorCalculator::scaleFactor(const int partonFlavor, const double jetPt, const double jetEta, const int btagger)
{
  // if parton flavor is c or b
  if( abs(partonFlavor)==4 || abs(partonFlavor)==5 )
  {
    if(btagger==1)      return scaleFactorBC_SSVHEM(jetPt, jetEta);
    else if(btagger==2) return scaleFactorBC_SSVHPT(jetPt, jetEta);
    else                return scaleFactorBC_TCHEM(jetPt, jetEta);
  }
  // if parton flavor is different from c or b (N.B.: this also includes parton flavors different from u,d,s,g which might not be the fully correct treatment)
  else
  {
    if(btagger==1)      return scaleFactorUDSG_SSVHEM(jetPt, jetEta);
    else if(btagger==2) return scaleFactorUDSG_SSVHPT(jetPt, jetEta);
    else                return scaleFactorUDSG_TCHEM(jetPt, jetEta);
  }
}

double
BTagScaleFactorCalculator::scaleFactorBC_TCHEM(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  double eta = fabs(jetEta);
  // for scale factor extrapolation
  if(Pt<20) Pt = 20.1;
  if(Pt>=240) Pt = 239.9;
  if(eta>=2.4) eta = 2.39;

  if(eta<1.2)
  {
    int bin = h1_SF_BC_TCHEM_lt1p2->GetXaxis()->FindBin(Pt);
    return h1_SF_BC_TCHEM_lt1p2->GetBinContent(bin);
  }
  else
  {
    int bin = h1_SF_BC_TCHEM_1p2to2p4->GetXaxis()->FindBin(Pt);
    return h1_SF_BC_TCHEM_1p2to2p4->GetBinContent(bin);
  }
}

double
BTagScaleFactorCalculator::scaleFactorBC_SSVHEM(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  double eta = fabs(jetEta);
  // for scale factor extrapolation
  if(Pt<20) Pt = 20.1;
  if(Pt>=240) Pt = 239.9;
  if(eta>=2.4) eta = 2.39;

  if(eta<1.2)
  {
    int bin = h1_SF_BC_SSVHEM_lt1p2->GetXaxis()->FindBin(Pt);
    return h1_SF_BC_SSVHEM_lt1p2->GetBinContent(bin);
  }
  else
  {
    int bin = h1_SF_BC_SSVHEM_1p2to2p4->GetXaxis()->FindBin(Pt);
    return h1_SF_BC_SSVHEM_1p2to2p4->GetBinContent(bin);
  }
}

double
BTagScaleFactorCalculator::scaleFactorBC_SSVHPT(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  double eta = fabs(jetEta);
  // for scale factor extrapolation
  if(Pt<20) Pt = 20.1;
  if(Pt>=240) Pt = 239.9;
  if(eta>=2.4) eta = 2.39;

  if(eta<1.2)
  {
    int bin = h1_SF_BC_SSVHPT_lt1p2->GetXaxis()->FindBin(Pt);
    return h1_SF_BC_SSVHPT_lt1p2->GetBinContent(bin);
  }
  else
  {
    int bin = h1_SF_BC_SSVHPT_1p2to2p4->GetXaxis()->FindBin(Pt);
    return h1_SF_BC_SSVHPT_1p2to2p4->GetBinContent(bin);
  }
}

double
BTagScaleFactorCalculator::scaleFactorUDSG_TCHEM(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  double eta = fabs(jetEta);
  // for scale factor extrapolation
  if(Pt<20) Pt = 20.1;
  if(Pt>=520) Pt = 519.9;
  if(eta>=2.4) eta = 2.39;

  if(eta<0.8)
  {
    int bin = h1_SF_UDSG_TCHEM_lt0p8->GetXaxis()->FindBin(Pt);
    return h1_SF_UDSG_TCHEM_lt0p8->GetBinContent(bin);
  }
  else if(eta>=0.8 && eta<1.6)
  {
    int bin = h1_SF_UDSG_TCHEM_0p8to1p6->GetXaxis()->FindBin(Pt);
    return h1_SF_UDSG_TCHEM_0p8to1p6->GetBinContent(bin);
  }
  else
  {
    int bin = h1_SF_UDSG_TCHEM_1p6to2p4->GetXaxis()->FindBin(Pt);
    return h1_SF_UDSG_TCHEM_1p6to2p4->GetBinContent(bin);
  }
}

double
BTagScaleFactorCalculator::scaleFactorUDSG_SSVHEM(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  double eta = fabs(jetEta);
  // for scale factor extrapolation
  if(Pt<20) Pt = 20.1;
  if(Pt>=520) Pt = 519.9;
  if(eta>=2.4) eta = 2.39;

  if(eta<0.8)
  {
    int bin = h1_SF_UDSG_SSVHEM_lt0p8->GetXaxis()->FindBin(Pt);
    return h1_SF_UDSG_SSVHEM_lt0p8->GetBinContent(bin);
  }
  else if(eta>=0.8 && eta<1.6)
  {
    int bin = h1_SF_UDSG_SSVHEM_0p8to1p6->GetXaxis()->FindBin(Pt);
    return h1_SF_UDSG_SSVHEM_0p8to1p6->GetBinContent(bin);
  }
  else
  {
    int bin = h1_SF_UDSG_SSVHEM_1p6to2p4->GetXaxis()->FindBin(Pt);
    return h1_SF_UDSG_SSVHEM_1p6to2p4->GetBinContent(bin);
  }
}

double
BTagScaleFactorCalculator::scaleFactorUDSG_SSVHPT(const double jetPt, const double jetEta)
{
  double Pt = jetPt;
  // for scale factor extrapolation
  if(Pt<20) Pt = 20.1;
  if(Pt>=520) Pt = 519.9;

  int bin = h1_SF_UDSG_SSVHPT_lt2p4->GetXaxis()->FindBin(Pt);
  return h1_SF_UDSG_SSVHPT_lt2p4->GetBinContent(bin);
}

BTagScaleFactorCalculator::BTagScaleFactorCalculator()
{
   h1_SF_BC_TCHEM_lt1p2 = new TH1D("h1_SF_BC_TCHEM_lt1p2","TCHEM_BTagScaleFactors",1,20,240);
   h1_SF_BC_TCHEM_lt1p2->SetBinContent(1,0.94);
   h1_SF_BC_TCHEM_lt1p2->SetBinError(1,0.094);
   
   h1_SF_BC_TCHEM_1p2to2p4 = new TH1D("h1_SF_BC_TCHEM_1p2to2p4","TCHEM_BTagScaleFactors",1,20,240);
   h1_SF_BC_TCHEM_1p2to2p4->SetBinContent(1,0.93);
   h1_SF_BC_TCHEM_1p2to2p4->SetBinError(1,0.093);

   h1_SF_BC_SSVHEM_lt1p2 = new TH1D("h1_SF_BC_SSVHEM_lt1p2","SSVHEM_BTagScaleFactors",1,20,240);
   h1_SF_BC_SSVHEM_lt1p2->SetBinContent(1,0.95);
   h1_SF_BC_SSVHEM_lt1p2->SetBinError(1,0.095);
   
   h1_SF_BC_SSVHEM_1p2to2p4 = new TH1D("h1_SF_BC_SSVHEM_1p2to2p4","SSVHEM_BTagScaleFactors",1,20,240);
   h1_SF_BC_SSVHEM_1p2to2p4->SetBinContent(1,0.93);
   h1_SF_BC_SSVHEM_1p2to2p4->SetBinError(1,0.093);

   h1_SF_BC_SSVHPT_lt1p2 = new TH1D("h1_SF_BC_SSVHPT_lt1p2","SSVHPT_BTagScaleFactors",1,20,240);
   h1_SF_BC_SSVHPT_lt1p2->SetBinContent(1,0.89);
   h1_SF_BC_SSVHPT_lt1p2->SetBinError(1,0.089);

   h1_SF_BC_SSVHPT_1p2to2p4 = new TH1D("h1_SF_BC_SSVHPT_1p2to2p4","SSVHPT_BTagScaleFactors",1,20,240);
   h1_SF_BC_SSVHPT_1p2to2p4->SetBinContent(1,0.9);
   h1_SF_BC_SSVHPT_1p2to2p4->SetBinError(1,0.09);
   
   h1_SF_UDSG_TCHEM_lt0p8 = new TH1D("h1_SF_UDSG_TCHEM_lt0p8","TCHEM_MistagScaleFactors",50,20,520);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(1,1.27465);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(2,1.26837);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(3,1.262178);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(4,1.256079);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(5,1.250078);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(6,1.24418);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(7,1.23839);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(8,1.232713);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(9,1.227155);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(10,1.22172);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(11,1.216414);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(12,1.211241);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(13,1.206207);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(14,1.201316);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(15,1.196574);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(16,1.191987);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(17,1.187558);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(18,1.183293);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(19,1.179197);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(20,1.175275);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(21,1.171533);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(22,1.167975);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(23,1.164606);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(24,1.161432);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(25,1.158458);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(26,1.155688);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(27,1.153128);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(28,1.150782);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(29,1.148657);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(30,1.146756);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(31,1.145086);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(32,1.14365);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(33,1.142455);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(34,1.141505);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(35,1.140805);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(36,1.140361);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(37,1.140177);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(38,1.14026);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(39,1.140612);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(40,1.14124);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(41,1.14215);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(42,1.143345);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(43,1.144831);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(44,1.146613);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(45,1.148696);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(46,1.151085);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(47,1.153786);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(48,1.156803);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(49,1.160142);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinContent(50,1.163807);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(1,0.2199212);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(2,0.2144088);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(3,0.2091886);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(4,0.2042528);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(5,0.1995932);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(6,0.1952018);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(7,0.1910705);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(8,0.1871914);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(9,0.1835564);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(10,0.1801575);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(11,0.1769865);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(12,0.1740356);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(13,0.1712966);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(14,0.1687616);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(15,0.1664224);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(16,0.164271);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(17,0.1622995);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(18,0.1604997);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(19,0.1588637);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(20,0.1573833);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(21,0.1560507);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(22,0.1548576);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(23,0.1537961);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(24,0.1528582);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(25,0.1520358);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(26,0.1513209);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(27,0.1507054);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(28,0.1501813);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(29,0.1497406);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(30,0.1493752);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(31,0.1490771);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(32,0.1488383);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(33,0.1486507);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(34,0.1485063);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(35,0.148397);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(36,0.1483149);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(37,0.1482518);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(38,0.1481998);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(39,0.1481507);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(40,0.1480967);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(41,0.1480295);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(42,0.1479413);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(43,0.1478239);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(44,0.1476693);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(45,0.1474696);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(46,0.1472165);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(47,0.1469022);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(48,0.1465186);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(49,0.1460576);
   h1_SF_UDSG_TCHEM_lt0p8->SetBinError(50,0.1455112);
   
   h1_SF_UDSG_TCHEM_0p8to1p6 = new TH1D("h1_SF_UDSG_TCHEM_0p8to1p6","TCHEM_MistagScaleFactors",50,20,520);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(1,1.269763);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(2,1.255674);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(3,1.242437);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(4,1.230034);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(5,1.218451);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(6,1.207673);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(7,1.197683);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(8,1.188466);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(9,1.180006);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(10,1.172288);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(11,1.165296);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(12,1.159015);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(13,1.153429);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(14,1.148523);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(15,1.14428);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(16,1.140685);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(17,1.137724);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(18,1.135379);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(19,1.133636);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(20,1.132478);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(21,1.131891);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(22,1.131859);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(23,1.132366);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(24,1.133396);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(25,1.134935);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(26,1.136966);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(27,1.139473);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(28,1.142442);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(29,1.145857);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(30,1.149701);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(31,1.15396);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(32,1.158618);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(33,1.163659);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(34,1.169068);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(35,1.174828);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(36,1.180926);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(37,1.187344);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(38,1.194068);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(39,1.201081);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(40,1.208369);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(41,1.215915);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(42,1.223704);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(43,1.23172);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(44,1.239948);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(45,1.248373);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(46,1.256978);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(47,1.265748);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(48,1.274668);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(49,1.283721);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinContent(50,1.292893);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(1,0.1895142);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(2,0.1843697);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(3,0.1795465);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(4,0.1750363);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(5,0.1708301);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(6,0.1669195);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(7,0.1632958);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(8,0.1599502);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(9,0.1568742);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(10,0.1540591);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(11,0.1514962);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(12,0.1491768);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(13,0.1470924);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(14,0.1452343);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(15,0.1435937);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(16,0.1421621);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(17,0.1409308);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(18,0.1398911);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(19,0.1390344);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(20,0.138352);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(21,0.1378353);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(22,0.1374755);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(23,0.1372642);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(24,0.1371925);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(25,0.1372519);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(26,0.1374336);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(27,0.137729);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(28,0.1381296);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(29,0.1386265);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(30,0.1392112);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(31,0.139875);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(32,0.1406093);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(33,0.1414053);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(34,0.1422544);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(35,0.1431481);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(36,0.1440775);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(37,0.1450341);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(38,0.1460092);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(39,0.1469942);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(40,0.1479803);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(41,0.148959);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(42,0.1499215);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(43,0.1508593);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(44,0.1517636);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(45,0.1526259);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(46,0.1534374);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(47,0.1541894);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(48,0.1548735);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(49,0.1554808);
   h1_SF_UDSG_TCHEM_0p8to1p6->SetBinError(50,0.1560027);
   
   h1_SF_UDSG_TCHEM_1p6to2p4 = new TH1D("h1_SF_UDSG_TCHEM_1p6to2p4","TCHEM_MistagScaleFactors",50,20,520);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(1,1.183647);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(2,1.168223);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(3,1.15477);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(4,1.143198);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(5,1.133419);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(6,1.125344);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(7,1.118883);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(8,1.113949);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(9,1.110452);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(10,1.108303);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(11,1.107414);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(12,1.107695);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(13,1.109057);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(14,1.111412);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(15,1.114671);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(16,1.118745);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(17,1.123544);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(18,1.128981);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(19,1.134966);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(20,1.141411);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(21,1.148225);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(22,1.155322);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(23,1.162611);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(24,1.170003);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(25,1.177411);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(26,1.184744);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(27,1.191914);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(28,1.198833);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(29,1.205411);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(30,1.211559);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(31,1.217189);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(32,1.222211);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(33,1.226538);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(34,1.230078);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(35,1.232745);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(36,1.234449);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(37,1.235101);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(38,1.234612);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(39,1.232894);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(40,1.229857);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(41,1.225413);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(42,1.219472);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(43,1.211946);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(44,1.202746);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(45,1.191783);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(46,1.178968);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(47,1.164212);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(48,1.147427);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(49,1.128523);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinContent(50,1.107411);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(1,0.1551869);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(2,0.1503253);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(3,0.1459676);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(4,0.1420933);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(5,0.1386818);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(6,0.1357125);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(7,0.1331648);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(8,0.1310183);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(9,0.1292522);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(10,0.1278459);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(11,0.126779);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(12,0.1260308);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(13,0.1255808);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(14,0.1254084);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(15,0.1254929);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(16,0.1258138);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(17,0.1263505);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(18,0.1270825);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(19,0.1279891);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(20,0.1290498);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(21,0.130244);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(22,0.131551);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(23,0.1329505);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(24,0.1344216);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(25,0.1359439);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(26,0.1374968);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(27,0.1390597);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(28,0.140612);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(29,0.1421331);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(30,0.1436025);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(31,0.1449995);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(32,0.1463036);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(33,0.1474942);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(34,0.1485508);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(35,0.1494527);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(36,0.1501793);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(37,0.1507101);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(38,0.1510245);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(39,0.1511019);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(40,0.1509217);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(41,0.1504634);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(42,0.1497064);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(43,0.14863);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(44,0.1472137);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(45,0.1454369);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(46,0.1432791);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(47,0.1407196);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(48,0.1377379);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(49,0.1343133);
   h1_SF_UDSG_TCHEM_1p6to2p4->SetBinError(50,0.1304254);

   h1_SF_UDSG_SSVHEM_lt0p8 = new TH1D("h1_SF_UDSG_SSVHEM_lt0p8","SSVHEM_MistagScaleFactors",50,20,520);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(1,0.8849058);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(2,0.8874509);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(3,0.8902115);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(4,0.8931665);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(5,0.8962955);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(6,0.8995775);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(7,0.9029919);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(8,0.9065179);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(9,0.9101347);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(10,0.9138217);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(11,0.917558);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(12,0.921323);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(13,0.9250959);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(14,0.9288558);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(15,0.9325822);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(16,0.9362542);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(17,0.9398511);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(18,0.9433522);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(19,0.9467367);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(20,0.9499838);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(21,0.953073);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(22,0.9559832);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(23,0.9586939);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(24,0.9611843);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(25,0.9634336);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(26,0.9654211);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(27,0.9671261);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(28,0.9685279);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(29,0.9696055);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(30,0.9703384);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(31,0.9707057);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(32,0.9706869);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(33,0.9702609);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(34,0.9694073);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(35,0.9681051);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(36,0.9663336);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(37,0.9640722);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(38,0.9613001);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(39,0.9579964);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(40,0.9541405);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(41,0.9497117);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(42,0.9446892);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(43,0.9390521);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(44,0.9327799);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(45,0.9258517);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(46,0.9182469);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(47,0.9099445);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(48,0.900924);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(49,0.8911646);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinContent(50,0.8806455);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(1,0.1054525);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(2,0.1043027);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(3,0.1033287);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(4,0.1025217);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(5,0.1018734);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(6,0.101375);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(7,0.1010181);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(8,0.1007939);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(9,0.1006941);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(10,0.1007099);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(11,0.1008328);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(12,0.1010543);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(13,0.1013657);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(14,0.1017585);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(15,0.102224);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(16,0.1027539);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(17,0.1033393);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(18,0.1039718);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(19,0.1046429);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(20,0.1053438);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(21,0.1060661);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(22,0.1068012);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(23,0.1075404);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(24,0.1082752);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(25,0.1089971);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(26,0.1096975);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(27,0.1103677);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(28,0.1109993);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(29,0.1115836);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(30,0.112112);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(31,0.112576);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(32,0.112967);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(33,0.1132765);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(34,0.1134958);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(35,0.1136163);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(36,0.1136296);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(37,0.113527);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(38,0.1132999);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(39,0.1129398);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(40,0.1124381);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(41,0.1117862);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(42,0.1109756);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(43,0.1099976);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(44,0.1088437);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(45,0.1075053);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(46,0.1059738);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(47,0.1042407);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(48,0.1022974);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(49,0.1001353);
   h1_SF_UDSG_SSVHEM_lt0p8->SetBinError(50,0.09774583);
   
   h1_SF_UDSG_SSVHEM_0p8to1p6 = new TH1D("h1_SF_UDSG_SSVHEM_0p8to1p6","SSVHEM_MistagScaleFactors",50,20,520);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(1,0.9400575);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(2,0.9359007);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(3,0.932184);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(4,0.9288992);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(5,0.9260378);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(6,0.9235915);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(7,0.9215519);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(8,0.9199107);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(9,0.9186596);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(10,0.9177901);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(11,0.9172938);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(12,0.9171626);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(13,0.917388);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(14,0.9179615);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(15,0.918875);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(16,0.92012);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(17,0.9216881);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(18,0.9235712);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(19,0.9257606);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(20,0.9282481);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(21,0.9310254);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(22,0.9340841);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(23,0.9374158);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(24,0.9410122);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(25,0.9448649);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(26,0.9489656);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(27,0.953306);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(28,0.9578775);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(29,0.9626719);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(30,0.9676809);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(31,0.9728961);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(32,0.9783091);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(33,0.9839116);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(34,0.9896952);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(35,0.9956515);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(36,1.001772);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(37,1.008049);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(38,1.014474);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(39,1.021037);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(40,1.027732);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(41,1.03455);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(42,1.041481);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(43,1.048519);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(44,1.055654);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(45,1.062879);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(46,1.070184);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(47,1.077561);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(48,1.085003);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(49,1.092501);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinContent(50,1.100046);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(1,0.113446);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(2,0.1108837);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(3,0.1086086);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(4,0.1066096);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(5,0.1048753);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(6,0.1033946);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(7,0.1021562);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(8,0.1011489);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(9,0.1003615);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(10,0.09978278);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(11,0.09940149);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(12,0.09920642);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(13,0.09918634);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(14,0.09933003);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(15,0.09962627);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(16,0.1000638);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(17,0.1006315);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(18,0.101318);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(19,0.1021122);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(20,0.1030028);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(21,0.1039786);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(22,0.1050284);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(23,0.106141);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(24,0.1073051);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(25,0.1085095);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(26,0.109743);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(27,0.1109944);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(28,0.1122524);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(29,0.1135058);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(30,0.1147434);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(31,0.115954);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(32,0.1171263);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(33,0.1182491);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(34,0.1193113);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(35,0.1203015);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(36,0.1212086);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(37,0.1220213);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(38,0.1227284);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(39,0.1233187);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(40,0.123781);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(41,0.124104);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(42,0.1242765);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(43,0.1242873);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(44,0.1241251);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(45,0.1237788);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(46,0.1232372);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(47,0.1224889);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(48,0.1215228);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(49,0.1203276);
   h1_SF_UDSG_SSVHEM_0p8to1p6->SetBinError(50,0.1188921);
   
   h1_SF_UDSG_SSVHEM_1p6to2p4 = new TH1D("h1_SF_UDSG_SSVHEM_1p6to2p4","SSVHEM_MistagScaleFactors",50,20,520);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(1,0.9725683);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(2,0.9583432);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(3,0.9456732);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(4,0.9344997);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(5,0.9247639);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(6,0.9164071);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(7,0.9093708);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(8,0.903596);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(9,0.8990244);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(10,0.8955969);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(11,0.8932552);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(12,0.8919403);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(13,0.8915936);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(14,0.8921565);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(15,0.8935703);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(16,0.8957762);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(17,0.8987156);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(18,0.9023298);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(19,0.9065601);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(20,0.9113477);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(21,0.9166341);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(22,0.9223606);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(23,0.9284683);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(24,0.9348987);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(25,0.9415931);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(26,0.9484927);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(27,0.9555389);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(28,0.962673);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(29,0.9698364);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(30,0.9769701);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(31,0.9840158);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(32,0.9909145);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(33,0.9976077);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(34,1.004037);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(35,1.010143);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(36,1.015867);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(37,1.021151);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(38,1.025936);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(39,1.030164);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(40,1.033774);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(41,1.03671);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(42,1.038913);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(43,1.040322);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(44,1.040881);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(45,1.040529);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(46,1.039209);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(47,1.036862);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(48,1.03343);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(49,1.028852);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinContent(50,1.023071);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(1,0.1131297);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(2,0.1103853);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(3,0.1078943);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(4,0.1056489);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(5,0.1036417);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(6,0.1018652);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(7,0.1003118);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(8,0.09897406);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(9,0.09784432);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(10,0.0969151);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(11,0.09617887);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(12,0.09562808);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(13,0.09525521);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(14,0.09505273);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(15,0.09501309);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(16,0.09512877);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(17,0.09539222);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(18,0.09579591);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(19,0.09633232);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(20,0.0969939);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(21,0.09777312);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(22,0.09866245);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(23,0.09965435);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(24,0.1007413);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(25,0.1019157);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(26,0.1031701);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(27,0.104497);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(28,0.1058887);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(29,0.1073378);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(30,0.1088367);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(31,0.110378);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(32,0.1119539);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(33,0.1135572);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(34,0.1151801);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(35,0.1168151);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(36,0.1184548);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(37,0.1200916);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(38,0.1217179);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(39,0.1233263);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(40,0.1249091);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(41,0.1264589);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(42,0.1279681);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(43,0.1294291);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(44,0.1308346);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(45,0.1321768);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(46,0.1334483);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(47,0.1346416);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(48,0.135749);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(49,0.1367632);
   h1_SF_UDSG_SSVHEM_1p6to2p4->SetBinError(50,0.1376765);

   h1_SF_UDSG_SSVHPT_lt2p4 = new TH1D("h1_SF_UDSG_SSVHPT_lt2p4","SSVHPT_MistagScaleFactors",50,20,520);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(1,0.8759499);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(2,0.8917812);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(3,0.9061931);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(4,0.9192517);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(5,0.9310229);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(6,0.9415727);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(7,0.9509671);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(8,0.9592722);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(9,0.9665538);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(10,0.972878);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(11,0.9783108);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(12,0.9829181);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(13,0.986766);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(14,0.9899204);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(15,0.9924474);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(16,0.9944128);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(17,0.9958827);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(18,0.9969232);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(19,0.9976001);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(20,0.9979795);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(21,0.9981274);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(22,0.9981097);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(23,0.9979925);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(24,0.9978417);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(25,0.9977232);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(26,0.9977033);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(27,0.9978476);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(28,0.9982224);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(29,0.9988936);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(30,0.999927);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(31,1.001389);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(32,1.003345);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(33,1.005862);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(34,1.009004);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(35,1.01284);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(36,1.017433);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(37,1.022851);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(38,1.029159);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(39,1.036423);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(40,1.04471);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(41,1.054085);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(42,1.064614);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(43,1.076363);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(44,1.089398);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(45,1.103786);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(46,1.119592);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(47,1.136882);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(48,1.155722);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(49,1.176178);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinContent(50,1.198317);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(1,0.1420702);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(2,0.1393674);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(3,0.136849);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(4,0.1345104);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(5,0.1323471);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(6,0.1303548);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(7,0.1285288);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(8,0.1268647);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(9,0.125358);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(10,0.1240042);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(11,0.1227989);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(12,0.1217374);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(13,0.1208155);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(14,0.1200285);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(15,0.1193719);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(16,0.1188414);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(17,0.1184323);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(18,0.1181403);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(19,0.1179608);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(20,0.1178893);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(21,0.1179213);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(22,0.1180524);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(23,0.1182781);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(24,0.1185938);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(25,0.1189952);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(26,0.1194776);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(27,0.1200366);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(28,0.1206678);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(29,0.1213666);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(30,0.1221286);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(31,0.1229492);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(32,0.123824);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(33,0.1247484);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(34,0.1257181);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(35,0.1267284);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(36,0.127775);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(37,0.1288534);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(38,0.1299589);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(39,0.1310872);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(40,0.1322338);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(41,0.1333941);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(42,0.1345637);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(43,0.1357382);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(44,0.1369129);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(45,0.1380835);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(46,0.1392454);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(47,0.1403941);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(48,0.1415252);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(49,0.1426342);
   h1_SF_UDSG_SSVHPT_lt2p4->SetBinError(50,0.1437165);
}


DEFINE_FWK_MODULE(MyAnalyzer);
