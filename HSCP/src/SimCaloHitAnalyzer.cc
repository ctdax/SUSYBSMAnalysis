// -*- C++ -*-
//
// Package:    HSCP/SimCaloHitAnalyzer
// Class:      SimCaloHitAnalyzer
//
/**\class SimCaloHitAnalyzer SimCaloHitAnalyzer.cc HSCP/SimCaloHitAnalyzer/plugins/SimCaloHitAnalyzer.cc

 Description: [Reads SIM ROOTs with gluino samples and outputs histograms for: calorimiter hit location, 4 momenta of R-Hadrons 1 and 2, ...]

 Implementation:
     [Github repository: https://github.com/ctdax/HSCP]
*/
//
// Original Author:  Colby Thompson
//         Created:  Thu, 22 Feb 2024 20:45:33 GMT
//
//


//System include files
#include <memory>
#include <cmath>

//Triggers and Handles
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//Tracker
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

//Calorimiter
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

//Muon Chamber
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTChamber.h"
#include "Geometry/DTGeometry/interface/DTSuperLayer.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

//Calculations
#include "DataFormats/Common/interface/HLTPathStatus.h"
#include "DataFormats/Common/interface/HLTenums.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"

//HSCP specific packages
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
#include "SUSYBSMAnalysis/Analyzer/interface/CommonFunction.h"
#include "SUSYBSMAnalysis/Analyzer/interface/DeDxUtility.h"
#include "SUSYBSMAnalysis/HSCP/interface/HSCPHelpers.h"

//ROOT
#include "TFile.h"
#include "TH1.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TStyle.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TMatrix.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TGraph.h"

//FWCORE
#define FWCORE

// Thresholds for candidate preselection
float globalMaxEta_ = 1.0;
float globalMinPt_ = 55.0;
unsigned int globalMinNOPH_ = 2;
float globalMinFOVH_ = 0.8;
unsigned int globalMinNOM_ = 10;
float globalMaxChi2_ = 5.0;
float globalMaxEoP_ = 0.3; //unsure
float globalMaxDZ_ = 0.1;
float globalMaxDXY_ = 0.02;
float globalMaxTIsol_, globalMiniRelIsoAll_; //unsure
float globalMinIh_ = 3.47; //unsure
float trackProbQCut_;
unsigned int minMuStations_;
float globalMinIs_ = 0.0;
float globalMinTOF_;
float GlobalMinNDOF = 8;            // cut on number of     DegreeOfFreedom used for muon TOF measurement
float GlobalMinNDOFDT = 6;          // cut on number of DT  DegreeOfFreedom used for muon TOF measurement
float GlobalMinNDOFCSC = 6;         // cut on number of CSC DegreeOfFreedom used for muon TOF measurement
float GlobalMaxTOFErr = 0.15;       //0.07;   // cut on error on muon TOF measurement
bool useClusterCleaning = true; //unsure

int evtcount = 0;
int evtcount00 = 0;
int evtcount01 = 0;
int evtcount01n = 0;
int evtcount01c = 0;
int evtcount11 = 0;
int badiqstat = 0;
int passevtcount0 = 0;
int passevtcount1 = 0;
int passevtcount00 = 0;
int passevtcount01 = 0;
int passevtcount01c = 0;
int passevtcount01n = 0;
int passevtcount11 = 0;

int ngencheta09 = 0;
int nsimmatch[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
int nsimnomatch[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
int nsimnomatchqsame[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
int nsimnomatchqdiff[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
double layerr[30] = {0, 3., 7., 11., 16., 26., 33.5, 42., 50., 61., 69., 78., 87., 97., 108., 440., 525., 630., 740., 0., 0., 0.};

struct simhitinfo {
  int ilayer;
  int pdgID;
//  PSimHit *hit;
};

using namespace edm;
using namespace reco;
using namespace std;
using namespace __gnu_cxx;
using namespace trigger;

class TupleMaker;
class MCWeight;

#define DRMATCH 0.002
#define DPTFRACMIN -0.5
#define DPTFRACMAX 1.0

#define TIBEVT1 2
#define TIBEVT2 3
#define TIBEVT3 5
#define TRKEVT1 11
#define TRKEVT2 13
#define TRKEVT3 15
#define BPIXEVT1 32
#define BPIXEVT2 34
#define BPIXEVT3 27

#define NTRLAYERS 19
#define BPIX0 0
#define TIB0 4
#define TOB0 8
#define DT0 14

#define TYPEGLUE 9 
#define TYPEMESON 111
#define TYPEBARYON 1111
#define TYPELEPTON 11

#define ETACUT 0.9

class SimCaloHitAnalyzer : public edm::EDAnalyzer {
public:
  explicit SimCaloHitAnalyzer (const edm::ParameterSet&);
  ~SimCaloHitAnalyzer();


private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  float deltaphi(float phi1, float phi2);
  float deltaphiabs(float phi1, float phi2);
  float deltaRta(float phi1, float eta1, float phi2, float eta2);
  int determineCharge(int pid);
  int determineType(int pid);
  bool passHSCPPreselection(int typeMode, susybsm::HSCParticle hscpCan, reco::VertexCollection vertColl, reco::VertexCollection inclusiveSecondaryVertices, reco::PFCandidateCollection pf, reco::DeDxData *dedxSObj, reco::DeDxData *dedxMObj, bool passedCutsArray[]);
  bool comparePID(pair<int, int> p1, pair<int, int> p2); 

  edm::EDGetTokenT<vector<reco::GenParticle>> genParticlesToken_;

  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tGeomEsToken_;

  // Tracker hits
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTIBLow_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTIBHigh_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTOBLow_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTOBHigh_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTIDLow_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTIDHigh_Token_;
  edm::EDGetTokenT<edm::PSimHitContainer> edmPSimHitContainer_siTECLow_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_siTECHigh_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_pxlBrlLow_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_pxlBrlHigh_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_pxlFwdLow_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_pxlFwdHigh_Token_;

  // Calorimiter hits
  edm::EDGetTokenT <edm::PSimHitContainer> edmCaloHitContainer_EcalHitsEB_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmCaloHitContainer_EcalHitsEE_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmCaloHitContainer_EcalHitsES_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmCaloHitContainer_HcalHits_Token_;


  // Muon hits
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_muonCSC_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_muonDT_Token_;

  // Miscellaneous
  edm::EDGetTokenT<edm::SimTrackContainer> edmSimTrackContainerToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerSummary_;
  edm::EDGetTokenT<std::vector<reco::PFMET>> pfmet_;
  edm::EDGetTokenT<std::vector<reco::CaloMET>> calomet_;
  edm::EDGetTokenT<VertexCollection> vertexToken_;
  edm::EDGetTokenT<VertexCollection> inclusiveSecondaryvertexToken_;
  edm::EDGetTokenT<PFCandidateCollection> pfToken_;
  edm::EDGetTokenT<std::vector<reco::Track>> gentrkToken_;
  edm::EDGetTokenT<vector<susybsm::HSCParticle>> hscpCandToken_;
  edm::EDGetTokenT<reco::DeDxHitInfoAss> dedxToken_;

  //declare variables here
#define TLEN 1625
  string trig_name[TLEN] = {};
  int trig_pass[TLEN]    = {};
  int trig_total[TLEN]   = {};
  int trig_pass_or = 0;
  int trig_pass_and = 0;
  //int trig_total_special = 0;

   std::vector<std::pair<int, int>> genrhadcounts;

  ///double met_values[];

  int typeMode_;
  int sampleType_;
  string period_;
  int isData;
  bool useTemplateLayer_;
  bool skipPixel_;
  string dEdxTemplate_;
  bool  enableDeDxCalibration_;
  string dEdxCalibration_;

  dedxGainCorrector trackerCorrector;
  TH3F* dEdxTemplates = nullptr;
  float dEdxSF_0_, dEdxSF_1_;
  float dEdxSF[2] = {dEdxSF_0_, dEdxSF_1_};

  TFile *outputFile_;
  TH1F *EBHit_ieta;
  TH1F *EBHit_iphi;
  TH1F *EBHit_x;
  TH1F *EBHit_y;
  TH1F *EBHit_z;
  TH1F *RHadron1_px;
  TH1F *RHadron1_py;
  TH1F *RHadron1_pz;
  TH1F *RHadron2_px;
  TH1F *RHadron2_py;
  TH1F *RHadron2_pz;

};

int SimCaloHitAnalyzer::determineCharge(int pdgID) {
  int charge = -99;
  int pid = abs(pdgID);
  if (pid==1000993) {charge = 0;} // ~g_glueball
  else if (pid==1009113) {charge = 0;} // ~g_rho0
  else if (pid==1009223) {charge = 0;} // ~g_omega
  else if (pid==1009213) {charge = 1;} // ~g_rho+
  else if (pid==1009313) {charge = 0;} // ~g_K*0
  else if (pid==1009323) {charge = 1;} // ~g_K*+
  else if (pid==1009333) {charge = 0;} // ~g_phi
  else if (pid==1091114) {charge = -1;} // ~g_Delta-
  else if (pid==1092114) {charge = 0;} // ~g_Delta0
  else if (pid==1092214) {charge = 1;} // ~g_Delta+
  else if (pid==1092224) {charge = 2;} // ~g_Delta++
  else if (pid==1093114) {charge = -1;} // ~g_Sigma*-
  else if (pid==1093214) {charge = 0;} // ~g_Sigma*0
  else if (pid==1093224) {charge = 1;} // ~g_Sigma*+
  else if (pid==1093314) {charge = -1;} // ~g_Xi*-
  else if (pid==1093324) {charge = 0;} // ~g_Xi*0
  else if (pid==1093334) {charge = -1;} // ~g_Omega-
  else if (pid==1000015) {charge = -1;} // ~stau-
  else if (pid==1000612) {charge = +1;} // ~T+
  else if (pid==1000622 ) {charge = 0;} // ~T0 
  else if (pid==1000632 ) {charge = +1;} // ~T_s+
  else if (pid==1000642 ) {charge = 0;} // ~T_c0
  else if (pid==1000652 ) {charge = +1;} // ~T_b+
  else if (pid==1006113 ) {charge = 0;} // ~T_dd10
  else if (pid==1006211 ) {charge = +1;} // ~T_ud0+
  else if (pid==1006213 ) {charge = +1;} // ~T_ud1+
  else if (pid==1006223 ) {charge = +2;} // ~T_uu1++ 
  else if (pid==1006311 ) {charge = 0;} // ~T_sd00  
  else if (pid==1006313 ) {charge = 0;} // ~T_sd10  
  else if (pid==1006321 ) {charge = +1;} // ~T_su0+  
  else if (pid==1006323 ) {charge = +1;} // ~T_su1+  
  else if (pid==1006333 ) {charge = 0;} // ~T_ss10  
/* 1000021   1500.0     # ~g
      1000993   1500.700   # ~g_glueball
      1009213   1500.650   # ~g_rho+
      1009313   1500.825   # ~g_K*0
      1009323   1500.825   # ~g_K*+
      1009113   1500.650   # ~g_rho0
      1009223   1500.650   # ~g_omega
      1009333   1501.800   # ~g_phi
      1091114   1500.975   # ~g_Delta-
      1092114   1500.975   # ~g_Delta0
      1092214   1500.975   # ~g_Delta+
      1092224   1500.975   # ~g_Delta++
      1093114   1501.150   # ~g_Sigma*-
      1093214   1501.150   # ~g_Sigma*0
      1093224   1501.150   # ~g_Sigma*+
      1093314   1501.300   # ~g_Xi*-
      1093324   1501.300   # ~g_Xi*0
      1093334   1501.600   # ~g_Omega- 
  stop
        1000006 800.000   # ~t_1
	1000612 800.325   # ~T+  
	1000622 800.325   # ~T0  
	1000632 800.500   # ~T_s+
	1000642 801.500   # ~T_c0
	1000652 804.800   # ~T_b+
	1006113 800.650   # ~T_dd10
	1006211 800.650   # ~T_ud0+
	1006213 800.650   # ~T_ud1+
	1006223 800.650   # ~T_uu1++ 
	1006311 800.825   # ~T_sd00  
	1006313 800.825   # ~T_sd10  
	1006321 800.825   # ~T_su0+  
	1006323 800.825   # ~T_su1+  
	1006333 801.000   # ~T_ss10  
        -1000006 800.000   # ~t_1bar
	-1000612 800.325   # ~Tbar-  
	-1000622 800.325   # ~Tbar0  
	-1000632 800.500   # ~Tbar_s-
	-1000642 801.500   # ~Tbar_c0
	-1000652 804.800   # ~Tbar_b-
	-1006113 800.650   # ~Tbar_dd10
	-1006211 800.650   # ~Tbar_ud0-
	-1006213 800.650   # ~Tbar_ud1-
	-1006223 800.650   # ~Tbar_uu1-- 
	-1006311 800.825   # ~Tbar_sd00  
	-1006313 800.825   # ~Tbar_sd10  
	-1006321 800.825   # ~Tbar_su0-  
	-1006323 800.825   # ~Tbar_su1-  
	-1006333 801.000   # ~Tbar_ss10  

*/
  //else {cout << "Unknown PID = " << pid << endl;}
  if ((charge>-98)&&(pdgID<0)) charge *= -1;
  return(charge);
}

int SimCaloHitAnalyzer::determineType(int pdgID) {
     int pidabs = abs(pdgID);
     if ((pidabs>1000900)&&(pidabs<1000999)) {
       return(TYPEGLUE);
     } else if ((pidabs>1000999)&&(pidabs<1009999)) {
       return(TYPEMESON);
     } else if ((pidabs>1009999)&&(pidabs<1099999)) {
       return(TYPEBARYON);
     } else if ((pidabs>1000009)&&(pidabs<1000020)) {
       return(TYPELEPTON);
     }
  return(-1);
}

// delta(phi1-phi2) - return value from -pi to +pi
float SimCaloHitAnalyzer::deltaphi(float phi1, float phi2) {
  float tdelta = phi1 - phi2;
  while (tdelta<(-1.0*M_PI)) {tdelta += (2.0*M_PI);}
  while (tdelta>(1.0*M_PI)) {tdelta -= (2.0*M_PI);}
  return(tdelta);
}

// |delta(phi1-phi2)| - return value from 0 to +pi
float SimCaloHitAnalyzer::deltaphiabs(float phi1, float phi2) {
  float tdphi = deltaphi(phi1,phi2);
  return(fabs(tdphi));
}

float SimCaloHitAnalyzer::deltaRta(float phi1, float eta1, float phi2, float eta2) {
  float diffphi1 = deltaphiabs(phi1,phi2);
  float diffeta1 = eta1 - eta2;
  float tmpdeltaR1 = sqrt(diffphi1*diffphi1 + diffeta1*diffeta1);
  return(tmpdeltaR1);
}

// preselection
bool SimCaloHitAnalyzer::passHSCPPreselection(int typeMode, susybsm::HSCParticle hscpCan, reco::VertexCollection vertexColl, reco::VertexCollection inclusiveSecondaryVertices, reco::PFCandidateCollection pf, reco::DeDxData *dedxSObj, reco::DeDxData *dedxMObj, bool passedCutsArray[]) {

  reco::TrackRef track = hscpCan.trackRef();
  if (track.isNull()) return(false);

    // find closest vertex
  int closestZGoodVertex = -1;
  int goodVerts = 0;
  float dzMin = 10000;
    // Loop on the vertices in the event
  for (unsigned int i = 0; i < vertexColl.size(); i++) {
    if (vertexColl[i].isFake() || fabs(vertexColl[i].z()) > 24 || vertexColl[i].position().rho() > 2 || vertexColl[i].ndof() <= 4) continue;  //only consider good vertex
    goodVerts++;
    
    if (fabs(track->dz(vertexColl[i].position())) < fabs(dzMin)) {
      dzMin = fabs(track->dz(vertexColl[i].position()));
      closestZGoodVertex = i;
    }
  } // End loop on the vertices in the event
  
  if (closestZGoodVertex < 0) {
    closestZGoodVertex = 0;
  }
  
  // Impact paramters dz and dxy
  float dz = track->dz(vertexColl[closestZGoodVertex].position());
  float dxy = track->dxy(vertexColl[closestZGoodVertex].position());

  bool isMaterialTrack = false;

  // Loop on the secondary vertices to find tracks that original from the pixel layers
  // i.e. are due to NI
  for (unsigned int i = 0; i < inclusiveSecondaryVertices.size(); i++) {
    if (inclusiveSecondaryVertices[i].isFake()) {
      continue;
    }
    auto rho = inclusiveSecondaryVertices[i].position().rho();
    if ( (( 2.80-0.075 ) < rho && rho < ( 3.10+0.075 )) || (( 6.60-0.075 ) < rho && rho < ( 7.00+0.075 ))
        || (( 10.9-0.075 ) < rho && rho < ( 10.9+0.075 )) || (( 16.0-0.075 ) < rho && rho < ( 16.0+0.075 )) ) {
      for( const auto& rf_track : inclusiveSecondaryVertices[i].refittedTracks() ) {
        const reco::Track& origTrk = *( inclusiveSecondaryVertices[i].originalTrack( rf_track ));
        if( track->pt() == origTrk.pt() ){
          isMaterialTrack = true;
          break;
        } 
      } 
    } else {
      continue;
    } 
  } 

    // Save PF informations and isolation
    float pfIsolation_DZ_ = 0.1;

    float track_PFIso005_sumCharHadPt = 0, track_PFIso005_sumNeutHadPt = 0, track_PFIso005_sumPhotonPt = 0, track_PFIso005_sumPUPt = 0;
    float track_PFIso01_sumCharHadPt = 0, track_PFIso01_sumNeutHadPt = 0, track_PFIso01_sumPhotonPt = 0, track_PFIso01_sumPUPt = 0;
    float track_PFIso03_sumCharHadPt = 0, track_PFIso03_sumNeutHadPt = 0, track_PFIso03_sumPhotonPt = 0, track_PFIso03_sumPUPt = 0;
    float track_PFIso05_sumCharHadPt = 0, track_PFIso05_sumNeutHadPt = 0, track_PFIso05_sumPhotonPt = 0, track_PFIso05_sumPUPt = 0;

    float track_PFMiniIso_sumCharHadPt = 0, track_PFMiniIso_sumNeutHadPt = 0, track_PFMiniIso_sumPhotonPt = 0, track_PFMiniIso_sumPUPt = 0, track_PFMiniIso_sumMuonPt = 0;
    float pf_energy=0;

    float RMin = 9999.;
    unsigned int idx_pf_RMin = 9999;

      for (unsigned int i = 0; i < pf.size(); i++){
          const reco::PFCandidate pfCand = pf[i];
          float dr = deltaR(pfCand.eta(),pfCand.phi(),track->eta(),track->phi());
          if(dr < RMin){
              RMin = dr;
              idx_pf_RMin = i;
           }
       }//end loop PFCandidates

      // https://urldefense.proofpoint.com/v2/url?u=https-3A__github.com_cms-2Dsw_cmssw_blob_72d0fc00976da53d1fb745eb7f37b2a4ad965d7e_&d=DwIGAg&c=gRgGjJ3BkIsb5y6s49QqsA&r=iYf-W5o_XDmJ-u42pi_qsGfw1LmudygYMQ_Az9XJWjc&m=oBFJTDv2gfKrlPynZ9fSMGuSOoJZcNBR8epoz7KKL7P7LrQKFQ2Uhpkq3B3sODWH&s=Y8xRnPaFLOTBqfLjih3mXkPcBp2SP-_OPNc8fWzBgOg&e= 
      // PhysicsTools/PatAlgos/plugins/PATIsolatedTrackProducer.cc#L555
//      for(unsigned int i=0;i<pf->size();i++){
      for (unsigned int i=0;i<pf.size();i++){
        const reco::PFCandidate* pfCand = &(pf)[i];

        pf_energy = pfCand->ecalEnergy() + pfCand->hcalEnergy();
        bool pf_isPhotonForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::gamma;
        bool pf_isChHadronForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h;
        bool pf_isNeutHadronForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h0;
        bool pf_isMuonForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::mu;

        if(i == idx_pf_RMin) continue; //don't count itself
        float dr = deltaR(pfCand->eta(),pfCand->phi(),track->eta(),track->phi());
        bool fromPV = (fabs(dz) < pfIsolation_DZ_);
        int id = std::abs(pfCand->pdgId());
        float pt = pfCand->p4().pt();
        if(dr<0.05){
            // charged cands from PV get added to trackIso
          if(id == 211 && fromPV) track_PFIso005_sumCharHadPt+=pt;
            // charged cands not from PV get added to pileup iso
          else if(id == 211) track_PFIso005_sumPUPt+=pt;
            // neutral hadron iso
          if(id == 130) track_PFIso005_sumNeutHadPt+=pt;
            // photon iso
          if(id == 22) track_PFIso005_sumPhotonPt+=pt;
        }if(dr<0.1){
            if(id == 211 && fromPV) track_PFIso01_sumCharHadPt+=pt;
            else if(id == 211) track_PFIso01_sumPUPt+=pt;
            if(id == 130) track_PFIso01_sumNeutHadPt+=pt;
            if(id == 22) track_PFIso01_sumPhotonPt+=pt;
        }if(dr<0.3){
            if(id == 211 && fromPV) track_PFIso03_sumCharHadPt+=pt;
            else if(id == 211) track_PFIso03_sumPUPt+=pt;
            if(id == 130) track_PFIso03_sumNeutHadPt+=pt;
            if(id == 22) track_PFIso03_sumPhotonPt+=pt;
        }if(dr<0.5){
            if(id == 211 && fromPV) track_PFIso05_sumCharHadPt+=pt;
            else if(id == 211) track_PFIso05_sumPUPt+=pt;
            if(id == 130) track_PFIso05_sumNeutHadPt+=pt;
            if(id == 22) track_PFIso05_sumPhotonPt+=pt;
        }

        float drForMiniIso = 0.0;
        if (track->pt() < 50 ) {
          drForMiniIso = 0.2;
        } else if (track->pt() < 200) {
          drForMiniIso = 10/track->pt();
        } else {
          drForMiniIso = 0.05;
        }
        if (dr<drForMiniIso) {
            // charged cands from PV get added to trackIso
          if(pf_isChHadronForIdx && fromPV) track_PFMiniIso_sumCharHadPt+=pt;
            // charged cands not from PV get added to pileup iso
          else if(pf_isChHadronForIdx) track_PFMiniIso_sumPUPt+=pt;
            // neutral hadron iso
          if(pf_isNeutHadronForIdx) track_PFMiniIso_sumNeutHadPt+=pt;
            // photon iso
          if(pf_isPhotonForIdx) track_PFMiniIso_sumPhotonPt+=pt;
            // muon iso
          if(pf_isMuonForIdx) track_PFMiniIso_sumMuonPt+=pt;
        }
      }//end loop PFCandidates
//    }
            

  // Calculate PF mini relative isolation
  float miniRelIsoAll = (track_PFMiniIso_sumCharHadPt + std::max(0.0, track_PFMiniIso_sumNeutHadPt + track_PFMiniIso_sumPhotonPt - 0.5* track_PFMiniIso_sumPUPt))/track->pt();

  float EoP = pf_energy / track->p();

// Number of DeDx hits
  unsigned int numDeDxHits = (dedxSObj) ? dedxSObj->numberOfMeasurements() : 0;

  float Ih = (dedxMObj) ?  dedxMObj->dEdx() : 0.0;

  // No cut, i.e. events after trigger
  passedCutsArray[0]  = true;
  // Check if eta is inside the max eta cut
  passedCutsArray[1]  = (fabs(track->eta()) < globalMaxEta_) ? true : false;
  // Cut on minimum track pT
  passedCutsArray[2]  = (track->pt() > globalMinPt_) ? true : false;
  // Check the number of pixel hits
  passedCutsArray[3]  = (typeMode != 3 && fabs(track->hitPattern().numberOfValidPixelHits()) > globalMinNOPH_) ? true : false;
  // Check the min fraction of valid hits
  passedCutsArray[4]  = (typeMode != 3 && track->validFraction() > globalMinFOVH_) ? true : false;
  // Cut for the number of dEdx hits
  passedCutsArray[5]  = (numDeDxHits >= globalMinNOM_)  ? true : false;
  // Select only high purity tracks
  passedCutsArray[6]  = (typeMode != 3 && track->quality(reco::TrackBase::highPurity)) ? true : false;
  // Cut on the chi2 / ndof
  passedCutsArray[7] = (typeMode != 3 && track->chi2() / track->ndof() < globalMaxChi2_) ? true : false;

  // Cut on the energy over momenta
  passedCutsArray[8] = (EoP < globalMaxEoP_) ? true : false;

 // Cut on the impact parameter
 // for typeMode_ 5 dz is supposed to come from the beamspot, TODO
  passedCutsArray[9] = (  (typeMode != 5 && fabs(dz) < globalMaxDZ_)
                       || (typeMode == 5 && fabs(dz) < 4)) ? true : false;
 // for typeMode_ 5 dxy is supposed to come from the beamspot, TODO
  passedCutsArray[10] = (  (typeMode != 5 && fabs(dxy) < globalMaxDXY_)
                         || (typeMode == 5 && fabs(dxy) < 4)) ? true : false;

 // Cut on the tracker based isolation
  passedCutsArray[12] = (!isMaterialTrack) ? true : false;

 // Cut on the PF based mini-isolation
  passedCutsArray[13] = ( miniRelIsoAll < globalMiniRelIsoAll_) ? true : false;
 // Cut on the PF electron ID

 // Cut on min Ih (or max for fractionally charged)
  passedCutsArray[15] = (  (typeMode != 5 &&  Ih > globalMinIh_)
                         || (typeMode == 5 && Ih < globalMinIh_)) ? true : false;
 
  return(true);
}

bool SimCaloHitAnalyzer::comparePID(pair<int, int> p1, pair<int, int> p2) {
  return (p1.first < p2.first);
}

//constructor
SimCaloHitAnalyzer::SimCaloHitAnalyzer(const edm::ParameterSet& iConfig) 
 :  tGeomEsToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>())


{
  // Tracker
  edmPSimHitContainer_siTIBLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTIBLowTof"));
  edmPSimHitContainer_siTIBHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTIBHighTof"));
  edmPSimHitContainer_siTOBLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTOBLowTof"));
  edmPSimHitContainer_siTOBHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTOBHighTof"));
  edmPSimHitContainer_siTIDLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTIDLowTof"));
  edmPSimHitContainer_siTIDHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTIDHighTof"));
  edmPSimHitContainer_siTECLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTECLowTof"));
  edmPSimHitContainer_siTECHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsTECHighTof"));
  edmPSimHitContainer_pxlBrlLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsPixelBarrelLowTof"));
  edmPSimHitContainer_pxlBrlHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsPixelBarrelHighTof"));
  edmPSimHitContainer_pxlFwdLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsPixelEndcapLowTof"));
  edmPSimHitContainer_pxlFwdHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("TrackerHitsPixelEndcapHighTof"));

  // Calorimiter
  edmCaloHitContainer_EcalHitsEB_Token_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("EcalHitsEB"));
  edmCaloHitContainer_EcalHitsEE_Token_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("EcalHitsEE"));
  edmCaloHitContainer_EcalHitsES_Token_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("EcalHitsES"));
  edmCaloHitContainer_HcalHits_Token_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("HcalHits"));

  // Muon Chamber
  edmPSimHitContainer_muonCSC_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("MuonCSCHits"));
  edmPSimHitContainer_muonDT_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("MuonDTHits"));

  // Miscellaneous
  edmSimTrackContainerToken_ = consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("G4TrkSrc"));
  triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"));
  triggerSummary_ = consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("trig_sum"));
  pfmet_ = consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("pfmet_reco"));
  calomet_ = consumes<std::vector<reco::CaloMET>>(iConfig.getParameter<edm::InputTag>("calomet_reco"));
  vertexToken_ = mayConsume<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex_reco"));
  inclusiveSecondaryvertexToken_ = mayConsume<VertexCollection>(iConfig.getParameter<edm::InputTag>("secondaryvertex_reco"));
  genParticlesToken_ = mayConsume<vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen_info"));
  pfToken_ = mayConsume<PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pf_reco"));
  gentrkToken_ = mayConsume<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("gentrk_reco"));
  hscpCandToken_ = mayConsume<vector<susybsm::HSCParticle>>(iConfig.getParameter<edm::InputTag>("hscp_cand"));
  dedxToken_ = consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("DedxCollection"));

  typeMode_ = iConfig.getUntrackedParameter<int>("TypeMode");
  sampleType_ = iConfig.getUntrackedParameter<int>("SampleType");
  period_ = iConfig.getUntrackedParameter<string>("Period");
  useTemplateLayer_ = iConfig.getUntrackedParameter<bool>("UseTemplateLayer");
  skipPixel_ = iConfig.getUntrackedParameter<bool>("SkipPixel");
  dEdxTemplate_ = iConfig.getUntrackedParameter<string>("DeDxTemplate");
  enableDeDxCalibration_ = iConfig.getUntrackedParameter<bool>("EnableDeDxCalibration");
  dEdxCalibration_ = iConfig.getUntrackedParameter<string>("DeDxCalibration");
  dEdxSF_0_ = iConfig.getUntrackedParameter<double>("DeDxSF_0");
  dEdxSF_1_ = iConfig.getUntrackedParameter<double>("DeDxSF_1");

  isData = (sampleType_ == 0);
  dEdxSF[0] = dEdxSF_0_;
  dEdxSF[1] = dEdxSF_1_;


  // Output File
  cout << " passevtcount01 = " << passevtcount01 << " evtcount01 = " << evtcount01 << " fraction = " << fr << endl;
  outputFile_ = new TFile("/uscms/home/cthompso/nobackup/CMSSW_10_6_30/src/HSCP/Saturation/SimCaloHitAnalysis/test/RHadron.root");  

  // Declare ROOT histograms
  EBHit_ieta = new TH1F("EBHit_ieta","EB_SimHits_ieta",100,-2.,2.)
  EBHit_iphi = new TH1F("EBHit_iphi","EB_SimHits_iphi",100,-3.5,3.5)
  EBHit_x = new TH1F("EBHit_x","EB_SimHits_x",100,-3.5,3.5)
  EBHit_y = new TH1F("EBHit_y","EB_SimHits_y",100,-3.5,3.5)
  EBHit_z = new TH1F("EBHit_z","EB_SimHits_z",100,-3.5,3.5)
  RHadron1_px = new TH1F("RHadron1_px","RHadron1_px",100,-100.,100.)
  RHadron1_py = new TH1F("RHadron1_py","RHadron1_py",100,-100.,100.)
  RHadron1_pz = new TH1F("RHadron1_pz","RHadron1_pz",100,-100.,100.)
  RHadron2_px = new TH1F("RHadron2_px","RHadron2_px",100,-100.,100.)
  RHadron2_py = new TH1F("RHadron2_py","RHadron2_py",100,-100.,100.)
  RHadron2_pz = new TH1F("RHadron2_pz","RHadron2_pz",100,-100.,100.)

  
  evtcount = 0;
  evtcount00 = 0;
  evtcount01 = 0;
  evtcount01n = 0;
  evtcount01c = 0;
  evtcount11 = 0;
  badiqstat = 0;
  passevtcount0 = 0;
  passevtcount1 = 0;
  passevtcount00 = 0;
  passevtcount01 = 0;
  passevtcount01n = 0;
  passevtcount01c = 0;
  passevtcount11 = 0;

  ngencheta09 = 0;

}


//destructor
SimCaloHitAnalyzer::~SimCaloHitAnalyzer() {

  float fr = (float)passevtcount0/(float)evtcount;
  //cout << " passevtcount0 = " << passevtcount0 << " evtcount = " << evtcount << " fraction = " << fr << endl;
  fr = (float)passevtcount1/(float)evtcount;
  //cout << " passevtcount1 = " << passevtcount1 << " evtcount = " << evtcount << " fraction = " << fr << endl;
  fr = (float)passevtcount00/(float)evtcount00;
  //cout << " passevtcount00 = " << passevtcount00 << " evtcount00 = " << evtcount00 << " fraction = " << fr << endl;
  fr = (float)passevtcount01/(float)evtcount01;
  //cout << " passevtcount01 = " << passevtcount01 << " evtcount01 = " << evtcount01 << " fraction = " << fr << endl;
  fr = (float)passevtcount11/(float)evtcount11;
  //cout << " passevtcount11 = " << passevtcount11 << " evtcount11 = " << evtcount11 << " fraction = " << fr << endl;
  //cout << " badiqstat count = " << badiqstat << endl;

  //cout << endl;
  //cout << " # of charged R-hadron |eta|<0.9 = " << ngencheta09 << endl;
  nsimmatch[0] = ngencheta09;
  nsimnomatch[0] = ngencheta09;
  nsimnomatchqsame[0] = ngencheta09;
  nsimnomatchqdiff[0] = ngencheta09;
  double nsimmatch2[NTRLAYERS];
  double nsimnomatch2[NTRLAYERS];
  double nsimnomatchqsame2[NTRLAYERS];
  double nsimnomatchqdiff2[NTRLAYERS];
  double nsimtotal2[NTRLAYERS];
  double nsimmatchfrac[NTRLAYERS];
  double nsimnomatchfrac[NTRLAYERS];
  double nsimnomatchqsamefrac[NTRLAYERS];
  double nsimnomatchqdifffrac[NTRLAYERS];
  double nsimtotalfrac[NTRLAYERS];
  for (int k=0; k<NTRLAYERS; k++) {
    //cout << "   # of matched hits in layer " << k << " = " << nsimmatch[k] << endl;
    nsimmatch2[k] = (double)nsimmatch[k];
    nsimnomatch2[k] = (double)nsimnomatch[k];
    nsimnomatchqsame2[k] = (double)nsimnomatchqsame[k];
    nsimnomatchqdiff2[k] = (double)nsimnomatchqdiff[k];
    if (k==0) {
      nsimtotal2[0] = (double)nsimmatch[0];
    } else {
      nsimtotal2[k] = (double)nsimmatch[k] + (double)nsimnomatch[k];
    }
    nsimmatchfrac[k] = (double)nsimmatch[k]/(double)ngencheta09;
    nsimnomatchfrac[k] = (double)nsimnomatch[k]/(double)ngencheta09;
    nsimnomatchqsamefrac[k] = (double)nsimnomatchqsame[k]/(double)ngencheta09;
    nsimnomatchqdifffrac[k] = (double)nsimnomatchqdiff[k]/(double)ngencheta09;
    nsimtotalfrac[k] = (double)nsimtotal2[k]/(double)ngencheta09;
  }
  //cout << endl;
   
   std::vector<std::pair<int, int>> sortgenpairs;
   int psize = genrhadcounts.size();
   int lastlow = -999999999;
   for (int i=0; i<psize; i++) {
     pair<int, int> plow(lastlow,0);
     int lowid = 999999999;
     for (std::vector<std::pair<int, int>>::iterator it = genrhadcounts.begin(); (it != genrhadcounts.end()); ++it) {
       int id = (*it).first;
       if ((id<lowid)&&(id>lastlow)) {
         plow = (*it);
         lowid = id;
       }
     }
     sortgenpairs.push_back(plow);
     lastlow = plow.first;

  // Write histograms and file
  outputFile_->cd();

  EBHit_ieta->Write();
  EBHit_iphi->Write();
  EBHit_x->Write();
  EBHit_y->Write();
  EBHit_z->Write();
  RHadron1_px->Write();
  RHadron1_py->Write();
  RHadron1_pz->Write();
  RHadron2_px->Write();
  RHadron2_py->Write();
  RHadron2_pz->Write();

  std::cout << "saving MET trigger histograms..." << std::endl;
  outputFile_->Write();
  outputFile_->Close();

}

void SimCaloHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   evtcount++;

   // Generator particles
   edm::Handle<vector<reco::GenParticle>> genColl;
   iEvent.getByToken(genParticlesToken_, genColl);
//   std::cout << " Gen particles: " << std::endl;

   vector<simhitinfo> linkedhits1;
   vector<simhitinfo> linkedhits2;

// simhit code inspired by
//   Validation/TrackerHits/src/TrackerHitAnalyzer.cc

   // PSimHits
   /////////////////////////////////
   // get Silicon TIB information
   //////////////////////////////////
   // extract TIB low container
   edm::Handle<edm::PSimHitContainer> SiTIBLowContainer;
   iEvent.getByToken(edmPSimHitContainer_siTIBLow_Token_, SiTIBLowContainer);
   if (!SiTIBLowContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTIBLowTof in event!";
     return;
   }
   //////////////////////////////////
   // extract TIB high container
   edm::Handle<edm::PSimHitContainer> SiTIBHighContainer;
   iEvent.getByToken(edmPSimHitContainer_siTIBHigh_Token_, SiTIBHighContainer);
   if (!SiTIBHighContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTIBHighTof in event!";
     return;
   }

   ///////////////////////////////////
   // get Silicon TOB information
   //////////////////////////////////
   // extract TOB low container
   edm::Handle<edm::PSimHitContainer> SiTOBLowContainer;
   iEvent.getByToken(edmPSimHitContainer_siTOBLow_Token_, SiTOBLowContainer);
   if (!SiTOBLowContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTOBLowTof in event!";
     return;
   }
   //////////////////////////////////
   // extract TOB high container
   edm::Handle<edm::PSimHitContainer> SiTOBHighContainer;
   iEvent.getByToken(edmPSimHitContainer_siTOBHigh_Token_, SiTOBHighContainer);
   if (!SiTOBHighContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTOBHighTof in event!";
     return;
   }

   //
   /////////////////////////////////////
   // get Silicon TID information
   //////////////////////////////////
   // extract TID low container
   edm::Handle<edm::PSimHitContainer> SiTIDLowContainer;
   iEvent.getByToken(edmPSimHitContainer_siTIDLow_Token_, SiTIDLowContainer);
   if (!SiTIDLowContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTIDLowTof in event!";
     return;
   }
   //////////////////////////////////
   // extract TID high container
   edm::Handle<edm::PSimHitContainer> SiTIDHighContainer;
   iEvent.getByToken(edmPSimHitContainer_siTIDHigh_Token_, SiTIDHighContainer);
   if (!SiTIDHighContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTIDHighTof in event!";
     return;
   }

   ///////////////////////////////////
   // get Silicon TEC information
   //////////////////////////////////
   // extract TEC low container
   edm::Handle<edm::PSimHitContainer> SiTECLowContainer;
   iEvent.getByToken(edmPSimHitContainer_siTECLow_Token_, SiTECLowContainer);
   if (!SiTECLowContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTECLowTof in event!";
     return;
   }
   //////////////////////////////////
   // extract TEC high container
   edm::Handle<edm::PSimHitContainer> SiTECHighContainer;
   iEvent.getByToken(edmPSimHitContainer_siTECHigh_Token_, SiTECHighContainer);
   if (!SiTECHighContainer.isValid()) {
     edm::LogError("TrackerHitProducer::analyze") << "Unable to find TrackerHitsTECHighTof in event!";
     return;
   }

   /////////////////////////////////
   // get Pixel Barrel information
   ////////////////////////////////
   // extract low container
   edm::Handle<edm::PSimHitContainer> PxlBrlLowContainer;
   iEvent.getByToken(edmPSimHitContainer_pxlBrlLow_Token_, PxlBrlLowContainer);
   if (!PxlBrlLowContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find TrackerHitsPixelBarrelLowTof in event!";
     return;
   }
   // extract high container
   edm::Handle<edm::PSimHitContainer> PxlBrlHighContainer;
   iEvent.getByToken(edmPSimHitContainer_pxlBrlHigh_Token_, PxlBrlHighContainer);
   if (!PxlBrlHighContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find TrackerHitsPixelBarrelHighTof in event!";
     return;
   }
   /////////////////////////////////
   // get Pixel Forward information
   ////////////////////////////////
   // extract low container
   edm::Handle<edm::PSimHitContainer> PxlFwdLowContainer;
   iEvent.getByToken(edmPSimHitContainer_pxlFwdLow_Token_, PxlFwdLowContainer);
   if (!PxlFwdLowContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find TrackerHitsPixelEndcapLowTof in event!";
     return;
   }
   // extract high container
   edm::Handle<edm::PSimHitContainer> PxlFwdHighContainer;
   iEvent.getByToken(edmPSimHitContainer_pxlFwdHigh_Token_, PxlFwdHighContainer);
   if (!PxlFwdHighContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find TrackerHitsPixelEndcapHighTof in event!";
     return;
   }
    ////////////////////////////////
    // Get calorimiter information//
    ////////////////////////////////
    //EcalEB
    edm::Handle<edm::PSimHitContainer> EcalEBContainer;
    iEvent.getByToken(edmCaloHitContainer_EcalHitsEB_Token_, EcalEBContainer);
    if (!EcalEBContainer.isValid()) {
        edm::LogError("TrackerHitProducer::analyze") << "Unable to find EcalHitsEB in event!";
        return;
   }
    //EcalEE
    edm::Handle<edm::PSimHitContainer> EcalEEContainer;
    iEvent.getByToken(edmCaloHitContainer_EcalHitsEE_Token_, EcalEEContainer);
    if (!EcalEEContainer.isValid()) {
        edm::LogError("TrackerHitProducer::analyze") << "Unable to find EcalHitsEE in event!";
        return;
   }
    //EcalES
    edm::Handle<edm::PSimHitContainer> EcalESContainer;
    iEvent.getByToken(edmCaloHitContainer_EcalHitsES_Token_, EcalESContainer);
    if (!EcalESContainer.isValid()) {
        edm::LogError("TrackerHitProducer::analyze") << "Unable to find EcalHitsES in event!";
        return;
   }
    //Hcal
    edm::Handle<edm::PSimHitContainer> HcalContainer;
    iEvent.getByToken(edmCaloHitContainer_HcalHits_Token_, HcalContainer);
    if (!HcalContainer.isValid()) {
        edm::LogError("TrackerHitProducer::analyze") << "Unable to find HcalHits in event!";
        return;
   }

   // Import tracker and calorimiter geometry
   edm::ESHandle<TrackerGeometry> tkGeometry;
   iSetup.get<TrackerDigiGeometryRecord>().get(tkGeometry);

   edm::ESHandle<CaloGeometry> caloGeometry;
   iSetup.get<IdealGeometryRecord>().get(caloGeometry); 

   // Find R-hadrons
   const reco::GenParticle *genrhad1=0;
   const reco::GenParticle *genrhad2=0;

   for (unsigned int ig=0;ig<genColl->size();ig++) {
     const reco::GenParticle & part = (*genColl)[ig];
     if ((abs(part.pdgId())>999999)&&(part.status()==1)) {
   //cout << "SUSY part = " << part.pdgId() << endl;
       if (genrhad1==0) {
         genrhad1 = &(*genColl)[ig];
       } else if (genrhad2==0) {
         genrhad2 = &(*genColl)[ig];
       } else {
         //std::cout << "WARNING - found more than two R-hadrons" << std::endl;
       }
         // find and increment count
       bool found = false;

       for (std::vector<std::pair<int, int>>::iterator it = genrhadcounts.begin(); (it != genrhadcounts.end())&&(!found); ++it) {
         if (part.pdgId()==(*it).first) {
           (*it).second++;
           found = true;
         }  
       }
       if (!found) {
         std::pair<int,int> pair1(part.pdgId(),1);
         genrhadcounts.push_back(pair1);
         //cout << " NEW ID = " << part.pdgId() << " " << endl;
       }
     }
   }

// Get R-Hadron properties from the tracker
  edm::PSimHitContainer::const_iterator itHit;

  int icount = 0;
  int nhitstib[20] = {0,0,0,0,0,0,0,0,0,0};
  int nhitstib1[20] = {0,0,0,0,0,0,0,0,0,0};
  int nhitstib2[20] = {0,0,0,0,0,0,0,0,0,0};

  for (itHit = SiTIBLowContainer->begin(); itHit != SiTIBLowContainer->end(); ++itHit) {
    DetId detid = DetId(itHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(itHit->localPosition());
    float r = sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y());

    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<999999999)) {
      float dphi1 = genrhad1->phi() - gpos.phi();
      float dphi2 = genrhad2->phi() - gpos.phi();
      float deta1 = genrhad1->eta() - gpos.eta();
      float deta2 = genrhad2->eta() - gpos.eta();
      float dR1 = sqrt(dphi1*dphi1+deta1*deta1);
      float dR2 = sqrt(dphi2*dphi2+deta2*deta2);

      // Fill R-Hadron momentum histograms
      RHadron1_px->Fill(genrhad1->p4().px());
      RHadron1_py->Fill(genrhad1->p4().py());
      RHadron1_pz->Fill(genrhad1->p4().pz());
      RHadron2_px->Fill(genrhad2->p4().px());
      RHadron2_py->Fill(genrhad2->p4().py());
      RHadron2_pz->Fill(genrhad2->p4().pz());
      
      const reco::GenParticle *genrhadclosest=0;
      int iclosest = -1;
      float dRclosest = 999999.;
      if (dR1<dR2) {
        //htibdeltaRclosest->Fill(dR1);
        dRclosest = dR1;
        iclosest = 1;
        genrhadclosest = genrhad1;
      } else {
        //htibdeltaRclosest->Fill(dR2);
        dRclosest = dR2;
        iclosest = 2;
        genrhadclosest = genrhad2;
      }
      int ilayer = -1;
      if (dRclosest<1.0) {
        if ((r>23)&&(r<29)) { // first layer
          ilayer = TIB0 + 1;
        } else if ((r>31)&&(r<36)) { // second layer
          ilayer = TIB0 + 2;
        } else if ((r>39)&&(r<45)) { // third layer
          ilayer = TIB0 + 3;
        } else if ((r>47)&&(r<53)) { // fourth layer
          ilayer = TIB0 + 4;
        } else {
          //cout << "didn't find layer " << r << " " << genrhadclosest->pdgId() << endl;
        }
      }
      if (ilayer>-1) {
        nhitstib[ilayer]++;
        if (dR1<0.4) nhitstib1[ilayer]++;
        if (dR2<0.4) nhitstib2[ilayer]++;
        if (iclosest==1) {
          linkedhits1.push_back({ilayer,itHit->particleType()});
        } else if (iclosest==2) {
          linkedhits2.push_back({ilayer,itHit->particleType()});
        }
      }
    }

  // Grab simulated calorimiter hits
  edm::PCaloHitContainer::const_iterator itHit;
  for (itHit = EcalEBContainer->begin(); itHit != EcalEBContainer->end(); ++itHit) {
      DetId detid = DetId(itHit->detUnitId());
      const GeomDetUnit *det = (const GeomDetUnit *)caloGeometry->idToDetUnit(detid);
      GlobalPoint gpos = det->toGlobal(itHit->localPosition());
      float r = sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y());

      // Fill calohit location histograms
      EBHit_ieta->Fill(gpos.eta());
      EBHit_iphi->Fill(gpos.phi());
      EBHit_x->Fill(gpos.x());
      EBHit_y->Fill(gpos.y());
      EBHit_z->Fill(gpos.z());
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimCaloHitAnalyzer);
