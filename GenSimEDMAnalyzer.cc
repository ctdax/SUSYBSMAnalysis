// system include files
#include <memory>
#include <cmath>

// user include files
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

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTChamber.h"
#include "Geometry/DTGeometry/interface/DTSuperLayer.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
 
#include "DataFormats/EcalDetId/interface/EBDetId.h"
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

#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
//#include "SUSYBSMAnalysis/Analyzer/interface/CommonFunction.h"
//#include "SUSYBSMAnalysis/Analyzer/interface/DeDxUtility.h"

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

#define FWCORE
#include "SUSYBSMAnalysis/Analyzer/interface/CommonFunction.h"
#include "SUSYBSMAnalysis/Analyzer/interface/DeDxUtility.h"
#include "SUSYBSMAnalysis/HSCP/interface/HSCPHelpers.h"
////#include "SUSYBSMAnalysis/Analyzer/interface/TOFUtility.h"
//#include "SUSYBSMAnalysis/Analyzer/interface/TupleMaker.h"

// Thresholds for candidate preselection
  float globalMaxEta_ = 1.0;
  float globalMinPt_ = 55.0;
  unsigned int globalMinNOPH_ = 1;
  float globalMinFOVH_ = 0.8;
  unsigned int globalMinNOM_ = 10;
  float globalMaxChi2_ = 5.0;
  float globalMaxEoP_ = 0.3;
  float globalMaxDZ_ = 0.1;
  float globalMaxDXY_ = 0.02;
  float globalMaxTIsol_, globalMiniRelIsoAll_;
  float globalMinIh_ = 3.47;
  float trackProbQCut_;
  unsigned int minMuStations_;
  float globalMinIs_ = 0.0;
  float globalMinTOF_;
  float GlobalMinNDOF = 8;            // cut on number of     DegreeOfFreedom used for muon TOF measurement
  float GlobalMinNDOFDT = 6;          // cut on number of DT  DegreeOfFreedom used for muon TOF measurement
  float GlobalMinNDOFCSC = 6;         // cut on number of CSC DegreeOfFreedom used for muon TOF measurement
  float GlobalMaxTOFErr = 0.15;       //0.07;   // cut on error on muon TOF measurement
  bool useClusterCleaning = true;

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

class GenSimEDMAnalyzer : public edm::EDAnalyzer {
public:
  explicit GenSimEDMAnalyzer (const edm::ParameterSet&);
  ~GenSimEDMAnalyzer();


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

  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_muonCSC_Token_;
  edm::EDGetTokenT <edm::PSimHitContainer> edmPSimHitContainer_muonDT_Token_;

  edm::EDGetTokenT<edm::SimTrackContainer> edmSimTrackContainerToken_;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerSummary_;
  edm::EDGetTokenT<std::vector<reco::PFMET>> pfmet_;
  edm::EDGetTokenT<std::vector<reco::CaloMET>> calomet_;
  //edm::EDGetTokenT<std::vector<reco::DeDxHitInfo>> dedxhitinfo_;
//TA  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
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
//  TH1F  *hltCalomet, *hltCalomet_muon, *hltCalomet_met120, *hltCalomet_met170;
  TH1F *hdeltaphicmet;
  TH1F *hdeltaphicmetgt100;
  TH1F *hdeltaphicmetlt100;
  TH2F *hdeltaphicmet2D;
  TH2F *hdeltaphicmet2Db;
  TH1F *hdeltaphimet;
  TH1F *hdeltaphimetgt100;
  TH1F *hdeltaphimetlt100;
  TH2F *hdeltaphimet2D;
  TH2F *hdeltaphimet2Db;
  TH1F *hdeltaphimetq0q0, *hdeltaphimetq0q1, *hdeltaphimetq1q1;
  TH1F *hdeltaphimetq0q2, *hdeltaphimetq1q2, *hdeltaphimetq2q2;
  TH1F *hdeltaphigen12;
  TH1F *hdRgenpf, *hdRgentrk, *hdRtrkpf;
  TH1F *hdRgenpf00, *hdRgenpf01, *hdRgenpf01n, *hdRgenpf01c, *hdRgenpf11;
//  TH1F *hdRgentrk00, *hdRgentrk01, *hdRgentrk01n, *hdRgentrk01c, *hdRgentrk11;
  TH1F *hdphigenpf, *hdphigentrk, *hdphitrkpf;
  TH1F *hdphigenpf00, *hdphigenpf01, *hdphigenpf01n, *hdphigenpf01c, *hdphigenpf11;
  TH1F *hdphigentrk00, *hdphigentrk01, *hdphigentrk01n, *hdphigentrk01c, *hdphigentrk11;
  TH1F *hdptgenpf, *hdptgentrk, *hdpttrkpf;
  TH1F *hdptgenpf00, *hdptgenpf01, *hdptgenpf01n, *hdptgenpf01c, *hdptgenpf11;
  TH1F *hdptgentrk00, *hdptgentrk01, *hdptgentrk01n, *hdptgentrk01c, *hdptgentrk11;
  TH1F *hdptgenpffrac, *hdptgentrkfrac, *hdpttrkpffrac;
  TH1F *hdptgenpffrac00, *hdptgenpffrac01, *hdptgenpffrac01n, *hdptgenpffrac01c, *hdptgenpffrac11;
//  TH1F *hdptgentrk00, *hdptgentrk01, *hdptgentrk01n, *hdptgentrk01c, *hdptgentrk11;
  TH1F *hdpttrkpf00, *hdpttrkpf01, *hdpttrkpf01n, *hdpttrkpf01c, *hdpttrkpf11;
  TH2F *h2pfptvsgenpt, *h2tkptvsgenpt, *h2tkptvspfpt;
  TH2F *h2pfptvsgenpt00, *h2pfptvsgenpt01, *h2pfptvsgenpt01n, *h2pfptvsgenpt01c, *h2pfptvsgenpt11;
  TH2F *h2tkptvsgenpt00, *h2tkptvsgenpt01, *h2tkptvsgenpt01n, *h2tkptvsgenpt01c, *h2tkptvsgenpt11;
  TH2F *h2tkptvspfpt00, *h2tkptvspfpt01, *h2tkptvspfpt01n, *h2tkptvspfpt01c, *h2tkptvspfpt11;

  TH1F *hgeneta, *hgenetaM;
  TH1F *hgenbeta, *hgenbetach, *hgenbetaneut;
  TH2F *h2genbeta;
  TH1F *hgendphi12;
  TH1F *hgendR12;

  TH1F *htibrdist,*htobrdist,*hbpixrdist,*hdtrdist,*hcscrdist;
  TH1F *htibdeltaphi,*htibdeltaphi2,*htibdeltaeta,*htibdeltaR,*htibdeltaRclosest;
  TH1F *htobdeltaphi,*htobdeltaphi2,*htobdeltaeta,*htobdeltaR,*htobdeltaRclosest;
  TH1F *hbpixdeltaphi,*hbpixdeltaphi2,*hbpixdeltaeta,*hbpixdeltaR,*hbpixdeltaRclosest;
  TH1F *hdtdeltaphi,*hdtdeltaphi2,*hdtdeltaeta,*hdtdeltaR,*hdtdeltaRclosest;
  TH2F *htracker0,*htracker1,*htracker2,*htracker3;
  TH2F *htib0,*htib1,*htib2,*htib3;
  TH2F *hbpix0,*hbpix1,*hbpix2,*hbpix3;
  TH2F *hrzbpix0,*hrzbpix1,*hrzbpix2,*hrzbpix3;

  TH1F *hlastlayersame, *hlastlayercharged;
  TH1F *hlastlayersame2, *hlastlayercharged2;
  TH1F *hhitlayersame, *hhitlayercharged;
  TH1F *hhitlayersame2, *hhitlayercharged2;

  //TH1F  *dedxnumbhits;

};

int GenSimEDMAnalyzer::determineCharge(int pdgID) {
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
  else {cout << "Unknown PID = " << pid << endl;}
  if ((charge>-98)&&(pdgID<0)) charge *= -1;
  return(charge);
}

int GenSimEDMAnalyzer::determineType(int pdgID) {
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
float GenSimEDMAnalyzer::deltaphi(float phi1, float phi2) {
  float tdelta = phi1 - phi2;
  while (tdelta<(-1.0*M_PI)) {tdelta += (2.0*M_PI);}
  while (tdelta>(1.0*M_PI)) {tdelta -= (2.0*M_PI);}
//  cout << " tdelta = " << tdelta << " phi1 = " << phi1 << " phi2 = " << phi2 << endl;
  return(tdelta);
}

// |delta(phi1-phi2)| - return value from 0 to +pi
float GenSimEDMAnalyzer::deltaphiabs(float phi1, float phi2) {
  float tdphi = deltaphi(phi1,phi2);
  return(fabs(tdphi));
}

float GenSimEDMAnalyzer::deltaRta(float phi1, float eta1, float phi2, float eta2) {
  float diffphi1 = deltaphiabs(phi1,phi2);
  float diffeta1 = eta1 - eta2;
  float tmpdeltaR1 = sqrt(diffphi1*diffphi1 + diffeta1*diffeta1);
  return(tmpdeltaR1);
}

// preselection
bool GenSimEDMAnalyzer::passHSCPPreselection(int typeMode, susybsm::HSCParticle hscpCan, reco::VertexCollection vertexColl, reco::VertexCollection inclusiveSecondaryVertices, reco::PFCandidateCollection pf, reco::DeDxData *dedxSObj, reco::DeDxData *dedxMObj, bool passedCutsArray[]) {

    //
//  reco::TrackRef track = (typeMode != 3) ? hscpCan.trackRef() : track = muon->standAloneMuon();
  reco::TrackRef track = hscpCan.trackRef();

  cout << "inside preselection" << endl;

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

//    bool pf_isMuon = false, pf_isElectron = false, pf_isChHadron = false, pf_isNeutHadron = false;
//    int pf_muon_selector = -1;
//    float pf_ecal_energy = 0, pf_hcal_energy = 0;

//    if(pfCandHandle.isValid() && !pfCandHandle->empty()) {
//      const reco::PFCandidateCollection* pf = pfCandHandle.product();
//      for (unsigned int i = 0; i < pf->size(); i++){
      for (unsigned int i = 0; i < pf.size(); i++){
//          const reco::PFCandidate* pfCand = &(*pf)[i];
          const reco::PFCandidate pfCand = pf[i];
//          float dr = deltaR(pfCand->eta(),pfCand->phi(),track->eta(),track->phi());
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
//        const reco::PFCandidate* pfCand = &(*pf)[i];
        const reco::PFCandidate* pfCand = &(pf)[i];
//        if(i == idx_pf_RMin) {
//            pf_isMuon = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::mu;
//            pf_isElectron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::e;
//            pf_isChHadron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h;
//            pf_isNeutHadron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h0;
//            pf_ecal_energy = pfCand->ecalEnergy();
//            pf_hcal_energy = pfCand->hcalEnergy();
//        }

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
//   unsigned int missingHitsTillLast =
//      track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) +
//      track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
//   float validFractionTillLast =
//            track->found() <= 0 ? -1 : track->found() / float(track->found() + missingHitsTillLast);

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

 // Cut on the uncertainty of the pt measurement
//Temp  passedCutsArray[11] = (typeMode != 3 && (track->ptError() / track->pt()) < pTerr_over_pT_etaBin(track->pt(), track->eta())) ? true : false;

 // Cut on the tracker based isolation
  passedCutsArray[12] = (!isMaterialTrack) ? true : false;
//  passedCutsArray[13] = ( IsoTK_SumEt < globalMaxTIsol_) ? true : false;
 // Cut on the PF based mini-isolation
  passedCutsArray[13] = ( miniRelIsoAll < globalMiniRelIsoAll_) ? true : false;
 // Cut on the PF electron ID
//TEMP  passedCutsArray[14] = ( !pf_isElectron  && !pf_isPhoton) ? true : false;
 // Cut on min Ih (or max for fractionally charged)
  passedCutsArray[15] = (  (typeMode != 5 &&  Ih > globalMinIh_)
                         || (typeMode == 5 && Ih < globalMinIh_)) ? true : false;
 // Cut away background events based on the probXY
//TEMP  passedCutsArray[16] = ((probXYonTrackNoLayer1 > 0.1 && probXYonTrackNoLayer1 < 1.0))  ? true : false;
 // Cut away background events based on the probQ
//TEMP  passedCutsArray[17] = (probQonTrackNoLayer1 < trackProbQCut_) ? true : false;
//passedCutsArray[17]  = (probQonTrack < trackProbQCut_ || probQonTrackNoLayer1 < trackProbQCut_) ? true : false;
 // TOF only cuts
//TEMP  passedCutsArray[18] = (typeMode != 3 || (typeMode == 3 && muonStations(track->hitPattern()) > minMuStations_)) ? true : false;
//TEMP  passedCutsArray[19] = (typeMode != 3 || (typeMode == 3 && fabs(track->phi()) > 1.2 && fabs(track->phi()) < 1.9)) ? true : false;
//TEMP  passedCutsArray[20] = (typeMode != 3 || (typeMode == 3 && fabs(minEta) > minSegEtaSep)) ? true : false;
  
 
  return(true);
}

bool GenSimEDMAnalyzer::comparePID(pair<int, int> p1, pair<int, int> p2) {
  return (p1.first < p2.first);
}

//constructor
GenSimEDMAnalyzer::GenSimEDMAnalyzer(const edm::ParameterSet& iConfig) 
 :  tGeomEsToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>())
  //triggerBits_(consumes(iConfig.getParameter("bits"))),
  //triggerObjects_(consumes(iConfig.getParameter("objects"))),
  //triggerPrescales_(consumes(iConfig.getParameter("prescales")))
{

  
  edmPSimHitContainer_siTIBLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("SiTIBLowSrc"));
  edmPSimHitContainer_siTIBHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("SiTIBHighSrc"));
  edmPSimHitContainer_siTOBLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("SiTOBLowSrc"));
  edmPSimHitContainer_siTOBHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("SiTOBHighSrc"));
  edmPSimHitContainer_siTIDLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("SiTIDLowSrc"));
  edmPSimHitContainer_siTIDHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("SiTIDHighSrc"));
  edmPSimHitContainer_siTECLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("SiTECLowSrc"));
  edmPSimHitContainer_siTECHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("SiTECHighSrc"));
  edmPSimHitContainer_pxlBrlLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("PxlBrlLowSrc"));
  edmPSimHitContainer_pxlBrlHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("PxlBrlHighSrc"));
  edmPSimHitContainer_pxlFwdLow_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("PxlFwdLowSrc"));
  edmPSimHitContainer_pxlFwdHigh_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("PxlFwdHighSrc"));
  edmPSimHitContainer_muonCSC_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("MuonCSCSrc"));
  edmPSimHitContainer_muonDT_Token_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("MuonDTSrc"));

  edmSimTrackContainerToken_ = consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("G4TrkSrc"));

  //triggerBits_ = consumes<edm::TriggerResults>(InputTag("TriggerResults","","HLT"));
  //triggerObjects_ = consumes<edm::TriggerResults>(InputTag("TriggerResults","","HLT"));
  //triggerPrescales_ = consumes<edm::TriggerResults>(InputTag("TriggerResults","","HLT"));

  triggerBits_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"));
  //triggerObjects_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("objects"));
  //triggerPrescales_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("prescales"));
  triggerSummary_ = consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("trig_sum"));
  pfmet_ = consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("pfmet_reco"));
  calomet_ = consumes<std::vector<reco::CaloMET>>(iConfig.getParameter<edm::InputTag>("calomet_reco"));
  //dedxhitinfo_ = consumes<std::vector<reco::DeDxHitInfo>>(iConfig.getParameter<edm::InputTag>("dedx_info"));

//TA   genParticlesToken_ = mayConsume<reco::GenParticleCollection>(InputTag("gen_info"));
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

  bool splitByModuleType = true;
  //dEdxTemplates = loadDeDxTemplate(dEdxTemplate_, splitByModuleType);
  //if (enableDeDxCalibration_)
  //  trackerCorrector.LoadDeDxCalibration(dEdxCalibration_);
  //else
  //  trackerCorrector.TrackerGains = nullptr;

  //outputFile_ = new TFile("/uscms/home/rkim/nobackup/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/met_histo_StopM800.root", "RECREATE" );  
  //outputFile_ = new TFile("/uscms/home/rkim/nobackup/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/met_histo_StopM800N.root", "RECREATE" );  
  //outputFile_ = new TFile("/uscms/home/rkim/nobackup/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/met_histo_GMStauM871.root", "RECREATE" );  
  //outputFile_ = new TFile("/uscms/home/rkim/nobackup/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/met_histo_GluinoM1800.root", "RECREATE" );  
  //outputFile_ = new TFile("/uscms/home/rkim/nobackup/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/met_histo_GluinoM1800N.root", "RECREATE" );  
  //outputFile_ = new TFile("/uscms/home/rkim/nobackup/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/met_histo_StopM800.root", "RECREATE" );  
  //outputFile_ = new TFile("/uscms/home/rkim/nobackup/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/met_histo_StopM800.root", "RECREATE" );  

//  outputFile_ = new TFile("/uscms/home/rkim/nobackup/HSCP_triggers/CMSSW_10_6_27/src/SUSYBSMAnalysis/HSCP/test/met_histo_GluinoM1400.root", "RECREATE" );  
//  outputFile_ = new TFile("/uscms_data/d2/tadams/hscp/CMSSW_10_6_27/src/SUSYBSMAnalysis/HSCP/test/sim_histo_gluinoM1800.root", "RECREATE" );  
  outputFile_ = new TFile("/uscms/home/cthompso/nobackup/CMSSW_10_6_30/src/HSCP/TrackerMuonHitAnalysis/test/test.root", "RECREATE");  


  hdRgenpf = new TH1F("hdRgenpf","delta(R) GEN - PF; delta(R)",120,0,0.02);
  hdRgenpf00 = new TH1F("hdRgenpf00","delta(R) GEN - PF q=00; delta(R)",120,0,0.02);
  hdRgenpf01 = new TH1F("hdRgenpf01","delta(R) GEN - PF q=01; delta(R)",120,0,0.02);
  hdRgenpf01n = new TH1F("hdRgenpf01n","delta(R) GEN - PF q=01n; delta(R)",120,0,0.02);
  hdRgenpf01c = new TH1F("hdRgenpf01c","delta(R) GEN - PF q=01c; delta(R)",120,0,0.02);
  hdRgenpf11 = new TH1F("hdRgenpf11","delta(R) GEN - PF q=11; delta(R)]",120,0,0.02);

  //1500 bins for TEfficiency plots
  //500 bins for MET distribution plots

/*
  hltCalomet = new TH1F("hltCalomet","HLT Calo MET; HLT Calo MET [GeV]",1500,0,20000);
  hltCalomet_muon = new TH1F("hltCalomet_muon","HLT Calo MET passing muon; HLT Calo MET [GeV]",1500,0,20000);
  hltCalomet_met120 = new TH1F("hltCalomet_met120","HLT Calo MET passing MET120; HLT Calo MET [GeV]",1500,0,20000);
  hltCalomet_met170 = new TH1F("hltCalomet_met170","HLT Calo MET passing MET170; HLT Calo MET [GeV]",1500,0,20000);
*/

  hdeltaphicmet = new TH1F("hdeltaphicmet","delta(phi) CALO MET - GEN; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphicmetgt100 = new TH1F("hdeltaphicmetgt100","delta(phi) CALO MET - GEN; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphicmetlt100 = new TH1F("hdeltaphicmetlt100","delta(phi) CALO MET - GEN; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphicmet2D = new TH2F("hdeltaphicmet2D","delta(phi) CALO MET - GEN; delta(phi) [rad]; delta(phi) [rad]",40,-0.5,3.5,40,-0.5,3.5);
  hdeltaphicmet2Db = new TH2F("hdeltaphicmet2Db","delta(phi) vs CALO MET; delta(phi) [rad]; CALO MET [GeV]",40,-0.5,3.5,40,0,2000);
  hdeltaphimet = new TH1F("hdeltaphimet","delta(phi) RECO PF MET - GEN; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphimetgt100 = new TH1F("hdeltaphimetgt100","delta(phi) RECO PF MET - GEN; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphimetlt100 = new TH1F("hdeltaphimetlt100","delta(phi) RECO PF MET - GEN; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphimet2D = new TH2F("hdeltaphimet2D","delta(phi) RECO PF MET - GEN; delta(phi) [rad]; delta(phi) [rad]",40,-0.5,3.5,40,-0.5,3.5);
  hdeltaphimet2Db = new TH2F("hdeltaphimet2Db","delta(phi) vs PF MET; delta(phi) [rad]; PF MET [GeV]",40,-0.5,3.5,40,0,2000);
  hdeltaphimetq0q0 = new TH1F("hdeltaphimetq0q0","delta(phi) RECO PF MET - GEN q0q0; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphimetq0q1 = new TH1F("hdeltaphimetq0q1","delta(phi) RECO PF MET - GEN q0q1; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphimetq1q1 = new TH1F("hdeltaphimetq1q1","delta(phi) RECO PF MET - GEN q1q1; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphimetq0q2 = new TH1F("hdeltaphimetq0q2","delta(phi) RECO PF MET - GEN q0q2; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphimetq1q2 = new TH1F("hdeltaphimetq1q2","delta(phi) RECO PF MET - GEN q1q2; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphimetq2q2 = new TH1F("hdeltaphimetq2q2","delta(phi) RECO PF MET - GEN q2q2; delta(phi) [rad]",70,-3.5,3.5);
  hdeltaphigen12 = new TH1F("hdeltaphigen12","delta(phi) GEN1 - GEN2; delta(phi) [rad]",70,-3.5,3.5);

//  hdRgentrk = new TH1F("hdRgentrk","delta(R) GEN - trk; delta(R)",120,0,0.1);
//  hdRgentrk00 = new TH1F("hdRgentrk00","delta(R) GEN - trk q=00; delta(R)",120,0,0.1);
//  hdRgentrk01 = new TH1F("hdRgentrk01","delta(R) GEN - trk q=01; delta(R)",120,0,0.1);
//  hdRgentrk01n = new TH1F("hdRgentrk01n","delta(R) GEN - trk q=01n; delta(R)",120,0,0.1);
//  hdRgentrk01c = new TH1F("hdRgentrk01c","delta(R) GEN - trkq=01c; delta(R)",120,0,0.1);
//  hdRgentrk11 = new TH1F("hdRgentrk11","delta(R) GEN - trk q=11; delta(R)",120,0,0.1);
  hdphigenpf = new TH1F("hdphigenpf","delta(phi) GEN - PF; delta(phi) [rad]",200,-0.06,0.06);
  hdphigenpf00 = new TH1F("hdphigenpf00","delta(phi) GEN - PF q=00; delta(phi) [rad]",200,-0.06,0.06);
  hdphigenpf01 = new TH1F("hdphigenpf01","delta(phi) GEN - PF q=01; delta(phi) [rad]",200,-0.06,0.06);
  hdphigenpf01n = new TH1F("hdphigenpf01n","delta(phi) GEN - PF q=01n; delta(phi) [rad]",200,-0.06,0.06);
  hdphigenpf01c = new TH1F("hdphigenpf01c","delta(phi) GEN - PF q=01c; delta(phi) [rad]",200,-0.06,0.06);
  hdphigenpf11 = new TH1F("hdphigenpf11","delta(phi) GEN - PF q=11; delta(phi) [rad]",200,-0.06,0.06);
  hdphigentrk = new TH1F("hdphigentrk","delta(phi) GEN - gen trk; delta(phi) [rad]",200,-0.06,0.06);
  hdphigentrk00 = new TH1F("hdphigentrk00","delta(phi) GEN - trk q=00; delta(phi) [rad]",200,-0.06,0.06);
  hdphigentrk01 = new TH1F("hdphigentrk01","delta(phi) GEN - trk q=01; delta(phi) [rad]",200,-0.06,0.06);
  hdphigentrk01n = new TH1F("hdphigentrk01n","delta(phi) GEN - trk q=01n; delta(phi) [rad]",200,-0.06,0.06);
  hdphigentrk01c = new TH1F("hdphigentrk01c","delta(phi) GEN - trk q=01c; delta(phi) [rad]",200,-0.06,0.06);
  hdphigentrk11 = new TH1F("hdphigentrk11","delta(phi) GEN - trk q=11; delta(phi) [rad]",200,-0.06,0.06);
  hdphitrkpf = new TH1F("hdphitrkpf","delta(phi) gen trk - trk; delta(phi) [rad]",200,-0.06,0.06);
  hdptgenpf = new TH1F("hdptgenpf","delta(pT) PF - GEN; delta(pT) [GeV]",100,-500.,500);
  hdptgenpf00 = new TH1F("hdptgenpf00","delta(pt) PF - GEN q=00; delta(pt) [GeV]",100,-500.,500.);
  hdptgenpf01 = new TH1F("hdptgenpf01","delta(pt) PF - GEN q=01; delta(pt) [GeV]",100,-500.,500.);
  hdptgenpf01n = new TH1F("hdptgenpf01n","delta(pt) PF - GEN q=01n; delta(pt) [GeV]",100,-500.,500.);
  hdptgenpf01c = new TH1F("hdptgenpf01c","delta(pt) PF - GEN q=01c; delta(pt) [GeV]",100,-500.,500.);
  hdptgenpf11 = new TH1F("hdptgenpf11","delta(pt) PF - GEN q=11; delta(pt) [GeV]",100,-500.,500.);
  hdptgenpffrac = new TH1F("hdptgenpffrac","(delta(pT) PF - GEN)/pT GEN; delta(pT)/pT [GeV]",100,-2.,5);
  hdptgenpffrac00 = new TH1F("hdptgenpffrac00","(delta(pt) PF - GEN)/pT GEN q=00; delta(pt)/pT [GeV]",100,-2.,5.);
  hdptgenpffrac01 = new TH1F("hdptgenpffrac01","(delta(pt) PF - GEN)/pT GEN q=01; delta(pt)/pT [GeV]",100,-2.,5.);
  hdptgenpffrac01n = new TH1F("hdptgenpffrac01n","(delta(pt) PF - GEN)/pT GEN q=01n; delta(pt)/pT [GeV]",100,-2.,5.);
  hdptgenpffrac01c = new TH1F("hdptgenpffrac01c","(delta(pt) PF - GEN)/pT GEN q=01c; delta(pt)/pT [GeV]",100,-2.,5.);
  hdptgenpffrac11 = new TH1F("hdptgenpffrac11","(delta(pt) PF - GEN)/pT GEN q=11; delta(pt)/pT [GeV]",100,-2.,5.);
  hdptgentrk = new TH1F("hdptgentrk","delta(pT) track - GEN; delta(pT) [GeV]",100,-500.,500);
  hdptgentrk00 = new TH1F("hdptgentrk00","delta(pt) track - GEN q=00; delta(pt) [GeV]",100,-500.,500.);
  hdptgentrk01 = new TH1F("hdptgentrk01","delta(pt) track - GEN q=01; delta(pt) [GeV]",100,-500.,500.);
  hdptgentrk01n = new TH1F("hdptgentrk01n","delta(pt) track - GEN q=01n; delta(pt) [GeV]",100,-500.,500.);
  hdptgentrk01c = new TH1F("hdptgentrk01c","delta(pt) track - GEN q=01c; delta(pt) [GeV]",100,-500.,500.);
  hdptgentrk11 = new TH1F("hdptgentrk11","delta(pt) track - GEN q=11; delta(pt) [GeV]",100,-500.,500.);
  hdpttrkpf = new TH1F("hdpttrkpf","delta(pT) PF - track; delta(pT) [GeV]",100,-500.,500);
  hdpttrkpf00 = new TH1F("hdpttrkpf00","delta(pt) PF - track q=00; delta(pt) [GeV]",100,-500.,500.);
  hdpttrkpf01 = new TH1F("hdpttrkpf01","delta(pt) PF - track q=01; delta(pt) [GeV]",100,-500.,500.);
  hdpttrkpf01n = new TH1F("hdpttrkpf01n","delta(pt) PF - track q=01n; delta(pt) [GeV]",100,-500.,500.);
  hdpttrkpf01c = new TH1F("hdpttrkpf01c","delta(pt) PF - track q=01c; delta(pt) [GeV]",100,-500.,500.);
  hdpttrkpf11 = new TH1F("hdpttrkpf11","delta(pt) PF - track q=11; delta(pt) [GeV]",100,-500.,500.);
  h2pfptvsgenpt = new TH2F("h2pgptvsgenpt","PF pT vs. GEN pT; gen particle pT [GeV]; PF particle pT [GeV]",200,0,4000,200,0,4000);
  h2pfptvsgenpt00 = new TH2F("h2pgptvsgenpt00","PF pT vs. GEN pT q=00; gen particle pT [GeV]; PF particle pT [GeV]",200,0,4000,200,0,4000);
  h2pfptvsgenpt01 = new TH2F("h2pgptvsgenpt01","PF pT vs. GEN pT q=01; gen particle pT [GeV]; PF particle pT [GeV]",200,0,4000,200,0,4000);
  h2pfptvsgenpt01n = new TH2F("h2pgptvsgenpt01n","PF pT vs. GEN pT q=01n; gen particle pT [GeV]; PF particle pT [GeV]",200,0,4000,200,0,4000);
  h2pfptvsgenpt01c = new TH2F("h2pgptvsgenpt01c","PF pT vs. GEN pT q=01c; gen particle pT [GeV]; PF particle pT [GeV]",200,0,4000,200,0,4000);
  h2pfptvsgenpt11 = new TH2F("h2pgptvsgenpt11","PF pT vs. GEN pT q=11; gen particle pT [GeV]; PF particle pT [GeV]",200,0,4000,200,0,4000);
  h2tkptvsgenpt = new TH2F("h2tkptvsgenpt","track pT vs. GEN pT; gen particle pT [GeV]; general track pT [GeV]",200,0,4000,200,0,4000);
  h2tkptvspfpt = new TH2F("h2tkptvspfpt","track pT vs. PF particle pT; PF particle pT [GeV]; general track pT [GeV]",200,0,4000,200,0,4000);
  hgeneta = new TH1F("hgeneta","generator eta; eta",100,-5.,5);
  hgenetaM = new TH1F("hgenetaM","generator eta - matchA; eta",100,-5.,5);

  hgenbeta = new TH1F("hgenbeta","R-hadron generator beta; beta",50,0.0,1.0);
  h2genbeta = new TH2F("h2genbeta","R-hadron generator beta (high) vs beta (low);beta (high);beta (low)",100,0.0,1.0,100,0.0,1.0);
  hgenbetach = new TH1F("hgenbetach","R-hadron (Q!=0) generator beta; beta",50,0.0,1.0);
  hgenbetaneut = new TH1F("hgenbetaneut","R-hadron (Q==0) generator beta; beta",50,0.0,1.0);
  hgendphi12 =  new TH1F("hgendphi12","R-hadron generator |phi1-phi2|; |delta(phi1-phi2)|",90,-0.5,4.0);
  hgendR12 =  new TH1F("hgendR12","R-hadron generator |Delta(R)|; |delta(R)|",160,0.,8.0);

  htibrdist = new TH1F("htibrdist","hit radius (TIB); r (units?)",160,20,60.0);
  htobrdist = new TH1F("htobrdist","hit radius (TOB); r (units?)",240,55,115.0);
  hbpixrdist = new TH1F("hbpixrdist","hit radius (pixel barrel); r (units?)",160,0,20.0);
  hdtrdist = new TH1F("hdtrdist","hit radius (DT); r (units?)",500,350.0,850.0);
  hcscrdist = new TH1F("hcscrdist","hit radius (CSC); r (units?)",160,150.0,3000.0);
  htibdeltaphi = new TH1F("htibdeltaphi","delta(phi) GEN - simhit(TIB); delta(phi)",320,-4,4.0);
  htibdeltaphi2 = new TH1F("htibdeltaphi2","delta(phi) GEN - simhit(TIB); delta(phi)",200,-1,1.0);
  htibdeltaeta = new TH1F("htibdeltaeta","delta(eta) GEN - simhit(TIB); delta(eta)",320,-4,4.0);
  htibdeltaR = new TH1F("htibdeltaR","delta(R) GEN - simhit(TIB); delta(R)",200,-0.5,4.5);
  htibdeltaRclosest = new TH1F("htibdeltaRclosest","min(delta(R) GEN - simhit(TIB)); delta(R)",200,-0.5,4.5);
  htobdeltaphi = new TH1F("htobdeltaphi","delta(phi) GEN - simhit(TOB); delta(phi)",320,-4,4.0);
  htobdeltaphi2 = new TH1F("htobdeltaphi2","delta(phi) GEN - simhit(TOB); delta(phi)",200,-1,1.0);
  htobdeltaeta = new TH1F("htobdeltaeta","delta(eta) GEN - simhit(TOB); delta(eta)",320,-4,4.0);
  htobdeltaR = new TH1F("htobdeltaR","delta(R) GEN - simhit(TOB); delta(R)",200,-0.5,4.5);
  htobdeltaRclosest = new TH1F("htobdeltaRclosest","min(delta(R) GEN - simhit(TOB)); delta(R)",200,-0.5,4.5);
  hbpixdeltaphi = new TH1F("hbpixdeltaphi","delta(phi) GEN - simhit(BPIX); delta(phi)",320,-4,4.0);
  hbpixdeltaphi2 = new TH1F("hbpixdeltaphi2","delta(phi) GEN - simhit(BPIX); delta(phi)",200,-1,1.0);
  hbpixdeltaeta = new TH1F("hbpixdeltaeta","delta(eta) GEN - simhit(BPIX); delta(eta)",320,-4,4.0);
  hbpixdeltaR = new TH1F("hbpixdeltaR","delta(R) GEN - simhit(BPIX); delta(R)",200,-0.5,4.5);
  hbpixdeltaRclosest = new TH1F("hbpixdeltaRclosest","min(delta(R) GEN - simhit(BPIX)); delta(R)",200,-0.5,4.5);
  hdtdeltaphi = new TH1F("hdtdeltaphi","delta(phi) GEN - simhit(DT); delta(phi)",1280,-4,4.0);
  hdtdeltaphi2 = new TH1F("hdtdeltaphi2","delta(phi) GEN - simhit(DT); delta(phi)",200,-1,1.0);
  hdtdeltaeta = new TH1F("hdtdeltaeta","delta(eta) GEN - simhit(DT); delta(eta)",1280,-4,4.0);
  hdtdeltaR = new TH1F("hdtdeltaR","delta(R) GEN - simhit(DT); delta(R)",200,-0.5,4.5);
  hdtdeltaRclosest = new TH1F("hdtdeltaRclosest","min(delta(R) GEN - simhit(DT)); delta(R)",200,-0.5,4.5);
  htracker0 = new TH2F("htracker0","event display;y;x",480,-120.0,120.0,480,-120.0,120.0);
  htracker1 = new TH2F("htracker1","event display;y;x",480,-120.0,120.0,480,-120.0,120.0);
  htracker2 = new TH2F("htracker2","event display;y;x",480,-120.0,120.0,480,-120.0,120.0);
  htracker3 = new TH2F("htracker3","event display;y;x",480,-120.0,120.0,480,-120.0,120.0);
  htib0 = new TH2F("htib0","event display;y;x",240,-60.0,60.0,240,-60.0,60.0);
  htib1 = new TH2F("htib1","event display;y;x",240,-60.0,60.0,240,-60.0,60.0);
  htib2 = new TH2F("htib2","event display;y;x",240,-60.0,60.0,240,-60.0,60.0);
  htib3 = new TH2F("htib3","event display;y;x",240,-60.0,60.0,240,-60.0,60.0);
  hbpix0 = new TH2F("hbpix0","event display;y;x",160,-20.0,20.0,160,-20.0,20.0);
  hbpix1 = new TH2F("hbpix1","event display;y;x",160,-20.0,20.0,160,-20.0,20.0);
  hbpix2 = new TH2F("hbpix2","event display;y;x",160,-20.0,20.0,160,-20.0,20.0);
  hbpix3 = new TH2F("hbpix3","event display;y;x",160,-20.0,20.0,160,-20.0,20.0);
  hrzbpix0 = new TH2F("hrzbpix0","event display;r;z",240,-30.0,30.0,160,-20.0,20.0);
  hrzbpix1 = new TH2F("hrzbpix1","event display;r;z",240,-30.0,30.0,160,-20.0,20.0);
  hrzbpix2 = new TH2F("hrzbpix2","event display;r;z",240,-30.0,30.0,160,-20.0,20.0);
  hrzbpix3 = new TH2F("hrzbpix3","event display;r;z",240,-30.0,30.0,160,-20.0,20.0);
  hlastlayersame = new TH1F("hlastlayersame","last layer matched to gen particle;layer number",22,-1.5,20.5);
  hlastlayersame2 = new TH1F("hlastlayersame2","last layer matched to gen particle (w/ gap);layer number",22,-1.5,20.5);
  hlastlayercharged = new TH1F("hlastlayercharged","last layer matched to charged particle;layer number",22,-1.5,20.5);
  hlastlayercharged2 = new TH1F("hlastlayercharged2","last layer matched to charged particle (w/ gap);layer number",22,-1.5,20.5);
  hhitlayersame = new TH1F("hhitlayersame","continuous layers matched to gen particle;layer number",22,-1.5,20.5);
  hhitlayersame2 = new TH1F("hhitlayersame2","continuous layers matched to gen particle (w/ gap);layer number",22,-1.5,20.5);
  hhitlayercharged = new TH1F("hhitlayercharged","continuous layers matched to charged particle;layer number",22,-1.5,20.5);
  hhitlayercharged2 = new TH1F("hhilayercharged2","continuous layers matched to charged particle (w/ gap);layer number",22,-1.5,20.5);

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
GenSimEDMAnalyzer::~GenSimEDMAnalyzer() {

  float fr = (float)passevtcount0/(float)evtcount;
  cout << " passevtcount0 = " << passevtcount0 << " evtcount = " << evtcount << " fraction = " << fr << endl;
  fr = (float)passevtcount1/(float)evtcount;
  cout << " passevtcount1 = " << passevtcount1 << " evtcount = " << evtcount << " fraction = " << fr << endl;
  fr = (float)passevtcount00/(float)evtcount00;
  cout << " passevtcount00 = " << passevtcount00 << " evtcount00 = " << evtcount00 << " fraction = " << fr << endl;
  fr = (float)passevtcount01/(float)evtcount01;
  cout << " passevtcount01 = " << passevtcount01 << " evtcount01 = " << evtcount01 << " fraction = " << fr << endl;
  fr = (float)passevtcount11/(float)evtcount11;
  cout << " passevtcount11 = " << passevtcount11 << " evtcount11 = " << evtcount11 << " fraction = " << fr << endl;
  cout << " badiqstat count = " << badiqstat << endl;

  cout << endl;
  cout << " # of charged R-hadron |eta|<0.9 = " << ngencheta09 << endl;
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
    cout << "   # of matched hits in layer " << k << " = " << nsimmatch[k] << endl;
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
  cout << endl;

  TGraph *layermatch = new TGraph(NTRLAYERS, layerr, nsimmatch2);
  layermatch->SetName("layermatch");
  layermatch->SetTitle("# of Gen R-hadrons |eta|<0.9 matched to simhit");
  layermatch->GetXaxis()->SetTitle("radius (cm)");
  layermatch->GetYaxis()->SetTitle("# of particles with matched simhit");
  TGraph *layerfrac = new TGraph(NTRLAYERS, layerr, nsimmatchfrac);
  layerfrac->SetName("layerfrac");
  layerfrac->SetTitle("Fraction of Gen R-hadrons |eta|<0.9 matched to simhit");
  layerfrac->GetXaxis()->SetTitle("radius (cm)");
  layerfrac->GetYaxis()->SetTitle("fraction (matched/gen) particles");
  layermatch->Write();
  layerfrac->Write();
  TGraph *layernomatch = new TGraph(NTRLAYERS, layerr, nsimnomatch2);
  layernomatch->SetName("layernomatch");
  layernomatch->SetTitle("# of Gen R-hadrons |eta|<0.9 not matched to simhit");
  layernomatch->GetXaxis()->SetTitle("radius (cm)");
  layernomatch->GetYaxis()->SetTitle("# of particles with not matched simhit");
  TGraph *layernomatchqsame = new TGraph(NTRLAYERS, layerr, nsimnomatchqsame2);
  layernomatchqsame->SetName("layernomatchqsame");
  layernomatchqsame->SetTitle("# of Gen R-hadrons |eta|<0.9 not matched to simhit");
  layernomatchqsame->GetXaxis()->SetTitle("radius (cm)");
  layernomatchqsame->GetYaxis()->SetTitle("# of particles with not matched simhit");
  TGraph *layernomatchqdiff = new TGraph(NTRLAYERS, layerr, nsimnomatchqdiff2);
  layernomatchqdiff->SetName("layernomatchqdiff");
  layernomatchqdiff->SetTitle("# of Gen R-hadrons |eta|<0.9 not matched to simhit");
  layernomatchqdiff->GetXaxis()->SetTitle("radius (cm)");
  layernomatchqdiff->GetYaxis()->SetTitle("# of particles with not matched simhit");
  TGraph *layernofrac = new TGraph(NTRLAYERS, layerr, nsimnomatchfrac);
  layernofrac->SetName("layernofrac");
  layernofrac->SetTitle("Fraction of Gen R-hadrons |eta|<0.9 not matched to simhit");
  layernofrac->GetXaxis()->SetTitle("radius (cm)");
  layernofrac->GetYaxis()->SetTitle("fraction (not matched/gen) particles");
  TGraph *layernoqsamefrac = new TGraph(NTRLAYERS, layerr, nsimnomatchqsamefrac);
  layernoqsamefrac->SetName("layernoqsamefrac");
  layernoqsamefrac->SetTitle("Fraction of Gen R-hadrons |eta|<0.9 not matched to simhit");
  layernoqsamefrac->GetXaxis()->SetTitle("radius (cm)");
  layernoqsamefrac->GetYaxis()->SetTitle("fraction (not matched/gen) particles");
  TGraph *layernoqdifffrac = new TGraph(NTRLAYERS, layerr, nsimnomatchqdifffrac);
  layernoqdifffrac->SetName("layernoqdifffrac");
  layernoqdifffrac->SetTitle("Fraction of Gen R-hadrons |eta|<0.9 not matched to simhit");
  layernoqdifffrac->GetXaxis()->SetTitle("radius (cm)");
  layernoqdifffrac->GetYaxis()->SetTitle("fraction (not matched/gen) particles");
  layernomatch->Write();
  layernomatchqsame->Write();
  layernomatchqdiff->Write();
  layernofrac->Write();
  layernoqsamefrac->Write();
  layernoqdifffrac->Write();
  TGraph *layertotal = new TGraph(NTRLAYERS, layerr, nsimtotal2);
  layertotal->SetName("layertotal");
  layertotal->SetTitle("# of Gen R-hadrons |eta|<0.9 match or not matched to simhit");
  layertotal->GetXaxis()->SetTitle("radius (cm)");
  layertotal->GetYaxis()->SetTitle("# of particles with matched or not matched simhit");
  TGraph *layertotalfrac = new TGraph(NTRLAYERS, layerr, nsimtotalfrac);
  layertotalfrac->SetName("layertotalfrac");
  layertotalfrac->SetTitle("Fraction of Gen R-hadrons |eta|<0.9 matched or not matched to simhit");
  layertotalfrac->GetXaxis()->SetTitle("radius (cm)");
  layertotalfrac->GetYaxis()->SetTitle("fraction (match + not matched/gen) particles");
  layertotal->Write();
  layertotalfrac->Write();

//   sort(genrhadcounts.begin(),genrhadcounts.end(),GenSimEDMAnalyzer::comparePID);
   
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
   }

     // dump out particle info
   int np2=0,np1=0,n0=0,nn1=0,nn2=0;
   int nglue=0,nmeson=0,nbaryon=0,nlepton=0;
   cout << "Dump particle count table: " << genrhadcounts.size() << " " << sortgenpairs.size() << endl;
   for (std::vector<std::pair<int, int>>::iterator it = sortgenpairs.begin(); (it != sortgenpairs.end()); ++it) {
     cout << "  id = " << (*it).first << " charge = " << determineCharge((*it).first) << " count = " << (*it).second << endl;
     int cha = determineCharge((*it).first);
     if (cha==1) {
       np1 += (*it).second;
     } else if (cha==0) {
       n0 += (*it).second;
     } else if (cha==-1) {
       nn1 += (*it).second;
     } else if (cha==2) {
       np2 += (*it).second;
     } else if (cha==-2) {
       nn2 += (*it).second;
     }
     int pidabs = abs((*it).first);
     if ((pidabs>1000900)&&(pidabs<1000999)) {
       nglue += (*it).second;
     } else if ((pidabs>1000999)&&(pidabs<1009999)) {
       nmeson += (*it).second;
     } else if ((pidabs>1009999)&&(pidabs<1099999)) {
       nbaryon += (*it).second;
     } else if ((pidabs>1000009)&&(pidabs<1000020)) {
       nlepton += (*it).second;
     }

   }
   cout << " # of q=+2:  " << np2 << endl;
   cout << " # of q=+1:  " << np1 << endl;
   cout << " # of q=0:  " << n0 << endl;
   cout << " # of q=-1:  " << nn1 << endl;
   cout << " # of q=-2:  " << nn2 << endl;
   cout << " # of charged = " << (np2+np1+nn1+nn2) << " # of neutral = " << n0 << endl;
   int totpart = np2+np1+n0+nn1+nn2;
   cout << " # of gen particles = " << totpart << endl;

   cout << endl;
   cout << " # of gluinoballs = " << nglue << " fraction = " << ((float)nglue/(float)totpart) << endl;
   cout << " # of R-mesons = " << nmeson << " fraction = " << ((float)nmeson/(float)totpart) << endl;
   cout << " # of R-baryons = " << nbaryon << " fraction = " << ((float)nbaryon/(float)totpart) << endl;
   cout << " # of leptons = " << nlepton << " fraction = " << ((float)nlepton/(float)totpart) << endl;

/*TA
  std::cout << "\n\n\n\n\n" << std::endl;

  for(int i = 0; i < TLEN; ++i) {

    if(trig_pass[i]>100){
      std::cout << i+1 << ". " << trig_name[i] << " Passed: " << trig_pass[i] << "/" << trig_total[i] << std::endl;
    }

}

  std::cout << "\n\n\n" << std::endl;
*/


  //printing the list by ascending efficiencies
  /*sort(trig_pass);

  for(int i = 0; i < TLEN; ++i) {

    std::cout << " Passed: " << trig_pass[i] << "/" << trig_total[i] << std::endl;

    }*/


  std::cout << "\nMu50 OR PFMET120: " << trig_pass_or << "/" << trig_total[100] << std::endl;
  std::cout << "Mu50 AND PFMET120: " << trig_pass_and << "/" << trig_total[100] << std::endl;

  outputFile_->cd();

  htibrdist->Write();
  htobrdist->Write();
  hbpixrdist->Write();
  hdtrdist->Write();
  hcscrdist->Write();
  htibdeltaphi->Write();
  htibdeltaphi2->Write();
  htibdeltaeta->Write();
  htibdeltaR->Write();
  htibdeltaRclosest->Write();
  htobdeltaphi->Write();
  htobdeltaphi2->Write();
  htobdeltaeta->Write();
  htobdeltaR->Write();
  htobdeltaRclosest->Write();
  hbpixdeltaphi->Write();
  hbpixdeltaphi2->Write();
  hbpixdeltaeta->Write();
  hbpixdeltaR->Write();
  hbpixdeltaRclosest->Write();
  hdtdeltaphi->Write();
  hdtdeltaphi2->Write();
  hdtdeltaeta->Write();
  hdtdeltaR->Write();
  hdtdeltaRclosest->Write();
  htracker0->Write();
  htracker1->Write();
  htracker2->Write();
  htracker3->Write();
  htib0->Write();
  htib1->Write();
  htib2->Write();
  htib3->Write();
  hbpix0->Write();
  hbpix1->Write();
  hbpix2->Write();
  hbpix3->Write();
  hrzbpix0->Write();
  hrzbpix1->Write();
  hrzbpix2->Write();
  hrzbpix3->Write();
  hlastlayersame->Write();
  hlastlayersame2->Write();
  hlastlayercharged->Write();
  hlastlayercharged2->Write();
  hhitlayersame->Write();
  hhitlayersame2->Write();
  hhitlayercharged->Write();
  hhitlayercharged2->Write();

  hgenbeta->Write();
  h2genbeta->Write();
  hgenbetach->Write();
  hgenbetaneut->Write();
  hgendphi12->Write();
  hgendR12->Write();

  std::cout << "saving MET trigger histograms..." << std::endl;

/*  
  hltCalomet->Write();
  hltCalomet_muon->Write();
  hltCalomet_met120->Write();
  hltCalomet_met170->Write();
*/

/*
  hdeltaphicmet->Write();
  hdeltaphicmetgt100->Write();
  hdeltaphicmetlt100->Write();
  hdeltaphicmet2D->Write();
  hdeltaphicmet2Db->Write();
  hdeltaphimet->Write();
  hdeltaphimetgt100->Write();
  hdeltaphimetlt100->Write();
  hdeltaphimet2D->Write();
  hdeltaphimet2Db->Write();
  hdeltaphimetq0q0->Write();
  hdeltaphimetq0q1->Write();
  hdeltaphimetq1q1->Write();
  hdeltaphimetq0q2->Write();
  hdeltaphimetq1q2->Write();
  hdeltaphimetq2q2->Write();
  hdeltaphigen12->Write();
  
  hdRgenpf->Write();
  hdRgenpf00->Write();
  hdRgenpf01->Write(); hdRgenpf01n->Write(); hdRgenpf01c->Write();
  hdRgenpf11->Write();
*/
/*
  hdRgentrk->Write();
  hdRgentrk00->Write();
  hdRgentrk01->Write(); hdRgentrk01n->Write(); hdRgentrk01c->Write();
  hdRgentrk11->Write();
*/
/*
  hdphigenpf->Write();
  hdphigenpf00->Write();
  hdphigenpf01->Write(); hdphigenpf01n->Write(); hdphigenpf01c->Write();
  hdphigenpf11->Write();
  hdphigentrk->Write();
  hdphigentrk00->Write();
  hdphigentrk01->Write(); hdphigentrk01n->Write(); hdphigentrk01c->Write();
  hdphigentrk11->Write();
  hdptgenpf->Write();
  hdptgenpf00->Write();
  hdptgenpf01->Write(); hdptgenpf01n->Write(); hdptgenpf01c->Write();
  hdptgenpf11->Write();
  hdptgenpffrac->Write();
  hdptgenpffrac00->Write();
  hdptgenpffrac01->Write(); hdptgenpffrac01n->Write(); hdptgenpffrac01c->Write();
  hdptgenpffrac11->Write();
  hdptgentrk->Write();
  hdptgentrk00->Write();
  hdptgentrk01->Write(); hdptgentrk01n->Write(); hdptgentrk01c->Write();
  hdptgentrk11->Write();
  hdpttrkpf->Write();
  hdpttrkpf00->Write();
  hdpttrkpf01->Write(); hdpttrkpf01n->Write(); hdpttrkpf01c->Write();
  hdpttrkpf11->Write();
  h2pfptvsgenpt->Write();
  h2pfptvsgenpt00->Write();
  h2pfptvsgenpt01->Write();
  h2pfptvsgenpt01n->Write();
  h2pfptvsgenpt01c->Write();
  h2pfptvsgenpt11->Write();
  h2tkptvsgenpt->Write();
  h2tkptvspfpt->Write();
*/

  hgeneta->Write();
  hgenetaM->Write();

  outputFile_->Write();
  outputFile_->Close();

}

void GenSimEDMAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

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

   /////////////////////////////////
   // get muon information
   ////////////////////////////////
   // extract DT container
   edm::Handle<edm::PSimHitContainer> MuonDTContainer;
   iEvent.getByToken(edmPSimHitContainer_muonDT_Token_, MuonDTContainer);
   if (!MuonDTContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find MuonDT simhits in event!";
     return;
   }
   // extract CSC container
   edm::Handle<edm::PSimHitContainer> MuonCSCContainer;
   iEvent.getByToken(edmPSimHitContainer_muonCSC_Token_, MuonCSCContainer);
   if (!MuonCSCContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find MuonCSC simhits in event!";
     return;
   }

   // Get G4SimTracks
   edm::Handle<edm::SimTrackContainer> G4TrkContainer;
   iEvent.getByToken(edmSimTrackContainerToken_, G4TrkContainer);
   if (!G4TrkContainer.isValid()) {
     edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find SimTrack in event!";
     return;
   }

// Get geometry information
//   const auto &tracker = &iSetup.getData(tGeomEsToken_);
   edm::ESHandle<TrackerGeometry> tkGeometry;
   iSetup.get<TrackerDigiGeometryRecord>().get(tkGeometry);

   edm::ESHandle<DTGeometry> dtGeom;
   iSetup.get<MuonGeometryRecord>().get(dtGeom);
   edm::ESHandle<CSCGeometry> cscGeom;
   iSetup.get<MuonGeometryRecord>().get(cscGeom);

// Find R-hadrons
   const reco::GenParticle *genrhad1=0;
   const reco::GenParticle *genrhad2=0;

   for (unsigned int ig=0;ig<genColl->size();ig++) {
     const reco::GenParticle & part = (*genColl)[ig];
     if ((abs(part.pdgId())>999999)&&(part.status()==1)) {
   cout << "SUSY part = " << part.pdgId() << endl;
       if (genrhad1==0) {
         genrhad1 = &(*genColl)[ig];
       } else if (genrhad2==0) {
         genrhad2 = &(*genColl)[ig];
       } else {
         std::cout << "WARNING - found more than two R-hadrons" << std::endl;
       }
         // find and increment count
       bool found = false;

       for (std::vector<std::pair<int, int>>::iterator it = genrhadcounts.begin(); (it != genrhadcounts.end())&&(!found); ++it) {
         if (part.pdgId()==(*it).first) {
//           cout << " found it " << part.pdgId() << " " << (*it).first << endl;
           (*it).second++;
           found = true;
         }  
    /* std::cout << *it; ... */
       }
//       for (unsigned int i=0; ((i<genrhadcounts.size())&&(!found)); i++) {
//         if (part.pdgId()==genrhadcounts[i].first) {
//           cout << " found it " << part.pdgId() << " " << genrhadcounts[i].first << endl;
//           genrhadcounts[i].second++;
//         }  
//       }
       if (!found) {
         std::pair<int,int> pair1(part.pdgId(),1);
         genrhadcounts.push_back(pair1);
         cout << " NEW ID = " << part.pdgId() << " " << endl;
       }
     }
   }

   cout << "after loop SUSY part1,2 = " << genrhad1->pdgId() << " " << genrhad2->pdgId() << endl;
// determine combination of charges for generator particles
   int charge1 = determineCharge(genrhad1->pdgId());
   int charge2 = determineCharge(genrhad2->pdgId());
   int iqstat = -1;  // charge of gen particles: 0=q0q0 1=q0q1 2=q1q1
   if ((abs(charge1)==0)&&(abs(charge2)==0)) {
     iqstat = 0;
   } else if ((abs(charge1)==1)&&(abs(charge2)==1)) {
     iqstat = 2;
   } else if ((abs(charge1)+abs(charge2))==1) {
     cout << " iqstat==1 q1,2 = " << charge1 << " " << charge2 << endl;
     iqstat = 1;
   } else {
     std::cout << "WARNING - unknown charge state q=" << charge1 << " q=" << charge2 << endl;
   }
    cout << "iqstat = " << iqstat << endl;

    cout << " gen pids = " << genrhad1->pdgId() << " " << genrhad1->pdgId() << " " << (charge1+charge2) << endl;
     // plot beta
   float beta1 = genrhad1->p4().Beta();
   float beta2 = genrhad2->p4().Beta();
   hgenbeta->Fill(beta1);
   hgenbeta->Fill(beta2);
   if (beta1>beta2) {
     h2genbeta->Fill(beta1,beta2);
   } else {
     h2genbeta->Fill(beta2,beta1);
   }
   if (charge1==0) {
     hgenbetaneut->Fill(beta1);
   } else {
     hgenbetach->Fill(beta1);
   }
   if (charge2==0) {
     hgenbetaneut->Fill(beta2);
   } else {
     hgenbetach->Fill(beta2);
   }
   hgendphi12->Fill(fabs(genrhad1->phi()-genrhad2->phi()));
   hgendR12->Fill(deltaRta(genrhad1->phi(), genrhad1->eta(), genrhad2->phi(), genrhad2->eta()));

   hgeneta->Fill(genrhad1->eta());
   hgeneta->Fill(genrhad2->eta());

   if ((charge1!=0)&&(fabs(genrhad1->eta())<ETACUT)) ngencheta09++;
   if ((charge2!=0)&&(fabs(genrhad2->eta())<ETACUT)) ngencheta09++;

  // analyze SimHits
  // iterator to access containers
     edm::PSimHitContainer::const_iterator itHit;

  int icount = 0;
  int nhitstib[20] = {0,0,0,0,0,0,0,0,0,0};
  int nhitstib1[20] = {0,0,0,0,0,0,0,0,0,0};
  int nhitstib2[20] = {0,0,0,0,0,0,0,0,0,0};
  // start with TIB information
  for (itHit = SiTIBLowContainer->begin(); itHit != SiTIBLowContainer->end(); ++itHit) {
    DetId detid = DetId(itHit->detUnitId());
//    detid = 0;
//    const GeomDetUnit *det = (const GeomDetUnit *)tracker->idToDetUnit(detid);
//    cout << endl << icount << " rawId = " << detid.rawId();
/*    if (detid.null()) {
      cout << " warning detid is null " << endl;
    } else {
//      cout << " detid OK ";
//      const GeomDetUnit det = *tkGeometry->idToDetUnit(detid);
    } */
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
//const GeomDetUnit& geomDet = *tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(itHit->localPosition());
    float r = sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y());
/*    if (det.isValid()) {
    GlobalPoint gpos = det.toGlobal(itHit->localPosition());
    } else {
      cout << "det was invalid!" << endl;
    }
*/

//TA    if ((abs(itHit->particleType())>999999)||(abs(itHit->particleType())==13)) {
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<999999999)) {
//TA    if ((abs(itHit->particleType())>999999)) {
//      cout << "TIBlow simhit pid = " << itHit->particleType() << " " << itHit->phiAtEntry() << " " << itHit->energyLoss() << " eta = " << gpos.eta() << " phi = " << gpos.phi() << endl;
//      cout << " gpos x,y,z = " << gpos.x() << " " << gpos.y() << " " << gpos.z() << " r = " << sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y()) << endl;
      htibrdist->Fill(r);
//      if (evtcount<2) htracker1->Fill(gpos.x(),gpos.y());
      htracker0->Fill(gpos.x(),gpos.y());
      if (evtcount==TRKEVT1) htracker1->Fill(gpos.x(),gpos.y());
      if (evtcount==TRKEVT2) htracker2->Fill(gpos.x(),gpos.y());
      if (evtcount==TRKEVT3) htracker3->Fill(gpos.x(),gpos.y());
      htib0->Fill(gpos.x(),gpos.y());
      if (evtcount==TIBEVT1) htib1->Fill(gpos.x(),gpos.y());
      if (evtcount==TIBEVT2) htib2->Fill(gpos.x(),gpos.y());
      if (evtcount==TIBEVT3) htib3->Fill(gpos.x(),gpos.y());
      float dphi1 = genrhad1->phi() - gpos.phi();
      float dphi2 = genrhad2->phi() - gpos.phi();
      float deta1 = genrhad1->eta() - gpos.eta();
      float deta2 = genrhad2->eta() - gpos.eta();
      float dR1 = sqrt(dphi1*dphi1+deta1*deta1);
      float dR2 = sqrt(dphi2*dphi2+deta2*deta2);
      htibdeltaphi->Fill(dphi1);
      htibdeltaphi->Fill(dphi2);
      htibdeltaphi2->Fill(dphi1);
      htibdeltaphi2->Fill(dphi2);
      htibdeltaeta->Fill(deta1);
      htibdeltaeta->Fill(deta2);
      htibdeltaR->Fill(dR1);
      htibdeltaR->Fill(dR2);
      
      const reco::GenParticle *genrhadclosest=0;
      int iclosest = -1;
      float dRclosest = 999999.;
      if (dR1<dR2) {
        htibdeltaRclosest->Fill(dR1);
        dRclosest = dR1;
        iclosest = 1;
        genrhadclosest = genrhad1;
      } else {
        htibdeltaRclosest->Fill(dR2);
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
          cout << "didn't find layer " << r << " " << genrhadclosest->pdgId() << endl;
        }
/*        if ((r>23)&&(r<25)) { // first layer
          ilayer = TIB0 + 1;
        } else if ((r>26.2)&&(r<28.2)) { // second layer
          ilayer = TIB0 + 2;
        } else if ((r>31.5)&&(r<33.5)) { // third layer
          ilayer = TIB0 + 3;
        } else if ((r>34.5)&&(r<36.5)) { // fourth layer
          ilayer = TIB0 + 4;
        } else if ((r>39.5)&&(r<41.5)) { // fifth layer
          ilayer = TIB0 + 5;
        } else if ((r>42.5)&&(r<44.5)) { // sixth layer
          ilayer = TIB0 + 6;
        } else if ((r>47.5)&&(r<49.5)) { // seventh layer
          ilayer = TIB0 + 7;
        } else if ((r>50.5)&&(r<52.5)) { // eighth layer
          ilayer = TIB0 + 8;
        } else {
          cout << "didn't find layer " << r << " " << genrhadclosest->pdgId() << endl;
        }
*/
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
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<1000000000)&&(itHit->particleType()!=genrhad1->pdgId())&&(itHit->particleType()!=genrhad2->pdgId())) {
      cout << "simhit pid " << itHit->particleType() << "(" << determineCharge(itHit->particleType()) << ") does not match " << genrhad1->pdgId() << "(" << charge1 <<") or " << genrhad2->pdgId() << "(" << charge2 << ")" << endl;
    }
    icount++;
  }
//  cout << "nhitstib = " << nhitstib[0] << " " << nhitstib[1] << " " << nhitstib[2] << " " << nhitstib[3] << " " << endl;
//  cout << "nhitstib1 = " << nhitstib1[0] << " " << nhitstib1[1] << " " << nhitstib1[2] << " " << nhitstib1[3] << " " << endl;
//  cout << "nhitstib2 = " << nhitstib2[0] << " " << nhitstib2[1] << " " << nhitstib2[2] << " " << nhitstib2[3] << " " << endl;
  // go to TIB high information
  for (itHit = SiTIBHighContainer->begin(); itHit != SiTIBHighContainer->end(); ++itHit) {
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<999999999)) {
      cout << "TIBhigh simhit pid = " << itHit->particleType() << " " << itHit->phiAtEntry() << endl;
    }
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<1000000000)&&(itHit->particleType()!=genrhad1->pdgId())&&(itHit->particleType()!=genrhad2->pdgId())) {
      cout << "simhit pid " << itHit->particleType() << "(" << determineCharge(itHit->particleType()) << ") does not match " << genrhad1->pdgId() << "(" << charge1 <<") or " << genrhad2->pdgId() << "(" << charge2 << ")" << endl;
    }
  }

  // move on to T0B information
  for (itHit = SiTOBLowContainer->begin(); itHit != SiTOBLowContainer->end(); ++itHit) {
    DetId detid = DetId(itHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(itHit->localPosition());
    float r = sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y());
//    DetId detid = DetId(itHit->detUnitId());
//    const GeomDetUnit *det = (const GeomDetUnit *)tracker->idToDetUnit(detid);
//    GlobalPoint gpos = det->toGlobal(itHit->localPosition());
//    cout << "gpos = " << gpos.eta() << endl;
//TA    if ((abs(itHit->particleType())>999999)||(abs(itHit->particleType())==13)) {
//TA    if ((abs(itHit->particleType())>999999)) {
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<999999999)) {
//      cout << "TOBlow simhit pid = " << itHit->particleType() << " " << itHit->phiAtEntry() << endl;
      htobrdist->Fill(r);
      htracker0->Fill(gpos.x(),gpos.y());
      if (evtcount==TRKEVT1) htracker1->Fill(gpos.x(),gpos.y());
      if (evtcount==TRKEVT2) htracker2->Fill(gpos.x(),gpos.y());
      if (evtcount==TRKEVT3) htracker3->Fill(gpos.x(),gpos.y());

      float dphi1 = genrhad1->phi() - gpos.phi();
      float dphi2 = genrhad2->phi() - gpos.phi();
      float deta1 = genrhad1->eta() - gpos.eta();
      float deta2 = genrhad2->eta() - gpos.eta();
      float dR1 = sqrt(dphi1*dphi1+deta1*deta1);
      float dR2 = sqrt(dphi2*dphi2+deta2*deta2);
      htobdeltaphi->Fill(dphi1);
      htobdeltaphi->Fill(dphi2);
      htobdeltaphi2->Fill(dphi1);
      htobdeltaphi2->Fill(dphi2);
      htobdeltaeta->Fill(deta1);
      htobdeltaeta->Fill(deta2);
      htobdeltaR->Fill(dR1);
      htobdeltaR->Fill(dR2);
      
      const reco::GenParticle *genrhadclosest=0;
      float dRclosest = 999999.;
      int iclosest = -1;
      if (dR1<dR2) {
        htobdeltaRclosest->Fill(dR1);
        dRclosest = dR1;
        iclosest = 1;
        genrhadclosest = genrhad1;
      } else {
        htobdeltaRclosest->Fill(dR2);
        dRclosest = dR2;
        iclosest = 2;
        genrhadclosest = genrhad2;
      }
      int ilayer = -1;
      if (dRclosest<1.0) {
        if ((r>58)&&(r<64)) { // first layer
          ilayer = TOB0 + 1;
        } else if ((r>66)&&(r<72)) { // second layer
          ilayer = TOB0 + 2;
        } else if ((r>75)&&(r<81)) { // third layer
          ilayer = TOB0 + 3;
        } else if ((r>84)&&(r<90)) { // fourth layer
          ilayer = TOB0 + 4;
        } else if ((r>94)&&(r<100)) { // fifth layer
          ilayer = TOB0 + 5;
        } else if ((r>105)&&(r<111)) { // sixth layer
          ilayer = TOB0 + 6;
        } else {
          cout << "didn't find TOB layer " << r << " " << genrhadclosest->pdgId() << " " << ilayer << endl;
        }
/*        if ((r>58)&&(r<60.5)) { // first layer
          ilayer = TOB0 + 1;
        } else if ((r>61.2)&&(r<63.7)) { // second layer
          ilayer = TOB0 + 2;
        } else if ((r>66.5)&&(r<68.7)) { // third layer
          ilayer = TOB0 + 3;
        } else if ((r>69.8)&&(r<72.)) { // fourth layer
          ilayer = TOB0 + 4;
        } else if ((r>75.8)&&(r<77.2)) { // fifth layer
          ilayer = TOB0 + 5;
        } else if ((r>79.)&&(r<80.5)) { // sixth layer
          ilayer = TOB0 + 6;
        } else if ((r>84.5)&&(r<86.)) { // seventh layer
          ilayer = TOB0 + 7;
        } else if ((r>87.7)&&(r<89.2)) { // eighth layer
          ilayer = TOB0 + 8;
        } else if ((r>94.2)&&(r<95.8)) { // ninth layer
          ilayer = TOB0 + 9;
        } else if ((r>97.5)&&(r<99.)) { // tenth layer
          ilayer = TOB0 + 10;
        } else if ((r>105.5)&&(r<107.2)) { // eleventh layer
          ilayer = TOB0 + 11;
        } else if ((r>109.)&&(r<110.5)) { // twelvth layer
          ilayer = TOB0 + 12;
        } else {
          cout << "didn't find TOB layer " << r << " " << genrhadclosest->pdgId() << " " << ilayer << endl;
        }
*/
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
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<1000000000)&&(itHit->particleType()!=genrhad1->pdgId())&&(itHit->particleType()!=genrhad2->pdgId())) {
      cout << "simhit pid " << itHit->particleType() << "(" << determineCharge(itHit->particleType()) << ") does not match " << genrhad1->pdgId() << "(" << charge1 <<") or " << genrhad2->pdgId() << "(" << charge2 << ")" << endl;
    }
  }
  // go to TOB high information
  for (itHit = SiTOBHighContainer->begin(); itHit != SiTOBHighContainer->end(); ++itHit) {
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<999999999)) {
      cout << "TOBhigh simhit pid = " << itHit->particleType() << " " << itHit->phiAtEntry() << endl;
    }
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<1000000000)&&(itHit->particleType()!=genrhad1->pdgId())&&(itHit->particleType()!=genrhad2->pdgId())) {
      cout << "simhit pid " << itHit->particleType() << "(" << determineCharge(itHit->particleType()) << ") does not match " << genrhad1->pdgId() << "(" << charge1 <<") or " << genrhad2->pdgId() << "(" << charge2 << ")" << endl;
    }
  }

  // go back to pixel barrel information BPIX
  for (itHit = PxlBrlLowContainer->begin(); itHit != PxlBrlLowContainer->end(); ++itHit) {
    DetId detid = DetId(itHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(itHit->localPosition());
    float r = sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y());
    float rsign = r;
    if (gpos.y()<0) rsign = -1.0*r;
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<999999999)) {
//TA    if ((abs(itHit->particleType())>9)&&(abs(itHit->particleType())<19)) {  // lepton
      hbpixrdist->Fill(r);
      htracker0->Fill(gpos.x(),gpos.y());
      if (evtcount==TRKEVT1) htracker1->Fill(gpos.x(),gpos.y());
      if (evtcount==TRKEVT2) htracker2->Fill(gpos.x(),gpos.y());
      if (evtcount==TRKEVT3) htracker3->Fill(gpos.x(),gpos.y());
      hbpix0->Fill(gpos.x(),gpos.y());
      hrzbpix0->Fill(gpos.z(),rsign);
      if (evtcount==BPIXEVT1) hbpix1->Fill(gpos.x(),gpos.y());
      if (evtcount==BPIXEVT2) hbpix2->Fill(gpos.x(),gpos.y());
      if (evtcount==BPIXEVT3) hbpix3->Fill(gpos.x(),gpos.y());
      if (evtcount==BPIXEVT1) hrzbpix1->Fill(gpos.z(),rsign);
      if (evtcount==BPIXEVT2) hrzbpix2->Fill(gpos.z(),rsign);
      if (evtcount==BPIXEVT3) hrzbpix3->Fill(gpos.z(),rsign);

      float dphi1 = genrhad1->phi() - gpos.phi();
      float dphi2 = genrhad2->phi() - gpos.phi();
      float deta1 = genrhad1->eta() - gpos.eta();
      float deta2 = genrhad2->eta() - gpos.eta();
      float dR1 = sqrt(dphi1*dphi1+deta1*deta1);
      float dR2 = sqrt(dphi2*dphi2+deta2*deta2);
      hbpixdeltaphi->Fill(dphi1);
      hbpixdeltaphi->Fill(dphi2);
      hbpixdeltaphi2->Fill(dphi1);
      hbpixdeltaphi2->Fill(dphi2);
      hbpixdeltaeta->Fill(deta1);
      hbpixdeltaeta->Fill(deta2);
      hbpixdeltaR->Fill(dR1);
      hbpixdeltaR->Fill(dR2);
      
      const reco::GenParticle *genrhadclosest=0;
      float dRclosest = 999999.;
      float dphiclosest = 999999.;
      int iclosest = -1;
      if (dR1<dR2) {
        hbpixdeltaRclosest->Fill(dR1);
        dRclosest = dR1;
        dphiclosest = dphi1;
        iclosest = 1;
        genrhadclosest = genrhad1;
      } else {
        hbpixdeltaRclosest->Fill(dR2);
        dRclosest = dR2;
        dphiclosest = dphi2;
        iclosest = 2;
        genrhadclosest = genrhad2;
      }
    if (dRclosest>1.5) {
      cout << " BPIX dR = " << dR1 << " " << dR2 << endl;
      cout << " BPIX deta = " << deta1 << " " << deta2 << " " << gpos.eta() << " " << genrhad1->eta() << " " << genrhad2->eta() << endl;
      cout << " BPIX dphi = " << dphi1 << " " << dphi2 << endl;
    }
      int ilayer = -1;
//TA      if (dRclosest<1.0) {
//TA      if (dRclosest<2.0) {
//TA      if ((fabs(dphiclosest)<0.1)&&(dRclosest<3.0)) {
      if ((fabs(dphiclosest)<0.5)&&(dRclosest<2.5)) {
        if ((r>2.)&&(r<4.)) { // first layer
          ilayer = BPIX0 + 1;
        } else if ((r>6.)&&(r<8.)) { // second layer
          ilayer = BPIX0 + 2;
        } else if ((r>10.)&&(r<12.)) { // third layer
          ilayer = BPIX0 + 3;
        } else if ((r>15.)&&(r<17.)) { // fourth layer
          ilayer = BPIX0 + 4;
        } else {
          cout << "didn't find BPIX layer " << r << " " << genrhadclosest->pdgId() << " " << ilayer << endl;
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
  }
    // check out BPIX high
  for (itHit = PxlBrlHighContainer->begin(); itHit != PxlBrlHighContainer->end(); ++itHit) {
    DetId detid = DetId(itHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)tkGeometry->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(itHit->localPosition());
    float r = sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y());
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<999999999)) {
      cout << "PxlBrlhigh simhit pid = " << itHit->particleType() << " " << itHit->phiAtEntry() << endl;
    }
  }
 
  // get DT hits 
  for (itHit = MuonDTContainer->begin(); itHit != MuonDTContainer->end(); ++itHit) {
    DetId detid = DetId(itHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)dtGeom->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(itHit->localPosition());
    float r = sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y());
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<999999999)) {
      hdtrdist->Fill(r);
//      htracker0->Fill(gpos.x(),gpos.y());

      float dphi1 = genrhad1->phi() - gpos.phi();
      float dphi2 = genrhad2->phi() - gpos.phi();
      float deta1 = genrhad1->eta() - gpos.eta();
      float deta2 = genrhad2->eta() - gpos.eta();
      float dR1 = sqrt(dphi1*dphi1+deta1*deta1);
      float dR2 = sqrt(dphi2*dphi2+deta2*deta2);
      hdtdeltaphi->Fill(dphi1);
      hdtdeltaphi->Fill(dphi2);
      hdtdeltaphi2->Fill(dphi1);
      hdtdeltaphi2->Fill(dphi2);
      hdtdeltaeta->Fill(deta1);
      hdtdeltaeta->Fill(deta2);
      hdtdeltaR->Fill(dR1);
      hdtdeltaR->Fill(dR2);
      
      const reco::GenParticle *genrhadclosest=0;
      float dRclosest = 999999.;
      int iclosest = -1;
      if (dR1<dR2) {
        hdtdeltaRclosest->Fill(dR1);
        dRclosest = dR1;
        iclosest = 1;
        genrhadclosest = genrhad1;
      } else {
        hdtdeltaRclosest->Fill(dR2);
        dRclosest = dR2;
        iclosest = 2;
        genrhadclosest = genrhad2;
      }
      int ilayer = -1;
      if (dRclosest<1.0) {
        if ((r>410.)&&(r<475.)) { // first layer
          ilayer = DT0 + 1;
        } else if ((r>490.)&&(r<550.)) { // second layer
          ilayer = DT0 + 2;
        } else if ((r>600.)&&(r<660.)) { // third layer
          ilayer = DT0 + 3;
        } else if ((r>700.)&&(r<800.)) { // fourth layer
          ilayer = DT0 + 4;
        } else {
          cout << "didn't find DT layer " << r << " " << genrhadclosest->pdgId() << " " << ilayer << endl;
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
  }
  // get CSC hits
  for (itHit = MuonCSCContainer->begin(); itHit != MuonCSCContainer->end(); ++itHit) {
    DetId detid = DetId(itHit->detUnitId());
    const GeomDetUnit *det = (const GeomDetUnit *)cscGeom->idToDetUnit(detid);
    GlobalPoint gpos = det->toGlobal(itHit->localPosition());
    float r = sqrt(gpos.x()*gpos.x()+gpos.y()*gpos.y());
    if ((abs(itHit->particleType())>999999)&&(abs(itHit->particleType())<999999999)) {
      hcscrdist->Fill(r);
//      htracker0->Fill(gpos.x(),gpos.y());
    }
  }

    // sort linkedhits1 by layer number
  int count1 = 0;
  for (std::vector<simhitinfo>::iterator it1 = linkedhits1.begin(); it1 != linkedhits1.end(); ++it1) {
    int count2 = 0;
    for (std::vector<simhitinfo>::iterator it2 = linkedhits1.begin(); it2 != linkedhits1.end(); ++it2) {
      if ((count2>count1)&&(it2->ilayer < it1->ilayer)) {
        simhitinfo tmpsimhit = {it1->ilayer,it1->pdgID};
        it1->ilayer = it2->ilayer;
        it1->pdgID = it2->pdgID;
        it2->ilayer = tmpsimhit.ilayer;
        it2->pdgID = tmpsimhit.pdgID;
      }
      count2++;
    }
    count1++;
  }
    // sort linkedhits2 by layer number
  count1 = 0;
  for (std::vector<simhitinfo>::iterator it1 = linkedhits2.begin(); it1 != linkedhits2.end(); ++it1) {
    int count2 = 0;
    for (std::vector<simhitinfo>::iterator it2 = linkedhits2.begin(); it2 != linkedhits2.end(); ++it2) {
      if ((count2>count1)&&(it2->ilayer < it1->ilayer)) {
        simhitinfo tmpsimhit = {it1->ilayer,it1->pdgID};
        it1->ilayer = it2->ilayer;
        it1->pdgID = it2->pdgID;
        it2->ilayer = tmpsimhit.ilayer;
        it2->pdgID = tmpsimhit.pdgID;
      }
      count2++;
    }
    count1++;
  }

  if (linkedhits1.size()>0) {
    cout << "linkedhits1 list: " << endl;
    for (std::vector<simhitinfo>::iterator it = linkedhits1.begin(); it != linkedhits1.end(); ++it) {
      cout << " ilayer = " << it->ilayer << " " << it->pdgID << " " << genrhad1->pdgId() << endl;
    }
  }

  // analyze simhits by layer
  int lastlayerm = 0;
  int lastlayerc = 0;
  int lastlayerm2 = 0;
  int lastlayerc2 = 0;
  int layerhit[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int layernomatch[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int layernomatchqsame[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int layernomatchqdiff[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (std::vector<simhitinfo>::iterator it = linkedhits1.begin(); it != linkedhits1.end(); ++it) {
    if (fabs(genrhad1->eta())<ETACUT) {
      if (it->pdgID==genrhad1->pdgId()) {
        layerhit[it->ilayer] = 1;
        int difflayerm = it->ilayer - lastlayerm;
  cout << " last layer try = " << difflayerm << " " << it->ilayer << " " << lastlayerm << endl;
        if ((difflayerm==0)||(difflayerm==1)) lastlayerm = it->ilayer;
        if ((difflayerm==0)||(difflayerm==1)||(difflayerm==2)) lastlayerm2 = it->ilayer;
      } else {
        layernomatch[it->ilayer] = 1;
        if (determineCharge(it->pdgID)==charge1) {
          layernomatchqsame[it->ilayer] = 1;
        } else {
          layernomatchqdiff[it->ilayer] = 1;
        }
      }
      if (determineCharge(it->pdgID)!=0) {
        int difflayerc = it->ilayer - lastlayerc;
//  cout << " last layer tryc = " << difflayerc << " " << it->ilayer << " " << lastlayerc << endl;
        if ((difflayerc==0)||(difflayerc==1)) lastlayerc = it->ilayer;
        if ((difflayerc==0)||(difflayerc==1)||(difflayerc==2)) lastlayerc2 = it->ilayer;
      }
    }
  }
  for (int j=0;j<NTRLAYERS;j++) {
    nsimmatch[j] += layerhit[j];
    nsimnomatch[j] += layernomatch[j];
    nsimnomatchqsame[j] += layernomatchqsame[j];
    nsimnomatchqdiff[j] += layernomatchqdiff[j];
    if (j<=lastlayerm) hhitlayersame->Fill(j);
    if (j<=lastlayerm2) hhitlayersame2->Fill(j);
    if (j<=lastlayerc) hhitlayercharged->Fill(j);
    if (j<=lastlayerc2) hhitlayercharged2->Fill(j);
  }
  hlastlayersame->Fill(lastlayerm);
  hlastlayersame2->Fill(lastlayerm2);
  hlastlayercharged->Fill(lastlayerc);
  hlastlayercharged2->Fill(lastlayerc2);

  lastlayerm = 0;
  lastlayerc = 0;
  lastlayerm2 = 0;
  lastlayerc2 = 0;
  int layerhit2[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int layernomatch2[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int layernomatchqsame2[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int layernomatchqdiff2[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (std::vector<simhitinfo>::iterator it = linkedhits2.begin(); it != linkedhits2.end(); ++it) {
    if (fabs(genrhad2->eta())<ETACUT) {
      if (it->pdgID==genrhad2->pdgId()) {
        layerhit2[it->ilayer] = 1;
        int difflayerm = it->ilayer - lastlayerm;
  cout << " last layer try2 = " << difflayerm << " " << it->ilayer << " " << lastlayerm << endl;
        if ((difflayerm==0)||(difflayerm==1)) lastlayerm = it->ilayer;
        if ((difflayerm==0)||(difflayerm==1)||(difflayerm==2)) lastlayerm2 = it->ilayer;
      } else {
        layernomatch2[it->ilayer] = 1;
        if (determineCharge(it->pdgID)==charge2) {
          layernomatchqsame2[it->ilayer] = 1;
        } else {
          layernomatchqdiff2[it->ilayer] = 1;
        }
      }
      if (determineCharge(it->pdgID)!=0) {
        int difflayerc = it->ilayer - lastlayerc;
//  cout << " last layer tryc = " << difflayerc << " " << it->ilayer << " " << lastlayerc << endl;
        if ((difflayerc==0)||(difflayerc==1)) lastlayerc = it->ilayer;
        if ((difflayerc==0)||(difflayerc==1)||(difflayerc==2)) lastlayerc2 = it->ilayer;
      }
    }
  }
  for (int j=0;j<NTRLAYERS;j++) {
    nsimmatch[j] += layerhit2[j];
    nsimnomatch[j] += layernomatch2[j];
    nsimnomatchqsame[j] += layernomatchqsame2[j];
    nsimnomatchqdiff[j] += layernomatchqdiff2[j];
    if (j<=lastlayerm) hhitlayersame->Fill(j);
    if (j<=lastlayerm2) hhitlayersame2->Fill(j);
    if (j<=lastlayerc) hhitlayercharged->Fill(j);
    if (j<=lastlayerc2) hhitlayercharged2->Fill(j);
  }
  hlastlayersame->Fill(lastlayerm);
  hlastlayersame2->Fill(lastlayerm2);
  hlastlayercharged->Fill(lastlayerc);
  hlastlayercharged2->Fill(lastlayerc2);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenSimEDMAnalyzer);
