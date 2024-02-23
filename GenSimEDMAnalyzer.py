import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = "106X_mcRun3_2021_realistic_v3"

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/cmst3/user/gpetrucc/miniAOD/v1/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
        #'/store/user/lpchscp/noreplica/Run2016_274315.root'
        #'/store/user/lpchscp/noreplica/DY_13TeV_M800_Q3.root'
        #'/store/user/lpchscp/noreplica/GMStau_13TeV_M871.root'
        #'/store/user/lpchscp/noreplica/Gluino_13TeV_M1800.root'
#        '/store/user/lpchscp/noreplica/Gluino_13TeV_M1800N.root'
        #'/store/user/lpchscp/noreplica/PPStau_13TeV_M871.root'
        #'/store/user/lpchscp/noreplica/Stop_13TeV_M800.root'
        #'/store/user/lpchscp/noreplica/Stop_13TeV_M800N.root'
        # 'file:/uscms/home/rkim/gensimq36m300_tot.root'
        # 'file:/uscms/home/rkim/nobackup/CMSSW_8_0_30/src/recoq36m300_tot.root'
        # 'file:/uscms/home/rkim/nobackup/CMSSW_8_0_30/src/recoq3m300_tot.root'

#old        'file:/uscms_data/d2/tadams/hscp/fall22a/CMSSW_10_6_30/src/EXO-RunIISummer20UL18GENSIM-00010.root',
#        'file:/uscms_data/d2/tadams/hscp/fall22a/CMSSW_10_6_30/src/EXO-RunIISummer20UL18GENSIM-00010-v3.root',
#        'file:/uscms_data/d2/tadams/hscp/fall22a/CMSSW_10_6_30/src/EXO-RunIISummer20UL18GENSIM-00073-stau.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1800/NewHSCP_1.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1800/NewHSCP_2.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1800/NewHSCP_3.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1800/NewHSCP_4.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1800/NewHSCP_5.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1800/NewHSCP_6.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1800/NewHSCP_7.root',
         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1800/NewHSCP_8.root',

#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/HSCP_1.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/HSCP_2.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/HSCP_3.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/HSCP_4.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/HSCP_5.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_1.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_2.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_3.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_4.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_5.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_6.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_7.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_8.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_9.root',
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_10.root'
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_l1.root'
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_12.root'
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_13.root'
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_14.root'
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_15.root'
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_16.root'
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_17.root'
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_18.root'
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_19.root'
#         'file:/uscms/home/tadams/hscpdata/newsignal/gluinoM1400/RECO_20.root'
#        'root://cmsxrootd.fnal.gov///store/user/kazana/HSCP/MC2017/UL17_mc_noPU/HSCPgluino_M_1400_TuneCP5_13TeV_pythia8/UL17_mc_noPU/210324_153720/0000/RECO_1.root'
#        'root://cms-xrd-global.cern.ch///store/user/kazana/HSCP/MC2017/UL17_mc_noPU/HSCPgluino_M_1400_TuneCP5_13TeV_pythia8/UL17_mc_noPU/210324_153720/0000/RECO_1.root'
#        'root://cms-xrd-global.cern.ch///store/user/kazana/HSCP/MC2017/UL17_mc_noPU/HSCPgluino_M_1400_TuneCP5_13TeV_pythia8/UL17_mc_noPU/210324_153720/0000/RECO_2.root'

    )
)

process.demo = cms.EDAnalyzer("GenSimEDMAnalyzer",

#    gen_info = cms.InputTag("genParticles","","DIGI2RAW"),
    gen_info = cms.InputTag("genParticles","","SIM"),

    G4TrkSrc = cms.InputTag("g4SimHits"),

    PxlBrlLowSrc = cms.InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"),
    PxlBrlHighSrc = cms.InputTag("g4SimHits","TrackerHitsPixelBarrelHighTof"),
    PxlFwdLowSrc = cms.InputTag("g4SimHits","TrackerHitsPixelEndcapLowTof"),
    PxlFwdHighSrc = cms.InputTag("g4SimHits","TrackerHitsPixelEndcapHighTof"),
    SiTIBLowSrc = cms.InputTag("g4SimHits","TrackerHitsTIBLowTof"),
    SiTIBHighSrc = cms.InputTag("g4SimHits","TrackerHitsTIBHighTof"),
    SiTOBLowSrc = cms.InputTag("g4SimHits","TrackerHitsTOBLowTof"),
    SiTOBHighSrc = cms.InputTag("g4SimHits","TrackerHitsTOBHighTof"),
    SiTECLowSrc = cms.InputTag("g4SimHits","TrackerHitsTECLowTof"),
    SiTECHighSrc = cms.InputTag("g4SimHits","TrackerHitsTECHighTof"),
    SiTIDLowSrc = cms.InputTag("g4SimHits","TrackerHitsTIDLowTof"),
    SiTIDHighSrc = cms.InputTag("g4SimHits","TrackerHitsTIDHighTof"),

    MuonCSCSrc = cms.InputTag("g4SimHits","MuonCSCHits"),
    MuonDTSrc = cms.InputTag("g4SimHits","MuonDTHits"),

    bits = cms.InputTag("TriggerResults","","HLT"),
#    objects = cms.InputTag("selectedPatTrigger"),
#    prescales = cms.InputTag("patTrigger"),
    trig_sum = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    pfmet_reco = cms.InputTag("pfMet","","RECO"),
    calomet_reco = cms.InputTag("caloMet","","RECO"),
    dedx_info = cms.InputTag("dedxHitInfo","","RECO"),
    pf_reco = cms.InputTag("particleFlow"),
    vertex_reco = cms.InputTag("offlinePrimaryVertices"),
    secondaryvertex_reco = cms.InputTag("inclusiveSecondaryVertices"),
#did not work    vertex_reco = cms.InputTag("inclusiveSecondaryVertices"),
    gentrk_reco = cms.InputTag("generalTracks","","RECO"),
    hscp_cand = cms.InputTag("HSCParticleProducer","","HSCPAnalysis"),
    DedxCollection = cms.InputTag("dedxHitInfo","","RECO"),
#    TypeMode = cms.InputTag("dedxHitInfo","","RECO"),
    TypeMode = cms.untracked.int32(0),
    SampleType = cms.untracked.int32(0),
    Period = cms.untracked.string("2018"),
    UseTemplateLayer = cms.untracked.bool(False),
    SkipPixel = cms.untracked.bool(False),
    DeDxTemplate = cms.untracked.string("data/dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root"),
    EnableDeDxCalibration = cms.untracked.bool(False),
    DeDxCalibration = cms.untracked.string("Data13TeVGains_v2.root"),
    DeDxSF_0 = cms.untracked.double(1.0),
    DeDxSF_1 = cms.untracked.double(1.0325),


)

#process.demo.TypeMode=0
#process.analyzer.SampleType=0
#process.analyzer.saveTree=0 #all saved
#process.analyzer.saveGenTree=0
#process.analyzer.DeDxCalibration="Data13TeVGains_v2.root"
#process.analyzer.DeDxTemplate="dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root"

process.p = cms.Path(process.demo)
