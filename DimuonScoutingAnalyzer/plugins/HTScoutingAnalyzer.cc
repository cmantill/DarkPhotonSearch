// -*- C++ -*-
//
// Package:    Analyzer/HTScoutingAnalyzer
// Class:      HTScoutingAnalyzer
// 
/**\class HTScoutingAnalyzer HTScoutingAnalyzer.cc Analyzer/HTScoutingAnalyzer/plugins/HTScoutingAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Swagata Mukherjee
//         Created:  Mon, 26 Jun 2017 09:38:37 GMT
//
//


#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Scouting/interface/ScoutingMuon.h"
#include "DataFormats/Scouting/interface/ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/ScoutingCaloJet.h"
#include "DataFormats/Scouting/interface/ScoutingParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DarkPhotonSearch/DimuonScoutingAnalyzer/interface/EnergyCorrelations.h"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include <fastjet/JetDefinition.hh>
#include <TLorentzVector.h>
#include "TH1.h"
#include "TH2.h"
#include "fastjet/PseudoJet.hh"
#include <vector>
#include <map>
#include "TMath.h"
#include "TString.h"
#include "TTree.h"

using namespace std;
using namespace edm;

class HTScoutingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit HTScoutingAnalyzer(const edm::ParameterSet&);
      ~HTScoutingAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      static bool orderPseudoJet(fastjet::PseudoJet j1, fastjet::PseudoJet j2); 
      static std::vector<math::XYZTLorentzVector> makeP4s(const std::vector<fastjet::PseudoJet> &jets);
      const reco::BasicJet* match( const pat::Jet *iJet,const reco::BasicJetCollection *jets );

   private:
  
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::TriggerResults>            trgResultsLabel_;
  edm::EDGetTokenT<ScoutingCaloJetCollection>      caloJetLabel_;
  //edm::EDGetTokenT<double>                         caloRhoLabel_;
  edm::EDGetTokenT<ScoutingPFJetCollection>        pfJetLabel_;
  edm::EDGetTokenT<double>                         pfMetLabel_;
  edm::EDGetTokenT<pat::METCollection>             recoMetLabel_;
  edm::EDGetTokenT<ScoutingParticleCollection>     particleLabel_;
  edm::EDGetTokenT<pat::JetCollection>             recoJetLabel_;
  edm::EDGetTokenT<pat::JetCollection>             recoAK8JetLabel_;
  //edm::EDGetTokenT<reco::BasicJetCollection>       recoAK8SoftDropJetLabel_;
  //edm::EDGetTokenT<ScoutingVertexCollection>       pfVertexLabel_;
  edm::EDGetTokenT<ScoutingMuonCollection>         pfMuonLabel_;
  edm::EDGetTokenT<reco::VertexCollection>         recoVertexLabel_;
  edm::EDGetTokenT<pat::MuonCollection>            recoMuonLabel_;
  edm::Service<TFileService> fs;
  TTree *outTree_;

  int passrecoTT;
  int passrecoTTN2p;
  int passrecoTTN2f;
  int passrecoTTN2pl;
  int passrecoTTN2fl;
  int passpfTT;
  int passpfTTN2p;
  int passpfTTN2f;

  int passNominalHT250Trig;
  int passNominalHT410Trig;
  int passMonitoringTrig;

  // float pfAK8Msd_,pfAK8Pt_,recoAK8Msd_,recoAK8Pt_;

  double pfAK8HT;
  double pfAK8M;
  double pfAK8Msd;
  double pfAK8N2sdb1;
  double pfAK8Pt;
  double pfAK8dRMu;

  double caloRho;
  double caloHT;
  double caloMjj;
  double caloDeltaEtajj;
  double caloMjjWide; 
  double caloDeltaEtajjWide;
  double caloDeltaPhijjWide;

  double pfMet;
  double pfHT;
  double pfPt;
  double pfMjj;
  double pfDeltaEtajj;
  double pfMsdWide;
  double pfMjjWide; 
  double pfDeltaEtajjWide;
  double pfDeltaPhijjWide;
  double pfN2sdb1;
  double pfcsv;
  double pfdRAK8;
  double pfdRMu;

  double pfmuPt;
  double pfmuEta;
  double pfmuPtSel;
  double recomuPt;
  double recomuEta;
  double recomuPtSel;

  double recoMet;
  double recoHT;
  double recoPt;
  double recoMjj;
  double recoDeltaEtajj;
  double recoMjjWide; 
  double recoDeltaEtajjWide;
  double recoDeltaPhijjWide;
  double recocsv;
  double recodRAK8;
  double recodRMu;
  double recoMupt;
  double recoMueta;
  double recoAK8HT;
  double recoAK8M;
  double recoAK8Msd;
  double recoAK8N2sdb1;
  double recoAK8Pt;
  double recoAK8dRMu;
  
  TH1F *h1_pfAK8HT_nominalHT250_monitoring;
  TH1F *h1_pfAK8HT_nominalHT250;
  TH1F *h1_pfAK8HT_nominalHT410_monitoring;
  TH1F *h1_pfAK8HT_nominalHT410;
  TH1F *h1_pfAK8HT_monitoring;
  TH1F *h1_pfAK8M_nominalHT250_monitoring;
  TH1F *h1_pfAK8M_nominalHT250;
  TH1F *h1_pfAK8M_nominalHT410_monitoring;
  TH1F *h1_pfAK8M_nominalHT410;
  TH1F *h1_pfAK8M_monitoring;
  TH1F *h1_pfAK8Msd_nominalHT250_monitoring;
  TH1F *h1_pfAK8Msd_nominalHT250;
  TH1F *h1_pfAK8Msd_nominalHT410_monitoring;
  TH1F *h1_pfAK8Msd_nominalHT410;
  TH1F *h1_pfAK8Msd_monitoring;
  TH1F *h1_pfAK8N2sdb1_nominalHT250_monitoring;
  TH1F *h1_pfAK8N2sdb1_nominalHT250;
  TH1F *h1_pfAK8N2sdb1_nominalHT410_monitoring;
  TH1F *h1_pfAK8N2sdb1_nominalHT410;
  TH1F *h1_pfAK8N2sdb1_monitoring;
  TH1F *h1_pfAK8Pt_nominalHT250_monitoring;
  TH1F *h1_pfAK8Pt_nominalHT250;
  TH1F *h1_pfAK8Pt_nominalHT410_monitoring;
  TH1F *h1_pfAK8Pt_nominalHT410;
  TH1F *h1_pfAK8Pt_monitoring;
  TH2F *h2_pfAK8MsdPt_nominalHT250_monitoring;
  TH2F *h2_pfAK8MsdPt_nominalHT250;
  TH2F *h2_pfAK8MsdPt_nominalHT410_monitoring;
  TH2F *h2_pfAK8MsdPt_nominalHT410;
  TH2F *h2_pfAK8MsdPt_monitoring;

  TH1F *h1_pfAK8M_monitoring_passTT;
  TH1F *h1_pfAK8M_monitoring_passTTN2p;
  TH1F *h1_pfAK8M_monitoring_passTTN2f;
  TH1F *h1_pfAK8N2sdb1_monitoring_passTT;
  TH1F *h1_pfAK8N2sdb1_monitoring_passTTN2p;
  TH1F *h1_pfAK8N2sdb1_monitoring_passTTN2f;
  TH1F *h1_pfAK8Msd_monitoring_passTT;
  TH1F *h1_pfAK8Msd_monitoring_passTTN2p;
  TH1F *h1_pfAK8Msd_monitoring_passTTN2f;
  TH1F *h1_pfAK8Pt_monitoring_passTT;
  TH1F *h1_pfAK8Pt_monitoring_passTTN2p;
  TH1F *h1_pfAK8Pt_monitoring_passTTN2f;

  TH1F *h1_caloHT_nominalHT250_monitoring;
  TH1F *h1_caloHT_nominalHT410_monitoring;
  TH1F *h1_caloHT_monitoring;
  TH1F *h1_caloHT_nominalHT250;
  TH1F *h1_caloHT_nominalHT410;
  TH1F *h1_caloMjj_nominalHT250_monitoring;
  TH1F *h1_caloMjj_nominalHT410_monitoring;
  TH1F *h1_caloMjj_monitoring;
  TH1F *h1_caloMjj_nominalHT250;
  TH1F *h1_caloMjj_nominalHT410;
  TH1F *h1_caloDeltaEtajj_nominalHT250_monitoring;
  TH1F *h1_caloDeltaEtajj_nominalHT410_monitoring;
  TH1F *h1_caloDeltaEtajj_monitoring;
  TH1F *h1_caloDeltaEtajj_nominalHT250;
  TH1F *h1_caloDeltaEtajj_nominalHT410;
  TH1F *h1_caloMjjWide_nominalHT250_monitoring;
  TH1F *h1_caloMjjWide_nominalHT410_monitoring;
  TH1F *h1_caloMjjWide_monitoring;
  TH1F *h1_caloMjjWide_nominalHT250;
  TH1F *h1_caloMjjWide_nominalHT410;
  TH1F *h1_caloDeltaEtajjWide_nominalHT250_monitoring;
  TH1F *h1_caloDeltaEtajjWide_nominalHT410_monitoring;
  TH1F *h1_caloDeltaEtajjWide_monitoring;
  TH1F *h1_caloDeltaEtajjWide_nominalHT250;
  TH1F *h1_caloDeltaEtajjWide_nominalHT410;

  TH1F *h1_recoAK8HT_nominalHT250_monitoring;
  TH1F *h1_recoAK8HT_nominalHT250;
  TH1F *h1_recoAK8HT_nominalHT410_monitoring;
  TH1F *h1_recoAK8HT_nominalHT410;
  TH1F *h1_recoAK8HT_monitoring;
  TH1F *h1_recoAK8M_nominalHT250_monitoring;
  TH1F *h1_recoAK8M_nominalHT250;
  TH1F *h1_recoAK8M_nominalHT410_monitoring;
  TH1F *h1_recoAK8M_nominalHT410;
  TH1F *h1_recoAK8M_monitoring;
  TH1F *h1_recoAK8Msd_nominalHT250_monitoring;
  TH1F *h1_recoAK8Msd_nominalHT250;
  TH1F *h1_recoAK8Msd_nominalHT410_monitoring;
  TH1F *h1_recoAK8Msd_nominalHT410;
  TH1F *h1_recoAK8Msd_monitoring;
  TH1F *h1_recoAK8N2sdb1_nominalHT250_monitoring;
  TH1F *h1_recoAK8N2sdb1_nominalHT250;
  TH1F *h1_recoAK8N2sdb1_nominalHT410_monitoring;
  TH1F *h1_recoAK8N2sdb1_nominalHT410;
  TH1F *h1_recoAK8N2sdb1_monitoring;
  TH1F *h1_recoAK8Pt_nominalHT250_monitoring;
  TH1F *h1_recoAK8Pt_nominalHT250;
  TH1F *h1_recoAK8Pt_nominalHT410_monitoring;
  TH1F *h1_recoAK8Pt_nominalHT410;
  TH1F *h1_recoAK8Pt_monitoring;
  TH2F *h2_recoAK8MsdPt_nominalHT250_monitoring;
  TH2F *h2_recoAK8MsdPt_nominalHT250;
  TH2F *h2_recoAK8MsdPt_nominalHT410_monitoring;
  TH2F *h2_recoAK8MsdPt_nominalHT410;
  TH2F *h2_recoAK8MsdPt_monitoring;

  TH1F *h1_recoAK8M_monitoring_passTT;
  TH1F *h1_recoAK8M_monitoring_passTTN2p;
  TH1F *h1_recoAK8M_monitoring_passTTN2f;
  TH1F *h1_recoAK8M_monitoring_passTTN2pl;
  TH1F *h1_recoAK8M_monitoring_passTTN2fl;
  TH1F *h1_recoAK8N2sdb1_monitoring_passTT;
  TH1F *h1_recoAK8N2sdb1_monitoring_passTTN2p;
  TH1F *h1_recoAK8N2sdb1_monitoring_passTTN2f;
  TH1F *h1_recoAK8N2sdb1_monitoring_passTTN2pl;
  TH1F *h1_recoAK8N2sdb1_monitoring_passTTN2fl;
  TH1F *h1_recoAK8Msd_monitoring_passTT;
  TH1F *h1_recoAK8Msd_monitoring_passTTN2p;
  TH1F *h1_recoAK8Msd_monitoring_passTTN2f;
  TH1F *h1_recoAK8Msd_monitoring_passTTN2pl;
  TH1F *h1_recoAK8Msd_monitoring_passTTN2fl;
  TH1F *h1_recoAK8Msdmatched_monitoring_passTT;
  TH1F *h1_recoAK8Msdmatched_monitoring_passTTN2p;
  TH1F *h1_recoAK8Msdmatched_monitoring_passTTN2f;
  TH1F *h1_recoAK8Msdmatched_monitoring_passTTN2pl;
  TH1F *h1_recoAK8Msdmatched_monitoring_passTTN2fl;
  TH1F *h1_recoAK8Pt_monitoring_passTT;
  TH1F *h1_recoAK8Pt_monitoring_passTTN2p;
  TH1F *h1_recoAK8Pt_monitoring_passTTN2f;
  TH1F *h1_recoAK8Pt_monitoring_passTTN2pl;
  TH1F *h1_recoAK8Pt_monitoring_passTTN2fl;
  TH1F *h1_pfrecoAK8M_monitoring_passTT;
  TH1F *h1_pfrecoAK8M_monitoring_passTTN2p;
  TH1F *h1_pfrecoAK8M_monitoring_passTTN2f;
  TH1F *h1_pfrecoAK8M_monitoring_passTTN2pl;
  TH1F *h1_pfrecoAK8M_monitoring_passTTN2fl;
  TH1F *h1_pfrecoAK8N2sdb1_monitoring_passTT;
  TH1F *h1_pfrecoAK8N2sdb1_monitoring_passTTN2p;
  TH1F *h1_pfrecoAK8N2sdb1_monitoring_passTTN2f;
  TH1F *h1_pfrecoAK8N2sdb1_monitoring_passTTN2pl;
  TH1F *h1_pfrecoAK8N2sdb1_monitoring_passTTN2fl;
  TH1F *h1_pfrecoAK8Msd_monitoring_passTT;
  TH1F *h1_pfrecoAK8Msd_monitoring_passTTN2p;
  TH1F *h1_pfrecoAK8Msd_monitoring_passTTN2f;
  TH1F *h1_pfrecoAK8Msd_monitoring_passTTN2pl;
  TH1F *h1_pfrecoAK8Msd_monitoring_passTTN2fl;
  TH1F *h1_pfrecoAK8Pt_monitoring_passTT;
  TH1F *h1_pfrecoAK8Pt_monitoring_passTTN2p;
  TH1F *h1_pfrecoAK8Pt_monitoring_passTTN2f;


  TH1F *h1_pfHT_nominalHT250_monitoring;
  TH1F *h1_pfHT_nominalHT410_monitoring;
  TH1F *h1_pfHT_monitoring;
  TH1F *h1_pfHT_nominalHT250;
  TH1F *h1_pfHT_nominalHT410;
  TH1F *h1_pfMjj_nominalHT250_monitoring;
  TH1F *h1_pfMjj_nominalHT410_monitoring;
  TH1F *h1_pfMjj_monitoring;
  TH1F *h1_pfMjj_nominalHT250;
  TH1F *h1_pfMjj_nominalHT410;
  TH1F *h1_pfDeltaEtajj_nominalHT250_monitoring;
  TH1F *h1_pfDeltaEtajj_nominalHT410_monitoring;
  TH1F *h1_pfDeltaEtajj_monitoring;
  TH1F *h1_pfDeltaEtajj_nominalHT250;
  TH1F *h1_pfDeltaEtajj_nominalHT410;
  TH1F *h1_pfMjjWide_nominalHT250_monitoring;
  TH1F *h1_pfMjjWide_nominalHT410_monitoring;
  TH1F *h1_pfMjjWide_monitoring;
  TH1F *h1_pfMjjWide_nominalHT250;
  TH1F *h1_pfMjjWide_nominalHT410;
  TH1F *h1_pfDeltaEtajjWide_nominalHT250_monitoring;
  TH1F *h1_pfDeltaEtajjWide_nominalHT410_monitoring;
  TH1F *h1_pfDeltaEtajjWide_monitoring;
  TH1F *h1_pfDeltaEtajjWide_nominalHT250;
  TH1F *h1_pfDeltaEtajjWide_nominalHT410;

  TH1F *h1_recoHT_nominalHT250_monitoring;
  TH1F *h1_recoHT_nominalHT410_monitoring;
  TH1F *h1_recoHT_monitoring;
  TH1F *h1_recoHT_nominalHT250;
  TH1F *h1_recoHT_nominalHT410;
  TH1F *h1_recoMjj_nominalHT250_monitoring;
  TH1F *h1_recoMjj_nominalHT410_monitoring;
  TH1F *h1_recoMjj_monitoring;
  TH1F *h1_recoMjj_nominalHT250;
  TH1F *h1_recoMjj_nominalHT410;
  TH1F *h1_recoDeltaEtajj_nominalHT250_monitoring;
  TH1F *h1_recoDeltaEtajj_nominalHT410_monitoring;
  TH1F *h1_recoDeltaEtajj_monitoring;
  TH1F *h1_recoDeltaEtajj_nominalHT250;
  TH1F *h1_recoDeltaEtajj_nominalHT410;
  TH1F *h1_recoMjjWide_nominalHT250_monitoring;
  TH1F *h1_recoMjjWide_nominalHT410_monitoring;
  TH1F *h1_recoMjjWide_monitoring;
  TH1F *h1_recoMjjWide_nominalHT250;
  TH1F *h1_recoMjjWide_nominalHT410;
  TH1F *h1_recoN2sdb1_nominalHT250_monitoring;
  TH1F *h1_recoN2sdb1_nominalHT410_monitoring;
  TH1F *h1_recoN2sdb1_monitoring;
  TH1F *h1_recoN2sdb1_nominalHT250;
  TH1F *h1_recoN2sdb1_nominalHT410;
  TH1F *h1_recoDeltaEtajjWide_nominalHT250_monitoring;
  TH1F *h1_recoDeltaEtajjWide_nominalHT410_monitoring;
  TH1F *h1_recoDeltaEtajjWide_monitoring;
  TH1F *h1_recoDeltaEtajjWide_nominalHT250;
  TH1F *h1_recoDeltaEtajjWide_nominalHT410;

  // for JECs
  //bool doJECs_;
  //edm::FileInPath L1corrAK4_DATA_, L2corrAK4_DATA_, L3corrAK4_DATA_,ResCorrAK4_DATA_;
  // JetCorrectorParameters *L1ParAK4_DATA;
  // JetCorrectorParameters *L2ParAK4_DATA;
  // JetCorrectorParameters *L3ParAK4_DATA;
  // JetCorrectorParameters *L2L3ResAK4_DATA;
  // FactorizedJetCorrector *JetCorrectorAK4_DATA;
};

bool HTScoutingAnalyzer::orderPseudoJet(fastjet::PseudoJet j1, fastjet::PseudoJet j2) {
  // to be used to order pseudojets in decreasing pT order                                                                                                   
  return j1.perp2() > j2.perp2();
}

const reco::BasicJet* HTScoutingAnalyzer::match( const pat::Jet *iJet,const reco::BasicJetCollection *jets ) { 
  int lId = -1;
  double dRmin = 999.;
  for ( unsigned int i=0; i< jets->size(); ++i ) {
    const reco::BasicJet* jet = &(*jets)[i];
    float dR = deltaR( iJet->eta(), iJet->phi(), jet->eta(), jet->phi() );
    if(dR > 0.8) continue;
    if ( dR < dRmin ) {
      dRmin = dR;
      lId = i;
    }
  }
  const reco::BasicJet* lJet = 0; 
  if(lId != -1) lJet = &((*jets)[lId]);
  return lJet;
}

HTScoutingAnalyzer::HTScoutingAnalyzer(const edm::ParameterSet& iConfig)

{
  //doJECs_                  = iConfig.getParameter<bool>("doJECs");
  //L1corrAK4_DATA_          = iConfig.getParameter<edm::FileInPath>("L1corrAK4_DATA");
  //L2corrAK4_DATA_          = iConfig.getParameter<edm::FileInPath>("L2corrAK4_DATA");
  //L3corrAK4_DATA_          = iConfig.getParameter<edm::FileInPath>("L3corrAK4_DATA");
  //caloRhoLabel_            = consumes<double>(iConfig.getParameter<edm::InputTag>("caloRho")),
  trgResultsLabel_         = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"));
  caloJetLabel_            = consumes<ScoutingCaloJetCollection>(iConfig.getParameter<edm::InputTag>("caloJets"));
  pfJetLabel_              = consumes<ScoutingPFJetCollection>(iConfig.getParameter<edm::InputTag>("pfJets"));
  pfMetLabel_              = consumes<double>(iConfig.getParameter<edm::InputTag>("pfMetPt"));
  particleLabel_           = consumes<ScoutingParticleCollection>(iConfig.getParameter<edm::InputTag>("candidates"));
  recoJetLabel_            = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("recoJets"));
  recoMetLabel_            = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metreco"));
  recoAK8JetLabel_         = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("recoAK8Jets"));
  //recoAK8SoftDropJetLabel_ = consumes<reco::BasicJetCollection>(iConfig.getParameter<edm::InputTag>("softdropAK8Jets"));
  //pfVertexLabel_           = consumes<ScoutingVertexCollection>(iConfig.getParameter<InputTag>("vtx"));
  pfMuonLabel_             = consumes<ScoutingMuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  recoVertexLabel_         = consumes<reco::VertexCollection>(iConfig.getParameter<InputTag>("vtxreco"));
  recoMuonLabel_           = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("recoMuons"));

  usesResource("TFileService");

  outTree_ = fs->make<TTree>("events", "events");  
  outTree_->Branch("passrecoTT", &passrecoTT, "passrecoTT/I");
  outTree_->Branch("passrecoTTN2p", &passrecoTTN2p, "passrecoTTN2p/I");
  outTree_->Branch("passrecoTTN2f", &passrecoTTN2f, "passrecoTTN2f/I");
  outTree_->Branch("passrecoTTN2pl", &passrecoTTN2pl, "passrecoTTN2pl/I");
  outTree_->Branch("passrecoTTN2fl", &passrecoTTN2fl, "passrecoTTN2fl/I");
  outTree_->Branch("passpfTT", &passpfTT, "passpfTT/I");
  outTree_->Branch("passpfTTN2p", &passpfTTN2p, "passpfTTN2p/I");
  outTree_->Branch("passpfTTN2f", &passpfTTN2f, "passpfTTN2f/I");

  outTree_->Branch("passNominalHT250Trig", &passNominalHT250Trig, "passNominalHT250Trig/I");
  outTree_->Branch("passNominalHT410Trig", &passNominalHT410Trig, "passNominalHT410Trig/I");
  outTree_->Branch("passMonitoringTrig", &passMonitoringTrig, "passMonitoringTrig/I");

  outTree_->Branch("pfAK8Msd", &pfAK8Msd,  "pfAK8Msd/D");
  outTree_->Branch("pfAK8Pt" , &pfAK8Pt, "pfAK8Pt/D");
  outTree_->Branch("pfAK8N2sdb1", &pfAK8N2sdb1, "pfAK8N2sdb1/D");
  outTree_->Branch("pfAK8dRMu", &pfAK8dRMu, "pfAK8dRMu/D");
  outTree_->Branch("pfMet", &pfMet, "pfMet/D");
  outTree_->Branch("pfPt", &pfPt, "pfPt/D");
  outTree_->Branch("pfcsv", &pfcsv, "pfcsv/D");
  outTree_->Branch("pfdRAK8", &pfdRAK8, "pfdRAK8/D");
  outTree_->Branch("pfdRMu", &pfdRMu, "pfdRMu/D");
  
  outTree_->Branch("pfmuPt", &pfmuPt, "pfmuPt/D");
  outTree_->Branch("pfmuEta", &pfmuEta, "pfmuEta/D");
  outTree_->Branch("pfmuPtSel", &pfmuPtSel, "pfmuPtSel/D");
  outTree_->Branch("recomuPt", &recomuPt, "recomuPt/D");
  outTree_->Branch("recomuEta", &recomuEta, "recomuEta/D");
  outTree_->Branch("recomuPtSel", &recomuPtSel, "recomuPtSel/D");

  outTree_->Branch("recoAK8Msd", &recoAK8Msd,  "recoAK8Msd/D");
  outTree_->Branch("recoAK8Pt" , &recoAK8Pt, "recoAK8Pt/D");
  outTree_->Branch("recoAK8N2sdb1", &recoAK8N2sdb1, "recoAK8N2sdb1/D");
  outTree_->Branch("recoAK8dRMu", &recoAK8dRMu, "recoAK8dRMu/D");
  outTree_->Branch("recoMet", &recoMet, "recoMet/D");
  outTree_->Branch("recoPt", &recoPt, "recoPt/D");
  outTree_->Branch("recocsv", &recocsv, "recocsv/D");
  outTree_->Branch("recodRAK8", &recodRAK8, "recodRAK8/D");
  outTree_->Branch("recodRMu", &recodRMu, "recodRMu/D");

  TFileDirectory histoDir = fs->mkdir("histoDir");

  h1_pfAK8HT_nominalHT250_monitoring = histoDir.make<TH1F>("pfAK8HT_nominalHT250_monitoring", "pfAK8HT_nominalHT250_monitoring", 14000, 0, 14000);
  h1_pfAK8HT_nominalHT250            = histoDir.make<TH1F>("pfAK8HT_nominalHT250", "pfAK8HT_nominalHT250", 14000, 0, 14000);
  h1_pfAK8HT_nominalHT410_monitoring = histoDir.make<TH1F>("pfAK8HT_nominalHT410_monitoring", "pfAK8HT_nominalHT410_monitoring", 14000, 0, 14000);
  h1_pfAK8HT_nominalHT410            = histoDir.make<TH1F>("pfAK8HT_nominalHT410", "pfAK8HT_nominalHT410", 14000, 0, 14000);
  h1_pfAK8HT_monitoring              = histoDir.make<TH1F>("pfAK8HT_monitoring", "pfAK8HT_monitoring", 14000, 0, 14000);

  h1_pfAK8M_nominalHT250_monitoring = histoDir.make<TH1F>("pfAK8M_nominalHT250_monitoring", "pfAK8M_nominalHT250_monitoring", 40, 0, 250);
  h1_pfAK8M_nominalHT250            = histoDir.make<TH1F>("pfAK8M_nominalHT250", "pfAK8M_nominalHT250", 40, 0, 250);
  h1_pfAK8M_nominalHT410_monitoring = histoDir.make<TH1F>("pfAK8M_nominalHT410_monitoring", "pfAK8M_nominalHT410_monitoring", 40, 0, 250);
  h1_pfAK8M_nominalHT410            = histoDir.make<TH1F>("pfAK8M_nominalHT410", "pfAK8M_nominalHT410", 40, 0, 250);
  h1_pfAK8M_monitoring              = histoDir.make<TH1F>("pfAK8M_monitoring", "pfAK8M_monitoring", 40, 0, 250);

  h1_pfAK8Msd_nominalHT250_monitoring = histoDir.make<TH1F>("pfAK8Msd_nominalHT250_monitoring", "pfAK8Msd_nominalHT250_monitoring", 40, 0, 250);
  h1_pfAK8Msd_nominalHT250            = histoDir.make<TH1F>("pfAK8Msd_nominalHT250", "pfAK8Msd_nominalHT250", 40, 0, 250);
  h1_pfAK8Msd_nominalHT410_monitoring = histoDir.make<TH1F>("pfAK8Msd_nominalHT410_monitoring", "pfAK8Msd_nominalHT410_monitoring", 40, 0, 250);
  h1_pfAK8Msd_nominalHT410            = histoDir.make<TH1F>("pfAK8Msd_nominalHT410", "pfAK8Msd_nominalHT410", 40, 0, 250);
  h1_pfAK8Msd_monitoring              = histoDir.make<TH1F>("pfAK8Msd_monitoring", "pfAK8Msd_monitoring", 40, 0, 250);

  h1_pfAK8N2sdb1_nominalHT250_monitoring = histoDir.make<TH1F>("pfAK8N2sdb1_nominalHT250_monitoring", "pfAK8N2sdb1_nominalHT250_monitoring", 25,0,0.5);
  h1_pfAK8N2sdb1_nominalHT250            = histoDir.make<TH1F>("pfAK8N2sdb1_nominalHT250", "pfAK8N2sdb1_nominalHT250", 25,0,0.5);
  h1_pfAK8N2sdb1_nominalHT410_monitoring = histoDir.make<TH1F>("pfAK8N2sdb1_nominalHT410_monitoring", "pfAK8N2sdb1_nominalHT410_monitoring", 25,0,0.5);
  h1_pfAK8N2sdb1_nominalHT410            = histoDir.make<TH1F>("pfAK8N2sdb1_nominalHT410", "pfAK8N2sdb1_nominalHT410", 25,0,0.5);
  h1_pfAK8N2sdb1_monitoring              = histoDir.make<TH1F>("pfAK8N2sdb1_monitoring", "pfAK8N2sdb1_monitoring", 25,0,0.5);

  h1_pfAK8Pt_nominalHT250_monitoring = histoDir.make<TH1F>("pfAK8Pt_nominalHT250_monitoring", "pfAK8Pt_nominalHT250_monitoring", 20, 200, 700);
  h1_pfAK8Pt_nominalHT250            = histoDir.make<TH1F>("pfAK8Pt_nominalHT250", "pfAK8Pt_nominalHT250", 20, 200, 700);
  h1_pfAK8Pt_nominalHT410_monitoring = histoDir.make<TH1F>("pfAK8Pt_nominalHT410_monitoring", "pfAK8Pt_nominalHT410_monitoring", 20, 200, 700);
  h1_pfAK8Pt_nominalHT410            = histoDir.make<TH1F>("pfAK8Pt_nominalHT410", "pfAK8Pt_nominalHT410", 20, 200, 700);
  h1_pfAK8Pt_monitoring              = histoDir.make<TH1F>("pfAK8Pt_monitoring", "pfAK8Pt_monitoring", 20, 200, 700);

  h2_pfAK8MsdPt_nominalHT250_monitoring = histoDir.make<TH2F>("pfAK8MsdPt_nominalHT250_monitoring", "pfAK8MsdPt_nominalHT250_monitoring", 40, 0, 250, 20, 200, 700);
  h2_pfAK8MsdPt_nominalHT250            = histoDir.make<TH2F>("pfAK8MsdPt_nominalHT250", "pfAK8MsdPt_nominalHT250", 40, 0, 250, 20, 200, 700);
  h2_pfAK8MsdPt_nominalHT410_monitoring = histoDir.make<TH2F>("pfAK8MsdPt_nominalHT410_monitoring", "pfAK8MsdPt_nominalHT410_monitoring", 40, 0, 250, 20, 200, 700);
  h2_pfAK8MsdPt_nominalHT410            = histoDir.make<TH2F>("pfAK8MsdPt_nominalHT410", "pfAK8MsdPt_nominalHT410", 40, 0, 250, 20, 200, 700);
  h2_pfAK8MsdPt_monitoring              = histoDir.make<TH2F>("pfAK8MsdPt_monitoring", "pfAK8MsdPt_monitoring", 40, 0, 250, 20, 200, 700);

  h1_pfAK8M_monitoring_passTT      = histoDir.make<TH1F>("pfAK8M_monitoring_passTT", "pfAK8M_monitoring_passTT", 40, 0, 250);
  h1_pfAK8Msd_monitoring_passTT    = histoDir.make<TH1F>("pfAK8Msd_monitoring_passTT", "pfAK8Msd_monitoring_passTT", 40, 0, 250);
  h1_pfAK8N2sdb1_monitoring_passTT = histoDir.make<TH1F>("pfAK8N2sdb1_monitoring_passTT", "pfAK8N2sdb1_monitoring_passTT", 25, 0, 0.5);
  h1_pfAK8Pt_monitoring_passTT     = histoDir.make<TH1F>("pfAK8Pt_monitoring_passTT", "pfAK8Pt_monitoring_passTT", 20, 200, 700);
  h1_pfAK8M_monitoring_passTTN2p      = histoDir.make<TH1F>("pfAK8M_monitoring_passTTN2p", "pfAK8M_monitoring_passTTN2p", 40, 0, 250);
  h1_pfAK8Msd_monitoring_passTTN2p    = histoDir.make<TH1F>("pfAK8Msd_monitoring_passTTN2p", "pfAK8Msd_monitoring_passTTN2p", 40, 0, 250);
  h1_pfAK8N2sdb1_monitoring_passTTN2p = histoDir.make<TH1F>("pfAK8N2sdb1_monitoring_passTTN2p", "pfAK8N2sdb1_monitoring_passTTN2p", 25, 0, 0.5);
  h1_pfAK8Pt_monitoring_passTTN2p     = histoDir.make<TH1F>("pfAK8Pt_monitoring_passTTN2p", "pfAK8Pt_monitoring_passTTN2p", 20, 200, 700);
  h1_pfAK8M_monitoring_passTTN2f      = histoDir.make<TH1F>("pfAK8M_monitoring_passTTN2f", "pfAK8M_monitoring_passTTN2f", 40, 0, 250);
  h1_pfAK8Msd_monitoring_passTTN2f    = histoDir.make<TH1F>("pfAK8Msd_monitoring_passTTN2f", "pfAK8Msd_monitoring_passTTN2f", 40, 0, 250);
  h1_pfAK8N2sdb1_monitoring_passTTN2f = histoDir.make<TH1F>("pfAK8N2sdb1_monitoring_passTTN2f", "pfAK8N2sdb1_monitoring_passTTN2f", 25, 0, 0.5);
  h1_pfAK8Pt_monitoring_passTTN2f     = histoDir.make<TH1F>("pfAK8Pt_monitoring_passTTN2f", "pfAK8Pt_monitoring_passTTN2f", 20, 200, 700);

  h1_caloHT_nominalHT250_monitoring = histoDir.make<TH1F>("caloHT_nominalHT250_monitoring", "caloHT_nominalHT250_monitoring", 14000, 0, 14000);
  h1_caloHT_nominalHT250            = histoDir.make<TH1F>("caloHT_nominalHT250", "caloHT_nominalHT250", 14000, 0, 14000);
  h1_caloHT_nominalHT410_monitoring = histoDir.make<TH1F>("caloHT_nominalHT410_monitoring", "caloHT_nominalHT410_monitoring", 14000, 0, 14000);
  h1_caloHT_nominalHT410            = histoDir.make<TH1F>("caloHT_nominalHT410", "caloHT_nominalHT410", 14000, 0, 14000);
  h1_caloHT_monitoring              = histoDir.make<TH1F>("caloHT_monitoring", "caloHT_monitoring", 14000, 0, 14000);

  h1_caloMjj_nominalHT250_monitoring = histoDir.make<TH1F>("caloMjj_nominalHT250_monitoring", "caloMjj_nominalHT250_monitoring", 14000, 0, 14000);
  h1_caloMjj_nominalHT250            = histoDir.make<TH1F>("caloMjj_nominalHT250", "caloMjj_nominalHT250", 14000, 0, 14000);
  h1_caloMjj_nominalHT410_monitoring = histoDir.make<TH1F>("caloMjj_nominalHT410_monitoring", "caloMjj_nominalHT410_monitoring", 14000, 0, 14000);
  h1_caloMjj_nominalHT410            = histoDir.make<TH1F>("caloMjj_nominalHT410", "caloMjj_nominalHT410", 14000, 0, 14000);
  h1_caloMjj_monitoring              = histoDir.make<TH1F>("caloMjj_monitoring", "caloMjj_monitoring", 14000, 0, 14000);
  
  h1_caloDeltaEtajj_nominalHT250_monitoring = histoDir.make<TH1F>("caloDeltaEtajj_nominalHT250_monitoring", "caloDeltaEtajj_nominalHT250_monitoring", 1000, 0, 5);
  h1_caloDeltaEtajj_nominalHT250            = histoDir.make<TH1F>("caloDeltaEtajj_nominalHT250", "caloDeltaEtajj_nominalHT250", 1000, 0, 5);
  h1_caloDeltaEtajj_nominalHT410_monitoring = histoDir.make<TH1F>("caloDeltaEtajj_nominalHT410_monitoring", "caloDeltaEtajj_nominalHT410_monitoring", 1000, 0, 5);
  h1_caloDeltaEtajj_nominalHT410            = histoDir.make<TH1F>("caloDeltaEtajj_nominalHT410", "caloDeltaEtajj_nominalHT410", 1000, 0, 5);
  h1_caloDeltaEtajj_monitoring              = histoDir.make<TH1F>("caloDeltaEtajj_monitoring", "caloDeltaEtajj_monitoring", 1000, 0, 5);
  
  h1_caloMjjWide_nominalHT250_monitoring = histoDir.make<TH1F>("caloMjjWide_nominalHT250_monitoring", "caloMjjWide_nominalHT250_monitoring", 14000, 0, 14000);
  h1_caloMjjWide_nominalHT250            = histoDir.make<TH1F>("caloMjjWide_nominalHT250", "caloMjjWide_nominalHT250", 14000, 0, 14000);
  h1_caloMjjWide_nominalHT410_monitoring = histoDir.make<TH1F>("caloMjjWide_nominalHT410_monitoring", "caloMjjWide_nominalHT410_monitoring", 14000, 0, 14000);
  h1_caloMjjWide_nominalHT410            = histoDir.make<TH1F>("caloMjjWide_nominalHT410", "caloMjjWide_nominalHT410", 14000, 0, 14000);
  h1_caloMjjWide_monitoring              = histoDir.make<TH1F>("caloMjjWide_monitoring", "caloMjjWide_monitoring", 14000, 0, 14000);
  
  h1_caloDeltaEtajjWide_nominalHT250_monitoring = histoDir.make<TH1F>("caloDeltaEtajjWide_nominalHT250_monitoring", "caloDeltaEtajjWide_nominalHT250_monitoring", 1000, 0, 5);
  h1_caloDeltaEtajjWide_nominalHT250            = histoDir.make<TH1F>("caloDeltaEtajjWide_nominalHT250", "caloDeltaEtajjWide_nominalHT250", 1000, 0, 5);
  h1_caloDeltaEtajjWide_nominalHT410_monitoring = histoDir.make<TH1F>("caloDeltaEtajjWide_nominalHT410_monitoring", "caloDeltaEtajjWide_nominalHT410_monitoring", 1000, 0, 5);
  h1_caloDeltaEtajjWide_nominalHT410            = histoDir.make<TH1F>("caloDeltaEtajjWide_nominalHT410", "caloDeltaEtajjWide_nominalHT410", 1000, 0, 5);
  h1_caloDeltaEtajjWide_monitoring              = histoDir.make<TH1F>("caloDeltaEtajjWide_monitoring", "caloDeltaEtajjWide_monitoring", 1000, 0, 5);

  h1_pfHT_nominalHT250_monitoring = histoDir.make<TH1F>("pfHT_nominalHT250_monitoring", "pfHT_nominalHT250_monitoring", 14000, 0, 14000);
  h1_pfHT_nominalHT250            = histoDir.make<TH1F>("pfHT_nominalHT250", "pfHT_nominalHT250", 14000, 0, 14000);
  h1_pfHT_nominalHT410_monitoring = histoDir.make<TH1F>("pfHT_nominalHT410_monitoring", "pfHT_nominalHT410_monitoring", 14000, 0, 14000);
  h1_pfHT_nominalHT410            = histoDir.make<TH1F>("pfHT_nominalHT410", "pfHT_nominalHT410", 14000, 0, 14000);
  h1_pfHT_monitoring              = histoDir.make<TH1F>("pfHT_monitoring", "pfHT_monitoring", 14000, 0, 14000);

  h1_pfMjj_nominalHT250_monitoring = histoDir.make<TH1F>("pfMjj_nominalHT250_monitoring", "pfMjj_nominalHT250_monitoring", 14000, 0, 14000);
  h1_pfMjj_nominalHT250            = histoDir.make<TH1F>("pfMjj_nominalHT250", "pfMjj_nominalHT250", 14000, 0, 14000);
  h1_pfMjj_nominalHT410_monitoring = histoDir.make<TH1F>("pfMjj_nominalHT410_monitoring", "pfMjj_nominalHT410_monitoring", 14000, 0, 14000);
  h1_pfMjj_nominalHT410            = histoDir.make<TH1F>("pfMjj_nominalHT410", "pfMjj_nominalHT410", 14000, 0, 14000);
  h1_pfMjj_monitoring              = histoDir.make<TH1F>("pfMjj_monitoring", "pfMjj_monitoring", 14000, 0, 14000);

  h1_pfDeltaEtajj_nominalHT250_monitoring = histoDir.make<TH1F>("pfDeltaEtajj_nominalHT250_monitoring", "pfDeltaEtajj_nominalHT250_monitoring", 1000, 0, 5);
  h1_pfDeltaEtajj_nominalHT250            = histoDir.make<TH1F>("pfDeltaEtajj_nominalHT250", "pfDeltaEtajj_nominalHT250", 1000, 0, 5);
  h1_pfDeltaEtajj_nominalHT410_monitoring = histoDir.make<TH1F>("pfDeltaEtajj_nominalHT410_monitoring", "pfDeltaEtajj_nominalHT410_monitoring", 1000, 0, 5);
  h1_pfDeltaEtajj_nominalHT410            = histoDir.make<TH1F>("pfDeltaEtajj_nominalHT410", "pfDeltaEtajj_nominalHT410", 1000, 0, 5);
  h1_pfDeltaEtajj_monitoring              = histoDir.make<TH1F>("pfDeltaEtajj_monitoring", "pfDeltaEtajj_monitoring", 1000, 0, 5);
  
  h1_pfMjjWide_nominalHT250_monitoring = histoDir.make<TH1F>("pfMjjWide_nominalHT250_monitoring", "pfMjjWide_nominalHT250_monitoring", 14000, 0, 14000);
  h1_pfMjjWide_nominalHT250            = histoDir.make<TH1F>("pfMjjWide_nominalHT250", "pfMjjWide_nominalHT250", 14000, 0, 14000);
  h1_pfMjjWide_nominalHT410_monitoring = histoDir.make<TH1F>("pfMjjWide_nominalHT410_monitoring", "pfMjjWide_nominalHT410_monitoring", 14000, 0, 14000);
  h1_pfMjjWide_nominalHT410            = histoDir.make<TH1F>("pfMjjWide_nominalHT410", "pfMjjWide_nominalHT410", 14000, 0, 14000);
  h1_pfMjjWide_monitoring              = histoDir.make<TH1F>("pfMjjWide_monitoring", "pfMjjWide_monitoring", 14000, 0, 14000);
  
  h1_pfDeltaEtajjWide_nominalHT250_monitoring = histoDir.make<TH1F>("pfDeltaEtajjWide_nominalHT250_monitoring", "pfDeltaEtajjWide_nominalHT250_monitoring", 1000, 0, 5);
  h1_pfDeltaEtajjWide_nominalHT250            = histoDir.make<TH1F>("pfDeltaEtajjWide_nominalHT250", "pfDeltaEtajjWide_nominalHT250", 1000, 0, 5);
  h1_pfDeltaEtajjWide_nominalHT410_monitoring = histoDir.make<TH1F>("pfDeltaEtajjWide_nominalHT410_monitoring", "pfDeltaEtajjWide_nominalHT410_monitoring", 1000, 0, 5);
  h1_pfDeltaEtajjWide_nominalHT410            = histoDir.make<TH1F>("pfDeltaEtajjWide_nominalHT410", "pfDeltaEtajjWide_nominalHT410", 1000, 0, 5);
  h1_pfDeltaEtajjWide_monitoring              = histoDir.make<TH1F>("pfDeltaEtajjWide_monitoring", "pfDeltaEtajjWide_monitoring", 1000, 0, 5);

  h1_recoHT_nominalHT250_monitoring = histoDir.make<TH1F>("recoHT_nominalHT250_monitoring", "recoHT_nominalHT250_monitoring", 14000, 0, 14000);
  h1_recoHT_nominalHT250            = histoDir.make<TH1F>("recoHT_nominalHT250", "recoHT_nominalHT250", 14000, 0, 14000);
  h1_recoHT_nominalHT410_monitoring = histoDir.make<TH1F>("recoHT_nominalHT410_monitoring", "recoHT_nominalHT410_monitoring", 14000, 0, 14000);
  h1_recoHT_nominalHT410            = histoDir.make<TH1F>("recoHT_nominalHT410", "recoHT_nominalHT410", 14000, 0, 14000);
  h1_recoHT_monitoring              = histoDir.make<TH1F>("recoHT_monitoring", "recoHT_monitoring", 14000, 0, 14000);

  h1_recoMjj_nominalHT250_monitoring = histoDir.make<TH1F>("recoMjj_nominalHT250_monitoring", "recoMjj_nominalHT250_monitoring", 14000, 0, 14000);
  h1_recoMjj_nominalHT250            = histoDir.make<TH1F>("recoMjj_nominalHT250", "recoMjj_nominalHT250", 14000, 0, 14000);
  h1_recoMjj_nominalHT410_monitoring = histoDir.make<TH1F>("recoMjj_nominalHT410_monitoring", "recoMjj_nominalHT410_monitoring", 14000, 0, 14000);
  h1_recoMjj_nominalHT410            = histoDir.make<TH1F>("recoMjj_nominalHT410", "recoMjj_nominalHT410", 14000, 0, 14000);
  h1_recoMjj_monitoring              = histoDir.make<TH1F>("recoMjj_monitoring", "recoMjj_monitoring", 14000, 0, 14000);
  
  h1_recoN2sdb1_nominalHT250_monitoring = histoDir.make<TH1F>("recoN2sdb1_nominalHT250_monitoring", "recoN2sdb1_nominalHT250_monitoring", 25, 0, 0.5);
  h1_recoN2sdb1_nominalHT250            = histoDir.make<TH1F>("recoN2sdb1_nominalHT250", "recoN2sdb1_nominalHT250", 25, 0, 0.5);
  h1_recoN2sdb1_nominalHT410_monitoring = histoDir.make<TH1F>("recoN2sdb1_nominalHT410_monitoring", "recoN2sdb1_nominalHT410_monitoring", 25, 0, 0.5);
  h1_recoN2sdb1_nominalHT410            = histoDir.make<TH1F>("recoN2sdb1_nominalHT410", "recoN2sdb1_nominalHT410", 25, 0, 0.5);
  h1_recoN2sdb1_monitoring              = histoDir.make<TH1F>("recoN2sdb1_monitoring", "recoN2sdb1_monitoring", 25, 0, 0.5);

  h1_recoDeltaEtajj_nominalHT250_monitoring = histoDir.make<TH1F>("recoDeltaEtajj_nominalHT250_monitoring", "recoDeltaEtajj_nominalHT250_monitoring", 1000, 0, 5);
  h1_recoDeltaEtajj_nominalHT250            = histoDir.make<TH1F>("recoDeltaEtajj_nominalHT250", "recoDeltaEtajj_nominalHT250", 1000, 0, 5);
  h1_recoDeltaEtajj_nominalHT410_monitoring = histoDir.make<TH1F>("recoDeltaEtajj_nominalHT410_monitoring", "recoDeltaEtajj_nominalHT410_monitoring", 1000, 0, 5);
  h1_recoDeltaEtajj_nominalHT410            = histoDir.make<TH1F>("recoDeltaEtajj_nominalHT410", "recoDeltaEtajj_nominalHT410", 1000, 0, 5);
  h1_recoDeltaEtajj_monitoring              = histoDir.make<TH1F>("recoDeltaEtajj_monitoring", "recoDeltaEtajj_monitoring", 1000, 0, 5);
  
  h1_recoMjjWide_nominalHT250_monitoring = histoDir.make<TH1F>("recoMjjWide_nominalHT250_monitoring", "recoMjjWide_nominalHT250_monitoring", 14000, 0, 14000);
  h1_recoMjjWide_nominalHT250            = histoDir.make<TH1F>("recoMjjWide_nominalHT250", "recoMjjWide_nominalHT250", 14000, 0, 14000);
  h1_recoMjjWide_nominalHT410_monitoring = histoDir.make<TH1F>("recoMjjWide_nominalHT410_monitoring", "recoMjjWide_nominalHT410_monitoring", 14000, 0, 14000);
  h1_recoMjjWide_nominalHT410            = histoDir.make<TH1F>("recoMjjWide_nominalHT410", "recoMjjWide_nominalHT410", 14000, 0, 14000);
  h1_recoMjjWide_monitoring              = histoDir.make<TH1F>("recoMjjWide_monitoring", "recoMjjWide_monitoring", 14000, 0, 14000);
  
  h1_recoDeltaEtajjWide_nominalHT250_monitoring = histoDir.make<TH1F>("recoDeltaEtajjWide_nominalHT250_monitoring", "recoDeltaEtajjWide_nominalHT250_monitoring", 1000, 0, 5);
  h1_recoDeltaEtajjWide_nominalHT250            = histoDir.make<TH1F>("recoDeltaEtajjWide_nominalHT250", "recoDeltaEtajjWide_nominalHT250", 1000, 0, 5);
  h1_recoDeltaEtajjWide_nominalHT410_monitoring = histoDir.make<TH1F>("recoDeltaEtajjWide_nominalHT410_monitoring", "recoDeltaEtajjWide_nominalHT410_monitoring", 1000, 0, 5);
  h1_recoDeltaEtajjWide_nominalHT410            = histoDir.make<TH1F>("recoDeltaEtajjWide_nominalHT410", "recoDeltaEtajjWide_nominalHT410", 1000, 0, 5);
  h1_recoDeltaEtajjWide_monitoring              = histoDir.make<TH1F>("recoDeltaEtajjWide_monitoring", "recoDeltaEtajjWide_monitoring", 1000, 0, 5);

  h1_recoAK8HT_nominalHT250_monitoring = histoDir.make<TH1F>("recoAK8HT_nominalHT250_monitoring", "recoAK8HT_nominalHT250_monitoring", 14000, 0, 14000);
  h1_recoAK8HT_nominalHT250            = histoDir.make<TH1F>("recoAK8HT_nominalHT250", "recoAK8HT_nominalHT250", 14000, 0, 14000);
  h1_recoAK8HT_nominalHT410_monitoring = histoDir.make<TH1F>("recoAK8HT_nominalHT410_monitoring", "recoAK8HT_nominalHT410_monitoring", 14000, 0, 14000);
  h1_recoAK8HT_nominalHT410            = histoDir.make<TH1F>("recoAK8HT_nominalHT410", "recoAK8HT_nominalHT410", 14000, 0, 14000);
  h1_recoAK8HT_monitoring              = histoDir.make<TH1F>("recoAK8HT_monitoring", "recoAK8HT_monitoring", 14000, 0, 14000);

  h1_recoAK8M_nominalHT250_monitoring = histoDir.make<TH1F>("recoAK8M_nominalHT250_monitoring", "recoAK8M_nominalHT250_monitoring", 40, 0, 250);
  h1_recoAK8M_nominalHT250            = histoDir.make<TH1F>("recoAK8M_nominalHT250", "recoAK8M_nominalHT250", 40, 0, 250);
  h1_recoAK8M_nominalHT410_monitoring = histoDir.make<TH1F>("recoAK8M_nominalHT410_monitoring", "recoAK8M_nominalHT410_monitoring", 40, 0, 250);
  h1_recoAK8M_nominalHT410            = histoDir.make<TH1F>("recoAK8M_nominalHT410", "recoAK8M_nominalHT410", 40, 0, 250);
  h1_recoAK8M_monitoring              = histoDir.make<TH1F>("recoAK8M_monitoring", "recoAK8M_monitoring", 40, 0, 250);

  h1_recoAK8Msd_nominalHT250_monitoring = histoDir.make<TH1F>("recoAK8Msd_nominalHT250_monitoring", "recoAK8Msd_nominalHT250_monitoring", 40, 0, 250);
  h1_recoAK8Msd_nominalHT250            = histoDir.make<TH1F>("recoAK8Msd_nominalHT250", "recoAK8Msd_nominalHT250", 40, 0, 250);
  h1_recoAK8Msd_nominalHT410_monitoring = histoDir.make<TH1F>("recoAK8Msd_nominalHT410_monitoring", "recoAK8Msd_nominalHT410_monitoring", 40, 0, 250);
  h1_recoAK8Msd_nominalHT410            = histoDir.make<TH1F>("recoAK8Msd_nominalHT410", "recoAK8Msd_nominalHT410", 40, 0, 250);
  h1_recoAK8Msd_monitoring              = histoDir.make<TH1F>("recoAK8Msd_monitoring", "recoAK8Msd_monitoring", 40, 0, 250);

  h1_recoAK8N2sdb1_nominalHT250_monitoring = histoDir.make<TH1F>("recoAK8N2sdb1_nominalHT250_monitoring", "recoAK8N2sdb1_nominalHT250_monitoring", 25,0,0.5);
  h1_recoAK8N2sdb1_nominalHT250            = histoDir.make<TH1F>("recoAK8N2sdb1_nominalHT250", "recoAK8N2sdb1_nominalHT250", 25,0,0.5);
  h1_recoAK8N2sdb1_nominalHT410_monitoring = histoDir.make<TH1F>("recoAK8N2sdb1_nominalHT410_monitoring", "recoAK8N2sdb1_nominalHT410_monitoring", 25,0,0.5);
  h1_recoAK8N2sdb1_nominalHT410            = histoDir.make<TH1F>("recoAK8N2sdb1_nominalHT410", "recoAK8N2sdb1_nominalHT410", 25,0,0.5);
  h1_recoAK8N2sdb1_monitoring              = histoDir.make<TH1F>("recoAK8N2sdb1_monitoring", "recoAK8N2sdb1_monitoring", 25,0,0.5);

  h1_recoAK8Pt_nominalHT250_monitoring = histoDir.make<TH1F>("recoAK8Pt_nominalHT250_monitoring", "recoAK8Pt_nominalHT250_monitoring", 20, 200, 700);
  h1_recoAK8Pt_nominalHT250            = histoDir.make<TH1F>("recoAK8Pt_nominalHT250", "recoAK8Pt_nominalHT250", 20, 200, 700);
  h1_recoAK8Pt_nominalHT410_monitoring = histoDir.make<TH1F>("recoAK8Pt_nominalHT410_monitoring", "recoAK8Pt_nominalHT410_monitoring", 20, 200, 700);
  h1_recoAK8Pt_nominalHT410            = histoDir.make<TH1F>("recoAK8Pt_nominalHT410", "recoAK8Pt_nominalHT410", 20, 200, 700);
  h1_recoAK8Pt_monitoring              = histoDir.make<TH1F>("recoAK8Pt_monitoring", "recoAK8Pt_monitoring", 20, 200, 700);

  h2_recoAK8MsdPt_nominalHT250_monitoring = histoDir.make<TH2F>("recoAK8MsdPt_nominalHT250_monitoring", "recoAK8MsdPt_nominalHT250_monitoring", 40, 0, 250, 20, 200, 700);
  h2_recoAK8MsdPt_nominalHT250            = histoDir.make<TH2F>("recoAK8MsdPt_nominalHT250", "recoAK8MsdPt_nominalHT250", 40, 0, 250, 20, 200, 700);
  h2_recoAK8MsdPt_nominalHT410_monitoring = histoDir.make<TH2F>("recoAK8MsdPt_nominalHT410_monitoring", "recoAK8MsdPt_nominalHT410_monitoring", 40, 0, 250, 20, 200, 700);
  h2_recoAK8MsdPt_nominalHT410            = histoDir.make<TH2F>("recoAK8MsdPt_nominalHT410", "recoAK8MsdPt_nominalHT410", 40, 0, 250, 20, 200, 700);
  h2_recoAK8MsdPt_monitoring              = histoDir.make<TH2F>("recoAK8MsdPt_monitoring", "recoAK8MsdPt_monitoring", 40, 0, 250, 20, 200, 700);
  
  h1_recoAK8M_monitoring_passTT      = histoDir.make<TH1F>("recoAK8M_monitoring_passTT", "recoAK8M_monitoring_passTT", 40, 0, 250);
  h1_recoAK8Msd_monitoring_passTT    = histoDir.make<TH1F>("recoAK8Msd_monitoring_passTT", "recoAK8Msd_monitoring_passTT", 40, 0, 250);
  h1_recoAK8N2sdb1_monitoring_passTT = histoDir.make<TH1F>("recoAK8N2sdb1_monitoring_passTT", "recoAK8N2sdb1_monitoring_passTT", 25, 0, 0.5);
  h1_recoAK8Pt_monitoring_passTT     = histoDir.make<TH1F>("recoAK8Pt_monitoring_passTT", "recoAK8Pt_monitoring_passTT", 20, 200, 700);
  h1_recoAK8M_monitoring_passTTN2p      = histoDir.make<TH1F>("recoAK8M_monitoring_passTTN2p", "recoAK8M_monitoring_passTTN2p", 40, 0, 250);
  h1_recoAK8Msd_monitoring_passTTN2p    = histoDir.make<TH1F>("recoAK8Msd_monitoring_passTTN2p", "recoAK8Msd_monitoring_passTTN2p", 40, 0, 250);
  h1_recoAK8N2sdb1_monitoring_passTTN2p = histoDir.make<TH1F>("recoAK8N2sdb1_monitoring_passTTN2p", "recoAK8N2sdb1_monitoring_passTTN2p", 25, 0, 0.5);
  h1_recoAK8Pt_monitoring_passTTN2p     = histoDir.make<TH1F>("recoAK8Pt_monitoring_passTTN2p", "recoAK8Pt_monitoring_passTTN2p", 20, 200, 700);
  h1_recoAK8M_monitoring_passTTN2f      = histoDir.make<TH1F>("recoAK8M_monitoring_passTTN2f", "recoAK8M_monitoring_passTTN2f", 40, 0, 250);
  h1_recoAK8Msd_monitoring_passTTN2f    = histoDir.make<TH1F>("recoAK8Msd_monitoring_passTTN2f", "recoAK8Msd_monitoring_passTTN2f", 40, 0, 250);
  h1_recoAK8N2sdb1_monitoring_passTTN2f = histoDir.make<TH1F>("recoAK8N2sdb1_monitoring_passTTN2f", "recoAK8N2sdb1_monitoring_passTTN2f", 25, 0, 0.5);
  h1_recoAK8Pt_monitoring_passTTN2f     = histoDir.make<TH1F>("recoAK8Pt_monitoring_passTTN2f", "recoAK8Pt_monitoring_passTTN2f", 20, 200, 700);

  /*
  if (doJECs_) {
    L1ParAK4_DATA = new JetCorrectorParameters(L1corrAK4_DATA_.fullPath());
    L2ParAK4_DATA = new JetCorrectorParameters(L2corrAK4_DATA_.fullPath());
    L3ParAK4_DATA = new JetCorrectorParameters(L3corrAK4_DATA_.fullPath());

    std::vector<JetCorrectorParameters> vParAK4_DATA;
    vParAK4_DATA.push_back(*L1ParAK4_DATA);
    vParAK4_DATA.push_back(*L2ParAK4_DATA);
    vParAK4_DATA.push_back(*L3ParAK4_DATA);

    JetCorrectorAK4_DATA = new FactorizedJetCorrector(vParAK4_DATA);
  }
  */
}



HTScoutingAnalyzer::~HTScoutingAnalyzer()
{
}

std::vector<math::XYZTLorentzVector> HTScoutingAnalyzer::makeP4s(const std::vector<fastjet::PseudoJet> &jets) {
  std::vector<math::XYZTLorentzVector> JetObjectsAll;
  for (const fastjet::PseudoJet & pj : jets) {
    JetObjectsAll.push_back( math::XYZTLorentzVector( pj.px(), pj.py(), pj.pz(), pj.e() ) );
  }
  return JetObjectsAll;
}

// ------------ method called for each event  ------------
void
HTScoutingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   std::cout << "\nEVT" << std::endl;
   using namespace edm;
   passpfTT=99;
   passpfTTN2p=99;
   passpfTTN2f=99;
   passrecoTT=99;
   passrecoTTN2p=99;
   passrecoTTN2f=99;
   passrecoTTN2pl=99;
   passrecoTTN2fl=99;
   passNominalHT250Trig=99;
   passNominalHT410Trig=99;
   passMonitoringTrig=99;
   caloHT=0.0;
   caloMjj=-1;
   caloDeltaEtajj = -1;
   caloMjjWide = -1; 
   caloDeltaEtajjWide = -1;
   caloDeltaPhijjWide = -1;
   
   pfmuPt = -1;
   pfmuPtSel = -1;

   recomuPt = -1;
   recomuPtSel = -1;

   pfMet = -1;
   recoMet = -1;

   pfHT=0.0;
   pfPt=-1;
   pfMjj=-1;
   pfDeltaEtajj = -1;
   pfMjjWide = -1; 
   pfMsdWide = -1;
   pfDeltaEtajjWide = -1;
   pfDeltaPhijjWide = -1;
   
   recoHT=0.0;
   recoPt=-1;
   recoMjj=-1;
   recoDeltaEtajj = -1;
   recoMjjWide = -1; 
   recoDeltaEtajjWide = -1;
   recoDeltaPhijjWide = -1;
   recocsv = -1;

   recoAK8HT=0.0;
   recoAK8M = -1;
   recoAK8Msd = -1;
   recoAK8N2sdb1 = -1;
   recoAK8Pt = -1;

   pfAK8HT=0.0;
   pfAK8M=-1;
   pfAK8Msd=-1;
   pfAK8N2sdb1=-1;
   pfAK8Pt= -1;
   pfAK8dRMu = -1;

   edm::Handle<edm::TriggerResults> trgResultsHandle;
   iEvent.getByToken(trgResultsLabel_, trgResultsHandle);
   
   const edm::TriggerNames &trgNames = iEvent.triggerNames(*trgResultsHandle);

   for (size_t i = 0; i < trgNames.size(); ++i) {
     const std::string &name = trgNames.triggerName(i);
   
     if ( (name.find("DST_HT250_CaloScouting_v") != std::string::npos )) {
       passNominalHT250Trig=trgResultsHandle->accept(i);
     }
     if ( (name.find("DST_HT410_PFScouting_v") != std::string::npos )) {
       passNominalHT410Trig=trgResultsHandle->accept(i);
     }
     if ( (name.find("HLT_Mu50_v") != std::string::npos )) {
       passMonitoringTrig=trgResultsHandle->accept(i);
     }
     
   }
   
   //std::cout << "HT250:" << passNominalHT250Trig << " Mu50:" << passMonitoringTrig << std::endl;
   //std::cout << "HT410:" << passNominalHT410Trig << " Mu50:" << passMonitoringTrig << std::endl;

   EnergyCorrelations* fECF;
   fECF = new EnergyCorrelations();

   //edm::Handle<ScoutingVertexCollection> pfVertexHandle;
   //iEvent.getByToken(pfVertexLabel_, pfVertexHandle);
   
   // But this selection is not done in data
   //GH = Handle("vector<reco::GenParticle>")
   //GL = ("prunedGenParticles", "")

   edm::Handle<ScoutingMuonCollection> pfMuonHandle;
   iEvent.getByToken(pfMuonLabel_, pfMuonHandle);
   
   edm::Handle<reco::VertexCollection> recoVertexHandle;
   iEvent.getByToken(recoVertexLabel_, recoVertexHandle);
   if (!recoVertexHandle.isValid()) {
     throw Exception(errors::ProductNotFound)
       << "Could not find reco::VertexCollection." << endl;
     return;
   }
   edm::Handle<pat::MuonCollection> recoMuonHandle;
   iEvent.getByToken(recoMuonLabel_, recoMuonHandle);

   edm::Handle<ScoutingCaloJetCollection> caloJetHandle;
   iEvent.getByToken(caloJetLabel_, caloJetHandle);
   /*
   edm::Handle<double> caloRhoHandle;
   iEvent.getByToken(caloRhoLabel_, caloRhoHandle);
   caloRho = *caloRhoHandle;
   */
   edm::Handle<double> pfMetHandle;
   iEvent.getByToken(pfMetLabel_, pfMetHandle);
   pfMet = *pfMetHandle;

   edm::Handle<pat::METCollection> recoMetHandle;
   iEvent.getByToken(recoMetLabel_, recoMetHandle);
   recoMet = (*recoMetHandle)[0].et();

   edm::Handle<ScoutingPFJetCollection> pfJetHandle;
   iEvent.getByToken(pfJetLabel_, pfJetHandle);

   edm::Handle<pat::JetCollection> recoJetHandle;
   iEvent.getByToken(recoJetLabel_, recoJetHandle);
   
   edm::Handle<pat::JetCollection> recoAK8JetHandle;
   iEvent.getByToken(recoAK8JetLabel_, recoAK8JetHandle);

   //edm::Handle<reco::BasicJetCollection> recoAK8SoftDropJetHandle;
   //iEvent.getByToken(recoAK8SoftDropJetLabel_, recoAK8SoftDropJetHandle);

   edm::Handle<ScoutingParticleCollection> particleHandle;
   iEvent.getByToken(particleLabel_, particleHandle);
   
   //JEC factors
   std::vector<double> jecFactorsAK4;
   // Sort AK4 jets by increasing pT
   std::vector<unsigned> sortedAK4JetIdx;
   std::multimap<double, unsigned> sortedAK4Jets;

   for (ScoutingCaloJetCollection::const_iterator iCj = caloJetHandle->begin(); iCj != caloJetHandle->end(); ++iCj) {	   
     TLorentzVector cj1;
     
     double correction = 1.0;
     /*
     if (doJECs_) {
       JetCorrectorAK4_DATA->setJetEta(iCj->eta());
       JetCorrectorAK4_DATA->setJetPt(iCj->pt());
       JetCorrectorAK4_DATA->setJetA(iCj->jetArea());
       JetCorrectorAK4_DATA->setRho(caloRho);
       correction = JetCorrectorAK4_DATA->getCorrection();
     }
     */
     jecFactorsAK4.push_back(correction);
     cj1.SetPtEtaPhiM(iCj->pt()*correction, iCj->eta(), iCj->phi(), iCj->m());
     sortedAK4Jets.insert(std::make_pair(cj1.Pt(), iCj - caloJetHandle->begin()));
     if (cj1.Pt() > 40. && fabs(cj1.Eta()) < 3.0) caloHT += cj1.Pt();
   }
   // Get jet indices in decreasing pT order
   for (std::multimap<double, unsigned>::const_reverse_iterator it=sortedAK4Jets.rbegin(); it!=sortedAK4Jets.rend(); ++it) {
     sortedAK4JetIdx.push_back(it->second);
   }
   

   for (ScoutingPFJetCollection::const_iterator iCj = pfJetHandle->begin(); iCj != pfJetHandle->end(); ++iCj) {	   
     TLorentzVector cj1;
     cj1.SetPtEtaPhiM(iCj->pt(), iCj->eta(), iCj->phi(), iCj->m());
     if (iCj->pt() > 40. && fabs(iCj->eta()) < 3.0) pfHT += iCj->pt();
   }

   for (pat::JetCollection::const_iterator iCj = recoJetHandle->begin(); iCj != recoJetHandle->end(); ++iCj) {	   
     TLorentzVector cj1;
     cj1.SetPtEtaPhiM(iCj->pt(), iCj->eta(), iCj->phi(), iCj->mass());
     if (iCj->pt() > 40. && fabs(iCj->eta()) < 3.0) recoHT += iCj->pt();
   }

   for (pat::JetCollection::const_iterator iCj = recoAK8JetHandle->begin(); iCj != recoAK8JetHandle->end(); ++iCj) {
     TLorentzVector cj1;
     cj1.SetPtEtaPhiM(iCj->pt(), iCj->eta(), iCj->phi(), iCj->mass());
     if (iCj->pt() > 40. && fabs(iCj->eta()) < 3.0) recoAK8HT += iCj->pt();
   }
      
   double wideJetDeltaR_ = 1.1;

   TLorentzVector pfmuon;
   // PFMuon Scouting
   if (pfMuonHandle->size() > 1){
     for(ScoutingMuonCollection::const_iterator itMu = pfMuonHandle->begin(); itMu!=pfMuonHandle->end(); ++itMu) {
       if(itMu->pt() < pfmuPt || itMu->pt() < 10) continue;
       if(fabs(itMu->eta()) >= 2.4) continue;
       pfmuPt = itMu->pt();
       if(itMu->pt() <= 53) continue;
       if(fabs(itMu->eta() >= 2.1)) continue;
       // sort of equivalent cuts to tight muon except isPFMuon
       if(!itMu->isGlobalMuon()) continue;
       if(itMu->nTrackerLayersWithMeasurement() <= 5 || itMu->nValidPixelHits() <= 0) continue;
       if(itMu->nMatchedStations() <= 1) continue;
       if(itMu->chi2() >= 10) continue;
       if(fabs(itMu->dxy()) >= 0.2 || fabs(itMu->dz()) >=0.5) continue;
       // tracker isolation 
       if(itMu->trackIso() >= 0.05*(itMu->pt())) continue;
       pfmuPtSel = itMu->pt();
       pfmuon.SetPtEtaPhiM(itMu->pt(),itMu->eta(),itMu->phi(),0.105658369);
     }
   }

   // RecoVertex
   const reco::VertexCollection *pvCol = recoVertexHandle.product();
   const reco::Vertex* pv = &(*pvCol->begin());
   int nvtx = 0;
   for(reco::VertexCollection::const_iterator itVtx = pvCol->begin(); itVtx!=pvCol->end(); ++itVtx) {    
     if(itVtx->chi2()==0 && itVtx->ndof()==0) continue; 
     //if(itVtx->tracksSize()     < 0) continue;
     if(itVtx->ndof()           < 4) continue;
     if(fabs(itVtx->z())        > 24) continue;
     if(itVtx->position().Rho() > 2) continue;
     // vertices are sorted by sum{pT^2}, so the first one passing cuts
     // is taken as the event primary vertex
     if(nvtx==0) {
       pv = &(*itVtx);
     }
     nvtx++;
   }
   // RecoMuon 
   TLorentzVector recomuon;
   if (recoMuonHandle->size() > 1){
     for(pat::MuonCollection::const_iterator itMu = recoMuonHandle->begin(); itMu!=recoMuonHandle->end(); ++itMu) {
       if(itMu->muonBestTrack()->pt() < recomuPt || itMu->muonBestTrack()->pt()<10) continue;
       if(fabs(itMu->muonBestTrack()->eta()) >=  2.4) continue;
       recomuPt = itMu->muonBestTrack()->pt();
       if(itMu->muonBestTrack()->pt() <= 53) continue;
       if(fabs(itMu->muonBestTrack()->eta() >= 2.1)) continue;
       if(!muon::isTightMuon( *itMu, *pv ) ) continue;
       double iso = itMu->pfIsolationR04().sumChargedHadronPt + TMath::Max(itMu->pfIsolationR04().sumNeutralHadronEt + itMu->pfIsolationR04().sumPhotonEt - 0.5*(itMu->pfIsolationR04().sumPUPt), double(0));
       if(iso >= 0.15*(itMu->muonBestTrack()->pt())) continue;
       recomuPtSel = itMu->muonBestTrack()->pt();
       recomuon.SetPtEtaPhiM(itMu->muonBestTrack()->pt(),itMu->muonBestTrack()->eta(),itMu->muonBestTrack()->phi(),0.105658369);
     }
   }

   // PFParticles Scouting
   // Get PF particles and clusterize in AK8 jets
   TLorentzVector pfak8jet;
   std::vector<fastjet::PseudoJet>  lClusterParticles;
   if(particleHandle->size() > 1){
     for (ScoutingParticleCollection::const_iterator kPj = particleHandle->begin(); kPj != particleHandle->end(); ++kPj) {
       TLorentzVector vC0; vC0.SetPtEtaPhiM(kPj->pt(),kPj->eta(),kPj->phi(),kPj->m());
       fastjet::PseudoJet pPart(vC0.Px(),vC0.Py(),vC0.Pz(),vC0.E());
       lClusterParticles.emplace_back(pPart);
     }
   }
   // Choose a Jet Definition and Sequence, (0) is ptMin
   fastjet::JetDefinition lCJet_def8(fastjet::antikt_algorithm, 0.8);
   fastjet::ClusterSequence lCClust_seq8(lClusterParticles, lCJet_def8);
   std::vector<fastjet::PseudoJet> ak8inclusive_jets = fastjet::sorted_by_pt(lCClust_seq8.inclusive_jets(0));
   std::vector<math::XYZTLorentzVector> ak8jets = makeP4s(ak8inclusive_jets);
   for(unsigned int i0 = 0; i0 < ak8inclusive_jets.size(); i0++) {
     fastjet::contrib::SoftDrop SD(0.,0.1,0.8);
     fastjet::PseudoJet SD_jet = SD(ak8inclusive_jets[i0]);
     if (ak8jets[i0].Pt() > 100. && fabs(ak8jets[i0].Eta()) < 3.0) pfAK8HT += ak8jets[i0].Pt();
     if(i0==0) {
       pfAK8M = ak8jets[i0].M();
       pfAK8Msd = SD_jet.m();
       pfAK8Pt = ak8jets[i0].Pt();
     }
     double beta=1;
     std::vector<fastjet::PseudoJet> lSDClusterParticles = SD_jet.constituents();
     std::sort(lSDClusterParticles.begin(),lSDClusterParticles.end(),orderPseudoJet);
     int nFilter = TMath::Min(100,(int)lSDClusterParticles.size());
     std::vector<fastjet::PseudoJet> lSDFilter(lSDClusterParticles.begin(),lSDClusterParticles.begin()+nFilter);
     fECF->calcECFN(beta,lSDFilter,true);
     float ak8pfe2_sdb1      = float(fECF->manager->ecfns["2_2"]);
     float ak8pfe3_v2_sdb1   = float(fECF->manager->ecfns["3_2"]);
     if(i0==0) {
       pfAK8N2sdb1 = ak8pfe3_v2_sdb1/(ak8pfe2_sdb1*ak8pfe2_sdb1);
     }
     if(i0==0 && pfmuPtSel > 53) {
       pfak8jet.SetPtEtaPhiM(ak8jets[i0].Pt(),ak8jets[i0].Eta(),ak8jets[i0].Phi(),SD_jet.m());
       pfAK8dRMu = pfak8jet.DeltaR(pfmuon);
     }
     //std::cout <<"pt " << ak8jets[i0].Pt() << " mass " << ak8jets[i0].M() << " Msd " << SD_jet.m() << " HT " << pfAK8HT << std::endl;
   }
   //std::cout << "n ak8 jets " << ak8inclusive_jets.size() << std::endl;
   //std::cout <<"ak8 pt " << pfAK8Pt << " mass " << pfAK8M << " Msd " << pfAK8Msd << " HT " << pfAK8HT << std::endl;

   // CaloScoutingAK4Jets and Wide
   if (caloJetHandle->size() > 2){     
     TLorentzVector wj1, wj2, wdijet;
     TLorentzVector wj1_tmp, wj2_tmp;
     const ScoutingCaloJet & iCj = (*caloJetHandle)[ sortedAK4JetIdx[0] ];
     const ScoutingCaloJet & jCj = (*caloJetHandle)[ sortedAK4JetIdx[1] ];
     TLorentzVector jet1, jet2;
     jet1.SetPtEtaPhiM(jecFactorsAK4[ sortedAK4JetIdx[0] ]*(iCj.pt()), iCj.eta(), iCj.phi(), iCj.m());
     jet2.SetPtEtaPhiM(jecFactorsAK4[ sortedAK4JetIdx[1] ]*(jCj.pt()), jCj.eta(), jCj.phi(), jCj.m());
     if (jet1.Pt() > 60. && fabs(jet1.Eta()) < 2.5) {	 
       if (jet2.Pt() > 30. && fabs(jet2.Eta()) < 2.5) {
	 caloMjj = (jet1+jet2).M();
	 caloDeltaEtajj = fabs(jet1.Eta()-jet2.Eta());
	 for(size_t ijet=0; ijet<caloJetHandle->size(); ijet++) {
	   const ScoutingCaloJet & kCj = (*caloJetHandle)[ sortedAK4JetIdx[ijet] ];
	   TLorentzVector currentJet;
	   currentJet.SetPtEtaPhiM(jecFactorsAK4[ sortedAK4JetIdx[ijet] ]*(kCj.pt()), kCj.eta(), kCj.phi(), kCj.m());	     
	   if (currentJet.Pt() > 30. && fabs(currentJet.Eta()) < 2.5) {
	     double DeltaR1 = currentJet.DeltaR(jet1);
	     double DeltaR2 = currentJet.DeltaR(jet2);			   
	     if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_)
	       {
		 wj1_tmp += currentJet;
	       }
	     else if(DeltaR2 < wideJetDeltaR_)
	       {
		 wj2_tmp += currentJet;
	       }			 
	   } // if AK4 jet passes fid and jetid.
	 } //end of ak4 jet loop		     
       } //fid, jet id, pt cut
     } //fid, jet id, pt cut
     
     // Re-order the wide jets in pt
     //std::cout << "wide " << wj1_tmp.Pt() << wj2_tmp.Pt() <<std::endl;
     if( wj1_tmp.Pt() > wj2_tmp.Pt())
       {
	 wj1 = wj1_tmp;
	 wj2 = wj2_tmp;
       }
     else
       {
	 wj1 = wj2_tmp;
	 wj2 = wj1_tmp;
       }
   
     if( wj1.Pt()>0 && wj2.Pt()>0 )
       {
	 // Create dijet system
	 wdijet = wj1 + wj2;
	 caloMjjWide = wdijet.M();
	 caloDeltaEtajjWide = fabs(wj1.Eta()-wj2.Eta());
	 caloDeltaPhijjWide = fabs(wj1.DeltaPhi(wj2));
       }

   } // end of two calo jets.

   // PfAK4 Scouting Jets and wide
   TLorentzVector pfak4jet;
   for (ScoutingPFJetCollection::const_iterator iCj = pfJetHandle->begin(); iCj != pfJetHandle->end(); ++iCj) {
     if(iCj->pt()>30. && fabs(iCj->eta()) < 2.4){
       pfPt = iCj->pt();
       pfcsv = iCj->csv();
       if(pfPt>30. && pfcsv>0.8484) {
	 pfak4jet.SetPtEtaPhiM(iCj->pt(),iCj->eta(),iCj->phi(),iCj->m());
	 if(pfmuPtSel > 53) {
	   pfdRMu = pfak4jet.DeltaR(pfmuon);
	 }
	 if(pfak8jet.Pt()>200) {
	   pfdRAK8 = pfak4jet.DeltaR(pfak8jet);
	 }
       }
       break;
     }
   }
   if (pfJetHandle->size() > 2){     
     TLorentzVector wj1, wj2, wdijet;
     TLorentzVector wj1_tmp, wj2_tmp;
     const ScoutingPFJet & iCj = (*pfJetHandle)[ 0 ];
     const ScoutingPFJet & jCj = (*pfJetHandle)[ 1 ];
     std::vector<fastjet::PseudoJet>  lClusterParticlesW;                                                                                                 
     for(unsigned int i0 = 0; i0 < iCj.constituents().size(); i0++) {
       const ScoutingParticle &iC0 = (*particleHandle)[i0];
       TLorentzVector vC0; vC0.SetPtEtaPhiM(iC0.pt(),iC0.eta(),iC0.phi(),iC0.m());
       fastjet::PseudoJet pPart(vC0.Px(),vC0.Py(),vC0.Pz(),vC0.E());
       lClusterParticlesW.emplace_back(pPart);  
     }

     if (iCj.pt() > 60. && fabs(iCj.eta()) < 2.5) {	 
       if (jCj.pt() > 30. && fabs(jCj.eta()) < 2.5) {
	 TLorentzVector jet1, jet2;
	 jet1.SetPtEtaPhiM(iCj.pt(), iCj.eta(), iCj.phi(), iCj.m());
	 jet2.SetPtEtaPhiM(jCj.pt(), jCj.eta(), jCj.phi(), jCj.m());
	 pfMjj = (jet1+jet2).M();
	 pfDeltaEtajj = fabs(jet1.Eta()-jet2.Eta());
	 for (ScoutingPFJetCollection::const_iterator kCj = pfJetHandle->begin(); kCj != pfJetHandle->end(); ++kCj) {
	   TLorentzVector currentJet;
	   std::vector<int> pfConstituents = kCj->constituents();
	   if (kCj->pt() > 30. && fabs(kCj->eta()) < 2.5) {
	     currentJet.SetPtEtaPhiM(kCj->pt(), kCj->eta(), kCj->phi(), kCj->m());			   
	     double DeltaR1 = currentJet.DeltaR(jet1);
	     double DeltaR2 = currentJet.DeltaR(jet2);			   
	     if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_)
	       {
		 wj1_tmp += currentJet;
		 for(unsigned int i0 = 0; i0 < pfConstituents.size(); i0++) {
		   const ScoutingParticle &kC0 = (*particleHandle)[i0];
		   TLorentzVector vC0; vC0.SetPtEtaPhiM(kC0.pt(),kC0.eta(),kC0.phi(),kC0.m());
		   fastjet::PseudoJet pPart(vC0.Px(),vC0.Py(),vC0.Pz(),vC0.E());
		   lClusterParticlesW.emplace_back(pPart);
		 }
	       }
	     else if(DeltaR2 < wideJetDeltaR_)
	       {
		 wj2_tmp += currentJet;
	       }			 
	   } // if AK4 jet passes fid and jetid.
	 } //end of ak4 jet loop		     
	
       } //fid, jet id, pt cut
     } //fid, jet id, pt cut

     // Recluster constituents of wide jet
     // fastjet::JetDefinition lCJet_def(fastjet::cambridge_algorithm, 0.8);
     // fastjet::ClusterSequence lCClust_seq(lClusterParticlesW, lCJet_def);
     // std::vector<fastjet::PseudoJet> inclusive_jets = lCClust_seq.inclusive_jets(0);
     // fastjet::contrib::SoftDrop sd(0.,0.1,0.8);
     // fastjet::PseudoJet sd_jet = sd(inclusive_jets[0]);
     // pfMsdWide = sd_jet.m();
     // double beta=1;
     // std::vector<fastjet::PseudoJet> lSDClusterParticlesW = sd_jet.constituents();
     // std::sort(lSDClusterParticlesW.begin(),lSDClusterParticlesW.end(),orderPseudoJet);
     // int nFilter = TMath::Min(100,(int)lSDClusterParticlesW.size());
     // std::vector<fastjet::PseudoJet> lSDFiltered(lSDClusterParticlesW.begin(),lSDClusterParticlesW.begin()+nFilter);
     // fECF->calcECFN(beta,lSDFiltered,true);
     // float pfe2_sdb1      = float(fECF->manager->ecfns["2_2"]);
     // //float pfe3_sdb1      = float(fECF->manager->ecfns["3_3"]);
     // //float pfe3_v1_sdb1   = float(fECF->manager->ecfns["3_1"]);
     // float pfe3_v2_sdb1   = float(fECF->manager->ecfns["3_2"]);
     // pfN2sdb1 = pfe3_v2_sdb1/(pfe2_sdb1*pfe2_sdb1);
     
     // Re-order the wide jets in pt
     if( wj1_tmp.Pt() > wj2_tmp.Pt())
       {
	 wj1 = wj1_tmp;
	 wj2 = wj2_tmp;
       }
     else
       {
	 wj1 = wj2_tmp;
	 wj2 = wj1_tmp;
       }
   
     if( wj1.Pt()>0 && wj2.Pt()>0 )
       {
	 // Create dijet system
	 wdijet = wj1 + wj2;
	 pfMjjWide = wdijet.M();
	 //std::cout << "mass " << wj1.M() << " inclusivejets reclustered " << inclusive_jets[0].m() << " size " << inclusive_jets.size() << " sd " << pfMsdWide << std::endl;
	 pfDeltaEtajjWide = fabs(wj1.Eta()-wj2.Eta());
	 pfDeltaPhijjWide = fabs(wj1.DeltaPhi(wj2));
       }


   } // end of two PF jets.

   // Reco AK4 jets
   TLorentzVector recoak4jet;
   for (pat::JetCollection::const_iterator iCj = recoJetHandle->begin(); iCj != recoJetHandle->end(); ++iCj) {
     if(iCj->pt()>30. && fabs(iCj->eta()) < 2.4){
       recoPt = iCj->pt();
       recocsv = iCj->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
       if(recoPt>30. && recocsv>0.8484) {
	 recoak4jet.SetPtEtaPhiM(iCj->pt(),iCj->eta(),iCj->phi(),iCj->mass());
         if(recomuPtSel > 53) {
           recodRMu = recoak4jet.DeltaR(recomuon);
         }
	 break;
       }
     }
   }
   if (recoJetHandle->size() > 2){     
     TLorentzVector wj1, wj2, wdijet;
     TLorentzVector wj1_tmp, wj2_tmp;
     const pat::Jet & iCj = (*recoJetHandle)[ 0 ];
     const pat::Jet & jCj = (*recoJetHandle)[ 1 ];
     if (iCj.pt() > 60. && fabs(iCj.eta()) < 2.5) {	 
       if (jCj.pt() > 30. && fabs(jCj.eta()) < 2.5) {
	 TLorentzVector jet1, jet2;
	 jet1.SetPtEtaPhiM(iCj.pt(), iCj.eta(), iCj.phi(), iCj.mass());
	 jet2.SetPtEtaPhiM(jCj.pt(), jCj.eta(), jCj.phi(), jCj.mass());
	 recoMjj = (jet1+jet2).M();
	 recoDeltaEtajj = fabs(jet1.Eta()-jet2.Eta());
	 for (pat::JetCollection::const_iterator kCj = recoJetHandle->begin(); kCj != recoJetHandle->end(); ++kCj) {
	   TLorentzVector currentJet;
	   if (kCj->pt() > 30. && fabs(kCj->eta()) < 2.5) {
	     currentJet.SetPtEtaPhiM(kCj->pt(), kCj->eta(), kCj->phi(), kCj->mass());			   
	     double DeltaR1 = currentJet.DeltaR(jet1);
	     double DeltaR2 = currentJet.DeltaR(jet2);			   
	     if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_)
	       {
		 wj1_tmp += currentJet;
	       }
	     else if(DeltaR2 < wideJetDeltaR_)
	       {
		 wj2_tmp += currentJet;
	       }			 
	   } // if AK4 jet passes fid and jetid.
	 } //end of ak4 jet loop		     
       } //fid, jet id, pt cut
     } //fid, jet id, pt cut
     
     // Re-order the wide jets in pt
     if( wj1_tmp.Pt() > wj2_tmp.Pt())
       {
	 wj1 = wj1_tmp;
	 wj2 = wj2_tmp;
       }
     else
       {
	 wj1 = wj2_tmp;
	 wj2 = wj1_tmp;
       }
   
     if( wj1.Pt()>0 && wj2.Pt()>0 )
       {
	 // Create dijet system
	 wdijet = wj1 + wj2;
	 recoMjjWide = wdijet.M();
	 recoDeltaEtajjWide = fabs(wj1.Eta()-wj2.Eta());
	 recoDeltaPhijjWide = fabs(wj1.DeltaPhi(wj2));
       }

   } // end of two RECO jets.
   
   // Reco AK8 Jets
   TLorentzVector recoak8jet;
   if (recoAK8JetHandle->size() > 0){
     TLorentzVector wj1;
     TLorentzVector wj1_tmp;
     //const reco::BasicJetCollection *softdropJetCol = recoAK8SoftDropJetHandle.product();
     const pat::Jet & iCj = (*recoAK8JetHandle)[ 0 ];
     //const reco::BasicJet* matchJet = 0;
     //matchJet = match(&iCj,softdropJetCol);
     //if(matchJet) {
     //  std::cout << "sd" << std::endl;
     //}
     std::vector<fastjet::PseudoJet>  lClusterParticles;
     std::vector<reco::CandidatePtr> pfConstituents = iCj.getJetConstituents();
     for(unsigned int ic=0; ic<pfConstituents.size(); ic++) {                                                                                  
       reco::CandidatePtr pfcand = pfConstituents[ic];
       fastjet::PseudoJet   pPart(pfcand->px(),pfcand->py(),pfcand->pz(),pfcand->energy());
       lClusterParticles.emplace_back(pPart);
     } 
     fastjet::JetDefinition lCJet_def(fastjet::antikt_algorithm, 0.8);
     fastjet::ClusterSequence lCClust_seq(lClusterParticles, lCJet_def);
     std::vector<fastjet::PseudoJet> inclusive_jets = lCClust_seq.inclusive_jets(0);
     fastjet::contrib::SoftDrop sd(0.,0.1,0.8);
     fastjet::PseudoJet sd_jet = sd(inclusive_jets[0]);
     recoAK8Pt = iCj.pt();
     recoAK8M = iCj.mass();
     recoAK8Msd = sd_jet.m();
     double beta=1;
     std::vector<fastjet::PseudoJet> lSDClusterParticles = sd_jet.constituents();
     std::sort(lSDClusterParticles.begin(),lSDClusterParticles.end(),orderPseudoJet);
     int nFilter = TMath::Min(100,(int)lSDClusterParticles.size());
     std::vector<fastjet::PseudoJet> lSDFiltered(lSDClusterParticles.begin(),lSDClusterParticles.begin()+nFilter);
     fECF->calcECFN(beta,lSDFiltered,true);
     float recoAK8e2_sdb1      = float(fECF->manager->ecfns["2_2"]);
     float recoAK8e3_v2_sdb1   = float(fECF->manager->ecfns["3_2"]);
     recoAK8N2sdb1 = recoAK8e3_v2_sdb1/(recoAK8e2_sdb1*recoAK8e2_sdb1);
     if(fabs(iCj.eta() < 2.4)){
       recoak8jet.SetPtEtaPhiM(iCj.pt(),iCj.eta(),iCj.phi(),recoAK8Msd);
     }
     if(recomuPtSel > 53) {
       recoAK8dRMu = recoak8jet.DeltaR(recomuon);
     }
     if(recoak8jet.Pt()>200) {
       recodRAK8 = recoak4jet.DeltaR(recoak8jet);
     }
   } // end of AK8 RECO jet. 

   //reco selecion
   std::cout << "recomuPtSel " << recomuPtSel << std::endl;
   std::cout << "recoMet " << recoMet << std::endl;
   std::cout << "recoAK8Pt " << recoak8jet.Pt() << std::endl;
   std::cout << "recoAK8Pt " << recoAK8Pt << std::endl;

   if(recomuPtSel>53&&recoMet>40&&recoak8jet.Pt()>200&&recoAK8dRMu>0.8&&recoak4jet.Pt()>30&&recodRMu>0.3&&recodRAK8>0.8) {
     std::cout << "is it passing" << std::endl;
     passrecoTT = 1;
     if(recoAK8N2sdb1 < 0.2) passrecoTTN2p=1;
     else passrecoTTN2f=1;
     if(recoAK8N2sdb1 < 0.3) passrecoTTN2pl=1;
     else passrecoTTN2fl=1;
   }
   // pf selection
   if(pfmuPtSel>53&&pfMet>40&&pfak8jet.Pt()>200&&pfAK8dRMu>0.8&&pfak4jet.Pt()>30&&pfdRMu>0.3&&pfdRAK8>0.8) {
     passpfTT= 1;
     if(pfAK8N2sdb1 < 0.2) passpfTTN2p=1;
     else passpfTTN2f=1;
   }

   std::cout << "filling hist" << std::endl;
   std::cout << "HT250" << std::endl;
   if (passNominalHT250Trig && passMonitoringTrig) {
     h1_pfAK8HT_nominalHT250_monitoring->Fill(pfAK8HT) ;
     if(pfAK8M > 0) h1_pfAK8M_nominalHT250_monitoring->Fill(pfAK8M) ;
     if(pfAK8Msd > 0) h1_pfAK8Msd_nominalHT250_monitoring->Fill(pfAK8Msd) ;
     if(pfAK8Pt > 0) {
       h1_pfAK8Pt_nominalHT250_monitoring->Fill(pfAK8Pt) ;
       h2_pfAK8MsdPt_nominalHT250_monitoring->Fill(pfAK8Msd,pfAK8Pt) ;
     }
     if(pfAK8N2sdb1 > 0) h1_pfAK8N2sdb1_nominalHT250_monitoring->Fill(pfAK8N2sdb1) ;
     
     h1_caloHT_nominalHT250_monitoring->Fill(caloHT) ;
     if (caloDeltaEtajj > -1 && caloDeltaEtajj < 1.3) h1_caloMjj_nominalHT250_monitoring->Fill(caloMjj) ;
     h1_caloDeltaEtajj_nominalHT250_monitoring->Fill(caloDeltaEtajj) ;
     if (caloDeltaEtajjWide > -1 && caloDeltaEtajjWide < 1.3) h1_caloMjjWide_nominalHT250_monitoring->Fill(caloMjjWide) ;
     h1_caloDeltaEtajjWide_nominalHT250_monitoring->Fill(caloDeltaEtajjWide) ;
     
     h1_pfHT_nominalHT250_monitoring->Fill(pfHT) ;
     if (pfDeltaEtajj > -1 && pfDeltaEtajj < 1.3) h1_pfMjj_nominalHT250_monitoring->Fill(pfMjj) ;
     h1_pfDeltaEtajj_nominalHT250_monitoring->Fill(pfDeltaEtajj) ;
     if (pfDeltaEtajjWide > -1 && pfDeltaEtajjWide < 1.3) h1_pfMjjWide_nominalHT250_monitoring->Fill(pfMjjWide) ;
     h1_pfDeltaEtajjWide_nominalHT250_monitoring->Fill(pfDeltaEtajjWide) ;
   
     h1_recoHT_nominalHT250_monitoring->Fill(recoHT) ;
     if (recoDeltaEtajj > -1 && recoDeltaEtajj < 1.3) h1_recoMjj_nominalHT250_monitoring->Fill(recoMjj) ;
     h1_recoDeltaEtajj_nominalHT250_monitoring->Fill(recoDeltaEtajj) ;
     if (recoDeltaEtajjWide > -1 && recoDeltaEtajjWide < 1.3) h1_recoMjjWide_nominalHT250_monitoring->Fill(recoMjjWide) ;
     h1_recoDeltaEtajjWide_nominalHT250_monitoring->Fill(recoDeltaEtajjWide) ;

     h1_recoAK8HT_nominalHT250_monitoring->Fill(recoAK8HT) ;
     if(recoAK8M > -1) h1_recoAK8M_nominalHT250_monitoring->Fill(recoAK8M) ;
     if(recoAK8Msd > 0) h1_recoAK8Msd_nominalHT250_monitoring->Fill(recoAK8Msd) ;
     if(recoAK8Pt > 0) {
       h1_recoAK8Pt_nominalHT250_monitoring->Fill(recoAK8Pt) ;
       h2_recoAK8MsdPt_nominalHT250_monitoring->Fill(recoAK8Msd,recoAK8Pt) ;
     }
     if(recoAK8N2sdb1 > 0) h1_recoAK8N2sdb1_nominalHT250_monitoring->Fill(recoAK8N2sdb1) ;
   }
   std::cout << "HT410"<< std::endl;

   if (passNominalHT410Trig && passMonitoringTrig) {
     h1_pfAK8HT_nominalHT410_monitoring->Fill(pfAK8HT) ;
     if(pfAK8M > 0) h1_pfAK8M_nominalHT410_monitoring->Fill(pfAK8M) ;
     if(pfAK8Msd > 0) h1_pfAK8Msd_nominalHT410_monitoring->Fill(pfAK8Msd) ;
     if(pfAK8Pt > 0) {
       h1_pfAK8Pt_nominalHT410_monitoring->Fill(pfAK8Pt) ;
       h2_pfAK8MsdPt_nominalHT410_monitoring->Fill(pfAK8Msd,pfAK8Pt) ;
     }
     if(pfAK8N2sdb1 > 0) h1_pfAK8N2sdb1_nominalHT410_monitoring->Fill(pfAK8N2sdb1) ;
     
     h1_caloHT_nominalHT410_monitoring->Fill(caloHT);
     if (caloDeltaEtajj > -1 && caloDeltaEtajj < 1.3) {
       h1_caloMjj_nominalHT410_monitoring->Fill(caloMjj) ;
     }
     h1_caloDeltaEtajj_nominalHT410_monitoring->Fill(caloDeltaEtajj) ;
     if (caloDeltaEtajjWide > -1 && caloDeltaEtajjWide < 1.3) h1_caloMjjWide_nominalHT410_monitoring->Fill(caloMjjWide) ;
     h1_caloDeltaEtajjWide_nominalHT410_monitoring->Fill(caloDeltaEtajjWide) ;

     h1_pfHT_nominalHT410_monitoring->Fill(pfHT) ;
     if (pfDeltaEtajj > -1 && pfDeltaEtajj < 1.3) h1_pfMjj_nominalHT410_monitoring->Fill(pfMjj) ;
     h1_pfDeltaEtajj_nominalHT410_monitoring->Fill(pfDeltaEtajj) ;
     if (pfDeltaEtajjWide > -1 && pfDeltaEtajjWide < 1.3) h1_pfMjjWide_nominalHT410_monitoring->Fill(pfMjjWide) ;
     h1_pfDeltaEtajjWide_nominalHT410_monitoring->Fill(pfDeltaEtajjWide) ;
 
     h1_recoHT_nominalHT410_monitoring->Fill(recoHT) ;
     if (recoDeltaEtajj > -1 && recoDeltaEtajj < 1.3) h1_recoMjj_nominalHT410_monitoring->Fill(recoMjj) ;
     h1_recoDeltaEtajj_nominalHT410_monitoring->Fill(recoDeltaEtajj) ;
     if (recoDeltaEtajjWide > -1 && recoDeltaEtajjWide < 1.3) h1_recoMjjWide_nominalHT410_monitoring->Fill(recoMjjWide) ;
     h1_recoDeltaEtajjWide_nominalHT410_monitoring->Fill(recoDeltaEtajjWide) ;
     
     h1_recoAK8HT_nominalHT410_monitoring->Fill(recoAK8HT) ;
     if(recoAK8M > 0) h1_recoAK8M_nominalHT410_monitoring->Fill(recoAK8M) ;
     if(recoAK8Msd > 0) h1_recoAK8Msd_nominalHT410_monitoring->Fill(recoAK8Msd) ;
     if(recoAK8Pt > 0) {
       h1_recoAK8Pt_nominalHT410_monitoring->Fill(recoAK8Pt) ;
       h2_recoAK8MsdPt_nominalHT410_monitoring->Fill(recoAK8Msd,recoAK8Pt) ;
     }
     if(recoAK8N2sdb1 > 0) h1_recoAK8N2sdb1_nominalHT410_monitoring->Fill(recoAK8N2sdb1) ;
   }
   std::cout << "HT250no"<< std::endl;

   if (passNominalHT250Trig) {
     h1_pfAK8HT_nominalHT250->Fill(pfAK8HT) ;
     if(pfAK8M > 0) h1_pfAK8M_nominalHT250->Fill(pfAK8M) ;
     if(pfAK8Msd > 0) h1_pfAK8Msd_nominalHT250->Fill(pfAK8Msd) ;
     if(pfAK8Pt > 0) {
       h1_pfAK8Pt_nominalHT250->Fill(pfAK8Pt) ;
       h2_pfAK8MsdPt_nominalHT250->Fill(pfAK8Msd,pfAK8Pt) ;
     }
     if(pfAK8N2sdb1 > 0) h1_pfAK8N2sdb1_nominalHT250->Fill(pfAK8N2sdb1) ;

     h1_caloHT_nominalHT250->Fill(caloHT) ;
     if (caloDeltaEtajj > -1 && caloDeltaEtajj < 1.3) h1_caloMjj_nominalHT250->Fill(caloMjj) ;
     h1_caloDeltaEtajj_nominalHT250->Fill(caloDeltaEtajj) ;
     if (caloDeltaEtajjWide > -1 && caloDeltaEtajjWide < 1.3) h1_caloMjjWide_nominalHT250->Fill(caloMjjWide) ;
     h1_caloDeltaEtajjWide_nominalHT250->Fill(caloDeltaEtajjWide) ;
     
     h1_pfHT_nominalHT250->Fill(pfHT) ;
     if (pfDeltaEtajj > -1 && pfDeltaEtajj < 1.3) h1_pfMjj_nominalHT250->Fill(pfMjj) ;
     h1_pfDeltaEtajj_nominalHT250->Fill(pfDeltaEtajj) ;
     if (pfDeltaEtajjWide > -1 && pfDeltaEtajjWide < 1.3) h1_pfMjjWide_nominalHT250->Fill(pfMjjWide) ;
     h1_pfDeltaEtajjWide_nominalHT250->Fill(pfDeltaEtajjWide) ;
     
     h1_recoHT_nominalHT250->Fill(recoHT) ;
     if (recoDeltaEtajj > -1 && recoDeltaEtajj < 1.3) h1_recoMjj_nominalHT250->Fill(recoMjj) ;
     h1_recoDeltaEtajj_nominalHT250->Fill(recoDeltaEtajj) ;
     if (recoDeltaEtajjWide > -1 && recoDeltaEtajjWide < 1.3) h1_recoMjjWide_nominalHT250->Fill(recoMjjWide) ;
     h1_recoDeltaEtajjWide_nominalHT250->Fill(recoDeltaEtajjWide) ;

     h1_recoAK8HT_nominalHT250->Fill(recoAK8HT) ;
     if(recoAK8M > 0) h1_recoAK8M_nominalHT250->Fill(recoAK8M) ;
     if(recoAK8Msd > 0) h1_recoAK8Msd_nominalHT250->Fill(recoAK8Msd) ;
     if(recoAK8Pt > 0) {
       h1_recoAK8Pt_nominalHT250->Fill(recoAK8Pt) ;
       h2_recoAK8MsdPt_nominalHT250->Fill(recoAK8Msd,recoAK8Pt) ;
     }
     if(recoAK8N2sdb1 > 0) h1_recoAK8N2sdb1_nominalHT250->Fill(recoAK8N2sdb1) ;
   }
   std::cout << "HT410no"<< std::endl;

   if (passNominalHT410Trig) {
     h1_pfAK8HT_nominalHT410->Fill(pfAK8HT) ;
     if(pfAK8M > 0) h1_pfAK8M_nominalHT410->Fill(pfAK8M) ;
     if(pfAK8Msd > 0) h1_pfAK8Msd_nominalHT410->Fill(pfAK8Msd) ;
     if(pfAK8Pt > 0) {
       h1_pfAK8Pt_nominalHT410->Fill(pfAK8Pt) ;
       h2_pfAK8MsdPt_nominalHT410->Fill(pfAK8Msd,pfAK8Pt) ;
     }
     if(pfAK8N2sdb1 > 0) h1_pfAK8N2sdb1_nominalHT410->Fill(pfAK8N2sdb1) ;

     h1_caloHT_nominalHT410->Fill(caloHT) ;
     if (caloDeltaEtajj > -1 && caloDeltaEtajj < 1.3) h1_caloMjj_nominalHT410->Fill(caloMjj) ;
     h1_caloDeltaEtajj_nominalHT410->Fill(caloDeltaEtajj) ;
     if (caloDeltaEtajjWide > -1 && caloDeltaEtajjWide < 1.3) h1_caloMjjWide_nominalHT410->Fill(caloMjjWide) ;
     h1_caloDeltaEtajjWide_nominalHT410->Fill(caloDeltaEtajjWide) ;
     
     h1_pfHT_nominalHT410->Fill(pfHT) ;
     if (pfDeltaEtajj > -1 && pfDeltaEtajj < 1.3) h1_pfMjj_nominalHT410->Fill(pfMjj) ;
     h1_pfDeltaEtajj_nominalHT410->Fill(pfDeltaEtajj) ;
     if (pfDeltaEtajjWide > -1 && pfDeltaEtajjWide < 1.3) h1_pfMjjWide_nominalHT410->Fill(pfMjjWide) ;
     h1_pfDeltaEtajjWide_nominalHT410->Fill(pfDeltaEtajjWide) ;
     
     h1_recoHT_nominalHT410->Fill(recoHT) ;
     if (recoDeltaEtajj > -1 && recoDeltaEtajj < 1.3) h1_recoMjj_nominalHT410->Fill(recoMjj) ;
     h1_recoDeltaEtajj_nominalHT410->Fill(recoDeltaEtajj) ;
     if (recoDeltaEtajjWide > -1 && recoDeltaEtajjWide < 1.3) h1_recoMjjWide_nominalHT410->Fill(recoMjjWide) ;
     h1_recoDeltaEtajjWide_nominalHT410->Fill(recoDeltaEtajjWide) ;

     h1_recoAK8HT_nominalHT410->Fill(recoAK8HT) ;
     if(recoAK8M > 0) h1_recoAK8M_nominalHT410->Fill(recoAK8M) ;
     if(recoAK8Msd > 0) h1_recoAK8Msd_nominalHT410->Fill(recoAK8Msd) ;
     if(recoAK8Pt > 0) {
       h1_recoAK8Pt_nominalHT410->Fill(recoAK8Pt) ;
       h2_recoAK8MsdPt_nominalHT410->Fill(recoAK8Msd,recoAK8Pt) ;
     }
     if(recoAK8N2sdb1 > 0) h1_recoAK8N2sdb1_nominalHT410->Fill(recoAK8N2sdb1) ;
   }
   std::cout << "HTno"<< std::endl;

   if (passMonitoringTrig) {
     h1_pfAK8HT_monitoring->Fill(pfAK8HT) ;
     if(pfAK8M > 0) h1_pfAK8M_monitoring->Fill(pfAK8M) ;
     if(pfAK8Msd > 0) {
       h1_pfAK8Msd_monitoring->Fill(pfAK8Msd) ;
     }
     if(pfAK8Pt > 0) {
       h1_pfAK8Pt_monitoring->Fill(pfAK8Pt) ;
       h2_pfAK8MsdPt_monitoring->Fill(pfAK8Msd,pfAK8Pt) ;
     }
     if(pfAK8N2sdb1 > 0) h1_pfAK8N2sdb1_monitoring->Fill(pfAK8N2sdb1) ;
     std::cout << "pfAK8HTmM,Msd,Pt,N2"<< std::endl;

     if(passpfTT == 1 && pfAK8N2sdb1>0 && pfAK8Msd >0 && pfAK8Pt>0){
       h1_pfAK8N2sdb1_monitoring_passTT->Fill(pfAK8N2sdb1);
       h1_pfAK8Msd_monitoring_passTT->Fill(pfAK8Msd);
       h1_pfAK8M_monitoring_passTT->Fill(pfAK8M);
       h1_pfAK8Pt_monitoring_passTT->Fill(pfAK8Pt);
     }
     std::cout << "passpfTT" << std::endl;

     if(passpfTTN2p == 1 && pfAK8N2sdb1>0 && pfAK8Msd >0 && pfAK8Pt>0){
       h1_pfAK8N2sdb1_monitoring_passTTN2p->Fill(pfAK8N2sdb1);
       h1_pfAK8Msd_monitoring_passTTN2p->Fill(pfAK8Msd);
       h1_pfAK8M_monitoring_passTTN2p->Fill(pfAK8M);
       h1_pfAK8Pt_monitoring_passTTN2p->Fill(pfAK8Pt);
     }
     std::cout << "passpfTTp" << std::endl;

     if(passpfTTN2f == 1){
       h1_pfAK8N2sdb1_monitoring_passTTN2f->Fill(pfAK8N2sdb1);
       h1_pfAK8Msd_monitoring_passTTN2f->Fill(pfAK8Msd);
       h1_pfAK8M_monitoring_passTTN2f->Fill(pfAK8M);
       h1_pfAK8Pt_monitoring_passTTN2f->Fill(pfAK8Pt);
     }   
     std::cout << "passpfTTf" << std::endl;

     h1_caloHT_monitoring->Fill(caloHT) ;
     if (caloDeltaEtajj > -1 && caloDeltaEtajj < 1.3) h1_caloMjj_monitoring->Fill(caloMjj) ;
     h1_caloDeltaEtajj_monitoring->Fill(caloDeltaEtajj) ;
     if (caloDeltaEtajjWide > -1 && caloDeltaEtajjWide < 1.3) h1_caloMjjWide_monitoring->Fill(caloMjjWide) ;
     h1_caloDeltaEtajjWide_monitoring->Fill(caloDeltaEtajjWide) ;
     std::cout << "calo" << std::endl;

     h1_pfHT_monitoring->Fill(pfHT) ;
     if (pfDeltaEtajj > -1 && pfDeltaEtajj < 1.3) h1_pfMjj_monitoring->Fill(pfMjj) ;
     h1_pfDeltaEtajj_monitoring->Fill(pfDeltaEtajj) ;
     if (pfDeltaEtajjWide > -1 && pfDeltaEtajjWide < 1.3) h1_pfMjjWide_monitoring->Fill(pfMjjWide) ;
     h1_pfDeltaEtajjWide_monitoring->Fill(pfDeltaEtajjWide) ;
     
     std::cout << "pfHT" << std::endl;

     h1_recoHT_monitoring->Fill(recoHT) ;
     if (recoDeltaEtajj > -1 && recoDeltaEtajj < 1.3) h1_recoMjj_monitoring->Fill(recoMjj) ;
     h1_recoDeltaEtajj_monitoring->Fill(recoDeltaEtajj) ;
     if (recoDeltaEtajjWide > -1 && recoDeltaEtajjWide < 1.3) h1_recoMjjWide_monitoring->Fill(recoMjjWide) ;
     h1_recoDeltaEtajjWide_monitoring->Fill(recoDeltaEtajjWide) ;

     std::cout << "recoWide" << std::endl;


     h1_recoAK8HT_monitoring->Fill(recoAK8HT) ;
     if(recoAK8M > 0) h1_recoAK8M_monitoring->Fill(recoAK8M) ;
     if(recoAK8Msd > 0) h1_recoAK8Msd_monitoring->Fill(recoAK8Msd) ;
     if(recoAK8Pt > 0) {
       h1_recoAK8Pt_monitoring->Fill(recoAK8Pt) ;
       h2_recoAK8MsdPt_monitoring->Fill(recoAK8Msd,recoAK8Pt) ;
     }
     if(recoAK8N2sdb1 > 0) h1_recoAK8N2sdb1_monitoring->Fill(recoAK8N2sdb1) ;
     std::cout << "recoAK8" << std::endl;

     if(passrecoTT == 1 &&recoAK8Msd>0 &&recoAK8M>0 &&recoAK8Pt>0 && recoAK8N2sdb1>0){
       h1_recoAK8N2sdb1_monitoring_passTT->Fill(recoAK8N2sdb1);
       h1_recoAK8Msd_monitoring_passTT->Fill(recoAK8Msd);
       h1_recoAK8M_monitoring_passTT->Fill(recoAK8M);
       h1_recoAK8Pt_monitoring_passTT->Fill(recoAK8Pt);
     }
     std::cout << "recoAK8 passrecoTT" << std::endl;

     if(passrecoTTN2p == 1 &&recoAK8Msd>0 &&recoAK8M>0 &&recoAK8Pt>0 && recoAK8N2sdb1>0){
       h1_recoAK8N2sdb1_monitoring_passTTN2p->Fill(recoAK8N2sdb1);
       h1_recoAK8Msd_monitoring_passTTN2p->Fill(recoAK8Msd);
       h1_recoAK8M_monitoring_passTTN2p->Fill(recoAK8M);
       h1_recoAK8Pt_monitoring_passTTN2p->Fill(recoAK8Pt);
     }
     std::cout << "recoAK8 passrecoTTN2p" << std::endl;

     if(passrecoTTN2f == 1 && recoAK8Msd>0 &&recoAK8M>0 &&recoAK8Pt>0 && recoAK8N2sdb1>0){
       h1_recoAK8N2sdb1_monitoring_passTTN2f->Fill(recoAK8N2sdb1);
       h1_recoAK8Msd_monitoring_passTTN2f->Fill(recoAK8Msd);
       h1_recoAK8M_monitoring_passTTN2f->Fill(recoAK8M);
       h1_recoAK8Pt_monitoring_passTTN2f->Fill(recoAK8Pt);
     }
     std::cout << "recoAK8 passrecoTTN2f" << std::endl;

     if(passrecoTTN2pl == 1 && recoAK8Msd>0 &&recoAK8M>0 &&recoAK8Pt>0 && recoAK8N2sdb1>0){
       h1_recoAK8N2sdb1_monitoring_passTTN2pl->Fill(recoAK8N2sdb1);
       h1_recoAK8Msd_monitoring_passTTN2pl->Fill(recoAK8Msd);
       h1_recoAK8M_monitoring_passTTN2pl->Fill(recoAK8M);
       h1_recoAK8Pt_monitoring_passTTN2pl->Fill(recoAK8Pt);
     }
     std::cout << "recoAK8 passrecoTTN2pl" << std::endl;
     
     /*
     if(passrecoTTN2fl == 1 && recoAK8Msd>0 &&recoAK8M>0 &&recoAK8Pt>0 && recoAK8N2sdb1>0){
       h1_recoAK8N2sdb1_monitoring_passTTN2fl->Fill(recoAK8N2sdb1);
       h1_recoAK8Msd_monitoring_passTTN2fl->Fill(recoAK8Msd);
       h1_recoAK8M_monitoring_passTTN2fl->Fill(recoAK8M);
       h1_recoAK8Pt_monitoring_passTTN2fl->Fill(recoAK8Pt);
     }
     std::cout << "recoAK8 passrecoTTN2fl" << std::endl;
     */

   }
   std::cout << "end hist" << std::endl;

   outTree_->Fill();
   
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
HTScoutingAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HTScoutingAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HTScoutingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HTScoutingAnalyzer);
