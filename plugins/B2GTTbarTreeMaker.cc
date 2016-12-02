// -*- C++ -*-
//
// Package:    Analysis/B2GTTbarTreeMaker
// Class:      B2GTTbarTreeMaker
//
/**\class B2GTTbarTreeMaker B2GTTbarTreeMaker.cc Analysis/B2GTTbarTreeMaker/plugins/B2GTTbarTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  James Dolen
//         Created:  Sat, 30 Apr 2016 17:40:42 GMT
//
//

//--------------------------
// To add:
// - NNPDF3weight
// - electron quality cuts
// - trigger
//
// Note. The following items should be applied at the tree reader level:
// - MU ID (HIP)
// -- TFile* f_muID = TFile::Open("MuonID_Z_RunBCD_prompt80X_7p65.root","read");
// -- TH1F* h_muID = (TH1F*) f_muID->Get("MC_NUM_MediumID_DEN_genTracks_PAR_eta/eta_ratio")->Clone();
// -- float SF_muID = h_muID->GetBinContent(h_muID->FindBin(eta););
// - BTagCalibrationReader
//--------------------------


// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <bitset>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// DataFormats
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Gen particle
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// JER
#include "JetMETCorrections/Modules/interface/JetResolution.h"

// Electron
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// LHE weights
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Utilities
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// root
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"

//RS gluon PDF weights
namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  int numberPDF(int nset);
  void usePDFMember(int nset, int member);
  double xfx(int nset, double x, double Q, int fl);
  double getXmin(int nset, int member);
  double getXmax(int nset, int member);
  double getQ2min(int nset, int member);
  double getQ2max(int nset, int member);
  void extrapolate(bool extrapolate=true);
}

//
// class declaration
//

class B2GTTbarTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit B2GTTbarTreeMaker(const edm::ParameterSet&);
      ~B2GTTbarTreeMaker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::JetCollection> ak4jetToken_;
      edm::EDGetTokenT<pat::JetCollection> ak8jetToken_;
      edm::EDGetTokenT<pat::JetCollection> puppijetToken_;
      edm::EDGetTokenT<pat::JetCollection> ak8CHSSoftDropSubjetsToken_;
      edm::EDGetTokenT<pat::JetCollection> ak8PuppiSoftDropSubjetsToken_;
      edm::EDGetTokenT<reco::GenJetCollection> ak4genjetToken_;
      edm::EDGetTokenT<reco::GenJetCollection> ak8genjetToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsMETFilterToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetTokenT<bool> badMuonFilterToken_;
      edm::EDGetTokenT<bool> badChargedCandidateFilterToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken_;
      edm::EDGetTokenT<LHEEventProduct> theSrc_;
      edm::EDGetTokenT<GenEventInfoProduct> pdfToken_;

      bool useToolbox_;
      bool verbose_;
      bool verboseGen_;
      bool runGenLoop_;
      bool isZprime_;
      bool isttbar_;
      bool isRSG_;
      std::vector<std::string>  jecPayloadsAK4chs_;
      std::vector<std::string>  jecPayloadsAK8chs_;
      std::vector<std::string>  jecPayloadsAK4pup_;
      std::vector<std::string>  jecPayloadsAK8pup_;
      boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK4chs;
      boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK8chs;
      boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK4pup;
      boost::shared_ptr<FactorizedJetCorrector> JetCorrectorAK8pup;
      boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK4chs;
      boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK8chs;
      boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK4pup;
      boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK8pup;

      std::string jerSFtext_;

      TFile* fPUweight;
      TH1D* hPUweight;
      TH1D* hPUweight_MBup;
      TH1D* hPUweight_MBdn;

      int count_GenTruth_Leptonic ;
      int count_nMu_gt1 ;
      int count_nEl_gt1 ;
      int count_nMu_e1 ;
      int count_nEl_e1 ;
      int count_nLep_e1 ;
      int count_JetPt300 ;
      int count_JetPt300Eta ;
      int count_JetPt300Eta_AK4 ;
      int count_JetPt300Eta_muPt40 ;
      int count_JetPt300Eta_muPt40_MET40 ;
      int count_JetPt300Eta_muPt40_MET40_AK4 ;


      TH1D * h_ak8puppi_softDropMass ;
      TH1D * h_ak8chs_softDropMass   ;
      TH1D * h_ak8chs_softDropMass_reweighted ;
      TH1D * h_ak8chs_pt   ;
      TH1D * h_ak8chs_pt_reweighted ;
      TH1D * h_NtrueIntPU ;
      TH1D * h_NPV           ;
      TH1D * h_NPVreweighted ;




/*

      Tree with double leptons

*/

      TTree *TreeLept;
      std::vector<std::string> *LeptTrigNames     = new std::vector<std::string>;
      std::vector<int> *LeptTrigPrescales = new std::vector<int>;
      std::vector<bool> *LeptTrigPass    = new std::vector<bool>;



      std::string LeptTrigAcceptBits;

      Float_t JetPtRaw                               ;
      Float_t JetEtaRaw                              ;
      Float_t JetPhiRaw                              ;
      Float_t JetMassRaw                             ;
      Float_t JetP                                   ;
      Float_t JetPt                                  ;
      Float_t JetEta                                 ;
      Float_t JetPhi                                 ;
      Float_t JetRap                                 ;
      Float_t JetEnergy                              ;
      Float_t JetMass                                ;
      Float_t JetArea                                ;
      Float_t JetSDmass                              ;
      Float_t JetSDmassRaw                           ;
      Float_t JetSDmassCorrL23                       ;
      Float_t JetSDmassCorrL23Up                     ;
      Float_t JetSDmassCorrL23Dn                     ;
      Float_t JetSDmassCorrL123                      ;
      Float_t JetSDmassCorrL123Up                    ;
      Float_t JetSDmassCorrL123Dn                    ;
      Float_t JetSDmassCorrL23Smear                  ;
      Float_t JetSDmassCorrL23SmearUp                ;
      Float_t JetSDmassCorrL23SmearDn                ;
      Float_t JetSDptRaw                             ;
      Float_t JetSDptCorrL23                         ;
      Float_t JetSDptCorrL23Up                       ;
      Float_t JetSDptCorrL23Dn                       ;
      Float_t JetSDptCorrL123                        ;
      Float_t JetSDptCorrL123Up                      ;
      Float_t JetSDptCorrL123Dn                      ;
      Float_t JetSDptCorrL23Smear                    ;
      Float_t JetSDptCorrL23SmearUp                  ;
      Float_t JetSDptCorrL23SmearDn                  ;
      Float_t JetSDetaRaw                            ;
      Float_t JetSDphiRaw                            ;
      Float_t JetMassPruned                          ;
      Float_t JetMassTrimmed                         ;
      Float_t JetTau1                                ;
      Float_t JetTau2                                ;
      Float_t JetTau3                                ;
      Float_t JetTau4                                ;
      Float_t JetTau32                               ;
      Float_t JetTau21                               ;
      Float_t JetSDsubjet0bdisc                      ;
      Float_t JetSDsubjet1bdisc                      ;
      Float_t JetSDmaxbdisc                          ;
      Float_t JetSDmaxbdiscflavHadron                ;
      Float_t JetSDmaxbdiscflavParton                ;
      Float_t JetSDsubjet0pt                         ;
      Float_t JetSDsubjet0mass                       ;
      Float_t JetSDsubjet0eta                        ;
      Float_t JetSDsubjet0phi                        ;
      Float_t JetSDsubjet0area                       ;
      Float_t JetSDsubjet0flavHadron                 ;
      Float_t JetSDsubjet0flavParton                 ;
      Float_t JetSDsubjet0tau1                       ;
      Float_t JetSDsubjet0tau2                       ;
      Float_t JetSDsubjet0tau3                       ;
      Float_t JetSDsubjet1pt                         ;
      Float_t JetSDsubjet1mass                       ;
      Float_t JetSDsubjet1eta                        ;
      Float_t JetSDsubjet1phi                        ;
      Float_t JetSDsubjet1area                       ;
      Float_t JetSDsubjet1flavHadron                 ;
      Float_t JetSDsubjet1flavParton                 ;
      Float_t JetSDsubjet1tau1                       ;
      Float_t JetSDsubjet1tau2                       ;
      Float_t JetSDsubjet1tau3                       ;
      Float_t JetPuppiP                              ;
      Float_t JetPuppiPt                             ;
      Float_t JetPuppiEta                            ;
      Float_t JetPuppiPhi                            ;
      Float_t JetPuppiMass                           ;
      Float_t JetPuppiSDmass                         ;
      Float_t JetPuppiSDmassCorr                     ;
      Float_t JetPuppiSDmassCorrUp                   ;
      Float_t JetPuppiSDmassCorrDn                   ;
      Float_t JetPuppiSDmassCorrL23Smear             ;
      Float_t JetPuppiSDmassCorrL23SmearUp           ;
      Float_t JetPuppiSDmassCorrL23SmearDn           ;
      Float_t JetPuppiSDpt                           ;
      Float_t JetPuppiSDptCorr                       ;
      Float_t JetPuppiSDptCorrUp                     ;
      Float_t JetPuppiSDptCorrDn                     ;
      Float_t JetPuppiSDptCorrL23Smear               ;
      Float_t JetPuppiSDptCorrL23SmearUp             ;
      Float_t JetPuppiSDptCorrL23SmearDn             ;
      Float_t JetPuppiSDeta                          ;
      Float_t JetPuppiSDphi                          ;
      Float_t JetPuppiTau1                           ;
      Float_t JetPuppiTau2                           ;
      Float_t JetPuppiTau3                           ;
      Float_t JetPuppiTau4                           ;
      Float_t JetPuppiTau32                          ;
      Float_t JetPuppiTau21                          ;
      Float_t JetPuppiSDsubjet0bdisc                 ;
      Float_t JetPuppiSDsubjet1bdisc                 ;
      Float_t JetPuppiSDmaxbdisc                     ;
      Float_t JetPuppiSDmaxbdiscflavHadron           ;
      Float_t JetPuppiSDmaxbdiscflavParton           ;
      Float_t JetPuppiSDsubjet0pt                    ;
      Float_t JetPuppiSDsubjet0mass                  ;
      Float_t JetPuppiSDsubjet0eta                   ;
      Float_t JetPuppiSDsubjet0phi                   ;
      Float_t JetPuppiSDsubjet0area                  ;
      Float_t JetPuppiSDsubjet0flavHadron            ;
      Float_t JetPuppiSDsubjet0flavParton            ;
      Float_t JetPuppiSDsubjet0tau1                  ;
      Float_t JetPuppiSDsubjet0tau2                  ;
      Float_t JetPuppiSDsubjet0tau3                  ;
      Float_t JetPuppiSDsubjet1pt                    ;
      Float_t JetPuppiSDsubjet1mass                  ;
      Float_t JetPuppiSDsubjet1eta                   ;
      Float_t JetPuppiSDsubjet1phi                   ;
      Float_t JetPuppiSDsubjet1area                  ;
      Float_t JetPuppiSDsubjet1flavHadron            ;
      Float_t JetPuppiSDsubjet1flavParton            ;
      Float_t JetPuppiSDsubjet1tau1                  ;
      Float_t JetPuppiSDsubjet1tau2                  ;
      Float_t JetPuppiSDsubjet1tau3                  ;
      Float_t JetCHF                                 ;
      Float_t JetNHF                                 ;
      Float_t JetCM                                  ;
      Float_t JetNM                                  ;
      Float_t JetNEF                                 ;
      Float_t JetCEF                                 ;
      Float_t JetMF                                  ;
      Float_t JetMult                                ;
      Float_t JetPuppiCHF                            ;
      Float_t JetPuppiNHF                            ;
      Float_t JetPuppiCM                             ;
      Float_t JetPuppiNM                             ;
      Float_t JetPuppiNEF                            ;
      Float_t JetPuppiCEF                            ;
      Float_t JetPuppiMF                             ;
      Float_t JetPuppiMult                           ;
      Float_t JetMassCorrFactor                      ;
      Float_t JetMassCorrFactorUp                    ;
      Float_t JetMassCorrFactorDn                    ;
      Float_t JetCorrFactor                          ;
      Float_t JetCorrFactorUp                        ;
      Float_t JetCorrFactorDn                        ;
      Float_t JetPtSmearFactor                       ;
      Float_t JetPtSmearFactorUp                     ;
      Float_t JetPtSmearFactorDn                     ;
      Float_t JetPuppiMassCorrFactor                 ;
      Float_t JetPuppiMassCorrFactorUp               ;
      Float_t JetPuppiMassCorrFactorDn               ;
      Float_t JetPuppiCorrFactor                     ;
      Float_t JetPuppiCorrFactorUp                   ;
      Float_t JetPuppiCorrFactorDn                   ;
      Float_t JetPuppiPtSmearFactor                  ;
      Float_t JetPuppiPtSmearFactorUp                ;
      Float_t JetPuppiPtSmearFactorDn                ;
      Float_t JetEtaScaleFactor                      ;
      Float_t JetPhiScaleFactor                      ;
      // Float_t JetMatchedGenJetDR                     ;
      Float_t JetMatchedGenJetPt                     ;
      Float_t JetMatchedGenJetMass                   ;
      Int_t   JetGenMatched_TopHadronic              ;
      Float_t JetGenMatched_TopPt                    ;
      Float_t JetGenMatched_TopEta                   ;
      Float_t JetGenMatched_TopPhi                   ;
      Float_t JetGenMatched_TopMass                  ;
      Float_t JetGenMatched_bPt                      ;
      Float_t JetGenMatched_WPt                      ;
      Float_t JetGenMatched_Wd1Pt                    ;
      Float_t JetGenMatched_Wd2Pt                    ;
      Float_t JetGenMatched_Wd1ID                    ;
      Float_t JetGenMatched_Wd2ID                    ;
      Float_t JetGenMatched_MaxDeltaRPartonTop       ;
      Float_t JetGenMatched_MaxDeltaRWPartonTop      ;
      Float_t JetGenMatched_MaxDeltaRWPartonW        ;
      Float_t JetGenMatched_DeltaR_t_b               ;
      Float_t JetGenMatched_DeltaR_t_W               ;
      Float_t JetGenMatched_DeltaR_t_Wd1             ;
      Float_t JetGenMatched_DeltaR_t_Wd2             ;
      Float_t JetGenMatched_DeltaR_W_b1              ;
      Float_t JetGenMatched_DeltaR_W_Wd1             ;
      Float_t JetGenMatched_DeltaR_W_Wd2             ;
      Float_t JetGenMatched_DeltaR_Wd1_Wd2           ;
      Float_t JetGenMatched_DeltaR_Wd1_b             ;
      Float_t JetGenMatched_DeltaR_Wd2_b             ;
      Float_t JetGenMatched_DeltaR_jet_t             ;
      Float_t JetGenMatched_DeltaR_jet_W             ;
      Float_t JetGenMatched_DeltaR_jet_b             ;
      Float_t JetGenMatched_DeltaR_jet_Wd1           ;
      Float_t JetGenMatched_DeltaR_jet_Wd2           ;
      Float_t JetGenMatched_DeltaR_pup0_b            ;
      Float_t JetGenMatched_DeltaR_pup0_Wd1          ;
      Float_t JetGenMatched_DeltaR_pup0_Wd2          ;
      Float_t JetGenMatched_DeltaR_pup1_b            ;
      Float_t JetGenMatched_DeltaR_pup1_Wd1          ;
      Float_t JetGenMatched_DeltaR_pup1_Wd2          ;
      Float_t JetGenMatched_partonPt                 ;
      Float_t JetGenMatched_partonEta                ;
      Float_t JetGenMatched_partonPhi                ;
      Float_t JetGenMatched_partonMass               ;
      Float_t JetGenMatched_partonID                 ;
      Float_t JetGenMatched_DeltaRjetParton          ;
      Float_t LeptMETpx                          ;
      Float_t LeptMETpy                          ;
      Float_t LeptMETpt                          ;
      Float_t LeptMETphi                         ;
      Float_t LeptMETsumET                       ;
      Float_t LeptNvtx                           ;
      Float_t LeptNPUtrue                        ;
      Float_t LeptRho                            ;
      Float_t LeptEventWeight                    ;
      Float_t LeptPUweight       ;
      Float_t LeptPUweight_MBup  ;
      Float_t LeptPUweight_MBdn  ;


      Float_t LeptGenTTmass                      ;

      Float_t HTlep0                                  ;
      Float_t ST0                                     ;
      Float_t ST0_CorrDn                              ;
      Float_t ST0_CorrUp                              ;
      Float_t ST0_PtSmearNom                          ;
      Float_t ST0_PtSmearUp                           ;
      Float_t ST0_PtSmearDn                           ;

      Float_t HTlep1                                  ;
      Float_t ST1                                     ;
      Float_t ST1_CorrDn                              ;
      Float_t ST1_CorrUp                              ;
      Float_t ST1_PtSmearNom                          ;
      Float_t ST1_PtSmearUp                           ;
      Float_t ST1_PtSmearDn                           ;


      Float_t LeptRunNum                         ;
      Float_t LeptLumiBlock                      ;
      Float_t LeptEventNum                       ;
      Int_t   LeptPassMETFilters                 ;

      Float_t AK4dRminPt                             ;
      Float_t AK4dRminEta                            ;
      Float_t AK4dRminPhi                            ;
      Float_t AK4dRminMass                           ;
      Float_t AK4dRminBdisc                          ;
      Float_t AK4dRminLep                            ;
      Float_t AK4BtagdRminPt                         ;
      Float_t AK4BtagdRminBdisc                      ;
      Float_t AK4BtagdRminLep                        ;
      Int_t   LepHemiContainsAK4BtagLoose            ;
      Int_t   LepHemiContainsAK4BtagMedium           ;
      Int_t   LepHemiContainsAK4BtagTight            ;


      Float_t Lepton0Phi                              ;
      Float_t Lepton0Pt                               ;
      Float_t Lepton0Eta                              ;
      Float_t Lepton0Mass                             ;
      Float_t PtRel                                  ;
      Int_t   LeptonIsMu                             ;
      Int_t   Mu0HighPt                                ;
      Int_t   Mu0Tight                                ;
      Int_t   Mu0Medium                               ;
      Float_t DeltaRJetLep0                           ;
      Float_t DeltaPhiJetLep0                         ;
      Float_t Mu0Iso                                  ;
      Float_t Mu1Iso                                  ;

      Float_t Elecron0_absiso                         ;
      Float_t Elecron0_relIsoWithDBeta                ;
      Float_t Elecron0_absiso_EA                      ;
      Float_t Elecron0_relIsoWithEA                   ;

      Int_t Electron0_isMedium   ;
      Int_t Electron0_isTight    ;

      Float_t Lepton1Phi                              ;
      Float_t Lepton1Pt                               ;
      Float_t Lepton1Eta                              ;
      Float_t Lepton1Mass                             ;
      Int_t   Mu1HighPt                                ;
      Int_t   Mu1Tight                                ;
      Int_t   Mu1Medium                               ;
      Float_t DeltaRJetLep1                           ;
      Float_t DeltaPhiJetLep1                         ;
      Float_t Elecron1_absiso                         ;
      Float_t Elecron1_relIsoWithDBeta                ;
      Float_t Elecron1_absiso_EA                      ;
      Float_t Elecron1_relIsoWithEA                   ;

      Int_t Electron1_isMedium   ;
      Int_t Electron1_isTight    ;


      /*

            Single Gamma



            TTree *TreeSingleGamma;
            std::vector<std::string> *SingleGammaTrigNames     = new std::vector<std::string>;
            std::vector<int> *SingleGammaTrigPrescales = new std::vector<int>;
            std::vector<bool> *SingleGammaTrigPass    = new std::vector<bool>;



            std::string SingleGammaTrigAcceptBits;

            Float_t JetPtRaw                               ;
            Float_t JetEtaRaw                              ;
            Float_t JetPhiRaw                              ;
            Float_t JetMassRaw                             ;
            Float_t JetP                                   ;
            Float_t JetPt                                  ;
            Float_t JetEta                                 ;
            Float_t JetPhi                                 ;
            Float_t JetRap                                 ;
            Float_t JetEnergy                              ;
            Float_t JetMass                                ;
            Float_t JetArea                                ;
            Float_t JetSDmass                              ;
            Float_t JetSDmassRaw                           ;
            Float_t JetSDmassCorrL23                       ;
            Float_t JetSDmassCorrL23Up                     ;
            Float_t JetSDmassCorrL23Dn                     ;
            Float_t JetSDmassCorrL123                      ;
            Float_t JetSDmassCorrL123Up                    ;
            Float_t JetSDmassCorrL123Dn                    ;
            Float_t JetSDmassCorrL23Smear                  ;
            Float_t JetSDmassCorrL23SmearUp                ;
            Float_t JetSDmassCorrL23SmearDn                ;
            Float_t JetSDptRaw                             ;
            Float_t JetSDptCorrL23                         ;
            Float_t JetSDptCorrL23Up                       ;
            Float_t JetSDptCorrL23Dn                       ;
            Float_t JetSDptCorrL123                        ;
            Float_t JetSDptCorrL123Up                      ;
            Float_t JetSDptCorrL123Dn                      ;
            Float_t JetSDptCorrL23Smear                    ;
            Float_t JetSDptCorrL23SmearUp                  ;
            Float_t JetSDptCorrL23SmearDn                  ;
            Float_t JetSDetaRaw                            ;
            Float_t JetSDphiRaw                            ;
            Float_t JetMassPruned                          ;
            Float_t JetMassTrimmed                         ;
            Float_t JetTau1                                ;
            Float_t JetTau2                                ;
            Float_t JetTau3                                ;
            Float_t JetTau4                                ;
            Float_t JetTau32                               ;
            Float_t JetTau21                               ;
            Float_t JetSDsubjet0bdisc                      ;
            Float_t JetSDsubjet1bdisc                      ;
            Float_t JetSDmaxbdisc                          ;
            Float_t JetSDmaxbdiscflavHadron                ;
            Float_t JetSDmaxbdiscflavParton                ;
            Float_t JetSDsubjet0pt                         ;
            Float_t JetSDsubjet0mass                       ;
            Float_t JetSDsubjet0eta                        ;
            Float_t JetSDsubjet0phi                        ;
            Float_t JetSDsubjet0area                       ;
            Float_t JetSDsubjet0flavHadron                 ;
            Float_t JetSDsubjet0flavParton                 ;
            Float_t JetSDsubjet0tau1                       ;
            Float_t JetSDsubjet0tau2                       ;
            Float_t JetSDsubjet0tau3                       ;
            Float_t JetSDsubjet1pt                         ;
            Float_t JetSDsubjet1mass                       ;
            Float_t JetSDsubjet1eta                        ;
            Float_t JetSDsubjet1phi                        ;
            Float_t JetSDsubjet1area                       ;
            Float_t JetSDsubjet1flavHadron                 ;
            Float_t JetSDsubjet1flavParton                 ;
            Float_t JetSDsubjet1tau1                       ;
            Float_t JetSDsubjet1tau2                       ;
            Float_t JetSDsubjet1tau3                       ;
            Float_t JetPuppiP                              ;
            Float_t JetPuppiPt                             ;
            Float_t JetPuppiEta                            ;
            Float_t JetPuppiPhi                            ;
            Float_t JetPuppiMass                           ;
            Float_t JetPuppiSDmass                         ;
            Float_t JetPuppiSDmassCorr                     ;
            Float_t JetPuppiSDmassCorrUp                   ;
            Float_t JetPuppiSDmassCorrDn                   ;
            Float_t JetPuppiSDmassCorrL23Smear             ;
            Float_t JetPuppiSDmassCorrL23SmearUp           ;
            Float_t JetPuppiSDmassCorrL23SmearDn           ;
            Float_t JetPuppiSDpt                           ;
            Float_t JetPuppiSDptCorr                       ;
            Float_t JetPuppiSDptCorrUp                     ;
            Float_t JetPuppiSDptCorrDn                     ;
            Float_t JetPuppiSDptCorrL23Smear               ;
            Float_t JetPuppiSDptCorrL23SmearUp             ;
            Float_t JetPuppiSDptCorrL23SmearDn             ;
            Float_t JetPuppiSDeta                          ;
            Float_t JetPuppiSDphi                          ;
            Float_t JetPuppiTau1                           ;
            Float_t JetPuppiTau2                           ;
            Float_t JetPuppiTau3                           ;
            Float_t JetPuppiTau4                           ;
            Float_t JetPuppiTau32                          ;
            Float_t JetPuppiTau21                          ;
            Float_t JetPuppiSDsubjet0bdisc                 ;
            Float_t JetPuppiSDsubjet1bdisc                 ;
            Float_t JetPuppiSDmaxbdisc                     ;
            Float_t JetPuppiSDmaxbdiscflavHadron           ;
            Float_t JetPuppiSDmaxbdiscflavParton           ;
            Float_t JetPuppiSDsubjet0pt                    ;
            Float_t JetPuppiSDsubjet0mass                  ;
            Float_t JetPuppiSDsubjet0eta                   ;
            Float_t JetPuppiSDsubjet0phi                   ;
            Float_t JetPuppiSDsubjet0area                  ;
            Float_t JetPuppiSDsubjet0flavHadron            ;
            Float_t JetPuppiSDsubjet0flavParton            ;
            Float_t JetPuppiSDsubjet0tau1                  ;
            Float_t JetPuppiSDsubjet0tau2                  ;
            Float_t JetPuppiSDsubjet0tau3                  ;
            Float_t JetPuppiSDsubjet1pt                    ;
            Float_t JetPuppiSDsubjet1mass                  ;
            Float_t JetPuppiSDsubjet1eta                   ;
            Float_t JetPuppiSDsubjet1phi                   ;
            Float_t JetPuppiSDsubjet1area                  ;
            Float_t JetPuppiSDsubjet1flavHadron            ;
            Float_t JetPuppiSDsubjet1flavParton            ;
            Float_t JetPuppiSDsubjet1tau1                  ;
            Float_t JetPuppiSDsubjet1tau2                  ;
            Float_t JetPuppiSDsubjet1tau3                  ;
            Float_t JetCHF                                 ;
            Float_t JetNHF                                 ;
            Float_t JetCM                                  ;
            Float_t JetNM                                  ;
            Float_t JetNEF                                 ;
            Float_t JetCEF                                 ;
            Float_t JetMF                                  ;
            Float_t JetMult                                ;
            Float_t JetPuppiCHF                            ;
            Float_t JetPuppiNHF                            ;
            Float_t JetPuppiCM                             ;
            Float_t JetPuppiNM                             ;
            Float_t JetPuppiNEF                            ;
            Float_t JetPuppiCEF                            ;
            Float_t JetPuppiMF                             ;
            Float_t JetPuppiMult                           ;
            Float_t JetMassCorrFactor                      ;
            Float_t JetMassCorrFactorUp                    ;
            Float_t JetMassCorrFactorDn                    ;
            Float_t JetCorrFactor                          ;
            Float_t JetCorrFactorUp                        ;
            Float_t JetCorrFactorDn                        ;
            Float_t JetPtSmearFactor                       ;
            Float_t JetPtSmearFactorUp                     ;
            Float_t JetPtSmearFactorDn                     ;
            Float_t JetPuppiMassCorrFactor                 ;
            Float_t JetPuppiMassCorrFactorUp               ;
            Float_t JetPuppiMassCorrFactorDn               ;
            Float_t JetPuppiCorrFactor                     ;
            Float_t JetPuppiCorrFactorUp                   ;
            Float_t JetPuppiCorrFactorDn                   ;
            Float_t JetPuppiPtSmearFactor                  ;
            Float_t JetPuppiPtSmearFactorUp                ;
            Float_t JetPuppiPtSmearFactorDn                ;
            Float_t JetEtaScaleFactor                      ;
            Float_t JetPhiScaleFactor                      ;
            // Float_t JetMatchedGenJetDR                     ;
            Float_t JetMatchedGenJetPt                     ;
            Float_t JetMatchedGenJetMass                   ;
            Int_t   JetGenMatched_TopHadronic              ;
            Float_t JetGenMatched_TopPt                    ;
            Float_t JetGenMatched_TopEta                   ;
            Float_t JetGenMatched_TopPhi                   ;
            Float_t JetGenMatched_TopMass                  ;
            Float_t JetGenMatched_bPt                      ;
            Float_t JetGenMatched_WPt                      ;
            Float_t JetGenMatched_Wd1Pt                    ;
            Float_t JetGenMatched_Wd2Pt                    ;
            Float_t JetGenMatched_Wd1ID                    ;
            Float_t JetGenMatched_Wd2ID                    ;
            Float_t JetGenMatched_MaxDeltaRPartonTop       ;
            Float_t JetGenMatched_MaxDeltaRWPartonTop      ;
            Float_t JetGenMatched_MaxDeltaRWPartonW        ;
            Float_t JetGenMatched_DeltaR_t_b               ;
            Float_t JetGenMatched_DeltaR_t_W               ;
            Float_t JetGenMatched_DeltaR_t_Wd1             ;
            Float_t JetGenMatched_DeltaR_t_Wd2             ;
            Float_t JetGenMatched_DeltaR_W_b1              ;
            Float_t JetGenMatched_DeltaR_W_Wd1             ;
            Float_t JetGenMatched_DeltaR_W_Wd2             ;
            Float_t JetGenMatched_DeltaR_Wd1_Wd2           ;
            Float_t JetGenMatched_DeltaR_Wd1_b             ;
            Float_t JetGenMatched_DeltaR_Wd2_b             ;
            Float_t JetGenMatched_DeltaR_jet_t             ;
            Float_t JetGenMatched_DeltaR_jet_W             ;
            Float_t JetGenMatched_DeltaR_jet_b             ;
            Float_t JetGenMatched_DeltaR_jet_Wd1           ;
            Float_t JetGenMatched_DeltaR_jet_Wd2           ;
            Float_t JetGenMatched_DeltaR_pup0_b            ;
            Float_t JetGenMatched_DeltaR_pup0_Wd1          ;
            Float_t JetGenMatched_DeltaR_pup0_Wd2          ;
            Float_t JetGenMatched_DeltaR_pup1_b            ;
            Float_t JetGenMatched_DeltaR_pup1_Wd1          ;
            Float_t JetGenMatched_DeltaR_pup1_Wd2          ;
            Float_t JetGenMatched_partonPt                 ;
            Float_t JetGenMatched_partonEta                ;
            Float_t JetGenMatched_partonPhi                ;
            Float_t JetGenMatched_partonMass               ;
            Float_t JetGenMatched_partonID                 ;
            Float_t JetGenMatched_DeltaRjetParton          ;
            Float_t GammaMETpx                          ;
            Float_t GammaMETpy                          ;
            Float_t GammaMETpt                          ;
            Float_t GammaMETphi                         ;
            Float_t GammaMETsumET                       ;
            Float_t GammaNvtx                           ;
            Float_t GammaNPUtrue                        ;
            Float_t GammaRho                            ;
            Float_t GammaEventWeight                    ;
            Float_t GammaPUweight       ;
            Float_t GammaPUweight_MBup  ;
            Float_t GammaPUweight_MBdn  ;


            Float_t ST                                     ;
            Float_t ST_CorrDn                              ;
            Float_t ST_CorrUp                              ;
            Float_t ST_PtSmearNom                          ;
            Float_t ST_PtSmearUp                           ;
            Float_t ST_PtSmearDn                           ;

            Float_t AK4dRminPt                             ;
            Float_t AK4dRminEta                            ;
            Float_t AK4dRminPhi                            ;
            Float_t AK4dRminMass                           ;
            Float_t AK4dRminBdisc                          ;
            Float_t AK4dRminLep                            ;
            Float_t AK4BtagdRminPt                         ;
            Float_t AK4BtagdRminBdisc                      ;
            Float_t AK4BtagdRminLep                        ;
  */

};

//
// constructors and destructor
//
B2GTTbarTreeMaker::B2GTTbarTreeMaker(const edm::ParameterSet& iConfig):
    ak4jetToken_(consumes<pat::JetCollection>(edm::InputTag("slimmedJets"))),
    ak8jetToken_(consumes<pat::JetCollection>(    iConfig.getParameter<edm::InputTag>("ak8chsInput"))),  //edm::InputTag("slimmedJetsAK8"))),
    puppijetToken_(consumes<pat::JetCollection>(  iConfig.getParameter<edm::InputTag>("ak8puppiInput"))),
    ak8CHSSoftDropSubjetsToken_(consumes<pat::JetCollection>(  iConfig.getParameter<edm::InputTag>("ak8chsSubjetsInput"))),
    ak8PuppiSoftDropSubjetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8puppiSubjetsInput"))),
    ak4genjetToken_(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJets"))),
    ak8genjetToken_(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJetsAK8"))),
    prunedGenToken_(consumes<edm::View<reco::GenParticle> >(edm::InputTag("prunedGenParticles"))),
    rhoToken_(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"))),
    vtxToken_(consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"))),
    triggerResultsMETFilterToken_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "RECO"))),  //"PAT"
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),//"TriggerResults", "", "HLT2"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("patTrigger"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),  //   selectedPatTrigger))),
    badMuonFilterToken_(consumes<bool>(edm::InputTag("BadPFMuonFilter",""))),
    badChargedCandidateFilterToken_(consumes<bool>(edm::InputTag("BadChargedCandidateFilter",""))),
    muonToken_(consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"))),
    electronToken_(consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"))),
    metToken_(consumes<pat::METCollection>(edm::InputTag("slimmedMETs"))),
    pileupInfoToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),
    theSrc_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("theSrc"))),
    pdfToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
    useToolbox_(iConfig.getParameter<bool>  ("useToolbox")),
    verbose_(iConfig.getParameter<bool>  ("verbose")),
    verboseGen_(iConfig.getParameter<bool>  ("verboseGen")),
    runGenLoop_(iConfig.getParameter<bool>  ("runGenLoop")),
    jecPayloadsAK4chs_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK4chs")),
    jecPayloadsAK8chs_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK8chs")),
    jecPayloadsAK4pup_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK4pup")),
    jecPayloadsAK8pup_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK8pup")),
    jerSFtext_ (iConfig.getParameter<std::string>  ("jerSFtext"))
{
  std::cout<<"B2GTTbarTreeMaker::B2GTTbarTreeMaker"<<std::endl;

  //RS gluon PDF weights
  LHAPDF::initPDFSet(1, "NNPDF30_lo_as_0130");

  usesResource("TFileService");

  edm::Service<TFileService> fs;

  h_ak8puppi_softDropMass            =  fs->make<TH1D>("h_ak8puppi_softDropMass"           ,"",200,0,400);
  h_ak8chs_softDropMass              =  fs->make<TH1D>("h_ak8chs_softDropMass"             ,"",200,0,400);
  h_ak8chs_softDropMass_reweighted   =  fs->make<TH1D>("h_ak8chs_softDropMass_reweighted"  ,"",200,0,400);
  h_ak8chs_pt                        =  fs->make<TH1D>("h_ak8chs_pt"                       ,"",200,0,4000);
  h_ak8chs_pt_reweighted             =  fs->make<TH1D>("h_ak8chs_pt_reweighted"            ,"",200,0,4000);
  h_NtrueIntPU                       =  fs->make<TH1D>("h_NtrueIntPU"                      ,"",200,0,200);
  h_NPV                              =  fs->make<TH1D>("h_NPV"                             ,"",200,0,200);
  h_NPVreweighted                    =  fs->make<TH1D>("h_NPVreweighted"                   ,"",200,0,200);


/*

    Lepton Tree

*/


  TreeLept = new TTree("TreeLept","TreeLept");

  TreeLept->Branch("LeptTrigNames"    , "vector<std::string>", &LeptTrigNames);
  TreeLept->Branch("LeptTrigPrescales" , "vector<int>",  &LeptTrigPrescales);
  TreeLept->Branch("LeptTrigPass"      , "vector<bool>", &LeptTrigPass);
  TreeLept->Branch("LeptTrigAcceptBits", &LeptTrigAcceptBits);



  TreeLept->Branch("JetPtRaw"                             , & JetPtRaw                          ,    "JetPtRaw/F"                               );
  TreeLept->Branch("JetEtaRaw"                            , & JetEtaRaw                         ,    "JetEtaRaw/F"                              );
  TreeLept->Branch("JetPhiRaw"                            , & JetPhiRaw                         ,    "JetPhiRaw/F"                              );
  TreeLept->Branch("JetMassRaw"                           , & JetMassRaw                        ,    "JetMassRaw/F"                             );
  TreeLept->Branch("JetP"                                 , & JetP                              ,    "JetP/F"                                   );
  TreeLept->Branch("JetPt"                                , & JetPt                             ,    "JetPt/F"                                  );
  TreeLept->Branch("JetEta"                               , & JetEta                            ,    "JetEta/F"                                 );
  TreeLept->Branch("JetPhi"                               , & JetPhi                            ,    "JetPhi/F"                                 );
  TreeLept->Branch("JetRap"                               , & JetRap                            ,    "JetRap/F"                                 );
  TreeLept->Branch("JetEnergy"                            , & JetEnergy                         ,    "JetEnergy/F"                              );
  TreeLept->Branch("JetMass"                              , & JetMass                           ,    "JetMass/F"                                );
  TreeLept->Branch("JetArea"                              , & JetArea                           ,    "JetArea/F"                                );

  TreeLept->Branch("JetSDmass"                            , & JetSDmass                         ,    "JetSDmass/F"                              );
  TreeLept->Branch("JetSDmassRaw"                         , & JetSDmassRaw                      ,    "JetSDmassRaw/F"                           );
  TreeLept->Branch("JetSDmassCorrL23"                     , & JetSDmassCorrL23                  ,    "JetSDmassCorrL23/F"                       );
  TreeLept->Branch("JetSDmassCorrL23Up"                   , & JetSDmassCorrL23Up                ,    "JetSDmassCorrL23Up/F"                     );
  TreeLept->Branch("JetSDmassCorrL23Dn"                   , & JetSDmassCorrL23Dn                ,    "JetSDmassCorrL23Dn/F"                     );
  TreeLept->Branch("JetSDmassCorrL123"                    , & JetSDmassCorrL123                 ,    "JetSDmassCorrL123/F"                      );
  TreeLept->Branch("JetSDmassCorrL123Up"                  , & JetSDmassCorrL123Up               ,    "JetSDmassCorrL123Up/F"                    );
  TreeLept->Branch("JetSDmassCorrL123Dn"                  , & JetSDmassCorrL123Dn               ,    "JetSDmassCorrL123Dn/F"                    );
  TreeLept->Branch("JetSDmassCorrL23Smear"                , & JetSDmassCorrL23Smear             ,    "JetSDmassCorrL23Smear/F"                     );
  TreeLept->Branch("JetSDmassCorrL23SmearUp"              , & JetSDmassCorrL23SmearUp           ,    "JetSDmassCorrL23SmearUp/F"                   );
  TreeLept->Branch("JetSDmassCorrL23SmearDn"              , & JetSDmassCorrL23SmearDn           ,    "JetSDmassCorrL23SmearDn/F"                   );
  TreeLept->Branch("JetSDptRaw"                           , & JetSDptRaw                        ,    "JetSDptRaw/F"                             );
  TreeLept->Branch("JetSDptCorrL23"                       , & JetSDptCorrL23                    ,    "JetSDptCorrL23/F"                         );
  TreeLept->Branch("JetSDptCorrL23Up"                     , & JetSDptCorrL23Up                  ,    "JetSDptCorrL23Up/F"                       );
  TreeLept->Branch("JetSDptCorrL23Dn"                     , & JetSDptCorrL23Dn                  ,    "JetSDptCorrL23Dn/F"                       );
  TreeLept->Branch("JetSDptCorrL123"                      , & JetSDptCorrL123                   ,    "JetSDptCorrL123/F"                        );
  TreeLept->Branch("JetSDptCorrL123Up"                    , & JetSDptCorrL123Up                 ,    "JetSDptCorrL123Up/F"                      );
  TreeLept->Branch("JetSDptCorrL123Dn"                    , & JetSDptCorrL123Dn                 ,    "JetSDptCorrL123Dn/F"                      );
  TreeLept->Branch("JetSDptCorrL23Smear"                  , & JetSDptCorrL23Smear               ,    "JetSDptCorrL23Smear/F"                       );
  TreeLept->Branch("JetSDptCorrL23SmearUp"                , & JetSDptCorrL23SmearUp             ,    "JetSDptCorrL23SmearUp/F"                     );
  TreeLept->Branch("JetSDptCorrL23SmearDn"                , & JetSDptCorrL23SmearDn             ,    "JetSDptCorrL23SmearDn/F"                     );
  TreeLept->Branch("JetSDetaRaw"                          , & JetSDetaRaw                       ,    "JetSDetaRaw/F"                            );
  TreeLept->Branch("JetSDphiRaw"                          , & JetSDphiRaw                       ,    "JetSDphiRaw/F"                            );

  TreeLept->Branch("JetMassPruned"                        , & JetMassPruned                     ,    "JetMassPruned/F"                          );
  TreeLept->Branch("JetMassTrimmed"                       , & JetMassTrimmed                    ,    "JetMassTrimmed/F"                         );
  TreeLept->Branch("JetTau1"                              , & JetTau1                           ,    "JetTau1/F"                                );
  TreeLept->Branch("JetTau2"                              , & JetTau2                           ,    "JetTau2/F"                                );
  TreeLept->Branch("JetTau3"                              , & JetTau3                           ,    "JetTau3/F"                                );
  TreeLept->Branch("JetTau4"                              , & JetTau4                           ,    "JetTau4/F"                                );
  TreeLept->Branch("JetTau32"                             , & JetTau32                          ,    "JetTau32/F"                               );
  TreeLept->Branch("JetTau21"                             , & JetTau21                          ,    "JetTau21/F"                               );
  TreeLept->Branch("JetSDmaxbdisc"                        , & JetSDmaxbdisc                     ,    "JetSDmaxbdisc/F"                          );
  TreeLept->Branch("JetSDmaxbdiscflavHadron"              , & JetSDmaxbdiscflavHadron           ,    "JetSDmaxbdiscflavHadron/F"                );
  TreeLept->Branch("JetSDmaxbdiscflavParton"              , & JetSDmaxbdiscflavParton           ,    "JetSDmaxbdiscflavParton/F"                );

  TreeLept->Branch("JetSDsubjet0pt"                       , & JetSDsubjet0pt                    ,    "JetSDsubjet0pt/F"                         );
  TreeLept->Branch("JetSDsubjet0mass"                     , & JetSDsubjet0mass                  ,    "JetSDsubjet0mass/F"                       );
  TreeLept->Branch("JetSDsubjet0eta"                      , & JetSDsubjet0eta                   ,    "JetSDsubjet0eta/F"                        );
  TreeLept->Branch("JetSDsubjet0phi"                      , & JetSDsubjet0phi                   ,    "JetSDsubjet0phi/F"                        );
  TreeLept->Branch("JetSDsubjet0area"                     , & JetSDsubjet0area                  ,    "JetSDsubjet0area/F"                       );
  TreeLept->Branch("JetSDsubjet0flavHadron"               , & JetSDsubjet0flavHadron            ,    "JetSDsubjet0flavHadron/F"                 );
  TreeLept->Branch("JetSDsubjet0flavParton"               , & JetSDsubjet0flavParton            ,    "JetSDsubjet0flavParton/F"                 );
  TreeLept->Branch("JetSDsubjet0tau1"                     , & JetSDsubjet0tau1                  ,    "JetSDsubjet0tau1/F"                       );
  TreeLept->Branch("JetSDsubjet0tau2"                     , & JetSDsubjet0tau2                  ,    "JetSDsubjet0tau2/F"                       );
  TreeLept->Branch("JetSDsubjet0tau3"                     , & JetSDsubjet0tau3                  ,    "JetSDsubjet0tau3/F"                       );
  TreeLept->Branch("JetSDsubjet0bdisc"                    , & JetSDsubjet0bdisc                 ,    "JetSDsubjet0bdisc/F"                      );
  TreeLept->Branch("JetSDsubjet1pt"                       , & JetSDsubjet1pt                    ,    "JetSDsubjet1pt/F"                         );
  TreeLept->Branch("JetSDsubjet1mass"                     , & JetSDsubjet1mass                  ,    "JetSDsubjet1mass/F"                       );
  TreeLept->Branch("JetSDsubjet1eta"                      , & JetSDsubjet1eta                   ,    "JetSDsubjet1eta/F"                        );
  TreeLept->Branch("JetSDsubjet1phi"                      , & JetSDsubjet1phi                   ,    "JetSDsubjet1phi/F"                        );
  TreeLept->Branch("JetSDsubjet1area"                     , & JetSDsubjet1area                  ,    "JetSDsubjet1area/F"                       );
  TreeLept->Branch("JetSDsubjet1flavHadron"               , & JetSDsubjet1flavHadron            ,    "JetSDsubjet1flavHadron/F"                 );
  TreeLept->Branch("JetSDsubjet1flavParton"               , & JetSDsubjet1flavParton            ,    "JetSDsubjet1flavParton/F"                 );
  TreeLept->Branch("JetSDsubjet1tau1"                     , & JetSDsubjet1tau1                  ,    "JetSDsubjet1tau1/F"                       );
  TreeLept->Branch("JetSDsubjet1tau2"                     , & JetSDsubjet1tau2                  ,    "JetSDsubjet1tau2/F"                       );
  TreeLept->Branch("JetSDsubjet1tau3"                     , & JetSDsubjet1tau3                  ,    "JetSDsubjet1tau3/F"                       );
  TreeLept->Branch("JetSDsubjet1bdisc"                    , & JetSDsubjet1bdisc                 ,    "JetSDsubjet1bdisc/F"                      );

  TreeLept->Branch("JetPuppiP"                            , & JetPuppiP                         ,    "JetPuppiP/F"                              );
  TreeLept->Branch("JetPuppiPt"                           , & JetPuppiPt                        ,    "JetPuppiPt/F"                             );
  TreeLept->Branch("JetPuppiEta"                          , & JetPuppiEta                       ,    "JetPuppiEta/F"                            );
  TreeLept->Branch("JetPuppiPhi"                          , & JetPuppiPhi                       ,    "JetPuppiPhi/F"                            );
  TreeLept->Branch("JetPuppiMass"                         , & JetPuppiMass                      ,    "JetPuppiMass/F"                           );


  TreeLept->Branch("JetPuppiSDmass"                         , & JetPuppiSDmass                    ,    "JetPuppiSDmass/F"                          );
  TreeLept->Branch("JetPuppiSDmassCorr"                     , & JetPuppiSDmassCorr                ,    "JetPuppiSDmassCorr/F"                      );
  TreeLept->Branch("JetPuppiSDmassCorrUp"                   , & JetPuppiSDmassCorrUp              ,    "JetPuppiSDmassCorrUp/F"                    );
  TreeLept->Branch("JetPuppiSDmassCorrDn"                   , & JetPuppiSDmassCorrDn              ,    "JetPuppiSDmassCorrDn/F"                    );
  TreeLept->Branch("JetPuppiSDmassCorrL23Smear"             , & JetPuppiSDmassCorrL23Smear        ,    "JetPuppiSDmassCorrL23Smear/F"              );
  TreeLept->Branch("JetPuppiSDmassCorrL23SmearUp"           , & JetPuppiSDmassCorrL23SmearUp      ,    "JetPuppiSDmassCorrL23SmearUp/F"            );
  TreeLept->Branch("JetPuppiSDmassCorrL23SmearDn"           , & JetPuppiSDmassCorrL23SmearDn      ,    "JetPuppiSDmassCorrL23SmearDn/F"            );
  TreeLept->Branch("JetPuppiSDpt"                           , & JetPuppiSDpt                      ,    "JetPuppiSDpt/F"                            );
  TreeLept->Branch("JetPuppiSDptCorr"                       , & JetPuppiSDptCorr                  ,    "JetPuppiSDptCorr/F"                        );
  TreeLept->Branch("JetPuppiSDptCorrUp"                     , & JetPuppiSDptCorrUp                ,    "JetPuppiSDptCorrUp/F"                      );
  TreeLept->Branch("JetPuppiSDptCorrDn"                     , & JetPuppiSDptCorrDn                ,    "JetPuppiSDptCorrDn/F"                      );
  TreeLept->Branch("JetPuppiSDptCorrL23Smear"               , & JetPuppiSDptCorrL23Smear          ,    "JetPuppiSDptCorrL23Smear/F"                );
  TreeLept->Branch("JetPuppiSDptCorrL23SmearUp"             , & JetPuppiSDptCorrL23SmearUp        ,    "JetPuppiSDptCorrL23SmearUp/F"              );
  TreeLept->Branch("JetPuppiSDptCorrL23SmearDn"             , & JetPuppiSDptCorrL23SmearDn        ,    "JetPuppiSDptCorrL23SmearDn/F"              );
  TreeLept->Branch("JetPuppiSDeta"                          , & JetPuppiSDeta                     ,    "JetPuppiSDeta/F"                           );
  TreeLept->Branch("JetPuppiSDphi"                          , & JetPuppiSDphi                     ,    "JetPuppiSDphi/F"                           );


  TreeLept->Branch("JetPuppiTau1"                         , & JetPuppiTau1                      ,    "JetPuppiTau1/F"                           );
  TreeLept->Branch("JetPuppiTau2"                         , & JetPuppiTau2                      ,    "JetPuppiTau2/F"                           );
  TreeLept->Branch("JetPuppiTau3"                         , & JetPuppiTau3                      ,    "JetPuppiTau3/F"                           );
  TreeLept->Branch("JetPuppiTau4"                         , & JetPuppiTau4                      ,    "JetPuppiTau4/F"                           );
  TreeLept->Branch("JetPuppiTau32"                        , & JetPuppiTau32                     ,    "JetPuppiTau32/F"                          );
  TreeLept->Branch("JetPuppiTau21"                        , & JetPuppiTau21                     ,    "JetPuppiTau21/F"                          );

  TreeLept->Branch("JetPuppiSDmaxbdisc"                   , & JetPuppiSDmaxbdisc                ,    "JetPuppiSDmaxbdisc/F"                     );
  TreeLept->Branch("JetPuppiSDmaxbdiscflavHadron"         , & JetPuppiSDmaxbdiscflavHadron      ,    "JetPuppiSDmaxbdiscflavHadron/F"           );
  TreeLept->Branch("JetPuppiSDmaxbdiscflavParton"         , & JetPuppiSDmaxbdiscflavParton      ,    "JetPuppiSDmaxbdiscflavParton/F"           );
  TreeLept->Branch("JetPuppiSDsubjet0pt"                  , & JetPuppiSDsubjet0pt               ,    "JetPuppiSDsubjet0pt/F"                    );
  TreeLept->Branch("JetPuppiSDsubjet0mass"                , & JetPuppiSDsubjet0mass             ,    "JetPuppiSDsubjet0mass/F"                  );
  TreeLept->Branch("JetPuppiSDsubjet0eta"                 , & JetPuppiSDsubjet0eta              ,    "JetPuppiSDsubjet0eta/F"                   );
  TreeLept->Branch("JetPuppiSDsubjet0phi"                 , & JetPuppiSDsubjet0phi              ,    "JetPuppiSDsubjet0phi/F"                   );
  TreeLept->Branch("JetPuppiSDsubjet0area"                , & JetPuppiSDsubjet0area             ,    "JetPuppiSDsubjet0area/F"                  );
  TreeLept->Branch("JetPuppiSDsubjet0flavHadron"          , & JetPuppiSDsubjet0flavHadron       ,    "JetPuppiSDsubjet0flavHadron/F"            );
  TreeLept->Branch("JetPuppiSDsubjet0flavParton"          , & JetPuppiSDsubjet0flavParton       ,    "JetPuppiSDsubjet0flavParton/F"            );
  TreeLept->Branch("JetPuppiSDsubjet0tau1"                , & JetPuppiSDsubjet0tau1             ,    "JetPuppiSDsubjet0tau1/F"                  );
  TreeLept->Branch("JetPuppiSDsubjet0tau2"                , & JetPuppiSDsubjet0tau2             ,    "JetPuppiSDsubjet0tau2/F"                  );
  TreeLept->Branch("JetPuppiSDsubjet0tau3"                , & JetPuppiSDsubjet0tau3             ,    "JetPuppiSDsubjet0tau3/F"                  );
  TreeLept->Branch("JetPuppiSDsubjet0bdisc"               , & JetPuppiSDsubjet0bdisc            ,    "JetPuppiSDsubjet0bdisc/F"                 );
  TreeLept->Branch("JetPuppiSDsubjet1pt"                  , & JetPuppiSDsubjet1pt               ,    "JetPuppiSDsubjet1pt/F"                    );
  TreeLept->Branch("JetPuppiSDsubjet1mass"                , & JetPuppiSDsubjet1mass             ,    "JetPuppiSDsubjet1mass/F"                  );
  TreeLept->Branch("JetPuppiSDsubjet1eta"                 , & JetPuppiSDsubjet1eta              ,    "JetPuppiSDsubjet1eta/F"                   );
  TreeLept->Branch("JetPuppiSDsubjet1phi"                 , & JetPuppiSDsubjet1phi              ,    "JetPuppiSDsubjet1phi/F"                   );
  TreeLept->Branch("JetPuppiSDsubjet1area"                , & JetPuppiSDsubjet1area             ,    "JetPuppiSDsubjet1area/F"                  );
  TreeLept->Branch("JetPuppiSDsubjet1flavHadron"          , & JetPuppiSDsubjet1flavHadron       ,    "JetPuppiSDsubjet1flavHadron/F"            );
  TreeLept->Branch("JetPuppiSDsubjet1flavParton"          , & JetPuppiSDsubjet1flavParton       ,    "JetPuppiSDsubjet1flavParton/F"            );
  TreeLept->Branch("JetPuppiSDsubjet1tau1"                , & JetPuppiSDsubjet1tau1             ,    "JetPuppiSDsubjet1tau1/F"                  );
  TreeLept->Branch("JetPuppiSDsubjet1tau2"                , & JetPuppiSDsubjet1tau2             ,    "JetPuppiSDsubjet1tau2/F"                  );
  TreeLept->Branch("JetPuppiSDsubjet1tau3"                , & JetPuppiSDsubjet1tau3             ,    "JetPuppiSDsubjet1tau3/F"                  );
  TreeLept->Branch("JetPuppiSDsubjet1bdisc"               , & JetPuppiSDsubjet1bdisc            ,    "JetPuppiSDsubjet1bdisc/F"                 );


  TreeLept->Branch("JetCHF"                               , & JetCHF                            ,    "JetCHF/F"                                 );
  TreeLept->Branch("JetNHF"                               , & JetNHF                            ,    "JetNHF/F"                                 );
  TreeLept->Branch("JetCM"                                , & JetCM                             ,    "JetCM/F"                                  );
  TreeLept->Branch("JetNM"                                , & JetNM                             ,    "JetNM/F"                                  );
  TreeLept->Branch("JetNEF"                               , & JetNEF                            ,    "JetNEF/F"                                 );
  TreeLept->Branch("JetCEF"                               , & JetCEF                            ,    "JetCEF/F"                                 );
  TreeLept->Branch("JetMF"                                , & JetMF                             ,    "JetMF/F"                                  );
  TreeLept->Branch("JetMult"                              , & JetMult                           ,    "JetMult/F"                                );
  TreeLept->Branch("JetPuppiCHF"                          , & JetPuppiCHF                       ,    "JetPuppiCHF/F"                            );
  TreeLept->Branch("JetPuppiNHF"                          , & JetPuppiNHF                       ,    "JetPuppiNHF/F"                            );
  TreeLept->Branch("JetPuppiCM"                           , & JetPuppiCM                        ,    "JetPuppiCM/F"                             );
  TreeLept->Branch("JetPuppiNM"                           , & JetPuppiNM                        ,    "JetPuppiNM/F"                             );
  TreeLept->Branch("JetPuppiNEF"                          , & JetPuppiNEF                       ,    "JetPuppiNEF/F"                            );
  TreeLept->Branch("JetPuppiCEF"                          , & JetPuppiCEF                       ,    "JetPuppiCEF/F"                            );
  TreeLept->Branch("JetPuppiMF"                           , & JetPuppiMF                        ,    "JetPuppiMF/F"                             );
  TreeLept->Branch("JetPuppiMult"                         , & JetPuppiMult                      ,    "JetPuppiMult/F"                           );
  TreeLept->Branch("JetMassCorrFactor"                    , & JetMassCorrFactor                 ,    "JetMassCorrFactor/F"                      );
  TreeLept->Branch("JetMassCorrFactorUp"                  , & JetMassCorrFactorUp               ,    "JetMassCorrFactorUp/F"                    );
  TreeLept->Branch("JetMassCorrFactorDn"                  , & JetMassCorrFactorDn               ,    "JetMassCorrFactorDn/F"                    );
  TreeLept->Branch("JetCorrFactor"                        , & JetCorrFactor                     ,    "JetCorrFactor/F"                          );
  TreeLept->Branch("JetCorrFactorUp"                      , & JetCorrFactorUp                   ,    "JetCorrFactorUp/F"                        );
  TreeLept->Branch("JetCorrFactorDn"                      , & JetCorrFactorDn                   ,    "JetCorrFactorDn/F"                        );
  TreeLept->Branch("JetPtSmearFactor"                     , & JetPtSmearFactor                  ,    "JetPtSmearFactor/F"                       );
  TreeLept->Branch("JetPtSmearFactorUp"                   , & JetPtSmearFactorUp                ,    "JetPtSmearFactorUp/F"                     );
  TreeLept->Branch("JetPtSmearFactorDn"                   , & JetPtSmearFactorDn                ,    "JetPtSmearFactorDn/F"                     );
  TreeLept->Branch("JetPuppiMassCorrFactor"               , & JetPuppiMassCorrFactor            ,    "JetPuppiMassCorrFactor/F"                 );
  TreeLept->Branch("JetPuppiMassCorrFactorUp"             , & JetPuppiMassCorrFactorUp          ,    "JetPuppiMassCorrFactorUp/F"               );
  TreeLept->Branch("JetPuppiMassCorrFactorDn"             , & JetPuppiMassCorrFactorDn          ,    "JetPuppiMassCorrFactorDn/F"               );
  TreeLept->Branch("JetPuppiCorrFactor"                   , & JetPuppiCorrFactor                ,    "JetPuppiCorrFactor/F"                     );
  TreeLept->Branch("JetPuppiCorrFactorUp"                 , & JetPuppiCorrFactorUp              ,    "JetPuppiCorrFactorUp/F"                   );
  TreeLept->Branch("JetPuppiCorrFactorDn"                 , & JetPuppiCorrFactorDn              ,    "JetPuppiCorrFactorDn/F"                   );
  TreeLept->Branch("JetPuppiPtSmearFactor"                , & JetPuppiPtSmearFactor             ,    "JetPuppiPtSmearFactor/F"                  );
  TreeLept->Branch("JetPuppiPtSmearFactorUp"              , & JetPuppiPtSmearFactorUp           ,    "JetPuppiPtSmearFactorUp/F"                );
  TreeLept->Branch("JetPuppiPtSmearFactorDn"              , & JetPuppiPtSmearFactorDn           ,    "JetPuppiPtSmearFactorDn/F"                );
  TreeLept->Branch("JetEtaScaleFactor"                    , & JetEtaScaleFactor                 ,    "JetEtaScaleFactor/F"                      );
  TreeLept->Branch("JetPhiScaleFactor"                    , & JetPhiScaleFactor                 ,    "JetPhiScaleFactor/F"                      );
  // TreeLept->Branch("JetMatchedGenJetDR"                   , & JetMatchedGenJetDR                ,    "JetMatchedGenJetDR/F"                     );
  TreeLept->Branch("JetMatchedGenJetPt"                   , & JetMatchedGenJetPt                ,    "JetMatchedGenJetPt/F"                     );
  TreeLept->Branch("JetMatchedGenJetMass"                 , & JetMatchedGenJetMass              ,    "JetMatchedGenJetMass/F"                   );

  std::cout<<"Setup semi-lept jets in tree"<<std::endl;

  TreeLept->Branch("LeptMETpx"                        , & LeptMETpx                     , "LeptMETpx/F"                  );
  TreeLept->Branch("LeptMETpy"                        , & LeptMETpy                     , "LeptMETpy/F"                  );
  TreeLept->Branch("LeptMETpt"                        , & LeptMETpt                     , "LeptMETpt/F"                  );
  TreeLept->Branch("LeptMETphi"                       , & LeptMETphi                    , "LeptMETphi/F"                 );
  TreeLept->Branch("LeptMETsumET"                     , & LeptMETsumET                  , "LeptMETsumET/F"               );
  TreeLept->Branch("LeptNvtx"                         , & LeptNvtx                      , "LeptNvtx/F"                   );
  TreeLept->Branch("LeptRho"                          , & LeptRho                       , "LeptRho/F"                    );
  TreeLept->Branch("LeptEventWeight"                  , & LeptEventWeight               , "LeptEventWeight/F"            );
  TreeLept->Branch("LeptPUweight"                  , & LeptPUweight               , "LeptPUweight/F"            );
  TreeLept->Branch("LeptPUweight_MBup"                  , & LeptPUweight_MBup               , "LeptPUweight_MBup/F"            );
  TreeLept->Branch("LeptPUweight_MBdn"                  , & LeptPUweight_MBdn               , "LeptPUweight_MBdn/F"            );


  TreeLept->Branch("HTlep0"                                , & HTlep0                             , "HTlep0/F"                  );
  TreeLept->Branch("ST0"                                   , & ST0                                , "ST0/F"                     );
  TreeLept->Branch("ST0_CorrDn"                            , & ST0_CorrDn                         , "ST0_CorrDn/F"              );
  TreeLept->Branch("ST0_CorrUp"                            , & ST0_CorrUp                         , "ST0_CorrUp/F"              );
  TreeLept->Branch("ST0_PtSmearNom"                        , & ST0_PtSmearNom                     , "ST0_PtSmearNom/F"          );
  TreeLept->Branch("ST0_PtSmearUp"                         , & ST0_PtSmearUp                      , "ST0_PtSmearUp/F"           );
  TreeLept->Branch("ST0_PtSmearDn"                         , & ST0_PtSmearDn                      , "ST0_PtSmearDn/F"           );

  TreeLept->Branch("HTlep1"                                , & HTlep1                             , "HTlep1/F"                  );
  TreeLept->Branch("ST1"                                   , & ST1                                , "ST1/F"                     );
  TreeLept->Branch("ST1_CorrDn"                            , & ST1_CorrDn                         , "ST1_CorrDn/F"              );
  TreeLept->Branch("ST1_CorrUp"                            , & ST1_CorrUp                         , "ST1_CorrUp/F"              );
  TreeLept->Branch("ST1_PtSmearNom"                        , & ST1_PtSmearNom                     , "ST1_PtSmearNom/F"          );
  TreeLept->Branch("ST1_PtSmearUp"                         , & ST1_PtSmearUp                      , "ST1_PtSmearUp/F"           );
  TreeLept->Branch("ST1_PtSmearDn"                         , & ST1_PtSmearDn                      , "ST1_PtSmearDn/F"           );


  TreeLept->Branch("LeptRunNum"                       , & LeptRunNum                    , "LeptRunNum/F"                 );
  TreeLept->Branch("LeptLumiBlock"                    , & LeptLumiBlock                 , "LeptLumiBlock/F"              );
  TreeLept->Branch("LeptEventNum"                     , & LeptEventNum                  , "LeptEventNum/F"               );
  TreeLept->Branch("LeptPassMETFilters"               , & LeptPassMETFilters            , "LeptPassMETFilters/I"         );

  TreeLept->Branch("LepHemiContainsAK4BtagLoose"          , & LepHemiContainsAK4BtagLoose       , "LepHemiContainsAK4BtagLoose/I"    );
  TreeLept->Branch("LepHemiContainsAK4BtagMedium"         , & LepHemiContainsAK4BtagMedium      , "LepHemiContainsAK4BtagMedium/I"   );
  TreeLept->Branch("LepHemiContainsAK4BtagTight"          , & LepHemiContainsAK4BtagTight       , "LepHemiContainsAK4BtagTight/I"    );

  TreeLept->Branch("Lepton0Phi"                            , &  Lepton0Phi                        , "Lepton0Phi/F"                      );
  TreeLept->Branch("Lepton0Pt"                             , &  Lepton0Pt                         , "Lepton0Pt/F"                       );
  TreeLept->Branch("Lepton0Eta"                            , &  Lepton0Eta                        , "Lepton0Eta/F"                      );
  TreeLept->Branch("Lepton0Mass"                           , &  Lepton0Mass                       , "Lepton0Mass/F"                     );
  TreeLept->Branch("PtRel"                                , &  PtRel                            , "PtRel/F"                          );
  TreeLept->Branch("Mu0Medium"                             , &  Mu0Medium                         , "Mu0Medium/I"                       );
  TreeLept->Branch("Mu0Tight"                              , &  Mu0Tight                          , "Mu0Tight/I"                        );
  TreeLept->Branch("Mu0HighPt"                             , &  Mu0HighPt                         , "Mu0HighPt/I"                       );
  TreeLept->Branch("DeltaRJetLep0"                         , &  DeltaRJetLep0                     , "DeltaRJetLep0/F"                   );
  TreeLept->Branch("DeltaPhiJetLep0"                       , &  DeltaPhiJetLep0                   , "DeltaPhiJetLep0/F"                 );

  TreeLept->Branch("Mu0Iso"                                , &  Mu0Iso                            , "Mu0Iso/F"                          );
  TreeLept->Branch("Elecron0_absiso"                       , &  Elecron0_absiso                   , "Elecron0_absiso/F"                 );
  TreeLept->Branch("Elecron0_relIsoWithDBeta"              , &  Elecron0_relIsoWithDBeta          , "Elecron0_relIsoWithDBeta/F"        );
  TreeLept->Branch("Elecron0_absiso_EA"                    , &  Elecron0_absiso_EA                , "Elecron0_absiso_EA/F"              );
  TreeLept->Branch("Elecron0_relIsoWithEA"                 , &  Elecron0_relIsoWithEA             , "Elecron0_relIsoWithEA/F"           );

  TreeLept->Branch("Electron0_isMedium"                 , &  Electron0_isMedium             , "Electron0_isMedium/I"           );
  TreeLept->Branch("Electron0_isTight"                  , &  Electron0_isTight              , "Electron0_isTight/I"            );

  TreeLept->Branch("Lepton1Phi"                            , &  Lepton1Phi                        , "Lepton1Phi/F"                      );
  TreeLept->Branch("Lepton1Pt"                             , &  Lepton1Pt                         , "Lepton1Pt/F"                       );
  TreeLept->Branch("Lepton1Eta"                            , &  Lepton1Eta                        , "Lepton1Eta/F"                      );
  TreeLept->Branch("Lepton1Mass"                           , &  Lepton1Mass                       , "Lepton1Mass/F"                     );
  TreeLept->Branch("PtRel"                                , &  PtRel                            , "PtRel/F"                          );
  TreeLept->Branch("LeptonIsMu"                           , &  LeptonIsMu                       , "LeptonIsMu/I"                     );
  TreeLept->Branch("Mu1Medium"                             , &  Mu1Medium                         , "Mu1Medium/I"                       );
  TreeLept->Branch("Mu1Tight"                              , &  Mu1Tight                          , "Mu1Tight/I"                        );
  TreeLept->Branch("Mu1HighPt"                             , &  Mu1HighPt                         , "Mu1HighPt/I"                       );
  TreeLept->Branch("DeltaRJetLep1"                         , &  DeltaRJetLep1                     , "DeltaRJetLep1/F"                   );
  TreeLept->Branch("DeltaPhiJetLep1"                       , &  DeltaPhiJetLep1                   , "DeltaPhiJetLep1/F"                 );


  TreeLept->Branch("Mu1Iso"                                , &  Mu1Iso                            , "Mu1Iso/F"                          );
  TreeLept->Branch("Elecron1_absiso"                       , &  Elecron1_absiso                   , "Elecron1_absiso/F"                 );
  TreeLept->Branch("Elecron1_relIsoWithDBeta"              , &  Elecron1_relIsoWithDBeta          , "Elecron1_relIsoWithDBeta/F"        );
  TreeLept->Branch("Elecron1_absiso_EA"                    , &  Elecron1_absiso_EA                , "Elecron1_absiso_EA/F"              );
  TreeLept->Branch("Elecron1_relIsoWithEA"                 , &  Elecron1_relIsoWithEA             , "Elecron1_relIsoWithEA/F"           );

  TreeLept->Branch("Electron1_isMedium"                 , &  Electron1_isMedium             , "Electron1_isMedium/I"           );
  TreeLept->Branch("Electron1_isTight"                  , &  Electron1_isTight              , "Electron1_isTight/I"            );









  std::cout<<"Setup Leptonic Tree"<<std::endl;

  std::cout<<"Finished constructor"<<std::endl;

}


B2GTTbarTreeMaker::~B2GTTbarTreeMaker()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
B2GTTbarTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  using namespace LHAPDF;

  if (verbose_) {
    cout<<"----------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Analyze event "<<iEvent.id().event()<<" run "<<iEvent.id().run()<<" lumiblock "<<iEvent.id().luminosityBlock()<<endl;
  }
  //
  //    .d8888b.  8888888888 888b    888     8888888b.                   888    d8b          888
  //   d88P  Y88b 888        8888b   888     888   Y88b                  888    Y8P          888
  //   888    888 888        88888b  888     888    888                  888                 888
  //   888        8888888    888Y88b 888     888   d88P  8888b.  888d888 888888 888  .d8888b 888  .d88b.  .d8888b
  //   888  88888 888        888 Y88b888     8888888P"      "88b 888P"   888    888 d88P"    888 d8P  Y8b 88K
  //   888    888 888        888  Y88888     888        .d888888 888     888    888 888      888 88888888 "Y8888b.
  //   Y88b  d88P 888        888   Y8888     888        888  888 888     Y88b.  888 Y88b.    888 Y8b.          X88
  //    "Y8888P88 8888888888 888    Y888     888        "Y888888 888      "Y888 888  "Y8888P 888  "Y8888   88888P'
  //

  bool GenTruth_Leptonic = false;

  if (!iEvent.isRealData() and runGenLoop_) {
    Handle<edm::View<reco::GenParticle> > genpart;
    iEvent.getByToken(prunedGenToken_,genpart);
  }

  //
  // 888    888 888      88888888888
  // 888    888 888          888
  // 888    888 888          888
  // 8888888888 888          888
  // 888    888 888          888
  // 888    888 888          888
  // 888    888 888          888
  // 888    888 88888888     888
  //

    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    iEvent.getByToken(triggerPrescales_, triggerPrescales);

    vector<string> trigsToRun;
    trigsToRun.push_back("HLT_PFHT300_v");
    trigsToRun.push_back("HLT_PFHT350_v");
    trigsToRun.push_back("HLT_PFHT400_v");
    trigsToRun.push_back("HLT_PFHT475_v");
    trigsToRun.push_back("HLT_PFHT600_v");
    trigsToRun.push_back("HLT_PFHT650_v");
    trigsToRun.push_back("HLT_PFHT800_v");
    trigsToRun.push_back("HLT_PFHT900_v");
    trigsToRun.push_back("HLT_PFJet320_v");
    trigsToRun.push_back("HLT_PFJet400_v");
    trigsToRun.push_back("HLT_PFJet450_v");
    trigsToRun.push_back("HLT_PFJet500_v");
    trigsToRun.push_back("HLT_AK8PFJet360_TrimMass30_v");
    trigsToRun.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v");
    trigsToRun.push_back("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v");
    trigsToRun.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v");
    trigsToRun.push_back("HLT_Mu45_eta2p1_v");
    trigsToRun.push_back("HLT_Mu50_v");
    trigsToRun.push_back("HLT_Mu55_v");
    trigsToRun.push_back("HLT_IsoMu22_eta2p1_v");
    trigsToRun.push_back("HLT_IsoMu24_v");
    trigsToRun.push_back("HLT_IsoMu27_v");
    trigsToRun.push_back("HLT_Mu30_eta2p1_PFJet150_PFJet50_v");
    trigsToRun.push_back("HLT_Mu40_eta2p1_PFJet200_PFJet50_v");
    trigsToRun.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");
    trigsToRun.push_back("HLT_Ele35_WPLoose_Gsf_v");
    trigsToRun.push_back("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v");
    trigsToRun.push_back("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140_v");
    trigsToRun.push_back("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v");
    trigsToRun.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");

    const int ntrigs = trigsToRun.size();
    if (verbose_) cout<<"trigsToRun size "<<ntrigs<<endl;

    // do the same thing two different ways ( to test)
    std::bitset<30> hltbit;
    vector<bool> trigAccept;

    LeptTrigNames     ->clear();
    LeptTrigPrescales ->clear();
    LeptTrigPass      ->clear();

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    if (verbose_) std::cout << "\n === TRIGGER PATHS === " << std::endl;
    int counttrigs =0;
    // for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    for (unsigned int j=0; j<trigsToRun.size(); j++){  // Running this loop first even though it is slow in order to preserve order (temporary solution)
      if (verbose_) cout<<"try to find "<< trigsToRun[j]<<endl;
      for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        string name = names.triggerName(i);
        //if (verbose_) cout<<" "<<name<<endl;
        std::size_t found = name.find( trigsToRun[j] );
        // cout<<" Check: "<<trigsToRun[j]  <<" = "<<name<< " ?" <<found<<endl;
        if ( found !=std::string::npos ) {
          int accept = triggerBits->accept(i) ;
          bool pass = false;
          if (accept ==1 ) pass = true;
          int prescale = triggerPrescales->getPrescaleForIndex(i)  ;

          if (verbose_) std::cout << "  Found Trigger " << names.triggerName(i) << ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") << std::endl;
          trigAccept.push_back(pass);
          // trigStrings.push_back(name);
          // trigPrescale.push_back(prescale);
          AllHadTrigNames       ->push_back(name);
          LeptTrigNames     ->push_back(name);
          AllHadTrigPrescales   ->push_back(prescale);
          LeptTrigPrescales ->push_back(prescale);
          AllHadTrigPass        ->push_back(pass);
          LeptTrigPass      ->push_back(pass);
          if (pass)  hltbit[counttrigs]=1;
          counttrigs++;
          break;
        }
      }
    }

    if (verbose_) {
      cout<<"trig accept size "<<trigAccept.size()<<endl;
      for (unsigned int i=0; i< trigAccept.size(); i++){
        cout<<trigAccept[trigAccept.size()-1-i];
      }
      cout<<endl;
      cout<<"hlt bit"<<endl;
      cout<<hltbit.to_string()<<endl;
    }

    LeptTrigAcceptBits = hltbit.to_string();


  //
  // 888b     d888 8888888888 88888888888     888b    888          d8b                       8888888888 d8b 888 888
  // 8888b   d8888 888            888         8888b   888          Y8P                       888        Y8P 888 888
  // 88888b.d88888 888            888         88888b  888                                    888            888 888
  // 888Y88888P888 8888888        888         888Y88b 888  .d88b.  888 .d8888b   .d88b.      8888888    888 888 888888  .d88b.  888d888 .d8888b
  // 888 Y888P 888 888            888         888 Y88b888 d88""88b 888 88K      d8P  Y8b     888        888 888 888    d8P  Y8b 888P"   88K
  // 888  Y8P  888 888            888         888  Y88888 888  888 888 "Y8888b. 88888888     888        888 888 888    88888888 888     "Y8888b.
  // 888   "   888 888            888         888   Y8888 Y88..88P 888      X88 Y8b.         888        888 888 Y88b.  Y8b.     888          X88
  // 888       888 8888888888     888         888    Y888  "Y88P"  888  88888P'  "Y8888      888        888 888  "Y888  "Y8888  888      88888P'

  bool passMETfilters=false;
  if ( iEvent.isRealData() ) {
    edm::Handle < edm::TriggerResults > metFilters;
    iEvent.getByToken(triggerResultsMETFilterToken_, metFilters);
    edm::TriggerNames const& filterTriggerNames = iEvent.triggerNames(*metFilters);

    int nMETfilters = metFilters->size();
    if (verbose_) cout<<"nMETfilters "<<nMETfilters<<endl;

    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#MiniAOD_8011_ICHEP_dataset
    // Flag_HBHENoiseFilter TO BE USED
    // Flag_HBHENoiseIsoFilter TO BE USED
    // Flag_EcalDeadCellTriggerPrimitiveFilter TO BE USED
    // Flag_goodVertices TO BE USED
    // Flag_eeBadScFilter TO BE USED
    // Flag_globalTightHalo2016Filter NEW TO BE USED
    // badMuon (NEW not available in 80X miniAOD v2, see snippet below) TO BE USED
    // badCharged hadron (NEW not available in 80X miniAOD v2, see snippet below) TO BE USED

    // Recomendation:
    // primary vertex filter (available in miniAOD v2)
    // beam halo filter (NEW available in miniAOD v2)
    // HBHE noise filter (available in miniAOD v2)
    // HBHEiso noise filter (available in miniAOD v2)
    // ee badSC noise filter (available in miniAOD v2)
    // ECAL TP filter (available in miniAOD v2)
    // badMuon (NEW not available in miniAOD v2, see snippet below)
    // badCharged hadron (NEW not available in miniAOD v2, see snippet below)

    vector <string> filterFlags;
    filterFlags.push_back("Flag_HBHENoiseFilter");
    filterFlags.push_back("Flag_HBHENoiseIsoFilter");
    filterFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
    filterFlags.push_back("Flag_goodVertices");
    filterFlags.push_back("Flag_eeBadScFilter");
    filterFlags.push_back("Flag_globalTightHalo2016Filter");

    unsigned int count_matched_accept = 0;
    for (int itrig = 0; itrig != nMETfilters; ++itrig){
     std::string trigName = filterTriggerNames.triggerName(itrig);
      bool accept = metFilters->accept(itrig);
      if (verbose_) cout<<trigName<<"  "<<accept;
      if ( std::find( filterFlags.begin(), filterFlags.end(), trigName ) != filterFlags.end() ) {
          if (verbose_) cout<<"  -> matches filterFlags list ("<<trigName<<")"<<endl;
          if (accept) count_matched_accept++;
      }
      else { if (verbose_) cout<<endl;}
    }
    if (verbose_) cout<<"filterFlags.size() "<<filterFlags.size()<<" count_matched_accept "<<count_matched_accept<<endl;
    if ( count_matched_accept==filterFlags.size() ){
      passMETfilters=true;
    }
    if (verbose_) cout<<"miniAOD Flags pass? "<<passMETfilters<<endl;
  }
  else{
    passMETfilters=true;
  }

  // RECO problemes -> apply to both data and MC
  Handle<bool> ifilterbadChCand;
  iEvent.getByToken(badChargedCandidateFilterToken_ , ifilterbadChCand);
  bool  filterbadChCandidate = *ifilterbadChCand;
  if (verbose_) cout <<"filterbadChCandidate "<<filterbadChCandidate<<endl;

  Handle<bool> ifilterbadPFMuon;
  iEvent.getByToken(badMuonFilterToken_, ifilterbadPFMuon);
  bool filterbadPFMuon = *ifilterbadPFMuon;
  if (verbose_) cout <<"filterbadPFMuon "<<filterbadPFMuon<<endl;

  passMETfilters = passMETfilters && filterbadChCandidate && filterbadPFMuon;
  if (verbose_) cout<<"passMETfilters = "<< passMETfilters <<endl;

  //
  //  888888 8888888888  .d8888b.      8888888b.                    888                        888
  //    "88b 888        d88P  Y88b     888   Y88b                   888                        888
  //     888 888        888    888     888    888                   888                        888
  //     888 8888888    888            888   d88P  8888b.  888  888 888  .d88b.   8888b.   .d88888 .d8888b
  //     888 888        888            8888888P"      "88b 888  888 888 d88""88b     "88b d88" 888 88K
  //     888 888        888    888     888        .d888888 888  888 888 888  888 .d888888 888  888 "Y8888b.
  //     88P 888        Y88b  d88P     888        888  888 Y88b 888 888 Y88..88P 888  888 Y88b 888      X88
  //     888 8888888888  "Y8888P"      888        "Y888888  "Y88888 888  "Y88P"  "Y888888  "Y88888  88888P'
  //   .d88P                                                    888
  // .d88P"                                                Y8b d88P
  //888P"                                                   "Y88P"
  //

  // AK4chs JEC
  std::vector<JetCorrectorParameters> vParAK4chs;
   for ( std::vector<std::string>::const_iterator ipayload = jecPayloadsAK4chs_.begin(),
     ipayloadEnd = jecPayloadsAK4chs_.end(); ipayload != ipayloadEnd - 1; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vParAK4chs.push_back(pars);
  }
  JetCorrectorAK4chs   = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4chs) );
  JetCorrUncertAK4chs  = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadsAK4chs_.back()));

  // AK8chs JEC
  std::vector<JetCorrectorParameters> vParAK8chs;
   for ( std::vector<std::string>::const_iterator ipayload = jecPayloadsAK8chs_.begin(),
     ipayloadEnd = jecPayloadsAK8chs_.end(); ipayload != ipayloadEnd - 1; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vParAK8chs.push_back(pars);
  }
  JetCorrectorAK8chs   = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK8chs) );
  JetCorrUncertAK8chs  = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadsAK8chs_.back()));

  // AK4pup JEC
  std::vector<JetCorrectorParameters> vParAK4pup;
   for ( std::vector<std::string>::const_iterator ipayload = jecPayloadsAK4pup_.begin(),
     ipayloadEnd = jecPayloadsAK4pup_.end(); ipayload != ipayloadEnd - 1; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vParAK4pup.push_back(pars);
  }
  JetCorrectorAK4pup   = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4pup) );
  JetCorrUncertAK4pup  = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadsAK4pup_.back()));

  // AK8pup JEC
  std::vector<JetCorrectorParameters> vParAK8pup;
   for ( std::vector<std::string>::const_iterator ipayload = jecPayloadsAK8pup_.begin(),
     ipayloadEnd = jecPayloadsAK8pup_.end(); ipayload != ipayloadEnd - 1; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vParAK8pup.push_back(pars);
  }
  JetCorrectorAK8pup   = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK8pup) );
  JetCorrUncertAK8pup  = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadsAK8pup_.back()));

  // jet resolution scale factor from text files
  JME::JetResolutionScaleFactor jer_scaler;
  jer_scaler = JME::JetResolutionScaleFactor(jerSFtext_);


  //
  //  888     888                  888    d8b
  //  888     888                  888    Y8P
  //  888     888                  888
  //  Y88b   d88P  .d88b.  888d888 888888 888  .d8888b  .d88b.  .d8888b
  //   Y88b d88P  d8P  Y8b 888P"   888    888 d88P"    d8P  Y8b 88K
  //    Y88o88P   88888888 888     888    888 888      88888888 "Y8888b.
  //     Y888P    Y8b.     888     Y88b.  888 Y88b.    Y8b.          X88
  //      Y8P      "Y8888  888      "Y888 888  "Y8888P  "Y8888   88888P'
  //
  //

  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(vtxToken_, vertices);
  int nvtx = vertices->size();
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();  // save PV for tight muon ID

  // int nvtxgood =0 ;
  // for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
  //   // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
  //   // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
  //   bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
  //   if ( !isFake &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
  //     nvtxgood++;
  //   }
  // }
  // cout<<"nvtx "<<nvtx<<" nvtxgood "<<nvtxgood<<endl;
  //
  //  8888888b.  888     888     888       888          d8b          888      888
  //  888   Y88b 888     888     888   o   888          Y8P          888      888
  //  888    888 888     888     888  d8b  888                       888      888
  //  888   d88P 888     888     888 d888b 888  .d88b.  888  .d88b.  88888b.  888888
  //  8888888P"  888     888     888d88888b888 d8P  Y8b 888 d88P"88b 888 "88b 888
  //  888        888     888     88888P Y88888 88888888 888 888  888 888  888 888
  //  888        Y88b. .d88P     8888P   Y8888 Y8b.     888 Y88b 888 888  888 Y88b.
  //  888         "Y88888P"      888P     Y888  "Y8888  888  "Y88888 888  888  "Y888
  //                                                             888
  //                                                        Y8b d88P
  //                                                         "Y88P"

  edm::Handle<std::vector<PileupSummaryInfo> > pileup;
  iEvent.getByToken(pileupInfoToken_, pileup);
  int nPU = 0;
  if(pileup.isValid()) { // protection for data
    for(std::vector<PileupSummaryInfo>::const_iterator iPV = pileup->begin(); iPV != pileup->end(); ++iPV) {
      if(iPV->getBunchCrossing() == 0) {
        nPU = iPV->getTrueNumInteractions();
        //  numGenPV = iPV->getPU_NumInteractions();
        break;
      }
    }
  }

  double puweight   = hPUweight     ->GetBinContent( hPUweight     ->GetXaxis()->FindBin( nPU ) )  ;
  double puweightUp = hPUweight_MBup->GetBinContent( hPUweight_MBup->GetXaxis()->FindBin( nPU ) )  ;
  double puweightDn = hPUweight_MBdn->GetBinContent( hPUweight_MBdn->GetXaxis()->FindBin( nPU ) )  ;

  if (verbose_) std::cout<<"PU weight: "<<puweight<<std::endl;

  h_NtrueIntPU ->Fill(nPU);
  h_NPV ->Fill(nvtx);
  h_NPVreweighted  ->Fill(nvtx,puweight);

  //  888      888    888 8888888888     888       888          d8b          888      888
  //  888      888    888 888            888   o   888          Y8P          888      888
  //  888      888    888 888            888  d8b  888                       888      888
  //  888      8888888888 8888888        888 d888b 888  .d88b.  888  .d88b.  88888b.  888888 .d8888b
  //  888      888    888 888            888d88888b888 d8P  Y8b 888 d88P"88b 888 "88b 888    88K
  //  888      888    888 888            88888P Y88888 88888888 888 888  888 888  888 888    "Y8888b.
  //  888      888    888 888            8888P   Y8888 Y8b.     888 Y88b 888 888  888 Y88b.       X88
  //  88888888 888    888 8888888888     888P     Y888  "Y8888  888  "Y88888 888  888  "Y888  88888P'
  //                                                                     888
  //                                                                Y8b d88P
  //                                                                 "Y88P"

/*

  NEED TO FIGURE THESE OUT FOR OUR EVENTS, PREVIOUS WAS FOR TTBAR RSG OR WPRIME
  double Q2wgt_up = -999;
  double Q2wgt_down = -999;

  double NNPDF3wgt_up = -999;
  double NNPDF3wgt_down = -999;

*/

  //
  // 8888888b.  888
  // 888   Y88b 888
  // 888    888 888
  // 888   d88P 88888b.   .d88b.
  // 8888888P"  888 "88b d88""88b
  // 888 T88b   888  888 888  888
  // 888  T88b  888  888 Y88..88P
  // 888   T88b 888  888  "Y88P"
  //

  Handle<double> rhoH;
  iEvent.getByToken(rhoToken_, rhoH);
  double rho = *rhoH;

  //
  // 888b     d888
  // 8888b   d8888
  // 88888b.d88888
  // 888Y88888P888 888  888  .d88b.  88888b.
  // 888 Y888P 888 888  888 d88""88b 888 "88b
  // 888  Y8P  888 888  888 888  888 888  888
  // 888   "   888 Y88b 888 Y88..88P 888  888
  // 888       888  "Y88888  "Y88P"  888  888
  //

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  TLorentzVector mu0_p4;
  bool mu0_isTight=false;
  bool mu0_isMedium=false;
  bool mu0_isHighPt = false;
  double mu0_iso04=0;

  TLorentzVector mu1_p4;
  bool mu1_isTight=false;
  bool mu1_isMedium=false;
  bool mu1_isHighPt = false;
  double mu1_iso04=0;

  int count_mu=0;


  std::vector<reco::CandidatePtr> muFootprint;

  for (const pat::Muon &mu : *muons) {
      if (mu.pt() < 30 || !mu.isLooseMuon() || fabs( mu.eta() ) > 2.1) continue;


      if (count_mu==0){
        mu0_p4.SetPtEtaPhiM( mu.pt(), mu.eta(), mu.phi(), mu.mass() );
        if ( mu.isTightMuon(PV) ) mu0_isTight = true;
        // ICHEP HIP MEDIUM MUON ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Short_Term_Instructions_for_ICHE
        bool goodGlob = mu.isGlobalMuon() &&
                      mu.globalTrack()->normalizedChi2() < 3 &&
                      mu.combinedQuality().chi2LocalPosition < 12 &&
                      mu.combinedQuality().trkKink < 20;
        bool isMedium = muon::isLooseMuon(mu) &&
                      mu.innerTrack()->validFraction() > 0.49 &&
                      muon::segmentCompatibility(mu) > (goodGlob ? 0.303 : 0.451);
        if ( isMedium ) mu0_isMedium = true;

        bool isHighPt = mu.isHighPtMuon(PV);
        if ( isHighPt ) mu0_isHighPt = true;

        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Accessing_PF_Isolation_from_reco
        double sumChargedHadronPt = mu.pfIsolationR04().sumChargedHadronPt;
        double sumNeutralHadronPt = mu.pfIsolationR04().sumNeutralHadronEt;
        double sumPhotonPt        = mu.pfIsolationR04().sumPhotonEt;
        double sumPUPt            = mu.pfIsolationR04().sumPUPt;
        double pt                 = mu.pt();
        double iso04 = (sumChargedHadronPt+TMath::Max(0.,sumNeutralHadronPt+sumPhotonPt-0.5*sumPUPt))/pt;
        mu0_iso04 = iso04;

        for (unsigned int i = 0, n = mu.numberOfSourceCandidatePtrs(); i < n; ++i) {
          muFootprint.push_back(mu.sourceCandidatePtr(i));
        }



        if (verbose_) cout<<"Muon pT "<<mu.pt()<<" iso04 "<<iso04<<endl;
      }

      if (count_mu==1){
        mu1_p4.SetPtEtaPhiM( mu.pt(), mu.eta(), mu.phi(), mu.mass() );
        if ( mu.isTightMuon(PV) ) mu1_isTight = true;
        // ICHEP HIP MEDIUM MUON ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Short_Term_Instructions_for_ICHE
        bool goodGlob = mu.isGlobalMuon() &&
                      mu.globalTrack()->normalizedChi2() < 3 &&
                      mu.combinedQuality().chi2LocalPosition < 12 &&
                      mu.combinedQuality().trkKink < 20;
        bool isMedium = muon::isLooseMuon(mu) &&
                      mu.innerTrack()->validFraction() > 0.49 &&
                      muon::segmentCompatibility(mu) > (goodGlob ? 0.303 : 0.451);
        if ( isMedium ) mu1_isMedium = true;

        bool isHighPt = mu.isHighPtMuon(PV);
        if ( isHighPt ) mu1_isHighPt = true;

        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Accessing_PF_Isolation_from_reco
        double sumChargedHadronPt = mu.pfIsolationR04().sumChargedHadronPt;
        double sumNeutralHadronPt = mu.pfIsolationR04().sumNeutralHadronEt;
        double sumPhotonPt        = mu.pfIsolationR04().sumPhotonEt;
        double sumPUPt            = mu.pfIsolationR04().sumPUPt;
        double pt                 = mu.pt();
        double iso04 = (sumChargedHadronPt+TMath::Max(0.,sumNeutralHadronPt+sumPhotonPt-0.5*sumPUPt))/pt;
        mu1_iso04 = iso04;

        for (unsigned int i = 0, n = mu.numberOfSourceCandidatePtrs(); i < n; ++i) {
          muFootprint.push_back(mu.sourceCandidatePtr(i));
        }



        if (verbose_) cout<<"Muon pT "<<mu.pt()<<" iso04 "<<iso04<<endl;
      }
      // printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
      // mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
      count_mu++;
  }

  //
  // 8888888888 888                   888
  // 888        888                   888
  // 888        888                   888
  // 8888888    888  .d88b.   .d8888b 888888 888d888  .d88b.  88888b.
  // 888        888 d8P  Y8b d88P"    888    888P"   d88""88b 888 "88b
  // 888        888 88888888 888      888    888     888  888 888  888
  // 888        888 Y8b.     Y88b.    Y88b.  888     Y88..88P 888  888
  // 8888888888 888  "Y8888   "Y8888P  "Y888 888      "Y88P"  888  888
  //
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);

  TLorentzVector el0_p4;
  Float_t el0_absiso           =0;
  Float_t el0_relIsoWithDBeta  =0;
  Float_t el0_absiso_EA        =0;
  Float_t el0_relIsoWithEA     =0;
  bool el0_isMedium     = false;
  bool el0_isTight      = false;

  TLorentzVector el1_p4;
  Float_t el1_absiso           =0;
  Float_t el1_relIsoWithDBeta  =0;
  Float_t el1_absiso_EA        =0;
  Float_t el1_relIsoWithEA     =0;
  bool el1_isMedium     = false;
  bool el1_isTight      = false;

  int count_el=0;

  std::vector<reco::CandidatePtr> elFootprint;



  for (const pat::Electron &el : *electrons) {
      if (el.pt() < 50 || fabs(el.eta())>2.4 ) continue;

      float eta = el.eta();
      // Implemation of loose Quality cuts (need to study this and improve)
      // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
      // Recomended isolation variables are stored separately and not included in the loose quality cut


      bool passLoose      = false;
      bool passMedium     = false;
      bool passTight      = false;

      float ooEmooP_;
      if( el.ecalEnergy() == 0 ){
        if(verbose_) printf("Electron energy is zero!\n");
        ooEmooP_ = 999;
      }else if( !std::isfinite(el.ecalEnergy())){
        if(verbose_) printf("Electron energy is not finite!\n");
        ooEmooP_ = 998;
      }else{
        ooEmooP_ = fabs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy() );
      }
      float missHits = el.gsfTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);


      GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
      float absiso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
      float relIsoWithDBeta = absiso/el.pt();

      float effArea = 0.;
      if(abs(eta)>0.0 && abs(eta)<=1.0) effArea = 0.1752;
      if(abs(eta)>1.0 && abs(eta)<=1.479) effArea = 0.1862;
      if(abs(eta)>1.479 && abs(eta)<=2.0) effArea = 0.1411;
      if(abs(eta)>2.0 && abs(eta)<=2.2) effArea = 0.1534;
      if(abs(eta)>2.2 && abs(eta)<=2.3) effArea = 0.1903;
      if(abs(eta)>2.3 && abs(eta)<=2.4) effArea = 0.2243;
      if(abs(eta)>2.4 && abs(eta)<=2.5) effArea = 0.2687;

      float absiso_EA = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * effArea );
      float relIsoWithEA = absiso_EA/el.pt();


      if (fabs(el.eta())<=1.479 ){
        if( el.full5x5_sigmaIetaIeta()                <      0.011      &&
            fabs(el.deltaEtaSuperClusterTrackAtVtx()) <      0.00477    &&
            fabs(el.deltaPhiSuperClusterTrackAtVtx()) <      0.222      &&
            el.hadronicOverEm()                       <      0.298      &&
            ooEmooP_                                  <      0.241      &&
            missHits                                  <=     1           ){
          passLoose =true ;
        }
        if( el.full5x5_sigmaIetaIeta()                <      0.00998     &&
            fabs(el.deltaEtaSuperClusterTrackAtVtx()) <      0.00311     &&
            fabs(el.deltaPhiSuperClusterTrackAtVtx()) <      0.103       &&
            el.hadronicOverEm()                       <      0.253       &&
            ooEmooP_                                  <      0.134       &&
            missHits                                  <=     1           ){
          passMedium =true ;
        }
        if( el.full5x5_sigmaIetaIeta()                <      0.00998     &&
            fabs(el.deltaEtaSuperClusterTrackAtVtx()) <      0.00308     &&
            fabs(el.deltaPhiSuperClusterTrackAtVtx()) <      0.0816      &&
            el.hadronicOverEm()                       <      0.0414      &&
            ooEmooP_                                  <      0.0129      &&
            missHits                                  <=     1           ){
          passTight =true ;
        }
      }
      else {
        if( el.full5x5_sigmaIetaIeta()                <      0.0314     &&
            fabs(el.deltaEtaSuperClusterTrackAtVtx()) <      0.00868    &&
            fabs(el.deltaPhiSuperClusterTrackAtVtx()) <      0.213      &&
            el.hadronicOverEm()                       <      0.101      &&
            ooEmooP_                                  <      0.14       &&
            missHits                                  <=     1             ){
          passLoose =true ;
        }
        if( el.full5x5_sigmaIetaIeta()                <      0.0298      &&
            fabs(el.deltaEtaSuperClusterTrackAtVtx()) <      0.00609     &&
            fabs(el.deltaPhiSuperClusterTrackAtVtx()) <      0.045       &&
            el.hadronicOverEm()                       <      0.0878      &&
            ooEmooP_                                  <      0.13        &&
            missHits                                  <=     1             ){
          passMedium =true ;
        }
        if( el.full5x5_sigmaIetaIeta()                <      0.0292      &&
            fabs(el.deltaEtaSuperClusterTrackAtVtx()) <      0.00605     &&
            fabs(el.deltaPhiSuperClusterTrackAtVtx()) <      0.0394      &&
            el.hadronicOverEm()                       <      0.0641      &&
            ooEmooP_                                  <      0.0129      &&
            missHits                                  <=     1             ){
          passTight =true ;
        }
      }



      if (!passLoose) continue;

      if (count_el==0){
        el0_p4.SetPtEtaPhiM( el.pt(), el.eta(), el.phi(), el.mass() );

        el0_absiso           = absiso;
        el0_relIsoWithDBeta  = relIsoWithDBeta;
        el0_absiso_EA        = absiso_EA;
        el0_relIsoWithEA     = relIsoWithEA;
        el0_isMedium         = passMedium;
        el0_isTight          = passTight ;

        for (unsigned int i = 0, n = el.numberOfSourceCandidatePtrs(); i < n; ++i) {
          elFootprint.push_back(el.sourceCandidatePtr(i));
        }


        if (verbose_) cout<<"Electron pT "<<el.pt()<<endl;
      }
      if (count_el==1){
        el1_p4.SetPtEtaPhiM( el.pt(), el.eta(), el.phi(), el.mass() );

        el1_absiso           = absiso;
        el1_relIsoWithDBeta  = relIsoWithDBeta;
        el1_absiso_EA        = absiso_EA;
        el1_relIsoWithEA     = relIsoWithEA;
        el1_isMedium         = passMedium;
        el1_isTight          = passTight ;

        for (unsigned int i = 0, n = el.numberOfSourceCandidatePtrs(); i < n; ++i) {
          elFootprint.push_back(el.sourceCandidatePtr(i));
        }

        if (verbose_) cout<<"Electron pT "<<el.pt()<<endl;
      }

      count_el++;

      //printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), lost hits %d, pass conv veto %d\n",
      //              el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(), el.passConversionVeto());
  }

  TLorentzVector lep0_p4;
  TLorentzVector lep1_p4;
  if (count_mu==2 &&  count_el==0){
    lep0_p4 = mu0_p4;
    lep1_p4 = mu1_p4;
  }
  else if (count_el==2 && count_mu ==0){
    lep0_p4 = el0_p4;
    lep1_p4 = el1_p4;
  }

  int count_lep = count_mu + count_el;

  if (verbose_){
    cout<<"count_mu  "<<count_mu<<endl;
    cout<<"count_el  "<<count_el<<endl;
    cout<<"count_lep "<<count_lep<<endl;
  }

  //
  // 888b     d888 8888888888 88888888888
  // 8888b   d8888 888            888
  // 88888b.d88888 888            888
  // 888Y88888P888 8888888        888
  // 888 Y888P 888 888            888
  // 888  Y8P  888 888            888
  // 888   "   888 888            888
  // 888       888 8888888888     888
  //

  if (verbose_) cout<<"debug: about to grab met"<<endl;
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();
  if (verbose_){
    cout<<"MET pt "<<met.pt()<<endl;
    cout<<"MET phi "<<met.phi()<<endl;
    cout<<"MET sumEt "<<met.sumEt()<<endl;
    if (!iEvent.isRealData() )  cout<<"genMET "<< met.genMET()->pt()<<endl;
  }

  //
  //         d8888 888    d8P   .d8888b.       .d8888b.  888    888  .d8888b.         d8b          888
  //        d88888 888   d8P   d88P  Y88b     d88P  Y88b 888    888 d88P  Y88b        Y8P          888
  //       d88P888 888  d8P    Y88b. d88P     888    888 888    888 Y88b.                          888
  //      d88P 888 888d88K      "Y88888"      888        8888888888  "Y888b.         8888  .d88b.  888888 .d8888b
  //     d88P  888 8888888b    .d8P""Y8b.     888        888    888     "Y88b.       "888 d8P  Y8b 888    88K
  //    d88P   888 888  Y88b   888    888     888    888 888    888       "888        888 88888888 888    "Y8888b.
  //   d8888888888 888   Y88b  Y88b  d88P     Y88b  d88P 888    888 Y88b  d88P        888 Y8b.     Y88b.       X88
  //  d88P     888 888    Y88b  "Y8888P"       "Y8888P"  888    888  "Y8888P"         888  "Y8888   "Y888  88888P'
  //                                                                                  888
  //                                                                               d88P
  //                                                                             888P"

  edm::Handle<pat::JetCollection> AK8CHS;
  iEvent.getByToken(ak8jetToken_, AK8CHS);

  edm::Handle<reco::GenJetCollection> AK8GENJET;
  iEvent.getByToken(ak8genjetToken_, AK8GENJET);

  edm::Handle<pat::JetCollection> AK8CHSsub;
  edm::Handle<pat::JetCollection> AK8PUPPI;
  edm::Handle<pat::JetCollection> AK8PUPPIsub;
  if (useToolbox_){
    iEvent.getByToken( ak8CHSSoftDropSubjetsToken_   , AK8CHSsub);
    iEvent.getByToken( puppijetToken_ , AK8PUPPI );
    iEvent.getByToken( ak8PuppiSoftDropSubjetsToken_ , AK8PUPPIsub);
  }

  int count_AK8CHS = 0;
  int count_fill_leptTree =0;

  TLorentzVector AK8jet_Lept_P4corr;
  TLorentzVector AK8jet0_P4corr;
  TLorentzVector AK8jet1_P4corr;
  TLorentzVector PUPPIjet0_P4;
  TLorentzVector PUPPIjet1_P4;
  TLorentzVector PUPPIjet0_P4corr;
  TLorentzVector PUPPIjet1_P4corr;

  TLorentzVector GenJetMatched0;
  TLorentzVector GenJetMatched1;

  if (verbose_) cout<<"debug: about to grab ak8 jets"<<endl;

  for (const pat::Jet &ijet : *AK8CHS) {
    if (count_AK8CHS>1) break;
    if (count_AK8CHS==0 && ijet.pt()<250) break;

    if (verbose_) cout<<"\nJet "<<count_AK8CHS<<" with pT "<<ijet.pt()<<" sdMass "<<ijet.userFloat("ak8PFJetsCHSSoftDropMass")<<endl;

    //------------------------------------
    // Noise jet ID
    //------------------------------------
    double NHF       = ijet.neutralHadronEnergyFraction();
    double NEMF      = ijet.neutralEmEnergyFraction();
    double CHF       = ijet.chargedHadronEnergyFraction();
    // double MUF       = ijet.muonEnergyFraction();
    double CEMF      = ijet.chargedEmEnergyFraction();
    double NumConst  = ijet.chargedMultiplicity()+ijet.neutralMultiplicity();
    double NM        = ijet.neutralMultiplicity();
    double CM        = ijet.chargedMultiplicity();
    double eta       = ijet.eta();

    bool goodJet_looseJetID =
         ( fabs(eta) <= 2.4 && NHF < 0.99 && NEMF < 0.99 && NumConst >1 && CHF > 0.0  && CM > 0 && CEMF < 0.99   )
      || ( fabs(eta) <= 2.7 && fabs(eta) > 2.4 && NHF < 0.99 && NEMF < 0.99 && NumConst >1 )
      || ( fabs(eta) <= 3.0 && fabs(eta) > 2.7 && NEMF < 0.9 && NM > 2 )
      || ( fabs(eta)  > 3.0 && NEMF < 0.9 && NM > 10 );
    if (verbose_) cout<<"  goodJet = "<<goodJet_looseJetID<<endl;

    if (!goodJet_looseJetID) {
      if(verbose_) cout<<"  bad AK8 jet. skip.  ( pt "<<ijet.pt()<<" eta "<<ijet.eta()<<" NumConst "<<NumConst<<" )"<<endl;
      continue;
    }

    //------------------------------------
    // AK8CHS JEC correction
    //------------------------------------
    reco::Candidate::LorentzVector uncorrJet = ijet.correctedP4(0);
    JetCorrectorAK8chs->setJetEta( uncorrJet.eta() );
    JetCorrectorAK8chs->setJetPt ( uncorrJet.pt() );
    JetCorrectorAK8chs->setJetE  ( uncorrJet.energy() );
    JetCorrectorAK8chs->setJetA  ( ijet.jetArea() );
    JetCorrectorAK8chs->setRho   ( rho );
    JetCorrectorAK8chs->setNPV   ( nvtx );
    double corr = JetCorrectorAK8chs->getCorrection();

    reco::Candidate::LorentzVector corrJet = corr * uncorrJet;
    if (verbose_) cout<<" uncorrected AK8 jet pt "<<uncorrJet.pt()<<" corrected jet pt "<<corrJet.pt()<<endl;

    if(count_AK8CHS==0) AK8jet0_P4corr.SetPtEtaPhiM( corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );
    if(count_AK8CHS==1) AK8jet1_P4corr.SetPtEtaPhiM( corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );

    //------------------------------------
    // AK8CHS JEC L23 correction
    //------------------------------------
    JetCorrectorAK8chs->setJetEta( uncorrJet.eta() );
    JetCorrectorAK8chs->setJetPt ( uncorrJet.pt() );
    JetCorrectorAK8chs->setJetE  ( uncorrJet.energy() );
    JetCorrectorAK8chs->setJetA  ( ijet.jetArea() );
    JetCorrectorAK8chs->setRho   ( rho );
    JetCorrectorAK8chs->setNPV   ( nvtx );
    // getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction.
    vector<float> factors = JetCorrectorAK8chs->getSubCorrections();
    float corr_factor_L1      = 1.0;
    float corr_factor_L12     = 1.0;
    float corr_factor_L123    = 1.0;
    float corr_factor_L123res = 1.0;
    if (factors.size() > 0) corr_factor_L1       = factors[0];
    if (factors.size() > 1) corr_factor_L12      = factors[1];
    if (factors.size() > 2) corr_factor_L123     = factors[2];
    if (factors.size() > 3) corr_factor_L123res  = factors[3];
    double corr_factor_L2 = corr_factor_L12/corr_factor_L1;
    double corr_factor_L3 = corr_factor_L123/corr_factor_L12;
    double corr_factor_res = corr_factor_L123res/corr_factor_L123;
    //double corr_factor_L23 = corr_factor_L2*corr_factor_L3;
    double corr_factor_L23res = corr_factor_L2*corr_factor_L3*corr_factor_res;

    //------------------------------------
    // AK8CHS JEC uncertainty
    //------------------------------------
    double corrDn_L23  = 1.0;
    double corrDn_L123 = 1.0;
    JetCorrUncertAK8chs->setJetPhi(  corrJet.phi()  );
    JetCorrUncertAK8chs->setJetEta(  corrJet.eta()  );
    JetCorrUncertAK8chs->setJetPt(   corrJet.pt()   );
    double corrDn_temp1 = JetCorrUncertAK8chs->getUncertainty(0);
    corrDn_L23   = corr_factor_L23res - corrDn_temp1;
    corrDn_L123 = corr - corrDn_temp1;
    double corrUp_L23  = 1.0;
    double corrUp_L123 = 1.0;
    JetCorrUncertAK8chs->setJetPhi(  corrJet.phi()  );
    JetCorrUncertAK8chs->setJetEta(  corrJet.eta()  );
    JetCorrUncertAK8chs->setJetPt(   corrJet.pt()   );
    double corrUp_temp1 = JetCorrUncertAK8chs->getUncertainty(1);
    corrUp_L23   = corr_factor_L23res + corrUp_temp1;
    corrUp_L123 = corr + corrUp_temp1;

    //------------------------------------
    // GenJet  matched
    //------------------------------------
    TLorentzVector GenJetMatched;
    if (!iEvent.isRealData()){
      const reco::GenJet* genJet = ijet.genJet();
      if (genJet) {
        GenJetMatched.SetPtEtaPhiM( genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass() );
        if (verbose_) cout<<"  ak8 genJet pt "<<genJet->pt()<<" mass "<<genJet->mass()<<endl;
      }
    }

    if ( count_AK8CHS==0 ) GenJetMatched0 = GenJetMatched;
    if ( count_AK8CHS==1 ) GenJetMatched1 = GenJetMatched;


    //------------------------------------
    // JER SF
    //------------------------------------
    double ptsmear   = 1;
    double ptsmearUp = 1;
    double ptsmearDn = 1;
    if (!iEvent.isRealData()) {
      double jer_sf    = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrJet.eta()}});
      double jer_sf_up = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrJet.eta()}}, Variation::UP);
      double jer_sf_dn = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrJet.eta()}}, Variation::DOWN);
      if (verbose_) std::cout << " JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn << std::endl;
      double recopt    = corrJet.pt();
      double genpt     = GenJetMatched.Perp();
      double deltapt   = (recopt-genpt)*(jer_sf-1.0);
      double deltaptUp = (recopt-genpt)*(jer_sf_up-1.0);
      double deltaptDn = (recopt-genpt)*(jer_sf_dn-1.0);
      ptsmear   = std::max((double)0.0, (recopt+deltapt)/recopt     );
      ptsmearUp = std::max((double)0.0, (recopt+deltaptUp)/recopt   );
      ptsmearDn = std::max((double)0.0, (recopt+deltaptDn)/recopt   );
    }

    //------------------------------------
    // AK8 variables from miniAOD
    //------------------------------------
    // double pt           = corrJet.pt();
    //double mass         = corrJet.mass();
    // double eta          = corrJet.eta();
    // double phi          = corrJet.phi();
    //double rapidity     = ijet.rapidity();
    //double ndau         = ijet.numberOfDaughters();

    double tau1         = 99;
    double tau2         = 99;
    double tau3         = 99;
    double tau4         = 99;
    double prunedMass   = ijet.userFloat("ak8PFJetsCHSPrunedMass");
    double softDropMass = ijet.userFloat("ak8PFJetsCHSSoftDropMass");
    double trimmedMass  = -1;

    if (useToolbox_){
      tau1         = ijet.userFloat("NjettinessAK8CHS:tau1");
      tau2         = ijet.userFloat("NjettinessAK8CHS:tau2");
      tau3         = ijet.userFloat("NjettinessAK8CHS:tau3");
      tau4         = ijet.userFloat("NjettinessAK8CHS:tau4");
      trimmedMass  = ijet.userFloat("ak8PFJetsCHSTrimmedMass");
    }
    else{
      tau1         = ijet.userFloat("NjettinessAK8:tau1");
      tau2         = ijet.userFloat("NjettinessAK8:tau2");
      tau3         = ijet.userFloat("NjettinessAK8:tau3");
    }
    double tau21        = 99;
    double tau32        = 99;

    double puppi_p              = -1;
    double puppi_pt             = -1;
    double puppi_mass           = -1;
    double puppi_eta            = -1;
    double puppi_phi            = -1;
    double puppi_tau1           = -1;
    double puppi_tau2           = -1;
    double puppi_tau3           = -1;
    double puppi_tau4           = -1;
    // double puppi_prunedMass     = -1;
    // double puppi_trimmedMass    = -1;
    // double puppi_softDropMass   = -1;

    double puppi_CHF    = -1;
    double puppi_NHF    = -1;
    double puppi_CM     = -1;
    double puppi_NM     = -1;
    double puppi_NEF    = -1;
    double puppi_CEF    = -1;
    double puppi_MF     = -1;
    double puppi_Mult   = -1;

    if (!useToolbox_){
      puppi_pt           = ijet.userFloat("ak8PFJetsPuppiValueMap:pt");
      puppi_mass         = ijet.userFloat("ak8PFJetsPuppiValueMap:mass");
      puppi_eta          = ijet.userFloat("ak8PFJetsPuppiValueMap:eta");
      puppi_phi          = ijet.userFloat("ak8PFJetsPuppiValueMap:phi");
      puppi_tau1         = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
      puppi_tau2         = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
      puppi_tau3         = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");
    }
    if (useToolbox_){
      for (const pat::Jet &ipup : *AK8PUPPI) {
        if (ipup.pt()<180) break;
        double deltaRpup = deltaR(ijet.eta(), ijet.phi(), ipup.eta(), ipup.phi() );
        if (deltaRpup<0.8){
          puppi_p            = ipup.p();
          puppi_pt           = ipup.pt();
          puppi_mass         = ipup.mass();
          puppi_eta          = ipup.eta();
          puppi_phi          = ipup.phi();
          // puppi_prunedMass   = ipup.userFloat("ak8PFJetsPuppiPrunedMass");
          // puppi_trimmedMass  = ipup.userFloat("ak8PFJetsPuppiTrimmedMass");
          // puppi_softDropMass = ipup.userFloat("ak8PFJetsPuppiSoftDropMass");
          puppi_tau1         = ipup.userFloat("NjettinessAK8Puppi:tau1");
          puppi_tau2         = ipup.userFloat("NjettinessAK8Puppi:tau2");
          puppi_tau3         = ipup.userFloat("NjettinessAK8Puppi:tau3");
          puppi_tau4         = ipup.userFloat("NjettinessAK8Puppi:tau4");

          puppi_CHF          = ipup.chargedHadronEnergy() / ipup.correctedP4(0).E()  ;
          puppi_NHF          = ipup.neutralHadronEnergy() / ipup.correctedP4(0).E()  ;
          puppi_CM           = ipup.chargedMultiplicity()  ;
          puppi_NM           = ipup.neutralMultiplicity()  ;
          puppi_NEF          = ipup.neutralEmEnergy() / ipup.correctedP4(0).E()  ;
          puppi_CEF          = ipup.chargedEmEnergy() / ipup.correctedP4(0).E()  ;
          puppi_MF           = ipup.muonEnergy() / ipup.correctedP4(0).E()  ;
          puppi_Mult         = ipup.numberOfDaughters() ;

        }
      }
    }
    double puppi_tau21        = 99;
    double puppi_tau32        = 99;


    if (tau1!=0) tau21 = tau2/tau1;
    if (tau2!=0) tau32 = tau3/tau2;

    if (puppi_tau1!=0) puppi_tau21 = puppi_tau2/puppi_tau1;
    if (puppi_tau2!=0) puppi_tau32 = puppi_tau3/puppi_tau2;

    TLorentzVector jet_p4;
    jet_p4.SetPtEtaPhiM( corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );


    // //------------------------------------
    // // AK8PUPPI JEC L23 correction
    // //------------------------------------

    // TLorentzVector jet_p4;
    // pupppijet_p4.SetPtEtaPhiM( puppi_pt, puppi_eta, puppi_phi, puppi_mass);


    // JetCorrectorAK8chs->setJetEta( puppi_eta );
    // JetCorrectorAK8chs->setJetPt ( puppi_pt );
    // JetCorrectorAK8chs->setJetE  ( puppi_e );
    // JetCorrectorAK8chs->setJetA  ( ijet.jetArea() );
    // JetCorrectorAK8chs->setRho   ( rho );
    // JetCorrectorAK8chs->setNPV   ( nvtx );
    // // getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction.
    // vector<float> factors = JetCorrectorAK8chs->getSubCorrections();
    // float corr_factor_L1      = 1.0;
    // float corr_factor_L12     = 1.0;
    // float corr_factor_L123    = 1.0;
    // float corr_factor_L123res = 1.0;
    // if (factors.size() > 0) corr_factor_L1       = factors[0];
    // if (factors.size() > 1) corr_factor_L12      = factors[1];
    // if (factors.size() > 2) corr_factor_L123     = factors[2];
    // if (factors.size() > 3) corr_factor_L123res  = factors[3];
    // double corr_factor_L2 = corr_factor_L12/corr_factor_L1;
    // double corr_factor_L3 = corr_factor_L123/corr_factor_L12;
    // double corr_factor_res = corr_factor_L123res/corr_factor_L123;
    // //double corr_factor_L23 = corr_factor_L2*corr_factor_L3;
    // double corr_factor_L23res = corr_factor_L2*corr_factor_L3*corr_factor_res;


    if(count_AK8CHS==0) PUPPIjet0_P4corr.SetPtEtaPhiM( puppi_pt, puppi_eta, puppi_phi, puppi_mass );
    if(count_AK8CHS==1) PUPPIjet1_P4corr.SetPtEtaPhiM( puppi_pt, puppi_eta, puppi_phi, puppi_mass );



    //------------------------------------
    // SoftDrop subjets
    //------------------------------------
    TLorentzVector sub0_P4_uncorr           ;
    TLorentzVector sub0_P4_L23res           ;
    TLorentzVector sub0_P4_L23resCorrUp     ;
    TLorentzVector sub0_P4_L23resCorrDn     ;
    TLorentzVector sub0_P4_L23resPtSmear    ;
    TLorentzVector sub0_P4_L23resPtSmearUp  ;
    TLorentzVector sub0_P4_L23resPtSmearDn  ;
    TLorentzVector sub0_P4_L123res          ;
    TLorentzVector sub0_P4_L123resCorrUp    ;
    TLorentzVector sub0_P4_L123resCorrDn    ;

    TLorentzVector sub1_P4_uncorr           ;
    TLorentzVector sub1_P4_L23res           ;
    TLorentzVector sub1_P4_L23resCorrUp     ;
    TLorentzVector sub1_P4_L23resCorrDn     ;
    TLorentzVector sub1_P4_L23resPtSmear    ;
    TLorentzVector sub1_P4_L23resPtSmearUp  ;
    TLorentzVector sub1_P4_L23resPtSmearDn  ;
    TLorentzVector sub1_P4_L123res          ;
    TLorentzVector sub1_P4_L123resCorrUp    ;
    TLorentzVector sub1_P4_L123resCorrDn    ;

    double sub0_area  = 0;
    double sub0_tau1  = 0;
    double sub0_tau2  = 0;
    double sub0_tau3  = 0;
    double sub0_flav_hadron  = 0;
    double sub0_flav_parton  = 0;
    double sub0_bdisc = 0;
    double sub1_area  = 0;
    double sub1_tau1  = 0;
    double sub1_tau2  = 0;
    double sub1_tau3  = 0;
    double sub1_flav_hadron  = 0;
    double sub1_flav_parton  = 0;
    double sub1_bdisc = 0;
    double mostMassiveSDsubjetMass = 0;
    int count_SD =0;

    if (!useToolbox_){
      auto const & sdSubjets = ijet.subjets("SoftDrop");
      for ( auto const & it : sdSubjets ) {
        double subjetPt       = it->correctedP4(0).pt();
        double subjetEta      = it->correctedP4(0).eta();
        double subjetPhi      = it->correctedP4(0).phi();
        double subjetMass     = it->correctedP4(0).mass();
        double subjetBdisc    = it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);

        //------------------------------------
        // subjet JEC
        //------------------------------------
        reco::Candidate::LorentzVector uncorrSubjet = it->correctedP4(0);
        JetCorrectorAK4chs -> setJetEta( uncorrSubjet.eta()    );
        JetCorrectorAK4chs -> setJetPt ( uncorrSubjet.pt()     );
        JetCorrectorAK4chs -> setJetE  ( uncorrSubjet.energy() );
        JetCorrectorAK4chs -> setJetA  ( it->jetArea()         );
        JetCorrectorAK4chs -> setRho   ( rho                   );
        JetCorrectorAK4chs -> setNPV   ( nvtx                  );
        double subjet_corr_factor_L123res_full = JetCorrectorAK4chs->getCorrection();
        reco::Candidate::LorentzVector corrSubjetL123res = subjet_corr_factor_L123res_full * uncorrSubjet;

        //------------------------------------
        // subjet L23 JEC
        //------------------------------------
        JetCorrectorAK4chs->setJetEta( uncorrSubjet.eta()    );
        JetCorrectorAK4chs->setJetPt ( uncorrSubjet.pt()     );
        JetCorrectorAK4chs->setJetE  ( uncorrSubjet.energy() );
        JetCorrectorAK4chs->setJetA  ( it->jetArea()         );
        JetCorrectorAK4chs->setRho   ( rho                   );
        JetCorrectorAK4chs->setNPV   ( nvtx                  );
        // getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction.
        vector<float> subjet_factors = JetCorrectorAK4chs->getSubCorrections();
        float subjet_corr_factor_L1      = 1.0;
        float subjet_corr_factor_L12     = 1.0;
        float subjet_corr_factor_L123    = 1.0;
        float subjet_corr_factor_L123res = 1.0;
        if (factors.size() > 0) subjet_corr_factor_L1      = subjet_factors[0];
        if (factors.size() > 1) subjet_corr_factor_L12     = subjet_factors[1];
        if (factors.size() > 2) subjet_corr_factor_L123    = subjet_factors[2];
        if (factors.size() > 3) subjet_corr_factor_L123res = subjet_factors[3];
        double subjet_corr_factor_L2     = subjet_corr_factor_L12     / subjet_corr_factor_L1     ;
        double subjet_corr_factor_L3     = subjet_corr_factor_L123    / subjet_corr_factor_L12    ;
        double subjet_corr_factor_res    = subjet_corr_factor_L123res / subjet_corr_factor_L123   ;
        double subjet_corr_factor_L23    = subjet_corr_factor_L2 * subjet_corr_factor_L3     ;
        double subjet_corr_factor_L23res = subjet_corr_factor_L2 * subjet_corr_factor_L3 * subjet_corr_factor_res    ;
        if (verbose_) cout<<"subjet corr: L1 "<<subjet_corr_factor_L1<<" L23 "<<subjet_corr_factor_L23<<" L23res "<<subjet_corr_factor_L23res<<" L123res"<<subjet_corr_factor_L123res<<endl;
        reco::Candidate::LorentzVector corrSubjetL23res   = subjet_corr_factor_L23res * uncorrSubjet;

        //------------------------------------
        // subjet JEC uncertainty
        //------------------------------------
        double subjet_corrDn_L23 =   1.0;
        double subjet_corrDn_L123 = 1.0;
        JetCorrUncertAK4chs->setJetPhi(  corrSubjetL123res.phi()  );
        JetCorrUncertAK4chs->setJetEta(  corrSubjetL123res.eta()  );
        JetCorrUncertAK4chs->setJetPt(   corrSubjetL123res.pt()   );
        double corrDn_temp2 = JetCorrUncertAK4chs->getUncertainty(0);
        subjet_corrDn_L23   = subjet_corr_factor_L23res - corrDn_temp2;
        subjet_corrDn_L123  = subjet_corr_factor_L123res_full - corrDn_temp2;

        double subjet_corrUp_L23   = 1.0;
        double subjet_corrUp_L123 = 1.0;
        JetCorrUncertAK4chs->setJetPhi(  corrSubjetL123res.phi()  );
        JetCorrUncertAK4chs->setJetEta(  corrSubjetL123res.eta()  );
        JetCorrUncertAK4chs->setJetPt(   corrSubjetL123res.pt()   );
        double corrUp_temp2 = JetCorrUncertAK4chs->getUncertainty(1);
        subjet_corrUp_L23   = subjet_corr_factor_L23res + corrUp_temp2;
        subjet_corrUp_L123  = subjet_corr_factor_L123res_full + corrUp_temp2;

        reco::Candidate::LorentzVector corrSubjetL123resCorrDn  = subjet_corrDn_L123  * uncorrSubjet;
        reco::Candidate::LorentzVector corrSubjetL123resCorrUp  = subjet_corrUp_L123  * uncorrSubjet;
        reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
        reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;


        //------------------------------------
        // subjet values for Tree
        //------------------------------------
        if (count_SD==0){
          sub0_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
          sub0_P4_L123res           .SetPtEtaPhiM( corrSubjetL123res.pt()   , corrSubjetL123res.eta()   , corrSubjetL123res.phi()   , corrSubjetL123res.mass()    );
          sub0_P4_L23res            .SetPtEtaPhiM( corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass()     );
          sub0_P4_L123resCorrUp    .SetPtEtaPhiM( corrSubjetL123resCorrUp.pt() , corrSubjetL123resCorrUp.eta() , corrSubjetL123resCorrUp.phi() , corrSubjetL123resCorrUp.mass()  );
          sub0_P4_L23resCorrUp     .SetPtEtaPhiM( corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass()   );
          sub0_P4_L123resCorrDn    .SetPtEtaPhiM( corrSubjetL123resCorrDn.pt() , corrSubjetL123resCorrDn.eta() , corrSubjetL123resCorrDn.phi() , corrSubjetL123resCorrUp.mass()  );
          sub0_P4_L23resCorrDn     .SetPtEtaPhiM( corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass()     );
          sub0_area   = it->jetArea() ;
          sub0_flav_parton   = it->partonFlavour();
          sub0_flav_hadron   = it->hadronFlavour();
          sub0_bdisc  = subjetBdisc;
        }
        if (count_SD==1) {
          sub1_P4_uncorr          .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
          sub1_P4_L123res         .SetPtEtaPhiM( corrSubjetL123res.pt()   , corrSubjetL123res.eta()   , corrSubjetL123res.phi()   , corrSubjetL123res.mass()    );
          sub1_P4_L23res          .SetPtEtaPhiM( corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass()     );
          sub1_P4_L123resCorrUp  .SetPtEtaPhiM( corrSubjetL123resCorrUp.pt() , corrSubjetL123resCorrUp.eta() , corrSubjetL123resCorrUp.phi() , corrSubjetL123resCorrUp.mass()  );
          sub1_P4_L23resCorrUp   .SetPtEtaPhiM( corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass()   );
          sub1_P4_L123resCorrDn  .SetPtEtaPhiM( corrSubjetL123resCorrDn.pt() , corrSubjetL123resCorrDn.eta() , corrSubjetL123resCorrDn.phi() , corrSubjetL123resCorrUp.mass()  );
          sub1_P4_L23resCorrDn   .SetPtEtaPhiM( corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass()     );
          sub1_area   = it->jetArea() ;
          sub1_flav_parton   = it->partonFlavour();
          sub1_flav_hadron   = it->hadronFlavour();
          sub1_bdisc  = subjetBdisc;
        }
        if (subjetMass > mostMassiveSDsubjetMass) mostMassiveSDsubjetMass = subjetMass;

        if (verbose_) cout<<" SD Subjet pt "<<subjetPt<<" Eta "<<subjetEta<<" deltaRsubjetJet "<<deltaRsubjetJet<<" Mass "<<subjetMass<<" Bdisc "<<subjetBdisc<<endl;
        count_SD++;
      }
    }
    if (useToolbox_){
      for (const pat::Jet &isub : *AK8CHSsub) {

        double subjetPt       = isub.correctedP4(0).pt();
        double subjetEta      = isub.correctedP4(0).eta();
        double subjetPhi      = isub.correctedP4(0).phi();
        double subjetMass     = isub.correctedP4(0).mass();
        double subjetBdisc    = isub.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

        double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);
        if (deltaRsubjetJet<0.8){
          if (verbose_) cout<<" matched subjet with mass "<< subjetMass<<endl;

          //------------------------------------
          // subjet JEC
          //------------------------------------
          reco::Candidate::LorentzVector uncorrSubjet = isub.correctedP4(0);
          JetCorrectorAK4chs -> setJetEta( uncorrSubjet.eta()    );
          JetCorrectorAK4chs -> setJetPt ( uncorrSubjet.pt()     );
          JetCorrectorAK4chs -> setJetE  ( uncorrSubjet.energy() );
          JetCorrectorAK4chs -> setJetA  ( isub.jetArea()        );
          JetCorrectorAK4chs -> setRho   ( rho                   );
          JetCorrectorAK4chs -> setNPV   ( nvtx                  );
          double subjet_corr_factor_L123res_full = JetCorrectorAK4chs->getCorrection();
          reco::Candidate::LorentzVector corrSubjetL123res = subjet_corr_factor_L123res_full * uncorrSubjet;

          //------------------------------------
          // subjet L23 JEC
          //------------------------------------
          JetCorrectorAK4chs->setJetEta( uncorrSubjet.eta()    );
          JetCorrectorAK4chs->setJetPt ( uncorrSubjet.pt()     );
          JetCorrectorAK4chs->setJetE  ( uncorrSubjet.energy() );
          JetCorrectorAK4chs->setJetA  ( isub.jetArea()         );
          JetCorrectorAK4chs->setRho   ( rho                   );
          JetCorrectorAK4chs->setNPV   ( nvtx                  );
          // getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction.
          vector<float> subjet_factors = JetCorrectorAK4chs->getSubCorrections();
          float subjet_corr_factor_L1      = 1.0;
          float subjet_corr_factor_L12     = 1.0;
          float subjet_corr_factor_L123    = 1.0;
          float subjet_corr_factor_L123res = 1.0;
          if (factors.size() > 0) subjet_corr_factor_L1      = subjet_factors[0];
          if (factors.size() > 1) subjet_corr_factor_L12     = subjet_factors[1];
          if (factors.size() > 2) subjet_corr_factor_L123    = subjet_factors[2];
          if (factors.size() > 3) subjet_corr_factor_L123res = subjet_factors[3];
          double subjet_corr_factor_L2     = subjet_corr_factor_L12     / subjet_corr_factor_L1     ;
          double subjet_corr_factor_L3     = subjet_corr_factor_L123    / subjet_corr_factor_L12    ;
          double subjet_corr_factor_res    = subjet_corr_factor_L123res / subjet_corr_factor_L123   ;
          double subjet_corr_factor_L23    = subjet_corr_factor_L2 * subjet_corr_factor_L3     ;
          double subjet_corr_factor_L23res = subjet_corr_factor_L2 * subjet_corr_factor_L3 * subjet_corr_factor_res    ;
          if (verbose_) cout<<"  subjet corr: L1 "<<subjet_corr_factor_L1<<" L23 "<<subjet_corr_factor_L23<<" L23res "<<subjet_corr_factor_L23res<<" L123res "<<subjet_corr_factor_L123res<<endl;
          reco::Candidate::LorentzVector corrSubjetL23res   = subjet_corr_factor_L23res * uncorrSubjet;

   cout<<"CHSsubjet area "<<isub.jetArea() <<endl;

          //------------------------------------
          // subjet JEC uncertainty
          //------------------------------------
          double subjet_corrDn_L23 =   1.0;
          double subjet_corrDn_L123 = 1.0;
          JetCorrUncertAK4chs->setJetPhi(  corrSubjetL123res.phi()  );
          JetCorrUncertAK4chs->setJetEta(  corrSubjetL123res.eta()  );
          JetCorrUncertAK4chs->setJetPt(   corrSubjetL123res.pt()   );
          double corrDn_temp2 = JetCorrUncertAK4chs->getUncertainty(0);
          subjet_corrDn_L23   = subjet_corr_factor_L23res - corrDn_temp2;
          subjet_corrDn_L123  = subjet_corr_factor_L123res_full - corrDn_temp2;

          double subjet_corrUp_L23   = 1.0;
          double subjet_corrUp_L123 = 1.0;
          JetCorrUncertAK4chs->setJetPhi(  corrSubjetL123res.phi()  );
          JetCorrUncertAK4chs->setJetEta(  corrSubjetL123res.eta()  );
          JetCorrUncertAK4chs->setJetPt(   corrSubjetL123res.pt()   );
          double corrUp_temp2 = JetCorrUncertAK4chs->getUncertainty(1);
          subjet_corrUp_L23   = subjet_corr_factor_L23res + corrUp_temp2;
          subjet_corrUp_L123  = subjet_corr_factor_L123res_full + corrUp_temp2;

          reco::Candidate::LorentzVector corrSubjetL123resCorrDn  = subjet_corrDn_L123  * uncorrSubjet;
          reco::Candidate::LorentzVector corrSubjetL123resCorrUp  = subjet_corrUp_L123  * uncorrSubjet;
          reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
          reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;

          //------------------------------------
          // subjet JER SF
          //------------------------------------
          TLorentzVector GenSubJet;
          double ptsmear   = 1;
          double ptsmearUp = 1;
          double ptsmearDn = 1;
          if (!iEvent.isRealData()){
            const reco::GenJet* genJet = isub.genJet();
            if (genJet) {
              GenSubJet.SetPtEtaPhiM( genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass() );
              if (verbose_) cout<<"  SD subjet genJet pt "<<genJet->pt()<<" mass "<<genJet->mass()<<" reco pt "<<subjetPt<<" reco mass "<<subjetMass<<endl;
            }
            double jer_sf    = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}});
            double jer_sf_up = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::UP);
            double jer_sf_dn = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::DOWN);
            if (verbose_) std::cout << " SD subjet JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn << std::endl;
            double recopt    = corrSubjetL23res.pt();
            double genpt     = GenJetMatched.Perp();
            double deltapt   = (recopt-genpt)*(jer_sf-1.0);
            double deltaptUp = (recopt-genpt)*(jer_sf_up-1.0);
            double deltaptDn = (recopt-genpt)*(jer_sf_dn-1.0);
            ptsmear   = std::max((double)0.0, (recopt+deltapt)/recopt     );
            ptsmearUp = std::max((double)0.0, (recopt+deltaptUp)/recopt   );
            ptsmearDn = std::max((double)0.0, (recopt+deltaptDn)/recopt   );
            if (verbose_) std::cout<<" SD subjet ptsmear "<<ptsmear<<" ptsmearUp "<<ptsmearUp<<" ptsmearDn "<<ptsmearDn<<endl;
          }
          reco::Candidate::LorentzVector corrSubjetL23resPtSmear   = ptsmear * corrSubjetL23res ;
          reco::Candidate::LorentzVector corrSubjetL23resPtSmearUp = ptsmearUp * corrSubjetL23res ;
          reco::Candidate::LorentzVector corrSubjetL23resPtSmearDn = ptsmearDn * corrSubjetL23res ;

          //------------------------------------
          // subjet values for Tree
          //------------------------------------
          if (count_SD==0){
            sub0_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
            sub0_P4_L23res            .SetPtEtaPhiM( corrSubjetL23res          .pt() , corrSubjetL23res          .eta()  , corrSubjetL23res          .phi()  , corrSubjetL23res          .mass()  );
            sub0_P4_L23resCorrUp      .SetPtEtaPhiM( corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta()  , corrSubjetL23resCorrUp    .phi()  , corrSubjetL23resCorrUp    .mass()  );
            sub0_P4_L23resCorrDn      .SetPtEtaPhiM( corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta()  , corrSubjetL23resCorrDn    .phi()  , corrSubjetL23resCorrDn    .mass()  );
            sub0_P4_L23resPtSmear     .SetPtEtaPhiM( corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta()  , corrSubjetL23resPtSmear   .phi()  , corrSubjetL23resPtSmear   .mass()  );
            sub0_P4_L23resPtSmearUp   .SetPtEtaPhiM( corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta()  , corrSubjetL23resPtSmearUp .phi()  , corrSubjetL23resPtSmearUp .mass()  );
            sub0_P4_L23resPtSmearDn   .SetPtEtaPhiM( corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta()  , corrSubjetL23resPtSmearDn .phi()  , corrSubjetL23resPtSmearDn .mass()  );
            sub0_P4_L123res           .SetPtEtaPhiM( corrSubjetL123res         .pt() , corrSubjetL123res         .eta()  , corrSubjetL123res         .phi()  , corrSubjetL123res         .mass()  );
            sub0_P4_L123resCorrUp     .SetPtEtaPhiM( corrSubjetL123resCorrUp   .pt() , corrSubjetL123resCorrUp   .eta()  , corrSubjetL123resCorrUp   .phi()  , corrSubjetL123resCorrUp   .mass()  );
            sub0_P4_L123resCorrDn     .SetPtEtaPhiM( corrSubjetL123resCorrDn   .pt() , corrSubjetL123resCorrDn   .eta()  , corrSubjetL123resCorrDn   .phi()  , corrSubjetL123resCorrDn   .mass()  );

            sub0_area          = isub.jetArea() ;
            sub0_flav_parton   = isub.partonFlavour();
            sub0_flav_hadron   = isub.hadronFlavour();
            sub0_bdisc         = subjetBdisc;
            // available from toolbox only (80X)
            sub0_tau1          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau1");
            sub0_tau2          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau2");
            sub0_tau3          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau3");
          }
          if (count_SD==1) {
            sub1_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
            sub1_P4_L23res            .SetPtEtaPhiM( corrSubjetL23res          .pt() , corrSubjetL23res          .eta()  , corrSubjetL23res          .phi()  , corrSubjetL23res          .mass()  );
            sub1_P4_L23resCorrUp      .SetPtEtaPhiM( corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta()  , corrSubjetL23resCorrUp    .phi()  , corrSubjetL23resCorrUp    .mass()  );
            sub1_P4_L23resCorrDn      .SetPtEtaPhiM( corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta()  , corrSubjetL23resCorrDn    .phi()  , corrSubjetL23resCorrDn    .mass()  );
            sub1_P4_L23resPtSmear     .SetPtEtaPhiM( corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta()  , corrSubjetL23resPtSmear   .phi()  , corrSubjetL23resPtSmear   .mass()  );
            sub1_P4_L23resPtSmearUp   .SetPtEtaPhiM( corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta()  , corrSubjetL23resPtSmearUp .phi()  , corrSubjetL23resPtSmearUp .mass()  );
            sub1_P4_L23resPtSmearDn   .SetPtEtaPhiM( corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta()  , corrSubjetL23resPtSmearDn .phi()  , corrSubjetL23resPtSmearDn .mass()  );
            sub1_P4_L123res           .SetPtEtaPhiM( corrSubjetL123res         .pt() , corrSubjetL123res         .eta()  , corrSubjetL123res         .phi()  , corrSubjetL123res         .mass()  );
            sub1_P4_L123resCorrUp     .SetPtEtaPhiM( corrSubjetL123resCorrUp   .pt() , corrSubjetL123resCorrUp   .eta()  , corrSubjetL123resCorrUp   .phi()  , corrSubjetL123resCorrUp   .mass()  );
            sub1_P4_L123resCorrDn     .SetPtEtaPhiM( corrSubjetL123resCorrDn   .pt() , corrSubjetL123resCorrDn   .eta()  , corrSubjetL123resCorrDn   .phi()  , corrSubjetL123resCorrDn   .mass()  );

            sub1_area          = isub.jetArea() ;
            sub1_flav_parton   = isub.partonFlavour();
            sub1_flav_hadron   = isub.hadronFlavour();
            sub1_bdisc         = subjetBdisc;
            // available from toolbox only (80X)
            sub1_tau1          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau1");
            sub1_tau2          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau2");
            sub1_tau3          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau3");
          }
          if (subjetMass > mostMassiveSDsubjetMass) mostMassiveSDsubjetMass = subjetMass;

          if (verbose_) cout<<"  SD Subjet pt "<<subjetPt<<" Eta "<<subjetEta<<" deltaRsubjetJet "<<deltaRsubjetJet<<" Mass "<<subjetMass<<" corrMass "<<corrSubjetL23res.mass() <<" Bdisc "<<subjetBdisc<<endl;
          if (verbose_) cout<<"    sub0_tau1 "<<sub0_tau1<<" sub0_tau2 "<<sub0_tau2<<" sub0_tau3 "<<sub0_tau3<<endl;
          count_SD++;

        }
      }
    }

    TLorentzVector sumSDsubjets_P4_uncorr           ;
    TLorentzVector sumSDsubjets_P4_L23res           ;
    TLorentzVector sumSDsubjets_P4_L23resCorrUp     ;
    TLorentzVector sumSDsubjets_P4_L23resCorrDn     ;
    TLorentzVector sumSDsubjets_P4_L23resPtSmear    ;
    TLorentzVector sumSDsubjets_P4_L23resPtSmearUp  ;
    TLorentzVector sumSDsubjets_P4_L23resPtSmearDn  ;
    TLorentzVector sumSDsubjets_P4_L123res          ;
    TLorentzVector sumSDsubjets_P4_L123resCorrDn    ;
    TLorentzVector sumSDsubjets_P4_L123resCorrUp    ;

    if (count_SD>1){
      sumSDsubjets_P4_uncorr             = sub0_P4_uncorr              + sub1_P4_uncorr            ;
      sumSDsubjets_P4_L23res             = sub0_P4_L23res              + sub1_P4_L23res            ;
      sumSDsubjets_P4_L23resCorrUp       = sub0_P4_L23resCorrUp        + sub1_P4_L23resCorrUp      ;
      sumSDsubjets_P4_L23resCorrDn       = sub0_P4_L23resCorrDn        + sub1_P4_L23resCorrDn      ;
      sumSDsubjets_P4_L23resPtSmear      = sub0_P4_L23resPtSmear       + sub1_P4_L23resPtSmear     ;
      sumSDsubjets_P4_L23resPtSmearUp    = sub0_P4_L23resPtSmearUp     + sub1_P4_L23resPtSmearUp   ;
      sumSDsubjets_P4_L23resPtSmearDn    = sub0_P4_L23resPtSmearDn     + sub1_P4_L23resPtSmearDn   ;
      sumSDsubjets_P4_L123res            = sub0_P4_L123res             + sub1_P4_L123res           ;
      sumSDsubjets_P4_L123resCorrUp      = sub0_P4_L123resCorrUp       + sub1_P4_L123resCorrUp     ;
      sumSDsubjets_P4_L123resCorrDn      = sub0_P4_L123resCorrDn       + sub1_P4_L123resCorrDn     ;
    }

    double maxbdisc = 0 ;
    double maxbdiscflav_hadron = 0 ;
    double maxbdiscflav_parton = 0 ;
    if (sub0_bdisc>=sub1_bdisc){
      maxbdisc = sub0_bdisc;
      maxbdiscflav_hadron = sub0_flav_hadron;
      maxbdiscflav_parton = sub0_flav_parton;
    }
    else if (sub1_bdisc>sub0_bdisc){
      maxbdisc = sub1_bdisc;
      maxbdiscflav_hadron = sub1_flav_hadron;
      maxbdiscflav_parton = sub1_flav_parton;
    }

    //------------------------------------
    // PUPPI SoftDrop subjets
    //------------------------------------
    TLorentzVector pup0_P4_uncorr           ;
    TLorentzVector pup0_P4_L23res           ;
    TLorentzVector pup0_P4_L23resCorrUp     ;
    TLorentzVector pup0_P4_L23resCorrDn     ;
    TLorentzVector pup0_P4_L23resPtSmear    ;
    TLorentzVector pup0_P4_L23resPtSmearUp  ;
    TLorentzVector pup0_P4_L23resPtSmearDn  ;

    TLorentzVector pup1_P4_uncorr           ;
    TLorentzVector pup1_P4_L23res           ;
    TLorentzVector pup1_P4_L23resCorrUp     ;
    TLorentzVector pup1_P4_L23resCorrDn     ;
    TLorentzVector pup1_P4_L23resPtSmear    ;
    TLorentzVector pup1_P4_L23resPtSmearUp  ;
    TLorentzVector pup1_P4_L23resPtSmearDn  ;


    double pup0_area  = 0;
    double pup0_tau1  = 0;
    double pup0_tau2  = 0;
    double pup0_tau3  = 0;
    double pup0_flav_hadron  = 0;
    double pup0_flav_parton  = 0;
    double pup0_bdisc = 0;
    double pup1_area  = 0;
    double pup1_tau1  = 0;
    double pup1_tau2  = 0;
    double pup1_tau3  = 0;
    double pup1_flav_hadron  = 0;
    double pup1_flav_parton  = 0;
    double pup1_bdisc = 0;
    double mostMassiveSDPUPPIsubjetMass = 0;
    int count_pup=0;

    if (!useToolbox_){
      auto const & sdSubjetsPuppi = ijet.subjets("SoftDropPuppi");
      for ( auto const & it : sdSubjetsPuppi ) {
        double subjetPt       = it->correctedP4(0).pt();
        double subjetEta      = it->correctedP4(0).eta();
        double subjetPhi      = it->correctedP4(0).phi();
        double subjetMass     = it->correctedP4(0).mass();
        double subjetBdisc    = it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);
        if (verbose_) cout<<" SD Subjet pt "<<subjetPt<<" Eta "<<subjetEta<<" deltaRsubjetJet "<<deltaRsubjetJet<<" Mass "<<subjetMass<<" Bdisc "<<subjetBdisc<<endl;

        //------------------------------------
        // PUPPI subjet JEC
        //------------------------------------
        reco::Candidate::LorentzVector uncorrSubjet = it->correctedP4(0);
        JetCorrectorAK4pup -> setJetEta( uncorrSubjet.eta()    );
        JetCorrectorAK4pup -> setJetPt ( uncorrSubjet.pt()     );
        JetCorrectorAK4pup -> setJetE  ( uncorrSubjet.energy() );
        JetCorrectorAK4pup -> setJetA  ( it->jetArea()         );
        JetCorrectorAK4pup -> setRho   ( rho                   );
        JetCorrectorAK4pup -> setNPV   ( nvtx                  );
        double subjet_corr_factor_L23res_full = JetCorrectorAK4pup->getCorrection();
        reco::Candidate::LorentzVector corrSubjetL23res = subjet_corr_factor_L23res_full * uncorrSubjet;

        //------------------------------------
        // PUPPI subjet JEC uncertainty
        //------------------------------------
        double subjet_corrDn_L23 =   1.0;
        JetCorrUncertAK4pup->setJetPhi(  corrSubjetL23res.phi()  );
        JetCorrUncertAK4pup->setJetEta(  corrSubjetL23res.eta()  );
        JetCorrUncertAK4pup->setJetPt(   corrSubjetL23res.pt()   );
        subjet_corrDn_L23   = subjet_corr_factor_L23res_full - JetCorrUncertAK4pup->getUncertainty(0);
        double subjet_corrUp_L23 =   1.0;
        JetCorrUncertAK4pup->setJetPhi(  corrSubjetL23res.phi()  );
        JetCorrUncertAK4pup->setJetEta(  corrSubjetL23res.eta()  );
        JetCorrUncertAK4pup->setJetPt(   corrSubjetL23res.pt()   );
        subjet_corrUp_L23   = subjet_corr_factor_L23res_full + JetCorrUncertAK4pup->getUncertainty(1);

        reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
        reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;

        //------------------------------------
        // subjet values for Tree
        //------------------------------------

        if (count_pup==0){
          pup0_P4_uncorr           .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
          pup0_P4_L23res           .SetPtEtaPhiM( corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass()     );
          pup0_P4_L23resCorrUp    .SetPtEtaPhiM( corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass()   );
          pup0_P4_L23resCorrDn    .SetPtEtaPhiM( corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass()     );
          pup0_area   = it->jetArea() ;
          if (useToolbox_){
            pup0_tau1   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
            pup0_tau2   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
            pup0_tau3   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");
          }
          pup0_flav_parton   = it->partonFlavour();
          pup0_flav_hadron   = it->hadronFlavour();
          pup0_bdisc  = subjetBdisc;
        }
        if (count_pup==1) {
          pup1_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
          pup1_P4_L23res           .SetPtEtaPhiM( corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass()     );
          pup1_P4_L23resCorrUp    .SetPtEtaPhiM( corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass()   );
          pup1_P4_L23resCorrDn    .SetPtEtaPhiM( corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass()     );
          pup1_area   = it->jetArea() ;
          if (useToolbox_){
            pup1_tau1   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
            pup1_tau2   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
            pup1_tau3   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");
          }
          pup1_flav_parton   = it->partonFlavour();
          pup1_flav_hadron   = it->hadronFlavour();
          pup1_bdisc  = subjetBdisc;
        }

        if (subjetMass > mostMassiveSDPUPPIsubjetMass) mostMassiveSDPUPPIsubjetMass = subjetMass;
        count_pup++;
      }
    }
    if (useToolbox_){
      for (const pat::Jet &isub : *AK8PUPPIsub) {

        double subjetPt        = isub.correctedP4(0).pt();
        double subjetEta       = isub.correctedP4(0).eta();
        double subjetPhi       = isub.correctedP4(0).phi();
        double subjetMass      = isub.correctedP4(0).mass();
        double subjetBdisc     = isub.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);

        if (deltaRsubjetJet<0.8){
          if (verbose_) cout<<" matched puppi subjet with mass "<< subjetMass<<endl;

          //------------------------------------
          // PUPPI subjet JEC
          //------------------------------------
          reco::Candidate::LorentzVector uncorrSubjet = isub.correctedP4(0);
          JetCorrectorAK4pup -> setJetEta( uncorrSubjet.eta()    );
          JetCorrectorAK4pup -> setJetPt ( uncorrSubjet.pt()     );
          JetCorrectorAK4pup -> setJetE  ( uncorrSubjet.energy() );
          JetCorrectorAK4pup -> setJetA  ( isub.jetArea()         );
          JetCorrectorAK4pup -> setRho   ( rho                   );
          JetCorrectorAK4pup -> setNPV   ( nvtx                  );
          double subjet_corr_factor_L23res_full = JetCorrectorAK4pup->getCorrection();
          reco::Candidate::LorentzVector corrSubjetL23res = subjet_corr_factor_L23res_full * uncorrSubjet;

          //------------------------------------
          // PUPPI subjet JEC uncertainty
          //------------------------------------
          double subjet_corrDn_L23 =   1.0;
          JetCorrUncertAK4pup->setJetPhi(  corrSubjetL23res.phi()  );
          JetCorrUncertAK4pup->setJetEta(  corrSubjetL23res.eta()  );
          JetCorrUncertAK4pup->setJetPt(   corrSubjetL23res.pt()   );
          subjet_corrDn_L23   = subjet_corr_factor_L23res_full - JetCorrUncertAK4pup->getUncertainty(0);
          double subjet_corrUp_L23 =   1.0;
          JetCorrUncertAK4pup->setJetPhi(  corrSubjetL23res.phi()  );
          JetCorrUncertAK4pup->setJetEta(  corrSubjetL23res.eta()  );
          JetCorrUncertAK4pup->setJetPt(   corrSubjetL23res.pt()   );
          subjet_corrUp_L23   = subjet_corr_factor_L23res_full + JetCorrUncertAK4pup->getUncertainty(1);

          reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
          reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;

          //------------------------------------
          // subjet JER SF
          //------------------------------------
          TLorentzVector GenSubJet;
          double ptsmear   = 1;
          double ptsmearUp = 1;
          double ptsmearDn = 1;
          if (!iEvent.isRealData()){
            const reco::GenJet* genJet = isub.genJet();
            if (genJet) {
              GenSubJet.SetPtEtaPhiM( genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass() );
              if (verbose_) cout<<"  SD subjet genJet pt "<<genJet->pt()<<" mass "<<genJet->mass()<<" reco pt "<<subjetPt<<" reco mass "<<subjetMass<<endl;
            }
            double jer_sf    = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}});
            double jer_sf_up = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::UP);
            double jer_sf_dn = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::DOWN);
            if (verbose_) std::cout << " SD subjet JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn << std::endl;
            double recopt    = corrSubjetL23res.pt();
            double genpt     = GenJetMatched.Perp();
            double deltapt   = (recopt-genpt)*(jer_sf-1.0);
            double deltaptUp = (recopt-genpt)*(jer_sf_up-1.0);
            double deltaptDn = (recopt-genpt)*(jer_sf_dn-1.0);
            ptsmear   = std::max((double)0.0, (recopt+deltapt)/recopt     );
            ptsmearUp = std::max((double)0.0, (recopt+deltaptUp)/recopt   );
            ptsmearDn = std::max((double)0.0, (recopt+deltaptDn)/recopt   );
            if (verbose_) std::cout<<" SD subjet ptsmear "<<ptsmear<<" ptsmearUp "<<ptsmearUp<<" ptsmearDn "<<ptsmearDn<<endl;
          }
          reco::Candidate::LorentzVector corrSubjetL23resPtSmear   = ptsmear   * corrSubjetL23res ;
          reco::Candidate::LorentzVector corrSubjetL23resPtSmearUp = ptsmearUp * corrSubjetL23res ;
          reco::Candidate::LorentzVector corrSubjetL23resPtSmearDn = ptsmearDn * corrSubjetL23res ;

          //------------------------------------
          // subjet values for Tree
          //------------------------------------

          if (count_pup==0){
            pup0_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
            pup0_P4_L23res            .SetPtEtaPhiM( corrSubjetL23res          .pt() , corrSubjetL23res          .eta() , corrSubjetL23res          .phi() , corrSubjetL23res          .mass() );
            pup0_P4_L23resCorrUp      .SetPtEtaPhiM( corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta() , corrSubjetL23resCorrUp    .phi() , corrSubjetL23resCorrUp    .mass() );
            pup0_P4_L23resCorrDn      .SetPtEtaPhiM( corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta() , corrSubjetL23resCorrDn    .phi() , corrSubjetL23resCorrDn    .mass() );
            pup0_P4_L23resPtSmear     .SetPtEtaPhiM( corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta() , corrSubjetL23resPtSmear   .phi() , corrSubjetL23resPtSmear   .mass() );
            pup0_P4_L23resPtSmearUp   .SetPtEtaPhiM( corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta() , corrSubjetL23resPtSmearUp .phi() , corrSubjetL23resPtSmearUp .mass() );
            pup0_P4_L23resPtSmearDn   .SetPtEtaPhiM( corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta() , corrSubjetL23resPtSmearDn .phi() , corrSubjetL23resPtSmearDn .mass() );

            pup0_tau1   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau1");
            pup0_tau2   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau2");
            pup0_tau3   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau3");

            pup0_flav_parton   = isub.partonFlavour();
            pup0_flav_hadron   = isub.hadronFlavour();
            pup0_area          = isub.jetArea() ;
            pup0_bdisc         = subjetBdisc;
          }
          if (count_pup==1) {
            pup1_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
            pup1_P4_L23res            .SetPtEtaPhiM( corrSubjetL23res          .pt() , corrSubjetL23res          .eta() , corrSubjetL23res          .phi() , corrSubjetL23res          .mass() );
            pup1_P4_L23resCorrUp      .SetPtEtaPhiM( corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta() , corrSubjetL23resCorrUp    .phi() , corrSubjetL23resCorrUp    .mass() );
            pup1_P4_L23resCorrDn      .SetPtEtaPhiM( corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta() , corrSubjetL23resCorrDn    .phi() , corrSubjetL23resCorrDn    .mass() );
            pup1_P4_L23resPtSmear     .SetPtEtaPhiM( corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta() , corrSubjetL23resPtSmear   .phi() , corrSubjetL23resPtSmear   .mass() );
            pup1_P4_L23resPtSmearUp   .SetPtEtaPhiM( corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta() , corrSubjetL23resPtSmearUp .phi() , corrSubjetL23resPtSmearUp .mass() );
            pup1_P4_L23resPtSmearDn   .SetPtEtaPhiM( corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta() , corrSubjetL23resPtSmearDn .phi() , corrSubjetL23resPtSmearDn .mass() );

            pup1_tau1   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau1");
            pup1_tau2   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau2");
            pup1_tau3   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau3");

            pup1_flav_parton   = isub.partonFlavour();
            pup1_flav_hadron   = isub.hadronFlavour();
            pup1_area          = isub.jetArea() ;
            pup1_bdisc         = subjetBdisc;
          }

          if (subjetMass > mostMassiveSDPUPPIsubjetMass) mostMassiveSDPUPPIsubjetMass = subjetMass;
          count_pup++;
           if (verbose_) cout<<"  SD Subjet pt "<<subjetPt<<" Eta "<<subjetEta<<" deltaRsubjetJet "<<deltaRsubjetJet<<" Mass "<<subjetMass<<" Bdisc "<<subjetBdisc<<endl;


        }
      }
    }

    TLorentzVector sumPUPsubjets_P4_uncorr          ;
    TLorentzVector sumPUPsubjets_P4_L23res          ;
    TLorentzVector sumPUPsubjets_P4_L23resCorrUp    ;
    TLorentzVector sumPUPsubjets_P4_L23resCorrDn    ;
    TLorentzVector sumPUPsubjets_P4_L23resPtSmear   ;
    TLorentzVector sumPUPsubjets_P4_L23resPtSmearUp ;
    TLorentzVector sumPUPsubjets_P4_L23resPtSmearDn ;
    if (count_SD>1){
      sumPUPsubjets_P4_uncorr            = pup0_P4_uncorr            + pup1_P4_uncorr            ;
      sumPUPsubjets_P4_L23res            = pup0_P4_L23res            + pup1_P4_L23res            ;
      sumPUPsubjets_P4_L23resCorrUp      = pup0_P4_L23resCorrUp      + pup1_P4_L23resCorrUp      ;
      sumPUPsubjets_P4_L23resCorrDn      = pup0_P4_L23resCorrDn      + pup1_P4_L23resCorrDn      ;
      sumPUPsubjets_P4_L23resPtSmear     = pup0_P4_L23resPtSmear     + pup1_P4_L23resPtSmear     ;
      sumPUPsubjets_P4_L23resPtSmearUp   = pup0_P4_L23resPtSmearUp   + pup1_P4_L23resPtSmearUp   ;
      sumPUPsubjets_P4_L23resPtSmearDn   = pup0_P4_L23resPtSmearDn   + pup1_P4_L23resPtSmearDn   ;
    }


    double pup_maxbdisc = 0 ;
    double pup_maxbdiscflav_hadron = 0 ;
    double pup_maxbdiscflav_parton = 0 ;
    if (pup0_bdisc>=pup1_bdisc){
      pup_maxbdisc = pup0_bdisc;
      pup_maxbdiscflav_hadron = pup0_flav_hadron;
      pup_maxbdiscflav_parton = pup0_flav_parton;
    }
    else if (pup1_bdisc>pup0_bdisc){
      pup_maxbdisc = pup1_bdisc;
      pup_maxbdiscflav_hadron = pup1_flav_hadron;
      pup_maxbdiscflav_parton = pup1_flav_parton;
    }

    h_ak8chs_softDropMass    ->Fill( sumSDsubjets_P4_uncorr.M()   );
    h_ak8puppi_softDropMass  ->Fill( sumPUPsubjets_P4_uncorr.M()  );

    //------------------------------------
    // Gen particle info
    //------------------------------------

    double deltaR_jet_p1 =  jet_p4.DeltaR(hardest_parton_hardScatterOutgoing_p4);
    double deltaR_jet_p2 =  jet_p4.DeltaR(second_hardest_parton_hardScatterOutgoing_p4);
    bool jet_matched_p1 = false;
    bool jet_matched_p2 = false;
    if (deltaR_jet_p1<deltaR_jet_p2) jet_matched_p1 = true;
    if (deltaR_jet_p2<deltaR_jet_p1) jet_matched_p2 = true;

    //------------------------------------
    // Fill LeptTree
    //------------------------------------
    if (verbose_) cout<<"Fill LeptTree "<<endl;
    double deltaPhi_lep0_jet = fabs( deltaPhi(corrJet.phi(), lep0_p4.Phi() )) ;
    double deltaPhi_lep1_jet = fabs( deltaPhi(corrJet.phi(), lep1_p4.Phi() )) ;
    // AK8 jet should be in opposite hemisphere from lepton. If leading jet matches then use it. If it doensn't then check the second leading jet.
    if ( ((count_AK8CHS==0&& deltaPhi_lep_jet >=3.14/2) || (count_AK8CHS==1&&deltaPhi_lep_jet >=3.14/2)) && count_lep >=2 && count_fill_leptTree==0 ){
      count_fill_leptTree++;
      DeltaRJetLep0                          = deltaR(corrJet.eta(), corrJet.phi(), lep0_p4.Eta(), lep0_p4.Phi() );
      DeltaPhiJetLep0                        = deltaPhi_lep0_jet;
      DeltaRJetLep1                          = deltaR(corrJet.eta(), corrJet.phi(), lep1_p4.Eta(), lep1_p4.Phi() );
      DeltaPhiJetLep1                        = deltaPhi_lep1_jet;
      JetPtRaw                              = uncorrJet.pt()      ;
      JetEtaRaw                             = uncorrJet.eta()     ;
      JetPhiRaw                             = uncorrJet.phi()     ;
      JetMassRaw                            = uncorrJet.mass()    ;
      JetP                                  = corrJet.P()         ;
      JetPt                                 = corrJet.pt()        ;
      JetEta                                = corrJet.eta()       ;
      JetPhi                                = corrJet.phi()       ;
      JetRap                                = corrJet.Rapidity()  ;
      JetEnergy                             = corrJet.energy()    ;
      JetMass                               = corrJet.mass()      ;
      JetArea                               = ijet.jetArea()      ;
      JetSDmass                             = softDropMass;                           // Full softDrop uncorrected 4-vector
      JetSDmassRaw                          = sumSDsubjets_P4_uncorr          .M()    ; // Full softDrop uncorrected 4-vector
      JetSDmassCorrL23                      = sumSDsubjets_P4_L23res          .M()    ;
      JetSDmassCorrL23Up                    = sumSDsubjets_P4_L23resCorrUp   .M()    ;
      JetSDmassCorrL23Dn                    = sumSDsubjets_P4_L23resCorrDn   .M()    ;
      JetSDmassCorrL123                     = sumSDsubjets_P4_L123res         .M()    ;
      JetSDmassCorrL123Up                   = sumSDsubjets_P4_L123resCorrUp  .M()    ;
      JetSDmassCorrL123Dn                   = sumSDsubjets_P4_L123resCorrDn  .M()    ;
      JetSDmassCorrL23Smear                    = sumSDsubjets_P4_L23resPtSmear   .M()    ;
      JetSDmassCorrL23SmearUp                  = sumSDsubjets_P4_L23resPtSmearUp .M()    ;
      JetSDmassCorrL23SmearDn                  = sumSDsubjets_P4_L23resPtSmearDn .M()    ;
      JetSDptRaw                            = sumSDsubjets_P4_uncorr          .Perp() ;  // Full softDrop uncorrected 4-vector
      JetSDptCorrL23                        = sumSDsubjets_P4_L23res          .Perp() ;
      JetSDptCorrL23Up                      = sumSDsubjets_P4_L23resCorrUp   .Perp() ;
      JetSDptCorrL23Dn                      = sumSDsubjets_P4_L23resCorrDn   .Perp() ;
      JetSDptCorrL123                       = sumSDsubjets_P4_L123res         .Perp() ;
      JetSDptCorrL123Up                     = sumSDsubjets_P4_L123resCorrUp  .Perp() ;
      JetSDptCorrL123Dn                     = sumSDsubjets_P4_L123resCorrDn  .Perp() ;
      JetSDptCorrL23Smear                      = sumSDsubjets_P4_L23resPtSmear   .Perp() ;
      JetSDptCorrL23SmearUp                    = sumSDsubjets_P4_L23resPtSmearUp .Perp() ;
      JetSDptCorrL23SmearDn                    = sumSDsubjets_P4_L23resPtSmearDn .Perp() ;
      JetSDetaRaw                           = sumSDsubjets_P4_uncorr.Eta()     ;      // Full softDrop uncorrected 4-vector
      JetSDphiRaw                           = sumSDsubjets_P4_uncorr.Phi()     ;      // Full softDrop uncorrected 4-vector
      JetMassPruned                         = prunedMass ;
      JetMassTrimmed                        = trimmedMass ;
      JetTau1                               = tau1 ;
      JetTau2                               = tau2 ;
      JetTau3                               = tau3 ;
      JetTau4                               = tau4 ;
      JetTau32                              = tau32 ;
      JetTau21                              = tau21 ;
      JetSDsubjet0bdisc                     = sub0_bdisc ;
      JetSDsubjet1bdisc                     = sub1_bdisc ;
      JetSDmaxbdisc                         = maxbdisc;
      JetSDmaxbdiscflavHadron               = maxbdiscflav_hadron;
      JetSDmaxbdiscflavParton               = maxbdiscflav_parton;
      JetSDsubjet0pt                        = sub0_P4_uncorr.Pt() ;
      JetSDsubjet0mass                      = sub0_P4_uncorr.M()  ;
      JetSDsubjet0eta                       = sub0_P4_uncorr.Eta()  ;
      JetSDsubjet0phi                       = sub0_P4_uncorr.Phi()  ;
      JetSDsubjet0area                      = sub0_area ;
      JetSDsubjet0flavHadron                = sub0_flav_hadron ;
      JetSDsubjet0flavParton                = sub0_flav_parton ;
      JetSDsubjet0tau1                      = sub0_tau1                  ;
      JetSDsubjet0tau2                      = sub0_tau2                  ;
      JetSDsubjet0tau3                      = sub0_tau3                  ;
      JetSDsubjet1pt                        = sub1_P4_uncorr.Pt() ;
      JetSDsubjet1mass                      = sub1_P4_uncorr.M()  ;
      JetSDsubjet1eta                       = sub1_P4_uncorr.Eta()  ;
      JetSDsubjet1phi                       = sub1_P4_uncorr.Phi()  ;
      JetSDsubjet1area                      = sub1_area ;
      JetSDsubjet1flavHadron                = sub1_flav_hadron ;
      JetSDsubjet1flavParton                = sub1_flav_parton ;
      JetSDsubjet1tau1                      = sub1_tau1                  ;
      JetSDsubjet1tau2                      = sub1_tau2                  ;
      JetSDsubjet1tau3                      = sub1_tau3                  ;

      AK8jet_Lept_P4corr.SetPtEtaPhiM( corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );

      JetPuppiP                             = puppi_p    ;
      JetPuppiPt                            = puppi_pt   ;
      JetPuppiEta                           = puppi_eta  ;
      JetPuppiPhi                           = puppi_phi  ;
      JetPuppiMass                          = puppi_mass ;

      Jet1PuppiSDmass                       = sumPUPsubjets_P4_uncorr           .M()   ;
      Jet1PuppiSDmassCorr                   = sumPUPsubjets_P4_L23res           .M()   ;
      Jet1PuppiSDmassCorrUp                 = sumPUPsubjets_P4_L23resCorrUp     .M()   ;
      Jet1PuppiSDmassCorrDn                 = sumPUPsubjets_P4_L23resCorrDn     .M()   ;
      Jet1PuppiSDmassCorrL23Smear           = sumPUPsubjets_P4_L23resPtSmear    .M()   ;
      Jet1PuppiSDmassCorrL23SmearUp         = sumPUPsubjets_P4_L23resPtSmearUp  .M()   ;
      Jet1PuppiSDmassCorrL23SmearDn         = sumPUPsubjets_P4_L23resPtSmearDn  .M()   ;
      Jet1PuppiSDpt                         = sumPUPsubjets_P4_uncorr           .Perp();
      Jet1PuppiSDptCorr                     = sumPUPsubjets_P4_L23res           .Perp();
      Jet1PuppiSDptCorrUp                   = sumPUPsubjets_P4_L23resCorrUp     .Perp();
      Jet1PuppiSDptCorrDn                   = sumPUPsubjets_P4_L23resCorrDn     .Perp();
      Jet1PuppiSDptCorrL23Smear             = sumPUPsubjets_P4_L23resPtSmear    .Perp();
      Jet1PuppiSDptCorrL23SmearUp           = sumPUPsubjets_P4_L23resPtSmearUp  .Perp();
      Jet1PuppiSDptCorrL23SmearDn           = sumPUPsubjets_P4_L23resPtSmearDn  .Perp();
      Jet1PuppiSDeta                        = sumPUPsubjets_P4_uncorr           .Eta() ;
      Jet1PuppiSDphi                        = sumPUPsubjets_P4_uncorr           .Phi() ;

      JetPuppiTau1                          = puppi_tau1       ;
      JetPuppiTau2                          = puppi_tau2       ;
      JetPuppiTau3                          = puppi_tau3       ;
      JetPuppiTau4                          = puppi_tau4       ;
      JetPuppiTau32                         = puppi_tau32      ;
      JetPuppiTau21                         = puppi_tau21      ;
      JetPuppiSDsubjet0bdisc                = pup0_bdisc       ;
      JetPuppiSDsubjet1bdisc                = pup1_bdisc       ;
      JetPuppiSDmaxbdisc                    = pup_maxbdisc     ;
      JetPuppiSDmaxbdiscflavHadron          = pup_maxbdiscflav_hadron    ;
      JetPuppiSDmaxbdiscflavParton          = pup_maxbdiscflav_parton    ;
      JetPuppiSDsubjet0pt                   = pup0_P4_uncorr.Pt()        ;
      JetPuppiSDsubjet0mass                 = pup0_P4_uncorr.M()         ;
      JetPuppiSDsubjet0eta                  = pup0_P4_uncorr.Eta()         ;
      JetPuppiSDsubjet0phi                  = pup0_P4_uncorr.Phi()         ;
      JetPuppiSDsubjet0area                 = pup0_area                  ;
      JetPuppiSDsubjet0flavHadron           = pup0_flav_hadron           ;
      JetPuppiSDsubjet0flavParton           = pup0_flav_parton           ;
      JetPuppiSDsubjet0tau1                 = pup0_tau1                  ;
      JetPuppiSDsubjet0tau2                 = pup0_tau2                  ;
      JetPuppiSDsubjet0tau3                 = pup0_tau3                  ;
      JetPuppiSDsubjet1pt                   = pup1_P4_uncorr.Pt()        ;
      JetPuppiSDsubjet1mass                 = pup1_P4_uncorr.M()         ;
      JetPuppiSDsubjet1eta                  = pup1_P4_uncorr.Eta()       ;
      JetPuppiSDsubjet1phi                  = pup1_P4_uncorr.Phi()       ;
      JetPuppiSDsubjet1area                 = pup1_area                  ;
      JetPuppiSDsubjet1flavHadron           = pup1_flav_hadron           ;
      JetPuppiSDsubjet1flavParton           = pup1_flav_parton           ;
      JetPuppiSDsubjet1tau1                 = pup1_tau1                  ;
      JetPuppiSDsubjet1tau2                 = pup1_tau2                  ;
      JetPuppiSDsubjet1tau3                 = pup1_tau3                  ;

      JetCHF                                = ijet.chargedHadronEnergy() / uncorrJet.E()  ;
      JetNHF                                = ijet.neutralHadronEnergy() / uncorrJet.E()  ;
      JetCM                                 = ijet.chargedMultiplicity()  ;
      JetNM                                 = ijet.neutralMultiplicity()  ;
      JetNEF                                = ijet.neutralEmEnergy() / uncorrJet.E()  ;
      JetCEF                                = ijet.chargedEmEnergy() / uncorrJet.E()  ;
      JetMF                                 = ijet.muonEnergy() / uncorrJet.E()  ;
      JetMult                               = ijet.numberOfDaughters() ;


      JetPuppiCHF                                = puppi_CHF   ;
      JetPuppiNHF                                = puppi_NHF   ;
      JetPuppiCM                                 = puppi_CM    ;
      JetPuppiNM                                 = puppi_NM    ;
      JetPuppiNEF                                = puppi_NEF   ;
      JetPuppiCEF                                = puppi_CEF   ;
      JetPuppiMF                                 = puppi_MF    ;
      JetPuppiMult                               = puppi_Mult  ;

      JetMassCorrFactor                     = corr_factor_L23res ;
      JetMassCorrFactorUp                   = corrUp_L23 ;
      JetMassCorrFactorDn                   = corrDn_L23 ;
      JetCorrFactor                         = corr ;
      JetCorrFactorUp                       = corrUp_L123 ;
      JetCorrFactorDn                       = corrDn_L123;
      JetPtSmearFactor                      = ptsmear  ;
      JetPtSmearFactorUp                    = ptsmearUp;
      JetPtSmearFactorDn                    = ptsmearDn;
      JetPuppiMassCorrFactor                = 1;
      JetPuppiMassCorrFactorUp              = 1;
      JetPuppiMassCorrFactorDn              = 1;
      JetPuppiCorrFactor                    = 1;
      JetPuppiCorrFactorUp                  = 1;
      JetPuppiCorrFactorDn                  = 1;
      JetPuppiPtSmearFactor                 = 1;
      JetPuppiPtSmearFactorUp               = 1;
      JetPuppiPtSmearFactorDn               = 1;
      JetEtaScaleFactor                     = 1;
      JetPhiScaleFactor                     = 1;
      // JetMatchedGenJetDR                    = GenJetMatched_dRmin;
      JetMatchedGenJetPt                    = GenJetMatched.Perp();
      JetMatchedGenJetMass                  = GenJetMatched.M();
    }

    count_AK8CHS++;
  }

  //
  //  .d8888b.                         d8b        888                        888        88888888888
  // d88P  Y88b                        Y8P        888                        888            888
  // Y88b.                                        888                        888            888
  //  "Y888b.    .d88b.  88888b.d88b.  888        888       .d88b.  88888b.  888888         888     888d888  .d88b.   .d88b.
  //     "Y88b. d8P  Y8b 888 "888 "88b 888        888      d8P  Y8b 888 "88b 888            888     888P"   d8P  Y8b d8P  Y8b
  //       "888 88888888 888  888  888 888 888888 888      88888888 888  888 888            888     888     88888888 88888888
  // Y88b  d88P Y8b.     888  888  888 888        888      Y8b.     888 d88P Y88b.          888     888     Y8b.     Y8b.
  //  "Y8888P"   "Y8888  888  888  888 888        88888888  "Y8888  88888P"   "Y888         888     888      "Y8888   "Y8888
  //                                                                888
  //                                                                888
  //                                                                888

  LeptMETpx                = met.px();
  LeptMETpy                = met.py();
  LeptMETpt                = met.pt();
  LeptMETphi               = met.phi();
  LeptMETsumET             = met.sumEt();
  LeptNvtx                 = nvtx;
  LeptNPUtrue              = nPU;
  LeptRho                  = rho ;
  LeptEventWeight          = 1 ;
  LeptPUweight       = puweight  ;
  LeptPUweight_MBup  = puweightUp ;
  LeptPUweight_MBdn  = puweightDn  ;

  double htlep0 = lep0_p4.Perp() + met0.pt() ;
  double htlep1 = lep1_p4.Perp() + met1.pt() ;

  HTlep0                = htlep0 ;
  ST0                   = htlep0 + HT_AK4_pt30           ;
  ST0_CorrDn            = htlep0 + HT_AK4_pt30_corrDn    ;
  ST0_CorrUp            = htlep0 + HT_AK4_pt30_corrUp    ;
  ST0_PtSmearNom        = htlep0 + HT_AK4_pt30_smearNom  ;
  ST0_PtSmearUp         = htlep0 + HT_AK4_pt30_smearUp   ;
  ST0_PtSmearDn         = htlep0 + HT_AK4_pt30_smearDn   ;

  HTlep1                = htlep1 ;
  ST1                   = htlep1 + HT_AK4_pt30           ;
  ST1_CorrDn            = htlep1 + HT_AK4_pt30_corrDn    ;
  ST1_CorrUp            = htlep1 + HT_AK4_pt30_corrUp    ;
  ST1_PtSmearNom        = htlep1 + HT_AK4_pt30_smearNom  ;
  ST1_PtSmearUp         = htlep1 + HT_AK4_pt30_smearUp   ;
  ST1_PtSmearDn         = htlep1 + HT_AK4_pt30_smearDn   ;


  LeptRunNum               = iEvent.id().run() ;
  LeptLumiBlock            = iEvent.id().luminosityBlock() ;
  LeptEventNum             = iEvent.id().event() ;
  if(passMETfilters) LeptPassMETFilters  = 1;
  else LeptPassMETFilters  = 0;

  Lepton0Phi   = lep0_p4.Phi()  ;
  Lepton0Pt    = lep0_p4.Perp() ;
  Lepton0Eta   = lep0_p4.Eta()  ;
  Lepton0Mass  = lep0_p4.M() ;

  Lepton1Phi   = lep1_p4.Phi()  ;
  Lepton1Pt    = lep1_p4.Perp() ;
  Lepton1Eta   = lep1_p4.Eta()  ;
  Lepton1Mass  = lep1_p4.M() ;


  if (count_mu==2 && count_el==0)      LeptonIsMu  = 1  ;
  else if (count_mu==0 && count_el==2) LeptonIsMu  = 0  ;
  else LeptonIsMu =0;

  Mu0Iso  = mu0_iso04;
  Mu1Iso  = mu1_iso04;


  Elecron0_absiso            = el0_absiso           ;
  Elecron0_relIsoWithDBeta   = el0_relIsoWithDBeta  ;
  Elecron0_absiso_EA         = el0_absiso_EA        ;
  Elecron0_relIsoWithEA      = el0_relIsoWithEA     ;

  Elecron1_absiso            = el1_absiso           ;
  Elecron1_relIsoWithDBeta   = el1_relIsoWithDBeta  ;
  Elecron1_absiso_EA         = el1_absiso_EA        ;
  Elecron1_relIsoWithEA      = el1_relIsoWithEA     ;


  if (el0_isMedium) Electron0_isMedium  = 1         ;
  else              Electron0_isMedium  = 0;
  if (el0_isTight)  Electron0_isTight   = 1         ;
  else              Electron0_isTight   = 0;

  if (el1_isMedium) Electron1_isMedium  = 1         ;
  else              Electron1_isMedium  = 0;
  if (el1_isTight)  Electron1_isTight   = 1         ;
  else              Electron1_isTight   = 0;


  if(mu0_isMedium) Mu0Medium = 1   ;
  else             Mu0Medium = 0   ;

  if(mu0_isTight) Mu0Tight = 1   ;
  else            Mu0Tight = 0   ;

  if(mu0_isHighPt) Mu0HighPt = 1   ;
  else             Mu0HighPt = 0   ;


  if(mu1_isMedium) Mu1Medium = 1   ;
  else             Mu1Medium = 0   ;

  if(mu1_isTight) Mu1Tight = 1   ;
  else            Mu1Tight = 0   ;

  if(mu1_isHighPt) Mu1HighPt = 1   ;
  else             Mu1HighPt = 0   ;


  //------------------------------------
  // WRITE TREE WITH BASELINE PT CUT AND ETA CUT
  //------------------------------------


  if (GenTruth_Leptonic)  count_GenTruth_Leptonic ++;
  if (count_mu  >=1 )  count_nMu_gt1 ++;
  if (count_el  >=1 )  count_nEl_gt1 ++;
  if (count_mu  ==2 )  count_nMu_e1 ++;
  if (count_el  ==2 )  count_nEl_e1 ++;
  if (count_lep ==2 )  count_nLep_e1 ++;
  if (count_lep ==2  && AK8jet_Lept_P4corr.Perp()>300)  count_JetPt300 ++;
  if (count_lep ==2  && AK8jet_Lept_P4corr.Perp()>300 && fabs( AK8jet_Lept_P4corr.Rapidity() ) <2.4 )  count_JetPt300Eta ++;
  if (count_lep ==2  && AK8jet_Lept_P4corr.Perp()>300 && fabs( AK8jet_Lept_P4corr.Rapidity() ) <2.4 &&  mu0_p4.Perp()>40 && mu1_p4.Perp() > 40)  count_JetPt300Eta_muPt40 ++;
  if (count_lep ==2  && AK8jet_Lept_P4corr.Perp()>300 && fabs( AK8jet_Lept_P4corr.Rapidity() ) <2.4 &&  mu0_p4.Perp()>40 && met.pt() > 40 && mu1_p4.Perp() > 40)  count_JetPt300Eta_muPt40_MET40 ++;

  if (count_lep == 2  && verbose_){
    cout<<" ak8pt "<<AK8jet_Lept_P4corr.Perp()<<endl;
    cout<<" mu pt "<<mu0_p4.Perp()<<endl;
    cout<<" el pt "<<el0_p4.Perp()<<endl;
    cout<<" met "<<met.pt() <<endl;
    cout<<" ak4 pt "<<AK4_dRMinLep_p4.Perp() <<endl;
  }
  if (count_lep == 2 && AK8jet_Lept_P4corr.Perp()>200 && fabs( AK8jet_Lept_P4corr.Rapidity() ) <2.4
    &&  lep0_p4.Perp()>30 && met.pt() > 30 ){
    TreeLept -> Fill();
  }


}


// ------------ method called once each job just before starting event loop  ------------
void
B2GTTbarTreeMaker::beginJob()
{
  fPUweight = new TFile("PUweight20160908.root") ;
  hPUweight      = (TH1D*) fPUweight->Get("PUweight_true");
  hPUweight_MBup = (TH1D*) fPUweight->Get("PUweight_true_MBup");
  hPUweight_MBdn = (TH1D*) fPUweight->Get("PUweight_true_MBdn");

  std::cout<<"Test PU reweight file: "<<hPUweight->GetBinContent( hPUweight->GetXaxis()->FindBin( 30 ) )<<std::endl;


  count_GenTruth_Leptonic =0;
  count_nMu_gt1 =0;
  count_nEl_gt1 =0;
  count_nMu_e1 =0;
  count_nEl_e1 =0;
  count_nLep_e1 =0;
  count_JetPt300 =0;
  count_JetPt300Eta =0;
  count_JetPt300Eta_AK4 =0;
  count_JetPt300Eta_muPt40 =0;
  count_JetPt300Eta_muPt40_MET40 =0;
  count_JetPt300Eta_muPt40_MET40_AK4 =0;

}

// ------------ method called once each job just after ending the event loop  ------------
void
B2GTTbarTreeMaker::endJob()
{

  std::cout<<" nEvents GenTruth Leptonic  :" <<count_GenTruth_Leptonic<<std::endl;
  std::cout<<" nEvents nMu =1   :" <<count_nMu_e1 <<std::endl;
  std::cout<<" nEvents nEl =1   :" <<count_nEl_e1 <<std::endl;
  std::cout<<" nEvents nLepton =1   :" <<count_nLep_e1 <<std::endl;
  std::cout<<" nEvents nLepton =1 && JetPt300   :" <<count_JetPt300 <<std::endl;
  std::cout<<" nEvents nLepton =1 && JetPt300Eta   :" <<count_JetPt300Eta <<std::endl;
  std::cout<<" nEvents nLepton =1 && JetPt300Eta && AK4pt>20   :" <<count_JetPt300Eta_AK4 <<std::endl;
  std::cout<<" nEvents nLepton =1 && JetPt300Eta && muPt40   :" <<count_JetPt300Eta_muPt40 <<std::endl;
  std::cout<<" nEvents nLepton =1 && JetPt300Eta && muPt40 && MET40   :" <<count_JetPt300Eta_muPt40_MET40 <<std::endl;
  std::cout<<" nEvents nLepton =1 && JetPt300Eta && muPt40 && MET40 && AK4pt>20  :" <<count_JetPt300Eta_muPt40_MET40_AK4 <<std::endl;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
B2GTTbarTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(B2GTTbarTreeMaker);

//  LocalWords:  NNPDF
