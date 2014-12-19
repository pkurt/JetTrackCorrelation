#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <TF1.h>

//#include "class_def/JetAna.h"
//#include "class_def/Tracks.h"
//#include "class_def/HLT.h"
//#include "class_def/HiTree.h"
//#include "class_def/Skim.h"

#include "class_def/mixing_tree.h"
#include "class_def/mixing_tree_MB.h"


#include <algorithm>


#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TVector3.h"

bool is_data = false; 
bool test_it = false;
bool is_selected = true;
bool is_dijet = true;
bool is_leadingjet = false; 
bool is_subleadingjet = true;

#define nCBins 4
#define nPtBins 1
//#define nPtBins 
//#define nTrkPtBins 4
//#define nTrkPtBins 1
#define nTrkPtBins 5


#define nPtBins_sub 10

//#include "TrackingCorrections2012.h"
//#include "TrackingCorrectionsv6.h"

bool doTrackCorrections = true;
//bool doTrackCorrections = false;
//vector<TrackingCorrections*> trackCorrections;

float trkPtCut=1.0;

enum enum_dataset_types {e_Data2011,e_HydJet30,e_HydJet50,e_HydJet80, e_HydJet100, e_HydJet120,e_HydJet170,e_HydJet200,e_HydJet250, e_HydJet300, e_n_dataset_types};
TString dataset_type_strs[e_n_dataset_types] = { "Data2011","HydJet30","HydJet50","HydJet80", "HydJet100", "HydJet120","HydJet170","HydJet200","HydJet250","HydJet300" };
int dataset_type_code = -999;


///float trkPtCut=1.;
/*
enum enum_npart_types {e_0_25, e_25_50, e_50_75, e_75_100, e_100_125, e_125_150, e_150_175, e_175_200,e_200_225,e_225_250,e_250_275, e_275_300,  e_300_325, e_325_350, e_350_375, e_375_400,  e_400_425, e_425_450, e_450_475, e_475_500, e_500_525, e_525_550, e_550_575, e_575_600, e_600_625, e_625_650, e_650_675, e_675_700, e_700_725, e_725_750, e_750_775, e_775_800,  e_AllParts, e_n_npart_types};

TString npart_type_strs[e_n_npart_types] = {"part1", "part2", "part3", "part4", "part5", "part6", "part7", "part8", "part9", "part10","part11", "part12","part13","part14","part15","part16","part17","part18","part19","part20","part21","part22","part23","part24","part25","part26","part27", "part28", "part29", "part30", "part31", "part32",  "AllParts" };
*/




enum enum_npart_types {
  e_0_10, e_10_20, e_20_30, e_30_40, e_40_50, e_50_60, e_60_70, e_70_80,e_80_90,e_90_100,
  e_100_110, e_110_120, e_120_130, e_130_140, e_140_150, e_150_160, e_160_170, e_170_180,e_180_190,e_190_200,
  e_200_210, e_210_220, e_220_230, e_230_240, e_240_250, e_250_260, e_260_270, e_270_280,e_280_290,e_290_300,
  e_300_310, e_310_320, e_320_330, e_330_340, e_340_350, e_350_360, e_360_370, e_370_380,e_380_390,e_390_400,
  e_400_410, e_410_420, e_420_430, e_430_440, e_440_450, e_450_460, e_460_70, e_470_480,e_480_490,e_490_500,
  e_500_510, e_510_520, e_520_530, e_530_540, e_540_550, e_550_560, e_560_570, e_570_580,e_580_590,e_590_600,
  e_600_610, e_610_620, e_620_630, e_630_640, e_640_650, e_650_660, e_660_670, e_670_680,e_680_690,e_690_700,
  e_700_710, e_710_720, e_720_730, e_730_740, e_740_750, e_750_760, e_760_770, e_770_780,e_780_790,e_790_800,
  e_AllParts, e_n_npart_types};

TString npart_type_strs[e_n_npart_types] = {
  "part1",  "part2", "part3", "part4", "part5", "part6", "part7", "part8", "part9", "part10",
  "part11", "part12","part13","part14","part15","part16","part17","part18","part19","part20",
  "part21", "part22","part23","part24","part25","part26","part27", "part28", "part29", "part30",
  "part31", "part32", "part33", "part34", "part35", "part36", "part37", "part38", "part39", "part40",
  "part41",  "part42", "part43", "part44", "part45", "part46", "part47", "part48", "part49", "part50",
  "part51",  "part52", "part53", "part54", "part55", "part56", "part57", "part58", "part59", "part60",
  "part61",  "part62", "part63", "part64", "part65", "part66", "part67", "part68", "part69", "part70",
  "part71",  "part72", "part73", "part74", "part75", "part76", "part77", "part78", "part79", "part80",
  "AllParts" };




int npart = -999;

float CBins[nCBins+1] = {0, 20, 60, 100 , 200};
TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};


//float CBins[nCBins+1] = {0, 4, 12, 20, 40};
//TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};

//float PtBins[nPtBins+1] = {100, 110, 120, 130, 150, 200, 300};
//TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt110", "Pt120","Pt130", "Pt150", "Pt200","Pt300"};

float PtBins[nPtBins+1] = {50,300};
TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};

//float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 5, 10};
//TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1","TrkPt2", "TrkPt3", "TrkPt5", "TrkPt10" };

float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};
TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt999" };


//float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 5, 6, 7, 8, 10};
//TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1","TrkPt2", "TrkPt3", "TrkPt5", "TrkPt6", "TrkPt7", "TrkPt8","TrkPt10" };


//float TrkPtBins[nTrkPtBins+1] = {1, 999};
//TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt999" };


float PtBins_sub[nPtBins_sub+1] = {50, 60, 70, 80, 100, 120, 140, 160, 200, 300, 500};
TString PtBin_sub_strs[nPtBins_sub+1] = {"Pt50", "Pt60", "Pt70", "Pt80",  "Pt100", "Pt120", "Pt140","Pt160", "Pt200" , "Pt300", "Pt500"};

const int nPtBintrk = 27;
Double_t ptBintrk[nPtBintrk+1] = {0, 0.5, 1, 1.203915, 1.449412, 1.74497, 2.100796, 2.529181, 3.04492, 3.665826, 4.413344, 5.313293, 6.396755, 7.701152, 9.271536, 11.16214, 13.43828, 16.17855, 19.47761, 23.44939, 28.23108, 33.98783, 40.91848, 49.26238, 59.30774, 71.40151, 85.96137, 103.4902};


class hist_class {
  public:
    hist_class(TString the_desc, bool is_it_data);
    void Delete();
    void NormalizeByCounts();
    void Write();
    void AddHists(hist_class *more_hists, float wt);
    bool is_data;
   // bool is_selected;
    TString desc;
    int n_evt_raw ;
    TH1F* NumberOfMatches;
    TH1F* NEvents;
    TH1F* NEvents_test;
    TH1F* NEvents_after_noise;
    TH1F* NEvents_after_spike;
    TH1F* NEvents_after_dphi;
    TH1F* NEvents_before_dphi;
    TH1F* NEvents_dijets;
    TH1F* NEvents_after_trigger;
    
    TH1F* track_vz;
    TH1F* track_vz_weighted;
    
    TH1F* AllJetPt_raw_hist;
    TH1F* AllJetPhi_hist;
    TH1F* AllJetEta_hist;
    TH1F* AllJetPt_hist;
    TH1F* First_AllJetPhi_hist;
    TH1F* First_AllJetEta_hist;
    TH1F* First_AllJetPt_hist;
    TH1F* Sub_AllJetPhi_hist;
    TH1F* Sub_AllJetEta_hist;
    TH1F* Sub_AllJetPt_hist;
    TH1F* JetPt_fraction;
    TH1F* Centrality;
    
    
    TH1F* jet_pT_hist[nCBins];
    TH1F* jet_phi_hist[nCBins];
    TH1F* jet_eta_hist[nCBins];
    TH1F* jet_corrpT_hist[nCBins];
    TH1F* LeadingJetPt_hist[nCBins];
    TH1F* LeadingJetPhi_hist[nCBins];
    TH1F* LeadingJetEta_hist[nCBins];
    TH1F* SubJetPt_hist[nCBins];
    TH1F* SubJetPhi_hist[nCBins];
    TH1F* SubJetEta_hist[nCBins];
    
    TH1F* ThirdJetPt_hist[nCBins][nPtBins];
    TH1F* ThirdJetPhi_hist[nCBins][nPtBins];
    TH1F* ThirdJetEta_hist[nCBins][nPtBins];
    TH1F* all_cand_pT_hist[nCBins][nPtBins];
    TH1F* all_cand_phi_hist[nCBins][nPtBins];
    TH1F* all_cand_eta_hist[nCBins][nPtBins];
    
    TH1F* track_cand_pT_hist[nCBins][nPtBins];
    TH1F* track_cand_phi_hist[nCBins][nPtBins];
    TH1F* track_cand_eta_hist[nCBins][nPtBins];
    TH1F* track_cand_pT_hist_subleadingjet[nCBins][nPtBins];
    
    
    //TH1F* dPhi_hist[nCBins][nPtBins];
    //TH1F* dPhi_after_hist[nCBins][nPtBins];
    TH1F* dPhi_hist[nCBins];
    TH1F* dPhi_after_hist[nCBins];

    
    TH1F* dPhi_leadingjet_hist[nCBins];
    TH1F* dPhi_subleadingjet_hist[nCBins];
    
    TH1F* all_jets_corrpT[nCBins][nPtBins];
    TH1F* all_jets_phi[nCBins][nPtBins];
    TH1F* all_jets_eta[nCBins][nPtBins];
    TH1F* only_leadingjets_corrpT[nCBins][nPtBins];
    TH1F* only_leadingjets_phi[nCBins][nPtBins];
    TH1F* only_leadingjets_eta[nCBins][nPtBins];
    TH1F* only_subleadingjets_corrpT[nCBins][nPtBins];
    TH1F* only_subleadingjets_phi[nCBins][nPtBins];
    TH1F* only_subleadingjets_eta[nCBins][nPtBins];
    TH1F* neutral_cand_pT_hist[nCBins][nPtBins];
    TH1F* neutral_cand_phi_hist[nCBins][nPtBins];
    TH1F* neutral_cand_eta_hist[nCBins][nPtBins];
    TH1F* photons_cand_pT_hist[nCBins][nPtBins];
    TH1F* photons_cand_phi_hist[nCBins][nPtBins];
    TH1F* photons_cand_eta_hist[nCBins][nPtBins];
    
    TH1F* NumNeutral[nCBins][nPtBins];
    TH1F* NumPhotons[nCBins][nPtBins];
    TH1F* NumAll[nCBins][nPtBins];
    TH1F* NumCharged[nCBins][nPtBins];
    
    TH1F* NumAll_bkg[nCBins][nPtBins];
    TH1F* NumNeutral_bkg[nCBins][nPtBins];
    TH1F* NumPhotons_bkg[nCBins][nPtBins];
    TH1F* NumChargedHadrons_bkg[nCBins][nPtBins];
    TH1F* NumChargedParticles_bkg[nCBins][nPtBins];
    TH1F* NumElectrons_bkg[nCBins][nPtBins];
    TH1F* NumMuons_bkg[nCBins][nPtBins];
    
    TH1F* radius_hist[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist[nCBins][nPtBins];
    
    TProfile *JetShapeIntegratedParticles[nCBins][nPtBins];
    TProfile *JetShapeDiffParticles[nCBins][nPtBins];
    TProfile *JetShapeIntegratedParticles_bkgsub[nCBins][nPtBins];
    TProfile *JetShapeDiffParticles_bkgsub[nCBins][nPtBins];
    TProfile *JetShapeIntegratedParticles_bkg[nCBins][nPtBins];
    TProfile *JetShapeDiffParticles_bkg[nCBins][nPtBins];
    
    TH1F *JetShapeDiffParticles_bkg_1D[nCBins][nPtBins];
    //  TH1F *JetShapeDiffParticles_bkg_1D[nCBins][nPtBins][nTrkPtBins];
    
    
    TH1F *JetShapeDiffParticles_bkg_1D_EtaRef[nCBins][nPtBins];
    TH1F *JetShapeDiffParticles_1D[nCBins][nPtBins];
    //TH1F *JetShapeDiffParticles_1D[nCBins][nPtBins][nTrkPtBins];
    
    TH1F *JetShapeIntegratedParticles_bkg_1D[nCBins][nPtBins];
    TH1F *JetShapeIntegratedParticles_1D[nCBins][nPtBins];
    
    
    TH1F* NumNeutral_subleadingjet[nCBins][nPtBins];
    TH1F* NumPhotons_subleadingjet[nCBins][nPtBins];
    TH1F* NumAll_subleadingjet[nCBins][nPtBins];
    TH1F* NumChargedHadrons_subleadingjet[nCBins][nPtBins];
    
    
    TH1F* SumJetPtFraction_hist_subleadingjet[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_photons[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_neutrals[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_chargedhadrons[nCBins][nPtBins];
    
    TH1F* SumJetPtFraction_hist_leadingjet_photons[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_neutrals[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_charged[nCBins][nPtBins];
    
    TH1F* all_cand_pT_hist_subleadingjet[nCBins][nPtBins];
    TH1F* all_cand_phi_hist_subleadingjet[nCBins][nPtBins];
    TH1F* all_cand_eta_hist_subleadingjet[nCBins][nPtBins];
    TH1F* radius_hist_subleadingjet[nCBins][nPtBins];
    
    TH1F* Centrality_hist[nCBins][nPtBins];
    //TH1F* Aj[nCBins][nPtBins];
    TH1F* Aj[nCBins];

    
    TH1F* neutral_cand_pT_hist_subleadingjet[nCBins][nPtBins];
    TH1F* neutral_cand_phi_hist_subleadingjet[nCBins][nPtBins];
    TH1F* neutral_cand_eta_hist_subleadingjet[nCBins][nPtBins];
    
    TH1F* photons_cand_pT_hist_subleadingjet[nCBins][nPtBins];
    TH1F* photons_cand_phi_hist_subleadingjet[nCBins][nPtBins];
    TH1F* photons_cand_eta_hist_subleadingjet[nCBins][nPtBins];
    
    
    
    TH1F* SumJetPtFraction_hist_leadingjet_chargedhadrons[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_chargedparticles[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_electrons[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_muons[nCBins][nPtBins];
    
    TH1F* NumChargedHadrons[nCBins][nPtBins];
    TH1F* NumChargedParticles[nCBins][nPtBins];
    TH1F* NumElectrons[nCBins][nPtBins];
    TH1F* NumMuons[nCBins][nPtBins];
    
    TH1F* chargedhadrons_cand_pT_hist[nCBins][nPtBins];
    TH1F* chargedhadrons_cand_phi_hist[nCBins][nPtBins];
    TH1F* chargedhadrons_cand_eta_hist[nCBins][nPtBins];
    
    TH1F* chargedparticles_cand_pT_hist[nCBins][nPtBins];
    TH1F* chargedparticles_cand_phi_hist[nCBins][nPtBins];
    TH1F* chargedparticles_cand_eta_hist[nCBins][nPtBins];
    
    TH1F* electrons_cand_pT_hist[nCBins][nPtBins];
    TH1F* electrons_cand_phi_hist[nCBins][nPtBins];
    TH1F* electrons_cand_eta_hist[nCBins][nPtBins];
    
    TH1F* muons_cand_pT_hist[nCBins][nPtBins];
    TH1F* muons_cand_phi_hist[nCBins][nPtBins];
    TH1F* muons_cand_eta_hist[nCBins][nPtBins];
    
    
    
    TH1F* chargedhadrons_cand_pT_hist_bkg[nCBins][nPtBins];
    TH1F* chargedparticles_cand_pT_hist_bkg[nCBins][nPtBins];
    TH1F* electrons_cand_pT_hist_bkg[nCBins][nPtBins];
    TH1F* muons_cand_pT_hist_bkg[nCBins][nPtBins];
    TH1F* neutral_cand_pT_hist_bkg[nCBins][nPtBins];
    TH1F* photons_cand_pT_hist_bkg[nCBins][nPtBins];
    
    
    
    TH1F* dN_tracks[nCBins][nPtBins];
    TH1F* dN_chargedhadrons[nCBins][nPtBins];
    TH1F* dN_chargedparticles[nCBins][nPtBins];
    TH1F* dN_electrons[nCBins][nPtBins];
    TH1F* dN_muons[nCBins][nPtBins];
    TH1F* dN_neutrals[nCBins][nPtBins];
    TH1F* dN_photons[nCBins][nPtBins];
    
    TH1F* SumJetPtFraction_hist_leadingjet_photons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_neutrals_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_electrons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_leadingjet_muons_bkg[nCBins][nPtBins];
    
    
    
    /// subleading jet
    TH1F* chargedhadrons_cand_pT_hist_subleadingjet[nCBins][nPtBins];
    TH1F* chargedhadrons_cand_phi_hist_subleadingjet[nCBins][nPtBins];
    TH1F* chargedhadrons_cand_eta_hist_subleadingjet[nCBins][nPtBins];
    
    TH1F* chargedhadrons_cand_pT_hist_subleadingjet_bkg[nCBins][nPtBins];
    TH1F* neutral_cand_pT_hist_subleadingjet_bkg[nCBins][nPtBins];
    TH1F* photons_cand_pT_hist_subleadingjet_bkg[nCBins][nPtBins];
    
    TH1F* SumJetPtFraction_hist_subleadingjet_photons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_neutrals_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_electrons_bkg[nCBins][nPtBins];
    TH1F* SumJetPtFraction_hist_subleadingjet_muons_bkg[nCBins][nPtBins];
    
    TH1F* NumAll_subleadingjet_bkg[nCBins][nPtBins];
    TH1F* NumNeutral_subleadingjet_bkg[nCBins][nPtBins];
    TH1F* NumPhotons_subleadingjet_bkg[nCBins][nPtBins];
    TH1F* NumChargedHadrons_subleadingjet_bkg[nCBins][nPtBins];
    
    TH1F* JetPt_fraction_hist[nCBins][nPtBins];
    TH1F* dPhi_jet_track[nCBins][nPtBins];
    TH1F* dPhi_jet_track_ptweight[nCBins][nPtBins];
    
    
    TH1F* SumPt_in_a_radius[nCBins][nPtBins];
    TH1F* multiplicity_in_a_radius[nCBins][nPtBins];
    
    TH1F* JetShapeDiffParticles_bkg_1D_total[nCBins];
    TH1F* JetShapeIntegratedParticles_bkg_1D_total[nCBins];
    
    TH1F* track_cand_pT_hist_weighted[nCBins][nPtBins];
    
    TH1F* JetShapeDiffParticles_bkg_1DEtaRef[nCBins][nPtBins];
    TH1F* JetShapeDiffParticles_bkg_1DEtaRefSig[nCBins][nPtBins];
    
    
    TH1F* SumPt_bkg_only[nCBins][nPtBins];
    TH1F* SumPt_ratio[nCBins][nPtBins];
    TH1F* SumPt_only[nCBins][nPtBins];
    TH1F* SumPt_only_centbin[nCBins];
    
    TH1F* JetShapeDiffParticles_subleadingjet_bkg_1DEtaRef[nCBins][nPtBins_sub];
    TH1F* JetShapeIntegratedParticles_subleadingjet_1D[nCBins][nPtBins_sub];
    TH1F* JetShapeDiffParticles_subleadingjet_bkg_1D[nCBins][nPtBins_sub];
    TH1F* JetShapeDiffParticles_subleadingjet_1D[nCBins][nPtBins_sub];
    
    double sig_count[nCBins][nPtBins];
    double bkg_count[nCBins][nPtBins];
    
    TH1F* EtaRef_bkg_pt[nCBins][nPtBins];
    TH1F* EtaRef_bkg_pt_weighted[nCBins][nPtBins];
    
    
    TH2D* hJetTrackSignalBackground[nCBins][nPtBins][nTrkPtBins];
    TH2D* hJetTrackBackground[nCBins][nPtBins][nTrkPtBins];
    TH2D* hJetTrackBackground_notrkcorr[nCBins][nPtBins][nTrkPtBins];
    TH2D* hJetGenTrackBackground[nCBins][nPtBins][nTrkPtBins];
    TH2D* hJetGenTrackBackground_pythia[nCBins][nPtBins][nTrkPtBins];
    TH2D* hJetGenTrackBackground_hydjet[nCBins][nPtBins][nTrkPtBins];

    TH1F* Jet_phi[nCBins][nPtBins][nTrkPtBins];

    TH1D* Bkg_phi[nCBins][nPtBins][nTrkPtBins];
    TH1D* Bkg_eta[nCBins][nPtBins][nTrkPtBins];
    
    TH1D* Bkg_phi_weighted[nCBins][nPtBins][nTrkPtBins];
    TH1D* Bkg_eta_weighted[nCBins][nPtBins][nTrkPtBins];
    
    
    
    TH1D* EMix_bkg_pt[nCBins][nPtBins];
    TH1D* EMix_bkg_pt_weighted[nCBins][nPtBins];
    
    TH1D* EMix_bkg_phi[nCBins][nPtBins];
    TH1D* EMix_bkg_phi_weighted[nCBins][nPtBins];
    
    TH1D* EMix_bkg_eta[nCBins][nPtBins];
    TH1D* EMix_bkg_eta_weighted[nCBins][nPtBins];
    
    
    TH1D* track_bkg_pT_hist_weighted[nCBins][nPtBins];
    TH1D* track_bkg_pT_hist[nCBins][nPtBins];
    
    TH1D* track_bkg_phi_hist_weighted[nCBins][nPtBins];
    TH1D* track_bkg_phi_hist[nCBins][nPtBins];
    
    TH1D* track_bkg_eta_hist_weighted[nCBins][nPtBins];
    TH1D* track_bkg_eta_hist[nCBins][nPtBins];

    TH1D* GenTrkEta_bkg[nCBins][nPtBins][nTrkPtBins];
    TH1D* GenTrkPhi_bkg[nCBins][nPtBins][nTrkPtBins];

    TH1D* track_gen_bkg_pT_hist[nCBins][nPtBins];
    TH1D* track_gen_bkg_eta_hist[nCBins][nPtBins];
    TH1D* track_gen_bkg_phi_hist[nCBins][nPtBins];

    TH1F* emix_jetpt[nCBins][nPtBins];
    TH1F* emix_jetpt_test[nCBins][nPtBins];

    

};

hist_class::hist_class(TString the_desc, bool is_it_data)
{



  for(int centi = 0; centi < nCBins; centi++) {
    for(int pti = 0; pti < nPtBins; pti++) {
      sig_count[centi][pti] = 0;
      bkg_count[centi][pti] = 0;

    }
    /*
       for(int pti_sub = 0; pti_sub < nPtBins_sub; pti_sub++) {
       sig_count_sub[centi][pti_sub] = 0;
       bkg_count_sub[centi][pti_sub] = 0;

       }*/

  }




  n_evt_raw = 0;
  desc = the_desc;
  is_data = is_it_data;
  NumberOfMatches = new TH1F((TString) (desc + "_NumberOfMatches"), "", 100, 0., 100.);     NumberOfMatches->Sumw2();
  NEvents = new TH1F((TString) (desc + "_Nevents"), "", 100, 0., 100.);     NEvents->Sumw2();
  NEvents_test = new TH1F((TString) (desc + "_Nevents_test"), "", 100, 0., 100.);     NEvents_test->Sumw2();
  NEvents_after_noise = new TH1F((TString) (desc + "_Nevents_after_noise"), "", 100, 0., 100.);     NEvents_after_noise->Sumw2();
  NEvents_after_spike = new TH1F((TString) (desc + "_Nevents_after_spike"), "", 100, 0., 100.);     NEvents_after_spike->Sumw2();

  NEvents_dijets = new TH1F((TString) (desc + "_Nevents_dijets"), "", 100, 0., 100.);     NEvents_dijets->Sumw2();
  NEvents_after_dphi = new TH1F((TString) (desc + "_Nevents_after_dphi"), "", 100, 0., 100.);     NEvents_after_dphi->Sumw2();
  NEvents_before_dphi = new TH1F((TString) (desc + "_Nevents_before_dphi"), "", 100, 0., 100.);     NEvents_before_dphi->Sumw2();
  NEvents_after_trigger = new TH1F((TString) (desc + "_Nevents_after_trigger"), "", 100, 0., 100.);     NEvents_after_trigger->Sumw2();


  AllJetPt_raw_hist = new TH1F((TString) (desc + "_AllJetPt_raw_hist"), "", 100, 0., 500.);     AllJetPt_raw_hist->Sumw2();
  AllJetPt_hist = new TH1F((TString) (desc + "_AllJetPt_hist"), "", 100, 0., 500.);     AllJetPt_hist->Sumw2();
  AllJetPhi_hist = new TH1F((TString) (desc + "_AllJetPhi_hist"), "", 36, 0.,TMath::Pi());     AllJetPhi_hist->Sumw2();
  AllJetEta_hist = new TH1F((TString) (desc + "_AllJetEta_hist"), "", 50, -5., 5.);     AllJetEta_hist->Sumw2();

  First_AllJetPt_hist = new TH1F((TString) (desc + "_First_AllJetPt_hist"), "", 100, 0., 500.);     First_AllJetPt_hist->Sumw2();
  First_AllJetPhi_hist = new TH1F((TString) (desc + "_First_AllJetPhi_hist"), "", 36, 0.,TMath::Pi());     First_AllJetPhi_hist->Sumw2();
  First_AllJetEta_hist = new TH1F((TString) (desc + "_First_AllJetEta_hist"), "", 50, -5., 5.);     First_AllJetEta_hist->Sumw2();

  Sub_AllJetPt_hist = new TH1F((TString) (desc + "_Sub_AllJetPt_hist"), "", 100, 0., 500.);     Sub_AllJetPt_hist->Sumw2();
  Sub_AllJetPhi_hist = new TH1F((TString) (desc + "_Sub_AllJetPhi_hist"), "", 36, 0.,TMath::Pi());     Sub_AllJetPhi_hist->Sumw2();
  Sub_AllJetEta_hist = new TH1F((TString) (desc + "_Sub_AllJetEta_hist"), "", 50, -5., 5.);     Sub_AllJetEta_hist->Sumw2();
  JetPt_fraction = new TH1F((TString) (desc + "_JetPt_fraction"), "", 100,0.,5);     JetPt_fraction->Sumw2();

  Centrality = new TH1F((TString) (desc + "_Centrality"), "", 40,0.,40);     Centrality->Sumw2();

  track_vz = new TH1F((TString) (desc + "_track_vz"), "", 80, -20., 20.);
  track_vz_weighted = new TH1F((TString) (desc + "_track_vz_weighted"), "", 80, -20., 20.);


  for (int ibin=0;ibin<nCBins;ibin++){


 jet_pT_hist[ibin] = new TH1F((TString) (desc + "_jet_pT_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 25, 0., 500.); jet_pT_hist[ibin]->Sumw2();
jet_phi_hist[ibin] = new TH1F((TString) (desc + "_jet_phi_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ), "", 36, 0.,TMath::Pi());
jet_phi_hist[ibin]->Sumw2();
jet_eta_hist[ibin] = new TH1F((TString) (desc + "_jet_eta_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ), "", 50, -5., 5.);     jet_eta_hist[ibin]->Sumw2();
jet_corrpT_hist[ibin] = new TH1F((TString) (desc + "_jet_corrpT_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 25, 0., 500.);     jet_corrpT_hist[ibin]->Sumw2();
LeadingJetPt_hist[ibin] = new TH1F((TString) (desc + "_LeadingJetPt_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ), "", 25, 0., 500.);     LeadingJetPt_hist[ibin]->Sumw2();
LeadingJetPhi_hist[ibin] = new TH1F((TString) (desc + "_LeadingJetPhi_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 36, 0.,TMath::Pi());     LeadingJetPhi_hist[ibin]->Sumw2();
LeadingJetEta_hist[ibin] = new TH1F((TString) (desc + "_LeadingJetEta_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 50, -5., 5.);     LeadingJetEta_hist[ibin]->Sumw2();
SubJetPt_hist[ibin] = new TH1F((TString) (desc + "_SubJetPt_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 25, 0., 500.);     SubJetPt_hist[ibin]->Sumw2();
SubJetPhi_hist[ibin] = new TH1F((TString) (desc + "_SubJetPhi_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 36, 0.,TMath::Pi());     SubJetPhi_hist[ibin]->Sumw2();
SubJetEta_hist[ibin] = new TH1F((TString) (desc + "_SubJetEta_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 50, -5., 5.);     SubJetEta_hist[ibin]->Sumw2();

dPhi_leadingjet_hist[ibin] = new TH1F((TString) (desc + "_dPhi_leadingjet_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 30,0,3.14159);     dPhi_leadingjet_hist[ibin]->Sumw2();
dPhi_subleadingjet_hist[ibin] = new TH1F((TString) (desc + "_dPhi_subleadingjet_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 30,0,3.14159);     dPhi_subleadingjet_hist[ibin]->Sumw2();




 JetShapeDiffParticles_bkg_1D_total[ibin] = new TH1F((TString) (desc + "_JetShapeDiffParticles_bkg_1D_total"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_bkg_1D_total[ibin]->Sumw2();     

JetShapeIntegratedParticles_bkg_1D_total[ibin] = new TH1F((TString) (desc + "_JetShapeIntegratedParticles_bkg_1D_total"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ), "", 6,0.05,0.35);     JetShapeIntegratedParticles_bkg_1D_total[ibin]->Sumw2();

SumPt_only_centbin[ibin] = new TH1F((TString) (desc + "_SumPt_only_centbin"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 500, 0., 10.);     SumPt_only_centbin[ibin]->Sumw2();


      Aj[ibin] = new TH1F((TString) (desc + "_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 25,0.,1.25);     Aj[ibin]->Sumw2();


      dPhi_hist[ibin] = new TH1F((TString) (desc + "_dPhi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 30, 0.,TMath::Pi());     dPhi_hist[ibin]->Sumw2();
      dPhi_after_hist[ibin] = new TH1F((TString) (desc + "_dPhi_after_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] ), "", 30,0,3.14159);     dPhi_after_hist[ibin]->Sumw2();


    for (int ibin2=0;ibin2<nPtBins;ibin2++){

     
      ThirdJetPt_hist[ibin][ibin2] = new TH1F((TString) (desc + "_ThirdJetPt_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     ThirdJetPt_hist[ibin][ibin2]->Sumw2();
      ThirdJetPhi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_ThirdJetPhi_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     ThirdJetPhi_hist[ibin][ibin2]->Sumw2();
      ThirdJetEta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_ThirdJetEta_hist_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     ThirdJetEta_hist[ibin][ibin2]->Sumw2();
     
//     dPhi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 60,-3.14159,3.14159);     dPhi_hist[ibin][ibin2]->Sumw2();
   
  //  dPhi_after_hist[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_after_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 60,-3.14159,3.14159);     dPhi_after_hist[ibin][ibin2]->Sumw2();

     dPhi_jet_track[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_jet_track"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 30,0,3.14159);     dPhi_jet_track[ibin][ibin2]->Sumw2();


     dPhi_jet_track_ptweight[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_jet_track_ptweight"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 30,0,3.14159);     dPhi_jet_track_ptweight[ibin][ibin2]->Sumw2();





      // dPhi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 6);     // dPhi_hist[ibin][ibin2]->Sumw2();
      //dPhi_after_hist[ibin][ibin2] = new TH1F((TString) (desc + "_dPhi_after_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 6);     //dPhi_after_hist[ibin][ibin2]->Sumw2();


      all_jets_corrpT[ibin][ibin2] = new TH1F((TString) (desc + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     all_jets_corrpT[ibin][ibin2]->Sumw2();
      all_jets_phi[ibin][ibin2] = new TH1F((TString) (desc + "_all_jets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     all_jets_phi[ibin][ibin2]->Sumw2();
      all_jets_eta[ibin][ibin2] = new TH1F((TString) (desc + "_all_jets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     all_jets_eta[ibin][ibin2]->Sumw2();

      //// leading jet histograms
      only_leadingjets_corrpT[ibin][ibin2] = new TH1F((TString) (desc + "_only_leadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     only_leadingjets_corrpT[ibin][ibin2]->Sumw2();
      only_leadingjets_phi[ibin][ibin2] = new TH1F((TString) (desc + "_only_leadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     only_leadingjets_phi[ibin][ibin2]->Sumw2();
      only_leadingjets_eta[ibin][ibin2] = new TH1F((TString) (desc + "_only_leadingjets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     only_leadingjets_eta[ibin][ibin2]->Sumw2();
      all_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     all_cand_pT_hist[ibin][ibin2]->Sumw2();
      all_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     all_cand_phi_hist[ibin][ibin2]->Sumw2();
      all_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     all_cand_eta_hist[ibin][ibin2]->Sumw2();

      track_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_track_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     track_cand_pT_hist[ibin][ibin2]->Sumw2();
      track_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_track_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     track_cand_phi_hist[ibin][ibin2]->Sumw2();
      track_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_track_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     track_cand_eta_hist[ibin][ibin2]->Sumw2();
      track_cand_pT_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_track_cand_pT_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     track_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();




      NumNeutral[ibin][ibin2] = new TH1F((TString) (desc + "_NumNeutral"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumNeutral[ibin][ibin2]->Sumw2();
      NumPhotons[ibin][ibin2] = new TH1F((TString) (desc + "_NumPhotons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumPhotons[ibin][ibin2]->Sumw2();
      NumAll[ibin][ibin2] = new TH1F((TString) (desc + "_NumAll"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumAll[ibin][ibin2]->Sumw2();
      NumCharged[ibin][ibin2] = new TH1F((TString) (desc + "_NumCharged"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumCharged[ibin][ibin2]->Sumw2();
      NumChargedHadrons[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedHadrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedHadrons[ibin][ibin2]->Sumw2();
      NumChargedParticles[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedParticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedParticles[ibin][ibin2]->Sumw2();
      NumElectrons[ibin][ibin2] = new TH1F((TString) (desc + "_NumElectrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumElectrons[ibin][ibin2]->Sumw2();
      NumMuons[ibin][ibin2] = new TH1F((TString) (desc + "_NumMuons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumMuons[ibin][ibin2]->Sumw2();






      NumNeutral_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumNeutral_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumNeutral_bkg[ibin][ibin2]->Sumw2();
      NumPhotons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumPhotons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumPhotons_bkg[ibin][ibin2]->Sumw2();
      NumAll_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumAll_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumAll_bkg[ibin][ibin2]->Sumw2();
      NumChargedHadrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedHadrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedHadrons_bkg[ibin][ibin2]->Sumw2();
      NumChargedParticles_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedParticles_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedParticles_bkg[ibin][ibin2]->Sumw2();
      NumElectrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumElectrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumElectrons_bkg[ibin][ibin2]->Sumw2();
      NumMuons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumMuons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumMuons_bkg[ibin][ibin2]->Sumw2();



      NumNeutral_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumNeutral_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumNeutral_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      NumPhotons_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumPhotons_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumPhotons_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      NumAll_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumAll_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumAll_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      NumChargedHadrons_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedHadrons_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedHadrons_subleadingjet_bkg[ibin][ibin2]->Sumw2();





      ///===================////////======================//////////////===================//////////=================//////==========

      chargedhadrons_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedhadrons_cand_pT_hist[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     chargedhadrons_cand_phi_hist[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     chargedhadrons_cand_eta_hist[ibin][ibin2]->Sumw2();

      chargedparticles_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedparticles_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedparticles_cand_pT_hist[ibin][ibin2]->Sumw2();
      chargedparticles_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedparticles_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     chargedparticles_cand_phi_hist[ibin][ibin2]->Sumw2();
      chargedparticles_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_chargedparticles_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     chargedparticles_cand_eta_hist[ibin][ibin2]->Sumw2();

      electrons_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_electrons_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     electrons_cand_pT_hist[ibin][ibin2]->Sumw2();
      electrons_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_electrons_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     electrons_cand_phi_hist[ibin][ibin2]->Sumw2();
      electrons_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_electrons_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     electrons_cand_eta_hist[ibin][ibin2]->Sumw2();

      muons_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_muons_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     muons_cand_pT_hist[ibin][ibin2]->Sumw2();
      muons_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_muons_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     muons_cand_phi_hist[ibin][ibin2]->Sumw2();
      muons_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_muons_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     muons_cand_eta_hist[ibin][ibin2]->Sumw2();

      neutral_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     neutral_cand_pT_hist[ibin][ibin2]->Sumw2();
      neutral_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     neutral_cand_phi_hist[ibin][ibin2]->Sumw2();
      neutral_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     neutral_cand_eta_hist[ibin][ibin2]->Sumw2();

      photons_cand_pT_hist[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     photons_cand_pT_hist[ibin][ibin2]->Sumw2();
      photons_cand_phi_hist[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     photons_cand_phi_hist[ibin][ibin2]->Sumw2();
      photons_cand_eta_hist[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     photons_cand_eta_hist[ibin][ibin2]->Sumw2();

      ///===================////////======================//////////////===================//////////=================//////==========
      ///===================////////======================//////////////===================//////////=================//////==========
      ///===================////////======================//////////////===================//////////=================//////==========



      chargedhadrons_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedhadrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();

      chargedparticles_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_chargedparticles_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedparticles_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();

      electrons_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_electrons_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     electrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();

      muons_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_muons_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     muons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();

      neutral_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     neutral_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();

      photons_cand_pT_hist_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_pT_hist_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     photons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();


      /// subleading jet

      chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_pT_hist_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_pT_hist_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();
      photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_pT_hist_subleadingjet_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();

      ///===================////////======================//////////////===================//////////=================//////==========



      //SumJetPtFraction_hist[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     //SumJetPtFraction_hist[ibin][ibin2]->Sumw2();


      SumJetPtFraction_hist[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100,0.,500);     SumJetPtFraction_hist[ibin][ibin2]->Sumw2();





      JetShapeIntegratedParticles[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeIntegratedParticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3,0.,2.);     JetShapeIntegratedParticles[ibin][ibin2]->Sumw2();

      JetShapeIntegratedParticles_bkgsub[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeIntegratedParticles_bkgsub"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3,0.,2.); JetShapeIntegratedParticles_bkgsub[ibin][ibin2]->Sumw2();




      JetShapeDiffParticles[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeDiffParticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3,0.,100.);     JetShapeDiffParticles[ibin][ibin2]->Sumw2();


      JetShapeDiffParticles_bkgsub[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeDiffParticles_bkgsub"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3,0.,100.);     JetShapeDiffParticles_bkgsub[ibin][ibin2]->Sumw2();

      JetShapeIntegratedParticles_bkg[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeIntegratedParticles_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.05,0.35,0.,2.);     JetShapeIntegratedParticles_bkg[ibin][ibin2]->Sumw2();
      JetShapeDiffParticles_bkg[ibin][ibin2] = new TProfile((TString) (desc + "_JetShapeDiffParticles_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3,0.,100.);     JetShapeDiffParticles_bkg[ibin][ibin2]->Sumw2();


      radius_hist[ibin][ibin2] = new TH1F((TString) (desc + "_radius_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, 0., 5.);     radius_hist[ibin][ibin2]->Sumw2();

      JetShapeDiffParticles_bkg_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticles_bkg_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Sumw2();     
     
      JetShapeDiffParticles_bkg_1D_EtaRef[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticles_bkg_1D_EtaRef"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_bkg_1D_EtaRef[ibin][ibin2]->Sumw2();     
     
     
      JetShapeDiffParticles_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticles_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_1D[ibin][ibin2]->Sumw2();


      JetShapeIntegratedParticles_bkg_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeIntegratedParticles_bkg_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.05,0.35);     JetShapeIntegratedParticles_bkg_1D[ibin][ibin2]->Sumw2();
      JetShapeIntegratedParticles_1D[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeIntegratedParticles_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.05,0.35);     JetShapeIntegratedParticles_1D[ibin][ibin2]->Sumw2();



      SumPt_in_a_radius[ibin][ibin2] = new TH1F((TString) (desc + "_SumPt_in_a_radius"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     SumPt_in_a_radius[ibin][ibin2]->Sumw2();

      multiplicity_in_a_radius[ibin][ibin2] = new TH1F((TString) (desc + "_multiplicity_in_a_radius"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     multiplicity_in_a_radius[ibin][ibin2]->Sumw2();




      //// subleading jet histograms
      only_subleadingjets_corrpT[ibin][ibin2] = new TH1F((TString) (desc + "_only_subleadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     only_subleadingjets_corrpT[ibin][ibin2]->Sumw2();
      only_subleadingjets_phi[ibin][ibin2] = new TH1F((TString) (desc + "_only_subleadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     only_subleadingjets_phi[ibin][ibin2]->Sumw2();
      only_subleadingjets_eta[ibin][ibin2] = new TH1F((TString) (desc + "_only_subleadingjets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     only_subleadingjets_eta[ibin][ibin2]->Sumw2();
      all_cand_pT_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_pT_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     all_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      all_cand_phi_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_phi_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     all_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      all_cand_eta_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_all_cand_eta_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     all_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
      neutral_cand_pT_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_pT_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     neutral_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      neutral_cand_phi_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_phi_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     neutral_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      neutral_cand_eta_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_neutral_cand_eta_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     neutral_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
      photons_cand_pT_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_pT_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     photons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      photons_cand_phi_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_phi_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     photons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      photons_cand_eta_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_photons_cand_eta_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     photons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();

      chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_pT_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_phi_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, 0.,TMath::Pi());     chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
      chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_chargedhadrons_cand_eta_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();


      NumNeutral_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_NumNeutral_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumNeutral_subleadingjet[ibin][ibin2]->Sumw2();
      NumPhotons_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_NumPhotons_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumPhotons_subleadingjet[ibin][ibin2]->Sumw2();
      NumAll_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_NumAll_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumAll_subleadingjet[ibin][ibin2]->Sumw2();

      NumChargedHadrons_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_NumChargedHadrons_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 100.);     NumChargedHadrons_subleadingjet[ibin][ibin2]->Sumw2();



      //SumJetPtFraction_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100,0.,300);     //SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Sumw2();


      SumJetPtFraction_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_photons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_neutrals"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_chargedhadrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2]->Sumw2();

      SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_photons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_neutrals"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_charged"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2]->Sumw2();

      SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_chargedhadrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_chargedparticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_electrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_muons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2]->Sumw2();




      /// bkg -- sumpt fraction plots 

      SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_photons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_neutrals_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_chargedparticles_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_electrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_leadingjet_muons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2]->Sumw2();



      SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_photons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_neutrals_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_electrons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2]->Sumw2();
      SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2] = new TH1F((TString) (desc + "_SumJetPtFraction_hist_subleadingjet_muons_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50,0.,1.2);     SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2]->Sumw2();



      radius_hist_subleadingjet[ibin][ibin2] = new TH1F((TString) (desc + "_radius_hist_subleadingjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, 0., 5.);     radius_hist_subleadingjet[ibin][ibin2]->Sumw2();
      Centrality_hist[ibin][ibin2] = new TH1F((TString) (desc + "_Centrality_hist_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 200, 0., 200.);     Centrality_hist[ibin][ibin2]->Sumw2();
      //Aj[ibin][ibin2] = new TH1F((TString) (desc + "_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 25,0.,1.25);     Aj[ibin][ibin2]->Sumw2();

      dN_tracks[ibin][ibin2] = new TH1F((TString) (desc + "_dN_tracks"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_tracks[ibin][ibin2]->Sumw2();
      dN_chargedhadrons[ibin][ibin2] = new TH1F((TString) (desc + "_dN_chargedhadrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_chargedhadrons[ibin][ibin2]->Sumw2();
      dN_chargedparticles[ibin][ibin2] = new TH1F((TString) (desc + "_dN_chargedparticles"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_chargedparticles[ibin][ibin2]->Sumw2();
      dN_electrons[ibin][ibin2] = new TH1F((TString) (desc + "_dN_electrons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_electrons[ibin][ibin2]->Sumw2();
      dN_muons[ibin][ibin2] = new TH1F((TString) (desc + "_dN_muons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_muons[ibin][ibin2]->Sumw2();
      dN_neutrals[ibin][ibin2] = new TH1F((TString) (desc + "_dN_neutrals"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_neutrals[ibin][ibin2]->Sumw2();
      dN_photons[ibin][ibin2] = new TH1F((TString) (desc + "_dN_photons"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 300,0.,10.);     dN_photons[ibin][ibin2]->Sumw2();

      JetPt_fraction_hist[ibin][ibin2] = new TH1F((TString) (desc + "_JetPt_fraction_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100,0.,5);     JetPt_fraction_hist[ibin][ibin2]->Sumw2();


      track_cand_pT_hist_weighted[ibin][ibin2] = new TH1F((TString) (desc + "_track_cand_pT_hist_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     track_cand_pT_hist_weighted[ibin][ibin2]->Sumw2();
     


      JetShapeDiffParticles_bkg_1DEtaRef[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticles_bkg_1DEtaRef"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_bkg_1DEtaRef[ibin][ibin2]->Sumw2();     


      JetShapeDiffParticles_bkg_1DEtaRefSig[ibin][ibin2] = new TH1F((TString) (desc + "_JetShapeDiffParticles_bkg_1DEtaRefSig"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_bkg_1DEtaRefSig[ibin][ibin2]->Sumw2();     




      SumPt_ratio[ibin][ibin2] = new TH1F((TString) (desc + "_SumPt_ratio"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 500,0.,10);     SumPt_ratio[ibin][ibin2]->Sumw2();

      SumPt_only[ibin][ibin2] = new TH1F((TString) (desc + "_SumPt_only"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 250,0.0,250);     SumPt_only[ibin][ibin2]->Sumw2();

      SumPt_bkg_only[ibin][ibin2] = new TH1F((TString) (desc + "_SumPt_bkg_only"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 250,0.0,250);     SumPt_bkg_only[ibin][ibin2]->Sumw2();




EtaRef_bkg_pt[ibin][ibin2] = new TH1F((TString) (desc + "_EtaRef_bkg_pt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk); EtaRef_bkg_pt[ibin][ibin2]->Sumw2();

EtaRef_bkg_pt_weighted[ibin][ibin2] = new TH1F((TString) (desc + "_EtaRef_bkg_pt_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk); EtaRef_bkg_pt_weighted[ibin][ibin2]->Sumw2();




        
   
EMix_bkg_pt[ibin][ibin2] = new TH1D((TString) (desc + "_EMix_bkg_pt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] +
                                                       "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk); EMix_bkg_pt[ibin][ibin2]->Sumw2();

EMix_bkg_phi[ibin][ibin2] = new TH1D((TString) (desc + "_EMix_bkg_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi(),TMath::Pi());     EMix_bkg_phi[ibin][ibin2]->Sumw2();
        
     
EMix_bkg_eta[ibin][ibin2] = new TH1D((TString) (desc + "_EMix_bkg_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     EMix_bkg_eta[ibin][ibin2]->Sumw2();
        
     
        
        EMix_bkg_pt_weighted[ibin][ibin2] = new TH1D((TString) (desc + "_EMix_bkg_pt_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk); EMix_bkg_pt_weighted[ibin][ibin2]->Sumw2();

        
        EMix_bkg_phi_weighted[ibin][ibin2] = new TH1D((TString) (desc + "_EMix_bkg_phi_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi(),TMath::Pi());     EMix_bkg_phi_weighted[ibin][ibin2]->Sumw2();
        
        
        EMix_bkg_eta_weighted[ibin][ibin2] = new TH1D((TString) (desc + "_EMix_bkg_eta_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);     EMix_bkg_eta_weighted[ibin][ibin2]->Sumw2();
        

        
        

        
		track_bkg_pT_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_bkg_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     track_bkg_pT_hist[ibin][ibin2]->Sumw2();

        
        track_bkg_phi_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_bkg_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi(),TMath::Pi());     track_bkg_phi_hist[ibin][ibin2]->Sumw2();
        
        
        track_bkg_eta_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_bkg_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);  track_bkg_eta_hist[ibin][ibin2]->Sumw2();
        
    
        
        
        track_bkg_pT_hist_weighted[ibin][ibin2] = new TH1D((TString) (desc + "_track_bkg_pT_hist_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);     track_bkg_pT_hist_weighted[ibin][ibin2]->Sumw2();

        track_bkg_phi_hist_weighted[ibin][ibin2] = new TH1D((TString) (desc + "_track_bkg_phi_hist_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi(),TMath::Pi());     track_bkg_phi_hist_weighted[ibin][ibin2]->Sumw2();


        track_bkg_eta_hist_weighted[ibin][ibin2] = new TH1D((TString) (desc + "_track_bkg_eta_hist_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);  track_bkg_eta_hist_weighted[ibin][ibin2]->Sumw2();


        track_gen_bkg_pT_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_bkg_pT_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", nPtBintrk,ptBintrk);   track_gen_bkg_pT_hist[ibin][ibin2]->Sumw2();

        track_gen_bkg_eta_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_bkg_eta_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, -5., 5.);   track_gen_bkg_eta_hist[ibin][ibin2]->Sumw2();

        track_gen_bkg_phi_hist[ibin][ibin2] = new TH1D((TString) (desc + "_track_gen_bkg_phi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi(),TMath::Pi());   track_gen_bkg_phi_hist[ibin][ibin2]->Sumw2();

        emix_jetpt[ibin][ibin2] = new TH1F((TString) (desc + "_emix_jetpt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, 0., 500.);     emix_jetpt[ibin][ibin2]->Sumw2(); // pt1 bin loop

        emix_jetpt_test[ibin][ibin2] = new TH1F((TString) (desc + "_emix_jetpt_test"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, 0., 500.);     emix_jetpt_test[ibin][ibin2]->Sumw2(); // pt1 bin loop


        for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

          hJetGenTrackBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackBackground[ibin][ibin2][ibin3]->Sumw2();

          hJetGenTrackBackground_pythia[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackBackground_pythia"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Sumw2();


          hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetGenTrackBackground_hydjet"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Sumw2();

          hJetTrackSignalBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());     hJetTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();

          hJetTrackBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());  hJetTrackBackground[ibin][ibin2][ibin3]->Sumw2();

          hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,75,-TMath::Pi(),2*TMath::Pi());  hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();

          Bkg_phi[ibin][ibin2][ibin3] = new TH1D((TString) (desc + "_Bkg_phi_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 72, -TMath::Pi(),TMath::Pi());     Bkg_phi[ibin][ibin2][ibin3]->Sumw2();

          Bkg_eta[ibin][ibin2][ibin3] = new TH1D((TString) (desc + "_Bkg_eta_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 50, -5, 5);     Bkg_eta[ibin][ibin2][ibin3]->Sumw2();

          Jet_phi[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_Jet_phi_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 72, -TMath::Pi(),TMath::Pi());     Jet_phi[ibin][ibin2][ibin3]->Sumw2();

        Bkg_phi_weighted[ibin][ibin2][ibin3] = new TH1D((TString) (desc + "_Bkg_phi_weighted_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 72, -TMath::Pi(),TMath::Pi());     Bkg_phi_weighted[ibin][ibin2][ibin3]->Sumw2();
        
        Bkg_eta_weighted[ibin][ibin2][ibin3] = new TH1D((TString) (desc + "_Bkg_eta_weighted_" + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 50, -5, 5);     Bkg_eta_weighted[ibin][ibin2][ibin3]->Sumw2();

        GenTrkEta_bkg[ibin][ibin2][ibin3] = new TH1D((TString) (desc + "_GenTrkEta_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 50, -5., 5.);     GenTrkEta_bkg[ibin][ibin2][ibin3]->Sumw2();

        GenTrkPhi_bkg[ibin][ibin2][ibin3] = new TH1D((TString) (desc + "_GenTrkPhi_bkg"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 72, -TMath::Pi(), TMath::Pi());     GenTrkPhi_bkg[ibin][ibin2][ibin3]->Sumw2();        

 }  /// ibin3

 } // ibin2  pt bin loop leading jets


  for (int ibin3=0;ibin3<nPtBins_sub;ibin3++){  /// pti3 subjet plots

      JetShapeDiffParticles_subleadingjet_bkg_1DEtaRef[ibin][ibin3] = new TH1F((TString) (desc + "_JetShapeDiffParticles_subleadingjet_bkg_1DEtaRef"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_sub_strs[ibin3] + "_" + PtBin_sub_strs[ibin3+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_subleadingjet_bkg_1DEtaRef[ibin][ibin3]->Sumw2();     

      JetShapeIntegratedParticles_subleadingjet_1D[ibin][ibin3] = new TH1F((TString) (desc + "_JetShapeIntegratedParticles_subleadingjet_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_sub_strs[ibin3] + "_" + PtBin_sub_strs[ibin3+1]), "", 6,0.05,0.35);     JetShapeIntegratedParticles_subleadingjet_1D[ibin][ibin3]->Sumw2();

      JetShapeDiffParticles_subleadingjet_bkg_1D[ibin][ibin3] = new TH1F((TString) (desc + "_JetShapeDiffParticles_subleadingjet_bkg_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_sub_strs[ibin3] + "_" + PtBin_sub_strs[ibin3+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_subleadingjet_bkg_1D[ibin][ibin3]->Sumw2();    

      JetShapeDiffParticles_subleadingjet_1D[ibin][ibin3] = new TH1F((TString) (desc + "_JetShapeDiffParticles_subleadingjet_1D"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_sub_strs[ibin3] + "_" + PtBin_sub_strs[ibin3+1]), "", 6,0.0,0.3);     JetShapeDiffParticles_subleadingjet_1D[ibin][ibin3]->Sumw2();   
  }


  } // centrality bin loop
} // hist class loop




void hist_class::NormalizeByCounts()
{
  /*for(int centi = 0; centi < nCBins; centi++) {
    for(int pti = 0; pti < nPtBins; pti++) {
      float int_just_before_scale = JetShapeDiffParticles_1D[centi][pti]->Integral(0, JetShapeDiffParticles_1D[centi][pti]->GetNbinsX()+1);
      // std::cout << " cent : " << centi << "  pti: " << pti << "   integral_just_before_scale " << int_just_before_scale << std::endl;

      if( sig_count[centi][pti] < 0.5e-18 ) continue;
      JetShapeDiffParticles_1D[centi][pti]->Scale(1. / sig_count[centi][pti]);

      //Multiplicity_vs_radius[centi][pti]->Scale(1. / sig_count[centi][pti]);

     // std::cout << " cent : " << centi << "    pti: " <<  pti << "    sig_count:  " << sig_count[centi][pti] << std::endl;

      //float int_just_after_scale = JetShapeDiffParticles_1D[centi][pti]->Integral(0, JetShapeDiffParticles_1D[centi][pti]->GetNbinsX()+1);
      //std::cout << " cent : " << centi << "  pti: " << pti << "   integral_just_after_scale " << int_just_after_scale << std::endl;
     // std::cout << " cent : " << centi << "    pti: " <<  pti << "    sig_bkg:  " << bkg_count[centi][pti] << std::endl;

      if( bkg_count[centi][pti] < 0.5e-18  ) continue;
      JetShapeDiffParticles_bkg_1D[centi][pti]->Scale(1. / bkg_count[centi][pti]);
      //  Multiplicity_vs_radius_bkg[centi][pti]->Scale(1. / bkg_count[centi][pti]);

      //std::cout << "centi: " << centi << ", pti: " << pti << ", scale by sig_count: " << sig_count[centi][pti] << ", bkg_count: " << bkg_count[centi][pti] << "\n";

    }
  }*/
}







void hist_class::AddHists(hist_class *more_hists, float wt)
{



  for(int centi = 0; centi < nCBins; centi++) {

    for(int pti = 0; pti < nPtBins; pti++) {
      sig_count[centi][pti] += more_hists->sig_count[centi][pti]*wt;
      bkg_count[centi][pti] += more_hists->bkg_count[centi][pti]*wt;
    }
    
    
    /*for(int pti_sub = 0; pti_sub < nPtBins_sub; pti_sub++) {

      sig_count_sub[centi][pti_sub] += more_hists->sig_count_sub[centi][pti_sub]*wt;
      bkg_count_sub[centi][pti_sub] += more_hists->bkg_count_sub[centi][pti_sub]*wt;

      }*/

  }



  NumberOfMatches->Sumw2();   more_hists->NumberOfMatches->Sumw2();
  NumberOfMatches->Add(more_hists->NumberOfMatches, wt);


  NEvents->Sumw2();   more_hists->NEvents->Sumw2();
  NEvents->Add(more_hists->NEvents, wt);

    NEvents_test->Sumw2();   more_hists->NEvents_test->Sumw2();
    NEvents_test->Add(more_hists->NEvents_test, wt);

    NEvents_after_noise->Add(more_hists->NEvents_after_noise, wt);
    NEvents_after_spike->Add(more_hists->NEvents_after_spike, wt);

    NEvents_after_dphi->Sumw2();   more_hists->NEvents_after_dphi->Sumw2();
    NEvents_after_dphi->Add(more_hists->NEvents_after_dphi, wt);

    NEvents_before_dphi->Sumw2();   more_hists->NEvents_before_dphi->Sumw2();
    NEvents_before_dphi->Add(more_hists->NEvents_before_dphi, wt);

    NEvents_dijets->Sumw2();   more_hists->NEvents_dijets->Sumw2();
    NEvents_dijets->Add(more_hists->NEvents_dijets, wt);

    NEvents_after_trigger->Sumw2();   more_hists->NEvents_after_trigger->Sumw2();
    NEvents_after_trigger->Add(more_hists->NEvents_after_trigger, wt);

    AllJetPt_raw_hist->Sumw2();   more_hists->AllJetPt_raw_hist->Sumw2();
    AllJetPt_raw_hist->Add(more_hists->AllJetPt_raw_hist, wt);
    AllJetPt_hist->Sumw2();   more_hists->AllJetPt_hist->Sumw2();
    AllJetPt_hist->Add(more_hists->AllJetPt_hist, wt);
    AllJetPhi_hist->Sumw2();   more_hists->AllJetPhi_hist->Sumw2();
    AllJetPhi_hist->Add(more_hists->AllJetPhi_hist, wt);
    AllJetEta_hist->Sumw2();   more_hists->AllJetEta_hist->Sumw2();
    AllJetEta_hist->Add(more_hists->AllJetEta_hist, wt);
 
    First_AllJetPt_hist->Sumw2();   more_hists->First_AllJetPt_hist->Sumw2();
    First_AllJetPt_hist->Add(more_hists->First_AllJetPt_hist, wt);

    First_AllJetPhi_hist->Sumw2();   more_hists->First_AllJetPhi_hist->Sumw2();
    First_AllJetPhi_hist->Add(more_hists->First_AllJetPhi_hist, wt);

    First_AllJetEta_hist->Sumw2();   more_hists->First_AllJetEta_hist->Sumw2();
    First_AllJetEta_hist->Add(more_hists->First_AllJetEta_hist, wt);

    Sub_AllJetPt_hist->Sumw2();   more_hists->Sub_AllJetPt_hist->Sumw2();
    Sub_AllJetPt_hist->Add(more_hists->Sub_AllJetPt_hist, wt);

    Sub_AllJetPhi_hist->Sumw2();   more_hists->Sub_AllJetPhi_hist->Sumw2();
    Sub_AllJetPhi_hist->Add(more_hists->Sub_AllJetPhi_hist, wt);

    Sub_AllJetEta_hist->Sumw2();   more_hists->Sub_AllJetEta_hist->Sumw2();
    Sub_AllJetEta_hist->Add(more_hists->Sub_AllJetEta_hist, wt);

    JetPt_fraction->Sumw2();   more_hists->JetPt_fraction->Sumw2();
    JetPt_fraction->Add(more_hists->JetPt_fraction, wt);

    Centrality->Sumw2();   more_hists->Centrality->Sumw2();
    Centrality->Add(more_hists->Centrality, wt);


  track_vz->Sumw2();   more_hists->track_vz->Sumw2();
  track_vz->Add(more_hists->track_vz, wt);

  track_vz_weighted->Sumw2();   more_hists->track_vz_weighted->Sumw2();
  track_vz_weighted->Add(more_hists->track_vz_weighted, wt);


  for (int ibin=0;ibin<nCBins;ibin++){

    jet_pT_hist[ibin]->Sumw2();   more_hists->jet_pT_hist[ibin]->Sumw2();
    jet_pT_hist[ibin]->Add(more_hists->jet_pT_hist[ibin], wt);
    jet_phi_hist[ibin]->Sumw2();   more_hists->jet_phi_hist[ibin]->Sumw2();
    jet_phi_hist[ibin]->Add(more_hists->jet_phi_hist[ibin], wt);
    jet_eta_hist[ibin]->Sumw2();   more_hists->jet_eta_hist[ibin]->Sumw2();
    jet_eta_hist[ibin]->Add(more_hists->jet_eta_hist[ibin], wt);
    jet_corrpT_hist[ibin]->Sumw2();   more_hists->jet_corrpT_hist[ibin]->Sumw2();
    jet_corrpT_hist[ibin]->Add(more_hists->jet_corrpT_hist[ibin], wt);
    LeadingJetPt_hist[ibin]->Sumw2();   more_hists->LeadingJetPt_hist[ibin]->Sumw2();
    LeadingJetPt_hist[ibin]->Add(more_hists->LeadingJetPt_hist[ibin], wt);
    LeadingJetPhi_hist[ibin]->Sumw2();   more_hists->LeadingJetPhi_hist[ibin]->Sumw2();
    LeadingJetPhi_hist[ibin]->Add(more_hists->LeadingJetPhi_hist[ibin], wt);
    LeadingJetEta_hist[ibin]->Sumw2();   more_hists->LeadingJetEta_hist[ibin]->Sumw2();
    LeadingJetEta_hist[ibin]->Add(more_hists->LeadingJetEta_hist[ibin], wt);
    SubJetPt_hist[ibin]->Sumw2();   more_hists->SubJetPt_hist[ibin]->Sumw2();
    SubJetPt_hist[ibin]->Add(more_hists->SubJetPt_hist[ibin], wt);
    SubJetPhi_hist[ibin]->Sumw2();   more_hists->SubJetPhi_hist[ibin]->Sumw2();
    SubJetPhi_hist[ibin]->Add(more_hists->SubJetPhi_hist[ibin], wt);
    SubJetEta_hist[ibin]->Sumw2();   more_hists->SubJetEta_hist[ibin]->Sumw2();
    SubJetEta_hist[ibin]->Add(more_hists->SubJetEta_hist[ibin], wt);


    dPhi_leadingjet_hist[ibin]->Sumw2();   more_hists->dPhi_leadingjet_hist[ibin]->Sumw2();
    dPhi_leadingjet_hist[ibin]->Add(more_hists->dPhi_leadingjet_hist[ibin], wt);
    dPhi_subleadingjet_hist[ibin]->Sumw2();   more_hists->dPhi_subleadingjet_hist[ibin]->Sumw2();
    dPhi_subleadingjet_hist[ibin]->Add(more_hists->dPhi_subleadingjet_hist[ibin], wt);


    JetShapeDiffParticles_bkg_1D_total[ibin]->Sumw2();   more_hists->JetShapeDiffParticles_bkg_1D_total[ibin]->Sumw2();
    JetShapeDiffParticles_bkg_1D_total[ibin]->Add(more_hists->JetShapeDiffParticles_bkg_1D_total[ibin], wt);

    JetShapeIntegratedParticles_bkg_1D_total[ibin]->Sumw2();   more_hists->JetShapeIntegratedParticles_bkg_1D_total[ibin]->Sumw2();
    JetShapeIntegratedParticles_bkg_1D_total[ibin]->Add(more_hists->JetShapeIntegratedParticles_bkg_1D_total[ibin], wt);

    SumPt_only_centbin[ibin]->Sumw2();   more_hists->SumPt_only_centbin[ibin]->Sumw2();
    SumPt_only_centbin[ibin]->Add(more_hists->SumPt_only_centbin[ibin], wt);


    dPhi_hist[ibin]->Sumw2();   more_hists->dPhi_hist[ibin]->Sumw2();
    dPhi_hist[ibin]->Add(more_hists->dPhi_hist[ibin], wt);

    dPhi_after_hist[ibin]->Sumw2();   more_hists->dPhi_after_hist[ibin]->Sumw2();
    dPhi_after_hist[ibin]->Add(more_hists->dPhi_after_hist[ibin], wt);

    Aj[ibin]->Sumw2();   more_hists->Aj[ibin]->Sumw2();
    Aj[ibin]->Add(more_hists->Aj[ibin], wt);




      for (int ibin2=0;ibin2<nPtBins;ibin2++){

    
    ThirdJetPt_hist[ibin][ibin2]->Sumw2();   more_hists->ThirdJetPt_hist[ibin][ibin2]->Sumw2();
    ThirdJetPt_hist[ibin][ibin2]->Add(more_hists->ThirdJetPt_hist[ibin][ibin2], wt);
    ThirdJetPhi_hist[ibin][ibin2]->Sumw2();   more_hists->ThirdJetPhi_hist[ibin][ibin2]->Sumw2();
    ThirdJetPhi_hist[ibin][ibin2]->Add(more_hists->ThirdJetPhi_hist[ibin][ibin2], wt);
    ThirdJetEta_hist[ibin][ibin2]->Sumw2();   more_hists->ThirdJetEta_hist[ibin][ibin2]->Sumw2();
    ThirdJetEta_hist[ibin][ibin2]->Add(more_hists->ThirdJetEta_hist[ibin][ibin2], wt);
    //dPhi_hist[ibin][ibin2]->Sumw2();   more_hists->dPhi_hist[ibin][ibin2]->Sumw2();
    //dPhi_hist[ibin][ibin2]->Add(more_hists->dPhi_hist[ibin][ibin2], wt);
    //dPhi_after_hist[ibin][ibin2]->Sumw2();   more_hists->dPhi_after_hist[ibin][ibin2]->Sumw2();
    //dPhi_after_hist[ibin][ibin2]->Add(more_hists->dPhi_after_hist[ibin][ibin2], wt);



    all_jets_corrpT[ibin][ibin2]->Sumw2();   more_hists->all_jets_corrpT[ibin][ibin2]->Sumw2();
    all_jets_corrpT[ibin][ibin2]->Add(more_hists->all_jets_corrpT[ibin][ibin2], wt);
    all_jets_phi[ibin][ibin2]->Sumw2();   more_hists->all_jets_phi[ibin][ibin2]->Sumw2();
    all_jets_phi[ibin][ibin2]->Add(more_hists->all_jets_phi[ibin][ibin2], wt);
    all_jets_eta[ibin][ibin2]->Sumw2();   more_hists->all_jets_eta[ibin][ibin2]->Sumw2();
    all_jets_eta[ibin][ibin2]->Add(more_hists->all_jets_eta[ibin][ibin2], wt);

    //// leading jet histograms
    only_leadingjets_corrpT[ibin][ibin2]->Sumw2();   more_hists->only_leadingjets_corrpT[ibin][ibin2]->Sumw2();
    only_leadingjets_corrpT[ibin][ibin2]->Add(more_hists->only_leadingjets_corrpT[ibin][ibin2], wt);
    only_leadingjets_phi[ibin][ibin2]->Sumw2();   more_hists->only_leadingjets_phi[ibin][ibin2]->Sumw2();
    only_leadingjets_phi[ibin][ibin2]->Add(more_hists->only_leadingjets_phi[ibin][ibin2], wt);
    only_leadingjets_eta[ibin][ibin2]->Sumw2();   more_hists->only_leadingjets_eta[ibin][ibin2]->Sumw2();
    only_leadingjets_eta[ibin][ibin2]->Add(more_hists->only_leadingjets_eta[ibin][ibin2], wt);
    all_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->all_cand_pT_hist[ibin][ibin2]->Sumw2();
    all_cand_pT_hist[ibin][ibin2]->Add(more_hists->all_cand_pT_hist[ibin][ibin2], wt);
    all_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->all_cand_phi_hist[ibin][ibin2]->Sumw2();
    all_cand_phi_hist[ibin][ibin2]->Add(more_hists->all_cand_phi_hist[ibin][ibin2], wt);
    all_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->all_cand_eta_hist[ibin][ibin2]->Sumw2();
    all_cand_eta_hist[ibin][ibin2]->Add(more_hists->all_cand_eta_hist[ibin][ibin2], wt);

    track_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->track_cand_pT_hist[ibin][ibin2]->Sumw2();
    track_cand_pT_hist[ibin][ibin2]->Add(more_hists->track_cand_pT_hist[ibin][ibin2], wt);
    track_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->track_cand_phi_hist[ibin][ibin2]->Sumw2();
    track_cand_phi_hist[ibin][ibin2]->Add(more_hists->track_cand_phi_hist[ibin][ibin2], wt);
    track_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->track_cand_eta_hist[ibin][ibin2]->Sumw2();
    track_cand_eta_hist[ibin][ibin2]->Add(more_hists->track_cand_eta_hist[ibin][ibin2], wt);


    track_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->track_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
    track_cand_pT_hist_subleadingjet[ibin][ibin2]->Add(more_hists->track_cand_pT_hist_subleadingjet[ibin][ibin2], wt);

    neutral_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_pT_hist[ibin][ibin2]->Sumw2();
    neutral_cand_pT_hist[ibin][ibin2]->Add(more_hists->neutral_cand_pT_hist[ibin][ibin2], wt);
    neutral_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_phi_hist[ibin][ibin2]->Sumw2();
    neutral_cand_phi_hist[ibin][ibin2]->Add(more_hists->neutral_cand_phi_hist[ibin][ibin2], wt);
    neutral_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_eta_hist[ibin][ibin2]->Sumw2();
    neutral_cand_eta_hist[ibin][ibin2]->Add(more_hists->neutral_cand_eta_hist[ibin][ibin2], wt);
    photons_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->photons_cand_pT_hist[ibin][ibin2]->Sumw2();
    photons_cand_pT_hist[ibin][ibin2]->Add(more_hists->photons_cand_pT_hist[ibin][ibin2], wt);
    photons_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->photons_cand_phi_hist[ibin][ibin2]->Sumw2();
    photons_cand_phi_hist[ibin][ibin2]->Add(more_hists->photons_cand_phi_hist[ibin][ibin2], wt);
    photons_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->photons_cand_eta_hist[ibin][ibin2]->Sumw2();
    photons_cand_eta_hist[ibin][ibin2]->Add(more_hists->photons_cand_eta_hist[ibin][ibin2], wt);
   

    NumNeutral[ibin][ibin2]->Sumw2();   more_hists->NumNeutral[ibin][ibin2]->Sumw2();
    NumNeutral[ibin][ibin2]->Add(more_hists->NumNeutral[ibin][ibin2], wt);
    NumPhotons[ibin][ibin2]->Sumw2();   more_hists->NumPhotons[ibin][ibin2]->Sumw2();
    NumPhotons[ibin][ibin2]->Add(more_hists->NumPhotons[ibin][ibin2], wt);
    NumAll[ibin][ibin2]->Sumw2();   more_hists->NumAll[ibin][ibin2]->Sumw2();
    NumAll[ibin][ibin2]->Add(more_hists->NumAll[ibin][ibin2], wt);
    NumCharged[ibin][ibin2]->Sumw2();   more_hists->NumCharged[ibin][ibin2]->Sumw2();
    NumCharged[ibin][ibin2]->Add(more_hists->NumCharged[ibin][ibin2], wt);
    NumChargedHadrons[ibin][ibin2]->Sumw2();   more_hists->NumChargedHadrons[ibin][ibin2]->Sumw2();
    NumChargedHadrons[ibin][ibin2]->Add(more_hists->NumChargedHadrons[ibin][ibin2], wt);
    NumChargedParticles[ibin][ibin2]->Sumw2();   more_hists->NumChargedParticles[ibin][ibin2]->Sumw2();
    NumChargedParticles[ibin][ibin2]->Add(more_hists->NumChargedParticles[ibin][ibin2], wt);
    NumElectrons[ibin][ibin2]->Sumw2();   more_hists->NumElectrons[ibin][ibin2]->Sumw2();
    NumElectrons[ibin][ibin2]->Add(more_hists->NumElectrons[ibin][ibin2], wt);
    NumMuons[ibin][ibin2]->Sumw2();   more_hists->NumMuons[ibin][ibin2]->Sumw2();
    NumMuons[ibin][ibin2]->Add(more_hists->NumMuons[ibin][ibin2], wt);


    NumNeutral_bkg[ibin][ibin2]->Sumw2();   more_hists->NumNeutral_bkg[ibin][ibin2]->Sumw2();
    NumNeutral_bkg[ibin][ibin2]->Add(more_hists->NumNeutral_bkg[ibin][ibin2], wt);
    NumPhotons_bkg[ibin][ibin2]->Sumw2();   more_hists->NumPhotons_bkg[ibin][ibin2]->Sumw2();
    NumPhotons_bkg[ibin][ibin2]->Add(more_hists->NumPhotons_bkg[ibin][ibin2], wt);
    NumAll_bkg[ibin][ibin2]->Sumw2();   more_hists->NumAll_bkg[ibin][ibin2]->Sumw2();
    NumAll_bkg[ibin][ibin2]->Add(more_hists->NumAll_bkg[ibin][ibin2], wt);
    NumChargedHadrons_bkg[ibin][ibin2]->Sumw2();   more_hists->NumChargedHadrons_bkg[ibin][ibin2]->Sumw2();
    NumChargedHadrons_bkg[ibin][ibin2]->Add(more_hists->NumChargedHadrons_bkg[ibin][ibin2], wt);
    NumChargedParticles_bkg[ibin][ibin2]->Sumw2();   more_hists->NumChargedParticles_bkg[ibin][ibin2]->Sumw2();
    NumChargedParticles_bkg[ibin][ibin2]->Add(more_hists->NumChargedParticles_bkg[ibin][ibin2], wt);
    NumElectrons_bkg[ibin][ibin2]->Sumw2();   more_hists->NumElectrons_bkg[ibin][ibin2]->Sumw2();
    NumElectrons_bkg[ibin][ibin2]->Add(more_hists->NumElectrons_bkg[ibin][ibin2], wt);
    NumMuons_bkg[ibin][ibin2]->Sumw2();   more_hists->NumMuons_bkg[ibin][ibin2]->Sumw2();
    NumMuons_bkg[ibin][ibin2]->Add(more_hists->NumMuons_bkg[ibin][ibin2], wt);


    chargedhadrons_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_pT_hist[ibin][ibin2]->Sumw2();
    chargedhadrons_cand_pT_hist[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_pT_hist[ibin][ibin2], wt);
    chargedhadrons_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_phi_hist[ibin][ibin2]->Sumw2();
    chargedhadrons_cand_phi_hist[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_phi_hist[ibin][ibin2], wt);
    chargedhadrons_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_eta_hist[ibin][ibin2]->Sumw2();
    chargedhadrons_cand_eta_hist[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_eta_hist[ibin][ibin2], wt);

    chargedparticles_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->chargedparticles_cand_pT_hist[ibin][ibin2]->Sumw2();
    chargedparticles_cand_pT_hist[ibin][ibin2]->Add(more_hists->chargedparticles_cand_pT_hist[ibin][ibin2], wt);
    chargedparticles_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->chargedparticles_cand_phi_hist[ibin][ibin2]->Sumw2();
    chargedparticles_cand_phi_hist[ibin][ibin2]->Add(more_hists->chargedparticles_cand_phi_hist[ibin][ibin2], wt);
    chargedparticles_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->chargedparticles_cand_eta_hist[ibin][ibin2]->Sumw2();
    chargedparticles_cand_eta_hist[ibin][ibin2]->Add(more_hists->chargedparticles_cand_eta_hist[ibin][ibin2], wt);

    electrons_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->electrons_cand_pT_hist[ibin][ibin2]->Sumw2();
    electrons_cand_pT_hist[ibin][ibin2]->Add(more_hists->electrons_cand_pT_hist[ibin][ibin2], wt);
    electrons_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->electrons_cand_phi_hist[ibin][ibin2]->Sumw2();
    electrons_cand_phi_hist[ibin][ibin2]->Add(more_hists->electrons_cand_phi_hist[ibin][ibin2], wt);
    electrons_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->electrons_cand_eta_hist[ibin][ibin2]->Sumw2();
    electrons_cand_eta_hist[ibin][ibin2]->Add(more_hists->electrons_cand_eta_hist[ibin][ibin2], wt);

    muons_cand_pT_hist[ibin][ibin2]->Sumw2();   more_hists->muons_cand_pT_hist[ibin][ibin2]->Sumw2();
    muons_cand_pT_hist[ibin][ibin2]->Add(more_hists->muons_cand_pT_hist[ibin][ibin2], wt);
    muons_cand_phi_hist[ibin][ibin2]->Sumw2();   more_hists->muons_cand_phi_hist[ibin][ibin2]->Sumw2();
    muons_cand_phi_hist[ibin][ibin2]->Add(more_hists->muons_cand_phi_hist[ibin][ibin2], wt);
    muons_cand_eta_hist[ibin][ibin2]->Sumw2();   more_hists->muons_cand_eta_hist[ibin][ibin2]->Sumw2();
    muons_cand_eta_hist[ibin][ibin2]->Add(more_hists->muons_cand_eta_hist[ibin][ibin2], wt);
    /////


    muons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->muons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
    muons_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->muons_cand_pT_hist_bkg[ibin][ibin2], wt);
    electrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->electrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
    electrons_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->electrons_cand_pT_hist_bkg[ibin][ibin2], wt);
    chargedparticles_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->chargedparticles_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
    chargedparticles_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->chargedparticles_cand_pT_hist_bkg[ibin][ibin2], wt);
    chargedhadrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
    chargedhadrons_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_pT_hist_bkg[ibin][ibin2], wt);
    neutral_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
    neutral_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->neutral_cand_pT_hist_bkg[ibin][ibin2], wt);
    photons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();   more_hists->photons_cand_pT_hist_bkg[ibin][ibin2]->Sumw2();
    photons_cand_pT_hist_bkg[ibin][ibin2]->Add(more_hists->photons_cand_pT_hist_bkg[ibin][ibin2], wt);


    SumJetPtFraction_hist[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist[ibin][ibin2], wt);

    JetShapeIntegratedParticles[ibin][ibin2]->Sumw2();   more_hists->JetShapeIntegratedParticles[ibin][ibin2]->Sumw2();
    JetShapeIntegratedParticles[ibin][ibin2]->Add(more_hists->JetShapeIntegratedParticles[ibin][ibin2], wt);
    JetShapeDiffParticles[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles[ibin][ibin2]->Sumw2();
    JetShapeDiffParticles[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles[ibin][ibin2], wt);

    JetShapeIntegratedParticles_bkgsub[ibin][ibin2]->Sumw2();   more_hists->JetShapeIntegratedParticles_bkgsub[ibin][ibin2]->Sumw2();
    JetShapeIntegratedParticles_bkgsub[ibin][ibin2]->Add(more_hists->JetShapeIntegratedParticles_bkgsub[ibin][ibin2], wt);
 
    JetShapeDiffParticles_bkgsub[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_bkgsub[ibin][ibin2]->Sumw2();
    JetShapeDiffParticles_bkgsub[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_bkgsub[ibin][ibin2], wt);


    JetShapeIntegratedParticles_bkg[ibin][ibin2]->Sumw2();   more_hists->JetShapeIntegratedParticles_bkg[ibin][ibin2]->Sumw2();
    JetShapeIntegratedParticles_bkg[ibin][ibin2]->Add(more_hists->JetShapeIntegratedParticles_bkg[ibin][ibin2], wt);
 
    JetShapeDiffParticles_bkg[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_bkg[ibin][ibin2]->Sumw2();
    JetShapeDiffParticles_bkg[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_bkg[ibin][ibin2], wt);


    JetShapeIntegratedParticles_bkg_1D[ibin][ibin2]->Sumw2();   more_hists->JetShapeIntegratedParticles_bkg_1D[ibin][ibin2]->Sumw2();
    JetShapeIntegratedParticles_bkg_1D[ibin][ibin2]->Add(more_hists->JetShapeIntegratedParticles_bkg_1D[ibin][ibin2], wt);
    JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Sumw2();
    JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_bkg_1D[ibin][ibin2], wt);
    JetShapeDiffParticles_bkg_1D_EtaRef[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_bkg_1D_EtaRef[ibin][ibin2]->Sumw2();
    JetShapeDiffParticles_bkg_1D_EtaRef[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_bkg_1D_EtaRef[ibin][ibin2], wt);

    JetShapeIntegratedParticles_1D[ibin][ibin2]->Sumw2();   more_hists->JetShapeIntegratedParticles_1D[ibin][ibin2]->Sumw2();
    JetShapeIntegratedParticles_1D[ibin][ibin2]->Add(more_hists->JetShapeIntegratedParticles_1D[ibin][ibin2], wt);
    JetShapeDiffParticles_1D[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_1D[ibin][ibin2]->Sumw2();
    JetShapeDiffParticles_1D[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_1D[ibin][ibin2], wt);

    SumPt_in_a_radius[ibin][ibin2]->Sumw2();   more_hists->SumPt_in_a_radius[ibin][ibin2]->Sumw2();
    SumPt_in_a_radius[ibin][ibin2]->Add(more_hists->SumPt_in_a_radius[ibin][ibin2], wt);

    multiplicity_in_a_radius[ibin][ibin2]->Sumw2();   more_hists->multiplicity_in_a_radius[ibin][ibin2]->Sumw2();
    multiplicity_in_a_radius[ibin][ibin2]->Add(more_hists->multiplicity_in_a_radius[ibin][ibin2], wt);

    radius_hist[ibin][ibin2]->Sumw2();   more_hists->radius_hist[ibin][ibin2]->Sumw2();
    radius_hist[ibin][ibin2]->Add(more_hists->radius_hist[ibin][ibin2], wt);

    //// subleading jet histograms
    only_subleadingjets_corrpT[ibin][ibin2]->Sumw2();   more_hists->only_subleadingjets_corrpT[ibin][ibin2]->Sumw2();
    only_subleadingjets_corrpT[ibin][ibin2]->Add(more_hists->only_subleadingjets_corrpT[ibin][ibin2], wt);
    only_subleadingjets_phi[ibin][ibin2]->Sumw2();   more_hists->only_subleadingjets_phi[ibin][ibin2]->Sumw2();
    only_subleadingjets_phi[ibin][ibin2]->Add(more_hists->only_subleadingjets_phi[ibin][ibin2], wt);
    only_subleadingjets_eta[ibin][ibin2]->Sumw2();   more_hists->only_subleadingjets_eta[ibin][ibin2]->Sumw2();
    only_subleadingjets_eta[ibin][ibin2]->Add(more_hists->only_subleadingjets_eta[ibin][ibin2], wt);
    all_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->all_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
    all_cand_pT_hist_subleadingjet[ibin][ibin2]->Add(more_hists->all_cand_pT_hist_subleadingjet[ibin][ibin2], wt);
    all_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->all_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
    all_cand_phi_hist_subleadingjet[ibin][ibin2]->Add(more_hists->all_cand_phi_hist_subleadingjet[ibin][ibin2], wt);
    all_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->all_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
    all_cand_eta_hist_subleadingjet[ibin][ibin2]->Add(more_hists->all_cand_eta_hist_subleadingjet[ibin][ibin2], wt);
    neutral_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
    neutral_cand_pT_hist_subleadingjet[ibin][ibin2]->Add(more_hists->neutral_cand_pT_hist_subleadingjet[ibin][ibin2], wt);
    neutral_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
    neutral_cand_phi_hist_subleadingjet[ibin][ibin2]->Add(more_hists->neutral_cand_phi_hist_subleadingjet[ibin][ibin2], wt);
    neutral_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
    neutral_cand_eta_hist_subleadingjet[ibin][ibin2]->Add(more_hists->neutral_cand_eta_hist_subleadingjet[ibin][ibin2], wt);
   

    photons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->photons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
    photons_cand_pT_hist_subleadingjet[ibin][ibin2]->Add(more_hists->photons_cand_pT_hist_subleadingjet[ibin][ibin2], wt);
    photons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->photons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
    photons_cand_phi_hist_subleadingjet[ibin][ibin2]->Add(more_hists->photons_cand_phi_hist_subleadingjet[ibin][ibin2], wt);
    photons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->photons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
    photons_cand_eta_hist_subleadingjet[ibin][ibin2]->Add(more_hists->photons_cand_eta_hist_subleadingjet[ibin][ibin2], wt);
    
    chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2]->Sumw2();
    chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2], wt);
    chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2]->Sumw2();
    chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2], wt);
    chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2]->Sumw2();
    chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2], wt);
    

    chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();
    chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2], wt);

    neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();
    neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2], wt);

    photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Sumw2();
    photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2], wt);

    NumNeutral_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->NumNeutral_subleadingjet[ibin][ibin2]->Sumw2();
    NumNeutral_subleadingjet[ibin][ibin2]->Add(more_hists->NumNeutral_subleadingjet[ibin][ibin2], wt);
    NumPhotons_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->NumPhotons_subleadingjet[ibin][ibin2]->Sumw2();
    NumPhotons_subleadingjet[ibin][ibin2]->Add(more_hists->NumPhotons_subleadingjet[ibin][ibin2], wt);
    NumAll_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->NumAll_subleadingjet[ibin][ibin2]->Sumw2();
    NumAll_subleadingjet[ibin][ibin2]->Add(more_hists->NumAll_subleadingjet[ibin][ibin2], wt);
    NumChargedHadrons_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->NumChargedHadrons_subleadingjet[ibin][ibin2]->Sumw2();
    NumChargedHadrons_subleadingjet[ibin][ibin2]->Add(more_hists->NumChargedHadrons_subleadingjet[ibin][ibin2], wt);

    NumNeutral_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->NumNeutral_subleadingjet_bkg[ibin][ibin2]->Sumw2();
    NumNeutral_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->NumNeutral_subleadingjet_bkg[ibin][ibin2], wt);
    NumPhotons_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->NumPhotons_subleadingjet_bkg[ibin][ibin2]->Sumw2();
    NumPhotons_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->NumPhotons_subleadingjet_bkg[ibin][ibin2], wt);
    NumAll_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->NumAll_subleadingjet_bkg[ibin][ibin2]->Sumw2();
    NumAll_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->NumAll_subleadingjet_bkg[ibin][ibin2], wt);
    NumChargedHadrons_subleadingjet_bkg[ibin][ibin2]->Sumw2();   more_hists->NumChargedHadrons_subleadingjet_bkg[ibin][ibin2]->Sumw2();
    NumChargedHadrons_subleadingjet_bkg[ibin][ibin2]->Add(more_hists->NumChargedHadrons_subleadingjet_bkg[ibin][ibin2], wt);




    SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet[ibin][ibin2], wt);

    SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2], wt);

    SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2], wt);

    SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2], wt);



/// bkg..
    SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2], wt);

    SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2], wt);

    SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2], wt);

    SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2], wt);


    SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2], wt);

    SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2], wt);





    SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2], wt);
    SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2], wt);

    SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2], wt);

    ///

    SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2], wt);

    SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2], wt);

    SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2], wt);

    SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2], wt);
    ////

/// starts -- bkg sumpt fraction plots
    SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2], wt);

    SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2], wt);

    SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();more_hists->SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2], wt);

    SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2], wt);

    SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2], wt);

    SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2]->Sumw2();   more_hists->SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2]->Sumw2();
    SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2]->Add(more_hists->SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2], wt);
/// end -- bkg sumpt fraction plots




    radius_hist_subleadingjet[ibin][ibin2]->Sumw2();   more_hists->radius_hist_subleadingjet[ibin][ibin2]->Sumw2();
    radius_hist_subleadingjet[ibin][ibin2]->Add(more_hists->radius_hist_subleadingjet[ibin][ibin2], wt);
    Centrality_hist[ibin][ibin2]->Sumw2();   more_hists->Centrality_hist[ibin][ibin2]->Sumw2();
    Centrality_hist[ibin][ibin2]->Add(more_hists->Centrality_hist[ibin][ibin2], wt);

    
    //Aj[ibin][ibin2]->Sumw2();   more_hists->Aj[ibin][ibin2]->Sumw2();
   // Aj[ibin][ibin2]->Add(more_hists->Aj[ibin][ibin2], wt);


    dN_tracks[ibin][ibin2]->Sumw2();   more_hists->dN_tracks[ibin][ibin2]->Sumw2();
    dN_tracks[ibin][ibin2]->Add(more_hists->dN_tracks[ibin][ibin2], wt);
    dN_chargedhadrons[ibin][ibin2]->Sumw2();   more_hists->dN_chargedhadrons[ibin][ibin2]->Sumw2();
    dN_chargedhadrons[ibin][ibin2]->Add(more_hists->dN_chargedhadrons[ibin][ibin2], wt);
    dN_chargedparticles[ibin][ibin2]->Sumw2();   more_hists->dN_chargedparticles[ibin][ibin2]->Sumw2();
    dN_chargedparticles[ibin][ibin2]->Add(more_hists->dN_chargedparticles[ibin][ibin2], wt);
    dN_electrons[ibin][ibin2]->Sumw2();   more_hists->dN_electrons[ibin][ibin2]->Sumw2();
    dN_electrons[ibin][ibin2]->Add(more_hists->dN_electrons[ibin][ibin2], wt);
    dN_muons[ibin][ibin2]->Sumw2();   more_hists->dN_muons[ibin][ibin2]->Sumw2();
    dN_muons[ibin][ibin2]->Add(more_hists->dN_muons[ibin][ibin2], wt);
    dN_neutrals[ibin][ibin2]->Sumw2();   more_hists->dN_neutrals[ibin][ibin2]->Sumw2();
    dN_neutrals[ibin][ibin2]->Add(more_hists->dN_neutrals[ibin][ibin2], wt);
    dN_photons[ibin][ibin2]->Sumw2();   more_hists->dN_photons[ibin][ibin2]->Sumw2();
    dN_photons[ibin][ibin2]->Add(more_hists->dN_photons[ibin][ibin2], wt);

    JetPt_fraction_hist[ibin][ibin2]->Sumw2();   more_hists->JetPt_fraction_hist[ibin][ibin2]->Sumw2();
    JetPt_fraction_hist[ibin][ibin2]->Add(more_hists->JetPt_fraction_hist[ibin][ibin2], wt);

    dPhi_jet_track[ibin][ibin2]->Sumw2();   more_hists->dPhi_jet_track[ibin][ibin2]->Sumw2();
    dPhi_jet_track[ibin][ibin2]->Add(more_hists->dPhi_jet_track[ibin][ibin2], wt);

    dPhi_jet_track_ptweight[ibin][ibin2]->Sumw2();   more_hists->dPhi_jet_track_ptweight[ibin][ibin2]->Sumw2();
    dPhi_jet_track_ptweight[ibin][ibin2]->Add(more_hists->dPhi_jet_track_ptweight[ibin][ibin2], wt);


    track_cand_pT_hist_weighted[ibin][ibin2]->Sumw2();   more_hists->track_cand_pT_hist_weighted[ibin][ibin2]->Sumw2();
    track_cand_pT_hist_weighted[ibin][ibin2]->Add(more_hists->track_cand_pT_hist_weighted[ibin][ibin2], wt);


    JetShapeDiffParticles_bkg_1DEtaRef[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_bkg_1DEtaRef[ibin][ibin2]->Sumw2();
    JetShapeDiffParticles_bkg_1DEtaRef[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_bkg_1DEtaRef[ibin][ibin2], wt);

    JetShapeDiffParticles_bkg_1DEtaRefSig[ibin][ibin2]->Sumw2();   more_hists->JetShapeDiffParticles_bkg_1DEtaRefSig[ibin][ibin2]->Sumw2();
    JetShapeDiffParticles_bkg_1DEtaRefSig[ibin][ibin2]->Add(more_hists->JetShapeDiffParticles_bkg_1DEtaRefSig[ibin][ibin2], wt);




    SumPt_only[ibin][ibin2]->Sumw2();   more_hists->SumPt_only[ibin][ibin2]->Sumw2();
    SumPt_only[ibin][ibin2]->Add(more_hists->SumPt_only[ibin][ibin2], wt);

    SumPt_bkg_only[ibin][ibin2]->Sumw2();   more_hists->SumPt_bkg_only[ibin][ibin2]->Sumw2();
    SumPt_bkg_only[ibin][ibin2]->Add(more_hists->SumPt_bkg_only[ibin][ibin2], wt);

    SumPt_ratio[ibin][ibin2]->Sumw2();   more_hists->SumPt_ratio[ibin][ibin2]->Sumw2();
    SumPt_ratio[ibin][ibin2]->Add(more_hists->SumPt_ratio[ibin][ibin2], wt);


    EtaRef_bkg_pt[ibin][ibin2]->Sumw2(); more_hists->EtaRef_bkg_pt[ibin][ibin2]->Sumw2(); 
    EtaRef_bkg_pt[ibin][ibin2]->Add(more_hists->EtaRef_bkg_pt[ibin][ibin2], wt);

    EtaRef_bkg_pt_weighted[ibin][ibin2]->Sumw2(); more_hists->EtaRef_bkg_pt_weighted[ibin][ibin2]->Sumw2();
    EtaRef_bkg_pt_weighted[ibin][ibin2]->Add(more_hists->EtaRef_bkg_pt_weighted[ibin][ibin2], wt);

    EMix_bkg_pt[ibin][ibin2]->Sumw2(); more_hists->EMix_bkg_pt[ibin][ibin2]->Sumw2();
    EMix_bkg_pt[ibin][ibin2]->Add(more_hists->EMix_bkg_pt[ibin][ibin2], wt);

    EMix_bkg_pt_weighted[ibin][ibin2]->Sumw2(); more_hists->EMix_bkg_pt_weighted[ibin][ibin2]->Sumw2();
    EMix_bkg_pt_weighted[ibin][ibin2]->Add(more_hists->EMix_bkg_pt_weighted[ibin][ibin2], wt);

    /// phi
    EMix_bkg_phi[ibin][ibin2]->Sumw2(); more_hists->EMix_bkg_phi[ibin][ibin2]->Sumw2();
    EMix_bkg_phi[ibin][ibin2]->Add(more_hists->EMix_bkg_phi[ibin][ibin2], wt);
          
    EMix_bkg_phi_weighted[ibin][ibin2]->Sumw2(); more_hists->EMix_bkg_phi_weighted[ibin][ibin2]->Sumw2();
    EMix_bkg_phi_weighted[ibin][ibin2]->Add(more_hists->EMix_bkg_phi_weighted[ibin][ibin2], wt);

    /// eta
    EMix_bkg_eta[ibin][ibin2]->Sumw2(); more_hists->EMix_bkg_eta[ibin][ibin2]->Sumw2();
    EMix_bkg_eta[ibin][ibin2]->Add(more_hists->EMix_bkg_eta[ibin][ibin2], wt);
          
    EMix_bkg_eta_weighted[ibin][ibin2]->Sumw2(); more_hists->EMix_bkg_eta_weighted[ibin][ibin2]->Sumw2();
    EMix_bkg_eta_weighted[ibin][ibin2]->Add(more_hists->EMix_bkg_eta_weighted[ibin][ibin2], wt);

          
    track_bkg_pT_hist_weighted[ibin][ibin2]->Sumw2();   more_hists->track_bkg_pT_hist_weighted[ibin][ibin2]->Sumw2();
    track_bkg_pT_hist_weighted[ibin][ibin2]->Add(more_hists->track_bkg_pT_hist_weighted[ibin][ibin2], wt);
          
    track_bkg_pT_hist[ibin][ibin2]->Sumw2();   more_hists->track_bkg_pT_hist[ibin][ibin2]->Sumw2();
    track_bkg_pT_hist[ibin][ibin2]->Add(more_hists->track_bkg_pT_hist[ibin][ibin2], wt);
          
    /// phi
    track_bkg_phi_hist_weighted[ibin][ibin2]->Sumw2();   more_hists->track_bkg_phi_hist_weighted[ibin][ibin2]->Sumw2();
    track_bkg_phi_hist_weighted[ibin][ibin2]->Add(more_hists->track_bkg_phi_hist_weighted[ibin][ibin2], wt);
          
    track_bkg_phi_hist[ibin][ibin2]->Sumw2();   more_hists->track_bkg_phi_hist[ibin][ibin2]->Sumw2();
    track_bkg_phi_hist[ibin][ibin2]->Add(more_hists->track_bkg_phi_hist[ibin][ibin2], wt);
          
          
    /// eta
    track_bkg_eta_hist_weighted[ibin][ibin2]->Sumw2();   more_hists->track_bkg_eta_hist_weighted[ibin][ibin2]->Sumw2();
    track_bkg_eta_hist_weighted[ibin][ibin2]->Add(more_hists->track_bkg_eta_hist_weighted[ibin][ibin2], wt);
          
    track_bkg_eta_hist[ibin][ibin2]->Sumw2();   more_hists->track_bkg_eta_hist[ibin][ibin2]->Sumw2();
    track_bkg_eta_hist[ibin][ibin2]->Add(more_hists->track_bkg_eta_hist[ibin][ibin2], wt);
       

    track_gen_bkg_pT_hist[ibin][ibin2]->Sumw2();   more_hists->track_gen_bkg_pT_hist[ibin][ibin2]->Sumw2();
    track_gen_bkg_pT_hist[ibin][ibin2]->Add(more_hists->track_gen_bkg_pT_hist[ibin][ibin2], wt);

    track_gen_bkg_eta_hist[ibin][ibin2]->Sumw2();   more_hists->track_gen_bkg_eta_hist[ibin][ibin2]->Sumw2();
    track_gen_bkg_eta_hist[ibin][ibin2]->Add(more_hists->track_gen_bkg_eta_hist[ibin][ibin2], wt);

    track_gen_bkg_phi_hist[ibin][ibin2]->Sumw2();   more_hists->track_gen_bkg_phi_hist[ibin][ibin2]->Sumw2();
    track_gen_bkg_phi_hist[ibin][ibin2]->Add(more_hists->track_gen_bkg_phi_hist[ibin][ibin2], wt);


    emix_jetpt[ibin][ibin2]->Sumw2();   more_hists->emix_jetpt[ibin][ibin2]->Sumw2();
    emix_jetpt[ibin][ibin2]->Add(more_hists->emix_jetpt[ibin][ibin2], wt);

    emix_jetpt_test[ibin][ibin2]->Sumw2();   more_hists->emix_jetpt_test[ibin][ibin2]->Sumw2();
    emix_jetpt_test[ibin][ibin2]->Add(more_hists->emix_jetpt_test[ibin][ibin2], wt);




    for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

      hJetGenTrackBackground[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetGenTrackBackground[ibin][ibin2][ibin3]->Sumw2();
      hJetGenTrackBackground[ibin][ibin2][ibin3]->Add(more_hists->hJetGenTrackBackground[ibin][ibin2][ibin3], wt);

      hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Sumw2();
      hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Add(more_hists->hJetGenTrackBackground_pythia[ibin][ibin2][ibin3], wt);

      hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Sumw2();
      hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Add(more_hists->hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3], wt);


	  hJetTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();
	  hJetTrackSignalBackground[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3], wt);
	
	  hJetTrackBackground[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetTrackBackground[ibin][ibin2][ibin3]->Sumw2();
	  hJetTrackBackground[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackBackground[ibin][ibin2][ibin3], wt);

	  //JetShapeDiffParticles_bkg_1D[ibin][ibin2][ibin3]->Sumw2();   more_hists->JetShapeDiffParticles_bkg_1D[ibin][ibin2][ibin3]->Sumw2();
	  ///JetShapeDiffParticles_bkg_1D[ibin][ibin2][ibin3]->Add(more_hists->JetShapeDiffParticles_bkg_1D[ibin][ibin2][ibin3], wt);

	  //JetShapeDiffParticles_1D[ibin][ibin2][ibin3]->Sumw2();   more_hists->JetShapeDiffParticles_1D[ibin][ibin2][ibin3]->Sumw2();
	  //JetShapeDiffParticles_1D[ibin][ibin2][ibin3]->Add(more_hists->JetShapeDiffParticles_1D[ibin][ibin2][ibin3], wt);


	  Bkg_phi[ibin][ibin2][ibin3]->Sumw2();   more_hists->Bkg_phi[ibin][ibin2][ibin3]->Sumw2();
	  Bkg_phi[ibin][ibin2][ibin3]->Add(more_hists->Bkg_phi[ibin][ibin2][ibin3], wt);

    Bkg_eta[ibin][ibin2][ibin3]->Sumw2();   more_hists->Bkg_eta[ibin][ibin2][ibin3]->Sumw2();
    Bkg_eta[ibin][ibin2][ibin3]->Add(more_hists->Bkg_eta[ibin][ibin2][ibin3], wt);


    Bkg_phi_weighted[ibin][ibin2][ibin3]->Sumw2();   more_hists->Bkg_phi_weighted[ibin][ibin2][ibin3]->Sumw2();
    Bkg_phi_weighted[ibin][ibin2][ibin3]->Add(more_hists->Bkg_phi_weighted[ibin][ibin2][ibin3], wt);

    Bkg_eta_weighted[ibin][ibin2][ibin3]->Sumw2();   more_hists->Bkg_eta_weighted[ibin][ibin2][ibin3]->Sumw2();
    Bkg_eta_weighted[ibin][ibin2][ibin3]->Add(more_hists->Bkg_eta_weighted[ibin][ibin2][ibin3], wt);




    Jet_phi[ibin][ibin2][ibin3]->Sumw2();   more_hists->Jet_phi[ibin][ibin2][ibin3]->Sumw2();
    Jet_phi[ibin][ibin2][ibin3]->Add(more_hists->Jet_phi[ibin][ibin2][ibin3], wt);


    hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();   more_hists->hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();
    hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3], wt);

    GenTrkEta_bkg[ibin][ibin2][ibin3]->Sumw2();   more_hists->GenTrkEta_bkg[ibin][ibin2][ibin3]->Sumw2();
    GenTrkEta_bkg[ibin][ibin2][ibin3]->Add(more_hists->GenTrkEta_bkg[ibin][ibin2][ibin3], wt);

    GenTrkPhi_bkg[ibin][ibin2][ibin3]->Sumw2();   more_hists->GenTrkPhi_bkg[ibin][ibin2][ibin3]->Sumw2();
    GenTrkPhi_bkg[ibin][ibin2][ibin3]->Add(more_hists->GenTrkPhi_bkg[ibin][ibin2][ibin3], wt);



	} // ibin3
    } /// ibin2
  } /// ibin
}


void hist_class::Delete()
{

  delete NumberOfMatches;
  delete NEvents;
  delete NEvents_test;
  delete NEvents_after_noise;
  delete NEvents_after_spike;
  delete NEvents_after_dphi;
  delete NEvents_before_dphi;
  delete NEvents_dijets;
  delete NEvents_after_trigger;

  delete AllJetPt_raw_hist;
  delete AllJetPhi_hist;
  delete AllJetEta_hist;
  delete AllJetPt_hist;
  delete First_AllJetPhi_hist;
  delete First_AllJetEta_hist;
  delete First_AllJetPt_hist;
  delete Sub_AllJetPhi_hist;
  delete Sub_AllJetEta_hist;
  delete Sub_AllJetPt_hist;
  delete JetPt_fraction;
  delete Centrality;
  delete track_vz;
  delete track_vz_weighted;


  for (int ibin=0;ibin<nCBins;ibin++){

  delete jet_pT_hist[ibin];
  delete jet_phi_hist[ibin];
  delete jet_eta_hist[ibin];
  delete jet_corrpT_hist[ibin];
  delete LeadingJetPt_hist[ibin];
  delete LeadingJetPhi_hist[ibin];
  delete LeadingJetEta_hist[ibin];
  delete SubJetPt_hist[ibin];
  delete SubJetPhi_hist[ibin];
  delete SubJetEta_hist[ibin];

  delete dPhi_leadingjet_hist[ibin];
  delete dPhi_subleadingjet_hist[ibin];

  delete JetShapeDiffParticles_bkg_1D_total[ibin];
  delete JetShapeIntegratedParticles_bkg_1D_total[ibin];
  delete SumPt_only_centbin[ibin];

  delete dPhi_hist[ibin];
  delete dPhi_after_hist[ibin];
  delete Aj[ibin];


  for (int ibin2=0;ibin2<nPtBins;ibin2++){

  
  delete ThirdJetPt_hist[ibin][ibin2];
  delete ThirdJetPhi_hist[ibin][ibin2];
  delete ThirdJetEta_hist[ibin][ibin2];
  //delete dPhi_hist[ibin][ibin2];
  //delete dPhi_after_hist[ibin][ibin2];




  delete all_jets_corrpT[ibin][ibin2];
  delete all_jets_phi[ibin][ibin2];
  delete all_jets_eta[ibin][ibin2];

  delete only_leadingjets_corrpT[ibin][ibin2];
  delete only_leadingjets_phi[ibin][ibin2];
  delete only_leadingjets_eta[ibin][ibin2];
  delete neutral_cand_pT_hist[ibin][ibin2];
  delete neutral_cand_phi_hist[ibin][ibin2];
  delete neutral_cand_eta_hist[ibin][ibin2];
  delete photons_cand_pT_hist[ibin][ibin2];
  delete photons_cand_phi_hist[ibin][ibin2];
  delete photons_cand_eta_hist[ibin][ibin2];
  delete NumNeutral[ibin][ibin2];
  delete NumPhotons[ibin][ibin2];
  delete NumAll[ibin][ibin2];
  delete NumCharged[ibin][ibin2];

  delete NumAll_bkg[ibin][ibin2];
  delete NumNeutral_bkg[ibin][ibin2];
  delete NumPhotons_bkg[ibin][ibin2];
  delete NumChargedHadrons_bkg[ibin][ibin2];
  delete NumChargedParticles_bkg[ibin][ibin2];
  delete NumElectrons_bkg[ibin][ibin2];
  delete NumMuons_bkg[ibin][ibin2];

  delete radius_hist[ibin][ibin2];
  delete JetShapeIntegratedParticles[ibin][ibin2];
  delete JetShapeDiffParticles[ibin][ibin2];
  delete JetShapeIntegratedParticles_bkgsub[ibin][ibin2];
  delete JetShapeDiffParticles_bkgsub[ibin][ibin2];
  delete JetShapeIntegratedParticles_bkg[ibin][ibin2];
  delete JetShapeDiffParticles_bkg[ibin][ibin2];



  delete SumJetPtFraction_hist[ibin][ibin2];
  delete all_cand_pT_hist[ibin][ibin2];
  delete all_cand_phi_hist[ibin][ibin2];
  delete all_cand_eta_hist[ibin][ibin2];

  delete track_cand_pT_hist[ibin][ibin2];
  delete track_cand_phi_hist[ibin][ibin2];
  delete track_cand_eta_hist[ibin][ibin2];
  delete track_cand_pT_hist_subleadingjet[ibin][ibin2];

  delete only_subleadingjets_corrpT[ibin][ibin2];
  delete only_subleadingjets_phi[ibin][ibin2];
  delete only_subleadingjets_eta[ibin][ibin2];
  delete neutral_cand_pT_hist_subleadingjet[ibin][ibin2];
  delete neutral_cand_phi_hist_subleadingjet[ibin][ibin2];
  delete neutral_cand_eta_hist_subleadingjet[ibin][ibin2];
  delete photons_cand_pT_hist_subleadingjet[ibin][ibin2];
  delete photons_cand_phi_hist_subleadingjet[ibin][ibin2];
  delete photons_cand_eta_hist_subleadingjet[ibin][ibin2];
  delete NumNeutral_subleadingjet[ibin][ibin2];
  delete NumPhotons_subleadingjet[ibin][ibin2];
  delete NumAll_subleadingjet[ibin][ibin2];
  delete NumChargedHadrons_subleadingjet[ibin][ibin2];





  delete JetShapeDiffParticles_bkg_1D[ibin][ibin2];
  delete JetShapeDiffParticles_bkg_1D_EtaRef[ibin][ibin2];
  delete JetShapeDiffParticles_1D[ibin][ibin2];
  delete JetShapeIntegratedParticles_bkg_1D[ibin][ibin2];
  delete JetShapeIntegratedParticles_1D[ibin][ibin2];




  delete SumJetPtFraction_hist_subleadingjet[ibin][ibin2];
  delete SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2];
  delete SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2];
  delete SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2];

  delete SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2];
  delete SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2];
  delete SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2];


  delete all_cand_pT_hist_subleadingjet[ibin][ibin2];
  delete all_cand_phi_hist_subleadingjet[ibin][ibin2];
  delete all_cand_eta_hist_subleadingjet[ibin][ibin2];
  delete radius_hist_subleadingjet[ibin][ibin2];
  delete Centrality_hist[ibin][ibin2];
  //delete Aj[ibin][ibin2];

  delete SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2];
  delete SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2];
  delete SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2];
  delete SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2];

  delete NumChargedHadrons[ibin][ibin2];
  delete NumChargedParticles[ibin][ibin2];
  delete NumElectrons[ibin][ibin2];
  delete NumMuons[ibin][ibin2];

  delete chargedhadrons_cand_pT_hist[ibin][ibin2];
  delete chargedhadrons_cand_phi_hist[ibin][ibin2];
  delete chargedhadrons_cand_eta_hist[ibin][ibin2];

  delete chargedparticles_cand_pT_hist[ibin][ibin2];
  delete chargedparticles_cand_phi_hist[ibin][ibin2];
  delete chargedparticles_cand_eta_hist[ibin][ibin2];

  delete electrons_cand_pT_hist[ibin][ibin2];
  delete electrons_cand_phi_hist[ibin][ibin2];
  delete electrons_cand_eta_hist[ibin][ibin2];

  delete muons_cand_pT_hist[ibin][ibin2];
  delete muons_cand_phi_hist[ibin][ibin2];
  delete muons_cand_eta_hist[ibin][ibin2];

  delete dN_tracks[ibin][ibin2];
  delete dN_chargedhadrons[ibin][ibin2];
  delete dN_chargedparticles[ibin][ibin2];
  delete dN_electrons[ibin][ibin2];
  delete dN_muons[ibin][ibin2];
  delete dN_neutrals[ibin][ibin2];
  delete dN_photons[ibin][ibin2];

  delete chargedhadrons_cand_pT_hist_bkg[ibin][ibin2];
  delete chargedparticles_cand_pT_hist_bkg[ibin][ibin2];
  delete electrons_cand_pT_hist_bkg[ibin][ibin2];
  delete muons_cand_pT_hist_bkg[ibin][ibin2];
  delete neutral_cand_pT_hist_bkg[ibin][ibin2];
  delete photons_cand_pT_hist_bkg[ibin][ibin2];

  delete SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2];
  delete SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2];
  delete SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2];
  delete SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2];
  delete SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2];
  delete SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2];


  /// subleading jet
  delete chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2];
  delete chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2];
  delete chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2];

  delete chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2];
  delete neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2];
  delete photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2];


  delete SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2];
  delete SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2];
  delete SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2];
  delete SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2];
  delete SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2];
  delete SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2];


  delete NumAll_subleadingjet_bkg[ibin][ibin2];
  delete NumNeutral_subleadingjet_bkg[ibin][ibin2];
  delete NumPhotons_subleadingjet_bkg[ibin][ibin2];
  delete NumChargedHadrons_subleadingjet_bkg[ibin][ibin2];

  delete JetPt_fraction_hist[ibin][ibin2];
  delete dPhi_jet_track[ibin][ibin2];
  delete dPhi_jet_track_ptweight[ibin][ibin2];

  delete SumPt_in_a_radius[ibin][ibin2];
  delete multiplicity_in_a_radius[ibin][ibin2];

  delete track_cand_pT_hist_weighted[ibin][ibin2];
  delete SumPt_only[ibin][ibin2];
  delete JetShapeDiffParticles_bkg_1DEtaRef[ibin][ibin2];
  delete JetShapeDiffParticles_bkg_1DEtaRefSig[ibin][ibin2];

  delete SumPt_bkg_only[ibin][ibin2];
  delete SumPt_ratio[ibin][ibin2];


  delete EtaRef_bkg_pt[ibin][ibin2];
  delete EtaRef_bkg_pt_weighted[ibin][ibin2];
      
  delete EMix_bkg_pt[ibin][ibin2];
  delete EMix_bkg_pt_weighted[ibin][ibin2];

  delete EMix_bkg_phi[ibin][ibin2];
  delete EMix_bkg_phi_weighted[ibin][ibin2];
      
  delete EMix_bkg_eta[ibin][ibin2];
  delete EMix_bkg_eta_weighted[ibin][ibin2];
            
  delete track_bkg_pT_hist[ibin][ibin2];
  delete track_bkg_pT_hist_weighted[ibin][ibin2];

  delete track_bkg_phi_hist[ibin][ibin2];
  delete track_bkg_phi_hist_weighted[ibin][ibin2];
	
  delete track_bkg_eta_hist[ibin][ibin2];
  delete track_bkg_eta_hist_weighted[ibin][ibin2];

  delete track_gen_bkg_pT_hist[ibin][ibin2];
  delete track_gen_bkg_phi_hist[ibin][ibin2];
  delete track_gen_bkg_eta_hist[ibin][ibin2];
  delete emix_jetpt[ibin][ibin2];
  delete emix_jetpt_test[ibin][ibin2];


  for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){
    delete hJetTrackSignalBackground[ibin][ibin2][ibin3];
    delete hJetTrackBackground[ibin][ibin2][ibin3];
    delete hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3];
    delete hJetGenTrackBackground[ibin][ibin2][ibin3];
    delete hJetGenTrackBackground_pythia[ibin][ibin2][ibin3];
    delete hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3];

      //delete JetShapeDiffParticles_bkg_1D[ibin][ibin2][ibin3];
      //delete JetShapeDiffParticles_1D[ibin][ibin2][ibin3];
      delete Jet_phi[ibin][ibin2][ibin3];


      delete Bkg_phi[ibin][ibin2][ibin3];
      delete Bkg_eta[ibin][ibin2][ibin3];

      delete Bkg_phi_weighted[ibin][ibin2][ibin3];
      delete Bkg_eta_weighted[ibin][ibin2][ibin3];

      delete GenTrkEta_bkg[ibin][ibin2][ibin3];
      delete GenTrkPhi_bkg[ibin][ibin2][ibin3];


    } /// ibin3

   } /// ibin2
  } /// ibin
}

void hist_class::Write()
{
  //  TString out_name = (TString) ("root_output/" + desc + "_PbPbHLTJet35U.root");
  TString pT_str = "";
  if( trkPtCut >= 0.49 && trkPtCut < 1.5 ) pT_str = "trkPtCut1";
  else if( trkPtCut >= 1.5 && trkPtCut < 2.5 ) pT_str = "trkPtCut2";
  else if( trkPtCut >= 2.5 && trkPtCut < 3.5 ) pT_str = "trkPtCut3";
  else if( trkPtCut >= 3.5 && trkPtCut < 4.5 ) pT_str = "trkPtCut4";
  else assert(0);  
  
  
   TString out_name = (TString) ("root_output/" + dataset_type_strs[dataset_type_code] + "_" + npart_type_strs[npart] + "_" + pT_str + ".root");
      if( test_it ) out_name = (TString) ("root_output/" + dataset_type_strs[dataset_type_code] + "_" + npart_type_strs[npart] + "_" + pT_str + "_test.root");


  /*TString npart_name = "file_num";   npart_name += npart;
   TString out_name = (TString) ("root_output/" + dataset_type_strs[dataset_type_code] + "_" + npart_name + "_" + pT_str + ".root");
   if( test_it ) out_name = (TString) ("root_output/" + dataset_type_strs[dataset_type_code] + "_" + npart_name + "_" + pT_str + "_test.root");
*/
  TFile *out_file = new TFile(out_name, "RECREATE");

  NumberOfMatches->Write();
  NEvents->Write();
  NEvents_test->Write();
  NEvents_after_noise->Write();
  NEvents_after_spike->Write();
  NEvents_after_trigger->Write();
  NEvents_after_dphi->Write();
  NEvents_before_dphi->Write();
  NEvents_dijets->Write();

  AllJetPt_raw_hist->Write();
  AllJetPhi_hist->Write();
  AllJetEta_hist->Write();
  AllJetPt_hist->Write();

  First_AllJetPhi_hist->Write();
  First_AllJetEta_hist->Write();
  First_AllJetPt_hist->Write();

  Sub_AllJetPhi_hist->Write();
  Sub_AllJetEta_hist->Write();
  Sub_AllJetPt_hist->Write();
  JetPt_fraction->Write();
  Centrality->Write();
   track_vz->Write();
   track_vz_weighted->Write();

  for (int ibin=0;ibin<nCBins;ibin++){

  jet_pT_hist[ibin]->Write();
  jet_phi_hist[ibin]->Write();
  jet_eta_hist[ibin]->Write();
  jet_corrpT_hist[ibin]->Write();
  LeadingJetPt_hist[ibin]->Write();
  LeadingJetPhi_hist[ibin]->Write();
  LeadingJetEta_hist[ibin]->Write();
  SubJetPt_hist[ibin]->Write();
  SubJetPhi_hist[ibin]->Write();
  SubJetEta_hist[ibin]->Write();
  dPhi_leadingjet_hist[ibin]->Write();
  dPhi_subleadingjet_hist[ibin]->Write();
JetShapeDiffParticles_bkg_1D_total[ibin]->Write();
JetShapeDiffParticles_bkg_1D_total[ibin]->Write();
JetShapeIntegratedParticles_bkg_1D_total[ibin]->Write();

SumPt_only_centbin[ibin]->Write();


dPhi_hist[ibin]->Write();
dPhi_after_hist[ibin]->Write();
Aj[ibin]->Write();


  for (int ibin2=0;ibin2<nPtBins;ibin2++){



  ThirdJetPt_hist[ibin][ibin2]->Write();
  ThirdJetPhi_hist[ibin][ibin2]->Write();
  ThirdJetEta_hist[ibin][ibin2]->Write();

  radius_hist[ibin][ibin2]->Write();
  //dPhi_hist[ibin][ibin2]->Write();
  //dPhi_after_hist[ibin][ibin2]->Write();


  all_jets_corrpT[ibin][ibin2]->Write();
  all_jets_phi[ibin][ibin2]->Write();
  all_jets_eta[ibin][ibin2]->Write();
  only_leadingjets_corrpT[ibin][ibin2]->Write();
  only_leadingjets_phi[ibin][ibin2]->Write();
  only_leadingjets_eta[ibin][ibin2]->Write();
  only_subleadingjets_corrpT[ibin][ibin2]->Write();
  only_subleadingjets_phi[ibin][ibin2]->Write();
  only_subleadingjets_eta[ibin][ibin2]->Write();


  neutral_cand_pT_hist[ibin][ibin2]->Write();
  neutral_cand_phi_hist[ibin][ibin2]->Write();
  neutral_cand_eta_hist[ibin][ibin2]->Write();
  photons_cand_pT_hist[ibin][ibin2]->Write();
  photons_cand_phi_hist[ibin][ibin2]->Write();
  photons_cand_eta_hist[ibin][ibin2]->Write();
  NumNeutral[ibin][ibin2]->Write();
  NumPhotons[ibin][ibin2]->Write();
  NumAll[ibin][ibin2]->Write();
  NumCharged[ibin][ibin2]->Write();

  NumAll_bkg[ibin][ibin2]->Write();
  NumNeutral_bkg[ibin][ibin2]->Write();
  NumPhotons_bkg[ibin][ibin2]->Write();
  NumChargedHadrons_bkg[ibin][ibin2]->Write();
  NumChargedParticles_bkg[ibin][ibin2]->Write();
  NumElectrons_bkg[ibin][ibin2]->Write();
  NumMuons_bkg[ibin][ibin2]->Write();

  JetShapeIntegratedParticles[ibin][ibin2]->Write();
  JetShapeDiffParticles[ibin][ibin2]->Write();
  JetShapeIntegratedParticles_bkgsub[ibin][ibin2]->Write();
  JetShapeDiffParticles_bkgsub[ibin][ibin2]->Write();
  JetShapeIntegratedParticles_bkg[ibin][ibin2]->Write();
  JetShapeDiffParticles_bkg[ibin][ibin2]->Write();


  JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Write();
  JetShapeDiffParticles_bkg_1D_EtaRef[ibin][ibin2]->Write();
  JetShapeDiffParticles_1D[ibin][ibin2]->Write();
  JetShapeIntegratedParticles_bkg_1D[ibin][ibin2]->Write();
  JetShapeIntegratedParticles_1D[ibin][ibin2]->Write();



  SumJetPtFraction_hist[ibin][ibin2]->Write();
  all_cand_pT_hist[ibin][ibin2]->Write();
  all_cand_phi_hist[ibin][ibin2]->Write();
  all_cand_eta_hist[ibin][ibin2]->Write();

  track_cand_pT_hist[ibin][ibin2]->Write();
  track_cand_phi_hist[ibin][ibin2]->Write();
  track_cand_eta_hist[ibin][ibin2]->Write();
  track_cand_pT_hist_subleadingjet[ibin][ibin2]->Write();

  neutral_cand_pT_hist_subleadingjet[ibin][ibin2]->Write();
  neutral_cand_phi_hist_subleadingjet[ibin][ibin2]->Write();
  neutral_cand_eta_hist_subleadingjet[ibin][ibin2]->Write();
  photons_cand_pT_hist_subleadingjet[ibin][ibin2]->Write();
  photons_cand_phi_hist_subleadingjet[ibin][ibin2]->Write();
  photons_cand_eta_hist_subleadingjet[ibin][ibin2]->Write();
  NumNeutral_subleadingjet[ibin][ibin2]->Write();
  NumPhotons_subleadingjet[ibin][ibin2]->Write();
  NumAll_subleadingjet[ibin][ibin2]->Write();
  NumChargedHadrons_subleadingjet[ibin][ibin2]->Write();


  SumJetPtFraction_hist_subleadingjet[ibin][ibin2]->Write();
  SumJetPtFraction_hist_subleadingjet_photons[ibin][ibin2]->Write();
  SumJetPtFraction_hist_subleadingjet_neutrals[ibin][ibin2]->Write();
  SumJetPtFraction_hist_subleadingjet_chargedhadrons[ibin][ibin2]->Write();

  SumJetPtFraction_hist_leadingjet_photons[ibin][ibin2]->Write();
  SumJetPtFraction_hist_leadingjet_neutrals[ibin][ibin2]->Write();
  SumJetPtFraction_hist_leadingjet_charged[ibin][ibin2]->Write();

  SumJetPtFraction_hist_leadingjet_chargedhadrons[ibin][ibin2]->Write();
  SumJetPtFraction_hist_leadingjet_chargedparticles[ibin][ibin2]->Write();
  SumJetPtFraction_hist_leadingjet_electrons[ibin][ibin2]->Write();
  SumJetPtFraction_hist_leadingjet_muons[ibin][ibin2]->Write();

  NumChargedHadrons[ibin][ibin2]->Write();
  NumChargedParticles[ibin][ibin2]->Write();
  NumElectrons[ibin][ibin2]->Write();
  NumMuons[ibin][ibin2]->Write();

  chargedhadrons_cand_pT_hist[ibin][ibin2]->Write();
  chargedhadrons_cand_phi_hist[ibin][ibin2]->Write();
  chargedhadrons_cand_eta_hist[ibin][ibin2]->Write();

  chargedparticles_cand_pT_hist[ibin][ibin2]->Write();
  chargedparticles_cand_phi_hist[ibin][ibin2]->Write();
  chargedparticles_cand_eta_hist[ibin][ibin2]->Write();

  electrons_cand_pT_hist[ibin][ibin2]->Write();
  electrons_cand_phi_hist[ibin][ibin2]->Write();
  electrons_cand_eta_hist[ibin][ibin2]->Write();

  muons_cand_pT_hist[ibin][ibin2]->Write();
  muons_cand_phi_hist[ibin][ibin2]->Write();
  muons_cand_eta_hist[ibin][ibin2]->Write();


  all_cand_pT_hist_subleadingjet[ibin][ibin2]->Write();
  all_cand_phi_hist_subleadingjet[ibin][ibin2]->Write();
  all_cand_eta_hist_subleadingjet[ibin][ibin2]->Write();
  radius_hist_subleadingjet[ibin][ibin2]->Write();
  Centrality_hist[ibin][ibin2]->Write();
  //Aj[ibin][ibin2]->Write();


  dN_tracks[ibin][ibin2]->Write();
  dN_chargedhadrons[ibin][ibin2]->Write();
  dN_chargedparticles[ibin][ibin2]->Write();
  dN_electrons[ibin][ibin2]->Write();
  dN_muons[ibin][ibin2]->Write();
  dN_neutrals[ibin][ibin2]->Write();
  dN_photons[ibin][ibin2]->Write();

  chargedhadrons_cand_pT_hist_bkg[ibin][ibin2]->Write();
  chargedparticles_cand_pT_hist_bkg[ibin][ibin2]->Write();
  electrons_cand_pT_hist_bkg[ibin][ibin2]->Write();
  muons_cand_pT_hist_bkg[ibin][ibin2]->Write();
  neutral_cand_pT_hist_bkg[ibin][ibin2]->Write();
  photons_cand_pT_hist_bkg[ibin][ibin2]->Write();

  SumJetPtFraction_hist_leadingjet_photons_bkg[ibin][ibin2]->Write();
  SumJetPtFraction_hist_leadingjet_neutrals_bkg[ibin][ibin2]->Write();
  SumJetPtFraction_hist_leadingjet_chargedhadrons_bkg[ibin][ibin2]->Write();
  SumJetPtFraction_hist_leadingjet_chargedparticles_bkg[ibin][ibin2]->Write();
  SumJetPtFraction_hist_leadingjet_electrons_bkg[ibin][ibin2]->Write();
  SumJetPtFraction_hist_leadingjet_muons_bkg[ibin][ibin2]->Write();


  /// subleading jet plots
  chargedhadrons_cand_pT_hist_subleadingjet[ibin][ibin2]->Write();
  chargedhadrons_cand_phi_hist_subleadingjet[ibin][ibin2]->Write();
  chargedhadrons_cand_eta_hist_subleadingjet[ibin][ibin2]->Write();

  chargedhadrons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Write();
  neutral_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Write();
  photons_cand_pT_hist_subleadingjet_bkg[ibin][ibin2]->Write();


  SumJetPtFraction_hist_subleadingjet_photons_bkg[ibin][ibin2]->Write();
  SumJetPtFraction_hist_subleadingjet_neutrals_bkg[ibin][ibin2]->Write();
  SumJetPtFraction_hist_subleadingjet_chargedhadrons_bkg[ibin][ibin2]->Write();
  SumJetPtFraction_hist_subleadingjet_chargedparticles_bkg[ibin][ibin2]->Write();
  SumJetPtFraction_hist_subleadingjet_electrons_bkg[ibin][ibin2]->Write();
  SumJetPtFraction_hist_subleadingjet_muons_bkg[ibin][ibin2]->Write();


  NumAll_subleadingjet_bkg[ibin][ibin2]->Write();
  NumNeutral_subleadingjet_bkg[ibin][ibin2]->Write();
  NumPhotons_subleadingjet_bkg[ibin][ibin2]->Write();
  NumChargedHadrons_subleadingjet_bkg[ibin][ibin2]->Write();

  JetPt_fraction_hist[ibin][ibin2]->Write();
  dPhi_jet_track[ibin][ibin2]->Write();
  dPhi_jet_track_ptweight[ibin][ibin2]->Write();

  SumPt_in_a_radius[ibin][ibin2]->Write();
  multiplicity_in_a_radius[ibin][ibin2]->Write();
  track_cand_pT_hist_weighted[ibin][ibin2]->Write();


  JetShapeDiffParticles_bkg_1DEtaRef[ibin][ibin2]->Write();
  JetShapeDiffParticles_bkg_1DEtaRefSig[ibin][ibin2]->Write();


  SumPt_only[ibin][ibin2]->Write();
  SumPt_bkg_only[ibin][ibin2]->Write();
  SumPt_ratio[ibin][ibin2]->Write();

  EtaRef_bkg_pt[ibin][ibin2]->Write();
  EtaRef_bkg_pt_weighted[ibin][ibin2]->Write();

  EMix_bkg_pt[ibin][ibin2]->Write();
  EMix_bkg_pt_weighted[ibin][ibin2]->Write();
      
  EMix_bkg_phi[ibin][ibin2]->Write();
  EMix_bkg_phi_weighted[ibin][ibin2]->Write();
      
  EMix_bkg_eta[ibin][ibin2]->Write();
  EMix_bkg_eta_weighted[ibin][ibin2]->Write();

  track_bkg_pT_hist[ibin][ibin2]->Write();
  track_bkg_pT_hist_weighted[ibin][ibin2]->Write();
      
  track_bkg_phi_hist[ibin][ibin2]->Write();
  track_bkg_phi_hist_weighted[ibin][ibin2]->Write();
      
  track_bkg_eta_hist[ibin][ibin2]->Write();
  track_bkg_eta_hist_weighted[ibin][ibin2]->Write();

  track_gen_bkg_pT_hist[ibin][ibin2]->Write();
  track_gen_bkg_phi_hist[ibin][ibin2]->Write();
  track_gen_bkg_eta_hist[ibin][ibin2]->Write();

  emix_jetpt[ibin][ibin2]->Write();
  emix_jetpt_test[ibin][ibin2]->Write();

  for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){
    hJetTrackSignalBackground[ibin][ibin2][ibin3]->Write();
    hJetTrackBackground[ibin][ibin2][ibin3]->Write();
    hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Write();
    hJetGenTrackBackground[ibin][ibin2][ibin3]->Write();
    hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Write();
    hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Write();

    //JetShapeDiffParticles_bkg_1D[ibin][ibin2][ibin3]->Write();
    //JetShapeDiffParticles_1D[ibin][ibin2][ibin3]->Write();
    Jet_phi[ibin][ibin2][ibin3]->Write();

    Bkg_phi[ibin][ibin2][ibin3]->Write();
    Bkg_eta[ibin][ibin2][ibin3]->Write();
      
    Bkg_phi_weighted[ibin][ibin2][ibin3]->Write();
    Bkg_eta_weighted[ibin][ibin2][ibin3]->Write();

    GenTrkEta_bkg[ibin][ibin2][ibin3]->Write();
    GenTrkPhi_bkg[ibin][ibin2][ibin3]->Write();

  } /// ibin3

  } /// ptbin2

  }  //centralitybin
 out_file->Close();
}



double mydeltaR(double eta1, double phi1, double eta2, double phi2);
double mydeltaR(double eta1, double phi1, double eta2, double phi2){
	  	   double deltaEta = fabs(eta1-eta2);
	  	   double deltaPhi = fabs(phi1-phi2);
	  	   if (deltaPhi>TMath::TwoPi())   std::cout << "deltaR calculation. Why is deltaPhi larger than 2pi?" << std::endl;
	  	   else if (deltaPhi>TMath::Pi()) deltaPhi = TMath::TwoPi()-deltaPhi;
		   return pow(deltaEta*deltaEta+deltaPhi*deltaPhi,0.5);
}


void StudyFiles(std::vector<TString> file_names, std::vector<TString> MB_file_names, int foi, hist_class *my_hists);
//void GetFilesOfFileNames(std::vector<TString> &files_of_file_names, std::vector<float> &Xsections); ////getting the data/simulation based on whatever dataset you pick (arg)
void GetFilesOfFileNames(std::vector<TString> &files_of_file_names, std::vector<float> &Xsections, std::vector<double> &assumed_n_evt); ////getting the data/simulation based on whatever dataset you pick (arg);

void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug=false);
float GetDatasetWeight(double n_evt_raw, double Xsection);//
//void GetBkgShape(TVector3 highest_jet_vec, TVector3 second_highest_jet_vec, std::vector<TVector3>& bkg_particles);  
//double GetBkgShape(TVector3 highest_jet_vec, TVector3 second_highest_jet_vec, std::vector<TVector3>& bkg_particles);  
void GetBkgShape(TVector3 highest_jet_vec, TVector3 second_highest_jet_vec, TVector3& bkg_dir, TVector3& bkg_dir2);

int main(int argc, char *argv[])
{




  assert(argc == 4);
  dataset_type_code = atoi(argv[1]);    //// pick datasets you want to run over
  trkPtCut = atof(argv[2]);
  npart = atoi(argv[3]);
  assert(npart >= 0 && npart < e_n_npart_types);
 // assert(npart >= 0);
  assert(trkPtCut > 0. && trkPtCut < 5.);
  std::cout << "Running with trkPtCut " << trkPtCut << " and npart: " << npart << std::endl;
  
  
  // if( dataset_type_code == e_HydJet50 || dataset_type_code == e_HydJet80 || dataset_type_code == e_HydJet110 || dataset_type_code == e_HydJet170 || dataset_type_code == e_Pyquen15 || dataset_type_code == e_Pyquen30 ||dataset_type_code == e_Pyquen50 || dataset_type_code == e_Pyquen80 || dataset_type_code == e_Pyquen110 || dataset_type_code == e_Pyquen170) is_data = false;


  if(dataset_type_code == e_Data2011 ) is_data = true;

    else if( dataset_type_code == e_HydJet30 || dataset_type_code == e_HydJet50 || dataset_type_code == e_HydJet80|| dataset_type_code == e_HydJet100|| dataset_type_code == e_HydJet120|| dataset_type_code == e_HydJet170|| dataset_type_code == e_HydJet200 || dataset_type_code == e_HydJet250 || dataset_type_code == e_HydJet300) is_data =false;
      else assert(0);



      std::vector<TString> files_of_file_names;   files_of_file_names.clear();
      std::vector<TString> MB_files_of_file_names;   MB_files_of_file_names.clear();


  if(is_data) {
     MB_files_of_file_names.push_back("minbias.txt");
  }else{
     MB_files_of_file_names.push_back("minbias_mc.txt");
  }
  
  
  std::vector<float> Xsections;   Xsections.clear();
  std::vector<double> assumed_n_evt;   assumed_n_evt.clear();

  ////// Get list of files for each dataset ... where do I look for the dataset
  std::cout << "GetFilesOf" << std::endl;

  //GetFilesOfFileNames(files_of_file_names, Xsections);
  GetFilesOfFileNames(files_of_file_names, Xsections, assumed_n_evt);

  //// Total histograms, will add to this for each dataset
  hist_class *hists = new hist_class((TString) ("hists"), is_data);
  //// Loop over each dataset (dataset = collection of similar files of event types)
  for(int foi = 0; foi < (int) files_of_file_names.size(); foi++) {
    TString dataset_str = "sim_dataset_";   dataset_str += foi;
    hist_class *these_hists = new hist_class((TString) ("hists_" + dataset_str), is_data);
    std::cout << "Got hist class for foi " << foi << std::endl;
    std::vector<TString> file_names;   file_names.clear();

    ////// Collect the file names for this dataset
    ReadFileList( file_names, files_of_file_names.at(foi), true);
    std::vector<TString> MB_file_names;   MB_file_names.clear();
    ReadFileList( MB_file_names, MB_files_of_file_names.at(foi), true);

    ////// Pass these files to StudyFiles function to fill histograms
    StudyFiles(file_names, MB_file_names, foi, these_hists);
   

    float dataset_scale_val = 1.;
    if( !is_data ) {    ////// determine the proper weighting for this dataset (how much you've simulated and likelihood of event)
     // assert(foi < (int) Xsections.size() );
     //dataset_scale_val = GetDatasetWeight(these_hists->n_evt_raw, Xsections.at(foi));
     //dataset_scale_val = GetDatasetWeight(assumed_n_evt.at(foi), Xsections.at(foi));

     // std::cout << "foi: " << foi << ", xsec: " << Xsections.at(foi) << ", n_evt_raw: " << these_hists->n_evt_raw << ", scale: " << dataset_scale_val << "\n";
    }
    /////// Add the histograms for this dataset to the full histogram class that you care about
    std::cout << "Addhists" << std::endl;
    hists->AddHists(these_hists, dataset_scale_val);
    std::cout << "Delete" << std::endl;
    these_hists->Delete();
    std::cout << "Deleted" << std::endl;
    delete these_hists;
  }
  
  
  hists->NormalizeByCounts();
  hists->Write();
}



void StudyFiles(std::vector<TString> file_names, std::vector<TString> MB_file_names, int foi, hist_class *my_hists)
{

  //////////###### PTHAT SAMPLES ###########///////////////
  ///TFile * wtfile_vtx_mc;

  int pthat =0;
  int pthatmax =0;

  if( dataset_type_code == e_HydJet50 ) {
    pthat = 50;
    pthatmax=80;

  }
  if( dataset_type_code == e_HydJet80 ) {
    pthat = 80;
    //    pthatmax=100;
    pthatmax=120;

  }
  if( dataset_type_code == e_HydJet100 ) {
    pthat = 100;
    pthatmax=120;

  }
  if( dataset_type_code == e_HydJet120 ) {
    pthat = 120;
    //pthatmax=170;
    pthatmax=300;

  }
  if( dataset_type_code == e_HydJet170 ) {
    pthat = 170;
    pthatmax=200;
  }
  if( dataset_type_code == e_HydJet200 ) {
    pthat = 200;
    pthatmax=250;
  }
  if( dataset_type_code == e_HydJet250 ) {
    pthat = 250;
    pthatmax=300;
  }
  if( dataset_type_code == e_HydJet300 ) {
    pthat = 300;
    pthatmax=999;
  }


  ///========= vertex weight =======////

  TH1F* hWeight_vtx;
  TH1F* hWeight_MC_vtx;

  TH1F* hWeight_cent;
  TH1F* hWeight_MC_cent;

  ////===========  data vertex =========//////
  if( !is_data ) {
    TFile * wtfile_vtx = TFile::Open("pbpb_data_vertex_cent.root", "readonly");
    hWeight_vtx = (TH1F*) ((TH1F*)wtfile_vtx->Get("hists_track_vz"))->Clone("hists_track_vz_clone");
    hWeight_cent = (TH1F*) ((TH1F*)wtfile_vtx->Get("hists_Centrality"))->Clone("hists_centrality_clone");

    hWeight_vtx->Scale(1./hWeight_vtx->Integral());
  hWeight_cent->Scale(1./hWeight_cent->Integral());

  // wtfile_vtx->Close();

  TFile * wtfile_vtx_mc = TFile::Open("HydjetMerged_vertex_cent.root", "readonly");
  hWeight_MC_vtx = (TH1F*) ((TH1F*)wtfile_vtx_mc->Get("hists_track_vz"))->Clone("hists_track_vz_clone");
  hWeight_MC_cent = (TH1F*) ((TH1F*)wtfile_vtx_mc->Get("hists_Centrality"))->Clone("hists_centrality_clone");


  hWeight_MC_vtx->Scale(1./hWeight_MC_vtx->Integral());
  hWeight_MC_cent->Scale(1./hWeight_MC_cent->Integral());
}





  //////////###### centrality Bins ###########///////////////

  std::cout << "I am working\n";

 
    
  int mb_evi = 1;     //// minimum bias events
 // assert(npart >= 0 && npart < (int) file_names.size() );

  std::vector<int> MB_evt_nums_found;     MB_evt_nums_found.clear();
        
  for(int fi = 0; fi < (int) file_names.size(); fi++) {
    //if( fi != npart ) continue;
   
    TFile *my_file = TFile::Open(file_names.at(fi));
    std::cout << "file it" << ", file_name: " << file_names.at(fi) << ", number " << fi << " of " << file_names.size() << std::endl;
    if(my_file->IsZombie()) {
      std::cout << "Is zombie" << std::endl;
    }

    TFile *MB_file = TFile::Open(MB_file_names.at(0));
    std::cout << "file it MB" << ", file_name_MB: " << MB_file_names.at(0) << ", number_MB " << fi << " of " << MB_file_names.size() << std::endl;

    if(MB_file->IsZombie()) {
      std::cout << "Is zombie" << std::endl;
    }


    TTree *inp_tree;
    inp_tree = (TTree*)  my_file->Get("mixing_tree");
    mixing_tree *my_primary = new mixing_tree(inp_tree);

    int n_evt = my_primary->fChain->GetEntriesFast();
    std::cout << " testTree_NEvents : " << n_evt << std::endl;

    TTree *inp_tree_MB;
    inp_tree_MB = (TTree*)  MB_file->Get("mixing_tree");
    mixing_tree_MB *my_primary_MB = new mixing_tree_MB(inp_tree_MB);

    std::cout << "Got CT" << std::endl;
//    std::cout << "pointers to MV_file: " << MB_file << ", mixing tree: " << my_primary_MB << std::endl;
  //  std::cout << "pointer to fChain: " << my_primary_MB->fChain << std::endl;

    int n_MB_evt = my_primary_MB->fChain->GetEntriesFast();
    std::cout << " testTree_NEvents_MB : " << n_MB_evt << std::endl;

   // TF1 *fcen = new TF1("fcen","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x)",0,40);
   // TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);

    ///==========================   Event Loop starts ==================================================
    ///==========================   Event Loop starts ==================================================
    ///==========================   Event Loop starts ==================================================

      
      int min_ev_num = 0.;      int max_ev_num = 0;
      min_ev_num = npart*15000;
      max_ev_num = (npart+1)*15000;
      if( npart == e_AllParts ) {
        min_ev_num = 0;
        max_ev_num = my_primary->fChain->GetEntriesFast();
        std::cout << "Run over all events from " << min_ev_num << " to " << max_ev_num << "\n";
      }
      //   std::cout << "min_ev_num " << min_ev_num << ", max_ev_num: " << max_ev_num << ", out of possible " << n_evt << "\n";
      for(int evi = min_ev_num; evi < max_ev_num; evi++) {
        if( evi > n_evt ) break;

      //if( evi % 100 == 0 ) std::cout << "I am running on file " << fi << " of " << ((int) file_names.size()) << ", evi: " << evi << " of " << n_evt << std::endl;
     
          if (evi%1000==0) std::cout << " I am running on file " << fi << " of " << ((int) file_names.size()) << ", evi: " << evi << " of " << n_evt <<  "  pthat : " << pthat << "  pthatmax : " << pthatmax  << std::endl;
           
      
      my_primary->fChain->GetEntry(evi);
      my_hists->NEvents->Fill(my_primary->fChain->GetEntry(evi));
      my_hists->n_evt_raw++;



      if( !is_data ) {
        double pthat_new = my_primary->pthat;
        if (pthat_new > pthatmax) continue;
      }



      //! Centrality reweighting function
      //fcen->SetParameters(1.98261e-02,5.55963e+00,-1.34951e-01,1.70895e-03,-9.28386e-05);
      //! vertex z reqeighting
      //fVz->SetParameters(7.62788e-01,-1.13228e-02,5.85199e-03,-3.04550e-04,4.43440e-05);

      Int_t hiBin = my_primary->hiBin;
      Float_t vertex = 0;
      vertex = my_primary->vz->at(0);

      /*int vtx_bin = hWeight_vtx->GetXaxis()->FindBin(vz);
        float num = hWeight_vtx->GetBinContent(vtx_bin);  /// DATA
        float denom = hWeight_MC_vtx->GetBinContent(vtx_bin); /// MC

        double wvz = 1.;
        if( denom > 0.0001 ) wvz = num / denom;

        double wcen=1;
        int cent_bin = hWeight_cent->GetXaxis()->FindBin(hiBin);
        float num_cent = hWeight_cent->GetBinContent(cent_bin);  /// DATA
        float denom_cent = hWeight_MC_cent->GetBinContent(cent_bin); /// MC
        if( denom_cent > 0.0001 ) wcen = num_cent / denom_cent;
       */


      if(fabs(vertex) > 15.) continue;
      double wvz=1;
      double wcen=1;

      if(is_data) {
        int noise_event_selection = my_primary->pHBHENoiseFilter;
        if(noise_event_selection==0) continue;

        int event_selection = my_primary->pcollisionEventSelection;
        if(event_selection==0) continue;
        my_hists->NEvents_test->Fill(1.0);

        if (my_primary->HLT_HIJet80_v1==0) continue;

        wvz=1;
        wcen=1;

      }else{
        //wvz=fVz->Eval(vertex);
        //wcen = fcen->Eval(hiBin);
        wvz=1.;
        wcen=1.;
        my_hists->track_vz_weighted->Fill(vertex, wvz*wcen);
      }


      ///==========================   MinBias-Mixing starts ==================================///
      ///==========================   MinBias-Mixing starts ==================================///
      ///==========================   MinBias-Mixing starts ==================================///


      ///====  Here: Load MB tree ===////
      int n_MBevs_found = 0;
      bool n_matched_duplicates = 0;
      int n_MB_resets = 0;
     if (evi%1000==0) std::cout << "event " << evi << "  passes selection" << std::endl;
      while(n_MBevs_found < 40  && n_MB_resets < 10 ) {     
        mb_evi++;
        if( mb_evi >= n_MB_evt ) {     ///// reset
          mb_evi = 1; 
          n_MB_resets++;
          std::cout << "Hit end of MB, reset it for the  " << n_MBevs_found << " time\n";
        }
          
          my_primary_MB->fChain->GetEntry(mb_evi);

          if(is_data) {
            if(my_primary_MB->pcollisionEventSelection==0) continue;
            if(my_primary_MB->HLT_HIMinBiasHfOrBSC_v1==0) continue;
          }

          ///==========================   Matching starts ==================================///
          ///==========================   Matching starts ==================================///
          ///==========================   Matching starts ==================================///

          ///// 1.vertex matching     
          Float_t vertex_data = my_primary->vz->at(0);
          Float_t vertex_minbias= my_primary_MB->vz->at(0);

          if( fabs(vertex_data - vertex_minbias) > 1.0 ) continue;

          float minbias_cent = my_primary_MB->hiBin;  
          float jet_cent = my_primary->hiBin;

          //// 2.centrality matching
            if( fabs(jet_cent-minbias_cent) / jet_cent > 0.1 ) continue;

          //// 3.event-plane matching
          // if( fabs(my_primary_MB->hiEvtPlanes - my_primary->hiEvtPlanes) > 0.25 ) continue;
          // if( fabs(my_primary_MB->hiEvtPlanes - my_primary->hiEvtPlanes) > 0.5 ) continue;

          //std::cout << "evi: " << evi << ", mb_evi: " << mb_evi << " of " << n_MB_evt << ", n found: " << n_MBevs_found << ", cent match , mb: " << minbias_cent << ", jet: " << jet_cent << "\n";



/*
          float n_JetTrig_tracks = 0;
          for(int tracks =0; tracks < (int) my_primary->trkPt->size(); tracks++){
            if(fabs(my_primary->trkEta->at(tracks))>2.4) continue;
            if (my_primary->trkAlgo->at(tracks)<4 || my_primary->highPurity->at(tracks)==1) {
              if(my_primary->trkPt->at(tracks)>trkPtCut) {
                if(my_primary->trkPt->at(tracks)>10) continue; 
                n_JetTrig_tracks++;

              }
            }
          }

           float  n_MB_tracks = 0;
          for(int tracks =0; tracks < (int) my_primary_MB->trkPt->size(); tracks++){
            if(fabs(my_primary_MB->trkEta->at(tracks))>2.4) continue;
            if (my_primary_MB->trkAlgo->at(tracks)<4 || my_primary_MB->highPurity->at(tracks)==1) {
              if(my_primary_MB->trkPt->at(tracks)>trkPtCut) {
                if(my_primary_MB->trkPt->at(tracks)>10)continue ;
                n_MB_tracks++;

              }
            }
          }
               // n_MB_tracks = my_primary_MB->trkPt->size();


          //// 4. ntracks matching
          if( fabs(n_JetTrig_tracks-n_MB_tracks) / n_JetTrig_tracks > 0.1 ) continue;
*/
          //std::cout << " ntracks_jettriggered: " << n_JetTrig_tracks  << "   ntracks_minbias: " << n_MB_tracks <<  "   selected : " << fabs(n_JetTrig_tracks-n_MB_tracks) / n_JetTrig_tracks <<  std::endl;

         //std::cout << "passed all MB selection\n";



         ////// MB event is ok
          n_MBevs_found++;
          MB_evt_nums_found.push_back(mb_evi);


          ///====== do dijet selection here ====//////
          TVector3 highest_jet_vec;
          highest_jet_vec.SetPtEtaPhi(0, 0, 0);

          TVector3 second_highest_jet_vec;
          second_highest_jet_vec.SetPtEtaPhi(0, 0, 0);

          TVector3 mb_second_highest_jet_vec;
          mb_second_highest_jet_vec.SetPtEtaPhi(0, 0, 0);

     const double etacut = 2.0 ;
     const double dphicut = 5.*(TMath::Pi())/6. ;
     const double leadingjetcut = 120. ;
     const double subleadjetcut = 50. ;

     double lead_pt=0. ;
     int lead_index=-1 ;
     double sublead_pt=0. ;
     int second_highest_idx=-1 ;
     int highest_idx=-1 ;
     int mb_second_highest_idx=-1 ;

     if(is_dijet){
       //search for leading jet
       for(int j4i = 0; j4i < my_primary->jtpt->size() ; j4i++) {
         if(TMath::Abs(my_primary->jteta->at(j4i))>etacut) continue ;
         //   if(my_primary->trackMax->at(j4i)/my_primary->jtpt->at(j4i) <=0.01) continue ;
         double jet_pt= my_primary->jtpt->at(j4i);
         if(jet_pt<120) continue ;
         if(jet_pt >lead_pt){
           lead_pt=jet_pt;
           highest_idx=j4i;
           highest_jet_vec.SetPtEtaPhi(jet_pt, my_primary->jteta->at(j4i), my_primary->jtphi->at(j4i));

         }
       } //search for leading jet loop from jet-triggered sample

       //search for subleading jet in the jet-triggered sample
       for(int j4i = 0; j4i < my_primary->jtpt->size() ; j4i++) {
         if(j4i==highest_idx) continue ;
         if(TMath::Abs(my_primary->jteta->at(j4i))>etacut) continue ;
         //     if(my_primary->trackMax->at(j4i)/my_primary->jtpt->at(j4i) <=0.01) continue ;
         double jet_pt= my_primary->jtpt->at(j4i);
         if(my_primary->jtpt->at(j4i)< 50) continue ;
         if(my_primary->jtpt->at(j4i) > sublead_pt){
           sublead_pt=my_primary->jtpt->at(j4i);
           second_highest_idx=j4i;
           second_highest_jet_vec.SetPtEtaPhi(my_primary->jtpt->at(j4i), my_primary->jteta->at(j4i), my_primary->jtphi->at(j4i));
         }

       } //search for subleading jet loop from jet-triggered sample




  //search for subleading jet in the minbias sample
  /*for(int ijet = 0 ; ijet <my_primary_MB->jtpt->size(); ijet++){
  
  if(TMath::Abs(my_primary_MB->jteta->at(ijet))> etacut) continue ;
 // if(ijet>0) continue;
  if(my_primary_MB->trackMax->at(ijet)/my_primary_MB->jtpt->at(ijet) <=0.01) continue ;

  if(my_primary_MB->jtpt->at(ijet)<50) continue ;

    if(my_primary_MB->jtpt->at(ijet) > sublead_pt){
    sublead_pt=my_primary_MB->jtpt->at(ijet);
    mb_second_highest_idx=ijet;
    mb_second_highest_jet_vec.SetPtEtaPhi(my_primary_MB->jtpt->at(ijet), my_primary_MB->jteta->at(ijet), my_primary_MB->jtphi->at(ijet));

    }
  }  //end of subleading jet search
*/


  if( highest_idx < 0 ) continue;
  if( second_highest_idx < 0 ) continue;
 // if( mb_second_highest_idx < 0 ) continue;

  if(my_primary->jtpt->at(highest_idx)< 120. ) continue;
  if( my_primary->jtpt->at(second_highest_idx) < 50) continue;
  //if( my_primary_MB->jtpt->at(mb_second_highest_idx) < 50) continue;

 
  float dphi = highest_jet_vec.DeltaPhi(second_highest_jet_vec);
  if (dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  
    for (int ibin=0;ibin<nCBins;ibin++){
     if (my_primary->hiBin >=CBins[ibin] && my_primary->hiBin<CBins[ibin+1]){
            my_hists->dPhi_hist[ibin]->Fill(fabs(dphi));
     }
    }

    if(fabs(dphi)<=(5.*TMath::Pi())/6.) continue;
    my_hists->NEvents_after_dphi->Fill(1.0);

    if(TMath::Abs(my_primary->jteta->at(highest_idx)) > 1.6 ) continue;
    if(TMath::Abs(my_primary->jteta->at(second_highest_idx)) > 1.6 ) continue;


  for (int ibin=0;ibin<nCBins;ibin++){
    if (my_primary->hiBin >=CBins[ibin] && my_primary->hiBin<CBins[ibin+1]){

        my_hists->LeadingJetPt_hist[ibin]->Fill(my_primary->jtpt->at(highest_idx));
        my_hists->LeadingJetPhi_hist[ibin]->Fill(my_primary->jtphi->at(highest_idx));
        my_hists->LeadingJetEta_hist[ibin]->Fill(my_primary->jteta->at(highest_idx));

        my_hists->SubJetPt_hist[ibin]->Fill(my_primary->jtpt->at(second_highest_idx));
        my_hists->SubJetPhi_hist[ibin]->Fill(my_primary->jtphi->at(second_highest_idx));
        my_hists->SubJetEta_hist[ibin]->Fill(my_primary->jteta->at(second_highest_idx));

//       double  Aj = (my_primary->jtpt->at(highest_idx) - my_primary_MB->jtpt->at(mb_second_highest_idx))/(my_primary->jtpt->at(highest_idx) + my_primary_MB->jtpt->at(mb_second_highest_idx));
       
       double  Aj = (my_primary->jtpt->at(highest_idx) - my_primary->jtpt->at(second_highest_idx))/(my_primary->jtpt->at(highest_idx) + my_primary->jtpt->at(second_highest_idx));
        my_hists->Aj[ibin]->Fill(Aj);

    }
  }

} /// is_dijets loop
          ///====== do dijet selection ends here ====//////

          int njets=0;
          ///// Now use jet trigger jets
          //// now loop over jets in the original events
          for(int j4i = 0; j4i < (int) my_primary->jtpt->size(); j4i++) {

            if(is_dijet){
              if ( (j4i != highest_idx) && (j4i != second_highest_idx) ) continue;  /// select 2 leading jet
              if(fabs(my_primary->jteta->at(j4i)) > 1.6 ) continue;
              if( my_primary->jtpt->at(highest_idx) <120 ) continue;
              if( my_primary->jtpt->at(highest_idx) >300 ) continue;
              if( my_primary->jtpt->at(second_highest_idx) <50 ) continue;
              if( my_primary->jtpt->at(second_highest_idx) >300 ) continue;

            }else{
              if(fabs(my_primary->jteta->at(j4i)) > 1.6 ) continue;
              if(my_primary->trackMax->at(j4i)/my_primary->jtpt->at(j4i) <=0.01) continue ;
              if( my_primary->jtpt->at(j4i) <120 ) continue;
              if( my_primary->jtpt->at(j4i) >300 ) continue;
            }


            if(is_dijet){
              if ( is_subleadingjet &&  j4i != second_highest_idx ) continue;  /// select subleading jets
              if ( is_leadingjet && j4i != highest_idx ) continue;  /// select leading jet
            }

            int ibin = 0;      int ibin2 = 0;  int ibin3=0;
            for (int centi=0;centi<nCBins;centi++){
              if (my_primary->hiBin >=CBins[centi] && my_primary->hiBin <CBins[centi+1])  ibin = centi;
            }

            for(int pti = 0; pti < nPtBins; pti++) {
              if (my_primary->jtpt->at(j4i) >=PtBins[pti] && my_primary->jtpt->at(j4i) < PtBins[pti+1])  ibin2 = pti ;
            } 

            float SumpTLeadJet = 0.;      float SumpTSubLeadJet = 0.;
            float SumPtJet = 0; 
            //// Get signal shape, but only if this is the first time ....

            double correctionFactors[4] = {0,0,0,0};
            float trkweight = 1.;
            float trkweight_bkg_num = 1.;
            float trkweight_bkg_er = 1.;
            float trkweight_MB_num = 1.;
            float trkweight_MB_den = 1.;


            if( n_MBevs_found == 1 ) {     //// only the first time 
            }   ///  is first MB event 
                //  std::cout << "past first MB event, mbevi: " << mb_evi << ", jet loop" << std::endl;


          ///==========     Done getting signal shape   ============/////

          ///======    Now do MB background   ========////
          double dir_eta= fabs(my_primary->jteta->at(j4i));
          double dir_phi = my_primary->jtphi->at(j4i);
          double SumPtTracks(0.0);
          double nbkg(0.0);
          double SumPt_bkg_particles_cone(0.0);
          int ntracks(0);
          const int NRBIN=6;
          const float RBIN[NRBIN+1]={0.0,0.05,0.1,0.15,0.2,0.25,0.3};
          std::vector<double> SumPt_bkg_particles(NRBIN);
          for(int vi = 0; vi < (int) SumPt_bkg_particles.size(); vi++) {
            SumPt_bkg_particles.at(vi) = 0.;
          }

          std::vector<double> multiplicity_bkg_particles_cone(NRBIN);
          for(int vi = 0; vi < (int) multiplicity_bkg_particles_cone.size(); vi++) {
            multiplicity_bkg_particles_cone.at(vi) = 0.;
          }

          //// Get MB background from sim-tracks

          int ngen=0;
          if( !is_data ) {

            /// closure with gen particles instead of simTracks ///

            for(int ipart = 0 ; ipart < (int) my_primary_MB->pt->size(); ipart++){ 
              double gen_pt = my_primary_MB->pt->at(ipart);
              double gen_phi = my_primary_MB->phi->at(ipart);
              double gen_eta = my_primary_MB->eta->at(ipart);
              int chg = my_primary_MB->chg->at(ipart);
              int sube = my_primary_MB->sube->at(ipart);

              if(chg==0) continue ;
              if(gen_pt<trkPtCut)continue ;

              if(TMath::Abs(gen_eta)>2.4)continue ;

              for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
                if (gen_pt >=TrkPtBins[trkpti] && gen_pt < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
              } /// trkpti loop

              double deta_gen = my_primary->jteta->at(j4i) - gen_eta;
              double dphi_gen = my_primary->jtphi->at(j4i) - gen_phi;

              if(sube==0) {   // selects pythia
                my_hists->hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen, wvz*wcen);
                my_hists->hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen+ (2.*TMath::Pi()), wvz*wcen);
                my_hists->hJetGenTrackBackground_pythia[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen- (2.*TMath::Pi()), wvz*wcen);
              }else{    // selects hydjet part
                my_hists->hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen, wvz*wcen);
                my_hists->hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen+ (2.*TMath::Pi()), wvz*wcen);
                my_hists->hJetGenTrackBackground_hydjet[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen- (2.*TMath::Pi()), wvz*wcen);
              }
            } // genparticle loop

           
            for(int ipart = 0 ; ipart < (int) my_primary_MB->pt->size() ; ipart++){ //sim track loop 
              double gen_pt = my_primary_MB->pt->at(ipart);
              double gen_phi = my_primary_MB->phi->at(ipart);
              double gen_eta = my_primary_MB->eta->at(ipart);
              if(gen_pt<trkPtCut)continue ;

              if(TMath::Abs(gen_eta)>2.4)continue ;
              double rr = mydeltaR(my_primary->jteta->at(j4i),my_primary->jtphi->at(j4i), gen_eta,gen_phi);

              for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
                if (gen_pt >=TrkPtBins[trkpti] && gen_pt < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
              } /// trkpti loop

              double deta_gen = my_primary->jteta->at(j4i) - gen_eta;
              double dphi_gen = my_primary->jtphi->at(j4i) - gen_phi;

              my_hists->hJetGenTrackBackground[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen, wvz*wcen);
              my_hists->hJetGenTrackBackground[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen+ (2.*TMath::Pi()), wvz*wcen);
              my_hists->hJetGenTrackBackground[ibin][ibin2][ibin3]->Fill(deta_gen,dphi_gen- (2.*TMath::Pi()), wvz*wcen);

              my_hists->GenTrkEta_bkg[ibin][ibin2][ibin3]->Fill(gen_eta,wvz*wcen);
              my_hists->GenTrkPhi_bkg[ibin][ibin2][ibin3]->Fill(gen_phi,wvz*wcen);

              my_hists->track_gen_bkg_pT_hist[ibin][ibin2]->Fill(gen_pt, wvz*wcen);
              my_hists->track_gen_bkg_eta_hist[ibin][ibin2]->Fill(gen_eta, wvz*wcen);

              my_hists->track_gen_bkg_phi_hist[ibin][ibin2]->Fill(gen_phi, wvz*wcen);

              if(rr<0.3){
                ngen++;
              }
            }/// simtrack loop 
          } /// isdata!



          /// Get MB sum track pT
          for(int tracks =0; tracks < (int) my_primary_MB->trkPt->size(); tracks++){
            if(fabs(my_primary_MB->trkEta->at(tracks))>2.4) continue;
            if (my_primary_MB->highPurity->at(tracks)==1) {
              if(my_primary_MB->trkPt->at(tracks)>trkPtCut) {

                double rr = mydeltaR(my_primary->jteta->at(j4i),my_primary->jtphi->at(j4i), my_primary_MB->trkEta->at(tracks),my_primary_MB->trkPhi->at(tracks));

                double ptParticle=my_primary_MB->trkPt->at(tracks);
                double phiParticle=my_primary_MB->trkPhi->at(tracks);
                double etaParticle=my_primary_MB->trkEta->at(tracks);

               Float_t trkweight_MB = (1-my_primary_MB->fake->at(tracks))/my_primary_MB->eff->at(tracks);

                for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) { 
                  if (my_primary_MB->trkPt->at(tracks) >=TrkPtBins[trkpti] && my_primary_MB->trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
                } /// trkpti loop


                my_hists->Bkg_phi[ibin][ibin2][ibin3]->Fill(my_primary_MB->trkPhi->at(tracks),wvz*wcen);
                my_hists->Bkg_eta[ibin][ibin2][ibin3]->Fill(my_primary_MB->trkEta->at(tracks),wvz*wcen);

                my_hists->Bkg_phi_weighted[ibin][ibin2][ibin3]->Fill(my_primary_MB->trkPhi->at(tracks),trkweight_MB*wvz*wcen);
                my_hists->Bkg_eta_weighted[ibin][ibin2][ibin3]->Fill(my_primary_MB->trkEta->at(tracks),trkweight_MB*wvz*wcen);

                double dphi_bkg = my_primary->jtphi->at(j4i) - my_primary_MB->trkPhi->at(tracks);
                double deta_bkg = my_primary->jteta->at(j4i) - my_primary_MB->trkEta->at(tracks);

                my_hists->hJetTrackBackground[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg, trkweight_MB*wvz*wcen);
                my_hists->hJetTrackBackground[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg + 2.*TMath::Pi(), trkweight_MB*wvz*wcen);
                my_hists->hJetTrackBackground[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg - 2.*TMath::Pi(), trkweight_MB*wvz*wcen);

                my_hists->hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg,wvz*wcen);
                my_hists->hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg + 2.*TMath::Pi(), wvz*wcen);
                my_hists->hJetTrackBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta_bkg,dphi_bkg - 2.*TMath::Pi(), wvz*wcen);

                float dphi_jet_track = dir_phi - my_primary_MB->trkPhi->at(tracks);
                my_hists->dPhi_jet_track[ibin][ibin2]->Fill(fabs(dphi_jet_track), wvz*wcen);

                my_hists->EMix_bkg_pt[ibin][ibin2]->Fill(ptParticle);
                my_hists->EMix_bkg_pt_weighted[ibin][ibin2]->Fill(ptParticle, trkweight_MB*wvz*wcen);

                my_hists->EMix_bkg_phi[ibin][ibin2]->Fill(phiParticle, wvz*wcen);
                my_hists->EMix_bkg_phi_weighted[ibin][ibin2]->Fill(phiParticle, trkweight_MB*wvz*wcen);

                my_hists->EMix_bkg_eta[ibin][ibin2]->Fill(etaParticle, wvz*wcen);
                my_hists->EMix_bkg_eta_weighted[ibin][ibin2]->Fill(etaParticle, trkweight_MB*wvz*wcen);


                if (rr<0.3){ 

                  my_hists->track_bkg_pT_hist[ibin][ibin2]->Fill(ptParticle, wvz*wcen);
                  my_hists->track_bkg_pT_hist_weighted[ibin][ibin2]->Fill(ptParticle, trkweight_MB*wvz*wcen);

                  my_hists->track_bkg_phi_hist[ibin][ibin2]->Fill(phiParticle, wvz*wcen);
                  my_hists->track_bkg_phi_hist_weighted[ibin][ibin2]->Fill(phiParticle, trkweight_MB*wvz*wcen);

                  my_hists->track_bkg_eta_hist[ibin][ibin2]->Fill(etaParticle, wvz*wcen);
                  my_hists->track_bkg_eta_hist_weighted[ibin][ibin2]->Fill(etaParticle, trkweight_MB*wvz*wcen);

                  my_hists->dN_tracks[ibin][ibin2]->Fill(TMath::Log(my_primary->jtpt->at(j4i)/ptParticle), wvz*wcen);
                  ntracks++;
                  SumPtTracks += ptParticle*trkweight_MB_den;
                  SumPt_bkg_particles_cone += ptParticle*trkweight_MB_den;
                }
              }
            }   /// end MB track loop
          }


          for(int tracks =0; tracks < (int) my_primary_MB->trkPt->size(); tracks++){
            if(fabs(my_primary_MB->trkEta->at(tracks))>2.4) continue;
            if (my_primary_MB->highPurity->at(tracks)==1) {
              if(my_primary_MB->trkPt->at(tracks)>trkPtCut) {

                double rr = mydeltaR(my_primary->jteta->at(j4i),my_primary->jtphi->at(j4i), my_primary_MB->trkEta->at(tracks),my_primary_MB->trkPhi->at(tracks));

                double ptParticle=my_primary_MB->trkPt->at(tracks);
               Float_t trkweight_MB_num = (1-my_primary_MB->fake->at(tracks))/my_primary_MB->eff->at(tracks);

                for(int ir=0;ir<NRBIN;ir++){
                  if(RBIN[ir]<rr && rr <= RBIN[ir+1]){
                    SumPt_bkg_particles[ir] += ptParticle*trkweight_MB_num;
                    multiplicity_bkg_particles_cone[ir] += ntracks;
                  }
                }
              }
            }  /// end MB track loop
          }


          ///=================================================================================================
          ///------------------------ Filling Background shape histograms  -----------------------------------
          ///=================================================================================================
          std::vector<double> DifferentialPt_bkg(NRBIN);
          std::vector<double> tot_sumpt_rbin_bkg(NRBIN);
          for(int ir = 0; ir < NRBIN; ir++) {
            DifferentialPt_bkg[ir] = 0;
            tot_sumpt_rbin_bkg[ir] = 0;
          }

          for(int ir=0;ir<NRBIN;ir++){
            DifferentialPt_bkg[ir] = SumPt_bkg_particles[ir]; 
          }

          float LeadSumPt_bkg_particles_cone = 0.;

          if(my_primary->jtpt->at(j4i)>0.0001){
            my_hists->bkg_count[ibin][ibin2] += wvz;
            nbkg += wvz;
          }

          for(int ir=0;ir<NRBIN;ir++){
            float fill_val = 0.05*(ir+1)-0.025;

            LeadSumPt_bkg_particles_cone += DifferentialPt_bkg[ir];

            if(my_primary->jtpt->at(j4i)>0.0001){
              my_hists->SumPt_bkg_only[ibin][ibin2]->Fill(LeadSumPt_bkg_particles_cone, wvz*wcen);
              my_hists->JetShapeDiffParticles_bkg_1D[ibin][ibin2]->Fill(fill_val,(DifferentialPt_bkg[ir] / my_primary->jtpt->at(j4i))*wvz*wcen);
              my_hists->SumPt_only[ibin][ibin2]->Fill(double(SumpTLeadJet), wvz*wcen);

            }
          }

          my_hists->emix_jetpt_test[ibin][ibin2]->Fill(my_primary->jtpt->at(j4i),wvz*wcen);
          my_hists->emix_jetpt[ibin][ibin2]->Fill(double(nbkg), wvz*wcen);


          ///=================================================================================================
          ///-------------- sanity check for LeadingJetPt and sumPt of candidates ----------------------------
          ///=================================================================================================

          } /// main jet loop
      }  /// looping over MB events until you find 10

      my_hists->NumberOfMatches->Fill(n_MBevs_found);

      if( evi % 500 == 0 ) my_hists->Write();
      }  ///event loop

      std::cout << "done with file " << fi << std::endl;

  }  /// file loop
}



void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug)
{
  ifstream file_stream(file_of_names);
  std::string line;
  my_file_names.clear();
  if( debug ) std::cout << "Open file " << file_of_names << " to extract files to run over" << std::endl;
  if( file_stream.is_open() ) {
    if( debug ) std::cout << "Opened " << file_of_names << " for reading" << std::endl;
    int line_num = 0;
    while( !file_stream.eof() ) {
      getline(file_stream, line);
      if( debug ) std::cout << line_num << ": " << line << std::endl;
      TString tstring_line(line);
      if( tstring_line.CompareTo("", TString::kExact) != 0 ) my_file_names.push_back(tstring_line);
      line_num++;
    }
  } else {
    std::cout << "Error, could not open " << file_of_names << " for reading" << std::endl;
    assert(0);
  }
}



//void GetFilesOfFileNames(std::vector<TString> &files_of_file_names, std::vector<float> &Xsections) ////getting the data/simulation based on whatever dataset you pick (arg)
  void GetFilesOfFileNames(std::vector<TString> &files_of_file_names, std::vector<float> &Xsections, std::vector<double> &assumed_n_evt) ////getting the data/simulation based on whatever dataset you pick (arg)
{
  files_of_file_names.clear();   Xsections.clear();

  if( dataset_type_code == e_HydJet30 ) {
    files_of_file_names.push_back("Hydjet_Pt30.txt");
    Xsections.push_back(1.079e-02);
    assumed_n_evt.push_back(109918);

  }else if( dataset_type_code == e_HydJet50 ) {
    files_of_file_names.push_back("Hydjet_Pt50.txt");
    Xsections.push_back(1.021e-03 - 9.913E-05);
    assumed_n_evt.push_back(112166);

  }else if( dataset_type_code == e_HydJet80 ) {
    files_of_file_names.push_back("Hydjet_Pt80.txt");
    Xsections.push_back(9.913E-05 - 3.069E-05);
    assumed_n_evt.push_back(146162);

  }else if( dataset_type_code == e_HydJet100 ) {
    files_of_file_names.push_back("Hydjet_Pt100.txt");
    Xsections.push_back(3.069E-05 - 1.128E-05);
    assumed_n_evt.push_back(593463);

  }else if( dataset_type_code == e_HydJet120 ) {
    files_of_file_names.push_back("Hydjet_Pt120.txt");
    Xsections.push_back(1.128E-05 - 1.470E-06);
    assumed_n_evt.push_back(141263);

  }else if( dataset_type_code == e_HydJet170 ) {
    files_of_file_names.push_back("Hydjet_Pt170.txt");
    Xsections.push_back(1.470E-06 - 5.310E-07);
    assumed_n_evt.push_back(62944);

  }else if( dataset_type_code == e_HydJet200 ) {
    files_of_file_names.push_back("Hydjet_Pt200.txt");
    Xsections.push_back(5.310E-07 - 1.192E-07);
    assumed_n_evt.push_back(54952);

  }else if( dataset_type_code == e_HydJet250 ) {
    files_of_file_names.push_back("Hydjet_Pt250.txt");
    Xsections.push_back(1.192E-07 - 3.176E-08);
  }else if( dataset_type_code == e_HydJet300 ) {
    files_of_file_names.push_back("Hydjet_Pt300.txt");
    Xsections.push_back(3.176E-08);
    assumed_n_evt.push_back(53009);

  }else if( dataset_type_code == e_Data2011 ) {
    files_of_file_names.push_back("Data2011.txt");
  }else {
    std::cout << "I don't understand the fileset" << std::endl;
    assert(0);
  }
}

float GetDatasetWeight(double n_evt_raw, double Xsection)
{
  assert( n_evt_raw > 0.5 );
 double wt  = Xsection / n_evt_raw;
 
  return wt;
}

