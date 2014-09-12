//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 30 11:51:56 2013 by ROOT version 5.32/00
// from TTree mixing_tree_MB/
// found on file: ../input_ntuple/minbias_mini_ntuple_p0.root
//////////////////////////////////////////////////////////

#ifndef mixing_tree_MB_h
#define mixing_tree_MB_h








#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <iostream>
using namespace std;


// Fixed size dimensions of array or collections stored in the TTree if any.

class mixing_tree_MB {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           mult;
   vector<float>   *pt;
   vector<float>   *phi;
   vector<float>   *eta;
   vector<int>     *chg;
   Int_t           nTrk;
   Int_t           nParticle;
   vector<float>   *pPt;
   vector<float>   *pPhi;
   vector<float>   *pEta;
   vector<float>   *pNRec;
   vector<float>   *mtrkAlgo;
   vector<int>     *mtrkQual;
   vector<float>   *trkEta;
   vector<float>   *trkPhi;
   vector<float>   *trkPt;
   vector<float>   *trkAlgo;
   vector<float>   *highPurity;
   vector<float>   *vz;
   vector<float>   *trkweight;
   Int_t           pHBHENoiseFilter;
   Int_t           pcollisionEventSelection;
   Int_t           HLT_HIMinBiasHfOrBSC_v1;
   Int_t           hiBin;
   vector<float>   *jteta;
   vector<float>   *jtphi;
   vector<float>   *jtpt;
   Float_t         pthat;
   vector<float>   *refeta;
   vector<float>   *refphi;
   vector<float>   *refpt;
   vector<float>   *trackMax;
   Float_t         hiEvtPlanes;

   // List of branches
   TBranch        *b_mult;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_chg;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_nParticle;   //!
   TBranch        *b_pPt;   //!
   TBranch        *b_pPhi;   //!
   TBranch        *b_pEta;   //!
   TBranch        *b_pNRec;   //!
   TBranch        *b_mtrkAlgo;   //!
   TBranch        *b_mtrkQual;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkAlgo;   //!
   TBranch        *b_highPurity;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_trkweight;   //!
   TBranch        *b_pHBHENoiseFilter;   //!
   TBranch        *b_pcollisionEventSelection;   //!
   TBranch        *b_HLT_HIMinBiasHfOrBSC_v1;   //!
   TBranch        *b_hiBin;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_refeta;   //!
   TBranch        *b_refphi;   //!
   TBranch        *b_refpt;   //!
   TBranch        *b_trackMax;   //!
   TBranch        *b_hiEvtPlanes;   //!

   mixing_tree_MB(TTree *tree=0);
   virtual ~mixing_tree_MB();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mixing_tree_MB_cxx
mixing_tree_MB::mixing_tree_MB(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../input_ntuple/minbias_mini_ntuple_p0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../input_ntuple/minbias_mini_ntuple_p0.root");
      }
      f->GetObject("mixing_tree_MB",tree);

   }
   Init(tree);
}

mixing_tree_MB::~mixing_tree_MB()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mixing_tree_MB::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mixing_tree_MB::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void mixing_tree_MB::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pt = 0;
   phi = 0;
   eta = 0;
   chg = 0;
   pPt = 0;
   pPhi = 0;
   pEta = 0;
   pNRec = 0;
   mtrkAlgo = 0;
   mtrkQual = 0;
   trkEta = 0;
   trkPhi = 0;
   trkPt = 0;
   trkAlgo = 0;
   highPurity = 0;
   vz = 0;
   trkweight = 0;
   jteta = 0;
   jtphi = 0;
   jtpt = 0;
   refeta = 0;
   refphi = 0;
   refpt = 0;
   trackMax = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mult", &mult, &b_mult);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("chg", &chg, &b_chg);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("nParticle", &nParticle, &b_nParticle);
   fChain->SetBranchAddress("pPt", &pPt, &b_pPt);
   fChain->SetBranchAddress("pPhi", &pPhi, &b_pPhi);
   fChain->SetBranchAddress("pEta", &pEta, &b_pEta);
   fChain->SetBranchAddress("pNRec", &pNRec, &b_pNRec);
   fChain->SetBranchAddress("mtrkAlgo", &mtrkAlgo, &b_mtrkAlgo);
   fChain->SetBranchAddress("mtrkQual", &mtrkQual, &b_mtrkQual);
   fChain->SetBranchAddress("trkEta", &trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", &trkPhi, &b_trkPhi);
   fChain->SetBranchAddress("trkPt", &trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkAlgo", &trkAlgo, &b_trkAlgo);
   fChain->SetBranchAddress("highPurity", &highPurity, &b_highPurity);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("trkweight", &trkweight, &b_trkweight);
   fChain->SetBranchAddress("pHBHENoiseFilter", &pHBHENoiseFilter, &b_pHBHENoiseFilter);
   fChain->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection, &b_pcollisionEventSelection);
   fChain->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1", &HLT_HIMinBiasHfOrBSC_v1, &b_HLT_HIMinBiasHfOrBSC_v1);
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("jteta", &jteta, &b_jteta);
   fChain->SetBranchAddress("jtphi", &jtphi, &b_jtphi);
   fChain->SetBranchAddress("jtpt", &jtpt, &b_jtpt);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("refeta", &refeta, &b_refeta);
   fChain->SetBranchAddress("refphi", &refphi, &b_refphi);
   fChain->SetBranchAddress("refpt", &refpt, &b_refpt);
   fChain->SetBranchAddress("trackMax", &trackMax, &b_trackMax);
   fChain->SetBranchAddress("hiEvtPlanes", &hiEvtPlanes, &b_hiEvtPlanes);
   Notify();
}

Bool_t mixing_tree_MB::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mixing_tree_MB::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mixing_tree_MB::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mixing_tree_MB_cxx
