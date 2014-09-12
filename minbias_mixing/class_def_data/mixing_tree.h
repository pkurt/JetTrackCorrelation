//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun  5 01:36:02 2014 by ROOT version 5.32/00
// from TTree mixing_tree/
// found on file: /afs/cern.ch/work/h/htrauger/public/OfficialReco_PbPbData_forMixingOnly_p0.root
//////////////////////////////////////////////////////////

#ifndef mixing_tree_h
#define mixing_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <iostream>
using namespace std;

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class mixing_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<float>   *vz;
   Int_t           pHBHENoiseFilter;
   Int_t           pcollisionEventSelection;
   Int_t           HLT_HIJet80_v1;
   Int_t           hiBin;
   vector<float>   *jteta;
   vector<float>   *jtphi;
   vector<float>   *jtpt;
   Float_t         pthat;
   vector<float>   *trackMax;

   // List of branches
   TBranch        *b_vz;   //!
   TBranch        *b_pHBHENoiseFilter;   //!
   TBranch        *b_pcollisionEventSelection;   //!
   TBranch        *b_HLT_HIJet80_v1;   //!
   TBranch        *b_hiBin;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_trackMax;   //!

   mixing_tree(TTree *tree=0);
   virtual ~mixing_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mixing_tree_cxx
mixing_tree::mixing_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/work/h/htrauger/public/OfficialReco_PbPbData_forMixingOnly_p0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/afs/cern.ch/work/h/htrauger/public/OfficialReco_PbPbData_forMixingOnly_p0.root");
      }
      f->GetObject("mixing_tree",tree);

   }
   Init(tree);
}

mixing_tree::~mixing_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mixing_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mixing_tree::LoadTree(Long64_t entry)
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

void mixing_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vz = 0;
   jteta = 0;
   jtphi = 0;
   jtpt = 0;
   trackMax = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("pHBHENoiseFilter", &pHBHENoiseFilter, &b_pHBHENoiseFilter);
   fChain->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection, &b_pcollisionEventSelection);
   fChain->SetBranchAddress("HLT_HIJet80_v1", &HLT_HIJet80_v1, &b_HLT_HIJet80_v1);
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("jteta", &jteta, &b_jteta);
   fChain->SetBranchAddress("jtphi", &jtphi, &b_jtphi);
   fChain->SetBranchAddress("jtpt", &jtpt, &b_jtpt);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("trackMax", &trackMax, &b_trackMax);
   Notify();
}

Bool_t mixing_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mixing_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mixing_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef mixing_tree_cxx
