//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul 18 00:13:58 2010 by ROOT version 5.22/00d
// from TTree mytree/
// found on file: murestree_100_MC_36Y_V10.root
//////////////////////////////////////////////////////////

#ifndef plotResHeader_h
#define plotResHeader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream> 
#include <TString.h>
#include <TCanvas.h>

class plotResHeader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           evt;
   Int_t           ls;
   Int_t           nmu;
   Bool_t          tevvalid[6];   //[nmu]
   Float_t         tevpt[6];   //[nmu]
   Float_t         tevipt[6];   //[nmu]
   Int_t           tevhits[6];   //[nmu]
   Float_t         tevres[6];   //[nmu]
   Bool_t          anewvalid[6];   //[nmu]
   Float_t         anewpt[6];   //[nmu]
   Float_t         anewipt[6];   //[nmu]
   Int_t           anewhits[6];   //[nmu]
   Float_t         anewres[6];   //[nmu]
   Bool_t          oddvalid[6];   //[nmu]
   Float_t         oddpt[6];   //[nmu]
   Float_t         oddipt[6];   //[nmu]
   Int_t           oddhits[6];   //[nmu]
   Float_t         oddres[6];   //[nmu]
   Bool_t          evevalid[6];   //[nmu]
   Float_t         evept[6];   //[nmu]
   Float_t         eveipt[6];   //[nmu]
   Int_t           evehits[6];   //[nmu]
   Float_t         everes[6];   //[nmu]
   Float_t         mres[6];   //[nmu]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_nmu;   //!
   TBranch        *b_tevvalid;   //!
   TBranch        *b_tevpt;   //!
   TBranch        *b_tevipt;   //!
   TBranch        *b_tevhits;   //!
   TBranch        *b_tevres;   //!
   TBranch        *b_anewvalid;   //!
   TBranch        *b_anewpt;   //!
   TBranch        *b_anewipt;   //!
   TBranch        *b_anewhits;   //!
   TBranch        *b_anewres;   //!
   TBranch        *b_oddvalid;   //!
   TBranch        *b_oddpt;   //!
   TBranch        *b_oddipt;   //!
   TBranch        *b_oddhits;   //!
   TBranch        *b_oddres;   //!
   TBranch        *b_evevalid;   //!
   TBranch        *b_evept;   //!
   TBranch        *b_eveipt;   //!
   TBranch        *b_evehits;   //!
   TBranch        *b_everes;   //!
   TBranch        *b_mres;   //!

   plotResHeader (TString fileName,float ptmuon,TString gt);
   virtual ~plotResHeader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   float thept;
   TString thegt;
   TString strpt;
  
};

#endif

#ifdef plotResHeader_cxx
plotResHeader::plotResHeader(TString fileName,float ptmuon, TString gt)
{
    std::cout << "Opening " << fileName << std::endl;
    thept = ptmuon;
    thegt = gt;
    char buffer[50];
    sprintf(buffer,"%0.0f",thept);
    strpt = buffer;
    TFile *theFile  = TFile::Open (fileName);
    TTree *tree     = (TTree*)theFile->Get("mytree");
    Init(tree);

}

plotResHeader::~plotResHeader()
{
   if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t plotResHeader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t plotResHeader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void plotResHeader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("nmu", &nmu, &b_nmu);
   fChain->SetBranchAddress("tevvalid", tevvalid, &b_tevvalid);
   fChain->SetBranchAddress("tevpt", tevpt, &b_tevpt);
   fChain->SetBranchAddress("tevipt", tevipt, &b_tevipt);
   fChain->SetBranchAddress("tevhits", tevhits, &b_tevhits);
   fChain->SetBranchAddress("tevres", tevres, &b_tevres);
   fChain->SetBranchAddress("anewvalid", anewvalid, &b_anewvalid);
   fChain->SetBranchAddress("anewpt", anewpt, &b_anewpt);
   fChain->SetBranchAddress("anewipt", anewipt, &b_anewipt);
   fChain->SetBranchAddress("anewhits", anewhits, &b_anewhits);
   fChain->SetBranchAddress("anewres", anewres, &b_anewres);
   fChain->SetBranchAddress("oddvalid", oddvalid, &b_oddvalid);
   fChain->SetBranchAddress("oddpt", oddpt, &b_oddpt);
   fChain->SetBranchAddress("oddipt", oddipt, &b_oddipt);
   fChain->SetBranchAddress("oddhits", oddhits, &b_oddhits);
   fChain->SetBranchAddress("oddres", oddres, &b_oddres);
   fChain->SetBranchAddress("evevalid", evevalid, &b_evevalid);
   fChain->SetBranchAddress("evept", evept, &b_evept);
   fChain->SetBranchAddress("eveipt", eveipt, &b_eveipt);
   fChain->SetBranchAddress("evehits", evehits, &b_evehits);
   fChain->SetBranchAddress("everes", everes, &b_everes);
   fChain->SetBranchAddress("mres", mres, &b_mres);
   Notify();
}

Bool_t plotResHeader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void plotResHeader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t plotResHeader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef plotResHeader_cxx
