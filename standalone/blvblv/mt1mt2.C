{
gROOT->Reset();

gROOT->SetStyle("Plain");
gStyle->SetTitleBorderSize(2);
gStyle->SetFrameFillColor(0);
gStyle->SetCanvasColor(0);
gStyle->SetPalette(1);


#include <fstream>

double mt=173.,Gt=1.4911;

/** **/
gStyle->SetOptStat(1101);
gStyle->SetOptTitle(0);
gStyle->SetStatX(.315);
gStyle->SetStatY(.31);

int bin = 50;
double min=168, max=178;
//double min=0,max=1000;
//int bin = 5000;

double yoff = 1.4;
double lm=.115, rm = .115;

TCanvas *c1 = new TCanvas("c1","mt1 vs mt2",0,0,1000,500);
c1->Divide(2,1);
c1->cd(1);
gPad->SetLeftMargin(lm);
gPad->SetRightMargin(rm);

int k;
int NB_both[10]={},NB_eith[10]={},NG_both[10]={},NG_eith[10]={};

TH2F *h1 = new TH2F ("LO","LHC /s = 14 TeV",bin,min,max,bin,min,max);
ifstream fin("RBorn/370/fort.33");
//ifstream fin("RBorn/fort.33");
double mt1,mt2;
while(fin>>mt1>>mt2) {
  h1->Fill(mt1,mt2);
  k=1;
  if (fabs(mt1-mt)<=k*Gt && fabs(mt2-mt)<=k*Gt) NB_both[k]++;
  if (fabs(mt1-mt)<=k*Gt || fabs(mt2-mt)<=k*Gt) NB_eith[k]++;
  k=2;
  if (fabs(mt1-mt)<=k*Gt && fabs(mt2-mt)<=k*Gt) NB_both[k]++;
  if (fabs(mt1-mt)<=k*Gt || fabs(mt2-mt)<=k*Gt) NB_eith[k]++;
  k=3;
  if (fabs(mt1-mt)<=k*Gt && fabs(mt2-mt)<=k*Gt) NB_both[k]++;
  if (fabs(mt1-mt)<=k*Gt || fabs(mt2-mt)<=k*Gt) NB_eith[k]++;
  k=4;
  if (fabs(mt1-mt)<=k*Gt && fabs(mt2-mt)<=k*Gt) NB_both[k]++;
  if (fabs(mt1-mt)<=k*Gt || fabs(mt2-mt)<=k*Gt) NB_eith[k]++;
  k=5;
  if (fabs(mt1-mt)<=k*Gt && fabs(mt2-mt)<=k*Gt) NB_both[k]++;
  if (fabs(mt1-mt)<=k*Gt || fabs(mt2-mt)<=k*Gt) NB_eith[k]++;
}
fin.close();
h1->SetXTitle("m_{bW^{+}} [GeV]");
h1->SetYTitle("m_{#bar{b}W^{-}} [GeV]");
h1->GetYaxis()->SetTitleOffset(yoff);
h1->Scale(10000./h1->Integral());
h1->Draw("COLZ");
//h1->SetContour(15);

//TPaletteAxis *pal1 = (TPaletteAxis *)h1->GetListOfFunctions()->FindObject("palette");
//pal1->SetX1NDC(.85);
//pal1->SetX2NDC(.9);
//pal1->Draw();


//TCanvas *c2 = new TCanvas("c2","mt1 vs mt2",650,0,550,500);
c1->cd(2);
gPad->SetLeftMargin(lm);
gPad->SetRightMargin(rm);

TH2F *h2 = new TH2F ("Full result","", bin,min,max,bin,min,max);
ifstream fin("Rgrn/370/fort.34");
//ifstream fin("Rgrn/fort.34");
double mt1,mt2;
while(fin>>mt1>>mt2) {
//for (Int_t i=0;i<100000;i++){
//  fin>>mt1>>mt2;
  h2->Fill(mt1,mt2);
  k=1;
  if (fabs(mt1-mt)<=k*Gt && fabs(mt2-mt)<=k*Gt) NG_both[k]++;
  if (fabs(mt1-mt)<=k*Gt || fabs(mt2-mt)<=k*Gt) NG_eith[k]++;
  k=2;
  if (fabs(mt1-mt)<=k*Gt && fabs(mt2-mt)<=k*Gt) NG_both[k]++;
  if (fabs(mt1-mt)<=k*Gt || fabs(mt2-mt)<=k*Gt) NG_eith[k]++;
  k=3;
  if (fabs(mt1-mt)<=k*Gt && fabs(mt2-mt)<=k*Gt) NG_both[k]++;
  if (fabs(mt1-mt)<=k*Gt || fabs(mt2-mt)<=k*Gt) NG_eith[k]++;
  k=4;
  if (fabs(mt1-mt)<=k*Gt && fabs(mt2-mt)<=k*Gt) NG_both[k]++;
  if (fabs(mt1-mt)<=k*Gt || fabs(mt2-mt)<=k*Gt) NG_eith[k]++;
  k=5;
  if (fabs(mt1-mt)<=k*Gt && fabs(mt2-mt)<=k*Gt) NG_both[k]++;
  if (fabs(mt1-mt)<=k*Gt || fabs(mt2-mt)<=k*Gt) NG_eith[k]++;
}
fin.close();

h2->SetXTitle("m_{bW^{+}} [GeV]");
h2->SetYTitle("m_{#bar{b}W^{-}} [GeV]");
h2->GetYaxis()->SetTitleOffset(yoff);
h2->Scale(10000./h2->Integral());
h2->Draw("COLZ");
//h2->SetContour(15);

//double zmax = 2100;
double zmax = h1->GetMaximum();
//zmax=600.;
h1->SetMaximum(zmax);
h2->SetMaximum(zmax);
//c1_1->SetLogz();
//c1_2->SetLogz();

std::cout<<NG_both[1]<<", "<<NG_both[2]<<", "<<NG_both[3]<<", "
<<NG_both[4]<<", "<<NG_both[5]
<<std::endl;
std::cout<<NB_both[1]<<", "<<NB_both[2]<<", "<<NB_both[3]<<", "
<<NB_both[4]<<", "<<NB_both[5]
<<std::endl;
std::cout<<NG_eith[1]<<", "<<NG_eith[2]<<", "<<NG_eith[3]<<", "
<<NG_eith[4]<<", "<<NG_eith[5]
<<std::endl;
std::cout<<NB_eith[1]<<", "<<NB_eith[2]<<", "<<NB_eith[3]<<", "
<<NB_eith[4]<<", "<<NB_eith[5]
<<std::endl;

//h1->GetXaxis()->SetRangeUser(168.,176.);
//h1->GetYaxis()->SetRangeUser(168.,176.);
//h2->GetXaxis()->SetRangeUser(168.,176.);
//h2->GetYaxis()->SetRangeUser(168.,176.);

c1->Update();

//c1->Print("lt360_ggoct.eps");
//c1->Print("lt370.eps");
//c1->Print("lt360_7l.root");

}
