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
gStyle->SetOptTitle(1);
gStyle->SetStatX(.315);
gStyle->SetStatY(.31);

int bin = 50;
double min=166, max=176;
//double min=0,max=1000;
//int bin = 5000;

double yoff = 1.4;
double lm=.115, rm = .115;

TCanvas *c1 = new TCanvas("c1","mt1 vs mt2",0,0,500,500);
gPad->SetLeftMargin(lm);
gPad->SetRightMargin(rm);

int k;
int NB_both[10]={},NB_eith[10]={},NG_both[10]={},NG_eith[10]={};

TH2F *h1 = new TH2F ("Full result","m_{tt}#leq 350 GeV, LHC /s = 14 TeV",
		     bin,min,max,bin,min,max);
ifstream fin("./fort.34");
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

c1->Update();

}
