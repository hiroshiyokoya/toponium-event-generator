{
gROOT->Reset();

#include <fstream>

gROOT->SetStyle("Plain");
gStyle->SetTitleBorderSize(2);
gStyle->SetFrameFillColor(0);
gStyle->SetCanvasColor(0);
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

TCanvas *c1 = new TCanvas("c1","c1",0,0,1400,800);
c1->Divide(5,3);
int i = 1;

c1->cd(i++);
TFile *f150 = new TFile ("./mtt_mt150.root");
TH2D *h_mean_150 = (TH2D*)f150->Get("mean vs mtt");
h_mean_150->SetTitle("Mean M=150 GeV");
h_mean_150->Draw();
c1->cd(i++);
TFile *f160 = new TFile ("./mtt_mt160.root");
TH2D *h_mean_160 = (TH2D*)f160->Get("mean vs mtt");
h_mean_160->SetTitle("Mean M=160 GeV");
h_mean_160->Draw();
c1->cd(i++);
TFile *f173 = new TFile ("./mtt_mt173.root");
TH2D *h_mean_173 = (TH2D*)f173->Get("mean vs mtt");
h_mean_173->SetTitle("Mean M=173 GeV");
h_mean_173->Draw();
c1->cd(i++);
TFile *f180 = new TFile ("./mtt_mt180.root");
TH2D *h_mean_180 = (TH2D*)f180->Get("mean vs mtt");
h_mean_180->SetTitle("Mean M=180 GeV");
h_mean_180->Draw();
c1->cd(i++);
TFile *f190 = new TFile ("./mtt_mt190.root");
TH2D *h_mean_190 = (TH2D*)f190->Get("mean vs mtt");
h_mean_190->SetTitle("Mean M=190 GeV");
h_mean_190->Draw();
//c1->cd(i++);
//TFile *f200 = new TFile ("./mtt_mt200.root");
//TH2D *h_mean_200 = (TH2D*)f200->Get("mean vs mtt");
//h_mean_200->SetTitle("Mean M=200 GeV");
//h_mean_200->Draw();

c1->cd(i++);
TH2D *h_min_150 = (TH2D*)f150->Get("min vs mtt");
h_min_150->SetTitle("Min M=150 GeV");
h_min_150->Draw();
c1->cd(i++);
TH2D *h_min_160 = (TH2D*)f160->Get("min vs mtt");
h_min_160->SetTitle("Min M=160 GeV");
h_min_160->Draw();
c1->cd(i++);
TH2D *h_min_173 = (TH2D*)f173->Get("min vs mtt");
h_min_173->SetTitle("Min M=173 GeV");
h_min_173->Draw();
c1->cd(i++);
TH2D *h_min_180 = (TH2D*)f180->Get("min vs mtt");
h_min_180->SetTitle("Min M=180 GeV");
h_min_180->Draw();
c1->cd(i++);
TH2D *h_min_190 = (TH2D*)f190->Get("min vs mtt");
h_min_190->SetTitle("Min M=190 GeV");
h_min_190->Draw();
//c1->cd(i++);
//TH2D *h_min_200 = (TH2D*)f200->Get("min vs mtt");
//h_min_200->SetTitle("Min M=200 GeV");
//h_min_200->Draw();

c1->cd(i++);
TH2D *h_max_150 = (TH2D*)f150->Get("max vs mtt");
h_max_150->SetTitle("Max M=150 GeV");
h_max_150->Draw();
c1->cd(i++);
TH2D *h_max_160 = (TH2D*)f160->Get("max vs mtt");
h_max_160->SetTitle("Max M=160 GeV");
h_max_160->Draw();
c1->cd(i++);
TH2D *h_max_173 = (TH2D*)f173->Get("max vs mtt");
h_max_173->SetTitle("Max M=173 GeV");
h_max_173->Draw();
c1->cd(i++);
TH2D *h_max_180 = (TH2D*)f180->Get("max vs mtt");
h_max_180->SetTitle("Max M=180 GeV");
h_max_180->Draw();
c1->cd(i++);
TH2D *h_max_190 = (TH2D*)f190->Get("max vs mtt");
h_max_190->SetTitle("Max M=190 GeV");
h_max_190->Draw();
//c1->cd(i++);
//TH2D *h_max_200 = (TH2D*)f200->Get("max vs mtt");
//h_max_200->SetTitle("Max M=200 GeV");
//h_max_200->Draw();


}
