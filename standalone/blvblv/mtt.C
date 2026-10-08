{
gROOT->Reset();

#include <fstream>

//TCanvas *c1 = new TCanvas("c1","c1",0,0,500,500);
TFile *fout = new TFile("mtt.root","recreate");

gROOT->SetStyle("Plain");
gStyle->SetTitleBorderSize(2);
gStyle->SetFrameFillColor(0);
gStyle->SetCanvasColor(0);
gStyle->SetPalette(1);

int bin=50,bmin=300,bmax=800;
TH1D *h_mean = new TH1D ("mean","mean",bin,bmin,bmax);
TH1D *h_min  = new TH1D ( "min", "min",bin,bmin,bmax);
TH1D *h_max  = new TH1D ( "max", "max",bin,bmin,bmax);
TH1D *h_non  = new TH1D ( "non", "non",bin,bmin,bmax);
TH1D *h_mtt  = new TH1D ( "mtt", "mtt",bin,bmin,bmax);
TH1D *h_sol  = new TH1D ( "sol", "sol",bin,bmin,bmax);

TH2D *h_minmtt = new TH2D ("min vs mtt","min vs mtt",
			   bin,bmin,bmax,bin,bmin,bmax);
TH2D *h_maxmtt = new TH2D ("max vs mtt","max vs mtt",
			   bin,bmin,bmax,bin,bmin,bmax);
TH2D *h_meanmtt = new TH2D ("mean vs mtt","mean vs mtt",
			    bin,bmin,bmax,bin,bmin,bmax);
h_minmtt->SetOption("COLZ");
h_maxmtt->SetOption("COLZ");
h_meanmtt->SetOption("COLZ");

ifstream fin("fort.2");
double mtt,mean,min,max;
while(fin>>mtt>>mean>>min>>max) {
  h_min->Fill(min);
  h_max->Fill(max);
  h_mean->Fill(mean);
  h_minmtt->Fill(mtt,min);
  h_maxmtt->Fill(mtt,max);
  h_meanmtt->Fill(mtt,mean);
}

ifstream fin3("fort.3");
double non;
while(fin3>>non) {
  h_non->Fill(non);
}

ifstream fin4("fort.4");
while(fin4>>mtt) {
  h_mtt->Fill(mtt);
}

ifstream fin10("fort.10");
double sol;
while(fin10>>sol) {
  h_sol->Fill(sol);
}

//h_mtt->SetLineColor(1);
//h_mtt->SetLineWidth(3);
//h_mtt->Scale(1./h_mtt->GetEntries());
//h_mtt->Draw();
//h_min->SetLineColor(2);
//h_min->Scale(1./h_min->GetEntries());
//h_min->Draw("SAME");
//h_max->SetLineColor(3);
//h_max->Scale(1./h_max->GetEntries());
//h_max->Draw("SAME");
//h_sol->SetLineColor(4);
//h_sol->Scale(1./h_sol->GetEntries());
//h_sol->Draw("SAME");

TH1D *h_frac = new TH1D(*h_non);
h_frac->Divide(h_mtt);
h_frac->SetName("fraction");
h_frac->SetTitle("fraction");
h_frac->SetMinimum(0);

fout->Write();
fout->Close();
//c1->Update();

}
