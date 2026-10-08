{
gROOT->Reset();

#include <fstream>

TCanvas *c1 = new TCanvas("c1","c1",0,0,500,500);

gROOT->SetStyle("Plain");
gStyle->SetTitleBorderSize(2);
gStyle->SetFrameFillColor(0);
gStyle->SetCanvasColor(0);
gStyle->SetPalette(1);

TH2I *h2 = new TH2I ("NSOL vs NCOMB","NSOL vs NCOMB",6,0,6,6,0,6);

ifstream fin("fort.1");
int nsol,ncomb;
while(fin>>nsol>>ncomb) {
  h2->Fill(nsol,ncomb);
}
h2->SetOption("Text");
h2->Draw();

c1->Update();

}
