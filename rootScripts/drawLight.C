void drawLight(string fnm="../o_avgModel_ideal23v2_mirror_offset_0.00_Nmodels_7.root"){
  TFile *fin=TFile::Open(fnm.c_str(),"READ");

  TH1D *hl=(TH1D*)fin->Get("peL_pos_0");
  TH1D *hr=(TH1D*)fin->Get("peR_pos_0");
  TH1D *xpos=(TH1D*)fin->Get("x_pos");

  gStyle->SetOptStat(1);
  hl->Divide(xpos);
  hl->Scale(1/ hl->Integral(hl->GetXaxis()->FindBin(-99),
			    hl->GetXaxis()->FindBin(99)));
  hr->Divide(xpos);
  hr->Scale(1/ hr->Integral(hr->GetXaxis()->FindBin(-99),
			    hr->GetXaxis()->FindBin(99)));  
  
  hl->SetLineColor(2);
  hl->SetTitle("Left: Red | Right: blue");
  hl->GetXaxis()->SetRangeUser(-100,100);

  TCanvas *c1=new TCanvas("c1",fnm.c_str());
  c1->Divide(3);
  c1->cd(1);
  hl->DrawCopy();
  hr->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  c1->cd(2);
  hl->DrawCopy();
  TH1D *hrmir=mirror(hr);
  hrmir->DrawCopy("same");
  gPad->SetGridx(1);
  gPad->SetGridy(1);

  c1->cd(3);
  TH1D *ratio=(TH1D*)hl->Clone("ratio");
  ratio->SetTitle("Left/mirrorRight");
  ratio->Divide(hl,hrmir);
  ratio->DrawCopy();
  gPad->SetGridx(1);
  gPad->SetGridy(1);
    
  fin->Close();
}

TH1D *mirror(TH1D *a){
  TH1D *mir=(TH1D*)a->Clone(Form("%s_mirror",a->GetName()));
  int nbin=a->GetXaxis()->GetNbins();
  for(int i=1;i<=nbin;i++)
    mir->SetBinContent(i,a->GetBinContent(nbin+1-i));
  return mir;
}
  
