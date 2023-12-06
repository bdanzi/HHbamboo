
void plotPU(){
  TFile* fMCUL18 = TFile::Open("pileup_UL_2018.root");  
  //  TFile* fDataUL18 = TFile::Open("puWeight_Spring2016MC_to_2017Data_294927-301141.root");
  //  TFile* fDataUL18 = TFile::Open("puWeight_Spring2016MC_to_2017Data_294927-301997.root");
  //TFile* fDataUL18 = TFile::Open("DataPileupHistogram2018UL_69200_80bins.root");
  //TFile* fDataUL18 = TFile::Open("DataPileupHistogram2018UL_69200_80bins.root");
  TFile* fDataUL18 = TFile::Open("DataPileupHistogram2016UL_69200_80bins.root");   
  // TFile* fDataUL18 = TFile::Open("DataPileupHistogram2017UL_69200_80bins.root");   

  auto hMC_out_of_the_box = (TH1F*) fMCUL18->Get("MC_out_of_the_box");
  auto hMC_reweighted = (TH1F*) fMCUL18->Get("MC_reweighted");
  auto hDataUL18_original = (TH1F*) fDataUL18->Get("pileup");
  auto hDataUL18 = new TH1F("Data","",80,0,80);
  bool isDataPoint = true;
  auto canvas = new TCanvas("canvas","canvas",600,400,600,400);
  gStyle->SetOptTitle(false);
  canvas->cd();
  for(int i=1;i<81; i++){
    hDataUL18->SetBinContent(i,hDataUL18_original->GetBinContent(i));
  }
  hDataUL18->Scale(1./hDataUL18->Integral());
  hDataUL18->SetStats(false);
  hMC_out_of_the_box->SetMaximum(0.05);
  hMC_out_of_the_box->SetStats(false);
  hMC_out_of_the_box->SetLineColor(30);
  hMC_out_of_the_box->SetFillColor(30);
  hMC_out_of_the_box->GetXaxis()->SetTitle("True number of interactions");
  hMC_out_of_the_box->GetYaxis()->SetTitle("Normalized to unity");
  hMC_reweighted->SetLineColor(kBlue-10);
  hMC_reweighted->SetFillColor(kBlue-10);
  if(isDataPoint){
  hDataUL18->SetMarkerStyle(20);
  hDataUL18->SetMarkerSize(0.8);
  hDataUL18->SetMarkerColor(kBlack);
  
  }
  
  else{
  hDataUL18->SetLineColor(kRed-2);
  hDataUL18->SetFillColor(kRed-2);
  hDataUL18->SetFillStyle(3004);
  
  }

  if(isDataPoint){
    hMC_out_of_the_box->Draw();
    hMC_reweighted->Draw("samehisto");
    hDataUL18->Draw("sameP");
  }
  else{
    hMC_out_of_the_box->Draw();
    hDataUL18->Draw("samehisto");
    hMC_reweighted->Draw("samehisto");
  }
  auto leg = new TLegend(0.56,0.71,0.88,0.88);
  leg->AddEntry(hMC_reweighted, "MC reweighted", "f");
  leg->AddEntry(hMC_out_of_the_box, "MC out of the box", "f");
  //leg->AddEntry(hData2016,"Data 2016", "p");
  if(isDataPoint){
    //leg->AddEntry(hDataUL18,"Data 2018, 59.7 fb^{-1}","p");}
    leg->AddEntry(hDataUL18,"Data 2017, 59.7 fb^{-1}","p");}
  // leg->AddEntry(hDataUL18,"Data 2016, 59.7 fb^{-1}","p");}
  // else{leg->AddEntry(hDataUL18,"Data 2018, 59.7 fb^{-1}","f");}
  // else{leg->AddEntry(hDataUL18,"Data 2017, 59.7 fb^{-1}","f");}
  else{leg->AddEntry(hDataUL18,"Data 2016, 59.7 fb^{-1}","f");}
  leg->SetLineColor(kWhite);
  leg->Draw();
  canvas->SaveAs("Pileup2018.png");
  canvas->SaveAs("Pileup2018.pdf");
}
