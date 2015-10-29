#include "RooStats/HistFactory/Measurement.h"

using namespace RooStats;

struct results{
	int status;
	double MU,B1,B2;
};

results minmzr(){
	results tmp;
	const char* filename = "results/StatModel_combined_meas_model.root";
	TFile* file = TFile::Open(filename);
	
	RooWorkspace* ws = (RooWorkspace*) file->Get("combined");
	// ws->Print();
	ModelConfig* mc = (ModelConfig*) ws->obj("ModelConfig");
	RooAbsData* data = (RooAbsData*) ws->data("obsData");
	RooSimultaneous* simPdf = (RooSimultaneous*) mc->GetPdf();
	RooAbsReal* nll = simPdf->createNLL(*data);
	
	RooMinimizer m(*nll);
	m.setStrategy(0);
	int status=m.minimize("Minuit2","Migrad");
	m.hesse();
	
	RooAbsReal* POI = (RooAbsReal*) mc->GetParametersOfInterest()->first();
	double Nui1 = ws->var("gamma_B1_bin_0")->getVal();
	double Nui2 = ws->var("gamma_B2_bin_0")->getVal();
	// cout<<"Status="<<status<<"\nMU="<<MU->getVal()<<endl;
	
	tmp.status=status;
	tmp.MU=POI->getVal();
	tmp.B1=Nui1;
	tmp.B2=Nui2;
	cout<<"MU="<<tmp.MU<<"\nB1="<<tmp.B1<<"\nB2="<<tmp.B2<<endl;
	return tmp;
}

void buildWS(){

   RooStats::HistFactory::Measurement meas("meas", "meas");

   meas.SetOutputFilePrefix("results/StatModel");
   meas.SetExportOnly(1);
   meas.SetPOI("MU");
   meas.SetLumi(1.0);
   meas.AddConstantParam("Lumi");

// Channel 1
   RooStats::HistFactory::Channel chan1("chan1");
   
   chan1.SetData("n1", "data/2bin.root"); 

   RooStats::HistFactory::Sample b("b", "b", "data/2bin.root");
   b.SetNormalizeByTheory(0);
   b.AddShapeFactor("B1");
   chan1.AddSample(b);         

   meas.AddChannel(chan1);
   
// Channel 2
   RooStats::HistFactory::Channel chan2("chan2");

   chan2.SetData("m1", "data/2bin.root");

   RooStats::HistFactory::Sample s1("s1", "s1", "data/2bin.root");
   s1.SetNormalizeByTheory(0);
   s1.AddNormFactor("MU", 1.0, -15.0, 15.0);
   chan2.AddSample(s1);     

   b.SetNormalizeByTheory(0);
   b.AddShapeFactor("B1");
   chan2.AddSample(b);    

   meas.AddChannel(chan2);
   
// Channel 3
   RooStats::HistFactory::Channel chan3("chan3");
   
   chan3.SetData("n2", "data/2bin.root");    

   RooStats::HistFactory::Sample c("c", "c", "data/2bin.root");
   c.SetNormalizeByTheory(0);
   c.AddShapeFactor("B2");
   chan3.AddSample(c);         

   meas.AddChannel(chan3);
   
// Channel 4
   RooStats::HistFactory::Channel chan4("chan4");

   chan4.SetData("m2", "data/2bin.root");

   RooStats::HistFactory::Sample s2("s2", "s2", "data/2bin.root");
   s2.SetNormalizeByTheory(0);
   s2.AddNormFactor("MU", 1.0, -15.0, 15.0);
   chan4.AddSample(s2);     

   c.SetNormalizeByTheory(0);
   c.AddShapeFactor("B2");
   chan4.AddSample(c);         

   meas.AddChannel(chan4);

// I/O
   meas.CollectHistograms();
   meas.PrintTree();
   RooStats::HistFactory::MakeModelAndMeasurementFast(meas);
}

void main_2bin(){
	const UInt_t N=1;
	UInt_t nbins=1;
	UInt_t i;
	UInt_t n1[N],m1[N],n2[N],m2[N];
	Double_t mu,b1,b2,s1,s2,c;
	// TRandom3 r;
	
	TFile f1("results/2bin.root","recreate");
	
	//Output histograms
	TH1F *hstatus = new TH1F("status","status",20,-2,2);
	TH1F *hMU = new TH1F("MU","MU",30,-15,15);
	TH1F *hMU_bli = new TH1F("MU_bli","MU_bli",30,-15,15);
	TH1F *hB1 = new TH1F("B1","B1",15,0,15);
	TH1F *hB1_bli = new TH1F("B1_bli","B1_bli",15,0,15);
	TH1F *hB2 = new TH1F("B2","B2",15,0,15);
	TH1F *hB2_bli = new TH1F("B2_bli","B2_bli",15,0,15);
	
	//Print input and results in .txt file
	// ofstream txtfile ("input_output.txt");

	for(i=0;i<N;i++)
	{
		TFile f("data/2bin.root","recreate");
		
		//Input histograms
		TH1F *hs1 = new TH1F("s1","s1",nbins,0,1);
		TH1F *hs2 = new TH1F("s2","s2",nbins,0,1);
		TH1F *hb = new TH1F("b","b",nbins,0,1);
		TH1F *hc = new TH1F("c","c",nbins,0,1);
		
		TH1F *hn1 = new TH1F("n1","n1",nbins,0,1);
		TH1F *hm1 = new TH1F("m1","m1",nbins,0,1);
		TH1F *hn2 = new TH1F("n2","n2",nbins,0,1);
		TH1F *hm2 = new TH1F("m2","m2",nbins,0,1);
		
		mu=1;
		b1=0.5;
		b2=0.5;
		
		b=1;
		c=1;
		s1=0.3;
		s2=0.4;
		
		hs1->SetBinContent(1,s1); hs2->SetBinContent(1,s2); hb->SetBinContent(1,b); hc->SetBinContent(1,c);
		
		// n1[i]=r.Poisson(b1);
		// m1[i]=r.Poisson(b1+mu*s1);
		// n2[i]=r.Poisson(b2);
		// m2[i]=r.Poisson(b2+mu*s2);
		
		n1[i]=1;
		m1[i]=0;
		n2[i]=2;
		m2[i]=3;

		hn1->SetBinContent(1,n1[i]); hm1->SetBinContent(1,m1[i]); hn2->SetBinContent(1,n2[i]); hm2->SetBinContent(1,m2[i]);
		
		f.Write();
		buildWS();
		results output=minmzr();
		
		hstatus->Fill(output.status);
		hMU->Fill(output.MU);
		hB1->Fill(output.B1);
		hB2->Fill(output.B2);
		
		if(output.status==0){
			hMU_bli->Fill(output.MU);
			hB1_bli->Fill(output.B1);
			hB2_bli->Fill(output.B2);
		}
		
		// Print in txtfile
		// if(i==0) {
		// 	txtfile << s1 << " " <<  s2 << " " << b1 << " " << b2 << " " << mu;
		// }
		// txtfile << "\n" << n1[i] << " " << m1[i] <<" " << n2[i] <<" " << m2[i];
		// txtfile << " " << output.MU << " " << output.B1 << " " << output.B2 << " " << output.status;
				
		f.Close();
	}
	
	f1.Write();
	f1.Close();
	txtfile.close();
	
	return;
}