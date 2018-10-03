void SmouthHisto(TString filename, TString dir, TH1F *histo, int radius, int threshold, double& peak, double& sigma, double& norm) {

cout<<filename<<endl;
gROOT->SetBatch(kTRUE);
int bins;
bins=histo->GetSize()-2;
double channel[bins];
TH1F *h = new TH1F("data","smouth+fit",bins,0.,bins);

double max=0;
int maxPos=0;

int range=2*radius+1;
double val;
int j;

for(j=1;j<=bins;j++) {
	
	if(j>=radius && j<bins-radius) {
		val=0;
		for(int k=j-radius; k<=j+radius;k++){
			val+=histo->GetBinContent(k);

		}
		
		h->SetBinContent(j,val/range); 
			
		if(j>threshold && val/range>max) {
			maxPos=j;
			max=val/range;

		}
	}
	else {
		h->SetBinContent(j,histo->GetBinContent(j));
	}
//cout<<histo->GetBinContent(k)<<endl;	

//cout<<h->GetBinContent(j)<<endl;	
}
TCanvas* cSN = new TCanvas(filename,filename,700,500);
cSN->SetGrid();

//cout<<endl;
//cout<<maxPos<<endl;
//cout<<endl;

//h->Draw();
h->Fit("gaus","","fitFuncSmooth",maxPos*0.9,maxPos*1.2);
TF1* fitfunc = h->GetFunction("gaus");

norm	=fitfunc->GetParameter(0);
peak	=fitfunc->GetParameter(1);
sigma	=fitfunc->GetParameter(2);

/*
TH1F *hfit = new TH1F("dataFit","smouth+normalized",bins,0.,bins);
int i;
double x;
for(i=1;i<=bins;i++) {
	x=h->GetBinContent(i);
	hfit->SetBinContent(i,x/max);
}
hfit->Draw("");
hfit->Fit("fitfunc","","Smooth+Normalized",maxPos*0.9,maxPos*1.2);
*/

cSN->SaveAs(dir+"smoothData/"+filename+".png","");


//return hfit;


return 0;

}

