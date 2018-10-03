#include "GetParamsFromFilename.h"
#include "SmoothHisto.h"
#include "GetElChainCalib.h"

void Analysis(TString filename, TString dir, int Npars, double& vRO, double& vr, double& vM, double& vD, double& gAmp, double& gConv, double& P, int& Mesh, int& Gas,  double& ampliFactor, int& MCA, int& bins, double& gain, double& gainEr, double& resol, double& resolEr, double& peakPos, double peakPosEr) {

GetParamsFromFilename(filename, Npars, vRO, vr, vM, vD, gAmp, gConv, P, Mesh, Gas, ampliFactor, MCA, bins);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
int radius;
int threshold;
if 	(bins==1024) {	threshold=100;	radius=5; if(ampliFactor==25){ threshold+=30; } }
else if (bins==512)  {	threshold=50;	radius=2; if(ampliFactor==25){ threshold+=15; }}
//_______________________________________
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// Drow histogram
TH1F *histo = new TH1F("Original", filename, bins, 0., bins);
histo->GetXaxis()->SetTitle("Channel");
histo->GetYaxis()->SetTitle("Counts");
ifstream inp;

// "data/" is the name of the folder where there are the data files (.mca).
inp.open(dir+filename);
int x;
int i = 1;

string line;
bool found=false;
bool comments = false;
while(inp>>line){
	if( !found || comments) {
		int j=0;
		string word="";
		while(j<line.length()) {
			word+=line[j];
			j++;	
		}
		if(!word.compare("<<PMCA")){ comments=true;}
		if(!word.compare("<<DATA>>")){ found=true; continue; }
		
	}
	if((!comments || found) && (i<=bins)){
		x=atoi(line.c_str());
		histo->SetBinContent(i,x); 
//cout<<histo->GetBinContent(i)<<endl;
		i++;
	}
}	
inp.close();

gStyle->SetOptStat("ne");
//histo->Draw();



double peak, sigma, norm;
SmouthHisto(filename, dir, histo, radius, threshold, peak, sigma, norm);

gStyle->SetOptFit(0001);

// Normalized histogram
TH1F *histoNF = new TH1F("Normalized", filename, bins, 0., bins);
histoNF->GetXaxis()->SetTitle("Channel");
histoNF->GetYaxis()->SetTitle("Counts");
for(i=1;i<=bins;i++) {
	histoNF->SetBinContent(i,histo->GetBinContent(i)/norm);
}


TCanvas* c = new TCanvas(filename, filename, 700, 500);
c->SetGrid();

histo->Draw("");

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  I AM HERE ! ! !
//TF1* doubleGaus = new TF1("double Gaussian", "[0]*exp(-0.5*(( x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2) ", 0, bins);
TF1* doubleGaus = new TF1("double Gaussian", "[0]*exp(-0.5*(( x-[1])/[2])**2) ", 0, bins);

doubleGaus->SetParameter(0, 1);
doubleGaus->SetParameter(1, peak);
doubleGaus->SetParameter(2, sigma);
doubleGaus->SetParNames ("Amplitude", "Mean", "Sigma");
/*
doubleGaus->SetParLimits(3,0.,1.);
doubleGaus->SetParameter(4, (6.4/5.9)*peak);
doubleGaus->SetParLimits(5,0.1,100.);
*/

histo->Fit("double Gaussian","","fitFunc",peak-sigma*0.7,peak+sigma*2);
gPad->Update();
TPaveStats *params = (TPaveStats*)histo->FindObject("stats");
	params->SetX1NDC(0.69);
	params->SetX2NDC(0.99);
	
	params->SetY1NDC(0.6);
	params->SetY2NDC(0.9);
	

TF1* fDGaus = histo->GetFunction("double Gaussian");
fDGaus->SetLineColor(2);
//fDGaus->Draw("SAME");

TF1* gaus59 = new TF1("5.9keV_peak", "gaus(0)", 0, bins);
double amplitude59	= fDGaus->GetParameter(0);
double peak59	= fDGaus->GetParameter(1);
double peak59er = fDGaus->GetParError(1);
peakPos   = peak59;
peakPosEr = peak59er;
double sigma59	= fDGaus->GetParameter(2);
double sigma59er= fDGaus->GetParError(2);
gaus59->SetParameters(amplitude59, peak59, sigma59);
gaus59->SetLineColor(kGreen+2);
gaus59->SetLineStyle(2);
gaus59->Draw("SAME");
/*
TF1* peak64 = new TF1("5.9keV peak", "gaus(0)", 0, bins);
peak64->SetParameters(fDGaus->GetParameter(3), fDGaus->GetParameter(4), fDGaus->GetParameter(5));
peak64->SetLineColor(4);
peak64->Draw("SAME");
*/



// Calculate the GAIN ==============================

// The signal will finally arrive to the MCA and will sorted in bins
// peak59 (channel) -> vamp
//   1024 (channel) -> 10V 
double vamp =MCA*peak59/bins; 	// in Volts.
// The signal will be amplified by the amplifier by a factor of "ampliFactor".
double vpreamp = 1000*vamp/ampliFactor; // in mV.
// The charge will be translated into voltage in mV from the preamplifier.
// Creation of 1 electron-hole in Si at 300K -> needs 3.62eV.
// The "nel" electrons -> x
// 1MeV Si will give 45mV. I have "vpreamp".
// The final number of electrons on the anode will be "nel"
double nel = vpreamp* 1000000/(3.62*45);
// The factor "1000000/(3.62*20)" is coming from the manual of the ORTEC Pre-Amplifier 142B, which is the model we use. 

// The X-ray of 5.9keV will produce "prim" (primary) electrons.
int nprim=226;

//gain = nel/nprim;
//gain   = (1000.*(MCA*peak59/bins)/ampliFactor) * 1000000/(3.62*20)/nprim;
//gainEr = (1000.*(MCA*peak59er/bins)/ampliFactor) * 1000000/(3.62*20)/nprim;

double alpha, alphaErr;
double beta, betaErr;
alpha= 3.1e4; alphaErr= 0.0; beta= 1e5; betaErr= 0.0;
GetElChainCalib(dir, alpha, alphaErr, beta, betaErr);


gain = (alpha*peak59 + beta)/nprim;
gainEr = sqrt(  ( (alpha*peak59er)/nprim )*( (alpha*peak59er)/nprim ) + ( (alphaErr*peak59)/nprim )*( (alphaErr*peak59)/nprim ) + ( betaErr/nprim )*( betaErr/nprim )  );

TString gainString="";
gainString+=Form("%.1f", gain);
//gainString+=Form("%.2f", gain);
TString gainStringEr="";
gainStringEr += Form("%.1f", gainEr);

// Calculate the RESOLUTION ========================

double FWHM59 = sqrt(8*log(2))* sigma59;
double FWHM59er = sqrt(8*log(2))* sigma59er;
resol = FWHM59/peak59*100; // Resolution (FWHM) in %.
resolEr = 100* pow( pow( FWHM59er/peak59 , 2) + pow( FWHM59*peak59er/pow(peak59,2) , 2) , 0.5);
TString resolString="";
resolString += Form("%.2f", resol);
TString resolStringEr="";
resolStringEr += Form("%.2f", resolEr);


TLegend *leg;
leg = new TLegend(0.69, 0.54, 0.99, 0.6);
leg->SetHeader(" Gain = " + gainString + "#pm " + gainStringEr);
TLegendEntry *headerGain = (TLegendEntry*)leg->GetListOfPrimitives()->First();
headerGain->SetTextSize(.04);
headerGain->SetTextFont(12);
headerGain->SetTextColor(kOrange+7);
leg->Draw("SAME");

leg = new TLegend(0.69, 0.48, 0.99, 0.54);
leg->SetHeader(" Resol(FWHM) = "+resolString+"%");
TLegendEntry *headerResol = (TLegendEntry*)leg->GetListOfPrimitives()->First();
headerResol->SetTextSize(.04);
headerResol->SetTextFont(12);
headerResol->SetTextColor(kAzure);
leg->Draw("SAME");


// "data/" is the name of the folder where the images are going to be saved.
c->SaveAs(dir+filename + ".png","");


return 0;
}

