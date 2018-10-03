#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TF1.h"
#include "TLegendEntry.h"
#include "TLegend.h"
#include "TFile.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TObject.h"

#include "ScanDirectory.h"
#include "Analysis.h"
#include "PrintParams.h"
#include "GetTempPres.h"

double bouteilleFunc(double A, double B, double voltage, double gap, double pressure, double temperature) {
	// voltage->Volt, gap->cm, Pressure->bar, Temperature->K
	double c0=1;
	double c1=1.;
	//c0=1.6;
	double x0 =  c0*A*pressure*gap/(temperature);
	double x1 = -B*pressure*gap/(temperature*c1*voltage);
	return exp( x0 * exp(x1)  );
}

void codeX() {

gROOT->SetBatch(kTRUE); // DO NOT OPEN CANVAS BUT STILL SAVE IMAGE. To change use "kFALSE".
gStyle->SetOptFit(0000);// Do not show the fit parameters.
// Scan \\directory// for files ending with the given \\extension// and return (1) the names in the pointer \\*filenamesList// and (2) the number of those files \\Nfiles//. Then save the names into a matrix of strings \\filenames[Nfiles]//. 

bool plotGain     = true;	// Process and produce the Gain curves
bool transparency = !plotGain; 	// Do the analysis either for the gain or the transparency.	
bool plotMaP	  = true;	// Plot the Moulin a Poivre datasets
bool plotRD3	  = true;	// Plot the RD3 datasets.

int Nfolders_MaP = 8; 	// Number of Moulin a Poivre datasets
int Nfolders_RD3 = 4;	// Number of RD3 detector datasets
int NfoldersTransp = 8;
int Nfolders;
int ifolder;
TString datasetSaveExtention = "";

if(plotMaP==true && plotRD3==false) { Nfolders = Nfolders_MaP; datasetSaveExtention="_MaP"; }
if(plotMaP==true && plotRD3==true)  { Nfolders = Nfolders_MaP + Nfolders_RD3; }
if(plotMaP==false && plotRD3==true) { Nfolders = Nfolders_RD3; datasetSaveExtention="_RD3"; }

if(transparency==true) { Nfolders=NfoldersTransp; }
const char *directory[Nfolders];
const char *dir_RD3[Nfolders_RD3];

if(plotGain==true && transparency==false) {
directory[0] = "/home/gpapad/Desktop/Stage/MaP/study/gap50um/";
directory[1] = "/home/gpapad/Desktop/Stage/MaP/study/gap65um/";
directory[2] = "/home/gpapad/Desktop/Stage/MaP/study/gap75um/";
directory[3] = "/home/gpapad/Desktop/Stage/MaP/study/gap90um/";
directory[4] = "/home/gpapad/Desktop/Stage/MaP/study/gap100um/";
directory[5] = "/home/gpapad/Desktop/Stage/MaP/study/gap125um/";
directory[6] = "/home/gpapad/Desktop/Stage/MaP/study/gap150um/";
directory[7] = "/home/gpapad/Desktop/Stage/MaP/study/gap190um/";

dir_RD3[0] = "/home/gpapad/Desktop/Stage/MaP/study/RD3_13.06.2018/";
dir_RD3[1] = "/home/gpapad/Desktop/Stage/MaP/study/RD3_15.06.2018/";
dir_RD3[2] = "/home/gpapad/Desktop/Stage/MaP/study/RD3_22.06.2018/";
dir_RD3[3] = "/home/gpapad/Desktop/Stage/MaP/study/RD3_28.06.2018/";

if(plotMaP==true && plotRD3==true)  {
	for(int i=Nfolders_MaP; i<Nfolders; i++) { directory[i]=dir_RD3[i-Nfolders_MaP]; }
}

}

if(plotGain==false && transparency==true) {
directory[0] = "/home/gpapad/Desktop/Stage/MaP/study/gap50um/transparency50um/";
directory[1] = "/home/gpapad/Desktop/Stage/MaP/study/gap65um/transparency65um/";
directory[2] = "/home/gpapad/Desktop/Stage/MaP/study/gap75um/transparency75um/";
directory[3] = "/home/gpapad/Desktop/Stage/MaP/study/gap90um/transparency90um/";
directory[4] = "/home/gpapad/Desktop/Stage/MaP/study/gap100um/transparency100um/";
directory[5] = "/home/gpapad/Desktop/Stage/MaP/study/gap125um/transparency125um/";
directory[6] = "/home/gpapad/Desktop/Stage/MaP/study/gap150um/transparency150um/";
directory[7] = "/home/gpapad/Desktop/Stage/MaP/study/gap190um/transparency190um/";
//directory[6] = "/home/gpapad/Desktop/Stage/MaP/study/RD3_22.06.2018/transpRD3/";
}

double maxGain[Nfolders-3];
double maxGainEr[Nfolders-3];

double alpha[Nfolders], alphaEr[Nfolders];
double beta[Nfolders], betaEr[Nfolders];

int marker[15];
marker[0] = 24;
marker[1] = 33;
marker[2] = 23;
marker[3] = 34;
marker[4] = 29;
marker[5] = 39;
marker[6] = 41;
marker[7] = 22;
marker[8] = 47;
marker[9] = 45;
marker[10] = 49;
marker[11] = 20;
marker[12] = 28;
marker[13] = 36;
marker[14] = 48;

int color[15];
color[0] = 1;
color[1] = 2;
color[2] = 3;
color[3] = 4;
color[4] = 5;
color[5] = 6;
color[6] = 7;
color[7] = 8;
color[8] = 800;
color[9] = 860;
color[10] = 881;
color[11] = 861;
color[12] = 898;
color[13] = 610;
color[14] = 801;


const char *extension = ".mca";

// Fit func: gain as a function of the voltage
TF1* bout_V_Func = new TF1("bout_V_Func", "exp( ([0])* exp(-[1]/x ) )", 0, 1000); // Expressed as a function of x=voltage.
bout_V_Func->SetLineWidth(1);
bout_V_Func->SetLineStyle(3);
TF1* expoFitFunc = new TF1("expoFitFunc", "exp( [0] + [1]*x  )", 0, 300); // Fit gain as a function of the el field
expoFitFunc->SetLineWidth(1);
expoFitFunc->SetLineStyle(3);
TLegend* legGainALL = new TLegend(0.75, 0.5, 0.99, 0.98); 	// legGainALL is used for both GainALL and ResolALL.
TLegend* legTranspALL = new TLegend(0.75, 0.5, 0.98, 0.98);

// Start the for loop over all the folders.
TCanvas* cGain_Eamp_ALL = new TCanvas("cGain_Eamp_ALL", "Plot all gain curves", 900, 600);
cGain_Eamp_ALL->SetLogy();
cGain_Eamp_ALL->SetGrid();
TCanvas* cGain_Vamp_ALL = new TCanvas("cGain_Vamp_ALL", "Plot all gain vs. v_amp", 900, 600);
cGain_Vamp_ALL->SetLogy();
cGain_Vamp_ALL->SetGrid();
TCanvas* cResol_Eamp_ALL = new TCanvas("cResol_Eamp_ALL", "Plot all resolution curves", 900, 600);
cResol_Eamp_ALL->SetGrid();
TCanvas* cResol_Eratio_ALL = new TCanvas("cResol_Eratio_ALL", "Plot all resolution curves", 900, 600);
cResol_Eratio_ALL->SetGrid();
TCanvas* cTranspALL = new TCanvas("cTranspALL", "Plot all transparency curves", 900, 600);
cTranspALL->SetGrid();
double gapAmp[Nfolders];

//*************************************************************************************
// Start processing all the folders
for( ifolder=0; ifolder<Nfolders; ifolder++) {

TString dir;
if(plotMaP==false && plotRD3==true) { dir = dir_RD3[ifolder]; }
else { dir = directory[ifolder]; }

int Nfiles; // Nfiles : number of ".mca" filesin the given directory
TString *filenamesList = ScanDirectory(dir, extension, Nfiles);

// Get temperature and pressure
double temperature, temperatureEr, pressure;
temperature=300.;
temperatureEr=0.2;
pressure=1.;
TString tempStr;
TString tempErStr;
GetTempPres(dir, temperature, temperatureEr, pressure);

// Definition of parameters
// filename[Nfiles] = "vRO_vr_vM_vD_gAmp_gConv_P_Mesh_Gas_ampliFactor_MCA_bins_.mca"
TString filenames[Nfiles]; // Name of the file as "string"
TString filename;	// Dummy variable holding a single filename.
int Npars=12; 		// Number of parameters

double vRO, vr, vM, vD; // Voltages [V]
double vEr = 0.1;	// Error of the voltage. [V]
double gAmp, gConv; 	// Gaps [μm,mm]
double gAmpEr = 3;	// Error of the amplification gap. [μm]
double P;		// Pressure [bar]
int Mesh, Gas;		// Mesh type, Gas mixture [0, 1...]
double ampliFactor;	// Amplification factor from Amplifier
int MCA;		// MCA maximum input voltage [5/10 V]
int bins;		// Number of bins [512, 1024...]

double Edrift[Nfiles];		// El. field in the conversion gap. [kV/cm]
double Eamp[Nfiles];		// El. field in the amplification gap. [kV/cm]
double EampEr;			// Error of the el. field in the amplification gap. [kV/cm]
double Vdrift;			// Voltage difference between Drift and Mesh |Vdrift-Vmesh|
double Vamp;			// Voltage difference between Anode and Mesh |Vanode-Vmesh|

// Definition of variable to extract.
double gain, gainEr;	// Gain
double resol, resolEr;	// Resolution
double peak, peakEr;	// Peak position (MCA channel)
double maxPeakTransp=0; 
double peakM[Nfiles], peakErM[Nfiles];	// Peak position (MCA channel)

// Save/write data in ".txt" file.
ofstream txt(dir+"MotherTree.txt");
txt << "vReadOut     vRing     vMesh      vDrift   gAmp gConv   Pres  Mesh  Gas ampliFactor  MCA  bins   gain   gainEr   resol    resolEr"<<std::endl;
// Definition of the Tree
TFile *tfile 		= new TFile(dir+"Mother.root","RECREATE");
TTree *mothertree	= new TTree("MotherTree", "Mother Tree");

TBranch *vRO_Br 	= mothertree->Branch("vReadOut", &vRO, "vRO/D");
TBranch *vr_Br 		= mothertree->Branch("vring", &vr, "vr/D");
TBranch *vM_Br 		= mothertree->Branch("vMesh", &vM, "vM/D");
TBranch *vD_Br 		= mothertree->Branch("vDrift", &vD, "vD/D");
TBranch *gAmp_Br	= mothertree->Branch("Amp_gap", &gAmp, "gAmp/D");
TBranch *gConv_Br 	= mothertree->Branch("Conv_gap", &gConv, "gConv/D");
TBranch *P_Br 		= mothertree->Branch("Pressure", &P, "P/D");
TBranch *Mesh_Br 	= mothertree->Branch("Mesh_type", &Mesh, "Mesh/I");
TBranch *Gas_Br		= mothertree->Branch("Gas_mix", &Gas, "Gas/I");
TBranch *ampliFactor_Br	= mothertree->Branch("ampliFactor", &ampliFactor, "ampliFactor/D");
TBranch *MCA_Br		= mothertree->Branch("MCA", &MCA, "MCA/I");
TBranch *bins_Br	= mothertree->Branch("bins", &bins, "bins/I");

TBranch *gain_Br 	= mothertree->Branch("Gain", &gain, "gain/D");
TBranch *gainEr_Br 	= mothertree->Branch("GainEr", &gainEr, "gainEr/D");
TBranch *resol_Br 	= mothertree->Branch("Resol_FWHM", &resol, "resol/D");
TBranch *resolEr_Br 	= mothertree->Branch("ResolEr_FWHM", &resolEr, "resolEr/D");

// Transparency curve
TGraphErrors* transparencyGain = new TGraphErrors(Nfiles);
TGraphErrors* transparencyResol = new TGraphErrors(Nfiles);
// Gain/Resolution curve
TGraphErrors* gain_Vamp  = new TGraphErrors(Nfiles);
gain_Vamp->SetName("gain_Vamp");
gain_Vamp->SetMarkerSize(1);
TGraphErrors* gain_Eamp  = new TGraphErrors(Nfiles);
gain_Eamp->SetName("gain_Eamp");
gain_Eamp->SetMarkerSize(1);
TGraphErrors* resol_Eamp = new TGraphErrors(Nfiles);
TGraphErrors* resol_Eratio = new TGraphErrors(Nfiles);

maxGain[ifolder]=0;
// Take the \\filenames// and extract the set parameters for each file. 
for(int i=0; i<Nfiles; i++) {

	filenames[i] = filenamesList[i];
	filename = filenames[i];
	//cout<< filename <<endl;

	Analysis(filename, dir, Npars, vRO, vr, vM, vD, gAmp, gConv, P, Mesh, Gas, ampliFactor, MCA, bins, gain, gainEr, resol, resolEr, peak, peakEr);
	Vdrift     = abs(vD-vM);
	Vamp       = abs(vM-vRO);
	Edrift[i]  = 1e-3*Vdrift/(1e-1*gConv);
	Eamp[i]	   = 1e-3*Vamp/(1e-4*gAmp);	
	peakM[i]   = peak; 
	peakErM[i] = peakEr;
	if(gain>maxGain[ifolder]) { maxGain[ifolder]=gain; maxGainEr[ifolder]=gainEr; }
	//PrintParams(filename, Npars, vRO, vr, vM, vD, gAmp, gConv, P, Mesh, Gas, ampliFactor, MCA, bins, gain, resol);
	
	mothertree->Fill();
	txt << vRO<<"    "<< vr<<"     "<<vM<<"     "<< vD<<"     "<< gAmp <<"     "<< gConv <<"     "<< P <<"     "<< Mesh <<"     "<< Gas <<"     "<< ampliFactor <<"     "<< MCA <<"     "<< bins <<"     "<< Form("%.1f",gain) <<"     "<< Form("%.1f",gainEr) <<"     "<< Form("%.2f",resol) <<"     "<< Form("%.2f",resolEr)<<std::endl;

	if(plotGain==true) {
// if you want to plot Gain vs. gap, use this line and change the range of x-axis, line ~285.
		EampEr = sqrt( pow(1e-3*vEr/(1e-4*gAmp) , 2) + pow( 1e-3*Vamp*1e-4*gAmpEr/( pow(1e-4*gAmp,2) ),2) );
		
		gain_Vamp->SetPoint(i, Vamp, gain); 
		gain_Vamp->SetPointError(i, vEr, gainEr);

		gain_Eamp->SetPoint(i, Eamp[i], gain); 
//		gain_Eamp->SetPointError(i, EampEr, gainEr);
		gain_Eamp->SetPointError(i, 0.0, gainEr);

		resol_Eamp->SetPoint(i, Eamp[i], resol);
//		resol_Eamp->SetPointError(i, EampEr, resolEr);
		resol_Eamp->SetPointError(i, 0.0, resolEr);

		resol_Eratio->SetPoint(i, Edrift[i]/Eamp[i], resol);
		resol_Eratio->SetPointError(i, 0.0, resolEr);
		
	}

	if(transparency==true) {
		if(maxPeakTransp<peak) {
			maxPeakTransp=peak;
		}
		transparencyResol->SetPoint(i, Edrift[i]/Eamp[i], resol);
		transparencyResol->SetPointError(i, 0.0, resolEr);
	}

}
gapAmp[ifolder]=gAmp;

if(transparency==true) {
for(int i=0; i<Nfiles; i++) {
	transparencyGain->SetPoint(i, Edrift[i]/Eamp[i], peakM[i]/maxPeakTransp);
	transparencyGain->SetPointError(i, 0.0, peakErM[i]/maxPeakTransp);
}
}

mothertree->Write();
txt.close();
//mothertree->Scan("vReadOut:vMesh:Amp_gap:Gain:GainEr:Resol_FWHM:ResolEr_FWHM");

// Plot gain/resolution.
if(plotGain==true) {
gStyle->SetOptFit(0000);
	gain_Eamp->SetMarkerStyle(20);
	gain_Eamp->SetMarkerColor(4);
	gain_Eamp->SetTitle("Gain vs. E_{amp}");
	gain_Eamp->GetXaxis()->SetTitle("E_{amp} (kV/cm)");
	gain_Eamp->GetYaxis()->SetTitle("Gain");
	resol_Eamp->SetMarkerStyle(20);
	resol_Eamp->SetMarkerColor(4);
	resol_Eamp->SetTitle("Resolution (FWHM) vs. E_{amp}");
	resol_Eamp->GetXaxis()->SetTitle("E_{amp} (kV/cm)");
	resol_Eamp->GetYaxis()->SetTitle("Resolution FWHM (%)");
	
	TCanvas* cGainResol = new TCanvas("cGain", "Gain/Resolution", 900, 600);	
	cGainResol->Divide(2,1);
	cGainResol->cd(1);
	cGainResol->cd(1)->SetGrid();
	cGainResol->cd(1)->SetLeftMargin(0.15);
	cGainResol->cd(1)->SetRightMargin(0.05);
	gain_Eamp->Draw("AP");
	gain_Eamp->Fit("expo"); 

	
	cGainResol->cd(2);
	cGainResol->cd(2)->SetGrid();
	resol_Eamp->Draw("AP");

	TString gapStr =Form("%.0f",gAmp);
	cGainResol->cd(0);
	TLegend* legGR = new TLegend(0.42, 0.96 , 0.58, 1);
	legGR->SetHeader("Amp gap = "+gapStr+"#mum", "C");
	legGR->Draw();

	cGainResol->SaveAs(dir+"Gain_Resolution_"+gapStr+"um.png","");


// Plot ALL gain/resolution curves in the same canvas.
	gain_Eamp->SetMarkerStyle(marker[ifolder]);
	gain_Eamp->SetMarkerColor(color[ifolder]);
	
	gain_Vamp->SetTitle("Gain vs. V_{amp}");
	gain_Vamp->GetXaxis()->SetTitle("V_{amp} (V)");
	gain_Vamp->GetYaxis()->SetTitle("Gain");
	gain_Vamp->SetMarkerStyle(marker[ifolder]);
	gain_Vamp->SetMarkerColor(color[ifolder]);

	resol_Eamp->SetMarkerStyle(marker[ifolder]);
	resol_Eamp->SetMarkerColor(color[ifolder]);

	resol_Eratio->SetMarkerStyle(marker[ifolder]);
	resol_Eratio->SetMarkerColor(color[ifolder]);

if(ifolder==0) {
	cGain_Eamp_ALL->cd(0);
	TAxis *axis0 = gain_Eamp->GetXaxis();
	axis0->SetLimits(20.,80.); 
	if(plotMaP==false && plotRD3==true) { axis0->SetLimits(26.,38.); }
	gain_Eamp->GetYaxis()->SetRangeUser(3e3,50e3);
	gain_Eamp->Draw("AP");
	expoFitFunc->SetLineColor(color[ifolder]);
	gain_Eamp->Fit("expoFitFunc");

	cGain_Vamp_ALL->cd(0);
	TAxis *axis0Vamp = gain_Vamp->GetXaxis();
	axis0Vamp->SetLimits(350.,650.); 
	gain_Vamp->GetYaxis()->SetRangeUser(3e3,50e3);
	gain_Vamp->Draw("AP");
	bout_V_Func->SetLineColor(color[ifolder]);
	gain_Vamp->Fit("bout_V_Func");

	cResol_Eamp_ALL->cd(0);
	TAxis *axis0resol = resol_Eamp->GetXaxis();
	axis0resol->SetLimits(20.,90.);
	if(plotMaP==false && plotRD3==true) { axis0resol->SetLimits(26.,38.); }
	resol_Eamp->GetYaxis()->SetRangeUser(20,65);
	resol_Eamp->Draw("AP");

	resol_Eratio->SetTitle("Resolution (FWHM) vs. E_{conv}/E_{amp}");
	resol_Eratio->GetXaxis()->SetTitle("E_{conv}/E_{amp}");
	resol_Eratio->GetYaxis()->SetTitle("Resolution FWHM (%)");
	cResol_Eratio_ALL->cd(0);
	TAxis *axis0resolRatio = resol_Eratio->GetXaxis();
	axis0resolRatio->SetLimits(0.,0.03);
	resol_Eratio->GetYaxis()->SetRangeUser(20,65);
	resol_Eratio->Draw("AP");


}
else {	
	cGain_Eamp_ALL->cd(0);
	gain_Eamp->Draw("PSAME");
	expoFitFunc->SetLineColor(color[ifolder]);
	gain_Eamp->Fit("expoFitFunc");

	cGain_Vamp_ALL->cd(0);
	gain_Vamp->Draw("PSAME");
	bout_V_Func->SetLineColor(color[ifolder]);
	gain_Vamp->Fit("bout_V_Func");

	cResol_Eamp_ALL->cd(0);	
	resol_Eamp->Draw("PSAME");

	cResol_Eratio_ALL->cd(0);
	resol_Eratio->Draw("PSAME");
}

	alpha[ifolder]= bout_V_Func->GetParameter(0) * temperature / ( pressure * gAmp*1e-4 );
	// Error includes fit-error and temperature-error
	alphaEr[ifolder]= sqrt( pow(bout_V_Func->GetParError(0) * temperature / ( pressure * gAmp*1e-4 ) ,2) + pow(bout_V_Func->GetParameter(0) * temperatureEr / ( pressure * gAmp*1e-4 ) ,2) );
	//beta[ifolder] = bout_V_Func->GetParameter(1) * 1e3 * temperature / pressure ;
	beta[ifolder] = bout_V_Func->GetParameter(1) * temperature / ( pressure * gAmp*1e-4 );
	betaEr[ifolder] = sqrt( pow(bout_V_Func->GetParError(1) * 1e3 * temperature / pressure, 2) + pow(bout_V_Func->GetParameter(1) * 1e3 * temperatureEr / pressure, 2) ) ;
	tempStr=Form("%.2f",temperature-273.15);
	tempErStr=Form("%.2f",temperatureEr);
	legGainALL->AddEntry(gain_Eamp, gapStr+" #mum @ "+tempStr+"#pm"+tempErStr+" #circC","P");

// Plot resolution ALL as a function of the electric field ratio.


// END  ---  GainALL

}
// End of "plotGain==true"
//##############################################################################


if(transparency==true) {
gStyle->SetOptFit(0000);	
	transparencyGain->SetMarkerStyle(20);
	transparencyGain->SetMarkerColor(4);
	transparencyGain->SetTitle("Transparency");
	transparencyGain->GetXaxis()->SetTitle("E_{drift}/E_{amp}");
	transparencyGain->GetYaxis()->SetTitle("Peak (MCA channel)");
	transparencyResol->SetMarkerStyle(20);
	transparencyResol->SetMarkerColor(4);
	transparencyResol->SetTitle("Resolution vs. E_{drift}/E_{amp}");
	transparencyResol->GetXaxis()->SetTitle("E_{drift}/E_{amp}");
	transparencyResol->GetYaxis()->SetTitle("Resolution (%)");

	TCanvas* cTransp = new TCanvas("cTransp", "Transparency", 900, 600);
	cTransp->Divide(2,1);
	cTransp->cd(1);
	cTransp->cd(1)->SetGrid();
	cTransp->cd(1)->SetLeftMargin(0.15);
	cTransp->cd(1)->SetRightMargin(0.05);
	transparencyGain->Draw("AP");
	
	cTransp->cd(2);
	cTransp->cd(2)->SetGrid();
	transparencyResol->Draw("AP");

	TString gapStr =Form("%.0f",gAmp);
	cTransp->cd(0);
	TLegend* legTransp = new TLegend(0.4, 0.94 , 0.6, 1);
	legTransp->SetHeader("Amp gap = "+gapStr+"#mum", "C");
	legTransp->Draw();

	cTransp->SaveAs(dir+"Transparency_"+gapStr+"um.png","");

	transparencyGain->SetMarkerStyle(marker[ifolder]);
	transparencyGain->SetMarkerColor(color[ifolder]);
	transparencyResol->SetMarkerStyle(marker[ifolder]);
	transparencyResol->SetMarkerColor(color[ifolder]);
	cTranspALL->cd(0);

if(ifolder==0) {
	TAxis *axis0 = transparencyGain->GetXaxis();
	axis0->SetLimits(0.,0.06);
	transparencyGain->GetYaxis()->SetRangeUser(0.6,1.05);
	transparencyGain->Draw("AP");
}
else{
	transparencyGain->Draw("PSAME");
}
	tempStr=Form("%.2f",temperature-273.15);
	tempErStr=Form("%.2f",temperatureEr);
	legTranspALL->AddEntry(transparencyGain, gapStr+" #mum @ "+tempStr+"#pm"+tempErStr+" #circC","P");

}
// End of "transparency==true"____



//cout<<endl;
//cout<<"Number of files found: "<< Nfiles << endl;

} // END  ---  ifolder-loop.

//##############################################################################


if(plotGain==true && transparency==false) {
cGain_Eamp_ALL->cd(0);
legGainALL->Draw();

cGain_Eamp_ALL->SaveAs("/home/gpapad/Desktop/Stage/MaP/study/Gain_vs_Eamp_ALL"+datasetSaveExtention+".png","");

cGain_Vamp_ALL->cd(0);
legGainALL->Draw();
cGain_Vamp_ALL->SaveAs("/home/gpapad/Desktop/Stage/MaP/study/Gain_vs_Vamp_ALL"+datasetSaveExtention+".png","");

cResol_Eamp_ALL->cd(0);
legGainALL->Draw();
cResol_Eamp_ALL->SaveAs("/home/gpapad/Desktop/Stage/MaP/study/Resolution_vs_Eamp_ALL"+datasetSaveExtention+".png","");

cResol_Eratio_ALL->cd(0);
legGainALL->Draw();
cResol_Eratio_ALL->SaveAs("/home/gpapad/Desktop/Stage/MaP/study/Resolution_vs_Eratio_ALL"+datasetSaveExtention+".png","");
}
if(plotGain==false && transparency==true) {
cTranspALL->cd(0);
legTranspALL->Draw();
cTranspALL->SaveAs("/home/gpapad/Desktop/Stage/MaP/study/TransparencyALL"+datasetSaveExtention+".png","");
}

//###################################################################################

gROOT->SetBatch(kFALSE);
TF1* function[Nfolders]; // x: voltage
TGraph* funcPoints = new TGraph(Nfolders-2);
TCanvas* cMAX = new TCanvas("cMAX", "sdkj", 700,500);
cMAX->SetLogx();
cMAX->SetGrid();
double voltage = 550;
TString voltStr=Form("%.0f",voltage);
TLegend* legcol = new TLegend(0.8,0.5,0.95,0.95);
legcol->SetHeader("V = "+voltStr,"C");

TCanvas* cAlpha = new TCanvas("cAlpha", "Alpha param", 700, 500);
cAlpha->SetGrid();
TGraphErrors* alphaGraph = new TGraphErrors(Nfolders);
alphaGraph->GetXaxis()->SetTitle("gap (#mum)");
alphaGraph->GetYaxis()->SetTitle("A");
alphaGraph->SetMarkerStyle(20);
alphaGraph->SetMarkerColor(4);
TCanvas* cBeta = new TCanvas("cBeta", "Beta param", 700, 500);
cBeta->SetGrid();
TGraphErrors* betaGraph = new TGraphErrors(Nfolders);
betaGraph->GetXaxis()->SetTitle("gap (#mum)");
betaGraph->GetYaxis()->SetTitle("B");
betaGraph->SetMarkerStyle(20);
betaGraph->SetMarkerColor(4);

for(int i=0; i<Nfolders; i++) {
cout<< "@ gap = "<<gapAmp[i]<<": A = "<<alpha[i]<<" +- "<<alphaEr[i]<<", B = "<<beta[i]<<" +- "<<betaEr[i]<<endl;
//cout<< "gap = "<<gapAmp[i]<<": A = "<<alpha[i]<<", B = "<<beta[i]<<endl;
alphaGraph->SetPoint(i, gapAmp[i], alpha[i]);
alphaGraph->SetPointError(i, 0.0, alphaEr[i]);
betaGraph->SetPoint(i, gapAmp[i], beta[i]);
betaGraph->SetPointError(i, 0.0, betaEr[i]);

function[i] = new TF1("function", "exp( ([0]*x)* exp(-[1]*x ) )", 9.99, 200);
function[i]->SetParameters(alpha[i]*1*1e-4/300, beta[i]*1*1e-4/(300*voltage));
//cout<< bouteilleFunc(alpha[i], beta[i],  voltage, gapAmp[i], 1., 300) <<endl;	

	function[i]->SetLineColor(i+1);
if(i==0) { 
	function[i]->Draw(""); 
	function[i]->SetTitle("Gain vs. gap");
	function[i]->GetXaxis()->SetTitle("gap (um)");
	function[i]->GetYaxis()->SetTitle("Gain");
	function[i]->GetYaxis()->SetRangeUser(0,1e8);

	funcPoints->SetPoint(i, gapAmp[i], bouteilleFunc(alpha[i], beta[i],  voltage, 1e-4*gapAmp[i], 1., 300));
}
else { 
	function[i]->Draw("SAME"); 
	if(i<Nfolders-3){ funcPoints->SetPoint(i, gapAmp[i], bouteilleFunc(alpha[i], beta[i],  voltage, 1e-4*gapAmp[i], 1., 300)); }
}
TString distance=Form("%.0f",gapAmp[i]);
legcol->AddEntry(function[i], distance+" um","l");

}
legcol->Draw();
//cMAX->SaveAs("/home/gpapad/Desktop/Stage/MaP/study/CurvesALL_V"+ voltStr+".png","");

cAlpha->cd(0);
alphaGraph->Draw("AP");
cBeta->cd(0);
betaGraph->Draw("AP");


gStyle->SetOptFit(0001);
TCanvas* cMAXGraph = new TCanvas("cMAXGraph","Graph",700,500);
cMAXGraph->SetGrid();
funcPoints->SetMarkerStyle(20);
funcPoints->SetMarkerColor(2);
funcPoints->GetYaxis()->SetTitle("Gain");
funcPoints->GetXaxis()->SetTitle("gap (um)");
funcPoints->Draw("AP");
funcPoints->Fit("function");
//cMAXGraph->SaveAs("/home/gpapad/Desktop/Stage/MaP/study/GraphALL_V"+ voltStr+".png","");




/*
gROOT->SetBatch(kTRUE);
TCanvas* cGap = new TCanvas("cGap", "V=500V",700,500);
cGap->SetGrid();
cGap->SetLogy();
cGap->SetLogx();
int Npoints=1;
int ipoint=0;
int ivar=0;
double vMin=400;
TF1* myFitFunc = new TF1("myFitFunc","exp([0]*x *exp([1]*x) )",  0. , 200);
myFitFunc->SetParameters(2e6*1/300, 3e7*1/(300*vMin));
for(double volt=vMin; volt<vMin+Npoints*50; volt=volt+50) {
TGraph* gapCurve = new TGraph(Nfolders-2);
gapCurve->GetYaxis()->SetRangeUser(1000,1e6);
	cout<<volt<<endl;
for(int i=0; i<Nfolders-2; i++) {
	gapCurve->SetPoint(i, gapAmp[i], bouteilleFunc(alpha[i], beta[i],  volt, 1e-4*gapAmp[i], 1., 25.5+273.15) );
	cout<< bouteilleFunc(alpha[i], beta[i],  volt, 1e-4*gapAmp[i], 1., 25.5+273.15) <<endl;	
}
gapCurve->SetMarkerStyle(20);
gapCurve->SetMarkerColor(ipoint+2);
if(ipoint==0) { gapCurve->Draw("APL");  }
else {gapCurve->Draw("PLSAME"); }
ipoint++;

}
cGap->SaveAs("/home/gpapad/Desktop/Stage/MaP/study/gap.png","");

*/
/*
TCanvas* cMaxGain = new TCanvas("cMaxGain", "MaxGain",700,500);
cMaxGain->SetGrid();
//cGap->SetLogy();
TGraphErrors* maxG = new TGraphErrors(Nfolders-1);
double maximum=0;
for(int jvar=0; jvar<Nfolders-1; jvar++) {
	if(maximum<maxGain[jvar]) {
		maximum=maxGain[jvar];
	}
}
for(int jvar=0; jvar<Nfolders-1; jvar++) {
	maxG->SetPoint(jvar, gapAmp[jvar], maxGain[jvar]/maximum);
	maxG->SetPointError(jvar, 0.0, maxGainEr[jvar]/maximum);
}
maxG->SetTitle("Maximum gain vs. gap");
maxG->GetYaxis()->SetTitle("Maximum gain (normalized)");
maxG->GetXaxis()->SetTitle("Amplification gap (#mum)");
maxG->SetMarkerStyle(20);
maxG->SetMarkerColor(2);
maxG->Draw("AP");

cMaxGain->SaveAs("/home/gpapad/Desktop/Stage/MaP/study/MaxGain_gap.png","");
*/



}



