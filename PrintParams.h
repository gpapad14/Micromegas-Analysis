void PrintParams(TString filename, int Npars, double vRO, double vr, double vM, double vD, double gAmp, double gConv, double P, int Mesh, int Gas,  double ampliFactor, int MCA, int bins, double gain, double resol) {

cout<< "***********************************" << endl;
cout<<endl;

cout<< "Name of file:		= " << filename << endl;
cout<< "vReadOut (A):	vRO	= " << vRO <<endl;
cout<< "vRing (A):	vr	= " << vr <<endl;
cout<< "vMesh: 		vM	= " << vM <<endl;
cout<< "vDrift: 	vD	= " << vD <<endl;
cout<< "Amplific gap: 	gAmp	= " << gAmp <<endl;
cout<< "Drift/Conv gap: gConv	= " << gConv <<endl;
cout<< "Pressure: 	P	= " << P <<endl;
cout<< "Type of mesh: 	Mesh	= " << Mesh <<endl;
cout<< "Gas mixture: 	Gas	= " << Gas <<endl;
cout<< "Amplifier:	ampliFac= " << ampliFactor <<endl;
cout<< "MCA voltage: 	MCA	= " << MCA <<endl;
cout<< "Binning: 	bins	= " << bins <<endl;

cout<<endl;
cout<< "***********************************" << endl;

}
