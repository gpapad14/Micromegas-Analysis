void GetParamsFromFilename(TString word, int Npars, double& vRO, double& vr, double& vM, double& vD, double& gAmp, double& gConv, double& P, int& Mesh, int& Gas,  double& ampliFactor, int& MCA, int& bins) {
// Npars: the number of parameters

string nums[Npars];
string num="";
string letter="";

int i=0;
int ipar=0;
while(i<word.Length()) {
letter=word[i];

	if(letter.compare("_")) { // If {letter is not "_"}. It will continue if letter is NOT "_"
		num += letter;	
		//cout<<num<<" nope"<<endl;
	}
	else {
		nums[ipar]=num;
		ipar++;
		num="";
		//cout<<"F"<<endl;
	}
	i++;
}	

vRO	=atof(nums[0].c_str());
vr	=atof(nums[1].c_str());
vM	=atof(nums[2].c_str());
vD	=atof(nums[3].c_str());
gAmp	=atof(nums[4].c_str());
gConv	=atof(nums[5].c_str());
P	=atof(nums[6].c_str());
Mesh	=atoi(nums[7].c_str());
Gas	=atoi(nums[8].c_str());
ampliFactor  =atof(nums[9].c_str());
MCA	=atoi(nums[10].c_str());
bins	=atoi(nums[11].c_str());



return 0;
}
