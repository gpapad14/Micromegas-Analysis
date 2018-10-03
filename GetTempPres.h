void GetTempPres(TString directory, double& temperature, double& temperatureEr, double& pressure) {

// Moulin a Poivre conditions
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/gap50um/") {
	temperature   = 26.22 + 273.15;
	temperatureEr = 0.13;
	pressure      = 0.9998;
}
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/gap65um/") {
	temperature   = 25.41 + 273.15;
	temperatureEr = 0.11;
	pressure      = 1.0003;
}
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/gap75um/") {
	temperature   = 25.69 + 273.15;
	temperatureEr = 0.02;
	pressure      = 1.0017;
}
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/gap90um/") {
	temperature   = 25.48 + 273.15;
	temperatureEr = 0.04;
	pressure      = 1.0017;
}
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/gap100um/") {
	temperature   = 24.93 + 273.15;
	temperatureEr = 0.05;
	pressure      = 1.0009;
}
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/gap125um/") {
	temperature   = 26.87 + 273.15;
	temperatureEr = 0.06;
	pressure      = 0.9994;
}
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/gap150um/") {
	temperature   = 30.82 + 273.15;
	temperatureEr = 0.04;
	pressure      = 1.0015;
}
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/gap190um/") {
	temperature   = 31.04 + 273.15;
	temperatureEr = 0.03;
	pressure      = 1.0006;
}


// RD3 condtions
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/RD3_13.06.2018/") {
	temperature   = 30.4 + 273.15;
	temperatureEr = 0.05;
	pressure      = 1.000;
}
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/RD3_15.06.2018/") {
	temperature   = 24.7 + 273.15;
	temperatureEr = 0.05;
	pressure      = 1.000;
}
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/RD3_22.06.2018/") {
	temperature   = 28.4 + 273.15;
	temperatureEr = 0.03;
	pressure      = 1.0109;
}
if (directory=="/home/gpapad/Desktop/Stage/MaP/study/RD3_28.06.2018/") {
	temperature   = 30.9 + 273.15;
	temperatureEr = 0.02;
	pressure      = 1.0022;
}





}
