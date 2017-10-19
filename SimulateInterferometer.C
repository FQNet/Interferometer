#include "TMath.h"

double GetModulationShape(double time, double RFModulationPeriod, double RFModulationVoltage, double VPI ) {

  //parameters
  double RFModulationFraction = RFModulationVoltage / VPI;

  double result = 0;  
  double cycleNumber = int( time / RFModulationPeriod);
  double timePhase = time - cycleNumber*RFModulationPeriod;
  
  double OnOrOff = 0;
  if (time > -1*(RFModulationPeriod/2) && time < RFModulationPeriod/2) OnOrOff = 0;
  else OnOrOff = 1;

  result = RFModulationFraction * OnOrOff;
  return result;
}


void SimulateInterferometer() {

  //RF modulation parameters
  double RFModulationVoltage = 2.4; //in V
  double VPI = 2.4; //in V
  double RFModulationPeriod = 30; //in ns


  //input parameters
  double LightSpeed = 3e8; // in nm / ns 
  double wavelength = 1550; //in nm
  double w = 2*TMath::Pi()*LightSpeed/wavelength; // frequency of light
  double k = 2*TMath::Pi() / wavelength;
  double LightPeriod = wavelength / LightSpeed; // in ns
  double amp1 = 10; //arbitrary units
  double amp2 = 10; //arbitrary units
  double L = 1e9; //in nm
  double deltaL = 1e7; //in nm

  //photodetector parameters
  double PhotodetectorCharPeriod = 0.35 / 38; // in ns; it's about 10ps.


  double xmin = -1; //in ns
  double xmax = 1; //in ns
  double binwidth = 0.01; //in ns
  int nbins = (xmax - xmin) / binwidth;

  // TH1D *output = new TH1D( "pulse", ";time;amplitude", nbins, xmin, xmax);
  vector<double> xVector;
  vector<double> yVector;

  //time bins on scope
  for (int i = 0; i < nbins; i++) {
    double t = xmin + i*binwidth;

    double value = 0;
    //need to do an integral over period characteristic to the time scale of the photodetector.
    double samplesPerLightPeriod = 10; //should sample at least 10 
    double nIntegralBins = samplesPerLightPeriod * (PhotodetectorCharPeriod / LightPeriod ) ;
    for (int k=0 ; k < nIntegralBins ; k++) {	       
	double dummyT = t - 0.5*PhotodetectorCharPeriod + k*(LightPeriod / samplesPerLightPeriod);
	double Pulse = GetModulationShape(dummyT,RFModulationPeriod,RFModulationVoltage, VPI);

	double Efield = amp1 * sin( k*2*(L+deltaL) - w*dummyT + TMath::Pi()*0.5*Pulse) 
	  + amp2 * sin( k*2*L - w*dummyT + TMath::Pi()*0.5*Pulse);
	double Power = 0.5*Efield*Efield;

	//add to the intgral
	value += Power * LightPeriod; //power * dt
      }
    value = value / PhotodetectorCharPeriod; //divide by full period to get the average power.

    cout << "time: " << t << " : average Power = " << value << "\n";
    xVector.push_back(t);
    yVector.push_back(value);
  }


  //Make Graph
  double xCoordinates[nbins];
  double yCoordinates[nbins];
  double maxY = 0;
  for (int i=0; i<nbins; i++) {
    if (yCoordinates[i] > maxY) maxY = yVector[i];
  }

  for (int i=0; i<nbins; i++) {
    xCoordinates[i] = xVector[i];
    yCoordinates[i] = yVector[i] / maxY;
  }
  TGraph *graph = new TGraph( nbins, xCoordinates, yCoordinates);
  graph->GetYaxis()->SetRangeUser(0, 1.2);
  graph->SetMarkerStyle(20);
  TCanvas *cv = new TCanvas("cv","cv", 800,800);
  graph->Draw("AP");
  

}
