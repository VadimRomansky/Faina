#ifndef CONSTANTS_H
#define CONSTANTS_H

const int randomParameter = 32*1024;
const double atomicUnionMass = 1.66053904E-24;
const double massProton = 1.67262177E-24;
const double massAlpha = 6.644656E-24;
const double massElectron = 0.910938291E-27;
const double massDeuterium = 3.34449696893E-24;
const double massHelium3 = 5.00823792874E-24;
const double massOxygen = 26.5601801672E-24;
const double massSilicon = 46.4567787264E-24;

const double kBoltzman = 1.3806488E-16;
const double speed_of_light = 2.99792458E10;
const double speed_of_light2 = speed_of_light * speed_of_light;
const double speed_of_light4 = speed_of_light2 * speed_of_light2;
const double me_c2 = massElectron * speed_of_light2;
const double electron_charge = 4.803529695E-10;
const double re = (electron_charge * electron_charge / (massElectron * speed_of_light2));
const double re2 = re*re;
const double hplank = 6.626E-27;
const double pi = 4*atan2(1.0,1.0);
const double four_pi = 4*pi;
const double criticalNuCoef = 3 * electron_charge / (4 * pi * massElectron * massElectron * massElectron * speed_of_light * speed_of_light4);
const double emissivityCoef = sqrt(3.0) * electron_charge * electron_charge * electron_charge / (massElectron * speed_of_light2);

//SN2009bb
//const double distance = 40*3.08*1.0E24;
//AT2018
//const double distance = 60*3.08*1.0E24;
//CSS161010
const double distance = 150*3.08*1.0E24;

//const std::string outputfileName = "radiation.dat";
//const std::string fileNameP = "../../tristan-mp-pitp/Pe";
//const std::string fileNameF = "../../tristan-mp-pitp/Fe";
//const std::string fileNameP = "Ee";
//const std::string fileNameF = "Fs";
//const std::string logFileName = "log.dat";


#endif