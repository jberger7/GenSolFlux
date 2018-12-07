//____________________________________________________________________________
/*!

\program gsolflux

\brief   A tool to generate a Tree of GSimpleNtpFlux entries for monochromatic
         dark matter incident on a detector from the Sun.

	 Syntax:
	   gsolflux [-h]
	             -n nev
		     -e energy
		     -m mass
		    [-o outfile_name]
		    [-d date_range]
	 Options:
	   [] Denotes an optional argument.
	   -h
	      Prints out help on using gsolflux and exists.

	   -n
	      Specifies number of incident dark matter to generate.
	     
	   -e
	      Specifies the dark matter energy in GeV.

           -m
	      Specifies the dark matter mass in GeV.

           -o
              Specifies the name of the output ROOT file.  Default: flux.root

           -d
	      Specifies a range of dates over which to generate incident dark
	      matter, input as MMDDYYYY,MMDDYYY where the first item is the
	      starting date and the last is ending date.  Default: The year 
	      2018.

\author  Joshua Berger <josh.berger \at pitt.edu>
         University of Pittsburgh
         Robert Hatcher <rhatcher \at fnal.gov>
	 Fermilab

\created December 5, 2018

*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <string>
#include <ctime>

#include <TFile.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TDatime.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/GVersion.h"
#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Tools/Flux/GSimpleNtpFlux.h"

extern "C" {
#include "SolTrack.h"
}

using std::string;

using namespace genie;
using namespace genie::flux;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

//User-specified options:
int         gOptNevents = 10000;           // number of incident Dm to generate
double      gOptDMEnergy;                  // dark matter energy E
double      gOptDMMass;                    // dark matter mass M
int         gOptDateStart = 20180101;      // starting date
int         gOptDateEnd = 20181231;        // ending date
double      gOptDetLat = 44.35;            // detector latitude (degrees) -- default SURF
double      gOptDetLong = -103.75;         // detector longitude (degrees) -- default SURF
double      gOptBeamAngle = 277.64;        // beam angle with respect to North (degrees) -- default LBNF to SURF
string      gOptOutFileName = "flux.root"; // output file name

// main routine
int main(int argc, char *argv[])
{ 
  GetCommandLineArgs(argc,argv);

  // Store the DM kinematics conveniently
  double Ed = gOptDMEnergy;
  double Md = gOptDMMass;
  double pd = TMath::Sqrt(Ed*Ed - Md*Md);

  // Create the ROOT file
  TFile* file = TFile::Open(gOptOutFileName.c_str(),"RECREATE");
  TTree* fluxntp = new TTree("flux","a simple flux n-tuple");
  TTree* metantp = new TTree("meta","metadata for flux n-tuple");
  genie::flux::GSimpleNtpEntry* fentry = new genie::flux::GSimpleNtpEntry;
  genie::flux::GSimpleNtpMeta*  fmeta  = new genie::flux::GSimpleNtpMeta;
  fluxntp->Branch("entry",&fentry);
  metantp->Branch("meta",&fmeta);

  // Set the metadata
  UInt_t metakey = TString::Hash(gOptOutFileName.c_str(),strlen(gOptOutFileName.c_str()));
  // ensure that we don't get smashed by UInt_t vs Int_t
  metakey &= 0x7FFFFFFF;

  // Some soltrack options
  int useDegrees = 0, useNorthEqualsZero = 1, computerRefrEquatorial = 0, computeDistance = 0;

  // Set our detector location
  struct Location loc;

  // Convert detector position to radians
  loc.latitude  = gOptDetLat * PI / 180.;
  loc.longitude =  gOptDetLong * PI / 180.;  // MicroBooNE

  // SolTrack requires some atmospheric conditions--pick something reasonable, but shouldn't be an effect
  loc.pressure = 101.0;
  loc.temperature = 300.0;

  // Store the beam angle in radians, converted to the right coord system
  const double angle = (360. - gOptBeamAngle) * PI / 180.;

  // Determine the range of times to generate.  We go through the ROOT time format
  // because C time messes with the times in odd ways
  TDatime startTime(gOptDateStart, 0);
  TDatime endTime(gOptDateEnd, 235959);
  double timeRange = difftime(endTime.Convert(),startTime.Convert());
  
  int ngen = 0;
  while ( ngen < gOptNevents ) {
    struct Time time;
    struct Position pos;

    // Generate a random time in our range
    time_t ctime = startTime.Convert() + timeRange * (double)rand() / RAND_MAX;
    
    // Convert to SolTrack struct format
    struct tm *ctimeInfo = localtime(&ctime);    
    time.year = ctimeInfo->tm_year + 1900;
    time.month = ctimeInfo->tm_mon + 1;
    time.day = ctimeInfo->tm_mday;
    time.hour = ctimeInfo->tm_hour;
    time.minute = ctimeInfo->tm_min;
    time.second = ctimeInfo->tm_sec;
    
    // Calculate the position of the sun
    SolTrack(time, loc, &pos, useDegrees, useNorthEqualsZero, computerRefrEquatorial, computeDistance);

    // Set the current entry
    fentry->Reset();
    fentry->metakey = metakey;
    fentry->pdg     = 2000010000;
    fentry->wgt     = 1.0;
    fentry->vtxx    = 0.0;
    fentry->vtxy    = 0.0;
    fentry->vtxz    = 0.0;
    fentry->dist    = 0.0;
    fentry->px      = - pd*cos(pos.altitudeRefract)*sin(pos.azimuthRefract - angle);
    fentry->py      = - pd*sin(pos.altitudeRefract);
    fentry->pz      = - pd*cos(pos.altitudeRefract)*cos(pos.azimuthRefract - angle);
    fentry->E       = Ed;

    // Fill the tree
    fluxntp->Fill();

    // Set the metadata for the first entry
    if (ngen == 0) {
      fmeta->AddFlavor(fentry->pdg);      
      fmeta->maxEnergy = fentry->E;
      fmeta->minWgt    = fentry->wgt;
      fmeta->maxWgt    = fentry->wgt;
      fmeta->metakey   = metakey;
      metantp->Fill();
    }
    
    ngen++;
  }

  // Write to the file
  file->cd();
  fluxntp->Write();
  metantp->Write();
  file->Close();

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
      PrintSyntax();
      exit(0);
  }

  // number of events
  if( parser.OptionExists('n') ) {
    gOptNevents = parser.ArgAsInt('n');
  } else {
    std::cout
      << "Unspecified number of events to generate - Using default of " << gOptNevents << std::endl;
  }

  // Output file name
  if( parser.OptionExists('o') ) {
    gOptOutFileName = parser.ArgAsString('o');
  } else {
    std::cout
      << "Unspecified output file name - Using default of " << gOptOutFileName << std::endl;
  }
    
  // Dark matter energy
  if( parser.OptionExists('e') ) {
    gOptDMEnergy = parser.ArgAsDouble('e');
  } else {
    std::cout << "Unspecified dark matter energy - Exiting" << std::endl;
    PrintSyntax();
    exit(1);
  }

  // Dark matter mass
  if( parser.OptionExists('m') ) {
    gOptDMMass = parser.ArgAsDouble('m');
  } else {
    std::cout << "Unspecified dark matter mass - Exiting" << std::endl;
    PrintSyntax();
    exit(1);
  }

  // date range
  if( parser.OptionExists('d') ) {
    string strdatran = parser.ArgAsString('d');
    vector<string> drange = utils::str::Split(strdatran, ",");
    assert(drange.size() == 2);
    gOptDateStart = atoi(drange[0].c_str());
    gOptDateEnd = atoi(drange[1].c_str());
    assert(gOptDateEnd > gOptDateStart && gOptDateStart > 10000);
  } else {
    std::cout << "Unspecified date range - Using default of the calendar year 2018" << std::endl;
  }

  // detector coordinates
  if( parser.OptionExists('p') ) {
    string strdetpos = parser.ArgAsString('p');
    vector<string> dpos = utils::str::Split(strdetpos, ",");
    assert(dpos.size() == 2);
    gOptDetLat = atof(dpos[0].c_str());
    gOptDetLong = atof(dpos[1].c_str());
  } else {
    std::cout << "Unspecified detector location - Using default of DUNE" << std::endl;
  }

  // Beam angle
  if( parser.OptionExists('b') ) {
    gOptBeamAngle = parser.ArgAsDouble('b');
  } else {
    std::cout << "Unspecified beam angle - Using default of DUNE" << std::endl;
  }

  std::cout
    << std::endl
     << "Running with options:" << std::endl;
  std::cout
    << "Number of events:     " << gOptNevents << std::endl;
  std::cout
    << "Dark Matter energy:   " << gOptDMEnergy << std::endl;
  std::cout
    << "Dark Matter mass:     " << gOptDMMass << std::endl;
  std::cout
    << "Output filename:      " << gOptOutFileName << std::endl; 
  std::cout
    << "Date range:           " << gOptDateStart << " - " << gOptDateEnd << std::endl;
  std::cout
    << "Detector coordinates: " << "( " << gOptDetLat << " , " << gOptDetLong << " )" << std::endl;
  std::cout
    << "Beam angle:           " << gOptBeamAngle << std::endl;
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  std::cout
    << "\n\n" << "Syntax:" << "\n"
    << "\n      gevgen [-h]"
    << "\n              -n nev"
    << "\n              -e energy"
    << "\n              -m mass"
    << "\n             [-o outfile_name]"
    << "\n             [-d date_range]"
    << "\n             [-p detector_coords]"
    << "\n             [-b beam angle]"
    << "\n";
}
//____________________________________________________________________________
