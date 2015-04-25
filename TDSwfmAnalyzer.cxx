#include "TStyle.h"
#include "TTree.h"
#include "TList.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1D.h"
#include "TKey.h"
#include "TNtuple.h"
#include "TGClient.h"
#include "TGFrame.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TColor.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TApplication.h"
#include "TRootBrowser.h"
#include "TBrowser.h"
#include "TVector.h"
#include "TProfile.h"
#include "TDSwfmAnalyzer.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <cctype>
#include <cmath>


/****************************************************************************
 ***      This function is used to print the header for each new segment  ***
 ****************************************************************************/
void printTitle()
{
  std::cout << "***********************************\n";
  std::cout << "***\tDPO Scope Analyzer\t***\n";
  std::cout << "***********************************\n";
  return;
}

/****************************************************************************
 ***  This function is used to create the nTuple from the raw data file   ***
 ****************************************************************************/
void nTupler(int numberOfArgs, char* arrayOfArgs[])
{
  system("clear");
  printTitle();
  
  std::cout << "\n\tNTupling the data...\n";
  TFile outputFile("STEP1_nTuples.root","RECREATE");
  
  /******************************************************************************************
   ***                      This loop runs over each input data file                      ***
   ******************************************************************************************/
  for (int whichArg = 1; whichArg != numberOfArgs; ++whichArg)
  {
    /******************************************************************************************
     ***                For handling data processed in another directory.                   ***
     ***        Turns input argument '/data/file/CosmicRun.xml' to 'CosmicRun.xml'          ***
     ******************************************************************************************/
    char *last_slash;
    bool isSlashed = 0;
    if (strrchr(arrayOfArgs[whichArg], '/'))
    {
      last_slash = strrchr(arrayOfArgs[whichArg], '/');
      isSlashed = 1;
    }
    else
    {
      last_slash = arrayOfArgs[whichArg];
    }
    
    if (isSlashed)
    {
      last_slash = (last_slash + 1);
    }
    /******************************************************************************************/

    
    std::cout << "\t\tProcessing: " << arrayOfArgs[whichArg] << std::endl;
    std::stringstream ss_dataName;  //!---------------------------------------------------------  "Title" of the Tree
    std::stringstream ss_ntupleName;  //!-------------------------------------------------------  "Name" of the Tree
    ss_dataName.clear();
    ss_ntupleName.clear();
    ss_dataName << last_slash;
    ss_ntupleName << "ntuple_" << whichArg;  //!------------------------------------------------  Each input file goes into own tree ntuple_1, ntuple_2...
    
    Double_t eachTimeSlice[4] = {0.2}; //!-------------------------------------------------------------------------------------  In nanoseconds
    Double_t verticalOffSet[4] = {0.0};
    std::vector<Float_t> timeSlice[4], amplitude[4];
    TTree theTree(ss_ntupleName.str().c_str(),ss_dataName.str().c_str());  //!------------------  Create the main output Tree, name ntuple_x, title CosmicRun.xml, e.g.
    theTree.SetAutoSave(10000000);  //!---------------------------------------------------------  This prevents duplicate entries in the output TFile.
    
    
    std::ifstream infile(arrayOfArgs[whichArg]);
    std::string headerString = arrayOfArgs[whichArg];
    std::string headerFile = headerString.substr(0, headerString.size()-3);
    headerFile.append("hea");
    
    std::ifstream headerInFile(headerFile);
    
    std::string currentLine;
    
    int numberOfEvents    = 0;  //!-------------------------------------------------------------  Container for "output every" calculation.
    int howManyGSamples   = 0;
    int howManySamples    = 0;
    int howManyChannels   = 0;
    int theFirstRun       = 0;
    /******************************************************************************************
     ***            This is the event processing loop, and loops over the data              ***
     ******************************************************************************************/
    while (std::getline(infile, currentLine).good())
    {
      if (numberOfEvents % 500 == 0)
      {
        std::cout << "Processing Event " << numberOfEvents << "." << std::endl;
      }
      ++numberOfEvents;
      
//      For Debug Purposes:
//      if (numberOfEvents > 1500)
//      {
//        break;
//      }
      
      std::istringstream parsedLine(currentLine);
      
      std::vector <std::string> entries;
      copy(std::istream_iterator<std::string>(parsedLine),std::istream_iterator<std::string>(),std::back_inserter<std::vector<std::string> >(entries));
      
      if (theFirstRun == 0)
      {
        for (int theChannels = 1; theChannels != 5; ++theChannels)
        {
          if ((entries.size() - 1032*theChannels) == 0)
          {
            howManySamples = 1032;
            howManyChannels = theChannels;
            break;
          }
          else if((entries.size() - 2032*theChannels) == 0)
          {
            howManySamples = 2032;
            howManyChannels = theChannels;
            break;
          }
          else if((entries.size() - 3032*theChannels) == 0)
          {
            howManySamples = 3032;
            howManyChannels = theChannels;
            break;
          }
          else if((entries.size() - 4032*theChannels) == 0)
          {
            howManySamples = 4032;
            howManyChannels = theChannels;
            break;
          }
          else if((entries.size() - 5032*theChannels) == 0)
          {
            howManySamples = 5032;
            howManyChannels = theChannels;
            break;
          }
          if (theChannels == 4)
          {
            std::cout << "\tCould not determine samples, gigasamples, or channels.  PLEASE HALP!" << std::endl;
            std::cout << "\t\tHow many channels are there?: ";
            std::cin >> howManyChannels;
            std::cout << "\n\t\tHow many gigasamples are there?: ";
            std::cin >> howManyGSamples;
            std::cout << "\n\t\tHow many samples re there?: ";
            std::cin >> howManySamples;
          }
        }
        
        for (int whichChannel = 0; whichChannel != howManyChannels;  ++whichChannel)
        {
          //Read the Header File
          int whichLineHasTimeSlice = whichChannel*13 + 10;
          int whichLineHasOffset = whichChannel*13 + 11;
          int lineCounter = 1;
          while (std::getline(headerInFile, currentLine).good())
          {
            if (whichLineHasTimeSlice == lineCounter)
            {
              eachTimeSlice[whichChannel] = atof((currentLine.substr(9,currentLine.size())).c_str())/1e-9; //Reads in the horizontal time slice and converts to ns from seconds.
            }
            if (whichLineHasOffset == lineCounter)
            {
              verticalOffSet[whichChannel] = atof((currentLine.substr(7,currentLine.size())).c_str())*100; //Reads in the vertical offset and converts to mV.
            }
            ++lineCounter;
          }
          headerInFile.clear();
          headerInFile.seekg (0, std::ios::beg);
          
          std::cout << "So many lines: " << lineCounter << std::endl;
          
          std::cout << "Channel " << whichChannel << " has offset of " << verticalOffSet[whichChannel] << " mV." << std::endl;
          std::cout << "Channel " << whichChannel << " has timeSlice of " << eachTimeSlice[whichChannel] << " ns." << std::endl;
          
        }
        
        for (int whichChannel = 0; whichChannel != howManyChannels;  ++whichChannel)
        {
          std::stringstream amplitudeBranchName;
          std::stringstream timeSliceBranchName;
          
          amplitudeBranchName << "Channel_" << whichChannel << "_Amplitude";
          timeSliceBranchName << "Channel_" << whichChannel << "_TimeSlice";
          
          theTree.Branch(amplitudeBranchName.str().c_str(),"std::vector<Float_t>",&amplitude[whichChannel]);
          theTree.Branch(timeSliceBranchName.str().c_str(),"std::vector<Float_t>",&timeSlice[whichChannel]);
          
          amplitudeBranchName.str("");
          timeSliceBranchName.str("");
        }
        theFirstRun = 1;
      }
      
//      std::cout << "\tFound: " << std::endl;
//      std::cout << "\t\t" << howManyChannels << " Channels." << std::endl;
//      std::cout << "\t\t" << howManySamples << " Samples." << std::endl;
//      std::cout << "\t\t" << howManyGSamples << " GS/s." << std::endl;
//      
//      if (entries[0]!="0") std::cout << "Found Channel 1. ";
//      if (entries[1032]!="0") std::cout << "Found Channel 2." << std::endl;
//      
//      
//      std::cout << "NEW EVENT, Channel 1:" << std::endl;

      
      for (int whichChannel = 0; whichChannel != howManyChannels;  ++whichChannel)
      {
        amplitude[whichChannel].clear();
        timeSlice[whichChannel].clear();
        for (int whichTimeSlice = 0; whichTimeSlice != howManySamples; ++whichTimeSlice)
        {
          amplitude[whichChannel].push_back(atof(entries[whichTimeSlice+whichChannel*howManySamples].c_str()) - verticalOffSet[whichChannel]);
          timeSlice[whichChannel].push_back(whichTimeSlice*eachTimeSlice[whichChannel]);
        }
      }
      theTree.Fill();
    }
    
    std::cout << "\t\tFinished Processing " << numberOfEvents << " Events.\n\n";
    
    if (whichArg == (numberOfArgs-1))
    {
      Int_t saveGains = 0;
      
      std::cout << "\tSave PMT Gains? 1 or 0: ";
      std::cin >> saveGains;
      
      std::cin.clear();
      std::cin.sync();
      std::cin.ignore();
      
      std::cout << std::endl;
      
      /******************************************************************************************************
       ***                  This is currently used in the analyzer for calculating P.E.s                  ***
       ******************************************************************************************************/
      if (saveGains == 1)
      {
        std::stringstream ss_vectorName;
        ss_vectorName << "PMTGains";
        
        Double_t pmtGains;
        std::string pmtSerials;
        
        std::vector<std::string> pmtSerialNumbers;
        std::vector<std::string> *branchpmtSerialNumbers = &pmtSerialNumbers;
        
        std::vector<Double_t> pmtGainNumbers;
        std::vector<Double_t> *branchpmtGainNumbers = &pmtGainNumbers;
        
        TTree *pmtTree = new TTree("PMTs","PMT Gains and Serials");
        pmtTree->Branch("pmtSerials", &branchpmtSerialNumbers);
        pmtTree->Branch("pmtGains",&branchpmtGainNumbers);

        for (int whichRun = 0; whichRun != howManyChannels; ++whichRun)
        {
          std::cout << "\tPMT Serial Channel " << whichRun+1 << ": ";
          std::cin >> pmtSerials;
          
          std::cin.sync();
          std::cin.ignore(256,'\n');
          
          std::cout << "\tPMT Gain Channel " << whichRun+1 << ": ";
          std::cin >> pmtGains;
          
          std::cout << std::endl;
          
          std::cin.ignore(256,'\n');
          
          pmtSerialNumbers.push_back(pmtSerials);
          pmtGainNumbers.push_back(pmtGains);
          
          pmtTree->Fill();
        }
        pmtTree->Write("PMTs");
      }
    }
    theTree.Write();
    std::cout << "\n\t\t* nTuple " << theTree.GetName() << " has been successfully created.\n\n\n";
  }
  
//  std::cin.clear();
//  std::cin.sync();
//  std::cin.ignore();
  
//  outputFile.Write();  //!----------------------------------------------------------------------  Sometimes you need to write, sometimes you don't.
  outputFile.Close();
}


/***************************************************************************
 ***                This function processes the ntuple                   ***
 ***************************************************************************/
void Analyze(int numberOfKeys, int numberOfWaveformsToPlot, double triggerLevel)
{
  system("clear");
  printTitle();
  
  TFile *inFile  = new TFile("STEP1_nTuples.root","READ");
  TFile outputFile("STEP2_Analyzed.root","RECREATE");

  Int_t howManyKeysToRunOver = 0;

  
  if ((TTree*)inFile->GetListOfKeys()->Contains("PMTGains"))
  {
    howManyKeysToRunOver = numberOfKeys-1;
  }
  else
  {
    howManyKeysToRunOver = numberOfKeys+1;
  }
  
  const Double_t chargeOfElectron = 1.602e-19;
  Double_t peCalculation = 0.0;
  Int_t numberOfTimeSlices = 2032;
  
  TH1F *eventNumberHistogram = new TH1F("EventsPerRun","Events per Run",howManyKeysToRunOver,1,howManyKeysToRunOver+1);
  
  for (int whichKey = 1; whichKey != howManyKeysToRunOver; ++whichKey)
  {
    std::cout << "********* THIS KEY IS " << whichKey << " OUT OF " << howManyKeysToRunOver << " ****************\n";
    
    std::stringstream whichKeyToGet;
    whichKeyToGet << "ntuple_" << whichKey;

    TTree *theTree = (TTree*)inFile->Get(whichKeyToGet.str().c_str());
    TTree *pmtTree = (TTree*)inFile->Get("PMTs");
    
    const char* dir = theTree->GetTitle();
    
    TDirectory *newDirectory = outputFile.mkdir(dir);
    newDirectory->cd();
    
    std::vector<Float_t> *timeSlice = new std::vector<Float_t>;
    std::vector<Float_t> *amplitude = new std::vector<Float_t>;

    std::vector<std::string> *pmtSerialNumber = 0;
    std::vector<Double_t> *thePMTGains = 0;
    
    TBranch * theTimeSliceBranch;
    TBranch * theAmplitudeBranch;
    
    Double_t thePMTGainsForCalcs[4];
    
    int numberOfBranches = 0;
    int hasPMTGains = 0;
    
    if (inFile->GetListOfKeys()->Contains("PMTs"))
    {
      pmtTree->SetBranchAddress("pmtSerials", &pmtSerialNumber);
      pmtTree->SetBranchAddress("pmtGains", &thePMTGains);
      pmtTree->GetEntries();
      
      for (int pmtEntry = 0; pmtEntry != pmtTree->GetEntries(); ++pmtEntry)
      {
        pmtTree->GetEntry(pmtEntry);
        thePMTGainsForCalcs[pmtEntry] = thePMTGains->at(pmtEntry);
      }
      
      numberOfBranches = (theTree->GetNbranches())/2;
      hasPMTGains = 1;
    }
    else
    {
      numberOfBranches = (theTree->GetNbranches()/2);
    }
    
    /****************************************************
     ***            Book the Histograms               ***
     ****************************************************/
    TH1F *amplitudeHistograms[numberOfBranches];
    TH1F *amplitudePedestalSubtractedHistograms[numberOfBranches];
    TH1F *integralHistograms[numberOfBranches];
    TH1F *waveFormHistograms[numberOfBranches];
    TH1F *waveFormHistogramsPedestalSubtracted[numberOfBranches];
    TH1F *waveFormPeakHistogram[numberOfBranches];
    TH1F *peHistogram[numberOfBranches];
    TProfile *waveformProfile[numberOfBranches];
    TProfile *waveformProfilePedestalSubtracted[numberOfBranches];
    
    /****************************************************
     ***            Book any fit functions            ***
     ****************************************************/
    TF1 *fitresultA[numberOfBranches];
    TF1 *fitresultB[numberOfBranches];
    TF1 *pedestalFits[numberOfBranches];
    
    Double_t CutOff = 0.0;
    
    Double_t pedestalValues[2][numberOfBranches];  //!--------------------------- 0 is Mean, 1 is StdDev, 2 is cutoff (par[0]+5*[1]) 5 sigma
    for (Int_t thisBranch = 0; thisBranch != numberOfBranches; ++thisBranch)
    {
      newDirectory->cd();
      std::stringstream channelDirectoryName;
      channelDirectoryName << "Channel_" << thisBranch << "_Analysis";
      TDirectory *channelAnalysisDirectory = gDirectory->mkdir(channelDirectoryName.str().c_str());
      channelAnalysisDirectory->cd();
      TDirectory *waveFormDirectory = gDirectory->mkdir("WaveForms");
      
      /*****************************************************************
       ***      Set the branches and histogram Names and Titles      ***
       *****************************************************************/
      std::stringstream timeSliceBranchName;
      timeSliceBranchName << "Channel_" << thisBranch << "_TimeSlice";
      
      std::stringstream amplitudeBranchName;
      amplitudeBranchName << "Channel_" << thisBranch << "_Amplitude";
      
      std::stringstream amplitudePedestalSubtractedHistName;
      amplitudePedestalSubtractedHistName << "Channel_" << thisBranch << "_Pedestal_Subtracted_Amplitude";
      
      std::stringstream pedestalFitName;
      pedestalFitName << "Channel_" << thisBranch << "_Pedestal";
      
      std::stringstream integralHistogramName;
      integralHistogramName << "Channel_" << thisBranch << "_Integral";
      
      std::stringstream waveformHistogramName;
//      waveformHistogramName << "Channel_" << thisBranch << "_WaveformHistogram";
      
      std::stringstream waveformHistogramTitle;
//      waveformHistogramTitle << "Channel " << thisBranch << " Waveform Histogram";
      
      std::stringstream waveFormHistogramPedSubName;
      waveFormHistogramPedSubName << "Channel_" << thisBranch << "_WaveformHistogramPedSub";
      
      std::stringstream waveFormHistogramPedSubTitle;
      waveFormHistogramPedSubTitle << "Channel " << thisBranch << " Waveform Histogram Pedestal Subtracted";
      
      std::stringstream waveformProfilePedSubName;
      waveformProfilePedSubName << "Channel_" << thisBranch << "WaveFormProfile_PedestalSubtracted";
      
      std::stringstream waveFormProfileName;
      waveFormProfileName << "Channel_" << thisBranch << "_WaveFormProfile";
      
      std::stringstream waveFormProfileTitle;
      waveFormProfileTitle << "Channel " << thisBranch << " WaveForm Profile";
      
      std::stringstream waveformProfilePedSubTitle;
      waveformProfilePedSubTitle << "Channel_" << thisBranch << "_WaveFormProfile_PedestalSubtracted";
      
      std::stringstream peakTitle;
      peakTitle << "Channel " << thisBranch << " Peaks";
      
      std::stringstream peakName;
      peakName << "Channel_" << thisBranch << "_Peaks";
      
      std::stringstream peName;
      if (hasPMTGains) peName << "Channel_" << thisBranch << "_PEs";
      
      std::stringstream peTitle;
      if (hasPMTGains) peTitle << "Channel " << thisBranch << " PEs";

      
      
      /****************************************************
       ***        Get the branches from the trees       ***
       ****************************************************/
      if (!(theTree->GetBranchStatus(timeSliceBranchName.str().c_str())))
      {
        continue;
      }
      theTimeSliceBranch = theTree->GetBranch(timeSliceBranchName.str().c_str());
      theAmplitudeBranch = theTree->GetBranch(amplitudeBranchName.str().c_str());
      
      
      /****************************************************
       ***          Set the branch addresses            ***
       ****************************************************/
      theTimeSliceBranch->SetAddress(&timeSlice);
      theAmplitudeBranch->SetAddress(&amplitude);
      
      
      /****************************************************
       ***        Containers for time calculations      ***
       ****************************************************/
      Double_t maximumValue = 0.0;  //!----------------------------------------------------------------------  Holds maximum amplitude of event
      Int_t maximumBin = 0;  //!-----------------------------------------------------------------------------  Holds bin location for max amplitude
      Int_t halfBin = 0;  //!--------------------------------------------------------------------------------  Holds bin location for last bin above half max amplitude
      Int_t eBin    = 0;  //!--------------------------------------------------------------------------------  Holds bin location for last bin above 1/e of max amplitude
      Int_t tenBin  = 0;  //!--------------------------------------------------------------------------------  Holds bin location for last bin above 1/10 of max amplitude
      Double_t nBins[numberOfTimeSlices+1];  //!-------------------------------------------------------------  There are 1032 samples out of DP0, nBins is 1032+1
      
      /********************************************
       ***          Loop to set bin edge        ***
       ********************************************/
      theAmplitudeBranch->GetEntry(1);
      theTimeSliceBranch->GetEntry(1);
      for (unsigned int amp = 0; amp != amplitude->size(); ++amp)
      {
        nBins[amp+1] = timeSlice->at(amp);
      }
      nBins[0] = 0.0;
      //***************************************//
      
      /****************************************************
       ***            Assign the Histograms             ***
       ****************************************************/
      amplitudeHistograms[thisBranch] = new TH1F(amplitudeBranchName.str().c_str(),amplitudeBranchName.str().c_str(),500,-100,300);
      amplitudeHistograms[thisBranch]->GetYaxis()->SetTitle("Events / 1 mV");
      amplitudeHistograms[thisBranch]->GetXaxis()->SetTitle("Amplitude (mV)");
      
      amplitudePedestalSubtractedHistograms[thisBranch] = new TH1F(amplitudePedestalSubtractedHistName.str().c_str(),amplitudePedestalSubtractedHistName.str().c_str(),500,-100,300);
      amplitudePedestalSubtractedHistograms[thisBranch]->GetYaxis()->SetTitle("Events / 1 mV");
      amplitudePedestalSubtractedHistograms[thisBranch]->GetXaxis()->SetTitle("Amplitude (mV)");
      
      integralHistograms[thisBranch] = new TH1F(integralHistogramName.str().c_str(), integralHistogramName.str().c_str(),1000,0,1000);
      integralHistograms[thisBranch]->GetYaxis()->SetTitle("Events / 200 fC");
      integralHistograms[thisBranch]->GetXaxis()->SetTitle("Charge (fC)");
      
      waveFormHistograms[thisBranch] = new TH1F(waveformHistogramName.str().c_str(), waveformHistogramTitle.str().c_str(), numberOfTimeSlices, nBins);
      waveFormHistograms[thisBranch]->GetYaxis()->SetTitle("Amplitude (mV)");
      waveFormHistograms[thisBranch]->GetXaxis()->SetTitle("Time (ns)");
      
      waveFormHistogramsPedestalSubtracted[thisBranch] = new TH1F(waveFormHistogramPedSubName.str().c_str(), waveFormHistogramPedSubTitle.str().c_str(), numberOfTimeSlices, nBins);
      waveFormHistogramsPedestalSubtracted[thisBranch]->GetYaxis()->SetTitle("Amplitude (mV)");
      waveFormHistogramsPedestalSubtracted[thisBranch]->GetXaxis()->SetTitle("Time (ns)");
      
      waveformProfile[thisBranch] = new TProfile(waveFormProfileName.str().c_str(), waveFormProfileTitle.str().c_str(), numberOfTimeSlices, nBins);
      waveformProfile[thisBranch]->GetYaxis()->SetTitle("Average mV");
      waveformProfile[thisBranch]->GetXaxis()->SetTitle("Time (ns)");
      
      waveformProfilePedestalSubtracted[thisBranch] = new TProfile(waveformProfilePedSubTitle.str().c_str(), waveformProfilePedSubTitle.str().c_str(), numberOfTimeSlices, nBins);
      waveformProfilePedestalSubtracted[thisBranch]->GetYaxis()->SetTitle("Average mV");
      waveformProfilePedestalSubtracted[thisBranch]->GetXaxis()->SetTitle("Time (ns)");
      
      waveFormPeakHistogram[thisBranch] = new TH1F(peakName.str().c_str(), peakTitle.str().c_str(),numberOfTimeSlices,nBins);
      
      if (hasPMTGains)
      {
        peHistogram[thisBranch] = new TH1F(peName.str().c_str(), peTitle.str().c_str(),1000,0,100);
        peHistogram[thisBranch]->GetYaxis()->SetTitle("Events / 1 pe");
        peHistogram[thisBranch]->GetXaxis()->SetTitle("# of P.E.s");
      }
      
      
      /****************************************************
       ***        Define the pedestal fit function      ***
       ****************************************************/
      pedestalFits[thisBranch] = new TF1(pedestalFitName.str().c_str(),"gaus(0)",-10,30);
      
      Long64_t numberOfEntries = (Int_t)theTree->GetEntries();
//      Long64_t numberOfEntries = 5000;  //!----------------------------------------------------------------  For debugging purposes

      if (thisBranch == 0)
      {
        eventNumberHistogram->Fill(whichKey, numberOfEntries);
      }
      
      TCanvas *forWaveForms = new TCanvas("Waveform","Waveforms", 800, 600);  //!----------------------------  Canvas for waveforms
      forWaveForms->SetRightMargin(0.2413793);  //!----------------------------------------------------------  This puts the stats box off to the side of the plot.
      
      int plottedWaveforms = 0;
      /*******************************************************************
       ***          First loop over events, get Pedestal Values        ***
       *******************************************************************/
    
      std::cout << "\tProcessing Pedestals...\n";
      for (Long64_t theEntry = 0; theEntry != numberOfEntries; ++theEntry)
      {
        if (!(theEntry % 1000))
        {
          std::cout << "\t\tProcessing Event: " << theEntry << " of " << numberOfEntries << std::endl;
        }
        
        theAmplitudeBranch->GetEntry(theEntry);
        theTimeSliceBranch->GetEntry(theEntry);

        for (unsigned int amp = 0; amp != amplitude->size(); ++amp)
        {
          nBins[amp+1] = timeSlice->at(amp);
        }
        nBins[0] = 0.0;

        int triggerLevelReached = 0;
        
        for (unsigned int amp = 0; amp != amplitude->size(); ++amp)
        {
          if (abs(amplitude->at(amp)) > triggerLevel && timeSlice->at(amp) > 0 && timeSlice->at(amp) < 200)
          {
            triggerLevelReached = 1;
          }
          amplitudeHistograms[thisBranch]->Fill(amplitude->at(amp));
          waveFormHistograms[thisBranch]->Fill(timeSlice->at(amp), amplitude->at(amp));
        }
        
        if (thisBranch == 0)
        {
          CutOff = 8.0;
        }
        
        if (thisBranch != 0)
        {
          CutOff = 100.0;
        }
        
        if (waveFormHistograms[thisBranch]->GetBinContent(waveFormHistograms[thisBranch]->GetMaximumBin()) > CutOff)
        {
          gStyle->SetOptStat(1111111);
          // Set stat options
          gStyle->SetStatY(0.9);
          // Set y-position (fraction of pad size)
          gStyle->SetStatX(0.98);
          // Set x-position (fraction of pad size)
          gStyle->SetStatW(0.2);
          // Set width of stat-box (fraction of pad size)
          gStyle->SetStatH(0.2);
          // Set height of stat-box (fraction of pad size)
          gStyle->SetPadRightMargin(0.5);
          
          waveformHistogramName << "Waveform_" << plottedWaveforms+1 << ".pdf";
          waveformHistogramTitle << "Waveform " << plottedWaveforms+1;
          waveFormHistograms[thisBranch]->SetTitle(waveformHistogramTitle.str().c_str());
          waveFormHistograms[thisBranch]->Draw();
          waveFormHistograms[thisBranch]->GetYaxis()->SetRangeUser(0,120);
          waveFormHistograms[thisBranch]->GetYaxis()->SetTitleOffset(1.3);
          
          maximumValue = waveFormHistograms[thisBranch]->GetBinContent(waveFormHistograms[thisBranch]->GetMaximumBin());
          maximumBin   = waveFormHistograms[thisBranch]->GetMaximumBin();
          halfBin = waveFormHistograms[thisBranch]->FindLastBinAbove(maximumValue/2.0);
          eBin    = waveFormHistograms[thisBranch]->FindLastBinAbove(maximumValue/2.71828);
          tenBin  = waveFormHistograms[thisBranch]->FindLastBinAbove(maximumValue/10);
          
          waveFormPeakHistogram[thisBranch]->Fill(nBins[maximumBin]);
          
          /*******************************************
           ***  Draw markers for Peak calculations ***
           *******************************************/
          //Peak Marker
          TMarker *peakMarker = new TMarker(nBins[maximumBin],maximumValue+2,23);
          peakMarker->SetMarkerColor(kRed-3);
          peakMarker->Draw();
          
          //Peak to 1/2 marker
          TMarker *halfMarker = new TMarker(nBins[halfBin]+2,waveFormHistograms[thisBranch]->GetBinContent(halfBin)+2,23);
          halfMarker->SetMarkerColor(833);
          halfMarker->Draw();
          
          //Peak to 1/e marker
          TMarker *eMarker = new TMarker(nBins[eBin]+2,waveFormHistograms[thisBranch]->GetBinContent(eBin)+2,23);
          eMarker->SetMarkerColor(kAzure-1);
          eMarker->Draw();
          
          //Peak to 1/10 marker
          TMarker *tenMarker = new TMarker(nBins[tenBin]+2,waveFormHistograms[thisBranch]->GetBinContent(tenBin)+2,23);
          tenMarker->SetMarkerColor(kOrange-3);
          tenMarker->Draw();
          
          Double_t halfDecay  = (nBins[halfBin] - nBins[maximumBin]);
          Double_t eDecay     = (nBins[eBin] - nBins[maximumBin]);
          Double_t tenDecay   = (nBins[tenBin] - nBins[maximumBin]);
          
          forWaveForms->Update();
          
          std::stringstream halfDecayText;
          std::stringstream eDecayText;
          std::stringstream tenDecayText;
          std::stringstream peakBinText;
          
          halfDecayText << "1/2 = " << halfDecay << " ns";
          eDecayText << "1/e = " << eDecay << " ns";
          tenDecayText << "1/10 = " << tenDecay << " ns";
          peakBinText << "Peak = " << nBins[maximumBin] << " ns";

          TPaveText *stats = (TPaveText*)gPad->GetPrimitive("stats");
          stats->SetName("Mystats");

          stats->AddText(peakBinText.str().c_str());
          stats->AddText(halfDecayText.str().c_str());
          stats->AddText(eDecayText.str().c_str());
          stats->AddText(tenDecayText.str().c_str());

          stats->GetLineWith(peakBinText.str().c_str())->SetTextColor(kRed-3);
          stats->GetLineWith(halfDecayText.str().c_str())->SetTextColor(833);
          stats->GetLineWith(eDecayText.str().c_str())->SetTextColor(kAzure-1);
          stats->GetLineWith(tenDecayText.str().c_str())->SetTextColor(kOrange-3);

          gPad->Modified();
          gPad->Update();
          
//          forWaveForms->SaveAs(waveformHistogramName.str().c_str());
          waveFormDirectory->cd();
          if (plottedWaveforms < 100)
          {
            waveFormHistograms[thisBranch]->Write(waveformHistogramTitle.str().c_str());
          }
          
          waveformHistogramName.str(std::string());
          waveformHistogramTitle.str(std::string());

          ++plottedWaveforms;
        }
        waveFormHistograms[thisBranch]->Reset();
        
      }
      std::cout << "\nFinished Processing Pedestals Successfully.\n\n";
      
      amplitudeHistograms[thisBranch]->Fit("gaus","q","",-15,25);
      pedestalFits[thisBranch] = amplitudeHistograms[thisBranch]->GetFunction("gaus");
      pedestalValues[0][thisBranch] = pedestalFits[thisBranch]->GetParameter(1);
      pedestalValues[1][thisBranch] = pedestalFits[thisBranch]->GetParameter(2);
      pedestalValues[2][thisBranch] = pedestalValues[0][thisBranch]+5*pedestalValues[1][thisBranch];
      
      std::cout << "Mean is: " << pedestalValues[0][thisBranch] << std::endl;
      std::cout << "Sigma is: " << pedestalValues[1][thisBranch] << std::endl;
      std::cout << "Cutoff is: " << pedestalValues[2][thisBranch] << std::endl;
      
      /*******************************************************************
       ***          Second loop over events, Do everything else        ***
       *******************************************************************/
      std::cout << "\tNow, Processing Waveforms...\n";
      for (Long64_t theEntry = 0; theEntry != numberOfEntries; ++theEntry)
      {
        int triggerLevelReached = 0;
        int numberOfPointsOverTrigger = 0;
        
        if (!(theEntry % 1000))
        {
          std::cout << "\t\tProcessing Event: " << theEntry << " of " << numberOfEntries << std::endl;
        }
        
        theTimeSliceBranch->GetEntry(theEntry);
        theAmplitudeBranch->GetEntry(theEntry);
        waveFormHistograms[thisBranch]->Reset();
        waveFormHistogramsPedestalSubtracted[thisBranch]->Reset();
        
        Double_t newAmplitude = 0.0;
        
        for (unsigned int amp = 0; amp != amplitude->size(); ++amp)
        {
          // Pedestal subtract and zero suppresss amplitude
          newAmplitude = ((amplitude->at(amp) - 6.0) > 0 ? (amplitude->at(amp) - 6.0) : 0.0);
          amplitudePedestalSubtractedHistograms[thisBranch]->Fill(newAmplitude);
          waveformProfile[thisBranch]->Fill(timeSlice->at(amp),amplitude->at(amp));
          waveformProfilePedestalSubtracted[thisBranch]->Fill(timeSlice->at(amp), newAmplitude);
          
          if (newAmplitude > (pedestalValues[2][thisBranch]+10))
          {
            triggerLevelReached = 1;
            ++numberOfPointsOverTrigger;
          }
          
          waveFormHistogramsPedestalSubtracted[thisBranch]->Fill(timeSlice->at(amp), newAmplitude);
        }
        
        if (waveFormHistogramsPedestalSubtracted[thisBranch]->GetBinContent(waveFormHistogramsPedestalSubtracted[thisBranch]->GetMaximumBin()) > CutOff)
        {
          integralHistograms[thisBranch]->Fill(waveFormHistogramsPedestalSubtracted[thisBranch]->Integral("width"));
        }
        
        if (hasPMTGains && (waveFormHistogramsPedestalSubtracted[thisBranch]->GetBinContent(waveFormHistogramsPedestalSubtracted[thisBranch]->GetMaximumBin()) > CutOff))
        {
          peCalculation = (waveFormHistogramsPedestalSubtracted[thisBranch]->Integral("width")*1e-15 / chargeOfElectron)/(thePMTGainsForCalcs[thisBranch]);
          peHistogram[thisBranch]->Fill(peCalculation);
        }
      }
      std::cout << "\nFinished Processing Waveforms Successfully.\n\nNow, fitting histograms...\n";

      /********************************************************************
       ***              Process Histograms, Draw/Save/Write             ***
       ********************************************************************/
      channelAnalysisDirectory->cd();
      gStyle->SetPadRightMargin(0.25);
      
      amplitudeHistograms[thisBranch]->Draw();
      amplitudeHistograms[thisBranch]->Write();
      
      amplitudePedestalSubtractedHistograms[thisBranch]->Draw();
      amplitudePedestalSubtractedHistograms[thisBranch]->Write();
      
      integralHistograms[thisBranch]->Draw();
      integralHistograms[thisBranch]->Write();
      
      waveformProfile[thisBranch]->Draw();
      waveformProfile[thisBranch]->Write();
      
      waveformProfilePedestalSubtracted[thisBranch]->Draw();
      waveformProfilePedestalSubtracted[thisBranch]->Write();
      
      if (hasPMTGains)
      {
        peHistogram[thisBranch]->Draw();
        peHistogram[thisBranch]->Write();
      }
      
      waveFormPeakHistogram[thisBranch]->GetYaxis()->SetTitle("Number of Events");
      waveFormPeakHistogram[thisBranch]->GetXaxis()->SetTitle("Peak arrival time (ns)");
      
      
      /********************************************************************
       ***    The waveFormPeakHistogram plots arrival time of pulses    ***
       ***    I am using two fitting functions because there were two   ***
       ***    Pulses arriving at CERN test beam, and I wanted to        ***
       ***    differentiate between them.                               ***
       ********************************************************************/
      TF1 *gausA = new TF1("gausA","gaus",0,100);
      TF1 *gausB = new TF1("gausB","gaus",80,100);
      
      waveFormPeakHistogram[thisBranch]->Rebin(8);
      
      waveFormPeakHistogram[thisBranch]->Draw();
      gPad->Update();
      waveFormPeakHistogram[thisBranch]->Fit(gausA,"R");
//      waveFormPeakHistogram[thisBranch]->Fit(gausB,"R+");

      fitresultA[thisBranch] = waveFormPeakHistogram[thisBranch]->GetFunction("gausA");
//      fitresultB[thisBranch] = waveFormPeakHistogram[thisBranch]->GetFunction("gausB");
      
      TPaveText *otherStats = (TPaveText*)gPad->GetPrimitive("stats");
      otherStats->SetName("MyOtherstats");
      waveFormPeakHistogram[thisBranch]->SetStats(0);
      
      if (gausA->GetParameter(1))
      {
        Double_t mean1 = gausA->GetParameter(1);
        
        std::stringstream mean1Text;
        
        mean1Text << "Mean 1 = " << mean1 << " ns";
        
        otherStats->AddText(mean1Text.str().c_str());
        
      }
      
      if (gausB->GetParameter(1))
      {
        Double_t mean2 = gausB->GetParameter(1);
        
        std::stringstream mean2Text;
        
        mean2Text << "Mean 2 = " << mean2 << " ns";
        
        otherStats->AddText(mean2Text.str().c_str());
      }
      
      gPad->Modified();
      gPad->Update();

      waveFormPeakHistogram[thisBranch]->Write();
      
      peakTitle.str(std::string());
      peakName.str(std::string());
      
      std::cout << "\nFinished Fitting Waveforms.\n\n";
    }
  }
  
  outputFile.cd();
  eventNumberHistogram->Draw();
  eventNumberHistogram->Write();
  inFile->Close();
  outputFile.Close();
  
  return;
}


/****************************************************************************
 ***                        Print Colors                                  ***
 ****************************************************************************/
void set_plot_style()
{
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}
/*==========================================================================*/


/****************************************************************************
 ***                        The MAIN function                             ***
 ****************************************************************************/
int main(int argc, char* argv[])
{
  int numberOfArguments = argc;
  char* arguments[numberOfArguments];
  for (int whichArg = 0; whichArg != numberOfArguments; ++whichArg)
  {
    arguments[whichArg] = argv[whichArg];
  }

  TApplication *theApp = new TApplication("Analysis", &argc, argv);
  gROOT->ProcessLine("#include <vector>");
  gROOT->SetBatch();
  
  system("clear");
  
  int theChoices = -1;
  int breakCondition = 0;
  
  printTitle();
  
  /****************************************************************************
   ***                Parse input Files / run nTupler                       ***
   ****************************************************************************/
  if (numberOfArguments > 1)
  {
    std::cout << "\n\tThese files have been loaded into memory:\n";
    
    for (int whichArg = 1; whichArg != numberOfArguments; ++whichArg)
    {
      std::cout << "\t\t" << arguments[whichArg] << std::endl;
    }
    
    std::string createNtuple = "-1";
    breakCondition = 0;
    while (1)
    {
      std::cout << "\n\tCreate nTuple? Enter 1 for yes, 0 for no, q to quit: ";
      std::cin >> createNtuple;
      
      if (createNtuple == "q" || createNtuple == "Q" || createNtuple == ".q")
      {
        std::cout << "\n\t\t\033[1;31mProgram terminated.\033[0m\n\n";
        return 0;
      }
      else if (createNtuple == "0")
      {
        theChoices = 0;
        break;
      }
      else if (createNtuple == "1")
      {
        theChoices = 1;
        break;
      }
      else if (breakCondition == 5)
      {
        std::cout << "\n\t\t\033[1;31m Too many tries.  Program terminated.\033[0m\n\n";
        return 0;
      }
      else
      {
        std::cout << "\t\t\033[1;31m ERROR: Please enter a 1 or 0.\033[0m\n";
      }
      ++breakCondition;
    }
    
    if (theChoices == 1)
    {
      nTupler(numberOfArguments, arguments);
    }
  }
  /*============================================================================*/
  
  
  TFile *STEP1_nTuples  = new TFile("STEP1_nTuples.root","READ");
  if (!STEP1_nTuples->IsOpen())
  {
    std::cout << "\t\033[1;31mERROR: No STEP1_nTuples.root file.\033[0m" << std::endl;
    std::cout << "\tFirst please create the nTuple: ./DRS4Analyzer File1.txt File2.txt..." << std::endl;
    std::cout << "\tExiting\n\n";
    return 0;
  }
  
  TList *listOfKeys     = new TList(STEP1_nTuples->GetListOfKeys());
  if (!listOfKeys)
  {
    std::cout << "\n\t\033[1;31m Error:  No keys found in file.  Please create nTuple First.  Exiting.\033[0m\n";
    return 0;
  }
  std::cout << "\nFound ROOT file \"STEP1_nTuples.root\".  Reading...\n";
  
  
  /****************************************************************************
   ***                Handle Data Maniuplation Functions                    ***
   ****************************************************************************/
  while (1)
  {
    theChoices = -1;
    std::string makeChoice = "-1";
    TIter nextKey(STEP1_nTuples->GetListOfKeys());
    TKey *key;
    int whichKey      = 1;
    int numberOfKeys  = STEP1_nTuples->GetNkeys();
    std::cout << "\n\n\tThere are " << numberOfKeys << " Keys:\n";
    if (STEP1_nTuples->GetListOfKeys()->Contains("PMTs"))
    {
      numberOfKeys = numberOfKeys - 1;
    }
    const char *keyNames[numberOfKeys];
    while ((key = (TKey*)nextKey()))
    {
      keyNames[whichKey] = key->GetName();
      std::cout << "\t\t" << "Key " << whichKey << ":\t" << keyNames[whichKey] << "\t" << key->GetTitle() << std::endl;
      ++whichKey;
    }
    
    int selectedDatasetKeys[2] = {0};
    std::cout << "\n\tThe following functions are available:\n";
    std::cout << "\t\t1 - Analyze Waveforms";
//    std::cout << "\n\t\t2 - Plot Waveforms";
//    std::cout << "\n\t\t3 - Divide two waveforms:\tKey A / Key B";
//    std::cout << "\n\t\t4 - Percent difference:\t\t|Key A - Key B| / Key A\n";
//    std::cout << "\n\t\t5 - Create Histograms of each nTuple\n";
    std::cout << "\n\t\t6 - Open a TBrowser";
    std::cout << "\n\t\tq - Quit \n";
    //    std::cout << "\t\t4 - Percent difference of two waveforms";
    std::cout << "\n\tEnter which function you would like to perform: ";
    std::cin >> makeChoice;
    
    if (makeChoice == "q" || makeChoice == "Q" || makeChoice == "exit" || makeChoice == "Exit")
    {
      std::cout << "\n\t\t\033[1;31m Exiting the program.\033[0m \n\n";
      break;
    }
    else if (makeChoice == "1")
    {
      std::cout << "\n\tHow many waveforms to plot?: ";
      std::cin >> selectedDatasetKeys[0];
      theChoices = 1;
      
      if (theChoices ==1 )
      {
        std::cout << "\n\tWhat threshold for plotting waveforms? (in mV): ";
        std::cin >> selectedDatasetKeys[1];
      }
      
    }
//    else if (makeChoice == "2")
//    {
//      std::cout << "How many waveforms to plot?: ";
//      std::cin >> selectedDatasetKeys[0];
//      theChoices = 2;
//      plotWaveForms(selectedDatasetKeys[0]);
//      continue;
//    }
//    else if (makeChoice == "3")
//    {
//      theChoices = 3;
//    }
//    else if (makeChoice == "4")
//    {
//      theChoices = 4;
//    }
//    else if (makeChoice == "5")
//    {
//      theChoices = 5;
//    }
    else if (makeChoice == "6")
    {
      gROOT->SetBatch(kFALSE);
      TFile *inFile  = new TFile("STEP1_nTuples.root","READ");
      TFile *outputFile = new TFile("STEP2_Analyzed.root","READ");
      new RBrowser();
      theApp->Run(kTRUE);
      outputFile->Close();
      outputFile->Delete();
      inFile->Close();
      inFile->Delete();
      gROOT->SetBatch();
      continue;
    }
    else
    {
      std::cout << "\t\t\033[1;31m ERROR: You selected: " << makeChoice << ".  Please make a valid selection.\033[0m\n";
      continue;
    }
    
    Analyze(numberOfKeys, selectedDatasetKeys[0], selectedDatasetKeys[1]);
    
  }

  return 0;
}



















































































