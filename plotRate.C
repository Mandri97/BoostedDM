#include <iostream>
#include <fstream>

#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TVectorT.h>

using namespace std;

void helper(){
    printf("Usage: ./plotRate listOfAnalyzedFiles.txt\n");
    exit(-1);
}

void FILE_NOT_FOUND(const char *filename){
    printf("Could not find %s, exiting the program.\n", filename);
    exit(-1);
}

void plotRate(char* filename){
    // Create a graph 
    TGraph *graph = new TGraph();
    
    int fileCounter = 1;

    // Open txt file
    ifstream _inFile;
    _inFile.open(filename);
    if (!_inFile.is_open()) FILE_NOT_FOUND(filename);

    string line;
    while (getline(_inFile, line)){
        // Open rootfile
        TFile* _rootfile = new TFile(line.c_str());
        if (_rootfile->IsZombie()) FILE_NOT_FOUND(line.c_str());

        // Get Tree
        auto t_selectedEvts = (TTree*) _rootfile->Get("SelectedEvents");

        // Create histogram
        auto hEnergySelectedEvents = new TH1D("hEnergySelectedEvents", "", 1000, 0, 100);
        
        // Fill histogram
        t_selectedEvts->Draw("E >> hEnergySelectedEvents", "", "goff");

        // Get runtime
        double runtime  = ((TVectorT<double> *) _rootfile->Get("runtime"))->Max();
        double binWidth = hEnergySelectedEvents->GetBinWidth(1);
        //double aliveSeg = 126;

        double scaling = 1 / (runtime * binWidth); //* aliveSeg * 24.2);

        hEnergySelectedEvents->Scale(scaling);

        // Integration start and end energy
        double startEnergy = 15;
        double endEnergy = 50;

        // Calculate the integration bin
        int binMin = hEnergySelectedEvents->GetXaxis()->FindBin(startEnergy);
        int binMax = hEnergySelectedEvents->GetXaxis()->FindBin(endEnergy);
        binMax = (binMax > hEnergySelectedEvents->GetNbinsX()) ? hEnergySelectedEvents->GetNbinsX() : binMax;

        // Integrate
        double rate = hEnergySelectedEvents->Integral(binMin, binMax, "width");

        graph->SetPoint(fileCounter -1, fileCounter, rate);
        
        fileCounter++;

        // Close file to avoid memory leak
        _rootfile->Close();
    }
    // Free memory
    _inFile.close();

    TCanvas *c = new TCanvas("c", "", 800, 500);
    graph->Draw();
    c->SaveAs("RateMuon.pdf");
}

int main(int argc, char* argv[]){
    // Check the number of argument
    if (argc != 2) helper();

    plotRate(argv[1]);

    return 0;
}
