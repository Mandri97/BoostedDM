#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

using namespace std;

auto hLiveSegment = new TH1F("hLiveSegment", "Live segment", 154, 0, 154);

void liveSegment(char *inFile, char *outFile){
    auto _inFile  = new TFile(inFile);
    auto _outFile = new TFile(outFile, "recreate");

    auto tree = (TTree*) _inFile->Get("PhysPulse");
    int segment;

    tree->SetBranchAddress("seg", &segment);

    auto nEntries = tree->GetEntries();

    for (long int iEntry = 0; iEntry < nEntries; iEntry++){
        tree->GetEntry(iEntry);

        hLiveSegment->Fill(segment) ;
    }

    hLiveSegment->Write();

    _inFile->Close();
    _outFile->Close();
}

void helper(){
    cout << "Not enough argurments\n";
}


int main(int argc, char *argv[]){
    if (argc != 3){
        helper();

        return 1;
    }

    liveSegment(argv[1], argv[2]);
}
