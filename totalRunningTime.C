#include <TFile.h>
#include <TVectorT.h>

void totalRunningTime(char *fileName){
    auto _file = new TFile(fileName);

    runingTime += ((TVectorT<double>*) _file->Get("runtime"))->Max();

    _file->Close();
}

int main(int argc, char *argv[]){
    double runingTime = 0.0;

    cout << "Total running time: " << runingTime << " s, " << runingTime / (60 * 60) << " h" << endl;    
}
