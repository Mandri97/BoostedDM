#include <TFile.h>
#include <TVectorT.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;


double runingTime = 0.0;


void totalRunningTime(string fileName){
    auto _file = new TFile(fileName);

    runingTime += ((TVectorT<double>*) _file->Get("runtime"))->Max();

    _file->Close();
}

int main(int argc, char *argv[]){
    ifstream _textFile(argv[1]);

    string file;

    while (_textFile >> file) totalRunningTime(file);

    cout << "Total running time: " << runingTime << " s, " << runingTime / (60 * 60) << " h" << endl;    
    _textFile.close();
}
