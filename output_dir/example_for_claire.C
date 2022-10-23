//-----------------------------------------------------------------------------------------------//
// Example macro to load a histogram from a ROOT file.
//-----------------------------------------------------------------------------------------------//
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>

//-----------------------------------------------------------------------------------------------//
/* ROOT main macro *must* always be a void and must always be the same as the filename
 *
 * To run the macro:
 * $ root example_for_claire.C
 */
//-----------------------------------------------------------------------------------------------//
void example_for_claire()
{
    // Create TFile pointer and initialize it by opening the root file (here in read mode only)
    auto tfile = TFile::Open("run7MuonEvents.root", "read");

    // Create histogram pointer by getting the TH1 object inside the ROOT file.
    // The Get("object_name") function gets any object inside the ROOT file. Since this is a
    // 1D histogram with doubles we cast the object as a TH1D (every object in ROOT starts with a
    // "T", "H" is for histogram, "1" is for 1D, "D" is for double, if it was a 2D histogram of
    // floats it would be "TH2F"). Hence the (TH1D*) before tfile.
    // If this were a TTree, we could cast it as (TTree*), and so on...
    auto histogram = (TH1D *)tfile->Get("hMuonEvent");

    // Create a canvas pointer with name/title, with size 700 (width) by 500 (height) pixels
    auto canvas = new TCanvas("Muon Events with energy > 15 MeV", "Muon Events with energy > 15MeV", 700, 500);

    // Set automatic scientific notation on axes with large numbers
    TGaxis::SetMaxDigits(3);

    histogram->SetTitle("Muon Events with energy > 15 MeV");                 // Set a title for the plot
    histogram->GetYaxis()->SetTitle("Counts"); // Set the title for the Y axis
    histogram->GetXaxis()->SetTitle("Energy (MeV)");    // Set the title for the X axis
    histogram->GetYaxis()->SetRangeUser(0, 3000e3);      // Set the limits in the Y axis
   // histogram->GetXaxis()->SetLimits(15,30);
    histogram->GetXaxis()->SetRangeUser(15,700);          // Set the limits in the X axis
    histogram->SetLineColor(kBlue-2);                      // Set a color for the histogram's line. You can find the name of the colors and color schemes at https://root.cern.ch/doc/master/classTColor.html
    histogram->SetLineWidth(2);                         // Set the line width. Higher the integer, thicker the line
    histogram->Draw();                                  // Draw histogram. You can choose different draw styles by including a code in the parenthesis. The options are at https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html#drawing-histograms in section 5.8.2
}
