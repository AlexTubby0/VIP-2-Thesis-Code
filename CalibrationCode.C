#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TFile.h"

// Function to plot histogram from multiple files and find peaks for a specific column
void plotHistogramAndFindPeaks(const std::vector<std::string>& filenames, int columnNumber, TH1F* combinedEnergyHist) {
    // Create a dynamic title for the histogram
    std::string histTitle = "Histogram of column " + std::to_string(columnNumber);

    // Create a ROOT histogram with a dynamic name and title
    int nBins = 4096; // Number of bins from 0 to 4096
    std::string histName = "hist_" + std::to_string(columnNumber);
    TH1F *hist = new TH1F(histName.c_str(), histTitle.c_str(), nBins, 0, 4096);

    // Process each file
    for (const auto& filename : filenames) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Unable to open file: " << filename << std::endl;
            continue;
        }

        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::vector<int> numbers;
            int number;
            while (iss >> number) {
                numbers.push_back(number);
            }
            if (numbers.size() >= 32) {
                int nthNumber = numbers[columnNumber - 1];
                int maxOfFirst32 = *std::max_element(numbers.begin(), numbers.begin() + 32);
                if (nthNumber == maxOfFirst32) {
                    hist->Fill(nthNumber);
                }
            }
        }
        file.close();
    }

    // Rebin the histogram by a factor of 8
    hist->Rebin(8);

    // Create a canvas and draw the histogram
    TCanvas *c1 = new TCanvas("c1", histTitle.c_str(), 800, 600);
    c1->SetLogy(); // Set Y-axis to logarithmic by default
    hist->GetXaxis()->SetTitle("Channels");
    hist->GetYaxis()->SetTitle("Count");
    hist->GetXaxis()->SetRangeUser(500, 1550); // Zoom into the region (500, 1550)
    hist->Draw("HIST"); // Draw histogram without red markers

    // Use TSpectrum to find peaks
    TSpectrum *spectrum = new TSpectrum();
    int nPeaks = spectrum->Search(hist, 2, "nodraw", 0.001); // Adjust parameters: sigma reduced, threshold lowered

    // Get peak positions and heights
    double *xPeaks = spectrum->GetPositionX();
    double *yPeaks = spectrum->GetPositionY();

    // Collect peaks in the desired range
    std::vector<std::pair<double, double>> peaks;
    for (int i = 0; i < nPeaks; i++) {
        double x = xPeaks[i];
        double y = yPeaks[i];
        if (x >= 790 && x <= 1550) { // Only consider peaks in the desired range
            peaks.emplace_back(x, y);
        }
    }

    // Sort peaks by height (y value) in descending order
    std::sort(peaks.begin(), peaks.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    // Select the most prominent peaks and ensure they are reasonably spaced apart
    std::vector<double> selectedPeaks;
    double minSpacing = 40.0; // Minimum spacing between peaks
    for (const auto& peak : peaks) {
        bool tooClose = false;
        for (const auto& selectedPeak : selectedPeaks) {
            if (std::abs(peak.first - selectedPeak) < minSpacing) {
                tooClose = true;
                break;
            }
        }
        if (!tooClose) {
            selectedPeaks.push_back(peak.first);
            if (selectedPeaks.size() == 5) break; // Only select up to 5 peaks
        }
    }

    // Sort the selected peaks by their x values (left to right)
    std::sort(selectedPeaks.begin(), selectedPeaks.end());

    // Define Gaussian fit functions and store their parameters
    std::vector<double> params;
    for (const auto& peak : selectedPeaks) {
        std::string funcName = "gaus_" + std::to_string(peak);
        TF1 *gaus = new TF1(funcName.c_str(), "gaus", peak - 20, peak + 20);

        hist->Fit(gaus, "RQ"); // 'Q' option for quiet mode

        // Get parameters: [0] = amplitude, [1] = mean, [2] = sigma
        double amp = gaus->GetParameter(0);
        double mean = gaus->GetParameter(1);
        double sigma = gaus->GetParameter(2);
        params.push_back(amp);
        params.push_back(mean);
        params.push_back(sigma);
    }

    // Define a multimodal Gaussian function with the new background and tail parameters
    TF1 *multiGaus = new TF1("multiGaus", [params](double *x, double *p) {
        double sum = 0;
        // Add the Gaussian components
        for (size_t i = 15; i < params.size() + 15; i += 3) {
            double amp = p[i];
            double mean = p[i+1];
            double sigma = p[i+2];
            sum += amp * exp(-0.5 * pow((x[0] - mean) / sigma, 2));
        }
        // Add the tail components for each Gaussian
        for (int j = 0; j < 5; ++j) {
            double amp = p[15 + j * 3];
            double mean = p[16 + j * 3];
            double sigma = p[17 + j * 3];
            double tailParam1 = p[4 + j * 2];
            double tailParam2 = p[5 + j * 2];
            sum += tailParam1 * amp * exp((x[0] - mean) / (tailParam2 * sigma) + 1 / (2 * tailParam2 * tailParam2)) *
                   TMath::Erfc((x[0] - mean) / (sqrt(2) * sigma) + 1 / (sqrt(2) * tailParam2));
        }
        // Add the new background function
        sum += p[0] * exp(-p[1] * x[0]) + p[2] * x[0] + p[3];
        return sum;
    }, 600, 1550, params.size() + 15);

    // Set the initial parameters for the new background function
    multiGaus->SetParameter(0, 75);    // A
    multiGaus->SetParameter(1, 0.005);  // B
    multiGaus->SetParameter(2, -0.0015); // C
    multiGaus->SetParameter(3, 11.5); // D
    // Set limits for the background parameter C to be negative
    multiGaus->SetParLimits(2, -10, -1e-5);

    // Set initial parameters for the new tail function parameters
    multiGaus->SetParameter(4, 0.04); // Initial value for tail amplitude factor for the first peak
    multiGaus->SetParameter(5, 3);     // Initial value for the 6 in the tail function for the first peak
    multiGaus->SetParameter(6, 0.03); // Initial value for tail amplitude factor for the second peak
    multiGaus->SetParameter(7, 3);     // Initial value for the 6 in the tail function for the second peak
    multiGaus->SetParameter(8, 0.05); // Initial value for tail amplitude factor for the third peak
    multiGaus->SetParameter(9, 3);     // Initial value for the 6 in the tail function for the third peak
    multiGaus->SetParameter(10, 0.04); // Initial value for tail amplitude factor for the fourth peak
    multiGaus->SetParameter(11, 3);     // Initial value for the 6 in the tail function for the fourth peak
    multiGaus->SetParameter(12, 0.04); // Initial value for tail amplitude factor for the fifth peak
    multiGaus->SetParameter(13, 3);     // Initial value for the 6 in the tail function for the fifth peak

    // Set limits for the tail function amplitude parameters to be between 0 and 0.08
    for (int i = 4; i <= 12; i += 2) {
        multiGaus->SetParLimits(i, 0, 0.08);
    }

    // Set limits for the tail function parameters to be between 0 and 5
    for (int i = 5; i <= 13; i += 2) {
        multiGaus->SetParLimits(i, 0, 5);
    }

    // Set the initial parameters for the multimodal Gaussian function
    for (size_t i = 0; i < params.size(); ++i) {
        multiGaus->SetParameter(i + 15, params[i]);
    }

    // Fit the global function to the histogram in the range 600-1550
    hist->Fit(multiGaus, "RQ", "", 600, 1550);

    // Draw the combined Gaussian function
    multiGaus->SetLineColor(kRed);
    multiGaus->SetNpx(2000); // Increase the number of points to make it smooth
    multiGaus->Draw("same");

    // Extract means from the global fit parameters
    std::vector<double> globalMeans;
    for (size_t i = 15; i < params.size() + 15; i += 3) {
        double mean = multiGaus->GetParameter(i + 1);
        globalMeans.push_back(mean);
    }

    // Create scatter plot of means vs. known energy values
    std::vector<double> energyValues = {4507.85, 4931.81, 5893.20, 6490.45, 8037.81};

    TCanvas *c2 = new TCanvas("c2", "Energy Calibration", 800, 600);
    TGraph *graph = new TGraph(globalMeans.size(), &globalMeans[0], &energyValues[0]);
    graph->SetTitle("Energy Calibration;Channels;Energy (eV)");
    graph->SetMarkerStyle(21);
    graph->Draw("AP");

    // Fit a second-order polynomial to the scatter plot
    TF1 *polyFit = new TF1("polyFit", "pol2", 600, 1550);
    graph->Fit(polyFit, "Q");

    // Draw the polynomial fit line
    polyFit->SetLineColor(kRed);
    polyFit->Draw("same");

    // Get the parameters of the polynomial fit
    double p0 = polyFit->GetParameter(0);
    double p1 = polyFit->GetParameter(1);
    double p2 = polyFit->GetParameter(2);

    // Print the equation of the polynomial fit line
    std::cout << "Equation of the polynomial fit line for column " << columnNumber << ": y = " << p0 << " + " << p1 << " * x + " << p2 << " * x^2" << std::endl;

    // Calculate chi-squared and NDF
    double chi2 = multiGaus->GetChisquare();
    int ndf = multiGaus->GetNDF();
    double chi2_ndf = chi2 / ndf;

    // Print chi-squared / NDF
    std::cout << "Chi-squared / NDF for column " << columnNumber << " = " << chi2_ndf << std::endl;

    // Create a new histogram for the energy values
    TH1F *energyHist = new TH1F(Form("energyHist_col_%d", columnNumber), "Energy Histogram", 800, 0, 20000); // Bin size of 25 (20000/800)

    // Process each file again to create the energy histogram
    for (const auto& filename : filenames) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Unable to open file: " << filename << std::endl;
            continue;
        }

        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::vector<int> numbers;
            int number;
            while (iss >> number) {
                numbers.push_back(number);
            }
            if (numbers.size() >= 32) {
                int nthNumber = numbers[columnNumber - 1];
                int maxOfFirst32 = *std::max_element(numbers.begin(), numbers.begin() + 32);
                if (nthNumber == maxOfFirst32) {
                    // Convert channel to energy using the polynomial fit
                    double energy = p0 + p1 * nthNumber + p2 * nthNumber * nthNumber;
                    if (energy >= 0 && energy <= 20000) {
                        energyHist->Fill(energy);
                        combinedEnergyHist->Fill(energy); // Add to the combined histogram
                    }
                }
            }
        }
        file.close();
    }

    // Draw the energy histogram for the current column
    TCanvas *c3 = new TCanvas(Form("c3_col_%d", columnNumber), "Energy Histogram", 800, 600);
    c3->SetLogy(); // Set Y-axis to logarithmic by default
    energyHist->GetXaxis()->SetTitle("Energy (eV)");
    energyHist->GetYaxis()->SetTitle("Count");
    energyHist->Draw("HIST");

    // Optionally save individual energy histograms
    // c3->SaveAs(Form("energy_histogram_col_%d.png", columnNumber));
}

int main(int argc, char **argv) {
    // Initialize the ROOT application
    TApplication app("app", &argc, argv);

    // Directory containing the data files
    std::string directory = "/mnt/c/Data_Background_December2023"; // Update with your actual directory path

    // Vector to hold file names
    std::vector<std::string> filenames;

    // Iterate through the directory and add text files to the filenames vector
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.path().extension() == ".txt") {
            filenames.push_back(entry.path());
        }
    }

    // Create a combined energy histogram
    TH1F *combinedEnergyHist = new TH1F("combinedEnergyHist", "Combined Energy Histogram", 800, 0, 20000);

    // Specify the columns to skip
    std::vector<int> skipColumns = {1, 13, 16, 26, 31};

    // Process each of the first 32 columns (skipping specified columns)
    for (int columnNumber = 1; columnNumber <= 32; ++columnNumber) {
        if (std::find(skipColumns.begin(), skipColumns.end(), columnNumber) == std::end(skipColumns)) {
            plotHistogramAndFindPeaks(filenames, columnNumber, combinedEnergyHist);
        }
    }

    // Draw the combined energy histogram
    TCanvas *c4 = new TCanvas("c4", "Combined Energy Histogram", 800, 600);
    c4->SetLogy(); // Set Y-axis to logarithmic by default
    combinedEnergyHist->GetXaxis()->SetTitle("Energy (eV)");
    combinedEnergyHist->GetYaxis()->SetTitle("Count");
    combinedEnergyHist->Draw("HIST");

    // Save the combined energy histogram as a PNG image and a ROOT file
    c4->SaveAs("combined_energy_histogram.png");
    TFile outFile("combined_energy_histogram.root", "RECREATE");
    combinedEnergyHist->Write();
    outFile.Close();

    // Run the ROOT application
    app.Run();

    return 0;
}
