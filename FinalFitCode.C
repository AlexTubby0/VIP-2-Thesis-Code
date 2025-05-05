#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TSystem.h>
#include <TSpectrum.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <TMath.h>
#include <TPaveText.h>
#include <TROOT.h>

void fitGaussianPeaks() {
    // Open the ROOT file
    TFile* file = TFile::Open("Dec_combined_energy_histogram.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // Retrieve the histogram from the file
    TH1F* hist = nullptr;
    file->GetObject("combinedEnergyHist", hist);
    if (!hist) {
        std::cerr << "Histogram not found in file" << std::endl;
        file->Close();
        return;
    }

    // Rename the histogram
    hist->SetTitle("Winter 2023 Parabolic Calibrated Data");

    // Check if the histogram has entries
    std::cout << "Number of entries in the histogram: " << hist->GetEntries() << std::endl;

    // Define the peak positions (approximate) and fit ranges
    std::vector<double> peakPositions = { 4500, 4900, 5850, 6400, 7400, 8000, 8750, 10500 };
    std::vector<double> fitRanges = { 70, 70, 70, 70, 70, 70, 70, 70 };

    // Create a canvas to draw the main histogram and fits
    TCanvas* c1 = new TCanvas("c1", "Summed Energy Histogram with Gaussian Fits", 800, 600);
    c1->cd();
    hist->GetXaxis()->SetTitle("Energy (eV)");
    hist->GetYaxis()->SetTitle("Counts / 8 eV");
    c1->SetLogy(); // Set logarithmic scale on y-axis

    // Set the range of the x-axis to zoom in
    hist->GetXaxis()->SetRangeUser(3500, 12500);

    gStyle->SetOptStat(0);
    hist->Draw("HIST");   // Draw the histogram
    hist->SetStats(0);    // Hide the stats box
    gPad->Modified();
    gPad->Update();       // Update the canvas

    // Use TSpectrum to refine peak positions within the specified ranges
    TSpectrum* s = new TSpectrum(10); // Increase peak buffer size to avoid "Peak buffer full" warnings

    // Initial parameter estimation
    std::vector<double> initialParams;
    for (size_t i = 0; i < peakPositions.size(); ++i) {
        double peak = peakPositions[i];
        double range = fitRanges[i];

        // Create a temporary histogram for the range around the peak
        double tempHistRange = 2 * range;
        int nBins = 2 * tempHistRange;
        TH1F* tempHist = new TH1F(Form("tempHist_%zu", i), "Temporary Histogram", nBins, peak - tempHistRange, peak + tempHistRange);

        // Copy bin contents to the temporary histogram
        for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
            double binCenter = hist->GetBinCenter(bin);
            if (binCenter >= peak - tempHistRange && binCenter <= peak + tempHistRange) {
                tempHist->SetBinContent(tempHist->FindBin(binCenter), hist->GetBinContent(bin));
            }
        }

        // Search for the peak in the temporary histogram
        int nPeaks = s->Search(tempHist, 2, "", 0.05);
        if (nPeaks == 0) {
            std::cerr << "No peak found in the range around " << peak << std::endl;
            delete tempHist;
            continue;
        }

        // Get the refined peak position
        double refinedPeak = s->GetPositionX()[0];

        // Define the Gaussian function with a tight initial range
        TF1* gaus = new TF1(Form("gaus_%zu", i), "gaus", refinedPeak - range, refinedPeak + range);

        // Fit the Gaussian to the temporary histogram first
        tempHist->Fit(gaus, "RQ");

        // Get refined parameters from the initial fit
        double amplitude = gaus->GetParameter(0);
        double mean = gaus->GetParameter(1);
        double sigma = gaus->GetParameter(2);

        // Ensure initial parameters are within reasonable bounds
        amplitude = std::max(amplitude, 0.0);
        sigma = std::min(std::max(sigma, 0.0), 200.0);

        // Add initial parameters to the vector
        initialParams.push_back(amplitude);
        initialParams.push_back(mean);
        initialParams.push_back(sigma);

        std::cout << "Initial parameters for Gaussian " << i + 1 << ":\n";
        std::cout << "  Amplitude: " << amplitude << "\n";
        std::cout << "  Mean: " << mean << "\n";
        std::cout << "  Sigma: " << sigma << "\n";

        // Draw the initial Gaussian fit
        gaus->SetLineColor(kMagenta);
        gaus->Draw("SAME");

        delete tempHist; // Clean up the temporary histogram
    }

    // Create the combined fit function with the new background and tail parameters
    TF1* multiGaus = new TF1("multiGaus", [](double* x, double* p) {
        double sum = 0;
        // Add the Gaussian components
        for (int i = 0; i < 8; ++i) {
            double amp = p[4 + 3 * i];
            double mean = p[5 + 3 * i];
            double sigma = p[6 + 3 * i];
            sum += amp * exp(-0.5 * pow((x[0] - mean) / sigma, 2));
        }
        // Add the tail components for the first four Gaussians
        for (int i = 0; i < 4; ++i) {
            double amp = p[4 + 3 * i];
            double mean = p[5 + 3 * i];
            double sigma = p[6 + 3 * i];
            double tailParam1 = p[28 + 2 * i];
            double tailParam2 = p[29 + 2 * i];
            sum += tailParam1 * amp * exp((x[0] - mean) / (tailParam2 * sigma) + 1 / (2 * tailParam2 * tailParam2)) *
                TMath::Erfc((x[0] - mean) / (sqrt(2) * sigma) + 1 / (sqrt(2) * tailParam2));
        }
        // Add the new background function
        sum += p[0] * exp(-p[1] * x[0]) + p[2] * x[0] + p[3];
        return sum > 0 ? sum : 1e-3; // Ensure the function value is always positive
        }, 4000, 12000, initialParams.size() + 12);

    // Set the initial parameters for the new background function
    multiGaus->SetParameter(0, 10000);   // A
    multiGaus->SetParameter(1, 0.0007);  // B
    multiGaus->SetParameter(2, -0.0001); // C
    multiGaus->SetParameter(3, 75);      // D

    // Set limits for the background parameter C to be negative
    multiGaus->SetParLimits(2, -0.01, 0);

    // Set the initial parameters for the Gaussians and tails
    for (size_t i = 0; i < initialParams.size(); ++i) {
        multiGaus->SetParameter(i + 4, initialParams[i]);
    }

    // Set initial parameters for the new tail function parameters
    for (int j = 0; j < 4; ++j) {
        multiGaus->SetParameter(28 + j * 2, 0.04); // Initial value for tail amplitude factor
        multiGaus->SetParameter(29 + j * 2, 3);    // Initial value for the 6 in the tail function
    }

    // Set limits for the tail function amplitude parameters to be between 0 and 0.08
    for (int i = 28; i < 36; i += 2) {
        multiGaus->SetParLimits(i, 0, 0.08);
    }

    // Set limits for the tail function parameters to be between 0 and 5
    for (int i = 29; i < 37; i += 2) {
        multiGaus->SetParLimits(i, 0, 5);
    }

    // Increase the number of points for the combined fit function to make it smooth
    multiGaus->SetNpx(2000);

    // Perform the combined fit
    hist->Fit(multiGaus, "RQ");

    // Draw the final individual Gaussian components using the parameters from the global fit
    for (int i = 0; i < 8; ++i) {
        TF1* gaus = new TF1(Form("gaus_%d", i), "gaus", 4000, 12000);
        gaus->SetParameters(multiGaus->GetParameter(4 + 3 * i), multiGaus->GetParameter(5 + 3 * i), multiGaus->GetParameter(6 + 3 * i));
        gaus->SetLineColor(kBlue);
        gaus->SetNpx(2000); // Increase the number of points to make it smooth
        gaus->Draw("SAME");
    }

    // Define the tail function for each peak separately to draw them
    for (int j = 0; j < 4; ++j) {
        TF1* tail = new TF1(Form("tail_%d", j), [initialParams, j](double* x, double* p) {
            double amp = initialParams[3 * j];
            double mean = initialParams[3 * j + 1];
            double sigma = initialParams[3 * j + 2];
            double tailParam1 = p[0];
            double tailParam2 = p[1];
            double val = tailParam1 * amp * exp((x[0] - mean) / (tailParam2 * sigma) + 1 / (2 * tailParam2 * tailParam2)) *
                TMath::Erfc((x[0] - mean) / (sqrt(2) * sigma) + 1 / (sqrt(2) * tailParam2));
            return val > 0 ? val : 1e-3; // Ensure the function value is always positive
            }, 4000, 12000, 2);

        tail->SetParameters(multiGaus->GetParameter(28 + j * 2), multiGaus->GetParameter(29 + j * 2));
        tail->SetLineColor(kGreen);
        tail->SetNpx(2000); // Increase the number of points to make it smooth
        tail->Draw("SAME");
    }

    // Define the background function separately to draw it
    TF1* background = new TF1("background", [](double* x, double* p) {
        double val = p[0] * exp(-p[1] * x[0]) + p[2] * x[0] + p[3];
        return val > 0 ? val : 1e-3; // Ensure the function value is always positive
        }, 4000, 12000, 4);
    background->SetParameters(multiGaus->GetParameter(0), multiGaus->GetParameter(1), multiGaus->GetParameter(2), multiGaus->GetParameter(3));
    background->SetLineColor(kMagenta); // Set the color to purple
    background->SetNpx(2000); // Increase the number of points to make it smooth
    background->Draw("SAME");

    // Create a TPaveText to display chi-squared, chi-squared/NDF, and p-value
    double chi2 = multiGaus->GetChisquare();
    int ndf = multiGaus->GetNDF();
    double chi2_ndf = chi2 / ndf;
    double pValue = TMath::Prob(chi2, ndf); // Calculate the p-value

    // Print chi-squared
    std::cout << "Chi-squared = " << chi2 << std::endl;


    TPaveText* pt = new TPaveText(0.65, 0.6, 0.85, 0.75, "NDC");
    pt->SetFillColor(0);
    pt->AddText(("Chi-squared = " + std::to_string(chi2)).c_str());
    pt->AddText(("Chi-squared / NDF = " + std::to_string(chi2_ndf)).c_str());
    pt->AddText(("p-value = " + std::to_string(pValue)).c_str());
    pt->Draw();

    // Print the final fit parameters with errors
    std::cout << "Global fit parameters:" << std::endl;
    for (int i = 0; i < 8; ++i) {
        std::cout << "Gaussian " << i + 1 << ": " << std::endl;
        std::cout << "  Amplitude = " << multiGaus->GetParameter(4 + 3 * i) << " +/- " << multiGaus->GetParError(4 + 3 * i) << std::endl;
        std::cout << "  Mean = " << multiGaus->GetParameter(5 + 3 * i) << " +/- " << multiGaus->GetParError(5 + 3 * i) << std::endl;
        std::cout << "  Sigma = " << multiGaus->GetParameter(6 + 3 * i) << " +/- " << multiGaus->GetParError(6 + 3 * i) << std::endl;
    }

    std::cout << "Tail function parameters:" << std::endl;
    for (int j = 0; j < 4; ++j) {
        std::cout << "Tail " << j + 1 << ": " << std::endl;
        std::cout << "  Tail amplitude factor = " << multiGaus->GetParameter(28 + j * 2) << " +/- " << multiGaus->GetParError(28 + j * 2) << std::endl;
        std::cout << "  Tail parameter = " << multiGaus->GetParameter(29 + j * 2) << " +/- " << multiGaus->GetParError(29 + j * 2) << std::endl;
    }

    std::cout << "Background function parameters:" << std::endl;
    std::cout << "  A = " << multiGaus->GetParameter(0) << " +/- " << multiGaus->GetParError(0) << std::endl;
    std::cout << "  B = " << multiGaus->GetParameter(1) << " +/- " << multiGaus->GetParError(1) << std::endl;
    std::cout << "  C = " << multiGaus->GetParameter(2) << " +/- " << multiGaus->GetParError(2) << std::endl;
    std::cout << "  D = " << multiGaus->GetParameter(3) << " +/- " << multiGaus->GetParError(3) << std::endl;
    std::cout << "Chi-squared / NDF = " << chi2_ndf << std::endl;
    std::cout << "p-value = " << pValue << std::endl;

    // Calculate the sum of counts in the range 7500 to 7900
    int binMin = hist->FindBin(7500);
    int binMax = hist->FindBin(7900);
    double sumCounts = 0;

    for (int bin = binMin; bin <= binMax; ++bin) {
        sumCounts += hist->GetBinContent(bin);
    }

    std::cout << "Number of counts in the range 7500 to 7900: " << sumCounts << std::endl;

    // Update the canvas to display the histogram and fits
    c1->Modified();
    c1->Update();
    gPad->Update(); // Ensure the pad is updated

    // Process ROOT events to keep the canvas open
    std::cout << "Press Ctrl+C to exit the application." << std::endl;
    while (!gROOT->IsBatch()) {
        gSystem->ProcessEvents();
        c1->Update(); // Ensure the canvas is updated continuously
        gSystem->Sleep(100);
    }

    // Close the file
    file->Close();
}

void FinalFitCode() {
    fitGaussianPeaks();
}