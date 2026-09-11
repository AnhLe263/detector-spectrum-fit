// ============================================================================
//  AnalyzeDetectorSpectrum_Calibrated.C
// ----------------------------------------------------------------------------
//  Usage   : 1) root -l AnalyzeDetectorSpectrum_Calibrated.C 
//            2) root -l 'AnalyzeDetectorSpectrum_Calibrated.C(true)'
//  The second usage will display the spectrum only, without fitting or area calculation.
//  Author  : Le Tuan Anh; See more: https://github.com/AnhLe263/detector-spectrum-fit.git 
//  Date    : 9/09/2026
// ============================================================================

#include "TCanvas.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <utility>
#include <tuple>

// ============================================================================
//  Forward declarations of functions
// ============================================================================
enum class ParamType { kAmp = 0, kMean = 1, kSigma = 2 }; // for fixed peaks
struct FixedParam {
    int       peakIndex;   // 1-based peak number (1 = first peak, 2 = second, ...)
    ParamType param;       // which parameter of that peak to fix
    double    value;       // the value to fix it at
};

void FitBackgroundSides(TH1* h,
                         double lowWinMin, double lowWinMax,
                         double highWinMin, double highWinMax,
                         bool useExpBkg,
                         double &p0, double &p1);

std::tuple<TF1*, TF1*, std::vector<TF1*>> FitGaussPeaks(TH1* h, int nPeaks,
                    double xmin, double xmax,
                    const std::vector<double>& peakGuess,
                    const std::vector<double>& sigmaGuess,
                    bool useExpBkg,
                    double lowWinMin, double lowWinMax,      
                    double highWinMin, double highWinMax,const std::vector<FixedParam>& fixedParams = {});


// ============================================================================
//  Main macro
// ============================================================================
void AnalyzeDetectorSpectrum_Calibrated(bool displayOnly = false)
{
    // ------------------------------------------------------------------
    // 0. Energy calibration: channel -> keV,  E = a*ch + b
    // Changes values of 4096, E0 and ch0 to match the new calibration points
    // ------------------------------------------------------------------
    TF1 *fEcal = new TF1("fEcal", "[0]*x+[1]", 0, 4096);
    double kE[4]  = {0, 3157, 5156.59, 5485};       // known energies (keV)
    double ch[4] = {80.50, 1008, 1592.2, 1691};    // corresponding channels
    TGraph *grCal = new TGraph(4, ch, kE);
    grCal->Fit(fEcal, "Q");   // Q = quiet, no printout

    std::cout << "Calibration: E(keV) = " << fEcal->GetParameter(0)
            << " * ch + " << fEcal->GetParameter(1) << std::endl;
    // ------------------------------------------------------------------
    // 1. User settings
    // ------------------------------------------------------------------
    const char* INPUT_FILE = "Histo_test.txt";  // path to the spectrum file
    std::vector<FixedParam> FIXED_PARAMS; //for fixed parameters, if any. 
    // Region of interest (ROI) to zoom in on, in channel units.
    // Edit these two values to look at a different part of the spectrum.
    double ROI_MIN = 2519.43;  
    double ROI_MAX = 3337.43; 
    // --- Peak-fitting settings for the ROI ---
    const int    N_PEAKS     = 5;
    const bool   USE_EXP_BKG = false; //Dùng Bkg exp thì bật true
    std::vector<double> PEAK_GUESS  = {2720.52, 2830, 2928.43, 3023.86, 3163.60};
    std::vector<double> SIGMA_GUESS = {40.90, 17.04, 20.45, 17.04, 11.93};
    // Peak-free side windows used to estimate the background
    const double BKG_LOW_MIN  = 2519.43, BKG_LOW_MAX  = 2638.72;   // just before peak 1
    const double BKG_HIGH_MIN = 3252.22, BKG_HIGH_MAX = 3337.43; // just after last peak 
    
    // Nếu muôn cố định các thông số đỉnh, hãy điền vector FixedParam :
    //FIXED_PARAMS = {
    //    { 3, ParamType::kMean, 2931 }   // fix peak #3's Mean = 2931
    //    ,{ 4, ParamType::kSigma, 15 }   // fix peak #4's Sigma = 15
    //};
    // ------------------------------------------------------------------
    // 2. Read the spectrum: one count value per line, channel = line index
    // ------------------------------------------------------------------
    std::ifstream inFile(INPUT_FILE);
    if (!inFile.is_open()) {
        std::cerr << "ERROR: cannot open input file: " << INPUT_FILE << std::endl;
        return;
    }

    std::vector<double> counts;
    double value;
    while (inFile >> value) {
        counts.push_back(value);
    }
    inFile.close();

    const int nChannels = static_cast<int>(counts.size());
    if (nChannels == 0) {
        std::cerr << "ERROR: no data read from " << INPUT_FILE << std::endl;
        return;
    }
    std::cout << "Read " << nChannels << " channels from " << INPUT_FILE << std::endl;
    // ------------------------------------------------------------------
    // 3. Build the histogram 
    // ------------------------------------------------------------------
    TH1F* hSpectrum = new TH1F("hSpectrum", "Detector spectrum;Energy (keV);Counts",
                                nChannels, fEcal->Eval(-0.5), fEcal->Eval(nChannels - 0.5));//-0.5 cho dich bin center: tâm bin giờ trùng khít với năng lượng calib đúng của kênh i
    for (int i = 0; i < nChannels; ++i) {
        hSpectrum->SetBinContent(i + 1, counts[i]);
    }
    hSpectrum->SetLineColor(kAzure + 2);
    hSpectrum->SetLineWidth(4);
    if (displayOnly) {
        ROI_MIN = hSpectrum->GetXaxis()->GetXmin();
        ROI_MAX = hSpectrum->GetXaxis()->GetXmax();
    }
    // Clone the histogram for the zoomed pad so the two pads keep
    // independent axis ranges.
    TH1F* hZoom = (TH1F*)hSpectrum->Clone("hZoom");
    hZoom->SetTitle(Form("Zoomed spectrum (ROI %.0f - %.0f keV);Energy (keV);Counts",
                      ROI_MIN, ROI_MAX));
    //hZoom->SetLineColor(kBlack + 1);
    hZoom->SetLineWidth(4);
    hZoom->SetMarkerStyle(20); hZoom->SetMarkerColor(kAzure + 2);
    hZoom->GetXaxis()->SetRangeUser(ROI_MIN, ROI_MAX);

    // ------------------------------------------------------------------
    // 4. Draw on a 2-pad canvas
    // ------------------------------------------------------------------
    gStyle->SetOptStat(0);   // Ẩn stats box

    TCanvas* c1 = new TCanvas("c1", "Detector spectrum", 1400, 600);
    c1->Divide(2, 1);

    // --- Pad 1: full spectrum, log Y ---
    c1->cd(1);
    gPad->SetLogy();
    //gPad->SetGrid();
    hSpectrum->Draw("HIST");

    // Mark the ROI on the full spectrum with two vertical dashed lines
    double yMax = hSpectrum->GetMaximum();
    TLine* lLow  = new TLine(ROI_MIN, 1, ROI_MIN, yMax);
    TLine* lHigh = new TLine(ROI_MAX, 1, ROI_MAX, yMax);
    lLow->SetLineColor(kRed + 1);
    lHigh->SetLineColor(kRed + 1);
    lLow->SetLineStyle(2);
    lHigh->SetLineStyle(2);
    lLow->Draw("SAME");
    lHigh->Draw("SAME");

    // --- Pad 2: zoomed region of interest, log Y ---
    c1->cd(2);
    gPad->SetLogy();
    // gPad->SetGrid();
    if (!displayOnly) {
        hZoom->Draw("P");
    } else {
        hZoom->Draw("HIST");
    }
    TLegend* leg = new TLegend(0.20, 0.7, 0.58, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hZoom,"Exp. Data","p");
    if (!displayOnly) {
        auto fitFuncs  = FitGaussPeaks(hZoom, N_PEAKS, ROI_MIN, ROI_MAX,
                               PEAK_GUESS, SIGMA_GUESS, USE_EXP_BKG,
                               BKG_LOW_MIN, BKG_LOW_MAX,
                               BKG_HIGH_MIN, BKG_HIGH_MAX
                               ,FIXED_PARAMS
                            );
        auto totalFit         = std::get<0>(fitFuncs);
        auto retrievedBkgFunc  = std::get<1>(fitFuncs);
        auto peakFuncList      = std::get<2>(fitFuncs);
        //if (totalFit) leg->AddEntry(totalFit, "Total fit","l");
        if (!peakFuncList.empty()) leg->AddEntry(peakFuncList[0], "Fitted peaks", "l");
        if (retrievedBkgFunc) leg->AddEntry(retrievedBkgFunc, "background","l");
    } else {
        std::cout << "\n[displayOnly = true] Skipping fit and area calculation.\n" << std::endl;
    }

    leg->Draw();

    c1->cd();
    c1->Update();

    c1->SaveAs("detector_spectrum.png");
}

std::tuple<TF1*, TF1*, std::vector<TF1*>> FitGaussPeaks(TH1* h, int nPeaks,
                    double xmin, double xmax,
                    const std::vector<double>& peakGuess,
                    const std::vector<double>& sigmaGuess,
                    bool useExpBkg,
                    double lowWinMin, double lowWinMax,      
                    double highWinMin, double highWinMax,const std::vector<FixedParam>& fixedParams)    
{
    if ((int)peakGuess.size() != nPeaks || (int)sigmaGuess.size() != nPeaks) {
        std::cerr << "ERROR: peakGuess/sigmaGuess size must match nPeaks" << std::endl;
        return std::make_tuple(nullptr, nullptr, std::vector<TF1*>());
    }

    // ------------------------------------------------------------
    // 1. Build the fit formula: sum of nPeaks Gaussian terms + background
    // ------------------------------------------------------------
    TString formula = "";
    for (int k = 0; k < nPeaks; ++k) {
        if (k > 0) formula += "+";
        formula += Form("gaus(%d)", 3 * k);
    }
    int bkgStart = 3 * nPeaks; // 3: 3 thông số của Gauss
    formula += useExpBkg ? Form("+expo(%d)", bkgStart)
                          : Form("+pol1(%d)", bkgStart);

    TF1* fitFunc = new TF1("fitFunc", formula, xmin, xmax);
    fitFunc->SetNpx(2000);
    fitFunc->SetLineColor(kBlack);
    fitFunc->SetLineWidth(2);

    // ------------------------------------------------------------
    // 2. Fit the background FIRST on the two peak-free side windows,
    //    
    // ------------------------------------------------------------
    double bkgP0, bkgP1;
    FitBackgroundSides(h, lowWinMin, lowWinMax, highWinMin, highWinMax,
                        useExpBkg, bkgP0, bkgP1);

    fitFunc->SetParameter(bkgStart,     bkgP0);
    fitFunc->SetParameter(bkgStart + 1, bkgP1);
    fitFunc->FixParameter(bkgStart,     bkgP0);   // background is now FIXED
    fitFunc->FixParameter(bkgStart + 1, bkgP1);   // not re-fitted with the peaks
    fitFunc->SetParName(bkgStart,     "Bkg_p0 (fixed)");
    fitFunc->SetParName(bkgStart + 1, "Bkg_p1 (fixed)");

    TF1* bkgFunc = new TF1("bkgFunc", useExpBkg ? "expo(0)" : "pol1(0)", xmin, xmax);
    bkgFunc->SetParameters(bkgP0, bkgP1);
    bkgFunc->SetLineColor(1);
    bkgFunc->SetLineStyle(1);
    bkgFunc->SetLineWidth(4);
    // ------------------------------------------------------------
    // 3. Set initial parameters and limits for each peak
    // ------------------------------------------------------------
    // Retrieve the bin indices first, then get their coordinate values
    int firstBin = h->GetXaxis()->GetFirst();
    int lastBin  = h->GetXaxis()->GetLast();
    double userMin = h->GetXaxis()->GetBinLowEdge(firstBin);
    double userMax = h->GetXaxis()->GetBinUpEdge(lastBin);
    double ampGuess = h->GetMaximum();
    for (int k = 0; k < nPeaks; ++k) {
        int ip = 3 * k;
        h->GetXaxis()->SetRangeUser(peakGuess[k] - 2 * sigmaGuess[k], peakGuess[k] + 2 * sigmaGuess[k]);
        ampGuess = h->GetMaximum() - bkgFunc->Eval(peakGuess[k]); // subtract background at peak position
        fitFunc->SetParameter(ip,     ampGuess);
        fitFunc->SetParameter(ip + 1, peakGuess[k]);
        fitFunc->SetParameter(ip + 2, sigmaGuess[k]);

        fitFunc->SetParLimits(ip,     0, 10 * ampGuess);
        fitFunc->SetParLimits(ip + 1, peakGuess[k] - 5 * sigmaGuess[k],
                                        peakGuess[k] + 5 * sigmaGuess[k]); // keep mean near guess, nên để
        fitFunc->SetParLimits(ip + 2, 1e-3, 10 * sigmaGuess[k]);

        fitFunc->SetParName(ip,     Form("Amp_%d",   k + 1));
        fitFunc->SetParName(ip + 1, Form("Mean_%d",  k + 1));
        fitFunc->SetParName(ip + 2, Form("Sigma_%d", k + 1));
    }
    h->GetXaxis()->SetRangeUser(userMin, userMax);// restore the original axis range after retrieving the bin edges
    // ------------------------------------------------------------
    // Apply any user-requested fixed parameters, overriding the
    //    Chú ý: free/limited setting from step 2 above for that one parameter.
    // ------------------------------------------------------------
    for (const auto& fp : fixedParams) {
        if (fp.peakIndex < 1 || fp.peakIndex > nPeaks) {
            std::cerr << "WARNING: FixedParam peakIndex " << fp.peakIndex
                       << " is out of range [1," << nPeaks << "], skipped." << std::endl;
            continue;
        }
        int ip  = 3 * (fp.peakIndex - 1);
        int idx = ip + static_cast<int>(fp.param);
        fitFunc->FixParameter(idx, fp.value);
        std::cout << "Fixed " << fitFunc->GetParName(idx)
                   << " = " << fp.value << std::endl;
    }
    

    // ------------------------------------------------------------
    // 4. Perform the fit (only the peak parameters are free now)
    // ------------------------------------------------------------
    TFitResultPtr fitResult = h->Fit(fitFunc, "RSM+NO");
    fitResult->Print();
    std::cout << "Chi2/NDF = " << fitFunc->GetChisquare() << " / " << fitFunc->GetNDF()
               << " = " << fitFunc->GetChisquare() / fitFunc->GetNDF() << std::endl;

    // ------------------------------------------------------------
    // 5. Draw each individual peak (on top of the fixed background)
    // ------------------------------------------------------------
    std::vector<TF1*> peakFuncs; 
    int colors[] = {kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet};
    for (int k = 0; k < nPeaks; ++k) {
        TF1* peakFunc = new TF1(Form("peak_%d", k + 1),
                                 useExpBkg ? "gaus(0)+expo(3)" : "gaus(0)+pol1(3)",
                                 xmin, xmax);
        peakFunc->SetParameters(fitFunc->GetParameter(3 * k),
                                 fitFunc->GetParameter(3 * k + 1),
                                 fitFunc->GetParameter(3 * k + 2),
                                 fitFunc->GetParameter(bkgStart),
                                 fitFunc->GetParameter(bkgStart + 1));
        peakFunc->SetLineColor(6);
        peakFunc->SetLineStyle(10);peakFunc->SetLineWidth(4);
        peakFunc->SetNpx(2000);
        peakFunc->Draw("SAME");
	peakFuncs.push_back(peakFunc);
    }

    bkgFunc->Draw("SAME");

    // ------------------------------------------------------------
    // 6. Compute the area (integral) of each Gaussian peak, with its
    //    uncertainty propagated from the fit's covariance matrix.
    //    For a Gaussian: Area = Amp * Sigma * sqrt(2*pi)/binWidth
    //    Amp and Sigma are correlated (both come from the same fit),
    //    so the correlation term MUST be included, not just the two
    //    individual errors added in quadrature.
    // ------------------------------------------------------------
    std::cout << "\n--- Peak areas (from fitted Gaussian parameters) ---" << std::endl;

    double binW = h->GetBinWidth(1); 
    for (int k = 0; k < nPeaks; ++k) {
        int ip = 3 * k;
        double amp      = fitFunc->GetParameter(ip);
        double mean     = fitFunc->GetParameter(ip + 1);
        double sigma    = fitFunc->GetParameter(ip + 2);
        double errAmp   = fitFunc->GetParError(ip);
        double errSigma = fitFunc->GetParError(ip + 2);
        double corr     = fitResult->Correlation(ip, ip + 2);

        const double sqrt2pi = std::sqrt(2.0 * TMath::Pi());
        double area = amp * sigma * sqrt2pi/binW;  

        // Error propagation: dA/dAmp = sigma*sqrt(2pi), dA/dSigma = amp*sqrt(2pi)
        double dA_damp   = sigma * sqrt2pi/binW;
        double dA_dsigma = amp   * sqrt2pi/binW;
        double covAmpSigma = corr * errAmp * errSigma;

        double areaErr2 = dA_damp * dA_damp * errAmp * errAmp
                         + dA_dsigma * dA_dsigma * errSigma * errSigma
                         + 2.0 * dA_damp * dA_dsigma * covAmpSigma;
        double areaErr = (areaErr2 > 0) ? std::sqrt(areaErr2) : 0.0;

        std::cout << "  Peak " << k + 1
                   << " : Mean = " << mean
                   << ",  Area = " << area
                   << " +/- " << areaErr
                   << "  (" << (area > 0 ? 100.0 * areaErr / area : 0.0) << " %)"
                   << std::endl;
    }
    
    return std::make_tuple(fitFunc, bkgFunc, peakFuncs);
}


// ============================================================================
//  FitBackgroundSides
// ----------------------------------------------------------------------------
//  Estimate the background (linear or exponential) using ONLY two peak-free
//  "side windows", one just below and one just above the peak region. Using
//  a TGraphErrors (not the histogram directly) avoids feeding zero-content
//  bins from the peak region into the fit, so the background is not biased
//  by the peaks themselves.
//
//    h                       : histogram containing the full spectrum
//    lowWinMin,  lowWinMax   : side window below the peaks (background only)
//    highWinMin, highWinMax  : side window above the peaks (background only)
//    useExpBkg               : true = exponential, false = linear
//    p0, p1 (output)         : fitted background parameters
// ============================================================================
void FitBackgroundSides(TH1* h,
                         double lowWinMin, double lowWinMax,
                         double highWinMin, double highWinMax,
                         bool useExpBkg,
                         double &p0, double &p1)
{
    // Collect only the bins inside the two peak-free side windows.
    // The empty region between them (where the peaks sit) is never seen
    // by the fit, so it cannot bias the background estimate.
    TGraphErrors* gBkg = new TGraphErrors();
    int nPts = 0;
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        double x = h->GetBinCenter(b);
        bool inLow  = (x >= lowWinMin  && x <= lowWinMax);
        bool inHigh = (x >= highWinMin && x <= highWinMax);
        if (!inLow && !inHigh) continue;

        double y  = h->GetBinContent(b);
        double ey = (y > 0) ? std::sqrt(y) : 1.0;   // Poisson error
        gBkg->SetPoint(nPts, x, y);
        gBkg->SetPointError(nPts, 0, ey);
        ++nPts;
    }

    if (nPts < 2) {
        std::cerr << "ERROR: not enough points in side windows to fit background" << std::endl;
        p0 = 0; p1 = 0;
        delete gBkg;
        return;
    }

    TF1* bkgOnly = new TF1("bkgOnly", useExpBkg ? "expo" : "pol1", lowWinMin, highWinMax);
    gBkg->Fit(bkgOnly, "RQ");   // R = respect range, Q = quiet (no printout)

    p0 = bkgOnly->GetParameter(0);
    p1 = bkgOnly->GetParameter(1);

    std::cout << "Background fit from side windows [" << lowWinMin << "," << lowWinMax
               << "] + [" << highWinMin << "," << highWinMax << "]:  "
               << "p0 = " << p0 << ",  p1 = " << p1 << std::endl;

    delete gBkg;
    delete bkgOnly;
}

