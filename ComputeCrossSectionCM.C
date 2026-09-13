// ============================================================================
//  ComputeCrossSectionCM.C
// ----------------------------------------------------------------------------
//  Purpose : Read yields measured at several Lab angles
//            for a series of incident particle energies, convert each yield
//            into a differential cross section dsigma/dOmega in the
//            Center-of-Mass (CM) frame, and:
//              - write a results table to cross_section_cm.txt
//              - plot dsigma/dOmega(CM) vs Ep for each angle (excitation
//                functions), as an inline ROOT canvas
//
//  Physics : theta_Lab -> theta_CM via the exact 2-body kinematic relation
//            theta_CM = theta_Lab + asin(gamma*sin(theta_Lab)), then the
//            Lab->CM solid-angle Jacobian J = dOmega_Lab/dOmega_CM converts
//            the Lab-frame cross section to the CM frame. Dead-time
//            correction and statistical + systematic errors are propagated
//            into the final uncertainty. The statistical error uses the
//            peak-fit uncertainty (dY) when available, falling back to
//            plain Poisson (1/sqrt(Y)) only when dY is missing ("x").
//
//  Input   : cross_section_input.txt -- see the comment header inside that file for the
//            full column layout. In short: lines starting with '#' are
//            skipped; the first non-comment line is the column header,
//            from which the Lab angle list is read dynamically (one Y/dY
//            column pair per angle); every following line is one data
//            point (Np, RunId, Ep, Nt, then a Y/dY pair per angle).
//
//  Usage   : root -l ComputeCrossSectionCM.C
//
//  Author  : Le Tuan Anh
//  Date    : 13/09/2026
// ============================================================================

#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>

#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

// ============================================================================
//  Physics constants for this specific reaction/setup
//  (elastic scattering: m1 = m3 = projectile, m2 = m4 = target, Q = 0)
// ============================================================================
const Double_t SOLID_ANGLE_SR = 0.01227;   // W: detector solid angle (sr), same for all detectors
const Double_t DEAD_TIME      = 0.05;      // 5% electronics dead time
// Systematic error budget (added in quadrature with the statistical error):
//   5% target thickness + 5% detector efficiency/solid angle + 2% beam charge
const Double_t SYS_ERR_SQ     = 0.05 * 0.05 + 0.05 * 0.05 + 0.02 * 0.02;

const Double_t M_PROTON = 1.007825;   // amu
const Double_t M_TARGET = 18.998403;  // amu (19F for example)
const Double_t Q_VALUE   = 0.0;        // MeV, elastic scattering

const Int_t    N_FIXED_COLUMNS = 4;    // Np, RunId, Ep, Nt, before the Y/dY pairs

// ============================================================================
//  Data structures and function declarations
// ============================================================================
struct CrossSectionResult {
    Double_t theta_CM_deg;
    Double_t dsigma_L;    // cm^2/sr
    Double_t dsigma_CM;   // cm^2/sr
    Double_t err_CM;      // cm^2/sr
};
CrossSectionResult ComputeDiffCrossSectionCM(Double_t yieldRaw, Double_t Ep_MeV, Double_t angleDeg,
                                              Double_t Np, Double_t Nt, Double_t errY = -1.0);
Bool_t ParseYield(const TString& tok, Double_t& val);
Int_t ParseHeaderAngles(const TString& headerLine, std::vector<Double_t>& labAngles);
Bool_t ParseDataLine(const TString& rawLine, Int_t nAngles,
                     Double_t& Np, Double_t& Ep, Double_t& Nt,
                     std::vector<TString>& yTokens, std::vector<TString>& dyTokens);

// ============================================================================
//  Main macro
// ============================================================================
void ComputeCrossSectionCM()
{
    const char* INPUT_FILE  = "cross_section_input.txt";
    const char* OUTPUT_FILE = "cross_section_cm.txt";

    std::ifstream inFile(INPUT_FILE);
    if (!inFile.is_open()) {
        std::cerr << "ERROR: cannot open input file: " << INPUT_FILE << std::endl;
        return;
    }

    // ------------------------------------------------------------------
    // 1. Skip comment lines ('#'), then read the header line to get the
    //    Lab angle list dynamically.
    // ------------------------------------------------------------------
    std::vector<Double_t> labAngles;
    Int_t nAngles = -1;
    std::string lineStd;

    while (std::getline(inFile, lineStd)) {
        TString line(lineStd.c_str());
        line.ReplaceAll("\r", "");
        if (line.IsWhitespace() || line.BeginsWith("#")) continue;   // skip comments/blank lines

        nAngles = ParseHeaderAngles(line, labAngles);
        break;   // this was the header line; data rows follow
    }

    if (nAngles <= 0) {
        std::cerr << "ERROR: could not read a valid header line from " << INPUT_FILE << std::endl;
        return;
    }

    std::cout << "Found " << nAngles << " angle(s): ";
    for (Double_t a : labAngles) std::cout << a << " ";
    std::cout << "(deg)" << std::endl;

    // ------------------------------------------------------------------
    // 2. Prepare the output table and per-angle TGraphErrors
    // ------------------------------------------------------------------
    std::ofstream outFile(OUTPUT_FILE);
    outFile << std::left << std::setw(10) << "Ep(MeV)" << std::setw(15) << "Np";
    for (Int_t i = 0; i < nAngles; ++i) {
        outFile << std::setw(15) << Form("CS_CM_%.0f", labAngles[i])
                 << std::setw(15) << Form("Err_%.0f",   labAngles[i]);
    }
    outFile << "\n" << std::string(15 * 2 * nAngles + 25, '-') << "\n";

    std::vector<TGraphErrors*> graphs(nAngles);
    for (Int_t i = 0; i < nAngles; ++i) {
        graphs[i] = new TGraphErrors();
        //graphs[i]->SetName(Form("gCS_%.0fdeg", labAngles[i]));
        graphs[i]->SetName(Form("%.0f#circ", labAngles[i]));
        //graphs[i]->SetTitle(Form("%.0f#circ", labAngles[i]));
    }

    // ------------------------------------------------------------------
    // 3. Read and process every remaining data line
    // ------------------------------------------------------------------
    while (std::getline(inFile, lineStd)) {
        TString line(lineStd.c_str());
        line.ReplaceAll("\r", "");
        if (line.IsWhitespace() || line.BeginsWith("#")) continue;   // skip comments/blank lines

        Double_t Np, Ep, Nt;
        std::vector<TString> yTokens, dyTokens;
        if (!ParseDataLine(line, nAngles, Np, Ep, Nt, yTokens, dyTokens)) continue;

        outFile << std::left << std::setw(10) << std::fixed << std::setprecision(3) << Ep
                 << std::setw(15) << std::scientific << std::setprecision(5) << Np;

        for (Int_t i = 0; i < nAngles; ++i) {
            Double_t yieldRaw, errY;
            Bool_t hasY  = ParseYield(yTokens[i], yieldRaw) && Np > 0 && Nt > 0 && yieldRaw > 0;
            Bool_t hasDY = hasY && ParseYield(dyTokens[i], errY);   // dY only meaningful if Y is valid

            if (!hasY) {
                outFile << std::left << std::setw(15) << "NaN" << std::setw(15) << "NaN";
                continue;
            }

            CrossSectionResult xs = ComputeDiffCrossSectionCM(
                yieldRaw, Ep, labAngles[i], Np, Nt, hasDY ? errY : -1.0);

            outFile << std::left << std::setw(15) << std::scientific << std::setprecision(5) << xs.dsigma_CM
                     << std::setw(15) << std::scientific << std::setprecision(5) << xs.err_CM;

            Int_t pt = graphs[i]->GetN();
            graphs[i]->SetPoint(pt, Ep, xs.dsigma_CM);
            graphs[i]->SetPointError(pt, 0.0, xs.err_CM);
        }
        outFile << "\n";
    }
    inFile.close();
    outFile.close();
    std::cout << "Done. Results written to " << OUTPUT_FILE << std::endl;

    //Done main part!!!
    // ------------------------------------------------------------------
    // 4. Plot the excitation functions dsigma/dOmega(CM) vs Ep
    // ------------------------------------------------------------------
    // Enable ticks on the opposite sides (top and right)
    gStyle->SetPadTickX(1); 
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(2);
    
    double global_xmin = 1e9;
    double global_xmax = -1e9;
    for (int i = 0; i < nAngles; i++) {
        double xmin, xmax, ymin, ymax;
        graphs[i]->ComputeRange(xmin, ymin, xmax, ymax); // Lấy range tự nhiên của đồ thị
        if (xmin < global_xmin) global_xmin = xmin;
        if (xmax > global_xmax) global_xmax = xmax;
    }
    // Thêm một chút khoảng thở (margin) 5% ở 2 đầu nếu muốn dữ liệu không chạm vạch
    double x_margin = (global_xmax - global_xmin) * 0.05;
    global_xmin -= x_margin;
    global_xmax += x_margin;
    
    // Khởi tạo Canvas (Tự động tăng chiều cao nếu có nhiều Pad để tránh bị nén)
    int canvas_width = 800;
    int canvas_height = 300 + (nAngles * 150); // Tự dãn chiều cao theo số lượng Pad
    TCanvas *c = new TCanvas("c", "Dynamic Pads", canvas_width, canvas_height);
    // padA bên trái chiếm 15% chiều rộng (0.0 đến 0.15) // padB bên phải chiếm 85% chiều rộng còn lại (0.15 đến 1.0)
    TPad *padA = new TPad("padA", "Pad chứa Latex", 0.0, 0.0, 0.05, 1.0);
    TPad *padB = new TPad("padB", "Pad chứa các pad con", 0.05, 0.0, 1.0, 1.0);
    padA->SetMargin(0,0,0,0);padB->SetMargin(0,0,0,0);
    padA->Draw();
    padB->Draw();
    
    //Vẽ và thêm nhãn Latex vào padA
    padA->cd();
    TLatex *tex = new TLatex(0.5, 0.5, "d#sigma/d#Omega_{CM} (cm^{2}/sr)");
    tex->SetTextAlign(22); // Căn giữa chữ
    tex->SetTextSize(0.5);
    tex->SetTextAngle(90);  // Xoay 90 độ vì pad này rất hẹp theo chiều ngang
    tex->Draw();

    // Định nghĩa các khoảng lề cố định bằng hằng số
    const double top_margin = 0.05;    // Lề trên cùng của Canvas
    const double bottom_margin = 0.05; // Lề dưới cùng để chứa nhãn trục X
    const double pad_margin = 0.001;  // Khoảng hở nhỏ giữa các pad để không đè vạch

    // Tính toán chiều cao thực tế của vùng vẽ (trừ đi lề trên và lề dưới)
    double available_height = 1.0 - top_margin - bottom_margin;
    double pad_height = available_height / nAngles; // Chiều cao bằng nhau cho mỗi ô vẽ

    // Vòng lặp tự động tạo Pad và vẽ đồ thị
    Int_t colorPalette[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
                             kCyan + 2, kViolet, kSpring + 4, kAzure + 6, kPink + 6};
    for (int i = 0; i < nAngles; i++) {
        padB->cd();
        // Tính toán tọa độ ylow và yup cho từng pad (chạy từ dưới lên trên)
        // i = 0 là pad dưới cùng, i = nAngles-1 là pad trên cùng
        double ylow = bottom_margin + i * pad_height;
        double yup = ylow + pad_height;

        // Mở rộng tọa độ của pad đầu và pad cuối để chứa lề bên ngoài
        if (i == 0)   ylow = 0.0;
        if (i == nAngles-1) yup = 1.0;

        TString pad_name = Form("pad_%d", i);
        TPad *pad = new TPad(pad_name, pad_name, 0.0, ylow, 1.0, yup);
        pad->SetLogy();
        pad->SetGrid();
        // Thiết lập lề cho từng vị trí pad
        pad->SetLeftMargin(0.1);
        pad->SetRightMargin(0.05);
        TLatex latex;
        latex.SetNDC();
        latex.SetTextAlign(22);
        if (i == nAngles-1) { // Pad trên cùng
            pad->SetTopMargin(top_margin/(top_margin + pad_height)); 
            pad->SetBottomMargin(pad_margin);
            latex.SetTextSize(0.07);
        } 
        else if (i == 0) { // Pad dưới cùng
            pad->SetTopMargin(pad_margin);
            pad->SetBottomMargin(bottom_margin / (bottom_margin + pad_height));
            latex.SetTextSize(0.07);
        } 
        else { // Các pad ở giữa
            pad->SetTopMargin(pad_margin);
            pad->SetBottomMargin(pad_margin);
            latex.SetTextSize(0.09);
        }

        pad->Draw();
        pad->cd(); // Chuyển vào pad hiện tại để vẽ
        
        // ối ưu hiển thị trục X và cỡ chữ chữ theo số lượng Pad
        graphs[i]->GetXaxis()->SetLimits(global_xmin, global_xmax);
        double xmin, xmax, ymin, ymax;
        graphs[i]->ComputeRange(xmin, ymin, xmax, ymax);
        if(ymin <= 0) ymin = 1e-5; // LogY không thể chứa giá trị <= 0
        graphs[i]->GetYaxis()->SetRangeUser(ymin * 0.5, ymax * 2.0);
        //graphs[i]->GetYaxis()->SetTitle("d#sigma/d#Omega_{CM} (cm^{2}/sr)");
        graphs[i]->SetMarkerStyle(20);
        graphs[i]->SetMarkerColor(colorPalette[i % 10]);
        graphs[i]->SetLineColor(colorPalette[i % 10]);
        graphs[i]->GetXaxis()->SetTitleFont(43);
        graphs[i]->GetYaxis()->SetTitleSize(0.07);
        graphs[i]->GetYaxis()->SetTitleOffset(0.6);graphs[i]->GetYaxis()->CenterTitle();
        graphs[i]->GetYaxis()->SetLabelFont(43);
        graphs[i]->GetYaxis()->SetLabelSize(16);
        graphs[i]->GetXaxis()->SetTickLength(0.1);
        graphs[i]->GetYaxis()->SetTickLength(0.05);
        graphs[i]->Draw("AP");
        if (i == 0) {
            // Chỉ hiện nhãn trục X ở ô dưới cùng
            graphs[i]->GetXaxis()->SetTitle("E_{p} (MeV)");
            graphs[i]->GetXaxis()->SetTitleFont(62);
            graphs[i]->GetXaxis()->SetTitleSize(0.08);
            graphs[i]->GetXaxis()->CenterTitle(); 
            graphs[i]->GetXaxis()->SetLabelFont(43);
            graphs[i]->GetXaxis()->SetLabelSize(16);
        } else {
            // Ẩn hoàn toàn trục X ở các ô phía trên
            graphs[i]->GetXaxis()->SetLabelSize(0);
            graphs[i]->GetXaxis()->SetTitleSize(0);
        }
        
        if (i == nAngles-1) { // Pad trên cùng
            latex.DrawLatex(0.5, 0.70, graphs[i]->GetName());
        } 
        else { 
            latex.DrawLatex(0.5, 0.85, graphs[i]->GetName());
        } 
         
    }
    
    c->Update();
    c->SaveAs("cross_section_cm.png");
   
}



// ============================================================================
//  CrossSectionResult / ComputeDiffCrossSectionCM
// ----------------------------------------------------------------------------
//  Convert one measured yield into a Lab and CM differential cross section
//  using standard 2-body elastic-scattering kinematics.
//
//    yieldRaw : peak yield (counts)
//    errY     : uncertainty on yieldRaw, e.g. from AnalyzeDetectorSpectrum.C's
//               ComputePeakAreas. Pass a NEGATIVE value (default -1) to mean
//               "not available" -- the function then falls back to plain
//               Poisson statistics (1/sqrt(yieldRaw)) for the statistical
//               error term instead.
// ============================================================================


CrossSectionResult ComputeDiffCrossSectionCM(Double_t yieldRaw, Double_t Ep_MeV, Double_t angleDeg,
                                              Double_t Np, Double_t Nt, Double_t errY)
{
    Double_t thetaLab = angleDeg * TMath::DegToRad();

    // gamma = ratio of CM-frame velocity to particle velocity;
    Double_t gamma = TMath::Sqrt( (M_PROTON * M_PROTON) / (M_TARGET * M_TARGET)
                                 * (Ep_MeV / (Ep_MeV + Q_VALUE * (1.0 + M_PROTON / M_TARGET))) );

    Double_t thetaCM = thetaLab + TMath::ASin(gamma * TMath::Sin(thetaLab));

    // Lab -> CM solid-angle Jacobian: dOmega_Lab / dOmega_CM
    Double_t jacobian = (1.0 + gamma * TMath::Cos(thetaCM))
                       / TMath::Power(1.0 + gamma * gamma + 2.0 * gamma * TMath::Cos(thetaCM), 1.5);

    Double_t yieldCorrected = yieldRaw / (1.0 - DEAD_TIME);          // dead-time correction
    Double_t dsigma_L  = yieldCorrected / (Nt * Np * SOLID_ANGLE_SR);
    Double_t dsigma_CM = dsigma_L * jacobian;

    // Statistical error: use the peak-fit uncertainty (errY) when available,
    // otherwise fall back to plain Poisson counting statistics.
    Double_t relErrStat_sq;
    if (errY >= 0.0 && yieldRaw > 0.0) {
        relErrStat_sq = (errY / yieldRaw) * (errY / yieldRaw);   // from AnalyzeDetectorSpectrum.C's fit
    } else {
        relErrStat_sq = 1.0 / yieldRaw;                           // Poisson: (sigma_Y/Y)^2 = 1/Y
    }

    Double_t relErrTotal = TMath::Sqrt(relErrStat_sq + SYS_ERR_SQ);
    Double_t err_CM       = dsigma_CM * relErrTotal;

    return { thetaCM * TMath::RadToDeg(), dsigma_L, dsigma_CM, err_CM };
}


// ============================================================================
//  ParseYield
// ----------------------------------------------------------------------------
//  Returns kFALSE for missing-data markers ("x", "X", "-"); otherwise parses
//  the token into 'val' and returns kTRUE.
// ============================================================================
Bool_t ParseYield(const TString& tok, Double_t& val)
{
    if (tok == "x" || tok == "X" || tok == "-") return kFALSE;
    val = tok.Atof();
    return kTRUE;
}


// ============================================================================
//  ParseHeaderAngles
// ----------------------------------------------------------------------------
//  Read the Lab angle list directly from the column header line, e.g.:
//     Np   RunId   Ep   Nt   Y60   dY60   Y80   dY80   ...
//  Angles are taken from every "Y<angle>" label (every 2nd column after the
//  4 fixed leading columns); the matching "dY<angle>" label is expected
//  immediately after it, but is not itself parsed for its value here.
//
//  Returns the number of angles found, or -1 on a malformed header.
// ============================================================================
Int_t ParseHeaderAngles(const TString& headerLine, std::vector<Double_t>& labAngles)
{
    TObjArray* tokens = headerLine.Tokenize(" \t");
    Int_t nTok = tokens->GetEntriesFast();

    if (nTok <= N_FIXED_COLUMNS || (nTok - N_FIXED_COLUMNS) % 2 != 0) {
        std::cerr << "ERROR: header must have " << N_FIXED_COLUMNS
                   << " fixed columns followed by Y/dY pairs (got " << nTok << " columns)" << std::endl;
        delete tokens;
        return -1;
    }

    Int_t nAngles = (nTok - N_FIXED_COLUMNS) / 2;
    labAngles.clear();
    for (Int_t i = 0; i < nAngles; ++i) {
        TString yLabel = ((TObjString*)tokens->At(N_FIXED_COLUMNS + 2 * i))->GetString();
        if (!yLabel.BeginsWith("Y")) {
            std::cerr << "WARNING: expected a 'Y<angle>' label at column "
                       << N_FIXED_COLUMNS + 2 * i + 1 << ", got '" << yLabel << "'" << std::endl;
        }
        TString angleStr(yLabel);
        angleStr.Remove(0, 1);   // strip the leading 'Y'
        labAngles.push_back(angleStr.Atof());
    }

    delete tokens;
    return nAngles;
}


// ============================================================================
//  ParseDataLine
// ----------------------------------------------------------------------------
//  Extract Np, Ep, Nt and the per-angle (Y, dY) token pairs from one data
//  line, given the number of angles already determined from the header.
//  Returns kFALSE if the line doesn't have the expected number of columns.
// ============================================================================
Bool_t ParseDataLine(const TString& rawLine, Int_t nAngles,
                     Double_t& Np, Double_t& Ep, Double_t& Nt,
                     std::vector<TString>& yTokens, std::vector<TString>& dyTokens)
{
    TObjArray* tokens = rawLine.Tokenize(" \t,;");   // also accept comma/semicolon separators
    Int_t nTok = tokens->GetEntriesFast();
    Int_t nExpected = N_FIXED_COLUMNS + 2 * nAngles;

    if (nTok < nExpected) {
        delete tokens;
        return kFALSE;
    }

    Np  = ((TObjString*)tokens->At(0))->GetString().Atof();
    // tokens->At(1) is RunId, intentionally not used in the calculation
    Ep  = ((TObjString*)tokens->At(2))->GetString().Atof();
    Nt = ((TObjString*)tokens->At(3))->GetString().Atof();

    yTokens.clear();
    dyTokens.clear();
    for (Int_t i = 0; i < nAngles; ++i) {
        yTokens.push_back(((TObjString*)tokens->At(N_FIXED_COLUMNS + 2 * i))->GetString());
        dyTokens.push_back(((TObjString*)tokens->At(N_FIXED_COLUMNS + 2 * i + 1))->GetString());
    }

    delete tokens;
    return kTRUE;
}



