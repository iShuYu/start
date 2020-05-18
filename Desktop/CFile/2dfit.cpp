// Code to Fit the mass of JPsi
// Author : zhaomr
// Date   : June 4 2015
// Input  : dirData
// Output : dirMass     : fit parameters
//          dirMassPlot : plot
//          dirMassFile : root file

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <RooFit.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooHistPdf.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussModel.h>
#include <RooAddModel.h>
#include <RooAddPdf.h>
#include <RooDecay.h>
#include <RooCBShape.h>
#include <RooExponential.h>
#include <RooProdPdf.h>
#include <RooStats/SPlot.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;
using namespace RooFit;

int main(int argc, char* argv[])
{

	//------------------------File-------------------
//	TCanvas *can = new TCanvas("","",800,600);
	ifstream config(argv[1]);
	string input_file_config, tail_file_config, splot_file_config, root_file_config, output_pyfile_config, mass_fit_result, tzbkg_fit_result;
	config >> input_file_config >> tail_file_config >> splot_file_config >> root_file_config >> output_pyfile_config >> mass_fit_result >> tzbkg_fit_result;

	//double pt_low, pt_hi, y_low, y_hi;
	//config >> pt_low >> pt_hi >> y_low >> y_hi;

	// Create DecayTree with the data after reconstruction
	TFile input_file(input_file_config.c_str());
	TTree* data_tree = (TTree*)input_file.Get("DecayTree");

	TFile tail_file(tail_file_config.c_str());
	TTree* tail_tree = (TTree*)tail_file.Get("DecayTree");

	TFile* root_file = new TFile(root_file_config.c_str(), "recreate");

	// Useful Variables
	RooRealVar psi_mass("psi_M", "psi_M", 3096.9 - 120, 3096.9 + 120);
	RooRealVar tz("tz", "tz", -10, 10);
	RooRealVar tzError("tzError", "tzError", 0.001, 0.3);
	RooRealVar bkgSweight("bkg_sw", "bkg_sw", -100, 100);

	// Data Set to Save the Fit Variables
	RooDataSet data_set("data_set", "data set", data_tree, RooArgSet(psi_mass, tz, tzError));
	RooDataSet tail_set("tail_set", "tail set", tail_tree, RooArgSet(tz));
	RooDataHist* tail_hist = tail_set.binnedClone("");

	double N = data_set.numEntries();
	cout << "Total Entries = " << N << endl;
	cout << "Wrong PV Entries = " << tail_set.numEntries() << endl;
//	cout << pt_low << " " << pt_hi << endl;
//	cout << y_low << " " << y_hi << endl;

	//---------------------pdf-----------------------
	double fmass_mean, fmass_sigma1, fmass_sigma2, fmass_ratio, fmass_p0;
	ifstream fmass(mass_fit_result.c_str());
	fmass >> fmass_mean >> fmass_sigma1 >> fmass_sigma2 >> fmass_ratio >> fmass_p0;
	fmass.close();
	cout << fmass_mean << " " << fmass_sigma1 << " " << fmass_sigma2 << " " << fmass_ratio << " " << fmass_p0 << endl;
	RooRealVar mass_mean("mass_mean", "mass mean", fmass_mean);
	RooRealVar mass_sigma1("mass_sigma1", "mass sigma1", fmass_sigma1);
	RooRealVar mass_sigma2("mass_sigma2", "mass sigma2", fmass_sigma2);
	RooRealVar mass_ratio("mass_ratio", "mass ratio", fmass_ratio);
	RooRealVar mass_p0("mass_p0", "mass p0", fmass_p0);
	RooRealVar a("a", "a", -0.00011);
	RooRealVar b("b", "b", 0.0085);
	RooRealVar c("c", "c", 2.066);
	RooRealVar n("n", "n", 1.0);
	RooFormulaVar alpha1("alpha1", "alpha1", "a*mass_sigma1*mass_sigma1 + b*mass_sigma1 + c", RooArgSet(a, b, c, mass_sigma1));
	RooFormulaVar alpha2("alpha2", "alpha2", "a*mass_sigma2*mass_sigma2 + b*mass_sigma2 + c", RooArgSet(a, b, c, mass_sigma2));
	RooCBShape mass_myCB1("mass_myCB1", "CBMass1", psi_mass, mass_mean, mass_sigma1, alpha1, n);
	RooCBShape mass_myCB2("mass_myCB2", "CBMass2", psi_mass, mass_mean, mass_sigma2, alpha2, n);
	RooAddPdf mass_signal("mass_signal", "mass signal", RooArgList(mass_myCB1, mass_myCB2), RooArgList(mass_ratio));
	RooExponential mass_bkg("mass_bkg", "mass bkg", psi_mass, mass_p0);

	double ftzbkg_bias, ftzbkg_sigma1, ftzbkg_sigma2, ftzbkg_beta, ftzbkg_tau1, ftzbkg_tau2, ftzbkg_tau3, ftzbkg_tau4, ftzbkg_tau5, ftzbkg_f1, ftzbkg_f2, ftzbkg_f3, ftzbkg_f4, ftzbkg_f5;
	ifstream ftzbkg(tzbkg_fit_result.c_str());
	ftzbkg >> ftzbkg_bias >> ftzbkg_sigma1 >> ftzbkg_sigma2 >> ftzbkg_beta >> ftzbkg_tau1 >> ftzbkg_tau2 >> ftzbkg_tau3 >> ftzbkg_tau4 >> ftzbkg_tau5 >> ftzbkg_f1 >> ftzbkg_f2 >> ftzbkg_f3 >> ftzbkg_f4 >> ftzbkg_f5;
	ftzbkg.close();
	cout << ftzbkg_bias << " " << ftzbkg_sigma1 << " " << ftzbkg_sigma2 << " " << ftzbkg_beta << " " << ftzbkg_tau1 << " " << ftzbkg_tau2 << " " << ftzbkg_tau3 << " " << ftzbkg_tau4 << " " << ftzbkg_tau5 << " " << ftzbkg_f1 << " " << ftzbkg_f2 << " " << ftzbkg_f3 << " " << ftzbkg_f4 << " " << ftzbkg_f5 << endl;
	RooRealVar tzbkg_bias("tzbkg_bias", "tzbkg bias", ftzbkg_bias);
	RooRealVar tzbkg_sigma1("tzbkg_sigma1", "tzbkg sigma1", ftzbkg_sigma1);
	RooRealVar tzbkg_sigma2("tzbkg_sigma2", "tzbkg sigma2", ftzbkg_sigma2);
	RooRealVar tzbkg_beta("tzbkg_beta", "tzbkg beta", ftzbkg_beta);
	RooRealVar tzbkg_tau1("tzbkg_tau1", "tzbkg tau1", ftzbkg_tau1);
	RooRealVar tzbkg_tau2("tzbkg_tau2", "tzbkg tau2", ftzbkg_tau2);
	RooRealVar tzbkg_tau3("tzbkg_tau3", "tzbkg tau3", ftzbkg_tau3);
	RooRealVar tzbkg_tau4("tzbkg_tau4", "tzbkg tau4", ftzbkg_tau4);
	RooRealVar tzbkg_tau5("tzbkg_tau5", "tzbkg tau5", ftzbkg_tau5);
	RooRealVar tzbkg_f1("tzbkg_f1", "tzbkg f1", ftzbkg_f1);
	RooRealVar tzbkg_f2("tzbkg_f2", "tzbkg f2", ftzbkg_f2);
	RooRealVar tzbkg_f3("tzbkg_f3", "tzbkg f3", ftzbkg_f3);
	RooRealVar tzbkg_f4("tzbkg_f4", "tzbkg f4", ftzbkg_f4);
	RooRealVar tzbkg_f5("tzbkg_f5", "tzbkg f5", ftzbkg_f5);
	RooRealVar zero("zero", "zero", 0);
	RooRealVar one("one", "one", 1);
	RooGaussModel tzbkg_res1("tzbkg_res1", "resolution model 1 for tzbkg", tz, tzbkg_bias, tzbkg_sigma1, one, tzError);
	RooGaussModel tzbkg_res2("tzbkg_res2", "resolution model 2 for tzbkg", tz, tzbkg_bias, tzbkg_sigma2, one, tzError);
	RooAddModel tzbkg_res("tzbkg_res", "resolution model for tzbkg", RooArgList(tzbkg_res1, tzbkg_res2), RooArgList(tzbkg_beta));
	RooDecay tzbkg_tzP1("tzbkg_tzP1", "tz decay function P1 for tzbkg", tz, tzbkg_tau1, tzbkg_res, RooDecay::SingleSided);
	RooDecay tzbkg_tzP2("tzbkg_tzP2", "tz decay function P2 for tzbkg", tz, tzbkg_tau2, tzbkg_res, RooDecay::SingleSided);
	RooDecay tzbkg_tzP3("tzbkg_tzP3", "tz decay function P3 for tzbkg", tz, tzbkg_tau5, tzbkg_res, RooDecay::SingleSided);
	RooDecay tzbkg_tzN("tzbkg_tzN", "tz decay function N for tzbkg", tz, tzbkg_tau3, tzbkg_res, RooDecay::Flipped);
	RooDecay tzbkg_tzPN("tzbkg_tzPN", "tz decay function PN for tzbkg", tz, tzbkg_tau4, tzbkg_res, RooDecay::DoubleSided);
	RooDecay tzbkg_tzDelta("tzbkg_tzDelta", "tz delta function for tzbkg", tz, zero, tzbkg_res, RooDecay::DoubleSided);
	RooAddPdf tzbkg_tzPDF("tzbkg_tzPDF", "tz pdf for tzbkg", RooArgList(tzbkg_tzP1, tzbkg_tzP2, tzbkg_tzN, tzbkg_tzPN, tzbkg_tzP3, tzbkg_tzDelta), RooArgList(tzbkg_f1, tzbkg_f2, tzbkg_f3, tzbkg_f4, tzbkg_f5)); 

	RooRealVar bias("bias", "bias", 0, -1, 1);
	RooRealVar sigma1("sigma1", "sigma1", 2.0, 0.01, 5);
	RooRealVar sigma2("sigma2", "sigma2", 1.0, 0.01, 5);
	RooRealVar beta("beta", "beta", 0.4, 0.0, 1.0);
	RooRealVar tau("tau", "tau", 0.4, 0, 20);
	RooGaussModel res1("res1", "resolution model 1 for tz", tz, bias, sigma1, one, tzError);
	RooGaussModel res2("res2", "resolution model 2 for tz", tz, bias, sigma2, one, tzError);
	RooAddModel res("res", "resolution model for tz", RooArgList(res1, res2), RooArgList(beta));
	RooDecay tzDelta("tzDelta", "tz delta function for tz", tz, zero, res, RooDecay::DoubleSided);
	RooDecay tzP("tzP", "tz decay function P for tz", tz, tau, res, RooDecay::SingleSided);
	RooHistPdf tzTail("tzTail", "tz tail for tz", tz, *tail_hist);

	RooRealVar Nprompt("Nprompt", "Nprompt", 0.7*N, 0, 1.2*N);
	RooRealVar Nfromb("Nfromb", "Nfromb", 0.1*N, 0, N);
	RooRealVar Ntail("Ntail", "Ntail", 0.01*N, 0, N);
	RooRealVar Nbkg("Nbkg", "Nbkg", 0.1*N, 0, N);

	RooProdPdf prompt("prompt", "pdf prompt", RooArgList(mass_signal, tzDelta));
	RooProdPdf fromb("fromb", "pdf fromb", RooArgList(mass_signal, tzP));
	RooProdPdf tail("tail", "pdf tail", RooArgList(mass_signal, tzTail));
	RooProdPdf bkg("bkg", "pdf bkg", RooArgList(mass_bkg, tzbkg_tzPDF));

	RooAddPdf tzPDF("tzPDF", "pdf tz", RooArgList(prompt, fromb, tail, bkg), RooArgList(Nprompt, Nfromb, Ntail, Nbkg));

	bool flag;
	config >> flag;
	if (flag) {
		beta.setVal(0);
		beta.setConstant(true);
		sigma1.setConstant(true);
	}
	config >> flag;
	if (flag) {
		Ntail.setVal(0);
		Ntail.setConstant(true);
	}

	// Fitting
	tzPDF.fitTo(data_set, Minos(true), Strategy(2), NumCPU(40), Save(false), ConditionalObservables(tzError));
	//tzPDF.getParameters(data_set);
	//tzPDF.fitTo(data_set, Minos(true), Strategy(2), NumCPU(40), Save(false), ConditionalObservables(tzError));
	data_set.Write();
	tzPDF.Write();
	root_file->Close();

	// Output
	ofstream fout(output_pyfile_config.c_str());

	fout << "bias = " << bias.getVal() << endl;
	fout << "biasError = " << bias.getError() << endl;
	fout << "sigma1 = " << sigma1.getVal() << endl;
	fout << "sigma1Error = " << sigma1.getError() << endl;
	fout << "sigma2 = " << sigma2.getVal() << endl;
	fout << "sigma2Error = " << sigma2.getError() << endl;
	fout << "beta = " << beta.getVal() << endl;
	fout << "betaError = " << beta.getError() << endl;
	fout << "tau = " << tau.getVal() << endl;
	fout << "tauError = " << tau.getError() << endl;
	fout << "Nprompt = " << Nprompt.getVal() << endl;
	fout << "NpromptError = " << Nprompt.getError() << endl;
	fout << "Nfromb = " << Nfromb.getVal() << endl;
	fout << "NfrombError = " << Nfromb.getError() << endl;
	fout << "Ntail = " << Ntail.getVal() << endl;
	fout << "NtailError = " << Ntail.getError() << endl;
	fout << "Nbkg = " << Nbkg.getVal() << endl;
	fout << "NbkgError = " << Nbkg.getError() << endl;

	fout.close();

	TFile* splot_file = new TFile(splot_file_config.c_str(), "recreate");
	RooStats::SPlot* splot = new RooStats::SPlot("splot", "splot", data_set, &tzPDF, RooArgList(Nprompt, Nfromb, Ntail, Nbkg));
	TTree* result = data_tree->CloneTree(0);
	double prompt_sw, fromb_sw, tail_sw, bkg_sw;
	result->Branch("prompt_sw", &prompt_sw, "prompt_sw/D");
	result->Branch("fromb_sw", &fromb_sw, "fromb_sw/D");
	result->Branch("tail_sw", &tail_sw, "tail_sw/D");
	result->Branch("bkg_sw", &bkg_sw, "bkg_sw/D");
	for (int i = 0; i < N; i++) {
		data_tree->GetEntry(i);
		prompt_sw = ((RooRealVar*)(data_set.get(i)->find("Nprompt_sw")))->getVal();
		fromb_sw = ((RooRealVar*)(data_set.get(i)->find("Nfromb_sw")))->getVal();
		tail_sw = ((RooRealVar*)(data_set.get(i)->find("Ntail_sw")))->getVal();
		bkg_sw = ((RooRealVar*)(data_set.get(i)->find("Nbkg_sw")))->getVal();
		result->Fill();
	}
	result->Write();
	splot_file->Close();

	return 0;
}
