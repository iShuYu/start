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
#include <RooFormulaVar.h>
#include <RooDataSet.h>
#include <RooGaussModel.h>
#include <RooAddModel.h>
#include <RooAddPdf.h>
#include <RooDecay.h>
#include <string>
#include <vector>
#include <fstream>

using namespace std;
using namespace RooFit;

int main(int argc, char* argv[])
{

	//------------------------File-------------------
	ifstream config(argv[1]);
	string input_file_config, root_file_config, output_pyfile_config, output_txtfile_config, plot_output;
	config >> input_file_config >> root_file_config >> output_pyfile_config >> output_txtfile_config;

	//double pt_low, pt_hi, y_low, y_hi;
	//config >> pt_low >> pt_hi >> y_low >> y_hi;

	// Create DecayTree with the data after reconstruction
	TFile input_file(input_file_config.c_str());
	TTree* data_tree = (TTree*)input_file.Get("DecayTree");

	TFile* root_file = new TFile(root_file_config.c_str(), "recreate");

	// Useful Variables
	RooRealVar tz("tz", "tz", -20, 20);
	RooRealVar tzError("tzError", "tzError", 0.001, 0.3);
	RooRealVar bkgSweight("bkg_sw", "bkg_sw", -100, 100);

	// Data Set to Save the Fit Variables
	RooDataSet data_set("data_set", "data set", data_tree, RooArgSet(tz, tzError));

	long N = data_set.numEntries();
	cout << "Total Entries = " << N << endl;
	cout << pt_low << " " << pt_hi << endl;
	cout << y_low << " " << y_hi << endl;

	//---------------------pdf-----------------------
	RooRealVar zero("zero", "zero", 0);
	RooRealVar one("one", "one", 1.0);
	RooRealVar bias("bias", "bias", 0, -1, 1);
	RooRealVar sigma1("sigma1", "sigma1", 1.0, 0.01, 5);
	RooRealVar sigma2("sigma2", "sigma2", 2.0, 0.1, 10);
	RooRealVar beta("beta", "beta", 0.4, 0.0, 1.0);
	RooRealVar tau1("tau1", "tau1", 0.6, 0, 5);
	RooRealVar tau2("tau2", "tau2", 1.5, 0, 5);
	RooRealVar tau3("tau3", "tau3", 1.6, 0, 10);
	RooRealVar tau4("tau4", "tau4", 10, 0.1, 150);
	RooRealVar tau5("tau5", "tau5", 3, 0, 20);
	RooRealVar f1("f1", "f1", 0.15, 0.0, 1.0);
	RooRealVar f2("f2", "f2", 0.04, 0.0, 1.0);
	RooRealVar f3("f3", "f3", 0.03, 0.0, 1.0);
	RooRealVar f4("f4", "f4", 0.02, 0.0, 1.0);
	RooRealVar f5("f5", "f5", 0.01, 0.0, 1.0);

	RooGaussModel res1("res1", "resolution model 1", tz, bias, sigma1, one, tzError);
	RooGaussModel res2("res2", "resolution model 2", tz, bias, sigma2, one, tzError);
	RooAddModel res("res", "resolution model", RooArgList(res1, res2), RooArgList(beta));
	RooDecay tzP1("tzP1", "tz decay function P1", tz, tau1, res, RooDecay::SingleSided);
	RooDecay tzP2("tzP2", "tz decay function P2", tz, tau2, res, RooDecay::SingleSided);
	RooDecay tzP3("tzP3", "tz decay function P3", tz, tau5, res, RooDecay::SingleSided);
	RooDecay tzN("tzN", "tz decay function N", tz, tau3, res, RooDecay::Flipped);
	RooDecay tzPN("tzPN", "tz decay function PN", tz, tau4, res, RooDecay::DoubleSided);
	RooDecay tzDelta("tzDelta", "tz delta function", tz, zero, res, RooDecay::DoubleSided);
	RooAddPdf tzPDF("tzPDF", "tz pdf", RooArgList(tzP1, tzP2, tzN, tzPN, tzP3, tzDelta), RooArgList(f1, f2, f3, f4, f5)); 
	
	bool flag;
	config >> flag;
	if (flag) {
		beta.setVal(0);
		beta.setConstant(true);
		sigma1.setConstant(true);
	}
	config >> flag;
	if (flag) {
		f1.setVal(0);
		f1.setConstant(true);
		tau1.setConstant(true);
	}
	config >> flag;
	if (flag) {
		f2.setVal(0);
		f2.setConstant(true);
		tau2.setConstant(true);
	}
	config >> flag;
	if (flag) {
		f3.setVal(0);
		f3.setConstant(true);
		tau3.setConstant(true);
	}
	config >> flag;
	if (flag) {
		f4.setVal(0);
		f4.setConstant(true);
		tau4.setConstant(true);
	}
	config >> flag;
	if (flag) {
		f5.setVal(0);
		f5.setConstant(true);
		tau5.setConstant(true);
	}
	
	// Fitting
	tzPDF.fitTo(data_set, Minos(true), Strategy(2), NumCPU(60), Save(false), ConditionalObservables(tzError));
	//tzPDF.getParameters(data_set);
	//tzPDF.fitTo(data_set, Minos(true), Strategy(2), NumCPU(60), Save(false), ConditionalObservables(tzError));

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
	fout << "tau1 = " << tau1.getVal() << endl;
	fout << "tau1Error = " << tau1.getError() << endl;
	fout << "tau2 = " << tau2.getVal() << endl;
	fout << "tau2Error = " << tau2.getError() << endl;
	fout << "tau3 = " << tau3.getVal() << endl;
	fout << "tau3Error = " << tau3.getError() << endl;
	fout << "tau4 = " << tau4.getVal() << endl;
	fout << "tau4Error = " << tau4.getError() << endl;
	fout << "tau5 = " << tau5.getVal() << endl;
	fout << "tau5Error = " << tau5.getError() << endl;
	fout << "f1 = " << f1.getVal() << endl;
	fout << "f1Error = " << f1.getError() << endl;
	fout << "f2 = " << f2.getVal() << endl;
	fout << "f2Error = " << f2.getError() << endl;
	fout << "f3 = " << f3.getVal() << endl;
	fout << "f3Error = " << f3.getError() << endl;
	fout << "f4 = " << f4.getVal() << endl;
	fout << "f4Error = " << f4.getError() << endl;
	fout << "f5 = " << f5.getVal() << endl;
	fout << "f5Error = " << f5.getError() << endl;

	fout.close();

	ofstream fresult(output_txtfile_config.c_str());

	fresult << bias.getVal() << endl;
	fresult << sigma1.getVal() << endl;
	fresult << sigma2.getVal() << endl;
	fresult << beta.getVal() << endl;
	fresult << tau1.getVal() << endl;
	fresult << tau2.getVal() << endl;
	fresult << tau3.getVal() << endl;
	fresult << tau4.getVal() << endl;
	fresult << tau5.getVal() << endl;
	fresult << f1.getVal() << endl;
	fresult << f2.getVal() << endl;
	fresult << f3.getVal() << endl;
	fresult << f4.getVal() << endl;
	fresult << f5.getVal() << endl;

	fresult.close();

	return 0;
}
