////////////////////////////////////////////
//                                        //
//      Date:2019-10-09                   // 
//                                        //
//      Author:KangYouen                  //
//                                        //
//      Organization:DEP in Tsinghua      //
//                                        //
//     (modified from ZhaoMR's）          //
////////////////////////////////////////////

#include <fstream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <RooFit.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooAddPdf.h>
#include <RooCBShape.h>
#include <RooExponential.h>
#include <string>
#include <RooStats/SPlot.h>
#include <vector>
#include "lhcbStyle.h"
using namespace std;
using namespace RooFit;

int main(int argc, char *argv[])
//生成的可执行文件可以添加参数（argv[1],argv[2]...）
{
	setLHCbStyle();
	//输入文件，通过输入流添加(一个接口)
	ifstream config(argv[1]);
	string inputfile, plotfile, resultfile, result_for_2dfit;
	config >> inputfile >> plotfile >> resultfile >> result_for_2dfit;
	string plotfile_format;//方便改画图的格式
	//进入root文件
	TFile input_file(inputfile.c_str());
	TTree *input_tree = (TTree*)input_file.Get("DecayTree");
	//需拟合的变量
	RooRealVar psi_mass("psi_M", "psi_mass", 3096.9 - 120, 3096.9 + 120);
	//新建一个储存拟合变量的数据集只储存psi_mass变量
	RooDataSet data_set("data_set", "data set", input_tree, RooArgSet(psi_mass));
	double N = data_set.numEntries();
	cout<<"total entries = "<<N<<endl;

	
	//拟合工作
	//1.构建拟合函数（两个水晶球加上指数衰减背景）
	//1.1相关参数
	RooRealVar a("a","a",-0.00011);
	RooRealVar b("b","b",0.0085);
	RooRealVar c("c","c",2.066);

	//1.2变量（两个水晶球函数有同一个均值所以只需一个mean，P0是背景信号的指数衰减常数）
	RooRealVar mean("mean","mean mass of Jpsi",3096.9,3096.9-10,3096.9+10);
	RooRealVar sigma1("sigma1","sigma of CB1",15,0,40);
	RooRealVar sigma2("sigma2","sigma of CB2",15,0,40);
	RooFormulaVar alpha1("alpha1","alpha1","a*sigma1*sigma1+b*sigma1+c",RooArgSet(a,b,c,sigma1));
	RooFormulaVar alpha2("alpha2","alpha2","a*sigma2*sigma2+b*sigma2+c",RooArgSet(a,b,c,sigma2));
	RooRealVar n("n","n",1.0);
	RooRealVar ratio("ratio","ratio of CB1 over CB1+CB2",0.1,0,1);
	RooRealVar p0("p0","p0(for background)",0.0,-2.0,2.0);
	//如果事例数少的话直接用一个水晶球函数就好不然会出现一个水晶球函数特别的窄，尖尖得戳在中间的情形
	if (N<30000){ratio.setVal(1); ratio.setConstant(true); sigma2.setConstant(true);}

	//1.3构造函数
	//水晶球函数（两个）
	RooCBShape CB1("CBMass1", "CB_Function 1", psi_mass, mean, sigma1, alpha1, n);
	RooCBShape CB2("CBMass2", "CB_Function 2", psi_mass, mean, sigma2, alpha2, n);
	//将两个函数以一定比例混合（自动归一化）
	RooAddPdf signal("signal", "signal_mass", RooArgList(CB1,CB2),RooArgList(ratio));
	//背景信号（指数衰减）
	RooExponential bkg("bkg", "background signal", psi_mass, p0);
	//将叠加后的水晶球函数加上背景信号（当然其实这里可以三个函数一起加，见RooAddPdf）
	//且相加有两种方式，一是给定各函数占的比例，这里只需n-1个，或者给出各个信号的事例数
	RooRealVar Nsig("Nsig", "number of signal", N*0.8, 0, N);
	RooRealVar Nbkg("Nbkg", "number of background", N*0.2, 0, N);
	RooAddPdf massPDF("massPDF", "total signal of mass", RooArgList(signal,bkg), RooArgList(Nsig,Nbkg));
	//至此一个总函数就构建完了，下面将函数拟合到这个数据上
	
	//2.拟合(见fitTo()函数)(这里用两种方法拟合，不知为何)
	massPDF.fitTo(data_set, Minos(true), Strategy(1), NumCPU(20), Save(false));
	massPDF.getParameters(data_set);
	massPDF.fitTo(data_set, Minos(true), Strategy(2), NumCPU(20), Save(false));


	//3.画图
	TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas for all", 1600, 1200);//画布
	RooPlot *frame1 = psi_mass.frame(Bins(60));//函数图像
	data_set.plotOn(frame1, Name("data_set"));//把原数据画上去
	//把拟合函数的各个组分画上去
	massPDF.plotOn(frame1, Name("signal"), Components("signal"), FillColor(kRed), FillStyle(3245), VLines(), DrawOption("F"), LineWidth(0));
	massPDF.plotOn(frame1, Name("bkg"), Components("bkg"), LineColor(kViolet), LineStyle(7));
	massPDF.plotOn(frame1, Name("total"), LineColor(kBlue), LineStyle(1));
	//图画设置
	frame1->SetLabelOffset(0.005, "Y");
	frame1->SetTitleOffset(0.8,"Y");
	frame1->SetXTitle("m_{#mu^{+}#mu^{-}}[MeV/#font[12]{c}^{2}]");
	frame1->SetYTitle("Candidates per 4 MeV/#font[12]{c}^{2}");
	//说明框(AddEntry:增加说明框内的说明条,“对象”“标签”“添加方式”)
	TLegend *leg = new TLegend(0.80, 0.80, 1.0, 1.0);
	leg->SetTextFont(132);
	leg->SetBorderSize(2);
	leg->SetTextSize(0.04);
	leg->AddEntry("data_set", "data", "lep");
	leg->AddEntry("massPDF", "total fit function", "l");
	leg->AddEntry("signal", "J/#psi signal", "f");
	leg->AddEntry("bkg", "background", "l");
	//划分区间(bin)分析每个地方是高了还是低了，用柱状图表示出来（沿途误差表）
	RooHist *hist1 = frame1->pullHist();
	RooPlot *plot1 = psi_mass.frame();
	for (int i=0; i<hist1->GetN(); i++) {hist1->SetPointError(i,0,0,0.0);}
	plot1->addPlotable(hist1, "B");
	plot1->GetYaxis()->SetTitleSize(0.20);
	plot1->GetYaxis()->SetLabelSize(0.20);
	plot1->GetYaxis()->SetNdivisions(5,0,0);
	plot1->GetYaxis()->SetTitle("#Delta/#sigma");
	plot1->GetYaxis()->SetRangeUser(-5,5);
	//将画布分为两部分
	gROOT->SetBatch(true);
	Canvas1->Divide(1,2,0,0,0);
	//进入画布第2部分
	Canvas1->cd(2);
	gPad->SetTopMargin(0);
	gPad->SetLeftMargin(0.12);
	//设置2部分上边和左边的空白
	gPad->SetPad(0.03,0.02,0.97,0.77);//位置,左下右上，下大半部分画frame1
	frame1->Draw();
	leg->Draw();
	//设置1部分上边和左边的空白,上面一小部分画沿途误差表
	Canvas1->cd(1);
	gPad->SetTopMargin(0);
        gPad->SetLeftMargin(0.12);
	gPad->SetPad(0.03,0.77,0.97,0.97);
	plot1->Draw();
	//改变画图的格式，这里用.pdf因为不知道为啥png和jpg在下载之后不能看了
	plotfile_format = plotfile + ".pdf";
	Canvas1->Print(plotfile_format.c_str(),"recreate");

	//4.输出（写到结果中）
	ofstream result(resultfile.c_str());
	//输出写入的内容：
	result << "Nsig = " << Nsig.getVal() << endl;
	result << "NsigError = " << Nsig.getError() << endl;
	result << "Nbkg = " << Nbkg.getVal() << endl;
	result << "NbkgError" << Nbkg.getError() << endl;
	result << "mean = " << mean.getVal() << endl;
	result << "meanError = " << mean.getError() << endl;
	result << "sigma1 = " << sigma1.getVal() << endl;
	result << "sigma1Error = " << sigma1.getError() <<endl;
	result << "sigma2 = " << sigma2.getVal() << endl;
    result << "sigma2Error = " << sigma2.getError() <<endl;
	result << "ratio of CB1 over CB2 = " << ratio.getVal() << endl;
	result << "ratioError = " << ratio.getError() << endl;
	result << "attenuation constant p0 = " << p0.getVal() << endl;
	result << "p0Error" << p0.getError() << endl;
	result.close();
	ofstream txtresult(result_for_2dfit.c_str());
	txtresult << mean.getVal() << endl;
	txtresult << sigma1.getVal() << endl;
	txtresult << sigma2.getVal() << endl;
	txtresult << ratio.getVal() << endl;
	txtresult << p0.getVal() << endl;
	txtresult.close();
	return 0;
}
