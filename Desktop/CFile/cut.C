#define DecayTree_cxx
#include "DecayTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;
int main(int argc, char* argv[])
{
	ifstream config(argv[1]);
	string outputfile;
	int SPD_low, SPD_high;
	config >> outputfile >> SPD_low >> SPD_high;
	//input file, using TFile or TChain,use pointer
	TChain *F_in = new TChain("tuple1/DecayTree");
	for (int i=0; i<=320; i++)
	F_in->Add(Form("../../../12/%d/output/Tuple.root",i));
	//output file
	TFile *F_out1 = new TFile(outputfile.c_str(),"recreate");
	TTree *T_out1 = F_in->CloneTree(0);
	Double_t tz, tzError;
	//要先声明这个变量，然后做个branch，在循环里每个entry给它们赋值
	T_out1->Branch("tz",&tz,"tz/D");
	T_out1->Branch("tzError",&tzError,"tzError/D");
	DecayTree t(F_in);	//define F_in as class tree(depend on which class you make in tree.h)
	
	//loop
	Long64_t Num_Entries = F_in->GetEntries();
	for(Long64_t jentry=0; jentry<=Num_Entries; jentry++)
	{
		F_in->GetEntry(jentry);	
		Double_t half_zeta_p = asin(t.mup_PT/t.mup_PE)/2;
		Double_t half_zeta_m = asin(t.mum_PT/t.mum_PE)/2;
		Double_t ETA_p = -log(tan(half_zeta_p));
		Double_t ETA_m = -log(tan(half_zeta_m));
		tz = (t.psi_ENDVERTEX_Z-t.psi_OWNPV_Z)/t.psi_PZ * 10.0 * 3096.9 / 2.9979;//c=2.9979 and 10.0 is for the unit
		tzError = t.psi_ENDVERTEX_ZERR  / t.psi_PZ * 10.0 * 3096.9 / 2.9979;
		cout << "jentry/nentries = " <<jentry << "/" << Num_Entries << endl;
		//cut
		if(!(
			t.mum_P>3000&&t.mup_P>3000&&	
			t.psi_ENDVERTEX_CHI2<8&&
			t.BCType == 3&&
			t.nSPDHits>=SPD_low&&t.nSPDHits<SPD_high&&
			t.mup_TRACK_GhostProb<0.3&&t.mum_TRACK_GhostProb<0.3&&
			t.mum_PT>650&&t.mup_PT>650&&
			tz>-10&&tz<10&&
			t.psi_MM>2976.9&&t.psi_MM<3216.9&&
			ETA_p>2.0&&ETA_p<4.9&&
			ETA_m>2.0&&ETA_m<4.9&&
			t.mup_isMuon==1&&t.mum_isMuon==1&&
			abs(t.mup_PIDmu)>0&&abs(t.mum_PIDmu)>0&&
			t.psi_L0DiMuonDecision_TOS==1&&t.nPVs!=0&&
			t.psi_Hlt1DiMuonHighMassDecision_TOS==1&&
			t.mup_isMuon == 1&&t.mum_isMuon == 1&&
			tzError <= 0.3
		))continue;
		T_out1->Fill();
	}
	T_out1->Write();
	F_out1->Close();
}
