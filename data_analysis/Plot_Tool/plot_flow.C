#include "1d2d_constants.h"
void plot_flow(){
	int mode=2;
	int SDclassifer_bin=0;
	double *SDclassifer_low=0;
	double *SDclassifer_high=0;
	if(mode==1){
    	SDclassifer_bin=sdZg_bin;
      	SDclassifer_low=&(sdZg_low[0]);
      	SDclassifer_high=&(sdZg_high[0]);

   	}
   	else if(mode==2){
        SDclassifer_bin=sdRg_bin;
      	SDclassifer_low=&(sdRg_low[0]);
      	SDclassifer_high=&(sdRg_high[0]);
   	}
   	else if(mode==3){
      	SDclassifer_bin=sdZgTgB_bin;
      	SDclassifer_low=&(sdZgTgB_low[0]);
      	SDclassifer_high=&(sdZgTgB_high[0]);
   	}
   	int PUbin=1;


	TFile *fin=TFile::Open("flow.root");
	TGraphErrors *gr[beta_bin][SDclassifer_bin][ptbin][PUbin];
	for(int ibeta=0;ibeta<beta_bin;ibeta++){
   		for(int iclass=0;iclass<SDclassifer_bin;iclass++){
   			for(int ipt=0;ipt<ptbin;ipt++){
   				int ipPU=0;
   				gr[ibeta][iclass][ipt][ipPU]=(TGraphErrors*)fin->Get(Form("beta_%d_iclass_%d_%d_%d",(int)beta_SD[ibeta],(int)(SDclassifer_low[iclass]*100),(int)(SDclassifer_high[iclass]*100),ipt ));
   			}
   		}
	}

	TCanvas *c1[beta_bin];
	int ipt=0;int ipPu=0;
	for(int ibeta=0;ibeta<beta_bin;ibeta++){
		c1[ibeta]=new TCanvas();
		c1[ibeta]->cd();
		TMultiGraph* mg1 = new TMultiGraph();
		for(int iclass=0;iclass<SDclassifer_bin;iclass++){
			gr[ibeta][iclass][ipt][ipPu]->SetMarkerSize(0.8);
			gr[ibeta][iclass][ipt][ipPu]->SetMarkerStyle(20);
			if(iclass==0) {gr[ibeta][iclass][ipt][ipPu]->SetLineColor(kRed);gr[ibeta][iclass][ipt][ipPu]->SetMarkerColor(kRed);}
			if(iclass==1) {gr[ibeta][iclass][ipt][ipPu]->SetLineColor(kBlue);gr[ibeta][iclass][ipt][ipPu]->SetMarkerColor(kBlue);}
			if(iclass==2) {gr[ibeta][iclass][ipt][ipPu]->SetLineColor(kBlack);gr[ibeta][iclass][ipt][ipPu]->SetMarkerColor(kBlack);}
			mg1->Add(gr[ibeta][iclass][ipt][ipPu],"P");
		}

		mg1->Draw("a");
		mg1->GetXaxis()->SetTitle("<N_{ch}^{j}>");
		mg1->GetYaxis()->SetTitle("V^{*}_{n#Delta} {2,|#Delta#eta^{*}| >2}");
		TLegend* leg1=new TLegend(0.4,0.5,0.85,0.7);
    	leg1->SetTextFont(22);
    	leg1->SetTextSize(0.05);
    	leg1->SetBorderSize(0);
    	for(int i=0;i<SDclassifer_bin;i++){
    		leg1->AddEntry(gr[ibeta][i][ipt][ipPu],Form("#beta=%d,R_{g}: %.2f ~ %.2f ",(int)beta_SD[ibeta],SDclassifer_low[i],SDclassifer_high[i]));
		}
		leg1->Draw("same");
		c1[ibeta]->SaveAs(Form("plot_flow/beta_%d.png",(int)beta_SD[ibeta]));
		c1[ibeta]->SaveAs(Form("plot_flow/beta_%d.pdf",(int)beta_SD[ibeta]));
	}

}