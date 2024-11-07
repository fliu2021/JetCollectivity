#include "1d2d_constants.h"
void DrawFlow(){
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
	TFile *fin=TFile::Open("h2DSigANDBkg.root");
	
    TH2D* hSignalShiftedCor[beta_bin][SDclassifer_bin][trackbin][ptbin][PUbin];
    TH2D* hBckrndShiftedCor[beta_bin][SDclassifer_bin][trackbin][ptbin][PUbin];
    TH1D* h1DFlow[beta_bin][SDclassifer_bin][trackbin][ptbin][PUbin];
    TGraphErrors *gr[beta_bin][SDclassifer_bin][ptbin][PUbin];
	for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        for(int ibeta=1; ibeta<beta_bin+1;ibeta++){
            for(int iclass=1;iclass<SDclassifer_bin+1;iclass++){
                 for(int wppt = 1; wppt<ptbin+1; wppt++){
                     for(int wpPU = 1; wpPU<PUbin+1; wpPU++){
                        hBckrndShiftedCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]   = (TH2D*)fin->Get(Form("hBckrndS_Cor_Nch_%d_%d_beta_%d_SDClass_%d_%d_pt_%d_%d_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)beta_SD[ibeta-1],(int)(100*SDclassifer_low[iclass-1]),(int)(100*SDclassifer_high[iclass-1]),(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ) );
                        hSignalShiftedCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]   = (TH2D*)fin->Get(Form("hSignalS_Cor_Nch_%d_%d_beta_%d_SDClass_%d_%d_pt_%d_%d_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)beta_SD[ibeta-1],(int)(100*SDclassifer_low[iclass-1]),(int)(100*SDclassifer_high[iclass-1]),(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ) );                   
                        
                    }
                }
            }
        }
   	}

   	TH2D* h2Jet_Pass[beta_bin]; 
   	for(int ibeta=0;ibeta<beta_bin;ibeta++){
       	h2Jet_Pass[ibeta]= (TH2D*)fin->Get( Form( "h2Jet_Pass_beta_%d",(int)beta_SD[ibeta] )  );
   	}
	TFile *fin2=TFile::Open("ana_Run3_2024.root");
	TH1D* hBinDist_unc[beta_bin][SDclassifer_bin][trackbin];
	 for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        for(int ibeta=1; ibeta<beta_bin+1;ibeta++){
            for(int iclass=1;iclass<SDclassifer_bin+1;iclass++){
                hBinDist_unc[ibeta-1][iclass-1][wtrk-1]=(TH1D*)fin2->Get(Form("hBinDist_unc_%d_%d_%d",ibeta,iclass,wtrk)) ;
            }     
        }    
    }
    

	int YPlo=28;
    int YPhi=34;
	for(int ibeta=0;ibeta<beta_bin;ibeta++){
		for(int iclass=0;iclass<SDclassifer_bin;iclass++){
			for(int ipt=0;ipt<ptbin;ipt++){
				int ipPu=0;
				double Nch_Ave[9];
				double Nch_err[9]={0.};
				double v2[9];
				double v2err[9];
				for(int iNch=0;iNch<9;iNch++){
					double eta_bw=hSignalShiftedCor[ibeta][iclass][iNch][ipt][ipPu]->GetXaxis()->GetBinWidth(1);//x:delta eta; y:delta phi
            		double phi_bw=hSignalShiftedCor[ibeta][iclass][iNch][ipt][ipPu]->GetYaxis()->GetBinWidth(1);

            		hSignalShiftedCor[ibeta][iclass][iNch][ipt][ipPu]->Scale(1.0/(h2Jet_Pass[ibeta]->GetBinContent(iNch+1,iclass+1)));
            		hSignalShiftedCor[ibeta][iclass][iNch][ipt][ipPu]->Scale(1./(phi_bw));

            		TH1D *histfit1 = (TH1D*) hSignalShiftedCor[ibeta][iclass][iNch][ipt][ipPu]->ProjectionY("",YPlo,YPhi)->Clone();
            		TH1D *histfit2 = (TH1D*) hBckrndShiftedCor[ibeta][iclass][iNch][ipt][ipPu]->ProjectionY("",YPlo,YPhi)->Clone();

            		histfit1->Divide(histfit2);
            		histfit1->Scale(hBckrndShiftedCor[ibeta][iclass][iNch][ipt][ipPu]->GetMaximum());

            		h1DFlow[ibeta][iclass][iNch][ipt][ipPu]=(TH1D*)histfit1->Clone(); 
            		std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)))";
            		TF1 func1("deltaPhi1", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
            		func1.SetParameter(0, histfit1->GetMaximum());
            		func1.SetParameter(1, 0.1);
            		func1.SetParameter(2, 0.1);
            		func1.SetParameter(3, 0.1);
            		func1.SetParameter(4, 0.1);
            		func1.SetParameter(5, 0.1);
            		h1DFlow[ibeta][iclass][iNch][ipt][ipPu]->Fit(&func1, "m E q");
            		v2[iNch]=func1.GetParameter(2);
                	v2err[iNch]=func1.GetParError(2)*TMath::Sqrt(2);
                	Nch_Ave[iNch]=hBinDist_unc[ibeta][iclass][iNch]->GetMean();
				}

				gr[ibeta][iclass][ipt][ipPu]=new TGraphErrors(9,Nch_Ave,v2,Nch_err,v2err);
				gr[ibeta][iclass][ipt][ipPu]->SetName(Form("beta_%d_iclass_%d_%d_%d",(int)beta_SD[ibeta],(int)(SDclassifer_low[iclass]*100),(int)(SDclassifer_high[iclass]*100),ipt ));






			}
		}
   	}
   	TFile *fout=TFile::Open("flow.root","RECREATE");
   	for(int ibeta=0;ibeta<beta_bin;ibeta++){
   		for(int iclass=0;iclass<SDclassifer_bin;iclass++){
   			for(int ipt=0;ipt<ptbin;ipt++){
   				int ipPU=0;
   				gr[ibeta][iclass][ipt][ipPU]->Write();
   			}
   		}
   	}
   	fout->Write();

  






}