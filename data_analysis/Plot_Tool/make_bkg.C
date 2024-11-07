
#include "1d2d_constants.h"


void GenBackGround(int  Npairs,TH2D *hraw,TH2D *hbkg){
    int backMult =10;

    long  NENT =  Npairs;
    long  XENT =  ((1+floor(sqrt(1+(4*2*backMult*NENT))))/2) ;
                
    //float A_ETA_Cor[XENT] = {0};
    //float A_PHI_Cor[XENT] = {0};
    std::vector<float> A_ETA_Cor(XENT);
    std::vector<float> A_PHI_Cor(XENT);

    for(int x = 0; x<XENT; x++){
        gRandom->SetSeed(0);

        double WEta1_Cor, WPhi1_Cor;//making the pseudoparticles
        //hEPDrawCor[wtrk-1][wppt-1][wpPU-1]->GetRandom2(WEta1_Cor, WPhi1_Cor);
        hraw->GetRandom2(WEta1_Cor, WPhi1_Cor);
        A_ETA_Cor[x] = WEta1_Cor;
        A_PHI_Cor[x] = WPhi1_Cor;
    }
    for(long int i = 0; i < (XENT-1); i++){
        for(long int j = (i+1); j < XENT; j++){

                double WdeltaEta_Cor = (A_ETA_Cor[i]-A_ETA_Cor[j]);
                double WdeltaPhi_Cor = (TMath::ACos(TMath::Cos(A_PHI_Cor[i]-A_PHI_Cor[j])));
                hbkg->Fill(WdeltaEta_Cor, WdeltaPhi_Cor,   1);//./XENT);
                hbkg->Fill(-WdeltaEta_Cor, WdeltaPhi_Cor,  1);//./XENT);
                hbkg->Fill(WdeltaEta_Cor, -WdeltaPhi_Cor,  1);//./XENT);
                hbkg->Fill(-WdeltaEta_Cor, -WdeltaPhi_Cor, 1);//./XENT);
                hbkg->Fill(WdeltaEta_Cor, 2*TMath::Pi() - WdeltaPhi_Cor, 1);//./XENT);
                hbkg->Fill(-WdeltaEta_Cor,2*TMath::Pi() - WdeltaPhi_Cor, 1);//./XENT);

        }
    }

}


void make_bkg(){
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

    for(int i=0;i<SDclassifer_bin;i++){
        std::cout<<SDclassifer_low[i]<<std::endl;
    }


   	int PUbin=1;

    TFile *fin=TFile::Open("ana_Run3_2024.root");

	TH2D* hEPDrawCor[beta_bin][SDclassifer_bin][trackbin][ptbin][PUbin];
    TH2D* hSignalShiftedCor[beta_bin][SDclassifer_bin][trackbin][ptbin][PUbin];
    TH2D* hBckrndShiftedCor[beta_bin][SDclassifer_bin][trackbin][ptbin][PUbin];



    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        for(int ibeta=1; ibeta<beta_bin+1;ibeta++){
            for(int iclass=1;iclass<SDclassifer_bin+1;iclass++){
                 for(int wppt = 1; wppt<ptbin+1; wppt++){
                     for(int wpPU = 1; wpPU<PUbin+1; wpPU++){
                        hBckrndShiftedCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]   = (TH2D*)fin->Get(Form("hBckrndS_Cor_Nch_%d_%d_beta_%d_SDClass_%d_%d_pt_%d_%d_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)beta_SD[ibeta-1],(int)(100*SDclassifer_low[iclass-1]),(int)(100*SDclassifer_high[iclass-1]),(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ) );
                        hSignalShiftedCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]   = (TH2D*)fin->Get(Form("hSignalS_Cor_Nch_%d_%d_beta_%d_SDClass_%d_%d_pt_%d_%d_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)beta_SD[ibeta-1],(int)(100*SDclassifer_low[iclass-1]),(int)(100*SDclassifer_high[iclass-1]),(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ) );
                               hEPDrawCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]   = (TH2D*)fin->Get(Form("hEPDraw_Cor_Nch_%d_%d_beta_%d_SDClass_%d_%d_pt_%d_%d_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)beta_SD[ibeta-1],(int)(100*SDclassifer_low[iclass-1]),(int)(100*SDclassifer_high[iclass-1]),(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ) );
                   
                        if(wtrk>6){ hBckrndShiftedCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]->Reset();}
                    }
                }
            }
        }
   }


   TH2D *hPairs[beta_bin][SDclassifer_bin];
   for(int ibeta=0;ibeta<beta_bin;ibeta++){
        for(int iclass=0;iclass<SDclassifer_bin;iclass++){
            hPairs[ibeta][iclass]= (TH2D*)fin->Get( Form("hPairs_beta_%d_SDclass_%d_%d",(int)beta_SD[ibeta],(int)(100*SDclassifer_low[iclass]),(int)(100*SDclassifer_high[iclass]) ));
        }
   }


   for(int ibeta=0;ibeta<beta_bin;ibeta++){
    std::cout<<ibeta<<std::endl;
        for(int iclass=0;iclass<SDclassifer_bin;iclass++){
            for(int iNch=6;iNch<trackbin;iNch++){
                //std::cout<<"iNch="<<iNch<<std::endl;
                for(int ipt=0;ipt<ptbin;ipt++){
                    int i_PU=0;
                    std::cout<<(int)hPairs[ibeta][iclass]->GetBinContent(iNch+1,ipt+1)<<std::endl;
                    GenBackGround((int)hPairs[ibeta][iclass]->GetBinContent(iNch+1,ipt+1),hEPDrawCor[ibeta][iclass][iNch][ipt][i_PU],hBckrndShiftedCor[ibeta][iclass][iNch][ipt][i_PU]);
                }
            }
        }
    }

    TH2D* h2Jet_Pass[beta_bin]; 
    for(int ibeta=0;ibeta<beta_bin;ibeta++){
       h2Jet_Pass[ibeta]= (TH2D*)fin->Get( Form( "h2Jet_Pass_beta_%d",(int)beta_SD[ibeta] )  );
    }


    TFile *fout=TFile::Open("h2DSigANDBkg.root","RECREATE");

    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        for(int ibeta=1; ibeta<beta_bin+1;ibeta++){
            for(int iclass=1;iclass<SDclassifer_bin+1;iclass++){
                 for(int wppt = 1; wppt<ptbin+1; wppt++){
                     for(int wpPU = 1; wpPU<PUbin+1; wpPU++){

                        hBckrndShiftedCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]->Write();
                        hSignalShiftedCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]->Write();
                   
                    }
                }
            }
        }
   }

   for(int ibeta=0;ibeta<beta_bin;ibeta++){
       h2Jet_Pass[ibeta]->Write();
    }
    fin->Close();
    fout->Close();














}