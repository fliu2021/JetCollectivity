#include "1d2d_constants.h"

void plot_substructure(){
	TFile *fin=TFile::Open("ana_Run3_2024.root");
	TH2D *hNch_Zg[beta_bin];                                                                                              
    TH2D *hNch_Rg[beta_bin];                                                                                                
    TH2D *hNch_ZgTgB[beta_bin];                                                                                            
    TH2D *hZg_Rg[beta_bin];

    TH1D *hZg_dis[beta_bin];
    TH1D *hRg_dis[beta_bin];
    TH1D *hZgTgB_dis[beta_bin];
    TH1D *hSDJetMass_dis[beta_bin];

                                                                                                                                                                                       //
    //for different Nch
    TH1D *hZgTgB_dis_Nch[beta_bin][trackbin];                                                                                       
    TH1D *hZg_dis_Nch[beta_bin][trackbin];                                                                                     
    TH1D *hRg_dis_Nch[beta_bin][trackbin];  
    TH1D *hSDJetMass_dis_Nch[beta_bin][trackbin]; 

    TH1D *hStoredSDJetMassNch=(TH1D*)fin->Get("hStoredSDJetMassNch"); 
    TH1D *hCalSDJetMassNch=(TH1D*)fin->Get("hCalSDJetMassNch"); 

    for(int ibeta=0;ibeta<beta_bin;ibeta++){




        hNch_ZgTgB[ibeta]=(TH2D*)fin->Get(Form("Nch_VS_ZgThetagBeta_beta_%d",(int)beta_SD[ibeta]) );
        hNch_Rg[ibeta]=(TH2D*)fin->Get(Form("Nch_VS_R_g_beta_%d",(int)beta_SD[ibeta])             );
        hNch_Zg[ibeta]=(TH2D*)fin->Get(Form("Nch_VS_Z_g_beta_%d",(int)beta_SD[ibeta])             );
        hZg_Rg[ibeta]=(TH2D*)fin->Get(Form("hZg_Vs_Rg_beta_%d",(int)beta_SD[ibeta])               ); 
        
        hZg_dis[ibeta]=(TH1D*)fin->Get(Form("hZg_dis_beta_%d",  (int)beta_SD[ibeta])              );
        hRg_dis[ibeta]=(TH1D*)fin->Get(Form("hRg_dis_beta_%d",  (int)beta_SD[ibeta])              );
        hZgTgB_dis[ibeta]=(TH1D*)fin->Get(Form("hZgTgB_dis_beta_%d",  (int)beta_SD[ibeta])        );
        hSDJetMass_dis[ibeta]=(TH1D*)fin->Get(Form("hSDJetMass_dis_beta_%d",(int)beta_SD[ibeta])  );

        for(int iNch=0;iNch<trackbin;iNch++){
            hZgTgB_dis_Nch[ibeta][iNch]=(TH1D*)fin->Get(Form("ZgThetagBeta_dis_Nch_%d_%d_beta_%d",trackbinbounds[iNch],trackbinboundsUpper[iNch],(int)beta_SD[ibeta]) );
            hZg_dis_Nch[ibeta][iNch]=(TH1D*)fin->Get(Form("ZgDis_Nch_%d_%d_beta_%d",trackbinbounds[iNch],trackbinboundsUpper[iNch],(int)beta_SD[ibeta])               );
            hRg_dis_Nch[ibeta][iNch]=(TH1D*)fin->Get(Form("RgDis_Nch_%d_%d_beta_%d",trackbinbounds[iNch],trackbinboundsUpper[iNch],(int)beta_SD[ibeta])               );
            hSDJetMass_dis_Nch[ibeta][iNch]=(TH1D*)fin->Get(Form("hSDJetMass_dis_beta_%d_Nch_%d_%d",(int)beta_SD[ibeta],trackbinbounds[iNch],trackbinboundsUpper[iNch]));

        }
    }
    /*
    TCanvas *c1[beta_bin];
    for(int i=0;i<beta_bin;i++){
    	c1[i]=new TCanvas();
    	c1[i]->cd();
    	hNch_Zg[i]->Draw("colz");
    }

    TCanvas *c2[beta_bin];
    for(int i=0;i<beta_bin;i++){
    	c2[i]=new TCanvas();
    	c2[i]->cd();
    	hNch_Rg[i]->Draw("colz");
    }
    TCanvas *c3[beta_bin];
    for(int i=0;i<beta_bin;i++){
    	c3[i]=new TCanvas();
    	c3[i]->cd();
    	hNch_ZgTgB[i]->Draw("colz");
    }

    TCanvas *c4[beta_bin];
    for(int i=0;i<beta_bin;i++){
    	c4[i]=new TCanvas();
    	c4[i]->cd();
    	hZg_Rg[i]->Draw("colz");
    }


    */

    TCanvas *c5=new TCanvas();
    c5->Divide(3,3);
    for(int i=0;i<9;i++){
    	c5->cd(i+1);
        if(i<3) gPad->SetLogy();
    	for(int j=beta_bin-1;j>=0;j--){
    		hZgTgB_dis_Nch[j][i]->GetXaxis()->SetTitle("Zg(#theta_{g})^{#beta}");
            hZgTgB_dis_Nch[j][i]->GetXaxis()->SetTitleSize(0.05);
            hZgTgB_dis_Nch[j][i]->GetXaxis()->SetRangeUser(0,0.7);
    		hZgTgB_dis_Nch[j][i]->GetYaxis()->SetTitle("");
    		hZgTgB_dis_Nch[j][i]->Scale(1./hZgTgB_dis_Nch[j][i]->Integral(),"width");
    		if(j==0)	hZgTgB_dis_Nch[j][i]->SetLineColor(kBlue);
    		if(j==1)    hZgTgB_dis_Nch[j][i]->SetLineColor(kGreen);
    		if(j==2)    hZgTgB_dis_Nch[j][i]->SetLineColor(kRed);
    		hZgTgB_dis_Nch[j][i]->Draw("same");
    	}
    }

    TLegend* leg5=new TLegend(0.6,0.5,0.85,0.7);
    leg5->SetTextFont(22);
    leg5->SetTextSize(0.05);
    leg5->SetBorderSize(0);
    for(int i=0;i<beta_bin;i++){
        leg5->AddEntry(hZgTgB_dis_Nch[i][0],Form("#beta=%d",(int)beta_SD[i]) );
    }
    c5->cd(1);
    leg5->Draw("same");
    c5->SaveAs("plot_substructure/ZgTgB_dis_Nch.pdf");
   	TCanvas *c6=new TCanvas();
    c6->Divide(3,3);
    for(int i=0;i<9;i++){
    	c6->cd(i+1);
    	for(int j=0;j<beta_bin;j++){
    		hZg_dis_Nch[j][i]->GetXaxis()->SetTitle("Zg");
            hZg_dis_Nch[j][i]->GetXaxis()->SetTitleSize(0.05);
    		hZg_dis_Nch[j][i]->GetYaxis()->SetTitle("");
    		hZg_dis_Nch[j][i]->Scale(1./hZg_dis_Nch[j][i]->Integral(),"width");
    		if(j==0)	hZg_dis_Nch[j][i]->SetLineColor(kBlue);
    		if(j==1)    hZg_dis_Nch[j][i]->SetLineColor(kGreen);
    		if(j==2)    hZg_dis_Nch[j][i]->SetLineColor(kRed);
    		hZg_dis_Nch[j][i]->Draw("same");
    	}
    }

    TLegend* leg6=new TLegend(0.6,0.5,0.85,0.7);
    leg6->SetTextFont(22);
    leg6->SetTextSize(0.05);
    leg6->SetBorderSize(0);
    for(int i=0;i<beta_bin;i++){
        leg6->AddEntry(hZg_dis_Nch[i][0],Form("#beta=%d",(int)beta_SD[i]) );
    }
    c6->cd(1);
    leg6->Draw("same");
    c6->SaveAs("plot_substructure/ZgB_dis_Nch.pdf");
    TCanvas *c7=new TCanvas();
    c7->Divide(3,3);
    for(int i=0;i<9;i++){
    	c7->cd(i+1);
    	for(int j=0;j<beta_bin;j++){
    		hRg_dis_Nch[j][i]->GetXaxis()->SetTitle("Rg");
            hRg_dis_Nch[j][i]->GetXaxis()->SetTitleSize(0.05);
    		hRg_dis_Nch[j][i]->GetYaxis()->SetTitle("");
    		hRg_dis_Nch[j][i]->Scale(1./hRg_dis_Nch[j][i]->Integral(),"width");
    		if(j==0)	hRg_dis_Nch[j][i]->SetLineColor(kBlue);
    		if(j==1)    hRg_dis_Nch[j][i]->SetLineColor(kGreen);
    		if(j==2)    hRg_dis_Nch[j][i]->SetLineColor(kRed);
    		hRg_dis_Nch[j][i]->Draw("same");
    	}
    }

    TLegend* leg7=new TLegend(0.6,0.5,0.85,0.7);
    leg7->SetTextFont(22);
    leg7->SetTextSize(0.05);
    leg7->SetBorderSize(0);
    for(int i=0;i<beta_bin;i++){
        leg7->AddEntry(hRg_dis_Nch[i][0],Form("#beta=%d",(int)beta_SD[i]) );
    }
    c7->cd(1);
    leg7->Draw("same");
    c7->SaveAs("plot_substructure/Rg_dis_Nch.pdf");


    TCanvas *c8=new TCanvas();
    
    c8->Divide(4,3);
    for(int ibeta=0;ibeta<beta_bin;ibeta++){
        c8->cd(ibeta*4+1);
        gPad->SetLogz();
        hNch_Zg[ibeta]->GetXaxis()->SetTitle("Nch");
        hNch_Zg[ibeta]->GetYaxis()->SetTitle("Zg");
        hNch_Zg[ibeta]->Draw("colz");
        c8->cd(ibeta*4+2);
        gPad->SetLogz();
        hNch_Rg[ibeta]->GetXaxis()->SetTitle("Nch");
        hNch_Rg[ibeta]->GetYaxis()->SetTitle("Rg");
        hNch_Rg[ibeta]->Draw("colz");
        c8->cd(ibeta*4+3);
        gPad->SetLogz();
        hNch_ZgTgB[ibeta]->GetXaxis()->SetTitle("Nch");
        hNch_ZgTgB[ibeta]->GetYaxis()->SetTitle("Zg(#theta g)^{#beta}");
        hNch_ZgTgB[ibeta]->Draw("colz");
        c8->cd(ibeta*4+4);
        gPad->SetLogz();
        hZg_Rg[ibeta]->GetXaxis()->SetTitle("Zg");
        hZg_Rg[ibeta]->GetYaxis()->SetTitle("Rg");
        hZg_Rg[ibeta]->Draw("colz");

       
    }
    c8->SaveAs("plot_substructure/Zg_Rg_Nch.pdf");
    TCanvas *c9=new TCanvas();
    c9->Divide(2,1);
    c9->cd(1);
    gPad->SetLogy();
    hStoredSDJetMassNch->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    hStoredSDJetMassNch->SetLineColor(kGreen);
    hStoredSDJetMassNch->Draw();
    hCalSDJetMassNch->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    hCalSDJetMassNch->SetLineColor(kBlue);
    hCalSDJetMassNch->Draw("same");
    c9->cd(2);
    hStoredSDJetMassNch->Rebin(5);
    hCalSDJetMassNch->Rebin(5);
    TH1D *h_ratio=(TH1D *)hCalSDJetMassNch->Clone("ratio_Calculated_Stored");
    h_ratio->Divide(hStoredSDJetMassNch);
    h_ratio->SetTitle("ratio_Calculated_Stored");
    h_ratio->SetLineColor(kViolet);
    h_ratio->Draw();

    TLegend* leg9=new TLegend(0.1,0.7,0.9,0.9);
    leg9->SetTextFont(22);
    leg9->SetTextSize(0.05);
    leg9->SetBorderSize(0);
    leg9->AddEntry(hStoredSDJetMassNch,"Stored_SoftDrop_Jet_Mass");
    leg9->AddEntry(hCalSDJetMassNch,"Calculated_SoftDrop_Jet_Mass");
    leg9->AddEntry(h_ratio,"ratio[Calculated/Stored]");
    c9->cd(1);
    leg9->Draw("same");
    c9->SaveAs("plot_substructure/SDJetMass_forCheck.pdf");















}