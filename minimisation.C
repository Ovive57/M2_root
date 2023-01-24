//#include <vector>


{

	gRandom = new TRandom3(2);

	double nb = 10000;
    double ns = 100;
    double b_slope = 1.0;
    double moy = 3.0;
    double sig = 0.2;
    double amp = 5.0;
    double ntot = ns+nb;
    double k = ns/ntot;
    double likelihood[100];
    //vector<double> Poisson(100);    

	TF1* fb = new TF1("background","[0]*exp(-[1]*x)",0,10);
	fb->SetParameters(1,b_slope);
	
	TF1* fs = new TF1("signal", "gaus", 0,10);
	fs->SetParameters(amp,moy,sig);
	
	

	TH1F* hb = new TH1F("histb","histb",100,0,10);
	hb->FillRandom("background",nb);
	hb->SetLineColor(kBlue);
	
	
	TH1F* hs = new TH1F("hists","hists",100,0,10);
	hs->FillRandom("signal",ns);
	hs->SetLineColor(kBlue);
    
    TH1F* htot = new TH1F("histot","histot",100,0,10);
    htot->Add(hs,hb);
    
    htot->SetTitle("Histogram title");
    htot->GetXaxis()->SetTitle("X axis title");
    htot->GetYaxis()->SetTitle("Y axis title");

    
	TCanvas* c = new TCanvas;
	c->cd(1);
	c->SetLogy(1);
	htot->Draw("e");

	TFile fout("gen.root","recreate");
	htot->Write();
	fout.Close();


    
    
    TF1* f = new TF1("likelihood","[0]*([1]*gaus(2)+(1-[1])*[5]*exp(-[6]*x))",0,10);

    
    double ki = 0;
    double x = 0;
    double val_att;
    double widthbin;
    widthbin = htot->GetBinWidth(1);
    int binmin=0;
    int binmax=99;
    int bin;
    double poisson;
    int i = 0;
    //double mult_poisson=1.;
    
    
    for(ki=0;ki<=0.02;ki=ki+0.005){
        f->SetParameters(ntot*widthbin,ki,amp,moy,sig,1,b_slope);
        //TH1F* hp = new TH1F("histpois","histpois",100,0,10);
        double mult_poisson=1.;
        for (bin = binmin; bin<=binmax; bin=bin+1){
            //bin=bin+1;
            x = htot->GetBinCenter(bin);
            val_att = f->Eval(x);
            double_t binContent = htot->GetBinContent(bin);
	        poisson = TMath::Poisson(binContent,val_att);
	        
            mult_poisson=mult_poisson*poisson;
            //Poisson[bin] = fp;
             
            }
            likelihood[i] = mult_poisson;
            i = i+1 ;
            
            //hp->Draw();
            //multiplications
            //garder valeurs
        }
    // plot pour chercher le minimum de k
    
    
    
    //f->Eval(1);
    //f->Draw();
    //s(double x){
    //for (i=0; i<=10 ; i = i+0.01) //start,stop,pas (si 100 bin en 10 x : 0.1)
    
    //}
    
    //void b(double x)
    
} 
    
   
