{
	gRandom = new TRandom3(2);
    
    // Variables signal/bruit
	double nb = 10000; //nombre d'évenements bruit
    double ns = 1000; //nombre d'évenements signal
    double ntot = ns+nb; //nombre d'évenements totales
    double k = ns/ntot; //paramètre signal sur bruit. Pour nous 0.09, on cherche s'approcher à ça
    
    // Variables Gauss et exponentielle pour créer signal et bruit et histogrames
    double b_slope = 1.0; // pente exponentielle (bruit)
    double moy = 3.0; // moyenne gaussienne (signal)
    double sig = 0.2; // déviation standard gaussienne (signal) 
    double amp = 1.0; // amplitude gaussienne (signal)
    int bins = 100; // Nombre de bins 
    double widthbin; // La largeur des bins, on la calcule après.
    
    // Variables pour les boucles 
    double ki; // Pour faire l'itération sur les différents k
    int bin; // Pour faire l'itération sur les bins
    
    
    // Variables nécésaires pour calculer poisson:
    double_t x; // Pour calculer le x au milieu du bin où on va évaluer la fonction f(x)
    double_t val_att; // Pour calculer la valeur attendu, i.e. f(x). (lamda)
    double_t val_mesure; // Pour calculer la valeur mesurée i.e. (n) 
    
    double_t poisson; // Fonction de poisson
    double tot_poisson; // Somme des logaritmes de toutes les fonctions poissons de chaque bin

    
    // Tableaux pour le plot
    double likelihood[100]; //likelihood des poisssons, axe y. Il faudrait voir la taille, même taille qu'iterations sur k.
    // Il manque un tableau de k pour l'axe x.
    double K[100];
    int i = 0; // Pour remplir chaque tableau pour le plot
    
    
    // Fonction du bruit : 
	TF1* fb = new TF1("background","[0]*exp(-[1]*x)",0,10);
	fb->SetParameters(1,b_slope);
	
	// Fonction du signal : 
	TF1* fs = new TF1("signal", "gaus", 0,10);
	fs->SetParameters(amp,moy,sig);
	
	// Histogramme du bruit :
	TH1F* hb = new TH1F("histb","histb",bins,0,10);
	hb->FillRandom("background",nb);
	hb->SetLineColor(kBlue);
	
	// Histogramme du signal :
	TH1F* hs = new TH1F("hists","hists",bins,0,10);
	hs->FillRandom("signal",ns);
	hs->SetLineColor(kBlue);
    
    // Les deux histogramme ensemble :
    TH1F* htot = new TH1F("histot","histot",bins,0,10);
    htot->Add(hs,hb);
    widthbin = htot->GetBinWidth(1); // On calcule la largeur des bins pour plus tard (les boucles) tous les bins sont pareils donc on choisi par exemple de mesurer la largeur du premier.
    
    // Le plot avec les titres des axes:
    
    htot->SetTitle("Signal sur bruit");
    htot->GetXaxis()->SetTitle("Temps"); // Pas sûre que ça soit temps xd
    htot->GetYaxis()->SetTitle("Nombre d'évènements");

    // Création de Canvas pour faire le plot, je n'ai aucune idée de comment ça marche, mais ce sont comme des figures.
    
	TCanvas* c = new TCanvas;
	c->cd(1);
	c->SetLogy(1);
	htot->Draw("e");

	TFile fout("gen.root","recreate");
	htot->Write();
	fout.Close();

 
    // Fonction f pour les valeurs attendues.
    TF1* f = new TF1("fontion f","[0]*([1]*gaus(2)+(1-[1])*[5]*exp(-[6]*x))",0,10);
    // Pour plot décommenter et commenter dans la boucle.
    //f->SetParameters(ntot*widthbin,k,1/(sqrt(2*3.1415)*sig),moy,sig,1,b_slope);
    //std::cout<<" ntot "<<ntot<<" width "<<widthbin<<"ki "<<k<<" amp "<<amp<<" moy "<<moy<<" sig "<<sig<<" slope "<<b_slope<<std::endl;
    //f->Draw("same");
   
    // Boucle pour chercher le maximum de vraisemblance:
    // Il faut chercher pour plusieurs k laquelle donne une probabilité la plus grande. On fait la probabilité pour tous les bins et on multiplie pour avoir la probabilité totale pour chaque k.

    for(ki=0.08;ki<=0.1;ki=ki+0.0005){
        f->SetParameters(ntot*widthbin,ki,1/(sqrt(2*3.1415)*sig),moy,sig,1,b_slope); //Commenter pour avoir le plot du fit.
        tot_poisson=0.; // On défini ici à 0 pour réinisializer à chaque fois qu'on fait un boucle sur k.
        for (bin = 0; bin<=bins-1; bin=bin+1){
            x = htot->GetBinCenter(bin+1);
            val_att = f->Eval(x);

            // Pour la valeur mesurée on prend la valeur du bin, combien de comptages il y a par bin.
            val_mesure = htot->GetBinContent(bin+1); //GetBinContent a minimum 1 bin, donc il peut pas avoir 0 comme argument.
	           	                    //std::cout<<val_mesure<<",,,,,,"<<val_att<<std::endl;
	        // Une fois on a la valeur attendue et la valeur mesuré on calcule avec Poisson la probabilité de trouver la valeur attendue en ayant la valeur mesurée P(mesurée,attendue):
	        poisson = TMath::Poisson(val_mesure,val_att);
    	                    //std::cout<<poisson<<std::endl;
	        // On fait le log : 
	        poisson = TMath::Log(poisson);
	                    //std::cout<<val_att<<std::endl;
	        // Comme on a fait le log on ne multiplie pas on somme:
            tot_poisson += poisson;
            }

    // On veut le -log du coup on stocke avec un - :
    likelihood[i] = - tot_poisson; //axe y, il faut aussi l'axe x pour le plot. L'axe x ce sont les k.
    K[i]=ki;
    //std::cout<<K[i]<<std::endl;
    i = i+1 ; // On actualise la valeur de i pour bouger dans les vecteurs likelihood et K
    // i en plus va être la quantité de points qu'on a, on va l'utiliser pour plot après !
     

        }
    

   TCanvas* c1 = new TCanvas;
   c1->cd(1);

   TGraph* g = new TGraph(i,K,likelihood); // ici on utilise i, parce que c'est le nombre de points qu'on a.
   //TAxis *axis = g->GetXaxis();
   //axis->SetLimits(0.084,0.1);
   g->SetTitle("Graph title;X title;Y title");
   //g->GetHistogram()->SetMaximum(236.);   // along          
   //g->GetHistogram()->SetMinimum(234.4);
   g->Draw("AC*");   //g->Draw("AC*");
   TF1* fit = new TF1("fit", "pol2");  
   g->Fit("fit"); // pour afficher: options->FitParameters
    // plot pour chercher le minimum de likelihood vs k, i.e. le meilleur k.
   
   double_t p0=fit->GetParameter(0); //c
   double_t p1=fit->GetParameter(1); //b
   double_t p2=fit->GetParameter(2); //a
   double_t chi2=fit->GetChisquare(); //chi2
   double_t ndl=fit->GetNDF(); // Number of Degrees of freedom
   
   // Le minimum d'une parabole es défini comme V = -b/(2a);
   double_t Kmin = -p1/(2*p2);
   
   
   std::cout<<"Valeur de k trouvée : "<<Kmin<<" notre k : "<<k<<std::endl;
   
    
} 
    
   
