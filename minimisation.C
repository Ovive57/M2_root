{
    /*************** Initialisation des variables ***************/
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
    
    double tot_poisson; // Somme des logaritmes de toutes les fonctions poissons de chaque bin

    
    // Tableaux pour le plot
    double likelihood[250]; //likelihood des poissons, axe y. Il faudrait voir la taille, même taille qu'iterations sur k.
    double K[250]; // l'axe x.
    int i = 0; // Compteur pour remplir chaque tableau pour le plot.
    
    // Parametres pour le fit du likelihood (y = a*x*x + b*x + c)
    double_t p0; //c
    double_t p1; //b
    double_t p2;//a
    double_t chi2; //chi2
    double_t ndl; // Number of Degrees of freedom
   
    double_t kmin; // K qui correspond au minimum de l'ajustement
    double_t ksig; // K qui correspond au sigma de l'ajustement
    double_t sigma; // sigma de l'ajustement
    
    
    /*************** Création pseudo-données ***************/
      
    /* Fonction du bruit */ 
	TF1* fb = new TF1("background","[0]*exp(-[1]*x)",0,10);
	fb->SetParameters(1,b_slope); // Amplitud exp, pente exp
	
	/* Fonction du signal */ 
	TF1* fs = new TF1("signal", "gaus", 0,10);
	fs->SetParameters(amp,moy,sig);
	
	/* Histogramme du bruit */
	TH1F* hb = new TH1F("histb","histb",bins,0,10);
	hb->FillRandom("background",nb);
	hb->SetLineColor(kBlue);
	
	/* Histogramme du signal */
	TH1F* hs = new TH1F("hists","hists",bins,0,10);
	hs->FillRandom("signal",ns);
	hs->SetLineColor(kBlue);
    
    /* Les deux histogrammes ensemble */
    TH1F* htot = new TH1F("histot","histot",bins,0,10);
    htot->Add(hs,hb);
        // On calcule la largeur des bins pour plus tard (les boucles) tous les bins sont pareils donc on choisi par exemple de mesurer la largeur du premier :
    widthbin = htot->GetBinWidth(1); 
    
    
    /*************** Plot pseudo-données avec les valeurs attendues ***************/
    
        // Création de Canvas pour faire le plot.
	TCanvas* c = new TCanvas;
	c->cd(1);
	c->SetLogy(1);
	
	    // Paramètres plot.
    htot->SetTitle("Histogramme d'evenements");
    htot->GetXaxis()->SetTitle("x"); // Pas sûre que ça soit temps xd
    htot->GetYaxis()->SetTitle("Nombre d'evenements");

	/* Plot des données */ 
	htot->Draw("e"); 

        // On souvegarde l'histogramme dans un fichier .root :
	TFile fout("gen.root","recreate");
	htot->Write();
	fout.Close();

 
    /* Fonction pour les valeurs attendues */
    TF1* ffit = new TF1("fontion f","[0]*([1]*gaus(2)+(1-[1])*[5]*exp(-[6]*x))",0,10);
    
        // Paramètres plot.
    ffit->SetTitle("Valeurs attendues");
    ffit->SetParameters(ntot*widthbin,k,1/(sqrt(2*3.1415)*sig),moy,sig,1,b_slope);
    
    /* Plot des valeurs attendues avec les données */
    ffit->Draw("same");
        // Pour ajouter une legende :
    c->BuildLegend();


    /*************** Méthode maximum de vraisemblance pour chercher le paramètre k, signal sur bruit ***************/

    /* Réinisalisation fonction pour les valeurs attendues */
    TF1* f = new TF1("fontion f","[0]*([1]*gaus(2)+(1-[1])*[5]*exp(-[6]*x))",0,10);
    f->SetParameters(ntot*widthbin,k,1/(sqrt(2*3.1415)*sig),moy,sig,1,b_slope);
    
    /*Boucle sur plusieurs k et tous les bins : On calcule la probabilité de poisson pour chaque bin et après on calcule la probabilité totale pour chaque k */
    for(ki=0.0775;ki<=0.11;ki=ki+0.0005){
        // On réinisialize le fit avec les diffèrents k.
        f->SetParameter(1,ki);
        // On réinisialize la somme des poissons à chaque fois qu'on fait un boucle sur k.
        tot_poisson=0.; 
        
        for (bin = 0; bin<=bins-1; bin++){
            // On cherche le x de chaque bin : 
            x = htot->GetBinCenter(bin+1);
            // On trouve la valeur attendue en chaque x : 
            val_att = f->Eval(x);

            // On trouve la valeur mesurée en prennant la valeur du bin : 
            val_mesure = htot->GetBinContent(bin+1); 
	         
	        // On calcule avec Poisson la probabilité de trouver la valeur attendue en ayant la valeur mesurée P(mesurée,attendue) et on fait la somme des logarithmes:
            tot_poisson += TMath::Log(TMath::Poisson(val_mesure,val_att));
            }

    // On sauvegarde les probabilités obtenues dans des tableaux
    likelihood[i] = - tot_poisson;
    K[i]=ki;
    // On actualise le compteur:
    i = i+1;
        }
    
    /*************** Plot vraisemblance ***************/
    
        // Création de Canvas pour faire le plot.
    TCanvas* c1 = new TCanvas;
    c1->cd(1);

    /* Plot avec TGraph */
    TGraph* g = new TGraph(i,K,likelihood); 
    	
    	// Paramètres plot.
    g->SetTitle("Likelihood;Signal sur bruit (k);Likelihood");
    g->Draw("AC*");
        
        // Fit avec un polynôme de second degré.
    TF1* fit = new TF1("fit", "pol2");  
    g->Fit("fit");

   
    /*************** Calcul du minimum, et donc de k ***************/
   
    /* Paramètres du fit y = a*x*x + b*x + c */
    p0=fit->GetParameter(0); //c
    p1=fit->GetParameter(1); //b
    p2=fit->GetParameter(2); //a
    chi2=fit->GetChisquare(); //chi2
    ndl=fit->GetNDF(); // Number of Degrees of freedom
   
        // Le minimum d'une parabole es défini comme V = -b/(2a);
    kmin = -p1/(2*p2);

    /* On cherche un sigma autour de notre k, +- 1/2 */
        // ymin = a*kmin*kmin + b*kmin + c; ymin + 0.5 = a*x*x + b*x + c ----> x = -b/2a +- sqrt(2a)/2a
    ksig = (-p1+TMath::Sqrt(2*p2))/(2*p2);

        // sigma = |x-kmin| 
    sigma = TMath::Abs(ksig-kmin);
    std::cout<<"Valeur de k trouvée : "<<kmin<<"+/-"<<sigma<<" notre k : "<<k<<std::endl;    
}  
