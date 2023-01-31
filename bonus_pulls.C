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
   
    float kmin; // K qui correspond au minimum de l'ajustement
    float ksig; // K qui correspond au sigma de l'ajustement
    float sigma; // sigma de l'ajustement
    
    // Parametres boucle et histogramme des pulls
    int j; // Compteur des boucles
    int N = 10000; // Nombre de boucles
    float pulls;

    
    /*************** Initialisation des fonctions ***************/
     
    /* Fonction du bruit */ 
	TF1* fb = new TF1("background","[0]*exp(-[1]*x)",0,10);
	fb->SetParameters(1,b_slope); // Amplitud exp, pente exp
	/* Histogramme du bruit */
    TH1F* hb = new TH1F("histb","histb",bins,0,10);
	
	/* Fonction du signal */ 
	TF1* fs = new TF1("signal", "gaus", 0,10);
	fs->SetParameters(amp,moy,sig);
	/* Histogramme du signal */	
    TH1F* hs = new TH1F("hists","hists",bins,0,10);
    
    /* Les deux histogrammes ensemble */
    TH1F* htot = new TH1F("histot","histot",bins,0,10);
	
	/* Fonction pour les valeurs attendues */   
    TF1* f = new TF1("fontion f","[0]*([1]*gaus(2)+(1-[1])*[5]*exp(-[6]*x))",0,10);
     
    /* Fonction pour le fit */   
    TF1* fit = new TF1("fit", "pol2"); 
    
	/* Histogramme des pulls */	
    TH1F* histpull = new TH1F("pulls","pulls",bins,-5,5);
    

    for(j=0;j<N;j++){
        /*************** Création pseudo-données ***************/
	        // Réinitialisation du compteur
	    i = 0; 
	    
	    /* Remplissage histogramme du bruit */
	    hb->FillRandom("background",nb);
	    hb->SetLineColor(kBlue);
	
	    /* Remplissage histogramme du signal */
	    hs->FillRandom("signal",ns);
	    hs->SetLineColor(kBlue);
    
        /* Les deux histogrammes ensemble */
        htot->Add(hs,hb);
        widthbin = htot->GetBinWidth(1); 
    


        /*************** Méthode maximum de vraisemblance pour chercher le paramètre k, signal sur bruit ***************/

        /* Inisalisation fonction pour les valeurs attendues */
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
    
        /*************** Fit vraisemblance ***************/
        TGraph* g = new TGraph(i,K,likelihood);  
       
            // Fit avec un polynôme de second degré.
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
        
            // Les résidus :
        pulls = (kmin-k)/sigma;
        std::cout<<"pull : "<<pulls<<std::endl;  
        
            // On réinitialise les histogrammes pour la prochaine boucle
        hb->Reset();
        hs->Reset();
        htot->Reset();
        
            // On remplit l'histogramme des pulls 
        histpull->Fill(pulls);
           
    
   }
   
    /*************** Plot histogramme des résidus ***************/

    	// Paramètres plot.
    histpull->SetTitle("Histogramme des residus;Pulls;Nombre");
    histpull->Draw();
        // Légende et fit
    gStyle->SetOptStat(0);
    histpull->Fit("gaus");
    gStyle->SetOptFit(0012); 
}  
