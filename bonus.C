{
// Ce fichier contient une macro qui génère un lot de N pseudo-données, la minimisation du likelihood correspondante ainsi que l'histogramme des pulls donnés en consigne.

//#include <iostream>

void generation(int N)
{

    /*
    après l'initialisation des variables :
    première boucle : génération de différents "signaux + bruit" tous associés à un snr choisi aléatoirement, il y en a 10
    partie génération des fonctions : c'est un copié collé, on reste sur le même type de signal
    génération des histogrammes : sur une boucle aussi : comme on a 10 signaux différents il faut que chacunes de leurs caractéristiques soient rangées dans un histogramme
                                  (au final c'est un peu comme si on a des histogrammes rangés dans un repère de 3D)
    deuxième boucle : la ça commence à puer du cul, j'ai pas trop réussi à afficher ce que je voulais, mais j'ai quand même modifié les objets puisqu'on utilise des tableaux mtn
    ensuite il y'a une partie tout en commentaire, je n'ai pas encore traité cette partie
    il y'a une dernière cellule de commentaire "création de l'histogramme des pulls" que j'ai voulu commencer à coder pour anticiper les objets que j'allais utiliser,
    elle est pas fonctionnelle mais si le reste en haut marche, c'est un peu comme ça que je pense les choses, sachant que forcément y'a des erreurs de type et de dimensions
    vu que j'ai pas pu vrmt la tester
                                  
    dernière chose : j'ai pensé à pleins de truc en terme de structures de fichier, de structure de la fonction mais vu que j'ai pas trop de connaissance en c++ et de root,
    c'est un peu l'un des seuls moyens que j'ai eu pour faire un truc un minimum cohérent et testable
    
    BISOUS 
    */
    
    
    // initialisation des variables 
    // gros copié collé, commentaire avec les nouvelles variables introduites
    
    gRandom = new TRandom3(2);
    
    double b_slope = 1.0;
    double moy = 3.0;
    double sig = 0.2;
    double amp = 1.0;
    int bins = 100; 
    double widthbin;
        
    double ki;
    int bin;
    
    double x[10];
    double val_att[10];
    double val_mesure[10];
    
    double poisson[10];
    double tot_poisson[10];

    double likelihood[250];
    double K[250];
    int i = 0;

    double_t p0;
    double_t p1;
    double_t p2;
    double_t chi2;
    double_t ndl;
   
    double_t kmin;
    double_t ksig;
    double_t sigma;
    
    double snr; // pour le calcul du snr
    double ns; // nombre de données du signal
    double nb; // nombre de données du bruit
    
    double snr_tab[10]; // tableau pour stocker les différents snr
    double ns_tab[10]; // tableau pour stocker les différents ns
    double nb_tab[10]; // tableau pour stocker les différents nb
    

    for (int i=0; i<10; i++) // ici le i correspond au nombre de signaux qu'on veut, on aurait pu faire une fonction mais les macros ça veut que 1 fonction c'est chiant
    {
        gRandom->SetSeed();
        snr = gRandom->Rndm(); // génération d'un snr random
        
        ns = snr*N; // données crspdt au signal   
        nb = N-ns; // données crspdt au bruit
        
    snr_tab[i] = snr;
    ns_tab[i] = ns;
    std::cout<<snr_tab[i]<<std::endl; // utile pour vérifier la cohérence des plots
    }
    
        
    // génération des fonctions
        
    TF1* fb = new TF1("background","[0]*exp(-[1]*x)",0,10);
    fb->SetParameters(1,b_slope);
	
    TF1* fs = new TF1("signal", "gaus", 0,10);
    fs->SetParameters(amp,moy,sig);
        
    // génération des histogrammes :
    
    for (int j=0; j<10; j++) // rappel : on a 10 signaux
    {
        TH1F* hb = new TH1F("histb","histb",bins,0,10);
        hb->FillRandom("background",nb_tab[j]);

    
        TH1F* hs = new TH1F("hists","hists",bins,0,10);
	    hs->FillRandom("signal",ns_tab[j]);


        TH1F* htot = new TH1F("histot","histot",bins,0,10);
        htot->Add(hs,hb);
        widthbin = htot->GetBinWidth(1);
        
        
     
    /* *************ATTENTION*************
     ici y'a 30 plots qui s'affiche si on décommente la cellule d'en dessous, donc ça prend un peu de temps
     le plus simple de c'est de plot les 10 plot qui représentent les différents htot et de vérifier leur pertinence avec la ligne 
     std::cout<<snr_tab[i]<<std::endl; , dans la première boucle for
     en effet : si fait generation(10), on aura très peu de données traitées et donc les plots ne sont pas évidents à analyser
    */
    
    
    /*   
    htot->SetTitle("Signal dans le bruit");
    htot->GetXaxis()->SetTitle("x");
    htot->GetYaxis()->SetTitle("Nombre d'evenements");
        
	TCanvas* c = new TCanvas;
    c->cd(1);
    c->SetLogy(1);
    htot->Draw("e");

    
    TFile fout("gen.root","recreate");
	htot->Write();
    fout.Close();
	
	TCanvas* a = new TCanvas;
    a->cd(1);
    hb->Draw("e");
    	
	TCanvas* b = new TCanvas;
    b->cd(1);
    hs->Draw("e");
	*/
	
	TF1* f = new TF1("fontion f","[0]*([1]*gaus(2)+(1-[1])*[5]*exp(-[6]*x))",0,10);
    
    for(ki=0.0775;ki<=0.11;ki=ki+0.0005)
    {
        f->SetParameters(N*widthbin,ki,1/(sqrt(2*3.1415)*sig),moy,sig,1,b_slope);
        tot_poisson[j]=0.;
    }
        for (bin = 0; bin<=bins-1; bin++)
        {
            for (int j=0; j<10; j++)
            {
                x[j] = htot->GetBinCenter(bin+1);
                
                val_att[j] = f->Eval(x[j]);
                val_mesure[j] = htot->GetBinContent(bin+1);
                poisson[j] = TMath::Poisson(val_mesure[j],val_att[j]);
                poisson[j] = TMath::Log(poisson[j]);
                tot_poisson[j] += poisson[j];
            }
        }
            
        likelihood[i] = - tot_poisson[i];
        K[i]=ki;
        i = i+1;
        
    TCanvas* c1 = new TCanvas;
    c1->cd(1);
	
    }
    
    /* *************ATTENTION*************
    pas encore traité
    */
    
    // recherche du maximum de vraisemblance
    

    /*TGraph* g = new TGraph(i,K,likelihood); // ici on utilise i, parce que c'est le nombre de points qu'on a.

    g->SetTitle("Likelihood;Signal sur bruit (k);Likelihood");
    g->Draw("AC*");   //g->Draw("AC*");
    TF1* fit = new TF1("fit", "pol2");  
    g->Fit("fit"); // pour afficher: options->FitParameters
    
   
    // plot pour chercher le minimum de likelihood vs k, i.e. le meilleur k.
   
    p0=fit->GetParameter(0); //c
    p1=fit->GetParameter(1); //b
    p2=fit->GetParameter(2); //a
    chi2=fit->GetChisquare(); //chi2
    ndl=fit->GetNDF(); // Number of Degrees of freedom
   
   // Le minimum d'une parabole es défini comme V = -b/(2a);
    kmin = -p1/(2*p2);

    std::cout<<kmin<<std::endl;
    
    ksig = (-p1+TMath::Sqrt(2*p2))/(2*p2);
   
   //std::cout<<xsig<<std::endl;
   sigma = TMath::Abs(ksig-kmin);
   std::cout<<"Valeur de k trouvée : "<<kmin<<"+-"<<sigma<<" notre k : "<<k<<std::endl;
  
  return kmin;
}
*/

// *************FIN DU ATTENTION*************



// création de l'histogramme des pulls, ***traduire le bon type des variables***

/*
double ns_exp = kmin * N // normalement on traite des histogrammes
double pull = (ns_exp - ns) / sigma

TH1F* h_pull = new TH1F("histp","histp",i,0,10);
h_pull->FillRandom("pulls", pull)

TCanvas* c2 = new TCanvas;
c2->cd(1);
h_pull->Draw("e")

*/
}
}
