/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/

#define prop_num 10
      
#define MASSES_INFO      
  /* Display information about mass spectrum  */

#define HIGGSBOUNDS "../Packages/HiggsBounds-4.2.0"                         // Initially undefined
#define HIGGSSIGNALS "../Packages/HiggsSignals-1.3.0"                       // Initially undefined


#define OMEGA            
  /* Calculate relic density and display contribution of  individual channels */
#define INDIRECT_DETECTION  
  /* Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation; 
     Calculate <sigma*v>;
     Integrate gamma signal over DM galactic squared density for given line 
     of sight; 
     Calculate galactic propagation of positrons and antiprotons.      
  */
      
#define RESET_FORMFACTORS                                       // Initially undefined
  /* Modify default nucleus form factors, 
    DM velocity distribution,
    A-dependence of Fermi-dencity
  */     
#define CDM_NUCLEON     
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */  

#define CDM_NUCLEUS                                              // Initially undefined  
  /* Calculate number of events for 1kg*day and recoil energy distibution
      for various nuclei
  */

#define NEUTRINO

#define DECAYS                                                 // Initially undefined
  
/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS*/
     /* Display  graphical plots on the screen */ 
#define CLEAN
/*===== End of DEFINE  settings ===== */


#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;
   
  ForceUG=0;  // to Force Unitary Gauge assign 1 
  
  if(argc==1)
  { 
      printf(" Correct usage:  ./main  <file with parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }
  
/*********************************************************************************************
Reading input file for properties
*********************************************************************************************/  
  int prop_bool[prop_num];
  int i = 0;
  
  
  FILE* prop = fopen(argv[2],"r");
  
  for(i=0; i<prop_num; i++)
  {
  	fscanf(prop,"%d",&prop_bool[i]);
  }
  	
//////////////////////////////////////////////////////////////////////////////////////////////
                    
                    
   printf(",\"");

                               
  err=readVar(argv[1]);
  
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}
           
  

  err=sortOddParticles(cdmName);



  if(err) { printf("Can't calculate %s\n",cdmName); printf("\"\n"); return 1;}
  
   qNumbers(cdmName, &spin2, &charge3, &cdim);
   
   printf("\""); 
  printf(",\"");
  
  
   printf("\nDark matter candidate is '%s' with spin=%d/2\n",
    cdmName,       spin2); 
   if(charge3) { printf("Dark Matter has electric charge %d*3\n",charge3); printf("\""); exit(1);}
   if(cdim!=1) { printf("Dark Matter ia a color particle\n"); printf("\""); exit(1);}
   
   printf("\"");
   
//#ifdef MASSES_INFO

if(prop_bool[0])
{
	
  printf(",\"");	
  printf("\n=== MASSES OF HIGG AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
  printf("\"");
  
}
//#endif


//#ifdef HIGGSBOUNDS

if(prop_bool[1])
{
   printf(",\"");
   if(access(HIGGSBOUNDS "/HiggsBounds",X_OK )) system( "cd " HIGGSBOUNDS "; ./configure; make ");
   HBblocks("HB.in");
   system(HIGGSBOUNDS "/HiggsBounds  LandH SLHA 1 0 HB.in HB.out > hb.stdout");
   slhaRead("HB.out",1+4);
    printf("HB result= %.0E  obsratio=%.2E\n",slhaValFormat("HiggsBoundsResults",0.,"1 2 %lf"), slhaValFormat("HiggsBoundsResults",0.,"1 3 %lf" )  );
   { char hbInfo[100];
    if(0==slhaSTRFormat("HiggsBoundsResults","1 5 ||%[^|]||",hbInfo)) printf("Channel: %s\n",hbInfo);
   }    
   printf("\"");
}
//#endif

//#ifdef HIGGSSIGNALS

if(prop_bool[2])
{
printf(",\"");

#define DataSet " latestresults "
//#define Method  " peak " 
//#define  Method " mass "
#define  Method " both "
#define PDF  " 2 "  // Gaussian
//#define PDF " 1 "  // box 
//#define PDF " 3 "  // box+Gaussia
#define dMh " 2 "
   printf("HiggsSignals:\n");
   if(access(HIGGSSIGNALS "/HiggsSignals",X_OK )) system( "cd " HIGGSSIGNALS "; ./configure; make ");
     system("rm -f HS.in HS.out");
     HBblocks("HS.in");
     system(HIGGSSIGNALS "/HiggsSignals" DataSet Method  PDF  " SLHA 1 0 HS.in > hs.stdout");
     system("grep -A 10000  HiggsSignalsResults HS.in > HS.out");
     slhaRead("HS.out",1+4);
     printf("  Number of observables %.0f\n",slhaVal("HiggsSignalsResults",0.,1,7));
     printf("  total chi^2= %.1E\n",slhaVal("HiggsSignalsResults",0.,1,12));
     printf("  HS p-value = %.1E\n", slhaVal("HiggsSignalsResults",0.,1,13));     
#undef dMh
#undef PDF
#undef Method
#undef DataSet


printf("\"");

}
//#endif                      // End of HIGGSSIGNALS


//#ifdef OMEGA

if(prop_bool[3])
{ int fast=1;
  double Beps=1.E-5, cut=0.01;
  double Omega,Xf;   
  int i;
  
  printf(",\"");
  
// to exclude processes with virtual W/Z in DM   annihilation      
    VZdecay=0; VWdecay=0; cleanDecayTable(); 

// to include processes with virtual W/Z  also  in co-annihilation 
//   VZdecay=2; VWdecay=2; cleanDecayTable(); 
  
  printf("\n==== Calculation of relic density =====\n");  
  Omega=darkOmega(&Xf,fast,Beps);
  printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
  printChannels(Xf,cut,Beps,1,stdout);   

//   VZdecay=1; VWdecay=1; cleanDecayTable();  // restore default

printf("\"");

}
//#endif


//#ifdef INDIRECT_DETECTION

if(prop_bool[4])
{ 
  
  printf(",\"");
  int err,i;
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum(4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */
  printf("sigmav=%.2E[cm^3/s]\n",sigmaV);  


  if(SpA)
  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);     
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{1} at f=%.2f[rad], cone angle %.2f[rad]",fi,2*dfi);
     displaySpectrum(FluxA,txt,Emin,Mcdm,1);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);       
  }

  if(SpE)
  { 
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
#ifdef SHOWPLOTS     
    displaySpectrum(FluxE,"positron flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm, 1);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
  if(SpP)
  { 
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
#ifdef SHOWPLOTS    
     displaySpectrum(FluxP,"antiproton flux [cm^2 s sr GeV]^{-1}" ,Emin, Mcdm,1);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
  }
  
  printf("\"");
  
}  
//#endif

//#ifdef RESET_FORMFACTORS

if(prop_bool[5])
{
/* 
   The user has approach to form factors  which specifies quark contents 
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma"; 
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigma0[MeV])  
   calculates and rewrites Scalar form factors
*/
  printf(",\"");
  
  printf("\n======== RESET_FORMFACTORS ======\n");
 
  printf("protonFF (default) d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//  To restore default form factors of  version 2  call 
//  calcScalarQuarkFF(0.553,18.9,55.,243.5);
 
  calcScalarFF(0.553,18.9,70.,35.);

  printf("protonFF (new)     d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

  printf("\"");

}
//#endif

//#ifdef CDM_NUCLEON

if(prop_bool[6])
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        
  int i;
  
  printf(",\"");
  
printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   
    nucleonAmplitudes(CDM1,NULL, pA0,pA5,nA0,nA5);
    printf("CDM[antiCDM]-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM[antiCDM]-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*pA0[0]*pA0[0],SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[0]*pA5[0],3*SCcoeff*pA5[1]*pA5[1]);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*nA0[0]*nA0[0],SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[0]*nA5[0],3*SCcoeff*nA5[1]*nA5[1]);
       
    printf("\"");   

}
//#endif
  
//#ifdef CDM_NUCLEUS

if(prop_bool[7])
{ double dNdE[300];
  double nEvents;
  
  printf(",\"");

printf("\n======== Direct Detection ========\n");    

  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,NULL,dNdE);

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
                                                                                                         
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 73Ge",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,NULL,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 131Xe",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,SxxNa23,NULL,dNdE);

  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 23Na",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,SxxI127,NULL,dNdE);

  printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 127I",0,199);
#endif

printf("\"");
  
}
//#endif 


//#ifdef NEUTRINO

if(prop_bool[8])
{ double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=0.01;
  
  printf(",\"");

 printf("\n===============Neutrino Telescope=======  for  ");
 if(forSun) printf("Sun\n"); else printf("Earth\n");

  err=neutrinoFlux(Maxwell,forSun, nu,nu_bar);
#ifdef SHOWPLOTS
  displaySpectrum(nu,"nu flux from Sun [1/Year/km^2/GeV]",Emin,Mcdm,1);
  displaySpectrum(nu_bar,"nu-bar from Sun [1/Year/km^2/GeV]",Emin,Mcdm,1);
#endif
{ double Ntot;
  double Emin=10; //GeV
  spectrInfo(Emin/Mcdm,nu, &Ntot,NULL);
    printf(" E>%.1E GeV neutrino flux       %.3E [1/Year/km^2] \n",Emin,Ntot);
  spectrInfo(Emin/Mcdm,nu_bar, &Ntot,NULL);
    printf(" E>%.1E GeV anti-neutrino flux  %.3E [1/Year/km^2]\n",Emin,Ntot);
}

/* Upward events */

  muonUpward(nu,nu_bar, mu);
#ifdef SHOWPLOTS  
  displaySpectrum(mu,"Upward muons[1/Year/km^2/GeV]",1,Mcdm/2,1);
#endif
  { double Ntot;
    double Emin=1; //GeV
    spectrInfo(Emin/Mcdm,mu, &Ntot,NULL);
    printf(" E>%.1E GeV Upward muon flux    %.3E [1/Year/km^2]\n",Emin,Ntot);
  }

/* Contained events */
  muonContained(nu,nu_bar,1., mu);
#ifdef SHOWPLOTS  
  displaySpectrum(mu,"Contained  muons[1/Year/km^3/GeV]",Emin,Mcdm,1);
#endif
  { double Ntot;
    double Emin=1; //GeV
    spectrInfo(Emin/Mcdm,mu, &Ntot,NULL);
    printf(" E>%.1E GeV Contained muon flux %.3E [1/Year/km^3]\n",Emin,Ntot);
  }
  
  printf("\"");
  
}
//#endif

//#ifdef DECAYS

if(prop_bool[9])
{ 

  printf(",\""); 
  txtList L;
   double width,br;
   char * pname;
   printf("\n================= Decays ==============\n");
      
   pname = "H";
   width=pWidth(pname,&L);
   printf("\n%s :   total width=%.2E \n and Branchings:\n",pname,width);
   printTxtList(L,stdout);

   pname = "~W+";
   width=pWidth(pname,&L);
   printf("\n%s :   total width=%.2E \n and Branchings:\n",pname,width);
   printTxtList(L,stdout);
   
   printf("\"");
}
//#endif

#ifdef CLEAN 
    system("rm -f HB.in HB.out HS.in HS.out hb.stdout hs.stdout debug_channels.txt debug_predratio.txt Key.dat");
    killPlots();   
#endif 
  return 0;
}
