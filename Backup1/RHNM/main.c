/*====== Modules ===============
   Keys to switch on
  various modules of micrOMEGAs
================================*/
      
#define MASSES_INFO      
#define OMEGA            
#define INDIRECT_DETECTION  
#define RESET_FORMFACTORS
#define CDM_NUCLEON     
#define CDM_NUCLEUS
#define NEUTRINO 
#define DECAYS
#define CROSS_SECTIONS

/*===== end of Modules  ======*/

/*===== Options ========*/
/* #define SHOWPLOTS */
     /* Display  graphical plots on the screen */ 

/*===== End of DEFINE  settings ===== */

#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;
   
   char s1,s2[3];
   s1='"';
   s2[0]=',';
   s2[1]='"';

  ForceUG=0;  /* to Force Unitary Gauge assign 1 */
  
  if(argc==1)
  { 
      printf(" Correct usage:  ./main <file with parameters> \n");
      exit(1);
  }
  
/*********************************************************************************************
Reading input file for properties
*********************************************************************************************/  
  int prop_num = 9,prop_bool[prop_num];
  int i = 0;
  
  
  FILE* prop = fopen(argv[2],"r");
  
  for(i=0; i<prop_num; i++)
  {
  	fscanf(prop,"%d",&prop_bool[i]);
  }
  	
//////////////////////////////////////////////////////////////////////////////////////////////
                               
printf("%s",s2);

/*  err=readVar(argv[1]);*/
   err=readVarRHNM(argv[1]);
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}
  

  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); printf("%c\n",s1); return 1;}
  qNumbers(cdmName, &spin2, &charge3, &cdim);
  
  printf("%c",s1);
  
   printf("%s",s2);
  
  printf("\nDark matter candidate is '%s' with spin=%d/2 \n",cdmName,spin2); 
  if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); printf("%c",s1); exit(1);}
  if(cdim!=1) { printf("Dark Matter is a color particle\n"); printf("%c",s1); exit(1);}
  if(strcmp(cdmName,"~n4")) printf(" ~n4 is not CDM\n"); 
  
  printf("%c",s1);
                               
if(prop_bool[0])      // MASSES_INFO
{
  printf("%s",s2);
  printf("\n=== MASSES OF PARTICLES OF ODD SECTOR: ===\n");
  printMasses(stdout,1);
  
  printf("%c",s1);
}


if(prop_bool[1])       //OMEGA
{ 
  printf("%s",s2);
  int fast=1;
  double Beps=1.E-5, cut=0.01;
  double Omega,Xf; 

// to exclude processes with virtual W/Z in DM   annihilation      
    VZdecay=0; VWdecay=0; cleanDecayTable(); 

// to include processes with virtual W/Z  also  in co-annihilation 
//   VZdecay=2; VWdecay=2; cleanDecayTable(); 
    
  printf("\n==== Calculation of relic density =====\n");  



  Omega=darkOmega(&Xf,fast,Beps);
  //printf("Xf=%.2e Omega=%.2e\n",Xf,Omega);
  printChannels(Xf,cut,Beps,1,stdout);
  
  VZdecay=1; VWdecay=1; cleanDecayTable();  // restore default
  
  printf("%c",s1);
  
  printf(",%lf,%lf",Xf,Omega);
}



if(prop_bool[2])         //INDIRECT_DETECTION
{ 
  printf("%s",s2);
  int err,i;
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum(2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
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
    posiFluxTab(Emin, sigmaV, SpE, FluxE);
#ifdef SHOWPLOTS     
    displaySpectrum(FluxE,"positron flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm,1);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
  if(SpP)
  { 
    pbarFluxTab(Emin, sigmaV, SpP, FluxP  ); 
#ifdef SHOWPLOTS    
     displaySpectrum(FluxP,"antiproton flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm,1);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
  }
  
printf("%c",s1);
}  


if(prop_bool[3])          //RESET_FORMFACTORS
{

printf("%s",s2);
/* 
   The user has approach to form factors  which specifies quark contents 
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma"; 
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigmaS[MeV])  
   calculates and rewrites Scalar form factors
*/

  printf("protonFF (default) d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

  calcScalarQuarkFF(0.46,27.5,34.,42.);

//  To restore default form factors of  version 2  call 
//  calcScalarQuarkFF(0.553,18.9,55.,243.5);
 
  printf("protonFF (new)     d %E, u %E, s %E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %E, u %E, s %E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

printf("%c",s1);

}


if(prop_bool[4])            //CDM_NUCLEON
{ 
  printf("%s",s2);

  double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        

printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

    nucleonAmplitudes(CDM1, FeScLoop, pA0,pA5,nA0,nA5);
    printf("CDM[antiCDM]-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM[antiCDM]-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*pA0[0]*pA0[0],SCcoeff*pA0[1]*pA0[1],3*SCcoeff*pA5[0]*pA5[0],3*SCcoeff*pA5[1]*pA5[1]);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n",
       SCcoeff*nA0[0]*nA0[0],SCcoeff*nA0[1]*nA0[1],3*SCcoeff*nA5[0]*nA5[0],3*SCcoeff*nA5[1]*nA5[1]);
printf("%c",s1);
}

  
if(prop_bool[5])              // CDM_NUCLEUS
{ 
  printf("%s",s2);

  double dNdE[300];
  double nEvents;

printf("\n======== Direct Detection ========\n");    

  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,FeScLoop,dNdE);

  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));
                                                                                                         
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 73Ge",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,FeScLoop,dNdE);

  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 131Xe",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,SxxNa23,FeScLoop,dNdE);

  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 23Na",0,199);
#endif

  nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,SxxI127,FeScLoop,dNdE);

  printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",
                                   cutRecoilResult(dNdE,10,50));                                   
#ifdef SHOWPLOTS
    displayRecoilPlot(dNdE,"Distribution of recoil energy of 127I",0,199);
#endif
  
printf("%c",s1);
}



if(prop_bool[6])  NEUTRINO
{ 
  printf("%s",s2);

  double nu[NZ], nu_bar[NZ],mu[NZ];
  double Ntot;
  int forSun=1;
  double Emin=0.01;

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
  displaySpectrum(mu,"Upward muons[1/Year/km^2/GeV]",1,Mcdm/2);
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
  printf("%c",s1);
}


if(prop_bool[7])             //DECAYS
{  
   printf("%s",s2);
   txtList L;
   double width,br;
   char * pname;
   
    if(!VZdecay || !VWdecay  ){ cleanDecayTable(); VZdecay=1; VWdecay=1;}

   printf("\n Calculation of particle decays\n");
   pname = "H";
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%.3E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);

   pname = "Zp";
    width=pWidth(pname,&L);
    printf("\n%s :   total width=%.3E \n and Branchings:\n",pname,width);
    printTxtList(L,stdout);
    
    printf("%c",s1);  
}


if(prop_bool[8])        //CROSS_SECTIONS
{
  printf("%s",s2);
  double Pcm=500, cosmin=-0.99, cosmax=0.99, cs;
  numout* cc;
printf("\n====== Calculation of cross section ====\n");  

printf(" e^+, e^- annihilation\n");
  Pcm=500.;
  Helicity[0]=0.5;    /* helicity : spin projection on direction of motion   */    
  Helicity[1]=-0.5;   /* helicities ={ 0.5, -0.5} corresponds to vector state */
  printf("Process e,E->2*x at Pcm=%.3E GeV\n",Pcm);
  cc=newProcess("e%,E%->2*x");
  if(cc)
  { int ntot,l;
    char * name[4];
    procInfo1(cc,&ntot,NULL,NULL);
    for(l=1;l<=ntot; l++)
    { int err;
      double cs;
      char txt[100];
      procInfo2(cc,l,name,NULL);
      sprintf(txt,"%3s,%3s -> %3s %3s  ",name[0],name[1],name[2],name[3]);
      cs= cs22(cc,l,Pcm,cosmin,cosmax,&err);
      if(err) printf("%-20.20s    Error\n",txt);
      else if(cs) printf("%-20.20s  %.2E [pb]\n",txt,cs); 
    }
  } 
printf("%c",s1);
}


  killPlots();
  
  printf("\n");
  
  return 0;
}
