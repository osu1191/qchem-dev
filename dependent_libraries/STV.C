#ifndef lint
static char vcid[] = "$Id: STV.C,v 1.32 2010/07/09 16:42:50 herbert Exp $";
#endif /* lint */

#include <stdio.h>
#include "qchem.h"
#include "STV.hh"
#include "BasisSet.hh"
#include "MultipoleField.h"
#include "PointCharges.h"
#include "efpman2/EFP2.h"
#include "rem_values.h"
#include <libscrf/pcm/pcm_class.h>
#include <libpeqs/peqs_class.h>
#include "fileman.h"
STV::STV(INTEGER fnum) :
  fNum(fnum),
  bCode(0),
  NOrb(0),
  NB2(0),
  LenX(0),
  MaxMultOrd(0),
  MaxK(0),
  hOffSet(0),
  xOffSet(0),
  sOffSet(0),
  tOffSet(0),
  vOffSet(0),
  eOffSet(0),
  mOffSet(0),
  h(0), 
  x(0), 
  sv(0), 
  sm(0), 
  t(0), 
  v(0), 
  e(0), 
  m(0),
  vemb(0),
  ExtCharges(0){}

STV::~STV()
{
  freeAll();
}
void STV::freeAll()
{
  freeH();
  freeX();
  freeS();
  freeT();
  freeV();
  freeE();
  freeM();
}

void STV::freeH() { if ( h != 0 ) QFree(h); h = 0; }
void STV::freeX() { if ( x != 0 ) QFree(x); x = 0; }
void STV::freeS() { if ( sv != 0 ) { QFree(sv); QFree(sm); } sv = 0; sm = 0; }
void STV::freeT() { if ( t != 0 ) QFree(t); t = 0; }
void STV::freeV() { if ( v != 0 ) QFree(v); v = 0; }
void STV::freeE() { if ( e != 0 ) QFree(e); e = 0; }
void STV::freeM() { if ( m != 0 ) QFree(m); m = 0; }

INTEGER STV::code() const { return bCode; }
INTEGER* STV::ftnCode() { return &bCode; }

INTEGER STV::getNOrb() const { return NOrb; }
INTEGER* STV::ftnGetNOrb() { return &NOrb; }

void STV::setOffSets(INTEGER crntOffSet)
{
  hOffSet = crntOffSet;
  xOffSet = hOffSet + NB2*sizeof(double);
  sOffSet = xOffSet + LenX*sizeof(double);
  tOffSet = sOffSet + (NB2+NBasis*NBasis)*sizeof(double);
  vOffSet = tOffSet + NB2*sizeof(double);
  if (rem_read(REM_PSEUDOPOTENTIAL)>0) {
    eOffSet = vOffSet + NB2*sizeof(double); //place ECP matrix just before m
    mOffSet = eOffSet + NB2*sizeof(double);
  }
  else {
    mOffSet = vOffSet + NB2*sizeof(double);
  }
}

static MultipoleField MField;
extern void compute_mess_e_energy(double* VPC);

void STV::makeH()
{
  if (h == 0) h = QAllocDouble(NB2);

  // Start with T + V

  VRadd(h,getT(),getV(),NB2);

  //Xing - embedding potential from input
  if (rem_read(REM_VEMB_DFET) == 1 ) VRadd(h,h,getVemb(), NB2);

  //RDA - COmbine the pseudopotential matrix, if required.
  if (rem_read(REM_PSEUDOPOTENTIAL) > 0) VRadd(h,h,getE(),NB2);

  // Next add any external multipole field contribution

  LOGICAL toFreeM = m == 0;

  LOGICAL doOnsager = (rem_read(REM_SOLVENT_METHOD) == ONSAGER) ? 1 : 0;
  if(!doOnsager){ //this messes up the E_sol in SCRF  KST 08/07
    for (INTEGER i = 0; i < MField.NComponents(); ++i)
      VRaxpy(h,-MField.Strength(i),getM(MField.Component(i)),h,NB2);
  }

  if ( toFreeM ) freeM();

  // ----------------------------------------------------------- //
  // ----- Next add any external point charge contribution ----- //
  // ----------------------------------------------------------- //

  if (rem_read(REM_SOLVE_PEQ) && rem_read(REM_PEQ_CORR) == 1){
    PEQSclass peqs;
    ShlPrs S2(code());
    INTEGER NB2car = S2.getNB2car();
    double *Vpeqs = QAllocDouble(NB2car);
    double *WSMuNu = QAllocDouble(NB2car);
    VRload(WSMuNu,NB2car,0.0);
    double Ediis=0.0;
    int NGrid=0,MaxPPer=peqs.Batch(),Nz=peqs.NZ();
    FileMan(FM_READ,FILE_PEQS_DATA,FM_INT,1,0,FM_BEG,&NGrid);
    FileMan(FM_READ,FILE_PEQS_DATA,FM_DP,1,0,FM_CUR,&Ediis); // gygi
    double *charges = QAllocDouble(NGrid);
    double *Carts = QAllocDouble(3*NGrid);
    cout << "NGrid = " << NGrid << endl;
    cout << "sccs_switch = " << peqs.switch1() << endl; //Gygi
    cout << "diis error = " << Ediis << endl; //Gygi
    double sccs_switch=peqs.switch1(); //Gygi
    FileMan(FM_READ,FILE_PEQS_DATA,FM_DP,NGrid,0,FM_CUR,charges);
    FileMan(FM_READ,FILE_PEQS_DATA,FM_DP,3*NGrid,0,FM_CUR,Carts);
    if (NGrid <= MaxPPer) MaxPPer = NGrid;
    else if ((MaxPPer % Nz) != 0) MaxPPer = MaxPPer - (MaxPPer % Nz);
    int NPPer = MaxPPer;
    double *tempWSMuNu = QAllocDouble(NB2car+MaxPPer);
    for (int ijk=0;ijk<NGrid;ijk+=MaxPPer){
      if (MaxPPer <= (NGrid-ijk)) NPPer = MaxPPer;
      else NPPer = (NGrid-ijk);
      VRload(tempWSMuNu,NB2car+NPPer,0.0);
      AOints(tempWSMuNu,NULL,NULL,NULL,NULL,&Carts[3*ijk],&charges[ijk],NULL,&NPPer,11,S2,S2);
      VRadd(WSMuNu,WSMuNu,tempWSMuNu,NB2car);
    }
    VRcopy(Vpeqs,WSMuNu,NB2car);
    VRadd(h,h,Vpeqs,NB2);
    if(Ediis < sccs_switch){
      cout << "Sccs will be added now" << endl;
      double *Vsccs = QAllocDouble(NB2car); // Gygi
      double *sc_potential = QAllocDouble(NBasis*NBasis); // Gygi
      FileMan(FM_READ,FILE_PEQS_DATA,FM_DP,NBasis*NBasis,0,FM_CUR,sc_potential); // Gygi
      ScaV2M(sc_potential,Vsccs,true,false); // Gygi
      VRadd(h,h,Vsccs,NB2);  // Gygi
      QFree(Vsccs);QFree(sc_potential); // Gygi
    }
    QFree(Vpeqs);QFree(WSMuNu);QFree(tempWSMuNu);
    QFree(Carts);QFree(charges);
  }

  INTEGER IGaussianBlur = rem_read(REM_GAUSSIAN_BLUR);
  INTEGER IpcmMeth      = rem_read(REM_PCM_TESSEL_METHOD);
  INTEGER Ipcmprint     = rem_read(REM_PCM_PRINT);
  PCMclass *thePCM = NULL;
  if(rem_read(REM_PCM_METHOD) > 0) thePCM = new PCMclass;
  //if(thePCM) cout << " Size of PCM object = " << sizeof(thePCM) << " bytes" << endl;

  // AWL -- PCM external point charges
  if(thePCM){
  if(thePCM->ready_for_Hcore())
  {
    // read PCM charges from PCM_SURFACE_CHARGES.
    int n_charges = 0;
    FileMan(FM_READ,FILE_PCM_SURFACE_CHARGES,FM_INT,1,0,FM_BEG,&n_charges);

    if(n_charges > 0){
      ShlPrs S2(code());
      INTEGER NB2car = S2.getNB2car();
      double* VpcmPC = QAllocDouble(NB2car);
      PointCharges pcmPC;

      if(IGaussianBlur != 1){
        pcmPC.Read(FILE_PCM_SURFACE_CHARGES);
        pcmPC.ansToBohr();
        if(Ipcmprint > 2) pcmPC.Print("PCM Charges being added to the core h");
        PointCharges* pnt_pcmPC = &pcmPC;
#ifdef DEVELOPMENT
        cout << "Updating external PCM charges" << endl;
#endif
        MakeV(VpcmPC,*pnt_pcmPC,S2);

	// Noneq SS-PCM: add the slow charge external potential of the reference state 
        if(thePCM->has_neqSS()){
          VRadd(h,h,VpcmPC,NB2);
          pcmPC.Read(FILE_PCM_FAST_SLOW_CHARGES);
          pcmPC.ansToBohr();
          if(Ipcmprint > 2) pcmPC.Print("PCM Slow Charges being added to the core h");
          MakeV(VpcmPC,*pnt_pcmPC,S2);
	}

      }else{
        if(Ipcmprint > 2){
          pcmPC.Read(FILE_PCM_SURFACE_CHARGES);	
          pcmPC.ansToBohr();
          pcmPC.Print("PCM Charges being added to the core h");
        }
        double* WSMuNu = QAllocDouble(NB2car + n_charges); //extra room for libint
        VRload(WSMuNu, NB2car + n_charges, 0.0);
        #ifdef DEVELOPMENT
           cout << "Evaluate Gaussian charge contribution to one-electron "
                << "Hamiltonian from PCM charges\n";
        #endif
        FileMan(FM_READ,FILE_PCM_SURFACE_CHARGES,FM_INT,1,0,FM_BEG,&n_charges);
        double *charges = QAllocDouble(n_charges);
        FileMan(FM_READ,FILE_PCM_SURFACE_CHARGES,FM_DP,n_charges,0,FM_CUR,charges);

	// Noneq SS-PCM: add the reference slow charge contribution
        if(thePCM->has_neqSS()){
	  // add up the external charges for nonequilbrium 
          FileMan(FM_READ,FILE_PCM_FAST_SLOW_CHARGES,FM_INT,1,0,FM_BEG,&n_charges);
	  double *qslow = QAllocDouble(n_charges);
          FileMan(FM_READ,FILE_PCM_FAST_SLOW_CHARGES,FM_DP,n_charges,0,FM_CUR,qslow);
          if(Ipcmprint > 2){
	    rem_write(0,REM_GAUSSIAN_BLUR);
            pcmPC.Read(FILE_PCM_FAST_SLOW_CHARGES);
            pcmPC.ansToBohr();
            pcmPC.Print("PCM Slow Charges being added to the core h");
	    rem_write(IGaussianBlur,REM_GAUSSIAN_BLUR);
          }
	  VRadd(charges,charges,qslow,n_charges);
	  QFree(qslow);
	}
	// form the reaction potential
        rem_write(1, REM_MODIFY_INT_JOB);
        AOints(WSMuNu,NULL,NULL,NULL,charges,NULL,NULL,NULL,NULL,15);
        rem_write(0, REM_MODIFY_INT_JOB);
        VRcopy(VpcmPC, WSMuNu, NB2car);
        QFree(WSMuNu);
        QFree(charges);
      }
      //Add to core hamiltonian
      VRadd(h,h,VpcmPC,NB2);
      QFree(VpcmPC);
    }
  }}

  // AWL -- In the case of a QM/MM/PCM calculation...
  // We have just taken care of adding the PCM surface charges, but we still need to add the
  // MM system point charges. So, read them in from disk and set. 
  if(thePCM && rem_read(REM_ONIOM_JOB) == 3 && rem_read(REM_ITESCF) > 1)
  {

     int n_mm_charges;
     FileMan(FM_READ,FILE_PCM_MM_CHARGES,FM_INT,1,0,FM_BEG,&n_mm_charges);
     double *mm_charges = QAllocDouble(5*n_mm_charges);
     FileMan(FM_READ,FILE_PCM_MM_CHARGES,FM_DP,5*n_mm_charges,0,FM_CUR,mm_charges);
     //if(rem_read(REM_FORCEMAN_PRINT) > 2)
     //  MatPrint(mm_charges,5,n_mm_charges,"mm_charges in stvman");
     // ...and then write to the external point charges file, overwriting PCM charges 
     FileMan(FM_WRITE,FILE_EXTERNAL_POINT_CHARGES,FM_INT,1,0,FM_BEG,&n_mm_charges);
     FileMan(FM_WRITE,FILE_EXTERNAL_POINT_CHARGES,FM_DP,5*n_mm_charges,0,FM_CUR,mm_charges);

     // use point charge stuff to set up
     PointCharges mmPC = GetExternalCharges();
     mmPC.Read();
     mmPC.ansToBohr();
     if(Ipcmprint > 2){
       cout << "MM Charges being added to the core h" << endl;
       mmPC.Print();
     }

     if(mmPC.NCharge() > 0){
        ShlPrs S2(code());
        INTEGER NB2car = S2.getNB2car();
        double* VmmPC = QAllocDouble(NB2car);
        // No gaussian blurring of MM system if doing gaussian smooth PCM for now
        if (IGaussianBlur != 1 || (IGaussianBlur == 1 && (IpcmMeth == SMOOTH || IpcmMeth == ISMOOTH))) { 
          PointCharges* pnt_mmPC = &mmPC;
          MakeV(VmmPC,*pnt_mmPC,S2);
        } else{
#ifdef DEVELOPMENT
          printf("evaluating GChg contribution to one-electron hamiltonian from MM charges\n");
#endif
          double* WSMuNu = QAllocDouble(n_mm_charges*NB2car);
          VRload(WSMuNu, NB2car*n_mm_charges, 0.0);
          rem_write(1, REM_MODIFY_INT_JOB);
          AOints(WSMuNu,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,15);
          rem_write(0, REM_MODIFY_INT_JOB);
          VRcopy(VmmPC, WSMuNu, NB2car);
          QFree(WSMuNu);
        }

        //Add to core hamiltonian
        VRadd(h,h,VmmPC,NB2);
        QFree(VmmPC);
     }
     QFree(mm_charges);
  }


  bool SpecialCase = thePCM && rem_read(REM_ONIOM_JOB) == 3 && rem_read(REM_ITESCF) < 1;

  if(SpecialCase){
     // make sure FILE_EXTERNAL_POINT_CHARGES is only the MM charges for the special case
     // we will reset afterwards
     int n_mm_charges;
     FileMan(FM_READ,FILE_PCM_MM_CHARGES,FM_INT,1,0,FM_BEG,&n_mm_charges);
     double *mm_charges = QAllocDouble(5*n_mm_charges);
     FileMan(FM_READ,FILE_PCM_MM_CHARGES,FM_DP,5*n_mm_charges,0,FM_CUR,mm_charges);
     // ...and then write to the external point charges file, overwriting PCM charges 
     FileMan(FM_WRITE,FILE_EXTERNAL_POINT_CHARGES,FM_INT,1,0,FM_BEG,&n_mm_charges);
     FileMan(FM_WRITE,FILE_EXTERNAL_POINT_CHARGES,FM_DP,5*n_mm_charges,0,FM_CUR,mm_charges);
     QFree(mm_charges);
  }


  //AWL -- Regular old way of doing the external charges
  if (ExtCharges == 0 ) QCrash("ExtCharges == 0."); 
  if (ExtCharges->NCharge() > 0 && (thePCM == NULL || SpecialCase)){
    ShlPrs S2(code());
    INTEGER NB2car = S2.getNB2car();
    double* VPC = QAllocDouble(NB2car);
    // Don't do gaussian charges for MM system if doing QM/MM/PCM-SWIG
    if ( IGaussianBlur != 1 ||
         (IGaussianBlur == 1 && (IpcmMeth == SMOOTH || IpcmMeth == ISMOOTH) && thePCM) )
    {
       //ExtCharges->Print();
       MakeV(VPC,*ExtCharges,S2);
       //MatPrint(VPC, -1, "PC Contribution");
       FileMan(FM_WRITE, FILE_PC_POTENTIAL, FM_DP, NB2, 0, FM_BEG, VPC);
       // MESS-E-QM/MM estimation for the polarization energy
       if (rem_read(REM_MESS_E_QMMM) > 0) compute_mess_e_energy(VPC);
       //cout << "MM charges have been added to hamiltonian in the regular way" << endl;
    } else {
       //printf("evaluating GChg contribution to one-electron hamiltonian\n");
       //ExtCharges->Print();
       INTEGER NChg = ExtCharges->NCharge();
       //double* WSMuNu = QAllocDouble(NChg*NB2car);
       // BJK and TAV changed from NChg*NB2car to just NB2car december 2009
       // The extra memory is not written to, and only the first part is
       // zeroed, and there was no change in output, so we believe it safe.
       //double* WSMuNu = QAllocDouble(NB2car);
       // AWL(12/2010): changed to NB2car + NChg 
       // Some QM/MM jobs screw up without some extra memory, so it's not safe!
       double* WSMuNu = QAllocDouble(NB2car + NChg);
       VRload(WSMuNu, NB2car, 0.0);
       rem_write(1, REM_MODIFY_INT_JOB);
       AOints(WSMuNu,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,15);
       //MatPrint(WSMuNu, NChg, NB2car, "SMuNu");
       rem_write(0, REM_MODIFY_INT_JOB);
       VRcopy(VPC, WSMuNu, NB2car);
       QFree(WSMuNu);
    }
    VRadd(h,h,VPC,NB2);
    //MatPrint(h,1,NB2,"Regular core hamiltonian");
    QFree(VPC);
  }

  if(thePCM != NULL) delete thePCM;

  if (EFP2::instance().initialized()) {
	  EFP2::instance().update_wf(h, code());
  }

  // Save updated core Hamiltonian to disk

  FileMan_Open_Write(fNum);
  FileMan_setCrntByteOff(fNum,hOffSet);
  FileMan(FM_WRITE,fNum,FM_DP,NB2,0,FM_CUR,h);
  FileMan_Close(fNum);
}

double* STV::getH()
{
  static LOGICAL Called = FALSE;

  /* We update automatically the first time we are called, and in other
     circumstances such as:

     The applied field doesn't match the last one used
     We are calculating solvation with point charges self-consistently */

  // See if field we're currently using is out of date
  MultipoleField CurField = GetMultipoleField();
  LOGICAL NewField = FALSE;
  LOGICAL doOnsager = (rem_read(REM_SOLVENT_METHOD) == ONSAGER) ? 1 : 0;
  if(!doOnsager){ 
    NewField = MField != CurField;
  }

  // -- See if the external point charges have changed -- //

  // GetExternalCharges converts from angstrom to bohr
  const PointCharges& CurCharges = GetExternalCharges();
  LOGICAL NewCharges = ExtCharges != &CurCharges; 
  bool test = false; // awl, temporary
  //force new charges for pcm and/or QMMM 
  if((rem_read(REM_PCM_METHOD) > 0 && rem_read(REM_ITESCF) >= 2) ||
      (rem_read(REM_QM_MM_INTERFACE) == 2 && rem_read(REM_ITESCF) < 2) ||
      (rem_read(REM_SOLVE_PEQ) && rem_read(REM_PEQ_CORR) == 1) || test)
  {
    NewCharges = true;
  }

  LOGICAL efp_update = EFP2::instance().initialized();

  LOGICAL Update = !Called || NewField || NewCharges || rem_read(REM_SVP) 
        || rem_read(REM_CHEMSOL) == 1 || efp_update
	|| rem_read(REM_CDFTCI_FRAGMENT);


  if (h == 0) h = QAllocDouble(NB2);
  if (Update) {  // Re-create H from components
    if (NewField) MField = CurField;
    if (NewCharges) ExtCharges = &CurCharges;
    makeH();
    Called =  TRUE;
  }
  else {  // Exists but is not loaded
    FileMan_Open_Read(fNum);
    FileMan_setCrntByteOff(fNum,hOffSet);
    FileMan(FM_READ,fNum,FM_DP,NB2,0,FM_CUR,h);
    FileMan_Close(fNum);
  }
  return h;
}

double* STV::getX()
{
  if ( x == 0 )
  {
    x = QAllocDouble(LenX);
    FileMan_Open_Read(fNum);
    FileMan_setCrntByteOff(fNum, xOffSet);
    FileMan(FM_READ,fNum,FM_DP,LenX,0,FM_CUR, x);
    FileMan_Close(fNum);
  }
  return x;
}

double* STV::getSv()
{
  if ( sv == 0 )
  {
    sv = QAllocDouble(NB2);
    FileMan_Open_Read(fNum);
    FileMan_setCrntByteOff(fNum, sOffSet);
    FileMan(FM_READ,fNum,FM_DP,NB2,0,FM_CUR, sv);
    FileMan_Close(fNum);
  }
  return sv;
}

double* STV::getSm()
{
  if ( sm == 0 )
  {
    sm = QAllocDouble(NBasis*NBasis);
    FileMan_Open_Read(fNum);
    FileMan_setCrntByteOff(fNum, sOffSet+NB2);
    FileMan(FM_READ,fNum,FM_DP,NBasis*NBasis,0,FM_CUR, sm);
    FileMan_Close(fNum);
  }
  return sm;
}

double* STV::getT()
{
  if ( t == 0 )
  {
    t = QAllocDouble(NB2);
    FileMan_Open_Read(fNum);
    FileMan_setCrntByteOff(fNum, tOffSet);
    FileMan(FM_READ,fNum,FM_DP,NB2,0,FM_CUR, t);
    FileMan_Close(fNum);
  }
  return t;
}

double* STV::getV()
{
  if ( v == 0 )
  {
    v = QAllocDouble2(NB2,shared);
    FileMan_Open_Read(fNum);
    FileMan_setCrntByteOff(fNum, vOffSet);
    FileMan(FM_READ,fNum,FM_DP,NB2,0,FM_CUR, v);
    FileMan_Close(fNum);
  }
  return v;
}

double* STV::getE()
{
  if (e == 0)
  {
    e = QAllocDouble2(NB2,shared);
    FileMan_Open_Read(fNum);
    FileMan_setCrntByteOff(fNum, eOffSet);
    FileMan(FM_READ,fNum,FM_DP,NB2,0,FM_CUR, e);
    FileMan_Close(fNum);
  }
  return e;
}

double* STV::getVemb()
{
  if(vemb == 0)
  {
    vemb = QAllocDouble2(NB2,shared);
    double *tmp = QAllocDouble(NBasis*NBasis);
    FileMan_Open_Read(FILE_VEMB_DFET);
    FileMan(FM_READ,FILE_VEMB_DFET,FM_DP,NBasis*NBasis,0,FM_BEG, tmp);
    FileMan_Close(FILE_VEMB_DFET);
    ScaV2M(tmp, vemb, true, false);
    QFree(tmp);
  }
  return vemb;
}

double* STV::getM(INTEGER k)
{
  if (k == 1)  // Special case of overlap matrix
    return getSv();
  else if (m == 0) {  // Load 'em up from disk (all at once for now)
    m = QAllocDouble(NB2*(MaxK-1));
    FileMan_Open_Read(fNum);
    FileMan_setCrntByteOff(fNum,mOffSet);
    FileMan(FM_READ,fNum,FM_DP,NB2*(MaxK-1),0,FM_CUR,m);
    FileMan_Close(fNum);
  }

  // Now that we have the mats loaded, return the one we want

  if (k < 2 || k > MaxK)
    QCrash("Requested multipole matrix not available in STV::getM");

  return m + NB2*(k-2);
}

double* STV::getM(INTEGER Lx, INTEGER Ly, INTEGER Lz)
{
  INTEGER k;
  KonL2K(&k,Lx,Ly,Lz);
  return getM(k);
}

void STV::dump()
{
  ShlPrs s2(bCode);
  PrtMat(getH(),-1,7,"Core Hamiltonian Matrix",s2);
  PrtMat(getX(),NBasis,NOrb,NBasis,NOrb,7,"Orthonormalization Matrix",s2);
  PrtMat(getSv(),-1,5,"Overlap Matrix",s2);
  PrtMat(getT(),-1,7,"Kinetic Energy Matrix",s2);
  PrtMat(getV(),-1,7,"Nuclear Attraction Matrix",s2);
  if (rem_read(REM_PSEUDOPOTENTIAL)>0)
     PrtMat(getE(),-1,4,"Effective Core Potential Matrix",s2);
  if(rem_read(REM_VEMB_DFET) == 1) 
        PrtMat(getVemb(),-1,5,"Vemb Matrix from DFET",s2);

  INTEGER Lx,Ly,Lz;
  char Message[80];
  for (INTEGER k = 2; k <= MaxK; ++k) {
    KonK2L(k,&Lx,&Ly,&Lz);
    sprintf(Message,"Multipole Matrix (%d,%d,%d)",Lx,Ly,Lz);
    PrtMat(getM(k),-1,7,Message,s2);
  }
}
