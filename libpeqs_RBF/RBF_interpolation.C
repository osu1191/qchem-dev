#include "peqs_class.h"

#include <cmath>
#include <Eigen/Core>
#include <fstream>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/Eigen>
#include <Eigen/IterativeLinearSolvers>

#include <random>
#include <memory>
#include <functional>
#include <algorithm>
#include <iostream>
#include <vector>

#include "qchem.h"
#include "rem_values.h"
#include "BasisSet.hh"
#include "OneEMtrx.hh"
#include "hirshfeld.h"
#include "$QCHEM_path/anlman/CM5.h"
#include "$QCHEM_path/anlman/makemesh.h"
#include "$QCHEM_path/libdftn/xcclass.h"
#include "$QCHEM_path/libgen/evlbasis.h"
#include "$QCHEM_path/liblas/liblas.h"
#include "ext_libs/fftw/include/fftw.h"
#include "ext_libs/fftw/include/rfftw.h"
#include <chrono>

using namespace std::chrono;
using namespace std;


using std::vector;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::FullPivLU;
using Eigen::PartialPivLU;
using Eigen::ConjugateGradient;
using Eigen::Map;
using Eigen::TriangularView;
using Eigen::Sparse;
using Eigen::Dense;


void PEQSclass::Make_Atom_Centered_Grid()
{
  INTEGER NAtoms = rem_read(REM_NATOMS);
    BasisSet BasisDen(bSetMgr.crntCode());
    BasisSet BasisOrb(bSetMgr.crntCode());
    int bCodeDen = BasisDen.code();
    int grdTyp   = rem_read(REM_IGRDTY);
    int IPrint   = rem_read(REM_PRINT_DFT_GRID);
    int IBCode   = bCodeDen;
    SET_QC_FCOMMON(&IBCode);
    XCFunctional xcFunc = XCFunctional(XCFUNC_SCF);     // read from input

    XCAtoms    *xcatom;
    XCJobPara  *xcpara;
    XCBasisSet *basDen;
    MoleGrid   *mgrid;

    xcatom = new XCAtoms;
    xcpara = new XCJobPara(xcFunc, *xcatom, 0);
    xcpara->nDen = 1;
    int nAtoms = xcatom->getNAtoms();
    double thresh = xcpara->thresh;
    basDen = new XCBasisSet(bCodeDen, nAtoms, thresh);
    xcatom->setSize(*basDen);
    mgrid = new MoleGrid(*xcatom, grdTyp, xcpara->nDrvNuc, thresh);


    int totalgrid(0),offset(0),offset2(0),offset3(0);
    double *jPts2, *jWts;
    int nBatch   = mgrid->getNBatch();
    int nDrvGrd = mgrid->getDrvGrd();
    double nThresh = mgrid->getThresh();
    int NPts = mgrid->getNPts();

    for (int ibat=0; ibat<nBatch; ibat++) {
      BatchGrid grid(*mgrid,ibat);
      nGrid  = grid.getNGrid();
      jPts2  = grid.getPts();
      jWts  = grid.getWts();


      offset = count;
      offset2 = totalgrid+offset;
      offset3 = offset2*3;
        FileMan(FM_WRITE,FILE_PEQS_GRID,FM_DP,nGrid*3,0,FM_CUR,jPts2);
        FileMan(FM_WRITE,FILE_ATM_BATCH_WTS,FM_DP,nGrid,0,FM_CUR,jWts);

      totalgrid += nGrid;
      if (ibat == (nBatch - 1)) count+= totalgrid;
    }
    delete xcatom, xcpara, basDen, mgrid;
}

void PEQSclass::Get_Electronic_Density_on_ACG(double *PAv)
{
  int i,j,k,ijk;
  double sum1=0.0,sum2=0.0,sum3=0.0,sumC1=0.0,sumC2=0.0;
  INTEGER NBasis = bSetMgr.crntShlsStats(STAT_NBASIS);
  INTEGER N2 = NBasis*NBasis;
  double *kPts2 = QAllocDouble(count*3);  // for ACG testing
  double *kWts = QAllocDouble(count);   // for ACG testing
  FileMan(FM_READ,FILE_PEQS_GRID,FM_DP,count*3,0,FM_BEG,kPts2);  // for ACG testing
  FileMan(FM_READ,FILE_ATM_BATCH_WTS,FM_DP,count,0,FM_BEG,kWts);   // for ACG testing

  qtime_t shackle08;
  double shackle09[3];
  shackle08 = QTimerOn();
  int cycles=MaxPPer;
  int acgPPer=cycles, offset=0;
  double *irregPts=QAllocDouble(3*cycles),*acg_phi=QAllocDouble(cycles);
  phi_irr=QAllocDouble(count);
  VRload(phi_irr,count,0.0);
  for(int l=0;l<count;l+=cycles){
  if (cycles <= (count-l)) acgPPer = cycles;
  else acgPPer = (count-l);
    VRload(acg_phi,acgPPer,0.0);
    VRcopy(irregPts,&kPts2[(3*l)],3*acgPPer);   // THIS IS THE DOUBTFUL PART... you have to copy all the ACG points using this line!!!
    AOints(acg_phi,NULL,NULL,NULL,PAv,irregPts,NULL,NULL,&acgPPer,12);
    VRcopy(&phi_irr[l],acg_phi,acgPPer); // Latter on scale the phi_elec by -1
  }
  cout << "Potential have been calculated" << endl;

   QTimerOff(shackle09,shackle08);
   printf(" Potential calculation in ACG points completed in %2.2f (cpu) %2.2f (wall) seconds. \n",shackle09[0],shackle09[2]);
   QFree(irregPts); QFree(acg_phi);
   QFree(kPts2); QFree(kWts);
}


void PEQSclass::Get_Neigh()
{
   int i,j,k,l,ijk;
   int tot_count;
   double xd,yd,zd;
   double *kPts2 = QAllocDouble(count*3);  // for ACG testing

   double *dist = QAllocDouble(count);
   double *tst_dist = QAllocDouble(count);

   INTEGER *ind = QAllocINTEGER(count);
   INTEGER *tst_ind = QAllocINTEGER(count);

   FileMan(FM_READ,FILE_PEQS_GRID,FM_DP,count*3,0,FM_BEG,kPts2);  // for ACG testing

   for(ijk=0; ijk<NTotPts; ijk++){
     get_3ind_from_multInd(i,j,k,ijk);
     double xx = valueOfX(i);
     double yy = valueOfY(j);
     double zz = valueOfZ(k);

     for(l=0; l<count; l++){
	xd = (xx-kPts2[3*l])*(xx-kPts2[3*l]);
     	yd = (yy-kPts2[3*l+1])*(yy-kPts2[3*l+1]);
     	zd = (zz-kPts2[3*l+2])*(zz-kPts2[3*l+2]);
     	dist[l] = sqrt(xd + yd + zd);
     	ind[l] = l;
	tst_dist[l] = dist[l];
	tst_ind[l] = ind[l];
     }
 
     int n = 0;
     double rad = domRad;
     int tst_count = count;
     double *near = QAllocDouble(count);
     INTEGER *near_ind = QAllocINTEGER(count);
     while(n < neigh_cutoff){
	int r = 0;

	for(l=0; l<tst_count; l++){
	  if(tst_dist[l] < rad){
	    near[n] = tst_dist[l];
	    near_ind[n] = tst_ind[l];
	    n++;
	  }
	  else{
	    tst_dist[r] = dist[l];
	    tst_ind[r] = ind[l];
	    r++;
	  }
	}
	tst_count = r;
	rad += (1.0*ang2bohr);
     }
     tot_count = n;
     mergeSort(near, near_ind, 0, tot_count-1);     
     for(l=0; l<500; l++){
	neighbors[ijk*500+l] = near_ind[l];
     }
     
    QFree(near), QFree(near_ind); 
   }

   QFree(kPts2), QFree(dist), QFree(tst_dist);
   QFree(ind), QFree(tst_ind);;
}


void PEQSclass::Local_RBF()
{
   double *hPts2=QAllocDouble(count*3);
   int ijk,ibox,jbox,kbox,i,j,k,l,m,n=0,sum_occ=0;
   int Nx=NPoints[0],Ny=NPoints[1],Nz=NPoints[2];

   FileMan(FM_READ,FILE_PEQS_GRID,FM_DP,count*3,0,FM_BEG,hPts2);  // ACG points
 
   qtime_t neigh02;
   double neigh03[3];
   neigh02 = QTimerOn();
   Get_Neigh();   // Neighbor Search
   QTimerOff(neigh03,neigh02);
   printf("\tPEqS:  Time taken for neighbor search in %2.2f (cpu) %2.2f (wall) sec. \n",
           neigh03[0],neigh03[2]);  
 
   VRload(rbf_itp,NTotPts,0.0);

   VectorXd xvec;
   #pragma omp parallel for
   for(ijk=0; ijk<NTotPts; ijk++){
     get_3ind_from_multInd(i,j,k,ijk);
     xvec = VectorXd::Zero(3);
     double xx = valueOfX(i);
     double yy = valueOfY(j);
     double zz = valueOfZ(k);
     xvec(0) = xx;
     xvec(1) = yy;
     xvec(2) = zz; 
 
     MatrixXd X_acg = MatrixXd::Zero(3, 500);
     VectorXd y_acg = VectorXd::Zero(500);
    
     for(l=0; l<500; l++){
	X_acg(0, l) = hPts2[3*neighbors[l]];
	X_acg(1, l) = hPts2[3*neighbors[l]+1];
	X_acg(2, l) = hPts2[3*neighbors[l]+2];	 
	y_acg(l) = phi_irr[neighbors[l]];
     }

     VectorXd phi = Map<VectorXd>(&y_acg[0], 500);
     GenerateMatrices(X_acg, phi);
  
     qtime_t weight02;
     double weight03[3];
     weight02 = QTimerOn();
     ComputeWeights();   // Computing Weights
     QTimerOff(weight03,weight02);
     printf("\tPEqS:  Time taken for computing weights in %2.2f (cpu) %2.2f (wall) sec. \n",
             weight03[0],weight03[2]);


     qtime_t interp02;
     double interp03[3];
     interp02 = QTimerOn();
     assert(xvec.rows() == m_X.rows());
     int num_data = m_w.rows();
     int dim = xvec.rows();
     Eigen::VectorXd norms = (m_X.colwise() - xvec).colwise().norm();
     VectorXd rbf_values{num_data};
     for (int i = 0; i < num_data; ++i){
        rbf_values(i) = GenerateRBF(norms(i));
     }

     double rbf_term = m_w.dot(rbf_values);
     if (m_use_polynomial_term){
       double polynomial_term = m_v(0) + xvec.transpose() * m_v.segment(1, dim);
       rbf_itp[ijk] = rbf_term + polynomial_term;
     }
     else{
       rbf_itp[ijk] = rbf_term;
     }
     QTimerOff(interp03,interp02);
     printf("\tPEqS:  Time taken for interp values in %2.2f (cpu) %2.2f (wall) sec. \n",
             interp03[0],interp03[2]);

   }

}


void PEQSclass::GenerateMatrices(const MatrixXd& X, const VectorXd& y)
{
    assert(y.rows() == X.cols());

    this->m_X = X;
    this->m_y = y;

}

void PEQSclass::ComputeWeights()
{
    const int num_data = m_y.rows();
    double lambda = 0.001;
    double value = 0.0;
 
    MatrixXd Phi = MatrixXd::Zero(num_data, num_data);
    for (int i = 0; i < num_data; ++i)
    {
        for (int j = i; j < num_data; ++j)
        {
            double r = (m_X.col(i) - m_X.col(j)).norm();
	    value = GenerateRBF(r);
 
            Phi(i, j) = value;
            Phi(j, i) = value;
        }
    }


    if (m_use_polynomial_term)
    {

    const int dim = m_X.rows();

    MatrixXd P{num_data, dim + 1};

    P.block(0, 0, num_data, 1)   = MatrixXd::Ones(num_data, 1);
    P.block(0, 1, num_data, dim) = m_X.transpose();

    MatrixXd A = MatrixXd::Zero(num_data + dim + 1, num_data + dim + 1);

    A.block(0, 0, num_data, num_data)       = Phi;
    A.block(0, num_data, num_data, dim + 1) = P;
    A.block(num_data, 0, dim + 1, num_data) = P.transpose();

    VectorXd b = VectorXd::Zero(num_data + dim + 1);

    b.segment(0, num_data) = m_y;

    VectorXd solution = VectorXd::Zero(num_data + dim + 1);

    if (use_regularization)
    {
    const auto I = MatrixXd::Identity(num_data + dim + 1, num_data + dim + 1);
    solution = ConjugateGradient<MatrixXd,Eigen::Lower|Eigen::Upper>(A.transpose() * A + lambda * I).solve(A.transpose() * b);

    }
    else
    {
	solution = ConjugateGradient<MatrixXd,Eigen::Lower|Eigen::Upper>(A).solve(b);	
    }

    m_w = solution.segment(0, num_data);
    m_v = solution.segment(num_data, dim + 1);


    } 
    else
    {
        const auto     I = MatrixXd::Identity(num_data, num_data);
        const MatrixXd A = use_regularization ? Phi.transpose() * Phi + lambda * I : Phi;
        const VectorXd b = use_regularization ? Phi.transpose() * m_y : m_y;

	int n = Eigen::nbThreads();
	omp_set_num_threads(8);
	Eigen::setNbThreads(8);
	m_w = ConjugateGradient<MatrixXd,Eigen::Lower|Eigen::Upper>(A).solve(b);

    }
}

double PEQSclass::GenerateRBF(const double r) const
{
    double result;

    if (rbf_kernel == RBF_KERNEL_GAUSSIAN){
        result = exp(- pow((nepsilon * r), 2.0));
    }
    else if (rbf_kernel == RBF_KERNEL_THINPLATESPLINE){
        if(r > 0.0)result = r * r * log(r);
	else result = 0.0;
    }
    else if (rbf_kernel == RBF_KERNEL_INVERSEQUADRATIC){
        result = 1.0 / (1.0 + pow((nepsilon * r), 2.0));
    }
    else if (rbf_kernel == RBF_KERNEL_BIHARMONICSPLINE){
        result = r;
    }

    return result;
}

VectorXd solveLinearSystem(const MatrixXd& A, const VectorXd& y)
{
    FullPivLU<MatrixXd> lu(A);
    return lu.solve(y);
}


