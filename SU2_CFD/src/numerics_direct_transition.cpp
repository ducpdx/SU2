/*!
 * \file numerics_direct_transition.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.9
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/numerics_structure.hpp"
#include <limits>

#include "c_routines_d/pow_d.c"
CUpwLin_TransLM::CUpwLin_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	grid_movement  = config->GetGrid_Movement();
	incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
  
}

CUpwLin_TransLM::~CUpwLin_TransLM(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
}

void CUpwLin_TransLM::ComputeResidual (double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {
  
  
	Density_i = U_i[0];
	Density_j = U_j[0];
  
	q_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1]/Density_i;
		Velocity_j[iDim] = U_j[iDim+1]/Density_j;
		q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
	}
  
	a0 = 0.5*(q_ij+fabs(q_ij));
	a1 = 0.5*(q_ij-fabs(q_ij));
	val_residual[0] = a0*TransVar_i[0]+a1*TransVar_j[0];
	val_residual[1] = a0*TransVar_i[1]+a1*TransVar_j[1];
  //	cout << "Velicity x: " << Velocity_i[0] << ", " << Velocity_j[0] << endl;
  //	cout << "Velicity y: " << Velocity_i[1] << ", " << Velocity_j[1] << endl;
  //	cout << "val_resid: " << val_residual[0] << ", " << val_residual[1] << endl;
  

	if (implicit) {
		val_Jacobian_i[0][0] = a0;
		val_Jacobian_i[1][1] = a0;
	}
}

CUpwSca_TransLM::CUpwSca_TransLM(unsigned short val_nDim, unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	grid_movement = config->GetGrid_Movement();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
}

CUpwSca_TransLM::~CUpwSca_TransLM(void) {
	delete [] Velocity_i;
	delete [] Velocity_j;
}

void CUpwSca_TransLM::ComputeResidual(double *val_residual, double **val_Jacobian_i, double **val_Jacobian_j, CConfig *config) {

	q_ij = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		q_ij += 0.5*(U_i[iDim+1]+U_j[iDim+1])*Normal[iDim];
	}

	a0 = 0.5*(q_ij+fabs(q_ij));
	a1 = 0.5*(q_ij-fabs(q_ij));
	val_residual[0] = a0*TransVar_i[0]+a1*TransVar_j[0];
	val_residual[1] = a0*TransVar_i[1]+a1*TransVar_j[1];

	if (implicit) {
		val_Jacobian_i[0][0] = a0;
		val_Jacobian_j[0][0] = a1;
		val_Jacobian_i[1][1] = a0;
		val_Jacobian_j[1][1] = a1;

		/* --- Zero out off-diagonal terms just in case ---*/
		val_Jacobian_i[0][1] = 0;
		val_Jacobian_j[0][1] = 0;
		val_Jacobian_i[1][0] = 0;
		val_Jacobian_j[1][0] = 0;
	}

}

CAvgGrad_TransLM::CAvgGrad_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	unsigned short iVar;
  
	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	sigma = 2./3.;
  
	Edge_Vector = new double [nDim];
	Proj_Mean_GradTransVar_Kappa = new double [nVar];
	Proj_Mean_GradTransVar_Edge = new double [nVar];
	Mean_GradTransVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradTransVar[iVar] = new double [nDim];
}

CAvgGrad_TransLM::~CAvgGrad_TransLM(void) {
  
	unsigned short iVar;
  
	delete [] Edge_Vector;
	delete [] Proj_Mean_GradTransVar_Kappa;
	delete [] Proj_Mean_GradTransVar_Edge;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradTransVar[iVar];
	delete [] Mean_GradTransVar;
}

void CAvgGrad_TransLM::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {

  //************************************************//
  // Please do not delete //SU2_CPP2C comment lines //
  //************************************************//

  double TransVar_id[nVar], TransVar_jd[nVar], val_residuald[nVar];
  //SU2_CPP2C START CAvgGrad_TransLM::SetResidual
  //SU2_CPP2C CALL_LIST START
  //SU2_CPP2C INVARS *TransVar_i *TransVar_j
  //SU2_CPP2C OUTVARS *val_residual
  //SU2_CPP2C VARS DOUBLE *U_i *U_j
  //SU2_CPP2C VARS DOUBLE Laminar_Viscosity_i Eddy_Viscosity_i
  //SU2_CPP2C VARS DOUBLE Laminar_Viscosity_j Eddy_Viscosity_j
  //SU2_CPP2C VARS DOUBLE *Coord_i *Coord_j *Normal
  //SU2_CPP2C VARS DOUBLE **ConsVar_Grad_i **ConsVar_Grad_j **TransVar_Grad_i **TransVar_Grad_j
  //SU2_CPP2C CALL_LIST END

  //SU2_CPP2C DEFINE nDim nVar

  //SU2_CPP2C DECL_LIST START
  //SU2_CPP2C VARS DOUBLE MATRIX SIZE=nDim Edge_Vector Normal Proj_Mean_GradTransVar_Kappa
  //SU2_CPP2C VARS DOUBLE MATRIX SIZE=nDim SIZE=nVar Mean_GradTransVar
  //SU2_CPP2C VARS DOUBLE SCALAR Density_i Density_j dist_ij_2 proj_vector_ij 
  //SU2_CPP2C VARS INT SCALAR iDim iVar
  //SU2_CPP2C DECL_LIST END
  double Density_Grad_i[nDim], Density_Grad_j[nDim], Conservative_Grad_i[nDim], Conservative_Grad_j[nDim];
  double Primitive_Grad_i[nDim], Primitive_Grad_j[nDim];

  /*--- Intermediate values for combining viscosities ---*/
  double Inter_Viscosity_i, Inter_Viscosity_j, REth_Viscosity_i, REth_Viscosity_j, Inter_Viscosity_Mean, REth_Viscosity_Mean;

  /*--- Model constants---*/
  double sigmaf       = 1.0;
  double sigma_thetat = 2.0;

  /*--- Get density ---*/
  Density_i = U_i[0]; 
  Density_j = U_j[0];

  /*--- Construct combinations of viscosity ---*/
  Inter_Viscosity_i    = (Laminar_Viscosity_i+Eddy_Viscosity_i/sigmaf);
  Inter_Viscosity_j    = (Laminar_Viscosity_j+Eddy_Viscosity_j/sigmaf);
  Inter_Viscosity_Mean = 0.5*(Inter_Viscosity_i+Inter_Viscosity_j);
  REth_Viscosity_i     = sigma_thetat*(Laminar_Viscosity_i+Eddy_Viscosity_i);
  REth_Viscosity_j     = sigma_thetat*(Laminar_Viscosity_j+Eddy_Viscosity_j);
  REth_Viscosity_Mean  = 0.5*(REth_Viscosity_i+REth_Viscosity_j);

  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  proj_vector_ij = proj_vector_ij/dist_ij_2; // to normalize vectors

  /*--- Mean gradient approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTransVar_Kappa[iVar] = 0.0;
    //Proj_Mean_GradTransVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {

      /* -- Compute primitive grad using chain rule -- */
      Density_Grad_i[iDim]      = ConsVar_Grad_i[0][iDim];
      Density_Grad_j[iDim]      = ConsVar_Grad_j[0][iDim];
      Conservative_Grad_i[iDim] = TransVar_Grad_i[iVar][iDim];
      Conservative_Grad_j[iDim] = TransVar_Grad_j[iVar][iDim];
      Primitive_Grad_i[iDim]    = 1./Density_i*(Conservative_Grad_i[iDim]-TransVar_i[iVar]*Density_Grad_i[iDim]);
      Primitive_Grad_j[iDim]    = 1./Density_j*(Conservative_Grad_j[iDim]-TransVar_j[iVar]*Density_Grad_j[iDim]);

      /*--- Compute the average primitive gradient and project it in the normal direction ---*/
      Mean_GradTransVar[iVar][iDim] = 0.5*(Primitive_Grad_i[iDim] + Primitive_Grad_j[iDim]);
      Proj_Mean_GradTransVar_Kappa[iVar] += Mean_GradTransVar[iVar][iDim]*Normal[iDim];
    }
  }

  val_residual[0] = Inter_Viscosity_Mean*Proj_Mean_GradTransVar_Kappa[0];
  val_residual[1] = REth_Viscosity_Mean*Proj_Mean_GradTransVar_Kappa[1];
  
  //SU2_CPP2C END CAvgGrad_TransLM::SetResidual

  if (implicit) {

    TransVar_id[0] = 1.0; TransVar_id[1] = 0.0;
    TransVar_jd[0] = 0.0; TransVar_jd[1] = 0.0;
    ComputeResidual_d(TransVar_i, TransVar_id,  TransVar_j,  TransVar_jd,  val_residual,  val_residuald);
    Jacobian_i[0][0] = val_residuald[0];
    Jacobian_i[1][0] = val_residuald[1];

    TransVar_id[0] = 0.0; TransVar_id[1] = 0.0;
    TransVar_jd[0] = 1.0; TransVar_jd[1] = 0.0;
    ComputeResidual_d(TransVar_i, TransVar_id,  TransVar_j,  TransVar_jd,  val_residual,  val_residuald);
    Jacobian_j[0][0] = val_residuald[0];
    Jacobian_j[1][0] = val_residuald[1];

    TransVar_id[0] = 0.0; TransVar_id[1] = 1.0;
    TransVar_jd[0] = 0.0; TransVar_jd[1] = 0.0;
    ComputeResidual_d(TransVar_i, TransVar_id,  TransVar_j,  TransVar_jd,  val_residual,  val_residuald);
    Jacobian_i[0][1] = val_residuald[0];
    Jacobian_i[1][1] = val_residuald[1];

    TransVar_id[0] = 0.0; TransVar_id[1] = 0.0;
    TransVar_jd[0] = 0.0; TransVar_jd[1] = 1.0;
    ComputeResidual_d(TransVar_i, TransVar_id,  TransVar_j,  TransVar_jd,  val_residual,  val_residuald);
    Jacobian_j[0][1] = val_residuald[0];
    Jacobian_j[1][1] = val_residuald[1];

  }

}

void CAvgGrad_TransLM::ComputeResidual_d(double *TransVar_i, double *TransVar_id, double *TransVar_j, double *TransVar_jd, double *val_residual, double *val_residuald)
{
//SU2_INSERT START
    double Density_Grad_i[nDim], Density_Grad_j[nDim], 
    Conservative_Grad_i[nDim], Conservative_Grad_j[nDim];
    double Density_Grad_id[nDim], Density_Grad_jd[nDim], 
    Conservative_Grad_id[nDim], Conservative_Grad_jd[nDim];
    double Primitive_Grad_i[nDim], Primitive_Grad_j[nDim];
    double Primitive_Grad_id[nDim], Primitive_Grad_jd[nDim];
    /*--- Intermediate values for combining viscosities ---*/
    double Inter_Viscosity_i, Inter_Viscosity_j, REth_Viscosity_i, 
    REth_Viscosity_j, Inter_Viscosity_Mean, REth_Viscosity_Mean;
    /*--- Model constants---*/
    double sigmaf = 1.0;
    double sigma_thetat = 2.0;
    double Edge_Vectord[nDim];
    double Proj_Mean_GradTransVar_Kappad[nDim];
    double Mean_GradTransVard[nDim][nVar];
    int ii2;
    int ii1;
    /*--- Get density ---*/
    Density_i = U_i[0];
    Density_j = U_j[0];
    /*--- Construct combinations of viscosity ---*/
    Inter_Viscosity_i = Laminar_Viscosity_i + Eddy_Viscosity_i/sigmaf;
    Inter_Viscosity_j = Laminar_Viscosity_j + Eddy_Viscosity_j/sigmaf;
    Inter_Viscosity_Mean = 0.5*(Inter_Viscosity_i+Inter_Viscosity_j);
    REth_Viscosity_i = sigma_thetat*(Laminar_Viscosity_i+Eddy_Viscosity_i);
    REth_Viscosity_j = sigma_thetat*(Laminar_Viscosity_j+Eddy_Viscosity_j);
    REth_Viscosity_Mean = 0.5*(REth_Viscosity_i+REth_Viscosity_j);
    /*--- Compute vector going from iPoint to jPoint ---*/
    dist_ij_2 = 0;
    proj_vector_ij = 0;
    for (iDim = 0; iDim < nDim; ++iDim) {
        Edge_Vectord[iDim] = 0.0;
        Edge_Vector[iDim] = Coord_j[iDim] - Coord_i[iDim];
        dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
        proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
    }
    proj_vector_ij = proj_vector_ij/dist_ij_2;
    for (ii1 = 0; ii1 < nDim; ++ii1)
        for (ii2 = 0; ii2 < nVar; ++ii2)
            Mean_GradTransVard[ii1][ii2] = 0.0;
    for (ii1 = 0; ii1 < nDim; ++ii1)
        Primitive_Grad_id[ii1] = 0.0;
    for (ii1 = 0; ii1 < nDim; ++ii1)
        Primitive_Grad_jd[ii1] = 0.0;
    for (ii1 = 0; ii1 < nDim; ++ii1)
        Proj_Mean_GradTransVar_Kappad[ii1] = 0.0;
    /*--- Mean gradient approximation ---*/
    // to normalize vectors
    for (iVar = 0; iVar < nVar; ++iVar) {
        Proj_Mean_GradTransVar_Kappad[iVar] = 0.0;
        Proj_Mean_GradTransVar_Kappa[iVar] = 0.0;
        // Proj_Mean_GradTransVar_Edge[iVar] = 0.0;
        for (iDim = 0; iDim < nDim; ++iDim) {
            /* -- Compute primitive grad using chain rule -- */
            Density_Grad_id[iDim] = 0.0;
            Density_Grad_i[iDim] = ConsVar_Grad_i[0][iDim];
            Density_Grad_jd[iDim] = 0.0;
            Density_Grad_j[iDim] = ConsVar_Grad_j[0][iDim];
            Conservative_Grad_id[iDim] = 0.0;
            Conservative_Grad_i[iDim] = TransVar_Grad_i[iVar][iDim];
            Conservative_Grad_jd[iDim] = 0.0;
            Conservative_Grad_j[iDim] = TransVar_Grad_j[iVar][iDim];
            Primitive_Grad_id[iDim] = -(Density_Grad_i[iDim]*TransVar_id[iVar]
                /Density_i);
            Primitive_Grad_i[iDim] = 1./Density_i*(Conservative_Grad_i[iDim]-
                TransVar_i[iVar]*Density_Grad_i[iDim]);
            Primitive_Grad_jd[iDim] = -(Density_Grad_j[iDim]*TransVar_jd[iVar]
                /Density_j);
            Primitive_Grad_j[iDim] = 1./Density_j*(Conservative_Grad_j[iDim]-
                TransVar_j[iVar]*Density_Grad_j[iDim]);
            /*--- Compute the average primitive gradient and project it in the normal direction ---
            */
            Mean_GradTransVard[iVar][iDim] = 0.5*(Primitive_Grad_id[iDim]+
                Primitive_Grad_jd[iDim]);
            Mean_GradTransVar[iVar][iDim] = 0.5*(Primitive_Grad_i[iDim]+
                Primitive_Grad_j[iDim]);
            Proj_Mean_GradTransVar_Kappad[iVar] = 
                Proj_Mean_GradTransVar_Kappad[iVar] + Normal[iDim]*
                Mean_GradTransVard[iVar][iDim];
            Proj_Mean_GradTransVar_Kappa[iVar] += Mean_GradTransVar[iVar][iDim
            ]*Normal[iDim];
        }
    }
    val_residuald[0] = Inter_Viscosity_Mean*Proj_Mean_GradTransVar_Kappad[0];
    val_residual[0] = Inter_Viscosity_Mean*Proj_Mean_GradTransVar_Kappa[0];
    val_residuald[1] = REth_Viscosity_Mean*Proj_Mean_GradTransVar_Kappad[1];
    val_residual[1] = REth_Viscosity_Mean*Proj_Mean_GradTransVar_Kappa[1];

//SU2_INSERT END
}

CAvgGradCorrected_TransLM::CAvgGradCorrected_TransLM(unsigned short val_nDim, unsigned short val_nVar,
                                                     CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	unsigned short iVar;
  
	implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
	incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	sigma = 2./3.;
  
	Edge_Vector = new double [nDim];
	Proj_Mean_GradTurbVar_Kappa = new double [nVar];
	Proj_Mean_GradTurbVar_Edge = new double [nVar];
	Proj_Mean_GradTurbVar_Corrected = new double [nVar];
	Mean_GradTurbVar = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		Mean_GradTurbVar[iVar] = new double [nDim];
}

CAvgGradCorrected_TransLM::~CAvgGradCorrected_TransLM(void) {
  
	unsigned short iVar;
  
	delete [] Edge_Vector;
	delete [] Proj_Mean_GradTurbVar_Kappa;
	delete [] Proj_Mean_GradTurbVar_Edge;
	delete [] Proj_Mean_GradTurbVar_Corrected;
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] Mean_GradTurbVar[iVar];
	delete [] Mean_GradTurbVar;
}

void CAvgGradCorrected_TransLM::ComputeResidual(double *val_residual, double **Jacobian_i, double **Jacobian_j, CConfig *config) {
  
  //	switch (config->GetKind_Turb_Model()) {
  //	case SA :
  //		/*--- Compute mean effective viscosity ---*/
  //		nu_i = Laminar_Viscosity_i/U_i[0];
  //		nu_j = Laminar_Viscosity_j/U_j[0];
  //		nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  //
  //		/*--- Compute vector going from iPoint to jPoint ---*/
  //		dist_ij_2 = 0; proj_vector_ij = 0;
  //		for (iDim = 0; iDim < nDim; iDim++) {
  //			Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
  //			dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  //			proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  //		}
  //		proj_vector_ij = proj_vector_ij/dist_ij_2;
  //
  //		/*--- Mean gradient approximation. Projection of the mean gradient
  //			 in the direction of the edge ---*/
  //		for (iVar = 0; iVar < nVar; iVar++) {
  //			Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
  //			Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
  //			for (iDim = 0; iDim < nDim; iDim++) {
  //				Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
  //				Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
  //				Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
  //			}
  //			Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
  //			Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
  //					(TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  //		}
  //
  //		val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  //
  //		/*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  //		if (implicit) {
  //			Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
  //			Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
  //		}
  //		break;
  //
  //	}
}

CSourcePieceWise_TransLM::CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	/*--- Spalart-Allmaras closure constants ---*/
	cv1_3 = pow(7.1,3.0);
	k2 = pow(0.41,2.0);
	cb1 = 0.1355;
	cw2 = 0.3;
	cw3_6 = pow(2.0,6.0);
	sigma = 2./3.;
	cb2 = 0.622;
	cw1 = cb1/k2+(1+cb2)/sigma;
  
	/*-- Gamma-theta closure constants --*/
	c_e1    = 1.0;
	c_a1    = 2.0;
	c_e2    = 50.0;
	c_a2    = 0.06;
	sigmaf  = 1.0;
	s1      = 2.0;
	c_theta = 0.03;
	sigmat  = 2.0;
  
	/*-- Correlation constants --*/
	flen_global  = 12.0;
	alpha_global = 0.85;

	/*-- Initialize member variables --*/
	strain   = 0.0;
	f_lambda = 0.0;

}

CSourcePieceWise_TransLM::~CSourcePieceWise_TransLM(void) { }

double CSourcePieceWise_TransLM::corr_func(double lambda_in) {

	double re_thetaRHS, re_thetaLHS;
	double f;

	/*-- compute f --*/
	re_thetaLHS = sqrt(U_i[0]*pow(Velocity_Mag,2)/(Laminar_Viscosity_i*du_ds)*lambda_in);

	if (lambda_in <= 0.0)
		f_lambda = 1 - (-12.986*lambda_in - 123.66*pow(lambda_in,2) - 405.689*pow(lambda_in,3))*exp(-pow(tu/1.5,1.5));
	else
		f_lambda = 1 + 0.275*(1-exp(-35.0*lambda_in))*exp(-tu/0.5);

	if (tu <= 1.3)
		re_thetaRHS = (1173.51 - 589.428*tu + 0.2196/pow(tu,2))*f_lambda;
	else
		re_thetaRHS = (331.5*pow(tu-0.5658,-0.671))*f_lambda;

	// cout << "U_i[0]=" << U_i[0] << ", Velocity_Mag = " << Velocity_Mag << ", Laminar_Viscosity_i = " << Laminar_Viscosity_i << ", du_ds = " << du_ds << endl;
	// cout << "LHS, RHS: " << re_thetaLHS << ", " << re_thetaRHS << endl;
	f = re_thetaLHS-re_thetaRHS;

	return f;

}


void CSourcePieceWise_TransLM::translm_helper(CConfig *config) {

	double dU_dx, dU_dy, dU_dz;

	double lambda_a, lambda_b, lambda_c, lambda;

	double f_a, f_b, f_c;
	double bracket_tol = 1e-10;
	double lambda_check;

	bool debug=false;


	rey  = config->GetReynolds();
	mach = config->GetMach_FreeStreamND();
	tu   = config->GetTurbulenceIntensity_FreeStream();

	/*--- Compute vorticity and strain ---*/
	Vorticity = (PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1])*(PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1]);
	if (nDim == 3) Vorticity += ( (PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2])*(PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2]) +
                                 (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) );

	/*-- Strain = sqrt(2*Sij*Sij) --*/
	strain = 0.0;
	for (int i=0; i<nDim; i++) {
		for(int j=0; j<nDim; j++) {
			strain += PrimVar_Grad_i[i+1][j]*PrimVar_Grad_i[i+1][j];
		}
	}
	strain = sqrt(2*strain);

	if (nDim==2) {
		Velocity_Mag = sqrt(U_i[1]*U_i[1]+U_i[2]*U_i[2])/U_i[0];
	} else if (nDim==3) {
		Velocity_Mag = sqrt(U_i[1]*U_i[1]+U_i[2]*U_i[2]+U_i[3]*U_i[3])/U_i[0];
	}

	/*-- Gradient of velocity magnitude ---*/
	dU_dx = 0.5*Velocity_Mag*( 2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][0]
	                          +2*U_i[2]/U_i[0]*PrimVar_Grad_i[2][0]);
	if (nDim==3)
		dU_dx += 0.5*Velocity_Mag*( 2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][0]);

	dU_dy = 0.5*Velocity_Mag*( 2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][1]
	                          +2*U_i[2]/U_i[0]*PrimVar_Grad_i[2][1]);
	if (nDim==3)
		dU_dy += 0.5*Velocity_Mag*( 2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][1]);

	if (nDim==3)
		dU_dz = 0.5*Velocity_Mag*( 2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][2]
		                          +2*U_i[2]/U_i[0]*PrimVar_Grad_i[2][2]
		                          +2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][2]);

	du_ds = U_i[1]/(U_i[0]*Velocity_Mag) * dU_dx +  // Streamwise velocity derivative
			U_i[2]/(U_i[0]*Velocity_Mag) * dU_dy;
	if (nDim==3)
		du_ds += U_i[3]/(U_i[0]*Velocity_Mag) * dU_dz;

	/*-- Begin iterations to solve REth correlation --*/
	f_lambda = 1.;
	tu = tu*100.;

	/*-- Initial bracket for lambda --*/
	if (du_ds>=0) {
		//cout << "Case A, du_ds>=0. " ;
		lambda_a = 0.0;
		lambda_b = 0.1;
	} else {
		//cout << "Case B, du_ds<0. " ;
		lambda_a = -0.1;
		lambda_b =  0.0;
	}

	if (abs(du_ds)<1e-12) {
		lambda = 0.0;
	} else {
		f_a = corr_func(lambda_a);
		f_b = corr_func(lambda_b);
		if (f_a*f_b > 0) {
			//cout << "Case I: f_a and f_b have same sign." << endl;
			lambda = 0.1*copysign(1.0,du_ds);
		} else {

			//cout << "Case II: f_a and f_b have same sign." << endl;
			/*-- Begin method of false position to solve for lambda --*/
			while(1) {

				/*-- Mathews and Fink, p.57 eq. 22 --*/
				lambda_c = lambda_b - (f_b*(lambda_b-lambda_a))/(f_b-f_a);
				f_c = corr_func(lambda_c);

				/*-- Monitor convergence --*/
				//cout << "\t lambda: [" << lambda_a << ", " << lambda_c << ", " << lambda_b << "] \t f: [" << f_a << ", " << f_c << ", " << f_b << "]" << endl;

				if (f_c==0.0) {
					lambda = lambda_c;
					//cout << "Exact solution found: lambda=" << lambda << endl;
					break;
				}
				if (lambda_c-lambda_a <= bracket_tol) {
					lambda = 0.5*(lambda_a+lambda_c);
					//cout << "A-C Bracket tolerance reached. lambda = " << lambda << endl;
					break;
				} else if (lambda_b-lambda_c <= bracket_tol) {
					lambda = 0.5*(lambda_b+lambda_c);
					//cout << "C-B Bracket tolerance reached. lambda = " << lambda << endl;
					break;
				}

				/*-- Modify bracket for next iteration --*/
				if (f_a*f_c <=0) {
					lambda_b = lambda_c;
					f_b      = f_c;
				} else {
					lambda_a = lambda_c;
					f_a      = f_c;
				}
			}
			debug = true;
		}
	}


  /*-- Evaluate correlation with solved value of lambda --*/
	if (lambda <= 0.0) f_lambda = 1 - (-12.986*lambda - 123.66*pow(lambda,2) - 405.689*pow(lambda,3))*exp(-pow(tu/1.5,1.5));
	else f_lambda = 1 + 0.275*(1-exp(-35.0*lambda))*exp(-tu/0.5);

	if (tu <= 1.3)
		re_theta_t = (1173.51 - 589.428*tu + 0.2196/pow(tu,2))*f_lambda;
	else
		re_theta_t = (331.5*pow(tu-0.5658,-0.671))*f_lambda;

	// Check lambda, TODO: Remove this
//	if (debug) {
//		lambda_check = pow(re_theta_t,2)*Laminar_Viscosity_i/(U_i[0]*pow(Velocity_Mag,2))*du_ds;
//		double check = abs(lambda-lambda_check);
//		if (check>1e-8)
//			cout << "WARNING: lambda_check does not match lambda. " << lambda << " " << lambda_check << endl;
//	}

	/*-- Restore tu to its regular value --*/
	tu /= 100;


	/*-- Calculate blending function f_theta --*/
	time_scale = 500.0*Laminar_Viscosity_i/(U_i[0]*Velocity_Mag*Velocity_Mag);

}


void CSourcePieceWise_TransLM::ComputeResidual_TransLM(double *val_residual, double **val_Jacobian_i, double &gamma_eff, CConfig *config, bool boundary, ofstream &sagt_debug) {

	//************************************************//
	// Please do not delete //SU2_CPP2C comment lines //
	//************************************************//

	//SU2_CPP2C START CSourcePieceWise_TransLM::ComputeResidual_TransLM
	//SU2_CPP2C CALL_LIST START
	//SU2_CPP2C INVARS *TransVar_i
	//SU2_CPP2C OUTVARS *val_residual
	//SU2_CPP2C VARS DOUBLE *U_i **PrimVar_Grad_i Laminar_Viscosity_i Eddy_Viscosity_i dist_i
	//SU2_CPP2C VARS DOUBLE SCALAR c_a1 c_e1 c_a2 c_e2 c_theta alpha_global flen_global Volume tu strain Velocity_Mag time_scale Vorticity re_theta_t
	//SU2_CPP2C CALL_LIST END

	//SU2_CPP2C DEFINE nDim

	//SU2_CPP2C DECL_LIST START
	//SU2_CPP2C VARS DOUBLE SCALAR Vorticity
	//SU2_CPP2C DECL_LIST END

	/*-- Local intermediate variables --*/
	double re_theta_c, flen, re_v, f_onset1, f_onset2, f_onset3, f_onset, f_turb;

	double prod, des;
	double r_t;
	double delta, theta, lambda, var1, f_theta;
	double theta_bl, f_reattach;
	double delta_bl, f_wake;
	double dU_dx, dU_dy, dU_dz;
	double gamma_sep;

	double fk[2], dx[2], J[2][2];
	double re_theta0, f_lambda0, dre_dlamb, normres, factor;
	double newton_tol = 1e-10;
	int iter;

	//SU2_CPP2C COMMENT START
	double val_residuald[2], TransVar_id[2];

	//SU2_CPP2C COMMENT END

	val_residual[0] = 0.0;
	val_residual[1] = 0.0;

	//SU2_CPP2C COMMENT START
	implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	if (implicit) {
		val_Jacobian_i[0][0] = 0.0;
		val_Jacobian_i[1][0] = 0.0;
		val_Jacobian_i[0][1] = 0.0;
		val_Jacobian_i[1][1] = 0.0;
	}

	/*-- Quit now if we're at the wall  --*/
	if (dist_i<1e-12) return;

	/* -- Compute intermediate correlations/expressions. These quantities which do not depend on TransVar are isolated to simplify the differentiated version of this routine.--*/
	translm_helper(config);
	//SU2_CPP2C COMMENT END


	/*-- Medida 2011, eq. 29-30 --*/
	re_theta_c = (4.45*pow(tu,3) - 5.7*pow(tu,2) + 1.37*tu + 0.585)*TransVar_i[1]/U_i[0];
	flen       = 0.171*pow(tu,2) - 0.0083*tu + 0.0306;

	re_v   = U_i[0]*pow(dist_i,2.)/Laminar_Viscosity_i*strain;  // Vorticity Reynolds number

	/*-- f_onset controls transition onset location --*/
	r_t      = Eddy_Viscosity_i/Laminar_Viscosity_i;
	f_onset1 = re_v / (2.193*re_theta_c);
	f_onset2 = min(max(f_onset1, pow(f_onset1,4.)), 2.);
	f_onset3 = max(1. - pow(0.4*r_t,3),0.);
	f_onset  = max(f_onset2 - f_onset3, 0.);

	f_turb = exp(-pow(0.25*r_t,4));  // Medida eq. 10

	prod = flen*c_a1*U_i[0]*strain*sqrt(f_onset*TransVar_i[0]/U_i[0]);
	prod = prod*(1. - c_e1*TransVar_i[0]/U_i[0]);

	des = c_a2*U_i[0]*Vorticity*TransVar_i[0]/U_i[0]*f_turb;
	des = des*(c_e2*TransVar_i[0]/U_i[0] - 1.);

	val_residual[0] = prod - des;
	val_residual[0] *= Volume;

	/*-- REtheta eq: --*/
	// Deactivated the f_wake parameter...
	theta_bl   = TransVar_i[1]/U_i[0]*Laminar_Viscosity_i / (U_i[0]*Velocity_Mag);
	delta_bl   = 7.5*theta_bl;
	delta      = 50.0*Vorticity*dist_i/Velocity_Mag*delta_bl + 1e-20;

	f_wake = 1.;

	var1 = (TransVar_i[0]/U_i[0]-1./c_e2)/(1.0-1./c_e2);
	var1 = 1. - pow(var1,2);

	f_theta = min(max(f_wake*exp(-pow(dist_i/delta,4)), var1),1.0);
	//f_theta = min(var1,1.0);

	val_residual[1] = c_theta*U_i[0]/time_scale *  (1.-f_theta) * (re_theta_t-TransVar_i[1]/U_i[0]);
	val_residual[1] *= Volume;

	//SU2_CPP2C COMMENT START
	sagt_debug << TransVar_i[0]/U_i[0] << " " << TransVar_i[1]/U_i[0] << " "
			<< re_theta_t << " " << flen << " " << re_theta_c << endl;

	/*-- Calculate term for separation correction --*/
	f_reattach = exp(-pow(0.05*r_t,4));
	gamma_sep = s1*max(0.,re_v/(3.235*re_theta_c)-1.)*f_reattach;
	gamma_sep = min(gamma_sep,2.0)*f_theta;
	gamma_eff = max(gamma_sep, TransVar_i[0]/U_i[0]);

	/*--- Implicit part ---*/
	TransVar_id[0] = 1.0; TransVar_id[1] = 0.0;
	CSourcePieceWise_TransLM__ComputeResidual_TransLM_d(TransVar_i, TransVar_id, val_residual, val_residuald, config, boundary);
	val_Jacobian_i[0][0] = val_residuald[0];
	val_Jacobian_i[1][0] = val_residuald[1];

	TransVar_id[0] = 0.0; TransVar_id[1] = 1.0;
	CSourcePieceWise_TransLM__ComputeResidual_TransLM_d(TransVar_i, TransVar_id, val_residual, val_residuald, config, boundary);
	val_Jacobian_i[0][1] = val_residuald[0];
	val_Jacobian_i[1][1] = val_residuald[1];

	//SU2_CPP2C COMMENT END

	//SU2_CPP2C END CSourcePieceWise_TransLM::ComputeResidual_TransLM
}


void CSourcePieceWise_TransLM::CSourcePieceWise_TransLM__ComputeResidual_TransLM_d(double *TransVar_i, double *TransVar_id, double *val_residual, double *val_residuald, CConfig *config, bool boundary)
{
	double re_theta_c, flen, re_v, f_onset1, f_onset2, f_onset3, f_onset,
	f_turb;
	double re_theta_cd, f_onset1d, f_onset2d, f_onsetd;
	double prod, des;
	double prodd, desd;
	double r_t;
	double delta, theta, lambda, var1, f_theta;
	double deltad, var1d, f_thetad;
	double theta_bl, f_reattach;
	double theta_bld;
	double delta_bl, f_wake;
	double delta_bld;
	double dU_dx, dU_dy, dU_dz;
	double fk[2], dx[2], J[2][2];
	double re_theta0, f_lambda0, dre_dlamb, normres, factor;
	double newton_tol = 1e-10;
	int iter;
	double result1;
	double result1d;
	double result2;
	double arg1;
	double arg1d;
	float result10;
	float result10d;
	double x4;
	double x3;
	double x2;
	double x1;
	double x1d;
	double x4d;
	double x3d;
	double y1;
	double y1d;
	val_residuald[0] = 0.0;
	val_residual[0] = 0.0;
	val_residuald[1] = 0.0;
	val_residual[1] = 0.0;
	/*-- Medida 2011, eq. 29-30 --*/
	result1 = pow(tu, 3);
	result2 = pow(tu, 2);
	re_theta_cd = (4.45*result1-5.7*result2+1.37*tu+0.585)*TransVar_id[1]/U_i[
	                                                                          0];
	re_theta_c = (4.45*result1-5.7*result2+1.37*tu+0.585)*TransVar_i[1]/U_i[0]
	                                                                        ;
	result1 = pow(tu, 2);
	flen = 0.171*result1 - 0.0083*tu + 0.0306;
	result1 = pow(dist_i, 2.);
	re_v = U_i[0]*result1/Laminar_Viscosity_i*strain;
	/*-- f_onset controls transition onset location --*/
	// Vorticity Reynolds number
	r_t = Eddy_Viscosity_i/Laminar_Viscosity_i;
	f_onset1d = -(re_v*2.193*re_theta_cd/(2.193*re_theta_c*(2.193*re_theta_c))
	);
	f_onset1 = re_v/(2.193*re_theta_c);
	y1d = pow_d(f_onset1, f_onset1d, 4., &y1);
	if (f_onset1 < y1) {
		x1d = y1d;
		x1 = y1;
	} else {
		x1d = f_onset1d;
		x1 = f_onset1;
	}
	if (x1 > 2.) {
		f_onset2 = 2.;
		f_onset2d = 0.0;
	} else {
		f_onset2d = x1d;
		f_onset2 = x1;
	}
	result1 = pow(0.4*r_t, 3);
	x2 = 1. - result1;
	if (x2 < 0.)
		f_onset3 = 0.;
	else
		f_onset3 = x2;
	if (f_onset2 - f_onset3 < 0.) {
		f_onset = 0.;
		f_onsetd = 0.0;
	} else {
		f_onsetd = f_onset2d;
		f_onset = f_onset2 - f_onset3;
	}
	result1 = pow(0.25*r_t, 4);
	f_turb = exp(-result1);
	// Medida eq. 10
	arg1d = (f_onsetd*TransVar_i[0]+f_onset*TransVar_id[0])/U_i[0];
	arg1 = f_onset*TransVar_i[0]/U_i[0];
	result10d = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1)));
	result10 = sqrt(arg1);
	prodd = flen*c_a1*U_i[0]*strain*result10d;
	prod = flen*c_a1*U_i[0]*strain*result10;
	prodd = prodd*(1.-c_e1*TransVar_i[0]/U_i[0]) - prod*c_e1*TransVar_id[0]/
			U_i[0];
	prod = prod*(1.-c_e1*TransVar_i[0]/U_i[0]);
	desd = f_turb*c_a2*Vorticity*TransVar_id[0];
	des = c_a2*U_i[0]*Vorticity*TransVar_i[0]/U_i[0]*f_turb;
	desd = desd*(c_e2*TransVar_i[0]/U_i[0]-1.) + des*c_e2*TransVar_id[0]/U_i[0
	                                                                         ];
	des = des*(c_e2*TransVar_i[0]/U_i[0]-1.);
	val_residuald[0] = prodd - desd;
	val_residual[0] = prod - des;
	val_residuald[0] = Volume*val_residuald[0];
	val_residual[0] *= Volume;
	/*-- REtheta eq: --*/
	// Deactivated the f_wake parameter...
	theta_bld = Laminar_Viscosity_i*TransVar_id[1]/U_i[0]/(U_i[0]*Velocity_Mag
	);
	theta_bl = TransVar_i[1]/U_i[0]*Laminar_Viscosity_i/(U_i[0]*Velocity_Mag);
	delta_bld = 7.5*theta_bld;
	delta_bl = 7.5*theta_bl;
	deltad = 50.0*Vorticity*dist_i*delta_bld/Velocity_Mag;
	delta = 50.0*Vorticity*dist_i/Velocity_Mag*delta_bl + 1e-20;
	f_wake = 1.;
	var1d = TransVar_id[0]/U_i[0]/(1.0-1./c_e2);
	var1 = (TransVar_i[0]/U_i[0]-1./c_e2)/(1.0-1./c_e2);
	result1d = pow_d(var1, var1d, 2, &result1);
	var1d = -result1d;
	var1 = 1. - result1;
	result1d = pow_d(dist_i/delta, -(dist_i*deltad/(delta*delta)), 4, &result1
	);
	x4d = -(f_wake*result1d*exp(-result1));
	x4 = f_wake*exp(-result1);
	if (x4 < var1) {
		x3d = var1d;
		x3 = var1;
	} else {
		x3d = x4d;
		x3 = x4;
	}
	if (x3 > 1.0) {
		f_theta = 1.0;
		f_thetad = 0.0;
	} else {
		f_thetad = x3d;
		f_theta = x3;
	}
	//f_theta = min(var1,1.0);
	val_residuald[1] = c_theta*U_i[0]*(-(f_thetad*(re_theta_t-TransVar_i[1]/U_i[
	                                                                          0]))-(1.-f_theta)*TransVar_id[1]/U_i[0])/time_scale;
	val_residual[1] = c_theta*U_i[0]/time_scale*(1.-f_theta)*(re_theta_t-
			TransVar_i[1]/U_i[0]);
	val_residuald[1] = Volume*val_residuald[1];
	val_residual[1] *= Volume;


}
