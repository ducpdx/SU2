/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

#include "../include/solver_structure.hpp"

CTransLMSolver::CTransLMSolver(void) : CTurbSolver() {}

CTransLMSolver::CTransLMSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CTurbSolver() {
	unsigned short iVar, iDim;
	unsigned long iPoint, index;
	double Viscosity_Inf, tu_Inf, dull_val, rey, mach;
	ifstream restart_file;
	char *cstr;
	string text_line;
	bool restart = (config->GetRestart() || config->GetRestart_Flow());
	bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
	bool freesurface = (config->GetKind_Regime() == FREESURFACE);
    bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
    bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
	
  cout << "Entered constructor for CTransLMSolver -AA\n";
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	
	/*--- Define geometry constans in the solver structure ---*/
	nDim = geometry->GetnDim();
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
	node = new CVariable*[geometry->GetnPoint()];
	
	/*--- Dimension of the problem --> 2 Transport equations (intermittency,Reth) ---*/
	nVar = 2;
	
	if (iMesh == MESH_0) {

		/*--- Define some auxillary vectors related to the residual ---*/
		Residual     = new double[nVar]; Residual_RMS = new double[nVar];
		Residual_i   = new double[nVar]; Residual_j   = new double[nVar];
    Residual_Max = new double[nVar]; Point_Max    = new unsigned long[nVar];

		/*--- Define some auxiliar vector related with the solution ---*/
		Solution   = new double[nVar];
		Solution_i = new double[nVar]; Solution_j = new double[nVar];
		
		/*--- Define some auxiliar vector related with the geometry ---*/
		Vector_i = new double[nDim]; Vector_j = new double[nDim];
				
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    
		/*--- Jacobians and vector structures for implicit computations ---*/
		if (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT) {

			/*--- Point to point Jacobians ---*/
			Jacobian_i = new double* [nVar];
			Jacobian_j = new double* [nVar];
			for (iVar = 0; iVar < nVar; iVar++) {
				Jacobian_i[iVar] = new double [nVar];
				Jacobian_j[iVar] = new double [nVar];
			}
			/*--- Initialization of the structure of the whole Jacobian ---*/
			Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry);
      
		}
	
	/*--- Computation of gradients by least squares ---*/
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		Smatrix = new double* [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Smatrix[iDim] = new double [nDim];
		/*--- c vector := transpose(WA)*(Wb) ---*/
		cvector = new double* [nVar];
		for (iVar = 0; iVar < nVar; iVar++)
			cvector[iVar] = new double [nDim];
	}

	/*--- Read farfield conditions from config ---*/
	Density_Inf       = config->GetDensity_FreeStreamND();
	Viscosity_Inf     = config->GetViscosity_FreeStreamND();
	Intermittency_Inf = config->GetIntermittency_FreeStream();
	tu_Inf            = config->GetTurbulenceIntensity_FreeStream()*100;
	
  /*-- Initialize REth from correlation --*/
  if (tu_Inf <= 1.3) {
    REth_Inf = (1173.51-589.428*tu_Inf+0.2196/(tu_Inf*tu_Inf));
  } else {
    REth_Inf = 331.5*pow(tu_Inf-0.5658,-0.671);
  }
  rey = config->GetReynolds();
  mach = config->GetMach_FreeStreamND();

//  REth_Inf *= mach/rey;
  cout << "REth_Inf = " << REth_Inf << ", rey: "<< rey << " -AA" << endl;
	
	/*--- Restart the solution from file information ---*/
	if (!restart) {
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint] = new CTransLMVariable(Density_Inf, Intermittency_Inf, REth_Inf, nDim, nVar, config);

	}
	else {

		/*--- Restart the solution from file information ---*/
		ifstream restart_file;
		string filename = config->GetSolution_FlowFileName();

        /*--- Modify file name for an unsteady restart ---*/
        if (dual_time) {
            int Unst_RestartIter;
            if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
                Unst_RestartIter = int(config->GetUnst_RestartIter())-1;
            else
                Unst_RestartIter = int(config->GetUnst_RestartIter())-2;
            filename = config->GetUnsteady_FileName(filename, Unst_RestartIter);
        }

        /*--- Open the restart file, throw an error if this fails. ---*/
		restart_file.open(filename.data(), ios::in);
		if (restart_file.fail()) {
			cout << "There is no turbulent restart file!!" << endl;
			exit(1);
		}

		/*--- In case this is a parallel simulation, we need to perform the
         Global2Local index transformation first. ---*/
		long *Global2Local;
		Global2Local = new long[geometry->GetGlobal_nPointDomain()];
		/*--- First, set all indices to a negative value by default ---*/
		for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
			Global2Local[iPoint] = -1;
		}
		/*--- Now fill array with the transform values only for local points ---*/
		for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
			Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
		}

		/*--- Read all lines in the restart file ---*/
		long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;

        /*--- The first line is the header ---*/
        getline (restart_file, text_line);

		while (getline (restart_file,text_line)) {
			istringstream point_line(text_line);

			/*--- Retrieve local index. If this node from the restart file lives
             on a different processor, the value of iPoint_Local will be -1.
             Otherwise, the local index for this node on the current processor
             will be returned and used to instantiate the vars. ---*/
			iPoint_Local = Global2Local[iPoint_Global];

			if (iPoint_Local >= 0) {

				double Density;
				if (compressible) {
					if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Density >> dull_val>> dull_val>> dull_val >> dull_val >> Solution[0] >> Solution[1];
					if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Density >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];

				}
				if (incompressible) {
					if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
					if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
				}

                if (freesurface) {
					if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1];
					if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >>  dull_val >> dull_val >> Solution[0] >> Solution[1];
				}

				/*--- Instantiate the solution at this node, note that the eddy viscosity should be recomputed ---*/
    			node[iPoint_Local] = new CTransLMVariable(Density, Solution[0]/Density, Solution[1]/Density, nDim, nVar, config);

    			/*--- Compute gamma_eff, which will be needed in first iteration of turbulence model ---*/
    			node[iPoint_Local]->SetGammaEff(1.0); // TODO: Need to recompute gamma_eff for restart
			}
			iPoint_Global++;
		}

		/*--- Instantiate the variable class with an arbitrary solution
         at any halo/periodic nodes. The initial solution can be arbitrary,
         because a send/recv is performed immediately in the solver. ---*/
		for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CTransLMVariable(Density_Inf, Intermittency_Inf, REth_Inf, nDim, nVar, config);
		}

		/*--- Close the restart file ---*/
		restart_file.close();

		/*--- Free memory needed for the transformation ---*/
		delete [] Global2Local;
	}

    /*--- MPI solution ---*/
    Set_MPI_Solution(geometry, config);

}

CTransLMSolver::~CTransLMSolver(void){
	unsigned short iVar, iDim;
	
	delete [] Residual; delete [] Residual_Max;
	delete [] Residual_i; delete [] Residual_j;
	delete [] Solution;
	delete [] Solution_i; delete [] Solution_j;
	delete [] Vector_i; delete [] Vector_j;
	
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] Jacobian_i[iVar];
		delete [] Jacobian_j[iVar];
	}
	delete [] Jacobian_i; delete [] Jacobian_j;
	
	for (iDim = 0; iDim < this->nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] Smatrix;
	
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
}

void CTransLMSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector, nBufferS_Scalar, nBufferR_Scalar;
	double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL, *Buffer_Receive_gamma_eff = NULL, *Buffer_Send_gamma_eff = NULL;
	int send_to, receive_from;

#ifndef NO_MPI
    MPI::Status status;
    MPI::Request send_request, recv_request;
#endif

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
            (config->GetMarker_All_SendRecv(iMarker) > 0)) {

			MarkerS = iMarker;  MarkerR = iMarker+1;

            send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
            nBufferS_Scalar = nVertexS;             nBufferR_Scalar = nVertexR;

            /*--- Allocate Receive and send buffers  ---*/
            Buffer_Receive_U = new double [nBufferR_Vector];
            Buffer_Send_U = new double[nBufferS_Vector];

            Buffer_Receive_gamma_eff = new double [nBufferR_Scalar];
            Buffer_Send_gamma_eff = new double[nBufferS_Scalar];

            /*--- Copy the solution that should be sended ---*/
            for (iVertex = 0; iVertex < nVertexS; iVertex++) {
                iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
                Buffer_Send_gamma_eff[iVertex] = node[iPoint]->GetGammaEff();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
            }

#ifndef NO_MPI

            /*--- Send/Receive information using Sendrecv ---*/
            MPI::COMM_WORLD.Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                                     Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);
            MPI::COMM_WORLD.Sendrecv(Buffer_Send_gamma_eff, nBufferS_Scalar, MPI::DOUBLE, send_to, 1,
                                     Buffer_Receive_gamma_eff, nBufferR_Scalar, MPI::DOUBLE, receive_from, 1);

#else

            /*--- Receive information without MPI ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                Buffer_Receive_gamma_eff[iVertex] = node[iPoint]->GetGammaEff();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
            }

#endif

            /*--- Deallocate send buffer ---*/
            delete [] Buffer_Send_U;
            delete [] Buffer_Send_gamma_eff;

            /*--- Do the coordinate transformation ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {

                /*--- Find point and its type of transformation ---*/
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

                /*--- Copy conservative variables. ---*/
                node[iPoint]->SetGammaEff(Buffer_Receive_gamma_eff[iVertex]);
                for (iVar = 0; iVar < nVar; iVar++)
                    node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);

            }

            /*--- Deallocate receive buffer ---*/
            delete [] Buffer_Receive_gamma_eff;
            delete [] Buffer_Receive_U;

        }

	}

}

void CTransLMSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
	int send_to, receive_from;

#ifndef NO_MPI
    MPI::Status status;
    MPI::Request send_request, recv_request;
#endif

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
            (config->GetMarker_All_SendRecv(iMarker) > 0)) {

			MarkerS = iMarker;  MarkerR = iMarker+1;

            send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

            /*--- Allocate Receive and send buffers  ---*/
            Buffer_Receive_U = new double [nBufferR_Vector];
            Buffer_Send_U = new double[nBufferS_Vector];

            /*--- Copy the solution old that should be sended ---*/
            for (iVertex = 0; iVertex < nVertexS; iVertex++) {
                iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
            }

#ifndef NO_MPI

            //      /*--- Send/Receive using non-blocking communications ---*/
            //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
            //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
            //      send_request.Wait(status);
            //      recv_request.Wait(status);

            /*--- Send/Receive information using Sendrecv ---*/
            MPI::COMM_WORLD.Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                                     Buffer_Receive_U, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);

#else

            /*--- Receive information without MPI ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
            }

#endif

            /*--- Deallocate send buffer ---*/
            delete [] Buffer_Send_U;

            /*--- Do the coordinate transformation ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {

                /*--- Find point and its type of transformation ---*/
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

                /*--- Copy transformed conserved variables back into buffer. ---*/
                for (iVar = 0; iVar < nVar; iVar++)
                    node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);

            }

            /*--- Deallocate receive buffer ---*/
            delete [] Buffer_Receive_U;

        }

	}
}

void CTransLMSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
    *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
	int send_to, receive_from;

#ifndef NO_MPI
    MPI::Status status;
    MPI::Request send_request, recv_request;
#endif

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
            (config->GetMarker_All_SendRecv(iMarker) > 0)) {

			MarkerS = iMarker;  MarkerR = iMarker+1;

            send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

            /*--- Allocate Receive and send buffers  ---*/
            Buffer_Receive_Limit = new double [nBufferR_Vector];
            Buffer_Send_Limit = new double[nBufferS_Vector];

            /*--- Copy the solution old that should be sended ---*/
            for (iVertex = 0; iVertex < nVertexS; iVertex++) {
                iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
            }

#ifndef NO_MPI

            //      /*--- Send/Receive using non-blocking communications ---*/
            //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Limit, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
            //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Limit, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
            //      send_request.Wait(status);
            //      recv_request.Wait(status);

            /*--- Send/Receive information using Sendrecv ---*/
            MPI::COMM_WORLD.Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                                     Buffer_Receive_Limit, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);

#else

            /*--- Receive information without MPI ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
            }

#endif

            /*--- Deallocate send buffer ---*/
            delete [] Buffer_Send_Limit;

            /*--- Do the coordinate transformation ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {

                /*--- Find point and its type of transformation ---*/
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

                /*--- Retrieve the supplied periodic information. ---*/
                angles = config->GetPeriodicRotation(iPeriodic_Index);

                /*--- Store angles separately for clarity. ---*/
                theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
                cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
                sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

                /*--- Compute the rotation matrix. Note that the implicit
                 ordering is rotation about the x-axis, y-axis,
                 then z-axis. Note that this is the transpose of the matrix
                 used during the preprocessing stage. ---*/
                rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
                rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
                rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

                /*--- Copy conserved variables before performing transformation. ---*/
                for (iVar = 0; iVar < nVar; iVar++)
                    Solution[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];

                /*--- Rotate the momentum components. ---*/
                if (nDim == 2) {
                    Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
                    rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
                    Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
                    rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
                }
                else {
                    Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
                    rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
                    rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
                    Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
                    rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
                    rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
                    Solution[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
                    rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
                    rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
                }

                /*--- Copy transformed conserved variables back into buffer. ---*/
                for (iVar = 0; iVar < nVar; iVar++)
                    node[iPoint]->SetLimiter(iVar, Solution[iVar]);

            }

            /*--- Deallocate receive buffer ---*/
            delete [] Buffer_Receive_Limit;

        }

	}
}

void CTransLMSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
    *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
	int send_to, receive_from;

#ifndef NO_MPI
    MPI::Status status;
    MPI::Request send_request, recv_request;
#endif

    double **Gradient = new double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
        Gradient[iVar] = new double[nDim];

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		if ((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) &&
            (config->GetMarker_All_SendRecv(iMarker) > 0)) {

			MarkerS = iMarker;  MarkerR = iMarker+1;

            send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
			receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;

			nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
			nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;

            /*--- Allocate Receive and send buffers  ---*/
            Buffer_Receive_Gradient = new double [nBufferR_Vector];
            Buffer_Send_Gradient = new double[nBufferS_Vector];

            /*--- Copy the solution old that should be sended ---*/
            for (iVertex = 0; iVertex < nVertexS; iVertex++) {
                iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    for (iDim = 0; iDim < nDim; iDim++)
                        Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
            }

#ifndef NO_MPI

            //      /*--- Send/Receive using non-blocking communications ---*/
            //      send_request = MPI::COMM_WORLD.Isend(Buffer_Send_Gradient, nBufferS_Vector, MPI::DOUBLE, 0, send_to);
            //      recv_request = MPI::COMM_WORLD.Irecv(Buffer_Receive_Gradient, nBufferR_Vector, MPI::DOUBLE, 0, receive_from);
            //      send_request.Wait(status);
            //      recv_request.Wait(status);

            /*--- Send/Receive information using Sendrecv ---*/
            MPI::COMM_WORLD.Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI::DOUBLE, send_to, 0,
                                     Buffer_Receive_Gradient, nBufferR_Vector, MPI::DOUBLE, receive_from, 0);

#else

            /*--- Receive information without MPI ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                for (iVar = 0; iVar < nVar; iVar++)
                    for (iDim = 0; iDim < nDim; iDim++)
                        Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
            }

#endif

            /*--- Deallocate send buffer ---*/
            delete [] Buffer_Send_Gradient;

            /*--- Do the coordinate transformation ---*/
            for (iVertex = 0; iVertex < nVertexR; iVertex++) {

                /*--- Find point and its type of transformation ---*/
                iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
                iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

                /*--- Retrieve the supplied periodic information. ---*/
                angles = config->GetPeriodicRotation(iPeriodic_Index);

                /*--- Store angles separately for clarity. ---*/
                theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
                cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
                sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

                /*--- Compute the rotation matrix. Note that the implicit
                 ordering is rotation about the x-axis, y-axis,
                 then z-axis. Note that this is the transpose of the matrix
                 used during the preprocessing stage. ---*/
                rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
                rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
                rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

                /*--- Copy conserved variables before performing transformation. ---*/
                for (iVar = 0; iVar < nVar; iVar++)
                    for (iDim = 0; iDim < nDim; iDim++)
                        Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];

                /*--- Need to rotate the gradients for all conserved variables. ---*/
                for (iVar = 0; iVar < nVar; iVar++) {
                    if (nDim == 2) {
                        Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
                        Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
                    }
                    else {
                        Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
                        Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
                        Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
                    }
                }

                /*--- Store the received information ---*/
                for (iVar = 0; iVar < nVar; iVar++)
                    for (iDim = 0; iDim < nDim; iDim++)
                        node[iPoint]->SetGradient(iVar, iDim, Gradient[iVar][iDim]);

            }

            /*--- Deallocate receive buffer ---*/
            delete [] Buffer_Receive_Gradient;

        }

	}

    for (iVar = 0; iVar < nVar; iVar++)
        delete [] Gradient[iVar];
    delete [] Gradient;

}

void CTransLMSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem) {
	unsigned long iPoint;

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		LinSysRes.SetBlock_Zero(iPoint);
  Jacobian.SetValZero();

	if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
	if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
}

void CTransLMSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

}

void CTransLMSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
	unsigned short iVar;
	unsigned long iPoint, total_index;
	double Delta, Delta_flow, Vol;
    
    
    /*--- Set maximum residual to zero ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
    }
    
	/*--- Build implicit system ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
        Vol = geometry->node[iPoint]->GetVolume();
        
        /*--- Modify matrix diagonal to assure diagonal dominance ---*/
        Delta_flow = Vol / (solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
        Delta = Delta_flow;
        Jacobian.AddVal2Diag(iPoint,Delta);
        
        for (iVar = 0; iVar < nVar; iVar++) {
            total_index = iPoint*nVar+iVar;
            
            /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
            LinSysRes[total_index] = -LinSysRes[total_index];
            LinSysSol[total_index] = 0.0;
            AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]*Vol);
            AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex());
        }
    }
    
    /*--- Initialize residual and solution at the ghost points ---*/
    for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
            total_index = iPoint*nVar + iVar;
            LinSysRes[total_index] = 0.0;
            LinSysSol[total_index] = 0.0;
        }
    }
	
	/*--- Solve the linear system (Krylov subspace methods) ---*/
  CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(Jacobian, geometry, config);
  
  CPreconditioner* precond = NULL;
  if (config->GetKind_Linear_Solver_Prec() == JACOBI) {
    Jacobian.BuildJacobiPreconditioner();
    precond = new CJacobiPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LU_SGS) {
    precond = new CLU_SGSPreconditioner(Jacobian, geometry, config);
  }
  else if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    Jacobian.BuildJacobiPreconditioner();
    Jacobian.BuildLineletPreconditioner(geometry, config);
    precond = new CLineletPreconditioner(Jacobian, geometry, config);
  }
  
  CSysSolve system;
  if (config->GetKind_Linear_Solver() == BCGSTAB)
    system.BCGSTAB(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                   config->GetLinear_Solver_Iter(), false);
  else if (config->GetKind_Linear_Solver() == FGMRES)
    system.FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, config->GetLinear_Solver_Error(),
                 config->GetLinear_Solver_Iter(), false);
  
  delete mat_vec;
  delete precond;
  
	/*--- Update solution (system written in terms of increments) ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++)
            node[iPoint]->AddSolution(iVar, config->GetLinear_Solver_Relax()*LinSysSol[iPoint*nVar+iVar]);
    }
    
    /*--- MPI solution ---*/
    Set_MPI_Solution(geometry, config);
    
    /*--- Compute the root mean square residual ---*/
    SetResidual_RMS(geometry, config);

}

void CTransLMSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {

  double *trans_var_i, *trans_var_j, *U_i, *U_j;
  double **Gradient_i, **Gradient_j, *Limiter_i = NULL, *Limiter_j = NULL, Project_Grad_i, Project_Grad_j;
  unsigned long iEdge, iPoint, jPoint, iDim, iVar;
	bool high_order_diss = (config->GetKind_Upwind_Turb() == SCALAR_UPWIND_2ND);
	bool limiter = (config->GetKind_SlopeLimit() != NONE);


	if (high_order_diss) {
		if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solver_container[FLOW_SOL]->SetSolution_Gradient_GG(geometry, config);
		if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solver_container[FLOW_SOL]->SetSolution_Gradient_LS(geometry, config);
	}

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge and normal vectors ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Conservative variables w/o reconstruction ---*/
    U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
    U_j = solver_container[FLOW_SOL]->node[jPoint]->GetSolution();
    numerics->SetConservative(U_i, U_j);

    /*--- Transition variables w/o reconstruction ---*/
    trans_var_i = node[iPoint]->GetSolution();
    trans_var_j = node[jPoint]->GetSolution();
    numerics->SetTransVar(trans_var_i,trans_var_j);

		if (high_order_diss) {
            
			/*--- Conservative solution using gradient reconstruction ---*/
			for (iDim = 0; iDim < nDim; iDim++) {
				Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
				Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
			}
			Gradient_i = solver_container[FLOW_SOL]->node[iPoint]->GetGradient();
			Gradient_j = solver_container[FLOW_SOL]->node[jPoint]->GetGradient();
            
			for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				FlowSolution_i[iVar] = U_i[iVar] + Project_Grad_i;
				FlowSolution_j[iVar] = U_j[iVar] + Project_Grad_j;
			}
			numerics->SetConservative(FlowSolution_i, FlowSolution_j);
            
			/*--- Transition variables using gradient reconstruction ---*/
			Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
			if (limiter) { Limiter_i = node[iPoint]->GetLimiter(); Limiter_j = node[jPoint]->GetLimiter(); }
            
			for (iVar = 0; iVar < nVar; iVar++) {
				Project_Grad_i = 0; Project_Grad_j = 0;
				for (iDim = 0; iDim < nDim; iDim++) {
					Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
					Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
				}
				if (limiter) {
					Solution_i[iVar] = trans_var_i[iVar] + Project_Grad_i*Limiter_i[iVar];
					Solution_j[iVar] = trans_var_j[iVar] + Project_Grad_j*Limiter_j[iVar];
				}
				else {
					Solution_i[iVar] = trans_var_i[iVar] + Project_Grad_i;
					Solution_j[iVar] = trans_var_j[iVar] + Project_Grad_j;
				}
			}
			numerics->SetTransVar(Solution_i, Solution_j);
		}

    /*--- Add and subtract Residual ---*/
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    LinSysRes.AddBlock(iPoint, Residual);
    LinSysRes.SubtractBlock(jPoint, Residual);

    /*--- Implicit part ---*/
    Jacobian.AddBlock(iPoint,iPoint,Jacobian_i);
    Jacobian.AddBlock(iPoint,jPoint,Jacobian_j);
    Jacobian.SubtractBlock(jPoint,iPoint,Jacobian_i);
    Jacobian.SubtractBlock(jPoint,jPoint,Jacobian_j);

  }

}


void CTransLMSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
																				CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
	unsigned long iEdge, iPoint, jPoint;
	
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Points coordinates, and normal vector ---*/
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                     geometry->node[jPoint]->GetCoord());

    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Conservative variables w/o reconstruction ---*/
    numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                            solver_container[FLOW_SOL]->node[jPoint]->GetSolution());

    /*--- Laminar Viscosity ---*/
    numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
    /*--- Eddy Viscosity ---*/
    numerics->SetEddyViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
                             solver_container[FLOW_SOL]->node[jPoint]->GetEddyViscosity());

    /*--- Transition variables w/o reconstruction, and its gradients ---*/
    numerics->SetTransVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    numerics->SetTransVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());

    numerics->SetConsVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient(),
                               solver_container[FLOW_SOL]->node[jPoint]->GetGradient());


    /*--- Compute residual, and Jacobians ---*/
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

    /*--- Add and subtract residual, and update Jacobians ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint,iPoint,Jacobian_i);
    Jacobian.SubtractBlock(iPoint,jPoint,Jacobian_j);
    Jacobian.AddBlock(jPoint,iPoint,Jacobian_i);
    Jacobian.AddBlock(jPoint,jPoint,Jacobian_j);
    
  }
}

void CTransLMSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {
  unsigned long iPoint, iVertex;
  unsigned short iMarker;
  double gamma_eff;
  bool boundary;

    //cout << "Setting Trans residual -AA " << endl;
    //cout << "\nBeginAA" << endl;

  // DEBUG
//  sagt_debug.open("sagt_debug.plt");
//  sagt_debug << "TITLE = \"SAGT (Langtry+Menter) Transition model debug file \" " << endl;
//  sagt_debug << "VARIABLES = \"itmc\" \"Re_th_bar\" \"re_theta_t\" \"flen\" \"re_theta_c\"" << endl;
//  sagt_debug << "ZONE DATAPACKING=POINT" << endl;

  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
	  //   cout << "\niPoint: " << iPoint << endl;

	  /*--- Conservative variables w/o reconstruction ---*/
	  numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(), NULL);

	  /*--- Gradient of the primitive and conservative variables ---*/
	  numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);

	  /*--- Laminar and eddy viscosity ---*/
	  numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
	  numerics->SetEddyViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),0.0);

	  /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
	  numerics->SetTransVar(node[iPoint]->GetSolution(), NULL);

	  /*--- Set volume ---*/
	  numerics->SetVolume(geometry->node[iPoint]->GetVolume());

	  /*--- Set distance to the surface ---*/
	  numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);

	  /*--- Set distance to the surface ---*/
	  boundary = geometry->node[iPoint]->GetBoundary();

	  /*--- Compute the source term ---*/
	  numerics->ComputeResidual_TransLM(Residual, Jacobian_i, gamma_eff, config, boundary, sagt_debug);

	  /*-- Store gamma_eff in variable class, where CTurbSASolver can access it --*/
	  node[iPoint]->SetGammaEff(gamma_eff);

	  /*--- Subtract residual and the jacobian ---*/
	  LinSysRes.SubtractBlock(iPoint, Residual);
	  Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

  }
//  sagt_debug.close();
  
}

void CTransLMSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
																			 CConfig *config, unsigned short iMesh) {
}

void CTransLMSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iPoint, iVertex;
  unsigned short iVar;
  int total_index;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /* --- Impose boundary values (Dirichlet) ---*/
      Solution[0] = 0.0;
      Solution[1] = 0.0;
      node[iPoint]->SetSolution_Old(Solution);
      LinSysRes.SetBlock_Zero(iPoint);

      /*--- includes 1 in the diagonal ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    }
  }


}

void CTransLMSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	unsigned long iPoint, iVertex;
  unsigned long total_index;
  unsigned short iVar;

	//cout << "Arrived in BC_Far_Field. -AA" << endl;
	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
		if (geometry->node[iPoint]->GetDomain()) {

			/* --- Impose boundary values (Dirichlet) ---*/
			Solution[0] = Intermittency_Inf*Density_Inf;
			Solution[1] = REth_Inf*Density_Inf;
			node[iPoint]->SetSolution_Old(Solution);
			LinSysRes.SetBlock_Zero(iPoint);

			/*--- includes 1 in the diagonal ---*/
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				Jacobian.DeleteValsRowi(total_index);
			}
		}
	}

}

void CTransLMSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

unsigned short iVar, iDim, Kind_Inlet = config->GetKind_Inlet();
  unsigned long iVertex, iPoint, Point_Normal;
  double P_Total, T_Total, Velocity[3];
  double Velocity2, H_Total, Temperature, Riemann;
  double Pressure, Density, Energy, *Flow_Dir, Mach2;
  double SoundSpeed2, SoundSpeed_Total2, Vel_Mag;
  double alpha, aa, bb, cc, dd;
  double Two_Gamma_M1 = 2.0/Gamma_Minus_One;
  double Gas_Constant = config->GetGas_ConstantND();
  double *Normal = new double[nDim];
  
  bool rotating_frame = config->GetRotating_Frame();
  bool grid_movement  = config->GetGrid_Movement();
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  string Marker_Tag = config->GetMarker_All_Tag(val_marker);
  
  double *Coord;
  unsigned long total_index;
  double density;

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      double Area = 0.0; double UnitaryNormal[3];
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitaryNormal[iDim] = Normal[iDim]/Area;

      /*--- Current conservative variables at this boundary node (U_domain) ---*/
      for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++)
        FlowSolution_i[iVar] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(iVar);

      /*--- Build the fictitious intlet state based on characteristics ---*/
      if (incompressible) {

        /*--- Pressure computation using the internal value ---*/
        FlowSolution_j[0] = solver_container[FLOW_SOL]->node[iPoint]->GetSolution(0);

        /*--- The velocity is computed from the interior, and normal
         derivative for the density ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          FlowSolution_j[iDim+1] = solver_container[FLOW_SOL]->GetVelocity_Inf(iDim)*solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc();

      }
      else {

        /*--- Subsonic inflow: there is one outgoing characteristic (u-c),
         therefore we can specify all but one state variable at the inlet.
         The outgoing Riemann invariant provides the final piece of info.
         Adapted from an original implementation in the Stanford University
         multi-block (SUmb) solver in the routine bcSubsonicInflow.f90
         written by Edwin van der Weide, last modified 04-20-2009. ---*/

        switch (Kind_Inlet) {

            /*--- Total properties have been specified at the inlet. ---*/
          case TOTAL_CONDITIONS:

            /*--- Retrieve the specified total conditions for this inlet. ---*/
            P_Total  = config->GetInlet_Ptotal(Marker_Tag);
            T_Total  = config->GetInlet_Ttotal(Marker_Tag);
            Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

            /*--- Non-dim. the inputs if necessary. ---*/
            P_Total /= config->GetPressure_Ref();
            T_Total /= config->GetTemperature_Ref();

            /*--- Store primitives and set some variables for clarity. ---*/
            Density = FlowSolution_i[0];
            Velocity2 = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              Velocity[iDim] = FlowSolution_i[iDim+1]/Density;
              Velocity2 += Velocity[iDim]*Velocity[iDim];
            }
            Energy      = FlowSolution_i[nVar-1]/Density;
            Pressure    = Gamma_Minus_One*Density*(Energy-0.5*Velocity2);
            H_Total     = (Gamma*Gas_Constant/Gamma_Minus_One)*T_Total;
            SoundSpeed2 = Gamma*Pressure/Density;

            /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/
            Riemann   = 2.0*sqrt(SoundSpeed2)/Gamma_Minus_One;
            for (iDim = 0; iDim < nDim; iDim++)
              Riemann += Velocity[iDim]*UnitaryNormal[iDim];

            /*--- Total speed of sound ---*/
            SoundSpeed_Total2 = Gamma_Minus_One*(H_Total - (Energy
                                                            + Pressure/Density)+0.5*Velocity2) + SoundSpeed2;

            /*--- Dot product of normal and flow direction. This should
             be negative due to outward facing boundary normal convention. ---*/
            alpha = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              alpha += UnitaryNormal[iDim]*Flow_Dir[iDim];

            /*--- Coefficients in the quadratic equation for the velocity ---*/
            aa =  1.0 + 0.5*Gamma_Minus_One*alpha*alpha;
            bb = -1.0*Gamma_Minus_One*alpha*Riemann;
            cc =  0.5*Gamma_Minus_One*Riemann*Riemann
            -2.0*SoundSpeed_Total2/Gamma_Minus_One;

            /*--- Solve quadratic equation for velocity magnitude. Value must
             be positive, so the choice of root is clear. ---*/
            dd = bb*bb - 4.0*aa*cc;
            dd = sqrt(max(0.0,dd));
            Vel_Mag   = (-bb + dd)/(2.0*aa);
            Vel_Mag   = max(0.0,Vel_Mag);
            Velocity2 = Vel_Mag*Vel_Mag;

            /*--- Compute speed of sound from total speed of sound eqn. ---*/
            SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

            /*--- Mach squared (cut between 0-1), use to adapt velocity ---*/
            Mach2 = Velocity2/SoundSpeed2;
            Mach2 = min(1.0,Mach2);
            Velocity2   = Mach2*SoundSpeed2;
            Vel_Mag     = sqrt(Velocity2);
            SoundSpeed2 = SoundSpeed_Total2 - 0.5*Gamma_Minus_One*Velocity2;

            /*--- Compute new velocity vector at the inlet ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              Velocity[iDim] = Vel_Mag*Flow_Dir[iDim];

            /*--- Static temperature from the speed of sound relation ---*/
            Temperature = SoundSpeed2/(Gamma*Gas_Constant);

            /*--- Static pressure using isentropic relation at a point ---*/
            Pressure = P_Total*pow((Temperature/T_Total),Gamma/Gamma_Minus_One);

            /*--- Density at the inlet from the gas law ---*/
            Density = Pressure/(Gas_Constant*Temperature);

            /*--- Using pressure, density, & velocity, compute the energy ---*/
            Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Velocity2;

            /*--- Conservative variables, using the derived quantities ---*/
            FlowSolution_j[0] = Density;
            FlowSolution_j[1] = Velocity[0]*Density;
            FlowSolution_j[2] = Velocity[1]*Density;
            FlowSolution_j[3] = Energy*Density;
            if (nDim == 3) {
              FlowSolution_j[3] = Velocity[2]*Density;
              FlowSolution_j[4] = Energy*Density;
            }

            break;

            /*--- Mass flow has been specified at the inlet. ---*/
          case MASS_FLOW:

            /*--- Retrieve the specified mass flow for the inlet. ---*/
            Density  = config->GetInlet_Ttotal(Marker_Tag);
            Vel_Mag  = config->GetInlet_Ptotal(Marker_Tag);
            Flow_Dir = config->GetInlet_FlowDir(Marker_Tag);

            /*--- Non-dim. the inputs if necessary. ---*/
            Density /= config->GetDensity_Ref();
            Vel_Mag /= config->GetVelocity_Ref();

            /*--- Get primitives from current inlet state. ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              Velocity[iDim] = node[iPoint]->GetVelocity(iDim, incompressible);
            Pressure    = node[iPoint]->GetPressure(incompressible);
            SoundSpeed2 = Gamma*Pressure/FlowSolution_i[0];

            /*--- Compute the acoustic Riemann invariant that is extrapolated
             from the domain interior. ---*/
            Riemann = Two_Gamma_M1*sqrt(SoundSpeed2);
            for (iDim = 0; iDim < nDim; iDim++)
              Riemann += Velocity[iDim]*UnitaryNormal[iDim];

            /*--- Speed of sound squared for fictitious inlet state ---*/
            SoundSpeed2 = Riemann;
            for (iDim = 0; iDim < nDim; iDim++)
              SoundSpeed2 -= Vel_Mag*Flow_Dir[iDim]*UnitaryNormal[iDim];

            SoundSpeed2 = max(0.0,0.5*Gamma_Minus_One*SoundSpeed2);
            SoundSpeed2 = SoundSpeed2*SoundSpeed2;

            /*--- Pressure for the fictitious inlet state ---*/
            Pressure = SoundSpeed2*Density/Gamma;

            /*--- Energy for the fictitious inlet state ---*/
            Energy = Pressure/(Density*Gamma_Minus_One)+0.5*Vel_Mag*Vel_Mag;

            /*--- Conservative variables, using the derived quantities ---*/
            FlowSolution_j[0] = Density;
            FlowSolution_j[1] = Vel_Mag*Flow_Dir[0]*Density;
            FlowSolution_j[2] = Vel_Mag*Flow_Dir[1]*Density;
            FlowSolution_j[3] = Energy*Density;
            if (nDim == 3) {
              FlowSolution_j[3] = Vel_Mag*Flow_Dir[2]*Density;
              FlowSolution_j[4] = Energy*Density;
            }

            break;
        }
      }

      /*--- Set the conservative variable states. Note that we really only
       need the density and momentum components so that we can compute
       the velocity for the convective part of the upwind residual. ---*/
      conv_numerics->SetConservative(FlowSolution_i, FlowSolution_i);

      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/
      Solution_i[0] = node[iPoint]->GetSolution(0);
      Solution_j[0] = Intermittency_Inf*Density_Inf;

      Solution_i[1] = node[iPoint]->GetSolution(1);
      Solution_j[1] = REth_Inf*Density_Inf;

      conv_numerics->SetTransVar(Solution_i, Solution_j);
      //cout << "Solution_i, Solution_j: " << Solution_i[1] << " " << Solution_j[1] << endl;

      /*--- Set various other quantities in the conv_numerics class ---*/
      conv_numerics->SetNormal(Normal);

      if (incompressible)
        conv_numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                                   solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution ---*/
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//      
//      /*--- Conservative variables w/o reconstruction ---*/
//      visc_numerics->SetConservative(FlowSolution_i, FlowSolution_j);
//      
//      /*--- Laminar Viscosity, and density (incompresible solver)  ---*/
//      if (incompressible) {
//        visc_numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
//                                         solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc());
//        visc_numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
//                                   solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
//      } else {
//        visc_numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
//                                         solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity());
//      }
//      
//      /*--- Eddy Viscosity ---*/
//      visc_numerics->SetEddyViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
//                                    solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity());
//      
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//      visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//      
//      /*--- Compute residual, and Jacobians ---*/
//      visc_numerics->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
//      
//      /*--- Subtract residual, and update Jacobians ---*/
//      SubtractResidual(iPoint, Residual);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  //cin.get();
}

void CTransLMSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                 CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex, Point_Normal;
  unsigned short iVar, iDim;
  
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool grid_movement  = config->GetGrid_Movement();
  
  double *Normal = new double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      
      /*--- Set the conservative variables, density & Velocity same as in
       the interior, don't need to modify pressure for the convection. ---*/
      conv_numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                                   solver_container[FLOW_SOL]->node[iPoint]->GetSolution());
      
      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> TurbVar_internal,
       Solution_j --> TurbVar_outlet ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
      }
      conv_numerics->SetTransVar(Solution_i, Solution_j);
      
      /*--- Set Normal (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Set various quantities in the solver class ---*/
      if (incompressible)
        conv_numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
                                   solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
//      /*--- Viscous contribution ---*/
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//      
//      /*--- Conservative variables w/o reconstruction ---*/
//      visc_numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(),
//                                   solver_container[FLOW_SOL]->node[iPoint]->GetSolution());
//      
//      /*--- Laminar Viscosity, and density (incompresible solver)  ---*/
//      if (incompressible) {
//        visc_numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc(),
//                                         solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosityInc());
//        visc_numerics->SetDensityInc(solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc(),
//                                   solver_container[FLOW_SOL]->node[iPoint]->GetDensityInc());
//      } else {
//        visc_numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
//                                         solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity());
//      }
//      
//      /*--- Eddy Viscosity ---*/
//      visc_numerics->SetEddyViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity(),
//                                    solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity());
//      
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//      visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//      
//      /*--- Compute residual, and Jacobians ---*/
//      visc_numerics->SetResidual(Residual, Jacobian_i, Jacobian_j, config);
//      
//      /*--- Subtract residual, and update Jacobians ---*/
//      SubtractResidual(iPoint, Residual);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

void CTransLMSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
																 CConfig *config, unsigned short val_marker) {

}
