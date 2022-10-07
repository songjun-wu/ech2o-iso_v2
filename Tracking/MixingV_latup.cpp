/*******************************************************************************
 * Ech2o, a spatially-distributed, ecohydrologic simulator
 * Copyright (c) 2016 Marco Maneta <marco.maneta@umontana.edu>
 *
 *     This file is part of ech2o, a hydrologic model developed at the 
 *     University of Montana.
 *
 *     Ech2o is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     Ech2o is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with Ech2o.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contributors:
 *    Marco Maneta, Sylvain Kuppel
 *******************************************************************************/
/*
 * MixingV_latup.cpp
 *
 *  Created on: Apr 25, 2018
 *      Author: Sylvain Kuppel
 */

#include "Basin.h"

void Tracking::MixingV_latup(Basin &bsn, Control &ctrl, 
			     double &d1, double &d2, double &d3, double &fc,
			     double &Qk1, double &dtdx, double &dx, int r, int c)
{
  int mixmod = ctrl.toggle_mix;

  // Soil state before routing
  double pond_old = bsn.getPondingWater()->matrix[r][c];  //without including exfiltration, so it's "old"
  double theta1_old = bsn.getSoilMoist1()->matrix[r][c];
  double theta2_old = bsn.getSoilMoist2()->matrix[r][c];
  double theta3_old = bsn.getSoilMoist3()->matrix[r][c];

  // Vertical fluxes
  double L1toSrf = bsn.getFluxExfilt()->matrix[r][c];
  double L2toL1 = bsn.getFluxL2toL1()->matrix[r][c];
  double L3toL2 = bsn.getFluxL3toL2()->matrix[r][c];


  // Lateral out
  double GWtoLat = bsn.getFluxGWtoLat()->matrix[r][c];
  double SrftoLat = bsn.getFluxSrftoLat()->matrix[r][c];
  // Lateral in
  double LattoGW = bsn.getFluxLattoGW()->matrix[r][c]; 
  //double LattoSrf = bsn.getFluxLattoSrf()->matrix[r][c];

  double d2Hin = 0;
  double d18Oin = 0;
  double Agein= 0;

  // Two-pore stuff
  double theta_MW1 = 0;
  double theta_MW2 = 0;
  double theta_r = 0;
  //double porosity = 0;
  double L3toTB2 = 0;
  double L3toMW2 = 0;
  double MW2toTB1 = 0;
  double MW2toMW1 = 0;
  double MW1toSrf = 0;

  double d_old = 0;

  if(ctrl.sw_TPD){
    theta_r = bsn.getSoilMoistR()->matrix[r][c];
    theta_MW1 = bsn.getMoistureMW1()->matrix[r][c];
    theta_MW2 = bsn.getMoistureMW2()->matrix[r][c];
    // Return flow to L2: weighted between TB2 (if there's deficit there) and MW2
    L3toTB2 = std::min<double>(L3toL2,//*(theta_MW-theta_r)/(porosity-theta_r),
				std::max<double>(0,d2*(theta_MW2-theta2_old)));
    L3toMW2 = std::max<double>(0, L3toL2 - L3toTB2);
    // Return flow to L1: : weighted between TB1 (if there's deficit there) and MW1
    MW2toTB1 = std::min<double>(L2toL1,//*(theta_MW-theta_r)/(porosity-theta_r),
				std::max<double>(0,d1*(theta_MW1-theta1_old)));
    MW2toMW1 = std::max<double>(0,L2toL1 - MW2toTB1);
    // Return flow to surface: : only from mobile water in L1
    MW1toSrf = L1toSrf;
  }

  // Layer 3 (GW included) --------------------------------------------------------------------

  if(LattoGW > RNDOFFERR){  
    if(ctrl.sw_2H){
      // Equivalent input signature for GW
      d2Hin = _Fd2HLattoGW->matrix[r][c] / LattoGW ;
      _d2Hsoil3->matrix[r][c] = InOutMix(theta3_old*d3, _d2Hsoil3->matrix[r][c],
					 LattoGW, d2Hin, L3toL2+GWtoLat, mixmod);
      _d2Hgroundwater->matrix[r][c] = _d2Hsoil3->matrix[r][c];
    }
    if(ctrl.sw_18O){
      d18Oin = _Fd18OLattoGW->matrix[r][c] / LattoGW ;
      _d18Osoil3->matrix[r][c] = InOutMix(theta3_old*d3, _d18Osoil3->matrix[r][c],
					  LattoGW, d18Oin, L3toL2+GWtoLat, mixmod);
      _d18Ogroundwater->matrix[r][c] = _d18Osoil3->matrix[r][c];
    }
    if(ctrl.sw_Age){
      Agein = _FAgeLattoGW->matrix[r][c] / LattoGW ;
      _Agesoil3->matrix[r][c] = InOutMix(theta3_old*d3, _Agesoil3->matrix[r][c],
					 LattoGW, Agein, L3toL2+GWtoLat, mixmod);
      _Agegroundwater->matrix[r][c] = _Agesoil3->matrix[r][c];
    }
  }

  // Layer 2 ------------------------------------------------------------------------

  // If two-pore domain activated
  if(ctrl.sw_TPD and L3toL2 > RNDOFFERR){
    // Tightly-bound
    if(ctrl.sw_2H)
      _d2H_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2, 
					_d2H_TB2->matrix[r][c],
					L3toTB2, _d2Hsoil3->matrix[r][c]);
    if(ctrl.sw_18O)
      _d18O_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2, 
					 _d18O_TB2->matrix[r][c],
					 L3toTB2, _d18Osoil3->matrix[r][c]);
    if(ctrl.sw_Age)
      _Age_TB2->matrix[r][c] = InputMix(std::min<double>(theta2_old,theta_MW2)*d2, 
					_Age_TB2->matrix[r][c],
					L3toTB2, _Agesoil3->matrix[r][c]);
    
    // Mobile water
    if(ctrl.sw_2H)
      _d2H_MW2->matrix[r][c] = std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta2_old-theta_MW2)*d2, _d2H_MW2->matrix[r][c], 
		 L3toMW2, _d2Hsoil3->matrix[r][c], L2toL1, mixmod) : _d2Hsoil3->matrix[r][c] ;

    if(ctrl.sw_18O)
      _d18O_MW2->matrix[r][c] = std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta2_old-theta_MW2)*d2, _d18O_MW2->matrix[r][c], 
		 L3toMW2, _d18Osoil3->matrix[r][c], L2toL1, mixmod) : _d18Osoil3->matrix[r][c] ;
    
    if(ctrl.sw_Age)
      _Age_MW2->matrix[r][c] = std::max<double>(0,theta2_old-theta_MW2) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta2_old-theta_MW2)*d2, _Age_MW2->matrix[r][c], 
		 L3toMW2, _Agesoil3->matrix[r][c], L2toL1, mixmod) : _Agesoil3->matrix[r][c] ;

  } else if (L3toL2 > RNDOFFERR) { // Soil-averaged values
    if(ctrl.sw_2H)
      _d2Hsoil2->matrix[r][c] = InOutMix(theta2_old*d2, _d2Hsoil2->matrix[r][c],
					 L3toL2, _d2Hsoil3->matrix[r][c], L2toL1, mixmod);
    if(ctrl.sw_18O)
      _d18Osoil2->matrix[r][c] = InOutMix(theta2_old*d2, _d18Osoil2->matrix[r][c],
					  L3toL2, _d18Osoil3->matrix[r][c], L2toL1, mixmod);
    
    if(ctrl.sw_Age)
      _Agesoil2->matrix[r][c] = InOutMix(theta2_old*d2, _Agesoil2->matrix[r][c],
					 L3toL2, _Agesoil3->matrix[r][c], L2toL1, mixmod);
  }

  // Layer 1 ------------------------------------------------------------------------

  // If two-pore domain activated: return flow only from MW2
  if(ctrl.sw_TPD and L2toL1 > RNDOFFERR){
    // Tightly-bound
    if(ctrl.sw_2H){
      d_old = _d2H_TB1->matrix[r][c];      
      _d2H_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1, 
					_d2H_TB1->matrix[r][c],	
					MW2toTB1, _d2H_MW2->matrix[r][c]);

      if(abs(_d2H_TB1->matrix[r][c])>100)
	cout << r << " " << c << "| d2H " << 
	  "| dTB1_new:" << _d2H_TB1->matrix[r][c] << "| dTB1_old:" << d_old << 
	  "| dMW2:" << _d2H_MW2->matrix[r][c] << "| MW2toTB1:" << MW2toTB1 << endl;
    }
    if(ctrl.sw_18O){
      d_old = _d18O_TB1->matrix[r][c];      
      _d18O_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1, 
					 _d18O_TB1->matrix[r][c],
					 MW2toTB1, _d18O_MW2->matrix[r][c]);
      
      if(abs(_d18O_TB1->matrix[r][c])>100)
	cout << r << " " << c << "| d18O " << 
	  "| dTB1_new:" << _d18O_TB1->matrix[r][c] << "| dTB1_old:" << d_old << 
	  "| dMW2:" << _d18O_MW2->matrix[r][c] << "| MW2toTB1:" << MW2toTB1 << endl;
    }
    if(ctrl.sw_Age)
      _Age_TB1->matrix[r][c] = InputMix(std::min<double>(theta1_old,theta_MW1)*d1, 
					_Age_TB1->matrix[r][c],
					MW2toTB1, _Age_MW2->matrix[r][c]);

    // Mobile water
    if(ctrl.sw_2H){
      d_old = _d2H_MW1->matrix[r][c];
      _d2H_MW1->matrix[r][c] = std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta1_old-theta_MW1)*d1, _d2H_MW1->matrix[r][c], 
		 MW2toMW1, _d2H_MW2->matrix[r][c], L1toSrf, mixmod) : _d2H_MW2->matrix[r][c] ;

      if(abs(_d2H_MW1->matrix[r][c])>100 ) //or (r==80 and c==119))
	cout << r << " " << c << "| d2H " << 
	  "| dMW1_new:" << _d2H_MW1->matrix[r][c] << "| dMW1_old:" << d_old << 
	  "| dMW2:" << _d2H_MW2->matrix[r][c] <<
	  "| MW2toMW1:" << MW2toMW1 << "| MW1toSrf:" << L1toSrf << endl;
    }
    if(ctrl.sw_18O){
      d_old = _d18O_MW1->matrix[r][c];
      _d18O_MW1->matrix[r][c] = std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta1_old-theta_MW1)*d1, _d18O_MW1->matrix[r][c], 
		 MW2toMW1, _d18O_MW2->matrix[r][c], L1toSrf, mixmod) : _d18O_MW2->matrix[r][c] ;

      if(abs(_d18O_MW1->matrix[r][c])>100)// or (r==80 and c==119))
	cout << r << " " << c << "| d2H " << 
	  "| dMW1_new:" << _d18O_MW1->matrix[r][c] << "| dMW1_old:" << d_old << 
	  "| dMW2:" << _d18O_MW2->matrix[r][c] <<
	  "| MW2toMW1:" << MW2toMW1 << "| MW1toSrf:" << L1toSrf << endl;
    }

    if(ctrl.sw_Age)
      _Age_MW1->matrix[r][c] = std::max<double>(0,theta1_old-theta_MW1) > RNDOFFERR ?
	InOutMix(std::max<double>(0,theta1_old-theta_MW1)*d1, _Age_MW1->matrix[r][c], 
		 MW2toMW1, _Age_MW2->matrix[r][c], L1toSrf, mixmod) : _Age_MW2->matrix[r][c] ;

  } else if (L2toL1 > RNDOFFERR) { // Soil-averaged values    
    if(ctrl.sw_2H)
      _d2Hsoil1->matrix[r][c] = InOutMix(theta1_old*d1, _d2Hsoil1->matrix[r][c],
					 L2toL1, _d2Hsoil2->matrix[r][c], L1toSrf, mixmod);
    if(ctrl.sw_18O)
      _d18Osoil1->matrix[r][c] = InOutMix(theta1_old*d1, _d18Osoil1->matrix[r][c],
					  L2toL1, _d18Osoil2->matrix[r][c], L1toSrf, mixmod);
    if(ctrl.sw_Age)
      _Agesoil1->matrix[r][c] = InOutMix(theta1_old*d1, _Agesoil1->matrix[r][c],
					 L2toL1, _Agesoil2->matrix[r][c], L1toSrf, mixmod);
  }

  // Surface --------------------------------------------------------------------------------
  //_d2Hsurface_r->matrix[r][c] = _d2Hsurface->matrix[r][c];
  //_d18Osurface_r->matrix[r][c] = _d18Osurface->matrix[r][c];
  //_Agesurface_r->matrix[r][c] = _Agesurface->matrix[r][c];
  // update tracer signals if there is return flow 
  if(L1toSrf > RNDOFFERR) {
    // If two-pore domain activated: return flow only from MW1
    if(ctrl.sw_TPD){
      if(ctrl.sw_2H){
          d2Hin = _d2H_MW1->matrix[r][c];
	      _d2Hsurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d2Hsurface->matrix[r][c], L1toSrf, d2Hin, 
		     SrftoLat, mixmod): d2Hin;
      }
      if(ctrl.sw_18O){
          d18Oin = _d18O_MW1->matrix[r][c];//d2h for L1toSrf
	      _d18Osurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d18Osurface->matrix[r][c], L1toSrf, d18Oin, 
		     SrftoLat, mixmod): d18Oin;
      }
      if(ctrl.sw_Age){
          Agein = _Age_MW1->matrix[r][c];//d2h for L1toSrf
	      _Agesurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _Agesurface->matrix[r][c], L1toSrf, Agein, 
		     SrftoLat, mixmod): Agein;
      }
    // *********************************
    } else { 		
      if(ctrl.sw_2H){
		//for ponding water, to update input exfiltration (L1toSrf), output SrftoLat  
          d2Hin =_d2Hsoil1->matrix[r][c];//d2h for L1toSrf
	      _d2Hsurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d2Hsurface->matrix[r][c], L1toSrf, d2Hin, 
		     SrftoLat, mixmod): d2Hin;		
      }     
      if(ctrl.sw_18O){
		//for ponding water, to update input exfiltration (L1toSrf), output SrftoLat 
          d18Oin = _d18Osoil1->matrix[r][c];//d2h for L1toSrf
	      _d18Osurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _d18Osurface->matrix[r][c], L1toSrf, d18Oin, 
		     SrftoLat, mixmod): d18Oin;
      }     
      if(ctrl.sw_Age) {
		//for ponding water, to update input exfiltration (L1toSrf), output SrftoLat 
          Agein = _Agesoil1->matrix[r][c];//d2h for L1toSrf
	      _Agesurface->matrix[r][c] = pond_old > RNDOFFERR ?
	         InOutMix(pond_old, _Agesurface->matrix[r][c], L1toSrf, Agein, 
		     SrftoLat, mixmod): Agein;
      }
    }
  }
  // -------------------------------------------------------------------------------
}
