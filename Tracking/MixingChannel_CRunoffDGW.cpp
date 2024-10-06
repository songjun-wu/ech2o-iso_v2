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
 * MixingChannel_CRunoffDGW.cpp
 *
 *  (1) A separate subroutine for updating tracer concentrations of all runoff components that
 *  seepage into channels
 
 *  Created on: May 5, 2022
 *      Author: Xiaoqiang Yang
 */

#include "Basin.h"

void Tracking::MixingChannel_CRunoffDGW(Basin &bsn, Control &ctrl,
			    double &poros1, int r, int c) // 
{
  int SrfOVFmix = ctrl.toggle_SrfOVFmix;
  int mixmod = ctrl.toggle_mix;
  double ExtraGWtoChn = bsn.getFluxExtraGWtoChn()->matrix[r][c];
  double GWtoChn = bsn.getFluxGWtoChn()->matrix[r][c];
  double SrftoChn = bsn.getFluxSrftoChn()->matrix[r][c];
  double LattoChn = bsn.getFluxLattoChn()->matrix[r][c];
  double ChntoLat = bsn.getFluxChntoLat()->matrix[r][c];// Qk1*dtdx/dx;
  //separate channel storage and ponding waters
  double FinSrf = SrftoChn + GWtoChn+ LattoChn + ExtraGWtoChn; //sum of channel inputs
  //double FinSrf2 = L1toSrf ; //sum of ponding inputs that needed to be mixed, the rest has already handled in mixingV_down.cpp
  double channel_old = bsn.getChannel_old()->matrix[r][c];//getChannelWater()

  double d2HinChn = 0;
  double d18OinChn = 0;
  double AgeinChn= 0;


  double theta1 = bsn.getSoilMoist1()->matrix[r][c];
  double theta_MW1 = 0;
  double ratio_theta;

  /*if ( SrfOVFmix == 0){
    //ratio_theta = std::max<double>(0.5,theta1/poros1); //the proportion of dsurface and dsoil1
    ratio_theta = 0;
  }else{
	ratio_theta = bsn.getRatioSrfOVF()->matrix[r][c];  
  }


  if(ctrl.sw_TPD){
    theta_MW1 = bsn.getMoistureMW1()->matrix[r][c];
	if ( SrfOVFmix == 0){
      //ratio_theta = std::max<double>(0.5,theta_MW1/poros1); //the proportion of dsurface and dsoil1
      ratio_theta = 0;
    }else{
	  ratio_theta = bsn.getRatioSrfOVF()->matrix[r][c];  
    }
  }
  //surface overland flow to channel
  //if (ctrl.sw_partial_reinf){
    //in this case, surface OVF conc gets variables ended with _r
    if(ctrl.sw_TPD ){
      if(ctrl.sw_2H){    
	    _d2HSrftoChn->matrix[r][c] = ratio_theta *_d2H_MW1->matrix[r][c] + 
		                             (1-ratio_theta) *_d2Hsurface->matrix[r][c];
      }
      if(ctrl.sw_18O){
	    _d18OSrftoChn->matrix[r][c] = ratio_theta * _d18O_MW1->matrix[r][c] +
	                                 (1-ratio_theta) * _d18Osurface->matrix[r][c];  
      }  
      if(ctrl.sw_Age){
	    _AgeSrftoChn->matrix[r][c] = ratio_theta * _Age_MW1->matrix[r][c] +
	                                 (1-ratio_theta) * _Agesurface->matrix[r][c]; 
      }
    }else { // Soil-averaged
      if(ctrl.sw_2H)
	    _d2HSrftoChn->matrix[r][c] = ratio_theta * _d2Hsoil1->matrix[r][c] +
	                                 (1-ratio_theta) * _d2Hsurface->matrix[r][c];
    
      if(ctrl.sw_18O)
	    _d18OSrftoChn->matrix[r][c] = ratio_theta * _d18Osoil1->matrix[r][c] +
	                                 (1-ratio_theta) * _d18Osurface->matrix[r][c];
    
      if(ctrl.sw_Age)
        _AgeSrftoChn->matrix[r][c] =  ratio_theta * _Agesoil1->matrix[r][c] +
	                                 (1-ratio_theta) * _Agesurface->matrix[r][c];
    }
    */

   /* else {
      if(ctrl.sw_2H)
	    _d2HSrftoChn->matrix[r][c] = _d2Hsurface->matrix[r][c];
    
      if(ctrl.sw_18O)
	    _d18OSrftoChn->matrix[r][c] = _d18Osurface->matrix[r][c];
    
      if(ctrl.sw_Age)
        _AgeSrftoChn->matrix[r][c] =  _Agesurface->matrix[r][c];
  } */

  // Channel mixing--------------------------------------------------------------------------------
 
  if(FinSrf > RNDOFFERR) {
    if(ctrl.sw_2H){
	//for channel water, to update inputs: GWtoChn, SrftoChn,  LattoChn and ExtraGWtoChn
	    d2HinChn = (GWtoChn *_d2Hgroundwater->matrix[r][c] + SrftoChn*_d2HSrftoChn->matrix[r][c] +
		    ExtraGWtoChn *_d2HExtraGWtoLat->matrix[r][c] + _Fd2HLattoChn->matrix[r][c]) / FinSrf ;
	    _d2Hchan->matrix[r][c] = channel_old > RNDOFFERR ?
	         InOutMix(channel_old, _d2Hchan->matrix[r][c], FinSrf, d2HinChn, 
		     ChntoLat, mixmod): d2HinChn;          
    }

    if(ctrl.sw_18O){
	    d18OinChn = (GWtoChn *_d18Ogroundwater->matrix[r][c] + SrftoChn*_d18OSrftoChn->matrix[r][c] +
		    ExtraGWtoChn *_d18OExtraGWtoLat->matrix[r][c] + _Fd18OLattoChn->matrix[r][c]) / FinSrf ;
	    _d18Ochan->matrix[r][c] = channel_old > RNDOFFERR ?
	         InOutMix(channel_old, _d18Ochan->matrix[r][c], FinSrf, d18OinChn, 
		     ChntoLat, mixmod): d18OinChn;          
    }

    if(ctrl.sw_Age){
	    AgeinChn = (GWtoChn *_Agegroundwater->matrix[r][c] + SrftoChn*_AgeSrftoChn->matrix[r][c] +
		    ExtraGWtoChn *_AgeExtraGWtoLat->matrix[r][c] + _FAgeLattoChn->matrix[r][c]) / FinSrf ;
	    _Agechan->matrix[r][c] = channel_old > RNDOFFERR ?
	         InOutMix(channel_old, _Agechan->matrix[r][c], FinSrf, AgeinChn, 
		     ChntoLat, mixmod): AgeinChn;          
        _AgeGWtoChn->matrix[r][c] = _Agegroundwater->matrix[r][c];
    }
  }      
}

