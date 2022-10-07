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
 * ReportAll2nc.cpp
 *
 *  Created on: Jan 27, 2022
 *      Author: Xiaoqiang Yang
 *      For coupling external modeling, this script stores necessary modeling state variables and fluxes as netcdf4 format file
 */

#include <netcdf>
#include <iostream>
#include <string>
#include "Sativa.h"
using namespace std;


int ReportAll2nc(){
  
  //create the nc file for water fluxes
  if(oControl->current_ts_count == 1){
    oReport->CreatALLOutputNC(oControl->path_ResultsFolder, "W");
  }  
//precipitation
	oReport->UpdateOutputNC(oAtmosphere->getPrecipitation_all(), "Pp", "W");
    //air temperature
	oReport->UpdateOutputNC(oAtmosphere->getTemperature(), "Tp", "W");
//above ground storage
    //interception storage
	oReport->UpdateOutputNC(oBasin->getCanopyStorage(), "Cs", "W");
  	//snow pack
	oReport->UpdateOutputNC(oBasin->getSnowWaterEquiv(), "SWE", "W");
	//ponding water storage
	oReport->UpdateOutputNC(oBasin->getPondingWater(), "Ponding_", "W");
//vertical fluxes
    //canopy evapration
	oReport->UpdateOutputNC(oBasin->getEvaporationI_all(), "EvapI", "W");
	//throughfall -rain
	oReport->UpdateOutputNC(oBasin->getFluxCnptoSrf(), "ThrghRn", "W"); //new
    //throughfall -snow
	oReport->UpdateOutputNC(oBasin->getFluxCnptoSnow(), "ThrghSn", "W"); //new  
	// snowmelt
	oReport->UpdateOutputNC(oBasin->getFluxSnowtoSrf(), "SnMelt", "W"); //new
	//infiltration
	oReport->UpdateOutputNC(oBasin->getFluxInfilt(), "Inf", "W");
	//L1toL2
	oReport->UpdateOutputNC(oBasin->getFluxL1toL2(), "L1toL2", "W"); //new
	//L2toL3
	oReport->UpdateOutputNC(oBasin->getFluxL2toL3(), "L2toL3", "W"); //new
	//reinfiltration	
    if (oControl->sw_reinfilt){
	  oReport->UpdateOutputNC(oBasin->getFluxSrftoL1R(), "InfR", "W");
	//re-L1toL2
	  oReport->UpdateOutputNC(oBasin->getFluxL1toL2R(), "L1toL2R", "W"); //new
	//re-L2toL3
	  oReport->UpdateOutputNC(oBasin->getFluxL2toL3R(), "L2toL3R", "W"); //new
    }
    //return flow
	//exfiltration
	oReport->UpdateOutputNC(oBasin->getFluxExfilt(), "RSrf", "W");
	//L2toL1
	oReport->UpdateOutputNC(oBasin->getFluxL2toL1(), "L2toL1", "W"); //new
	//L3toL2
	oReport->UpdateOutputNC(oBasin->getFluxL3toL2(), "L3toL2", "W"); //new
    //soil evaporation
	oReport->UpdateOutputNC(oBasin->getEvaporationS_all(), "EvapS", "W");
    //transpiration
	oReport->UpdateOutputNC(oBasin->getTranspiration_L1(), "EvapT1", "W");
	oReport->UpdateOutputNC(oBasin->getTranspiration_L2(), "EvapT2", "W");
	oReport->UpdateOutputNC(oBasin->getTranspiration_L3(), "EvapT3", "W");
//Soil water depth [m]
	oReport->UpdateOutputNC(oBasin->getSoilWaterDepthL1(), "SWD1_", "W");
	oReport->UpdateOutputNC(oBasin->getSoilWaterDepthL2(), "SWD2_", "W");
	oReport->UpdateOutputNC(oBasin->getSoilWaterDepthL3(), "SWD3_", "W");
	oReport->UpdateOutputNC(oBasin->getPorosityL1(), "satSM1", "W");
	oReport->UpdateOutputNC(oBasin->getPorosityL2(), "satSM2", "W");
	oReport->UpdateOutputNC(oBasin->getPorosityL3(), "satSM3", "W");
	oReport->UpdateOutputNC(oBasin->getWiltingPoint(), "Wiltpnt", "W"); //NEW
  //soil temperature	
	oReport->UpdateOutputNC(oBasin->getSoilTemp(), "Ts", "W");
	//conceptual grounwater storage for subsurface flow generation
	oReport->UpdateOutputNC(oBasin->getGrndWater(), "GW", "W");
//lateral exchange
    //surface outflow
	oReport->UpdateOutputNC(oBasin->getFluxSrftoLat(), "LSrfo", "W");
	//surface inflow
	oReport->UpdateOutputNC(oBasin->getFluxLattoSrf(), "LSrfi", "W");
    //subsurface outflow
	oReport->UpdateOutputNC(oBasin->getFluxGWtoLat(), "LGWo", "W");
	//subsurface inflow
	oReport->UpdateOutputNC(oBasin->getFluxLattoGW(), "LGWi", "W");
    //surface flow to channel
	oReport->UpdateOutputNC(oBasin->getFluxSrftoChn(), "SrfChn", "W");
    //surface flow tracer mixing ratio
    if (oControl->toggle_SrfOVFmix == 1){
	  oReport->UpdateOutputNC(oBasin->getRatioSrfOVF(), "RSrfmix", "W"); //added 2022-08
    }
	//subsurface flow to channel
	oReport->UpdateOutputNC(oBasin->getFluxGWtoChn(), "GWChn", "W");
    //channel storage [m3]
	oReport->UpdateOutputNC(oBasin->getChanstoreV(), "ChnS", "W"); //m^3
    //channel upstream input	
	oReport->UpdateOutputNC(oBasin->getChanUpin(), "LChQUi", "W"); //new FluxLattoChn --> m^3/s
	//channel output to down stream
	oReport->UpdateOutputNC(oBasin->getChanDownout(), "LChQDo", "W"); //new FluxChntoLat --> m^3/s
    //i.e., the final discharge
	oReport->UpdateOutputNC(oBasin->getStreamflow(), "Q", "W");
//percolation and deeper GW dynamics
	oReport->UpdateOutputNC(oBasin->getBedrockLeakage_Depth(), "LeakD", "W");	    
    //extra GW yangx 2020-05
    if(oControl->sw_extraGW){
	  oReport->UpdateOutputNC(oBasin->getExtraGW(), "ExtraGW", "W");
	  oReport->UpdateOutputNC(oBasin->getFluxExtraGWtoChn(), "ExGWtC", "W");
	  oReport->UpdateOutputNC(oBasin->getFluxLattoExtraGW(), "LtExGW", "W");
	  oReport->UpdateOutputNC(oBasin->getFluxExtraGWtoLat(), "ExGWtL", "W");
	  //oReport->UpdateOutputNC(oBasin->getAccExtraGWtoChn(), "ExGWtCA", "W");
	  //oReport->UpdateOutputNC(oBasin->getAccLattoExtraGW(), "LtExGWA", "W");
	  //oReport->UpdateOutputNC(oBasin->getAccExtraGWtoLat(), "ExGWtLA", "W");
    }

	//oReport->UpdateOutputNC(oBasin->getSoilTemp(), "Ts", "W");
    //oReport->UpdateOutputNC(oBasin->getTemp_w(), "Twater", "W");
    //oReport->UpdateOutputNC(oBasin->getChanEvap(), "EvapC", "W");
	//oReport->UpdateOutputNC(oBasin->getEvaporation(), "Evap", "W");

  return EXIT_SUCCESS;
}




