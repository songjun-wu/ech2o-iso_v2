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
 *  Created on: Oct 12, 2022
 *      Author: Songjun Wu
 *      For coupling external modeling, this script stores necessary modeling state variables and fluxes as binary format file
 */


#include <iostream>
#include <string>
#include "Sativa.h"
using namespace std;


int ReportAll2Bin(){
  
  //create the binary file for water fluxes
  if(oControl->current_ts_count == 1){
    oReport->CreatALLOutputBin(oControl->path_ResultstmpFolder);      
  }
  oReport->UpdateOutputBin(oAtmosphere->getPrecipitation_all(), "Pp", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oAtmosphere->getTemperature(), "Tp", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getCanopyStorage(), "Cs", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getSnowWaterEquiv(), "SWE", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getPondingWater_before_infiltration(), "Ponding_", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getEvaporationI_all(), "EvapI", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxCnptoSrf(), "ThrghRn", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxCnptoSnow(), "ThrghSn", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxSnowtoSrf(), "SnMelt", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxInfilt(), "Inf", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxL1toL2(), "L1toL2", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxL2toL3(), "L2toL3", oControl->path_ResultstmpFolder);

  if (oControl->sw_reinfilt){
  oReport->UpdateOutputBin(oBasin->getFluxSrftoL1R(), "InfR", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxL1toL2R(), "L1toL2R", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxL2toL3R(), "L2toL3R", oControl->path_ResultstmpFolder);
  }

  oReport->UpdateOutputBin(oBasin->getFluxExfilt(), "RSrf", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxL2toL1(), "L2toL1", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxL3toL2(), "L3toL2", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getEvaporationS_all(), "EvapS", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getTranspiration_L1(), "EvapT1", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getTranspiration_L2(), "EvapT2", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getTranspiration_L3(), "EvapT3", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getSoilWaterDepthL1(), "SWD1_", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getSoilWaterDepthL2(), "SWD2_", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getSoilWaterDepthL3(), "SWD3_", oControl->path_ResultstmpFolder);
  //oReport->UpdateOutputBin(oBasin->getPorosityL1(), "satSM1", oControl->path_ResultstmpFolder);
  //oReport->UpdateOutputBin(oBasin->getPorosityL2(), "satSM2", oControl->path_ResultstmpFolder);
  //oReport->UpdateOutputBin(oBasin->getPorosityL3(), "satSM3", oControl->path_ResultstmpFolder);
  //oReport->UpdateOutputBin(oBasin->getWiltingPoint(), "Wiltpnt", oControl->path_ResultstmpFolder);
  //oReport->UpdateOutputBin(oBasin->getSoilTemp(), "Ts", oControl->path_ResultstmpFolder);
  //oReport->UpdateOutputBin(oBasin->getGrndWater(), "GW", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxSrftoLat(), "LSrfo", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxLattoSrf(), "LSrfi", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxGWtoLat(), "LGWo", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxLattoGW(), "LGWi", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxSrftoChn(), "SrfChn", oControl->path_ResultstmpFolder);

  if (oControl->toggle_SrfOVFmix == 1){
    oReport->UpdateOutputBin(oBasin->getRatioSrfOVF(), "RSrfmix", oControl->path_ResultstmpFolder);
  }

  oReport->UpdateOutputBin(oBasin->getFluxGWtoChn(), "GWChn", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getChanstoreV(), "ChnS", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getChanUpin(), "LChQUi", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getChanDownout(), "LChQDo", oControl->path_ResultstmpFolder);
  //oReport->UpdateOutputBin(oBasin->getStreamflow(), "Q", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getBedrockLeakage_Depth(), "LeakD", oControl->path_ResultstmpFolder);

  if(oControl->sw_extraGW){
  oReport->UpdateOutputBin(oBasin->getExtraGW(), "ExtraGW", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxExtraGWtoChn(), "ExGWtC", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxLattoExtraGW(), "LtExGW", oControl->path_ResultstmpFolder);
  oReport->UpdateOutputBin(oBasin->getFluxExtraGWtoLat(), "ExGWtL", oControl->path_ResultstmpFolder);
  }

  return 0;

}