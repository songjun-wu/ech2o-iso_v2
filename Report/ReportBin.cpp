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
 *    Marco Maneta
 *******************************************************************************/
/*
 * ReportBin.cpp
 *  Created on: Oct 12, 2022
 */

#include <fstream>
#include "Report.h"
#include "Basin.h"
#include "Sativa.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int Report::CreatALLOutputBin(string filepath){
   REAL8 Rsol;
   UINT4  NLAT, NLON;
   REAL8 START_LAT,START_LON;
   REAL8 NODATA;
   std::vector<std::string> RprtNamesW;
   std::vector<std::string> outvars;
   string filename;
   
  //checking output variables
        RprtNamesW.push_back("Pp");
		RprtNamesW.push_back("Tp");
		RprtNamesW.push_back("Cs");
		RprtNamesW.push_back("SWE");
		RprtNamesW.push_back("Ponding_");
		RprtNamesW.push_back("EvapI");		
		RprtNamesW.push_back("ThrghRn");//NEW	
		RprtNamesW.push_back("ThrghSn");//NEW		
		RprtNamesW.push_back("SnMelt");//NEW
		RprtNamesW.push_back("Inf");
		RprtNamesW.push_back("L1toL2");//NEW
		RprtNamesW.push_back("L2toL3");//NEW
	 //reinfiltration	
        if (oControl->sw_reinfilt){
		  RprtNamesW.push_back("InfR");//NEW
		  RprtNamesW.push_back("L1toL2R");//NEW
		  RprtNamesW.push_back("L2toL3R");//NEW		
		}
		RprtNamesW.push_back("RSrf");
		RprtNamesW.push_back("L2toL1");
		RprtNamesW.push_back("L3toL2");
		RprtNamesW.push_back("EvapS");
		RprtNamesW.push_back("EvapT1");//NEW	
		RprtNamesW.push_back("EvapT2");//NEW	
		RprtNamesW.push_back("EvapT3");//NEW	
		RprtNamesW.push_back("SWD1_");//NEW
		RprtNamesW.push_back("SWD2_");//NEW
		RprtNamesW.push_back("SWD3_");//NEW
		RprtNamesW.push_back("satSM1");//NEW
		RprtNamesW.push_back("satSM2");//NEW
		RprtNamesW.push_back("satSM3");//NEW
		RprtNamesW.push_back("Wiltpnt");//NEW
		RprtNamesW.push_back("Ts");
		RprtNamesW.push_back("GW");
		RprtNamesW.push_back("LSrfo");
		RprtNamesW.push_back("LSrfi");
		RprtNamesW.push_back("LGWo");
		RprtNamesW.push_back("LGWi");
		RprtNamesW.push_back("SrfChn");
		if (oControl->toggle_SrfOVFmix == 1){
		  RprtNamesW.push_back("RSrfmix");  //NEW
		}
		RprtNamesW.push_back("GWChn");
 		RprtNamesW.push_back("ChnS"); //NEW       
		RprtNamesW.push_back("LChQUi");//NEW
		RprtNamesW.push_back("LChQDo");//NEW
		RprtNamesW.push_back("Q");
		RprtNamesW.push_back("LeakD");	
	    //extra GW yangx 2020-05
	    if(oControl->sw_extraGW){
		  RprtNamesW.push_back("ExtraGW");
		  RprtNamesW.push_back("ExGWtC");
		  RprtNamesW.push_back("LtExGW");
		  RprtNamesW.push_back("ExGWtL");
	    }

   
	  outvars = RprtNamesW;
      try{
		for (auto item : outvars) { 
          filename =  filepath + item.c_str() + ".bin";
          ofOutput.open(filename.c_str(), ios::binary);
          ofOutput.close();
		}
        return 0;
        } catch(const exception& e){
			    cerr << "Couldn't create binary file with message " << e.what() << endl;
			    exit(EXIT_FAILURE);
			  }
    }


int Report::UpdateOutputBin(const grid *input, string varname, string filepath){
   UINT4 NLAT, NLON;
   int TScount;
   string filename;
   //from Basin
   NLAT = oBasin->getNumRows();
   NLON = oBasin->getNumCols();  
   TScount = oControl->current_ts_count;
   REAL8 outdata[NLAT*NLON];   
   //get the outdata from input map
   for(unsigned i = 0; i < NLAT; i++){
     for(unsigned j = 0; j < NLON; j++){
       outdata[i*NLON + j] = input->matrix[i][j];
	 }
   }   
    filename =  filepath + varname.c_str() + ".bin";
    ofOutput.open(filename.c_str(), ios::binary|ios::app);

    ofOutput.write((char*)&outdata, sizeof(double)*NLAT*NLON);
 
    //ofOutput.write((char*)(&outdata), (sizeof(double)*NLAT*NLON));
    ofOutput.close();
    return 0;

}
