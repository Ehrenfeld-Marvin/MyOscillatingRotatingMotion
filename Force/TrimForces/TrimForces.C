/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "TrimForces.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include <fstream>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(TrimForces, 0);
    addToRunTimeSelectionTable(functionObject, TrimForces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

#include "checkIfExist.H"

void Foam::functionObjects::TrimForces::createFiles()
{
	
    if (writeToFile() && !Trim_Force_Ptr.valid())
    {
    		cout << "\n\n++++++++++CREATING TRIM_FORCES.DAT++++++++++\n\n";
	Trim_Force_Ptr = createFile("Trim_Forces", startTime);
	if(!restart_without_deleting) writeIntegratedHeaderNEW("Trim_Forces", Trim_Force_Ptr());
//	Force_Cycle = createFile(FileName);
//	writeIntegratedHeaderNEW("Force_Cycle", Force_Cycle());
    }
}




void Foam::functionObjects::TrimForces::writeIntegratedHeaderNEW
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeCommented(os, "Time");
    writeTabbed(os, "total_x");
    writeTabbed(os, "total_y");
    writeTabbed(os, "total_z");

    os  << endl;
}



/*DELETE FILE WIRD ERSTMAL AUSGELASSEN, KRÄFTE ALLE IN EINER DATEI SCHREIBEN
void Foam::functionObjects::TrimForces::DeleteFile
(

) const
{
		
	
	float run_T=time().value();
	float delta=time().deltaT().value();
//	float TimeInterval= (NuOfOsc*2*3.1416)/omega;
	float TimeInterval= 1;
	string Path=time().path();
	string PostProcessingPath="/postProcessing/TrimForces_Dir/0/Force_Cycle.dat";
	Path=Path.append(PostProcessingPath);
	
	writeTrimForces
	(
      	"Trim_Forces",
      	coordSys_.localVector(TrimForce_[0]),
      	coordSys_.localVector(TrimForce_[1]),
      	coordSys_.localVector(TrimForce_[2]),
	Force_Cycle
	);


	if(TimeToDelete==0) TimeToDelete = TimeInterval+delta;		//Werte sollen erst nach einem Cycle gelöscht werden
										//+delta wegen Zeitverzögerung zwischen schreiben und lesen

	if(run_T>=TimeToDelete)
	{

		TimeToDelete += (TimeInterval+delta);	
//		FileName += std::to_string(IT);
//		IT++;
		remove(Path.c_str());
		Force_Cycle= createFile(FileName);
		writeIntegratedHeaderNEW("Trim_Forces", Force_Cycle());    		

	}


}
*/

void Foam::functionObjects::TrimForces::writeIntegratedTrimForces
(
    const string& descriptor,
    const vectorField& fm0,
    const vectorField& fm1,
    const vectorField& fm2,
    autoPtr<OFstream>& osPtr
) const
{
	vector pressure = sum(fm0);
	vector viscous = sum(fm1);
	vector porous = sum(fm2);
	vector total = pressure + viscous + porous;
   	
   	float run_T=time().value();
   	
   	if (writeToFile())
	{
		Ostream& os = osPtr();
        	
//              writeCurrentTime(os);	

		os 	<< /*setprecision(5) <<*/ run_T;
						
		if(total[0]>0)	os << tab << setprecision(15) << "+" << total[0];
		else os << setprecision(15) << tab << total[0];
		if(total[1]>0)	os << tab << "+" << total[1];
		else os << tab << total[1];
		if(total[2]>0)	os << tab << "+" << total[2];
		else os << tab << total[2];
	
		os << endl;
    	}
}



void Foam::functionObjects::TrimForces::writeTrimForces()
{
    Log << type() << " " << name() << " write:" << nl;

    
    writeIntegratedTrimForces
    (
        "TrimForces",
        coordSys_.localVector(force_[0]),
        coordSys_.localVector(force_[1]),
        coordSys_.localVector(force_[2]),
        Trim_Force_Ptr
    );

    Log << endl;
}






// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::TrimForces::TrimForces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    forces(name, runTime, dict, false),
    Trim_Force_Ptr()
//    Force_Cycle(),
{
	if (readFields)
	{
  	      read(dict);
		setCoordinateSystem(dict);
		Log << endl;
	}

	if(startTime>0)
	{
			cout << "\n\n++++++++++WRTITIMG TEMP TRIM_FORCES.DAT++++++++++\n\n";
		restart_without_deleting=true;
		
  	  	string ReadTrimForce;
  	  	string Time_str;
  	  	float Time_float;
  	  	std::string TrimForceFile_line;
		std::fstream TrimForceFile; 
		TrimForceFile.open("Trim_Forces.dat");
		std::ofstream temp;
		temp.open("temp.dat");
    	
		while(getline(TrimForceFile, TrimForceFile_line))						//ließt die Force_Cycles Datei linie für linie
		{
			ReadTrimForce = TrimForceFile_line;
			
			if(ReadTrimForce[0]!='#')							//Überspringt Header
			{
      				
				std::size_t first_tab = ReadTrimForce.find('\t');			//Find place of first tab
				Time_str = ReadTrimForce.substr(0,first_tab);					//Save first Value in str1	: equals Time
				Time_float = std::stod(Time_str);					//Converting STRING to FLOAT

				if(Time_float<=startTime) temp << TrimForceFile_line << "\n";
			}
			else temp << TrimForceFile_line << "\n";
			
		}
		
		temp.close();
		TrimForceFile.close();
		remove("Trim_Forces.dat");
		rename("temp.dat","Trim_Forces.dat");
      		
	}
    
//	string Path=time().path();
//	string PostProcessingPath="/postProcessing/TrimForces_Dir/0/Trim_Forces.dat";
//	Path=Path.append(PostProcessingPath);
	
	if(startTime==0)
	{
		cout << "\n\n++++++++++REMOVING TRIM_FORCES.DAT++++++++++\n\n";
		remove("Trim_Forces.dat");
	}
	
	
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::TrimForces::read(const dictionary& dict)
{
    forces::read(dict);

        
    NuOfOsc = dict.get<int>("Oscillations");
    omega = dict.get<scalar>("omega");
    
    writeFields_ = dict.getOrDefault("writeFields", false);


    return true;
}


bool Foam::functionObjects::TrimForces::execute()
{
    forces::calcForcesMoment();

    if (Pstream::master())
    {
        createFiles();

        writeTrimForces();

//	DeleteFile();

        Log << endl;
    }

    return true;
}



// ************************************************************************* //
