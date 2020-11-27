/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "MyoscillatingRotatingMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include <fstream>
#include <iostream>



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(MyoscillatingRotatingMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        MyoscillatingRotatingMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::MyoscillatingRotatingMotion::
MyoscillatingRotatingMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
    if(omega_==0)omega_=1;
    float Time_Intervall_Value = (NuOfOsc*2*3.141592653)/omega_;
    if(startTime>0)
    {

    	while(Read_Time_Intervall<startTime) Read_Time_Intervall+=Time_Intervall_Value;

/*+++++++++++ READING PID BACK UP DATA ++++++*/
    	string ReadFile;
    	string FILE;
    	string Path=time_.path();
	string FileName="/Back_Up_Data_PID.dat";
	double Value;
	Path=Path.append(FileName);
    	std::ifstream Back_Up_Data_PID;
    	Back_Up_Data_PID.open(Path);
    	int i=0;
    	while(getline(Back_Up_Data_PID, FILE))
    	{
    		
    		ReadFile = FILE;
    		Value = std::stod(ReadFile);
    		if(i==0) IST_alt = Value; 
		if(i==1) e_t_SUMME = Value;
		if(i==2) e_t_alt = Value;
		i++;
    	}
    	
/*++++++++++++++++++++++++++++++++++++++++++++*/
    }	
    	if(startTime==0)
	{
		IST_alt=0;
		e_t_SUMME=0;
		Force_Average_X=0;
		Force_Average_Y=0;
		Force_Average_Z=0;
		Read_Time_Intervall=Time_Intervall_Value;
		
	}
	amplitude_.x()=0;
	amplitude_.y()=0;
	
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

#include "WriteFilePID.H"
#include "Amplitude_Member_Function.H"
#include "PID_BackUp_Data.H"


Foam::septernion
Foam::solidBodyMotionFunctions::MyoscillatingRotatingMotion::
transformation() const
{

    scalar t = time_.value();
    
	
//	vector eulerAngles = Amplitude()*sin(omega_*t);			
	vector eulerAngles = Amplitude()*(-1);



    // Convert the rotational motion from deg to rad
//    eulerAngles *= degToRad();

    quaternion R(quaternion::XYZ, eulerAngles);
    septernion TR(septernion(-origin_)*R*septernion(origin_));

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;
    
    


    return TR;
}





bool Foam::solidBodyMotionFunctions::MyoscillatingRotatingMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.readEntry("origin", origin_);
//    SBMFCoeffs_.readEntry("amplitude", amplitude_);
    SBMFCoeffs_.readEntry("omega", omega_);
    SBMFCoeffs_.readEntry("Oscillations", NuOfOsc);
    SBMFCoeffs_.readEntry("Target", Target);
    SBMFCoeffs_.readEntry("K_P", K_P);
    SBMFCoeffs_.readEntry("K_I", K_I);
    SBMFCoeffs_.readEntry("K_D", K_D);
    

    return true;
}





// ************************************************************************* //
