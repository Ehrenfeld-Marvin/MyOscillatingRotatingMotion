/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::solidBodyMotionFunctions::TrimMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Oscillating rotation.

SourceFiles
    TrimMotion.C

\*---------------------------------------------------------------------------*/

#ifndef TrimMotion_H
#define TrimMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "Function1.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class TrimMotion Declaration
\*---------------------------------------------------------------------------*/

class TrimMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of gravity
        point origin_;

        //- Amplitude
        mutable vector amplitude_;

        //- Radial velocity
        autoPtr<Function1<scalar>> omega_;

	//- Initial offset of the sine function
	scalar initialOffset_;

    // Private Member Functions

        //- No copy construct
        TrimMotion(const TrimMotion&) = delete;

        //- No copy assignment
        void operator=(const TrimMotion&) = delete;



	int NuOfOsc=1;
        
	mutable scalar Read_Time_Intervall;		//Intervall when the amplitude is supposed to change
	
	mutable float amplitude_target;			//Zielwert der Force
	mutable float Amplitude_K_P;			//PID Koeffizienten Amplitude
	mutable float Amplitude_K_I;
	mutable float Amplitude_K_D;
	mutable double IST_AMP;
	mutable double IST_alt_AMP;
	mutable double e_t_SUMME_AMP;
	mutable double e_t_AMP;
	mutable double e_t_alt_AMP;
	mutable double amplitude_y;
	mutable double P_AMP;
	mutable double I_AMP;
	mutable double D_AMP;
	
	mutable float amplitude_begin;
	
	mutable double phi_target;
	mutable float Phi_K_P;			//PID Koeffizienten Phasenverschiebung
	mutable float Phi_K_I;
	mutable float Phi_K_D;
	mutable double IST_PHI;
	mutable double IST_alt_PHI;
	mutable double e_t_SUMME_PHI;
	mutable double e_t_PHI;
	mutable double e_t_alt_PHI;
	mutable double phi_y;
	mutable double P_PHI;
	mutable double I_PHI;
	mutable double D_PHI;
	
	mutable scalar angle_F_X;
	

	scalar omega_PID;
	scalar startTime = time_.startTime().value();
	mutable scalar Time_Intervall_Value;
	mutable scalar Ramp_Time;
	mutable scalar Ramp_Time_Counter;
	mutable scalar amplitude_ramp;
	mutable scalar phi_ramp;
	mutable scalar phi_0;
	mutable scalar amplitude_0;
	mutable scalar delta_amplitude;
	mutable scalar delta_phi;
	mutable bool Ramp_change=false; 
	
	mutable scalar Amplitude_Value;
	mutable scalar phi_Value;

	
	mutable double Force_Average_X;
	mutable double Force_Average_Y;
	mutable double Force_Average_Z;

	mutable std::vector<double> I_vector;	

	mutable int I_iterator=0;
	mutable int I_intervall=10;


public:

    //- Runtime type information
    TypeName("TrimMotion");


    // Constructors

        //- Construct from components
        TrimMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new TrimMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~TrimMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
        
        void Amplitude(scalar* Amplitude_Value, scalar* phi_Value) const;
        
        void writePIDdata
	(    	
		string FileName,
		double Target,			
		float K_P,			
		float K_I,
		float K_D,
		double IST,
		double e_t_SUMME,
		double e_t,
		scalar y_t,
		float P,
		float I,
		float D,
		double Force_Average_X,
		double Force_Average_Y,
		double Force_Average_Z,
		float time
	) const;
  	
  	void Back_up_Data_PID
	(
		string FileName,
		scalar time,
		double IST_alt,
		double e_t_SUMME,
		double y_t
	) const;      
        
        scalar PID_REGLER
	(
		double Target,
		double Average,
		float K_P,
		float K_I,
		float K_D,
		double* e_t_SUMME,
		double* IST_alt,
		double* e_t_alt,
		double* P_PHI,
		double* I_PHI,
		double* D_PHI,
		double* IST,
		double* e_t
	) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
