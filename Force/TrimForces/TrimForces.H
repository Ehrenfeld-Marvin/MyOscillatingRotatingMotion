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

#ifndef functionObjects_TrimForces_H
#define functionObjects_TrimForces_H

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class TrimForces Declaration
\*---------------------------------------------------------------------------*/

class TrimForces
:
    public forces
{

protected:

    // Protected data


        // File streams

            //- TrimForces
	autoPtr<OFstream> Trim_Force_Ptr;
//	mutable autoPtr<OFstream> Force_Cycle;



        // Read from dictionary

            
	mutable scalar Amplitude;
	int NuOfOsc;
	scalar omega;
	scalar startTime=time().startTime().value();
	bool restart_without_deleting;

            //- Write fields flag
	bool writeFields_;



    // Protected Member Functions


        //- Create the output files
        void createFiles();

        //- Write header for integrated data
        void writeIntegratedHeaderNEW(const word& header, Ostream& os) const;


        void writeIntegratedTrimForces
        (
            	const string& descriptor,
              	const vectorField& fm0,
    		const vectorField& fm1,
    		const vectorField& fm2,
            	autoPtr<OFstream>& osPtr
        ) const;
        
        //- Write TrimForce data
        void writeTrimForces();

	
        //- No copy construct
        TrimForces(const TrimForces&) = delete;

        //- No copy assignment
        void operator=(const TrimForces&) = delete;


public:

    //- Runtime type information
    TypeName("TrimForces");
	
//	mutable float TimeToDelete=0;
//	mutable string test[10];
//	mutable int IT=0;
//	mutable std::string FileName="Force_Cycle";
    // Constructors

        //- Construct from Time and dictionary
        TrimForces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        TrimForces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );


    //- Destructor
    virtual ~TrimForces() = default;


    // Member Functions

        //- Read the TrimForces data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        
//        void DeleteFile() const;

       virtual autoPtr<OFstream> createFile
        (
            const word& name,
            scalar timeValue
        ) const;



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
