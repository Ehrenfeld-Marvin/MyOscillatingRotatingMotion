#include "IFstream.H"

Foam::autoPtr<Foam::OFstream> Foam::functionObjects::TrimForces::createFile
(
    const word& name,
    scalar timeValue
) const
{
	autoPtr<OFstream> osPtr;

	if (Pstream::master() && writeToFile_)
	{
//		if (useUserTime_)
//		{
//		timeValue = fileObr_.time().timeToUserTime(timeValue);
//		}

//		const word timeName = Time::timeName(timeValue);

//		fileName outputDir(baseFileDir()/prefix_/timeName);
//		mkDir(outputDir);

		word fName(name);
		
		
		// Check if file already exists
//		IFstream is(outputDir/(fName + ".dat"));
		IFstream is(fName + ".dat");
		
		
		
		
		if (is.good())
		{
//			fName = fName + "_" + timeName;
//			file = fopen("Trim_Forces.dat", "a");
			osPtr.set(new OFstream(fName + ".dat", IOstreamOption(), true));
		}

//		osPtr.reset(new OFstream(outputDir/(fName + ".dat")));
		osPtr.reset(new OFstream(fName + ".dat", IOstreamOption(), true));
		
		
        
		if (!osPtr->good())
		{
			FatalIOErrorInFunction(osPtr()) << "Cannot open file"
			<< exit(FatalIOError);
		}

			initStream(osPtr());
	}

		return osPtr;
}




/*
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeFile::writeFile(const writeFile& wf)
:
    fileObr_(wf.fileObr_),
    prefix_(wf.prefix_),
    fileName_(wf.fileName_),
    filePtr_(nullptr),
    writeToFile_(wf.writeToFile_),
    useUserTime_(wf.useUserTime_),
{}


Foam::functionObjects::writeFile::writeFile
(
    const objectRegistry& obr,
    const fileName& prefix,
    const word& name,
    const bool writeToFile
)
:
    fileObr_(obr),
    prefix_(prefix),
    fileName_(name),
    filePtr_(nullptr),
    writeToFile_(writeToFile),
    useUserTime_(true),
{}


Foam::functionObjects::writeFile::writeFile
(
    const objectRegistry& obr,
    const fileName& prefix,
    const word& name,
    const dictionary& dict,
    const bool writeToFile
)
:
    writeFile(obr, prefix, name, writeToFile)
{
    read(dict);

    if (writeToFile_)
    {
        filePtr_ = createFile(fileName_);
    }
}
*/


