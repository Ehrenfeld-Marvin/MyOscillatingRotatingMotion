Foam::vector Foam::solidBodyMotionFunctions::MyoscillatingRotatingMotion::Amplitude() const
{

//int ende=time_.endTime().value();

//ISstream TestFile("Test.dat");
/*if(!TestFile.is_open())
{
TestFile.open("Test.dat");
}
*/

float time=time_.value();
float Nullstelle = (5*2*3.1416)/omega_;
amplitude_.x()=0;
amplitude_.y()=0;
	
std::ifstream MyFile;
MyFile.open ("/home/marvin/OpenFOAM/marvin-v2006/applications/MyOscillatingRotatingMotion/Case/postProcessing/MyForces_Dir/0/Force_Array.dat");


	if(time>=changeTime)
	{
//		std::ifstream MyFile ("/home/marvin/OpenFOAM/marvin-v2006/applications/MyOscillatingRotatingMotion/Case/postProcessing/MyForces_Dir/0/Force_Array.dat");
    	
    		while(getline(MyFile, readFile))
      		{
      			readAmplitude = std::stoi(readFile);
      		}
		
		amplitude_.z() = readAmplitude;
		changeTime += Nullstelle;
		return amplitude_;
      	}

/*
	if(time>=changeTime)
	{
	
		amplitude_.z() += 3;
		changeTime += Nullstelle;
		return amplitude_;

	}
*/
	
	return amplitude_;
}

