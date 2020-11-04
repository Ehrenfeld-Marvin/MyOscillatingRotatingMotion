void Foam::functionObjects::Myforces::DeleteFile
(

) const
{
	float run_T=time().value();
	float delta=time().deltaT().value();
//	float TimeInterval= (NuOfOsc*2*3.1416)/omega;
	float TimeInterval= 3;
	
	writeForceArray
	(
      	"Force_Cycle",
      	coordSys_.localVector(Myforce_[0]),
      	coordSys_.localVector(Myforce_[1]),
      	coordSys_.localVector(Myforce_[2]),
	Force_Cycle
	);

	if(TimeToDelete==0) TimeToDelete = TimeInterval+delta;		//Werte sollen erst nach einem Cycle gelöscht werden
										//+delta wegen Zeitverzögerung zwischen schreiben und lesen

	if(run_T>=TimeToDelete)
	{

		TimeToDelete += (TimeInterval+delta);	
//		FileName += std::to_string(IT);
//		IT++;
		remove("/home/marvin/OpenFOAM/marvin-v2006/applications/MyOscillatingRotatingMotion/Case/postProcessing/MyForces_Dir/0/Force_Cycle.dat");
		Force_Cycle= createFile(FileName);
		writeIntegratedHeaderNEW("Force_Cycle", Force_Cycle());    		
		
		cout<< "\n\nFORCE: DELETING FORCE FILE. TIME:" << run_T << "\n\n";
	}


}
