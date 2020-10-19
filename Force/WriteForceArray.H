void Foam::functionObjects::Myforces::writeForceArray
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

	scalar run_T = time().value();
	float end_T = time().endTime().value();
	float start_T = time().startTime().value();
	float delta_T = time().deltaT().value();
	int Array_No = (end_T-start_T)/delta_T;		//Anzahl der Arrayfelder bestimmen
	vector Force_Array[Array_No];
	int iterator=(run_T/delta_T)-1;
	int Divide = iterator +1;
	
	vector Average = Zero;
	vector Force_Sum = Zero;
   	
   	Force_Array[iterator]=total;
   	
   	for(int m=iterator; m>=0; m--)
   	{
   		Force_Sum += Force_Array[m];
   	}   	
   	Average = Force_Sum/Divide;
   	
   	if (writeToFile())
	{
        	Ostream& os = osPtr();

        	writeCurrentTime(os);	

		os 	<< tab << total
			<< tab << Average;
	
	for(int j=iterator; j>=0; j--)
	{
		os 	<< tab << Force_Array[j] 
			<< tab << Force_Array[j][0] << endl << tab << tab;
			
	}

        os  << endl;
    }
}
