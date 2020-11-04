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
	
	float run_T = time().value();
	float end_T = time().endTime().value();
	float start_T = time().startTime().value();
	float delta_T = time().deltaT().value();
	int Array_No = (end_T-start_T)/delta_T;		//Anzahl der Arrayfelder bestimmen
	vector Force_Array[Array_No];
	int iterator=(run_T/delta_T)-1;


	
	
   	
   	Force_Array[iterator]=total;
   	
   	if (writeToFile())
	{
		Ostream& os = osPtr();
        	
//              writeCurrentTime(os);	

		os 	<< setprecision(5) << run_T;
					
	if(total[0]>0)	os << tab << setprecision(15) << "+" << total[0];
	else os << setprecision(15) << tab << total[0];
	if(total[1]>0)	os << tab << "+" << total[1];
	else os << tab << total[1];
	if(total[2]>0)	os << tab << "+" << total[2];
	else os << tab << total[2];

	cout <<"\n\nFORCE: Time to write Data at current Time:" << run_T << "\n\n";
	
	os << endl;
	
    	}
}
