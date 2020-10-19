	scalar run_T = time().value();
	float end_T = time().endTime().value();
	float start_T = time().startTime().value();
	float delta_T = time().deltaT().value();

/*	float Force_Average;
	float Force_Sum=0;
	float Force_Max;
	float Force_Min;
*/
	
	int Array_No = (end_T-start_T)/delta_T;		//Anzahl der Arrayfelder bestimmen
	float Force_Array[Array_No];
	
	int iterator=(run_T/delta_T)-1;
	
	
	Force_Array[iterator]=total[0];
	
	int Array_Test[Array_No];
	for(int i=0; i<=Array_No; i++)
	{
	Array_Test[i]=2*i;
	}
	

/*	for (int i=0; i<=Array_No; i++)			//Alle Force Werte in einem Array Speichern
	{
		Force_Array[i]=total[0];
	}
*/
	
/*	for (int j=0; j<=Array_No; j++)			//Alle Force Werte addieren und den Mittelwert bestimmen
	{
		Force_Sum = Force_Sum + Force_Array[j];
	}
	Force_Average = Force_Sum / Array_No;
	
	for (int m=0; m<=Array_No-1; m++)			//Maximum von Force_Array bestimmen
	{
		if(Force_Array[m] > Force_Array[m-1] && Force_Array [m] > Force_Array[m+1])
		{
		Force_Max = Force_Array[m];
		}
		else {}
	}
*/

