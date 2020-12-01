Foam::vector Foam::solidBodyMotionFunctions::TrimMotion::Amplitude(scalar omega) const
{
/*++++++++++++++++++++++++++++	VARIABLEN DEKLARIEREN		+++++++++++++++++++++++++++*/

	scalar time=time_.value();
	scalar delta_T=time_.deltaT().value();
//	float end_T=time_.endTime().value();
	scalar Time_Intervall_Value = (NuOfOsc*2*3.141592653)/omega_PID;
	string str1;
	string str2;
	string str3;
	string str4;
	string strReadForce;
	
	float TrimForce_File_Time;
	float Force_X_Value=0;
	float Force_Y_Value=0;
	float Force_Z_Value=0;

	
	int Average_Divider=0;
	int Average_Divider_Y=0;
	

	string Path=time_.path();
	string PostProcessingPath="/Trim_Forces.dat";
	Path=Path.append(PostProcessingPath);

//	amplitude_.x()=0;
//	amplitude_.y()=0;


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//	if(time==0) Read_Time_Intervall=Time_Intervall_Value; //+delta_T;				//Die Force Datei soll erst nach einem ganzen Cylce gelesen werden
										//+delta_T wegen Zeitverzögerung zwischen schreiben und lesen
/*##########################		FORCE_CYCLE LESEN	##################################*/
	
	

	if(time > Read_Time_Intervall)
	{

		std::string strForce;
		std::ifstream TrimForceFile; 
		TrimForceFile.open(Path);
    	
    		Force_Average_X = 0;							//Kräfte von einer Oszillation addieren
		Force_Average_Z = 0;
		Force_Average_Y = 0;
	
	    	while(getline(TrimForceFile, strForce))						//ließt die Force_Cycles Datei linie für linie
	      	{
      			strReadForce = strForce;
      				      		
      			if(strReadForce[0]!='#')							//Überspringt Header
      			{
      				
      				std::size_t first_tab = strReadForce.find('\t');			//Find place of first tab
      				std::size_t second_tab = strReadForce.find('\t', first_tab+1);	//Find place of second tab
				std::size_t third_tab = strReadForce.find('\t', second_tab+1);	//Find place of third tab
				
      				str1 = strReadForce.substr(0,first_tab);					//Save first Value in str1	: equals Time
      				
      				TrimForce_File_Time = std::stof(str1);					//Converting STRING to FLOAT
      				
      				if(TrimForce_File_Time>=(Read_Time_Intervall-Time_Intervall_Value))
      				{
      					
      					cout	<< "\nRead_Time_Intervall: " << Read_Time_Intervall
      						<< "\n" << strReadForce;
      					
					str2 = strReadForce.substr(first_tab+1,second_tab-first_tab-1);	//Save second Value in str2	: equals total[x]
					str3 = strReadForce.substr(second_tab+1,third_tab-second_tab-1);	//Save third Value in str3	: equals total[y]
					str4 = strReadForce.substr(third_tab+1);				//Save fourth Value in str4	: equals total[z]
      		
      					if(str2[0]=='+') str2[0]='0';						//ersetzt jedes '+' mit einer '0'
					if(str3[0]=='+') str3[0]='0';
					if(str4[0]=='+') str4[0]='0';
      		
				
					Force_X_Value = std::stod(str2);
					Force_Y_Value = std::stod(str3);
					Force_Z_Value = std::stod(str4);	

					cout	<<"\nForce_X_Value: " << Force_X_Value
						<<"\nForce_Y_Value: " << Force_Y_Value
						<<"\nForce_Z_Value: " << Force_Z_Value;

					Average_Divider++;							//Zählt wie viele Werte er gespeichert hat für den Mittelwert
		
					Force_Average_X += Force_X_Value;							//Kräfte von einer Oszillation addieren
					Force_Average_Z += Force_Z_Value;
					Force_Average_Y += Force_Y_Value;
					cout	<<"\nForce_Average_X: " << Force_Average_X
						<<"\nForce_Average_Y: " << Force_Average_Y
						<<"\nForce_Average_Z: " << Force_Average_Z;
				}
			}
      		}
      		
/*#################################################################################*/

/*++++++++++++++++++++++++++++++	MITTELWERTE BESTIMMEN	++++++++++++++++++++++++++++++*/

      		if(Average_Divider==0) Average_Divider=1;					//If Bedingung wird bei t=0 trotzdem aufgerugen, Divider wäre hier =0
   		if(Average_Divider_Y==0) Average_Divider_Y=1;									

      		Force_Average_X /= (Average_Divider);					//Mittelwerte von n-Cycles bestimmen
//    		Force_Average_Y /= (Average_Divider_Y);
		Force_Average_Y /= (Average_Divider);
     		Force_Average_Z /= (Average_Divider);
		
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/		

/*++##########################	INTEGRAL SUMME		##############################*/			//VERWENDEN WENN NICHT DAS GESAMTE INTEGRAL VERWENDET WERDEN SOLL
		
/*		if(I_iterator!=I_intervall) 						//Hier noch IST_SUMME korriegieren zu e_t_SUMME
		{
			I_vector.push_back(Force_Average_Y);
			I_iterator++;
		}
		else
		{
			I_vector.erase(I_vector.begin());
			I_vector.push_back(Force_Average_Y);
		}
		
		IST_SUMME=0;
		for(auto i=I_vector.begin(); i!=I_vector.end(); i++) e_t_SUMME+=*i;
*/		

/*++######################################################################################*/
		
		
/*++++++++++++++++++++++++++++	PID REGLER	+++++++++++++++++++++++++++++++++++++*/
		
		cout 	<<"TARGET " << Target << "\n"
		 	<<"FORCE_AVERAGE_Y: " << Force_Average_Y << "\n"
		 	<<"E_T_SUMME: " << e_t_SUMME << "\n"
		 	<<"IST_alt " << IST_alt << "\n"
		 	<<"e_t_SUMME " << e_t_SUMME << "\n"
		 	<<"y_t " << y_t << "\n";
		
		IST = Force_Average_Y;							
				 	cout <<"IST " << IST << "\n";
		e_t=Target-IST;
		e_t_alt=Target-IST_alt;
		e_t_SUMME+=e_t; 
		
		P=K_P*e_t;
		I=K_I*delta_T*e_t_SUMME;
		D=K_D*((e_t-e_t_alt)/delta_T);							//FEHLER
		
		y_t=P+I+D;			//STELLGRÖSSE		
		
		IST_alt=IST;

		
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
		
		
		amplitude_.z() = y_t;
		Read_Time_Intervall += Time_Intervall_Value; //+delta_T);					//+delta_T wegen Zeitverzögerung zwischen schreiben und lesen
	
		cout 	<<"TARGET " << Target << "\n"
		 	<<"FORCE_AVERAGE_Y: " << Force_Average_Y << "\n"
		 	<<"E_T_SUMME: " << e_t_SUMME << "\n"
		 	<<"IST " << IST << "\n"
		 	<<"IST_alt_neu " << IST_alt << "\n"
		 	<<"e_t " << e_t << "\n"
		 	<<"e_t_SUMME " << e_t_SUMME << "\n"
		 	<<"y_t " << y_t << "\n"
		 	<<"amplitude_ " << amplitude_.z() << "\n";

		writePIDdata
		(
		Target,			
		K_P,			
		K_I,
		K_D,
		IST,
		IST_alt,
		e_t_SUMME,
		e_t,
		e_t_alt,
		y_t,
		P,
		I,
		D,
		Force_Average_X,
		Force_Average_Y,
		Force_Average_Z,
		time
		);	
		
		Back_up_Data_PID
		(
			time,
			IST_alt,
			e_t_SUMME,
			y_t
		);
				
		return amplitude_;
	}
			cout 	<<"TARGET " << Target << "\n"
		 	<<"FORCE_AVERAGE_Y: " << Force_Average_Y << "\n"
		 	<<"E_T_SUMME: " << e_t_SUMME << "\n"
		 	<<"IST " << IST << "\n"
		 	<<"IST_alt_neu " << IST_alt << "\n"
		 	<<"e_t " << e_t << "\n"
		 	<<"e_t_SUMME " << e_t_SUMME << "\n"
		 	<<"y_t " << y_t << "\n"
		 	<<"amplitude_ " << amplitude_.z() << "\n";
	return amplitude_;
}

