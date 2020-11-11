Foam::vector Foam::solidBodyMotionFunctions::MyoscillatingRotatingMotion::Amplitude() const
{
/*++++++++++++++++++++++++++++	VARIABLEN DEKLARIEREN		+++++++++++++++++++++++++++*/

	float time=time_.value();
	float delta=time_.deltaT().value();
	float end_t=time_.endTime().value();
	
//	float Nullstelle = (NuOfOsc*2*3.141592653)/omega_;
	float Nullstelle = 4;
	string str1;
	string str2;
	string str3;
	string str4;
	string strReadForce;
	
	float forceTime;
	float total_X=0;
	float total_Y=0;
	float total_Z=0;

	int Average_Divider=0;
	int Average_Divider_Y=0;
	


//	float Sinus_Iterations = 2*3.141592653/(omega_*delta);

	amplitude_.x()=0;
	amplitude_.y()=0;

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

	if(time==0) changeTime=Nullstelle+delta;				//Die Force Datei soll erst nach einem ganzen Cylce gelesen werden
										//+delta wegen Zeitverzögerung zwischen schreiben und lesen
										
//	cout <<"\n\n\n Sinus Iterations " << Sinus_Iterations << "\n";
//	cout <<"Sinus Iterator " << Sinus_Iterator << "\n\n\n";
		
	if(time>=changeTime)							//Bedingung die Force Datei nur nach n-Zyklen des Sinus zu lesen
	{

/*##########################		FORCE_CYCLE LESEN	##################################*/

		std::string strForce;
		std::ifstream MyFile; 
		MyFile.open("/home/marvin/OpenFOAM/marvin-v2006/applications/MyOscillatingRotatingMotion/Case/postProcessing/MyForces_Dir/0/Force_Cycle.dat");
    	
	    	while(getline(MyFile, strForce))						//ließt die Force_Cycles Datei linie für linie
	      	{
	      		strReadForce = strForce;
      			
      			if(strReadForce[0]!='#')
      			{
	      			std::size_t first_tab = strReadForce.find('\t');			//Find place of first tab
				std::size_t second_tab = strReadForce.find('\t', first_tab+1);	//Find place of second tab
				std::size_t third_tab = strReadForce.find('\t', second_tab+1);	//Find place of third tab
				str1 = strReadForce.substr(0,first_tab);				//Save first Value in str1	: equals Time
				str2 = strReadForce.substr(first_tab+1,second_tab-first_tab-1);	//Save second Value in str2	: equals total[x]
				str3 = strReadForce.substr(second_tab+1,third_tab-second_tab-1);	//Save third Value in str3	: equals total[y]
				str4 = strReadForce.substr(third_tab+1);				//Save fourth Value in str4	: equals total[z]
      		
      				if(str2[0]=='+') str2[0]='0';						//ersetzt jedes '+' mit einer '0'
				if(str3[0]=='+') str3[0]='0';
				if(str4[0]=='+') str4[0]='0';
      		
				forceTime = std::stod(str1);						//Converting STRING to FLOAT
				total_X = std::stod(str2);
				total_Y = std::stod(str3);
				total_Z = std::stod(str4);	

//				Sinus_Iterator++;							//Zählt wieviele Linien die Schleife gelesen hat um später mit n-Linien zu rechnen 
				Average_Divider++;							//Zählt wie viele Werte er gespeichert hat für den Mittelwert
		
				Average_X += total_X;							//Kräfte von einer Oszillation addieren
				Average_Z += total_Z;

//				if(Sinus_Iterator>Sinus_Iterations/2) // geteilt durch 2		//Nur die erste Hälfte des Sinus wird addiert
//				{
					Average_Y += total_Y;
//					Average_Divider_Y++;
//				}
			}
      		}
/*#################################################################################*/

/*++++++++++++++++++++++++++++++	MITTELWERTE BESTIMMEN	++++++++++++++++++++++++++++++*/

      		if(Average_X==0) Average_Divider=1;					//If Bedingung wird bei t=0 trotzdem aufgerugen, Divider wäre hier =0
   		if(Average_Divider_Y==0) Average_Divider_Y=1;									

      		cout <<"\n\n\nDIVIDER_Y: " << Average_Divider_Y <<"\n\n\n";
      		
      		Average_X /= (Average_Divider);					//Mittelwerte von n-Cycles bestimmen
//    		Average_Y /= (Average_Divider_Y);
		Average_Y /= (Average_Divider);
     		Average_Z /= (Average_Divider);
		
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/		

/*++##########################	INTEGRAL SUMME		##############################*/		
		
/*		if(I_iterator!=I_intervall) 						//Hier noch IST_SUMME korriegieren zu e_t_SUMME
		{
			I_vector.push_back(Average_Y);
			I_iterator++;
		}
		else
		{
			I_vector.erase(I_vector.begin());
			I_vector.push_back(Average_Y);
		}
		
		IST_SUMME=0;
		for(auto i=I_vector.begin(); i!=I_vector.end(); i++) e_t_SUMME+=*i;
*/		

/*++######################################################################################*/
		
		
/*++++++++++++++++++++++++++++	PID REGLER	+++++++++++++++++++++++++++++++++++++*/
		
		cout 	<<"\n\n\nIST_alt " << IST_alt << "\n";
		
		IST = Average_Y;							
		
		e_t=Target-IST;
		e_t_alt=Target-IST_alt;
		e_t_SUMME+=e_t; 
		
		P=K_P*e_t;
		I=K_I*delta*e_t_SUMME;
		D=K_D*((e_t-e_t_alt)/delta)	;							//FEHLER
		
		y_t=P+I+D;			//STELLGRÖSSE		
		
//		if(y_t<0)y_t=0;								//Amplitude darf nicht negativ sein????

		IST_alt=IST;

		
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
		
		
		amplitude_.z() = y_t;
		changeTime += (Nullstelle+delta);					//+delta wegen Zeitverzögerung zwischen schreiben und lesen
		
		cout 	<<"TARGET " << Target << "\n"
		 	<<"AVERGAE Y: " << Average_Y << "\n"
		 	<<"IST_SUMME: " << e_t_SUMME << "\n"
//		 	<<"I_iterator: " << I_iterator << "\n"
		 	<<"IST " << IST << "\n"
		 	<<"IST_alt_neu " << IST << "\n"
		 	<<"e_t " << e_t << "\n"
		 	<<"y_t " << y_t << "\n"
		 	<<"K_P " << K_P << "\n"
		 	<<"K_I " << K_I << "\n"
		 	<<"K_D " << K_D << "\n"
		 	<<"P " << P << "\n"
		 	<<"I " << I << "\n"
		 	<<"D " << D << "\n"
		 	<<"amplitude_ " << amplitude_.z() << "\n";

					
		
//		if(Sinus_Iterator>Sinus_Iterations-1 )Sinus_Iterator=0;
		
		return amplitude_;
	}
	
 


	
	
	return amplitude_;
}

