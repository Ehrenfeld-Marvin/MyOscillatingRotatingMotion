Foam::vector Foam::solidBodyMotionFunctions::MyoscillatingRotatingMotion::Amplitude() const
{
/*++++++++++++++++++++++++++++	VARIABLEN DEKLARIEREN		+++++++++++++++++++++++++++*/

	float time=time_.value();
	float delta=time_.deltaT().value();
	float end_t=time_.endTime().value();
	
//	float Nullstelle = (NuOfOsc*2*3.141592653)/omega_;
	float Nullstelle = 3;
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
	


//	I_vector[I_intervall+1]=0;

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
/*	
//		I_array[I_iterator]=Average_Y;
		if(I_iterator<I_intervall) 
		{
			IST_SUMME+=Average_Y;
			I_vector[I_iterator]=Average_Y;
			cout <<"\n\n\nI_vector["<< I_iterator << "] " << I_vector[I_iterator] <<"\n";
			I_iterator++;

		}
		else 
		{
			IST_SUMME=0;
			for(int position=1; position <= I_intervall; position++)
			{
				
//				if(position < (I_intervall-1)) 
				I_vector[position]=I_vector[position+1];
//				else I_vector[position]=Average_Y;
//				IST_SUMME += I_vector[position];
			}
			
		}
*/		
		
/*############################	PID REGLER	###################################*/
		
		
		IST = Average_Y;							
		IST_SUMME+=Average_Y; 
		e_t=Target-IST;								//FEHLER
		y_t=K_P*e_t+K_I*delta*IST_SUMME+K_D*(IST-IST_alt)*delta;			//STELLGRÖSSE
//		if(y_t<0)y_t=0;								//Amplitude darf nicht negativ sein????
		IST_alt=IST;

		
/*#########################################################################################*/
		
		
		amplitude_.z() = y_t;
		changeTime += (Nullstelle+delta);					//+delta wegen Zeitverzögerung zwischen scheiben und lesen
		
		cout <<"\n\n\nTARGET " << Target << "\n";
		cout <<"AVERGAE Y: " << Average_Y << "\n";
		cout <<"IST_SUMME: " << IST_SUMME << "\n";
//		cout <<"I_iterator: " << I_iterator << "\n";		
		cout <<"IST " << IST << "\n";
		cout <<"e_t " << e_t << "\n";
		cout <<"y_t " << y_t << "\n";
		cout <<"amplitude_ " << amplitude_.z() << "\n";


					
		
//		if(Sinus_Iterator>Sinus_Iterations-1 )Sinus_Iterator=0;
		
		return amplitude_;
	}
	



	
	
	return amplitude_;
}

