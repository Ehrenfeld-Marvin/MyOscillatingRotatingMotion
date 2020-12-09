float Time_Intervall_Value = (NuOfOsc*2*3.141592653)/omega_PID;

	if(startTime>0)												//RESTART-BEDINGUNGEN FÜR MIT START =/= NULL
	{
		while(Read_Time_Intervall<startTime) Read_Time_Intervall+=Time_Intervall_Value;			//aktuellen Wert für Zeitintervall erstellen
		
		/*+++++++++++ READING PID BACK UP DATA - AMPLITUDE ++++++*/
		string ReadFile_AMP;
		string FILE_AMP;
    		string Path_AMP=time_.path();
		string FileName_AMP="/Amplitude-BackUp-Data.dat";
		scalar time;
		Path_AMP=Path_AMP.append(FileName_AMP);
    		std::ifstream AMP_BackUp;
    		AMP_BackUp.open(Path_AMP);
    		
    		string str1;
		string str2;
		string str3;
		string str4;
    		
    		while(getline(AMP_BackUp, FILE_AMP))
    		{
    			ReadFile_AMP = FILE_AMP;
    			
    			std::size_t first_tab = ReadFile_AMP.find('\t');			
      			std::size_t second_tab = ReadFile_AMP.find('\t', first_tab+1);	
			std::size_t third_tab = ReadFile_AMP.find('\t', second_tab+1);
    			
    			str1 = ReadFile_AMP.substr(0,first_tab);				
			str2 = ReadFile_AMP.substr(first_tab+1,second_tab-first_tab-1);	
			str3 = ReadFile_AMP.substr(second_tab+1,third_tab-second_tab-1);	
			str4 = ReadFile_AMP.substr(third_tab+1);	
    			
 	   		time = std::stof(str1);
 	   		
 	   		if(time<startTime)
 	   		{
 	   		
    				IST_alt_AMP = std::stod(str2); 
				e_t_SUMME_AMP = std::stod(str3);
				amplitude_y = std::stod(str4);
				amplitude_.z()=amplitude_y;
				break;
			}
    		}
    		
		string ReadFile_PHI;
		string FILE_PHI;
    		string Path_PHI=time_.path();
		string FileName_PHI="/Phi-BackUp-Data.dat";
		Path_PHI=Path_PHI.append(FileName_PHI);
    		std::ifstream PHI_BackUp;
    		PHI_BackUp.open(Path_PHI);
    		
    		while(getline(PHI_BackUp, FileName_PHI))
    		{
    			ReadFile_PHI = FileName_PHI;
    			
    			std::size_t first_tab = ReadFile_PHI.find('\t');			
      			std::size_t second_tab = ReadFile_PHI.find('\t', first_tab+1);	
			std::size_t third_tab = ReadFile_PHI.find('\t', second_tab+1);
    			
    			str1 = ReadFile_PHI.substr(0,first_tab);				
			str2 = ReadFile_PHI.substr(first_tab+1,second_tab-first_tab-1);	
			str3 = ReadFile_PHI.substr(second_tab+1,third_tab-second_tab-1);	
			str4 = ReadFile_PHI.substr(third_tab+1);	
    			
 	   		time = std::stof(str1);
 	   		
 	   		if(time<startTime)
 	   		{
 	   		
    				IST_alt_PHI = std::stod(str2); 
				e_t_SUMME_PHI = std::stod(str3);
				phi_y = std::stod(str4);
				break;
			}
    		}
    	
	/*++++++++++++++++++++++++++++++++++++++++++++*/
	}
	
	if(startTime==0)							//RESTART- UND START-BEDINUNGEN FÜR START = 0
	{
		amplitude_.z()=amplitude_begin;
		phi_Value=0;
		
		amplitude_ramp=amplitude_.z();
		phi_ramp=0;
		
		IST_alt_AMP=0;
		e_t_alt_AMP=0;
		e_t_SUMME_AMP=0;
		
		IST_alt_PHI=0;
		e_t_alt_PHI=0;
		e_t_SUMME_PHI=0;
		
		Force_Average_X=0;
		Force_Average_Y=0;
		Force_Average_Z=0;
		
		Read_Time_Intervall=Time_Intervall_Value;
		remove("Phi-BackUp-Data.dat");
		remove("Amplitude-BackUp-Data.dat");
	}
	amplitude_.x()=0;
	amplitude_.y()=0;
