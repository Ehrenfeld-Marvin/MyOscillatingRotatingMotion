float Time_Intervall_Value = (NuOfOsc*2*3.141592653)/omega_PID;

	if(startTime>0)												//RESTART-BEDINGUNGEN FÜR MIT START =/= NULL
	{
		while(Read_Time_Intervall<startTime) Read_Time_Intervall+=Time_Intervall_Value;			//aktuellen Wert für Zeitintervall erstellen
		
		/*+++++++++++ READING PID BACK UP DATA ++++++*/
		string ReadFile;
		string FILE;
    		string Path=time_.path();
		string FileName="/Back_Up_Data_PID.dat";
		scalar time;
		Path=Path.append(FileName);
    		std::ifstream Back_Up_Data_PID;
    		Back_Up_Data_PID.open(Path);
    		
    		string str1;
		string str2;
		string str3;
		string str4;
    		
    		while(getline(Back_Up_Data_PID, FILE))
    		{
    			ReadFile = FILE;
    			
    			std::size_t first_tab = ReadFile.find('\t');			
      			std::size_t second_tab = ReadFile.find('\t', first_tab+1);	
			std::size_t third_tab = ReadFile.find('\t', second_tab+1);
    			
    			str1 = ReadFile.substr(0,first_tab);				
			str2 = ReadFile.substr(first_tab+1,second_tab-first_tab-1);	
			str3 = ReadFile.substr(second_tab+1,third_tab-second_tab-1);	
			str4 = ReadFile.substr(third_tab+1);	
    			
 	   		time = std::stof(str1);
 	   		
 	   		if(time<startTime)
 	   		{
 	   		
    				IST_alt = std::stod(str2); 
				e_t_SUMME = std::stod(str3);
				y_t = std::stod(str4);
				amplitude_.z()=y_t;
				break;
			}
    		}
    	
	/*++++++++++++++++++++++++++++++++++++++++++++*/
	}
	
	if(startTime==0)							//RESTART- UND START-BEDINUNGEN FÜR START = 0
	{
		IST_alt=0;
		e_t_SUMME=0;
		Force_Average_X=0;
		Force_Average_Y=0;
		Force_Average_Z=0;
		Read_Time_Intervall=Time_Intervall_Value;
		amplitude_.z()=amplitude_begin;
		remove("Back_Up_Data_PID.dat");
	}
	amplitude_.x()=0;
	amplitude_.y()=0;
