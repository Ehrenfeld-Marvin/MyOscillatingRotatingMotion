		string str1;
		string str2;
		string str3;
		string str4;
		string strReadForce;
		std::string strForce;

		string Path=time_.path();
		string PostProcessingPath="/Trim_Forces.dat";
		Path=Path.append(PostProcessingPath);

		std::ifstream TrimForceFile; 
		TrimForceFile.open(Path);
    		float TrimForce_File_Time;
    	

	
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
      					
      					
					str2 = strReadForce.substr(first_tab+1,second_tab-first_tab-1);	//Save second Value in str2	: equals total[x]
					str3 = strReadForce.substr(second_tab+1,third_tab-second_tab-1);	//Save third Value in str3	: equals total[y]
					str4 = strReadForce.substr(third_tab+1);				//Save fourth Value in str4	: equals total[z]
      		
      					if(str2[0]=='+') str2[0]='0';						//ersetzt jedes '+' mit einer '0'
					if(str3[0]=='+') str3[0]='0';
					if(str4[0]=='+') str4[0]='0';
      		
				
					Force_X_Value = std::stod(str2);
					Force_Y_Value = std::stod(str3);
					Force_Z_Value = std::stod(str4);	

					Average_Divider++;							//Zählt wie viele Werte er gespeichert hat für den Mittelwert
		
					Force_Average_X += Force_X_Value;							//Kräfte von einer Oszillation addieren
					Force_Average_Z += Force_Z_Value;
					Force_Average_Y += Force_Y_Value;
				}
			}
      		}
