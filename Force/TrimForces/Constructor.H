	if(startTime>0)
	{
		restart_without_deleting=true;
		
  	  	string ReadTrimForce;
  	  	string Time_str;
  	  	float Time_float;
  	  	std::string TrimForceFile_line;
		std::fstream TrimForceFile; 
		TrimForceFile.open("Trim_Forces.dat");
		std::ofstream temp;
		temp.open("temp.dat");
    	
		while(getline(TrimForceFile, TrimForceFile_line))						//ließt die Force_Cycles Datei linie für linie
		{
			ReadTrimForce = TrimForceFile_line;
			
			if(ReadTrimForce[0]!='#')							//Überspringt Header
			{
      				
				std::size_t first_tab = ReadTrimForce.find('\t');			//Find place of first tab
				Time_str = ReadTrimForce.substr(0,first_tab);					//Save first Value in str1	: equals Time
				Time_float = std::stod(Time_str);					//Converting STRING to FLOAT

				if(Time_float<=startTime) temp << TrimForceFile_line << "\n";
			}
			else temp << TrimForceFile_line << "\n";
			
		}
		
		temp.close();
		TrimForceFile.close();
		remove("Trim_Forces.dat");
		rename("temp.dat","Trim_Forces.dat");
      		
	}
    
//	string Path=time().path();
//	string PostProcessingPath="/postProcessing/TrimForces_Dir/0/Trim_Forces.dat";
//	Path=Path.append(PostProcessingPath);
	
	if(startTime==0)
	{
		remove("Trim_Forces.dat");
	}
