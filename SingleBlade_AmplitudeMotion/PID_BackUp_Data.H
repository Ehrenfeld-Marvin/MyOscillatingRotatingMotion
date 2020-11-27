


void Foam::solidBodyMotionFunctions::MyoscillatingRotatingMotion::Back_up_Data_PID
(
	double IST_alt,
	double e_t_SUMME,
	double e_t_alt
) const
{
		string Path=time_.path();
		string FileName="/Back_Up_Data_PID.dat";
		Path=Path.append(FileName);
		std::ofstream Back_Up_Data_PID;
		Back_Up_Data_PID.open(Path, std::ios_base::trunc);
    	
		Back_Up_Data_PID	<< IST_alt << "\n" 
					<< e_t_SUMME << "\n"
					<< e_t_alt;
	
}



