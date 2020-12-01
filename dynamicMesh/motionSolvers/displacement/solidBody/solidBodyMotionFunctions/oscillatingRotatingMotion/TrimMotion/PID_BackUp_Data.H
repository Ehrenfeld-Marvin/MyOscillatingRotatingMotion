


void Foam::solidBodyMotionFunctions::TrimMotion::Back_up_Data_PID
(
	scalar time,
	double IST_alt,
	double e_t_SUMME,
	double y_t
) const
{
		string Path=time_.path();
		string FileName="/Back_Up_Data_PID.dat";
		Path=Path.append(FileName);
		std::ofstream Back_Up_Data_PID;
		Back_Up_Data_PID.open(Path, std::ios_base::app);
    	
		Back_Up_Data_PID	<< time << tab << IST_alt << tab << e_t_SUMME << tab << y_t <<"\n";
	
}


