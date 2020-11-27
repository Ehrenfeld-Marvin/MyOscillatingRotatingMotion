scalar Foam::solidBodyMotionFunctions::TrimMotion::PID-REGLER
(
	Average
	delta_T
) const
{


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
		

		
		IST = Average;							
		
		e_t=Target-IST;
		e_t_alt=Target-IST_alt;
		e_t_SUMME+=e_t; 
		
		P=K_P*e_t;
		I=K_I*delta_T*e_t_SUMME;
		D=K_D*((e_t-e_t_alt)/delta_T);							//FEHLER
		
		y_t=P+I+D;			//STELLGRÖSSE		
		
		IST_alt=IST;
		
		return y_t;
}
		

