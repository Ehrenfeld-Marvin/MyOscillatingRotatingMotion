Foam::scalar Foam::solidBodyMotionFunctions::TrimMotion::PID_REGLER
(
	double Target,
	double Average,
	float K_P,
	float K_I,
	float K_D,
	double* e_t_SUMME,
	double* IST_alt,
	double* e_t_alt,
	double* P,
	double* I,
	double* D,
	double* IST,
	double* e_t
) const
{
		double y_t;

		
		*IST = Average;							
		*e_t = Target - *IST;
		*e_t_SUMME += *e_t; 
		
		*P = K_P * *e_t;
		*I = K_I * Time_Intervall_Value * *e_t_SUMME;
		*D = K_D * ((*e_t - *e_t_alt)/Time_Intervall_Value);
		
		y_t= *P + *I + *D;					
		
		*IST_alt=*IST;
		*e_t_alt=Target-*IST;
		
		return y_t;
}
		

