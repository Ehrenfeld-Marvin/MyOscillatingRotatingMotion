vector ampli;
ampli.x()=0;
ampli.y()=0;

if(time_.value() < 2)
{	
	ampli.z()=5;
}	
else if (time_.value() >= 2 && time_.value() <= 5)
{
	ampli.z()=10;
}
else	
{
	ampli.z()=1;
}


/*ampli.z()=1;				//Hier bleibt dir Amplitude bei 1			
if(time_.value() == 2)		
{				
	ampli.z()=5;		
}				
else{}				
*/
