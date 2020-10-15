int Foam::solidBodyMotionFunctions::MyoscillatingRotatingMotion::Amplitude() const
{

//int ende=time_.endTime().value();
float time=time_.value();
int change;
int test;

amplitude_.x()=0;
amplitude_.y()=0;


if(time_.value()==2)
{
	change =2;
//	amplitude_.z()=1;
	test =10;
	
	return test;
}
/*else if (time==5)
	{
	test = 5;
	return test;
	}
*/

/*
if(time>=change)
{

//	amplitude_.z() += 1;
	test += 1;
	change += 2;
	return test;

}*/

	
	return test;
}

