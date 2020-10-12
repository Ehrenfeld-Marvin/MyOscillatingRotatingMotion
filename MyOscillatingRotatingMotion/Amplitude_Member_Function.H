void Foam::solidBodyMotionFunctions::MyoscillatingRotatingMotion::Amplitude() const
{
	if(time_.value() <= 2)
	{
	amplitude_.z()=10;
	}
	else if(time_.value() > 5)
	{
	amplitude_.z()=1;
	}
	
	return;
}

