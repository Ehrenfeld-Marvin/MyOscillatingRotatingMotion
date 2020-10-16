Foam::vector Foam::solidBodyMotionFunctions::MyoscillatingRotatingMotion::Amplitude() const
{

//int ende=time_.endTime().value();
float time=time_.value();
float Nullstelle = (5*2*3.1416)/omega_;
amplitude_.x()=0;
amplitude_.y()=0;


if(time>=changeTime)
{

	amplitude_.z() += 2;
	changeTime += Nullstelle;
	return amplitude_;

}

	
	return amplitude_;
}

