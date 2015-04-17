package Stat;

import java.util.Random;

public class ChiSquareRand
{
	public static float next(int df)
	{
		double 
			v = 0,
			x = 0;
		Random gen = new Random();
		for(int i = 0; i < df; i++)
		{
			 x = gen.nextGaussian();
			 v += x * x;
		}
		return (float)v;
	}
	
	public static double next2(double freedom)
	{
		     	double u = 0,v = 0,z = 0,zz = 0,r = 0, b = 0, vm = 0, vp = 0, vd = 0;
		     	double freedom_in = -1.0;
		      	//if( a < 1 )  return (-1.0); // Check for invalid input value
		      
		      	if (freedom == 1.0) {
		      		for(;;) {
		      			u = Math.random();
		      			v = Math.random() * 0.857763884960707;
		      			z = v / u;
		     			if (z < 0) continue;
		     			zz = z * z;
		     			r = 2.5 - zz;
		     			if (z < 0.0) r = r + zz * z / (3.0 * z);
		     			if (u < r * 0.3894003915) return(z*z);
		     			if (zz > (1.036961043 / u + 1.4)) continue;
		     			if (2.0 * Math.log(u) < (- zz * 0.5 )) return(z*z);
		     		}
		     	}
		     	else {
		     		if (freedom != freedom_in) {
		     			b = Math.sqrt(freedom - 1.0);
		     			vm = - 0.6065306597 * (1.0 - 0.25 / (b * b + 1.0));
		     			vm = (-b > vm) ? -b : vm;
		     			vp = 0.6065306597 * (0.7071067812 + b) / (0.5 + b);
		     			vd = vp - vm;
		     			freedom_in = freedom;
		     		}
		     		for(;;) 
		     		{
		     			u = Math.random();
		     			v = Math.random() * vd + vm;
		     			z = v / u;
		     			if (z < -b) continue;
		     			zz = z * z;
		     			r = 2.5 - zz;
		     			if (z < 0.0) r = r + zz * z / (3.0 * (z + b));
		     			if (u < r * 0.3894003915) return((z + b)*(z + b));
		     			if (zz > (1.036961043 / u + 1.4)) continue;
		     			if (2.0 * Math.log(u) < (Math.log(1.0 + z / b) * b * b - zz * 0.5 - z * b)) return((z + b)*(z + b));
		     		}
		     	}
	}

}
