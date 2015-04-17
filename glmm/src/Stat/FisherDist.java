/*
@(#)FisherDistribution.java
Copyright (C) 2001-2007  Kyle Siegrist, Dawn Duehring
Department of Mathematical Sciences
University of Alabama in Huntsville

This program is part of Virtual Laboratories in Probability and Statistics,
http://www.math.uah.edu/stat/, a project partially supported by the
National Science Foundation under grant number DUE-0089377.

This program is licensed under a Creative Commons License. Basically, you are free to copy,
distribute, and modify this program, and to make commercial use of the program.
However you must give proper attribution.
See http://creativecommons.org/licenses/by/2.0/ for more information.
*/
package Stat;
public class FisherDist
{
	private int numeratorDegrees, denominatorDegrees;
	private double num, den;
	/**
	* This general constructor creates a new Fisher distribution with a
	* specified number of degrees of freedom in the numerator and denominator.
	* @param n the numerator degrees of freedom
	* @param d the denominator degrees of freedom
	*/
	public FisherDist(int n, int d){
		setParameters(n, d);
	}
	/**
	* This default constructor creates a new Fisher distribution with 5
	* degrees of freedom in numerator and denominator.
	*/
	public FisherDist()
	{
		this(5, 5);
	}

	/**
	* This method sets the parameters, the degrees of freedom in the numerator
	* and denominator. Additionally, the normalizing constant and default domain
	* are computed.
	* @param n the numerator degrees of freedom
	* @param d the denominator degrees of freedom
	*/
	public void setParameters(int n, int d)
	{
		numeratorDegrees = n; denominatorDegrees = d;
		num = (double)numeratorDegrees; den = (double)denominatorDegrees;
	}

	/**
	* This method computes the cumulative distribution function in terms of
	* the beta CDF.
	* @param x a real number
	* @return the cumulative probability at x
	*/
	public double getCDF(double x)
	{
		double u = den / (den + num * x);
		if (x < 0) return 0;
		else return 1 - Probability.beta(0.5 * den, 0.5 * num, u);
	}

	/**
	* This method returns the numerator degrees of freedom.
	* @return the numerator degrees of freedom
	*/
	public double getNumeratorDegrees(){
		return numeratorDegrees;
	}

	/**
	* This method sets the numerator degrees of freedom.
	* @param n the numerator degrees of freedom
	*/
	public void setNumeratorDegrees(int n){
		setParameters(n, denominatorDegrees);
	}

	/**
	* This method gets the denominator degrees of freedom.
	* @return the denominator degrees of freedom
	*/
	public double getDenominatorDegrees(){
		return denominatorDegrees;
	}

	/**
	* This method sets the denominator degrees of freedom.
	* @param d the denominator degrees of freedom
	*/
	public void setDenominatorDegrees(int d){
		setParameters(numeratorDegrees, d);
	}
}


