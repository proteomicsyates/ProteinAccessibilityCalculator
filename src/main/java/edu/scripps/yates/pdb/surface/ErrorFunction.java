package edu.scripps.yates.pdb.surface;

import org.junit.Test;

public class ErrorFunction {
	private final static double a = -3.52;
	private final static double b = 17.64;

	public static double getErrorValue(double x, double y) {
		double y2 = a * x + b;
		double[] a = { -12, 60 };
		double[] b = { 5, 0 };
		double[] c = { x, y };
		double lineToPointDistance2D = LineToPointDistance2D(a, b, c, false);
		return lineToPointDistance2D;
	}

	@Test
	public void testing() {
		System.out.println(getErrorValue(0, 17));
		System.out.println(getErrorValue(5, 0.04));
		System.out.println(getErrorValue(0, 40));
	}

	// Compute the dot product AB . AC
	private static double DotProduct(double[] pointA, double[] pointB, double[] pointC) {
		double[] AB = new double[2];
		double[] BC = new double[2];
		AB[0] = pointB[0] - pointA[0];
		AB[1] = pointB[1] - pointA[1];
		BC[0] = pointC[0] - pointB[0];
		BC[1] = pointC[1] - pointB[1];
		double dot = AB[0] * BC[0] + AB[1] * BC[1];

		return dot;
	}

	// Compute the cross product AB x AC
	private static double CrossProduct(double[] pointA, double[] pointB, double[] pointC) {
		double[] AB = new double[2];
		double[] AC = new double[2];
		AB[0] = pointB[0] - pointA[0];
		AB[1] = pointB[1] - pointA[1];
		AC[0] = pointC[0] - pointA[0];
		AC[1] = pointC[1] - pointA[1];
		double cross = AB[0] * AC[1] - AB[1] * AC[0];

		return cross;
	}

	// Compute the distance from A to B
	private static double Distance(double[] pointA, double[] pointB) {
		double d1 = pointA[0] - pointB[0];
		double d2 = pointA[1] - pointB[1];

		return Math.sqrt(d1 * d1 + d2 * d2);
	}

	// Compute the distance from AB to C
	// if isSegment is true, AB is a segment, not a line.
	private static double LineToPointDistance2D(double[] pointA, double[] pointB, double[] pointC, boolean isSegment) {
		double dist = CrossProduct(pointA, pointB, pointC) / Distance(pointA, pointB);
		if (isSegment) {
			double dot1 = DotProduct(pointA, pointB, pointC);
			if (dot1 > 0)
				return Distance(pointB, pointC);

			double dot2 = DotProduct(pointB, pointA, pointC);
			if (dot2 > 0)
				return Distance(pointA, pointC);
		}
		return Math.abs(dist);
	}
}
