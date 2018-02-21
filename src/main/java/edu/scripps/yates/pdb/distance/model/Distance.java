package edu.scripps.yates.pdb.distance.model;

import edu.scripps.yates.pdb.JMolAtomReport;
import edu.scripps.yates.pdb.model.Atom3D;

public class Distance {
	private final Atom3D originAtom;
	private final Atom3D endAtom;
	private final Double distanceInAmstrongs;

	public Distance(Atom3D originalAtom, Atom3D endAtom, Double distanceInAmstrongs) {
		this.distanceInAmstrongs = distanceInAmstrongs;
		this.endAtom = endAtom;
		this.originAtom = originalAtom;
	}

	public Atom3D getOriginAtom() {
		return originAtom;
	}

	public Atom3D getEndAtom() {
		return endAtom;
	}

	public Double getDistanceInAmstrongs() {
		return distanceInAmstrongs;
	}

	@Override
	public String toString() {
		return distanceInAmstrongs + JMolAtomReport.sep + originAtom + JMolAtomReport.sep + endAtom;
	}

	public static String getStaticHeaders() {
		StringBuilder sb = new StringBuilder();
		sb.append("Distance (A)").append(JMolAtomReport.sep).append(Atom3D.getToStringHeaders())
				.append(JMolAtomReport.sep).append(Atom3D.getToStringHeaders());
		return sb.toString();
	}

	public static Distance fromString(String distanceString) {
		try {
			String[] split = distanceString.split(JMolAtomReport.sep);
			double distance = Double.valueOf(split[0]);
			StringBuilder sb = new StringBuilder();
			int index = 1;
			for (int i = 0; i < Atom3D.numElementsInPrint(); i++) {
				sb.append(split[index++]).append(JMolAtomReport.sep);
			}
			Atom3D atomOrigin = Atom3D.getFromString(sb.toString());
			StringBuilder sb2 = new StringBuilder();
			for (int i = 0; i < Atom3D.numElementsInPrint(); i++) {
				sb2.append(split[index++]).append(JMolAtomReport.sep);
			}
			Atom3D atomEnd = Atom3D.getFromString(sb2.toString());
			Distance ret = new Distance(atomOrigin, atomEnd, distance);
			return ret;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
}
