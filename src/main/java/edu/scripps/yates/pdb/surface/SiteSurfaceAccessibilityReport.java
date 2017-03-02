package edu.scripps.yates.pdb.surface;

import edu.scripps.yates.pdb.model.Atom;
import edu.scripps.yates.pdb.model.AtomType;

public class SiteSurfaceAccessibilityReport {
	private final double accesibility;
	private final String aa;
	private final AtomType atomType;
	private final String pdbID;
	private final int positionInPDB;
	private final String uniprotACC;
	private final int positionInUniprot;
	private final Atom atom;
	private final float resolution;
	private final static String sep = "\t";

	public SiteSurfaceAccessibilityReport(String pdbID, double accesibility, String aa, Atom atom, AtomType atomType,
			int positionInPDB, String uniprotACC, int positionInUniprot, float resolution) {
		super();
		this.pdbID = pdbID;
		this.accesibility = accesibility;
		this.aa = aa;
		this.atom = atom;
		this.atomType = atomType;
		this.positionInPDB = positionInPDB;
		this.uniprotACC = uniprotACC;
		this.positionInUniprot = positionInUniprot;
		this.resolution = resolution;
	}

	/**
	 * @return the accesibility
	 */
	public double getAccesibility() {
		return accesibility;
	}

	/**
	 * @return the aa
	 */
	public String getAa() {
		return aa;
	}

	/**
	 * @return the atomType
	 */
	public AtomType getAtomType() {
		return atomType;
	}

	/**
	 * @return the pdbID
	 */
	public String getPdbID() {
		return pdbID;
	}

	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append(accesibility).append(sep).append(aa).append(sep).append(pdbID).append(sep).append(positionInPDB)
				.append(sep).append(uniprotACC).append(sep).append(positionInUniprot).append(sep).append(resolution)
				.append(sep).append(atom.toString());
		return sb.toString();
	}

	public static String getToStringHeaders() {
		StringBuilder sb = new StringBuilder();
		sb.append("Surface_accessibility").append(sep).append("AA").append(sep).append("pdb_ID").append(sep)
				.append("position_in_PDB").append(sep).append("uniprot_ACC").append(sep).append("position_in_uniprot")
				.append(sep).append("resolution").append(sep).append(Atom.getToStringHeaders());
		return sb.toString();
	}

	/**
	 * Parse a string as oposite function of toString()
	 *
	 * @param string
	 * @return
	 */
	public static SiteSurfaceAccessibilityReport getFromString(String string) {
		try {
			final String[] split = string.split(sep);

			double accesibility = Double.valueOf(split[0]);
			String aa = split[1];
			AtomType atomType = AtomType.getByName(split[2]);
			String pdbID = split[3];
			int positionInPDB = Integer.valueOf(split[4]);
			String uniprotACC = split[5];
			int positionInUniprot = Integer.valueOf(split[6]);
			float resolution = Float.valueOf(split[7]);
			// get the rest of the splitted items to construct the Atom object
			StringBuilder sb = new StringBuilder();
			for (int i = 8; i < split.length; i++) {
				sb.append(split[i]).append(sep);
			}
			Atom atom = Atom.getFromString(sb.toString(), aa, positionInPDB);
			SiteSurfaceAccessibilityReport report = new SiteSurfaceAccessibilityReport(pdbID, accesibility, aa, atom,
					atomType, positionInPDB, uniprotACC, positionInUniprot, resolution);
			return report;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * @return the positionInPDB
	 */
	public int getPositionInPDB() {
		return positionInPDB;
	}

	public String getUniprotACC() {
		return uniprotACC;
	}

	public int getPositionInUniprot() {
		return positionInUniprot;
	}

	public Atom getAtom() {
		return atom;
	}

	/**
	 * @return the resolution
	 */
	public float getResolution() {
		return resolution;
	}
}
