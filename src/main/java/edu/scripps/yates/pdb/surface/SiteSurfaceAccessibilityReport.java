package edu.scripps.yates.pdb.surface;

import javax.json.Json;
import javax.json.JsonObject;
import javax.json.JsonObjectBuilder;

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
	private final Float resolution;
	private final static String sep = "\t";
	private final boolean otherChainsRemoved;
	private final boolean otherMoleculesRemoved;
	private final boolean mutation;
	private final String method;
	private boolean stored = false;

	public SiteSurfaceAccessibilityReport(String pdbID, double accesibility, String aa, Atom atom, AtomType atomType,
			int positionInPDB, String uniprotACC, int positionInUniprot, Float resolution, boolean otherChainsRemoved,
			boolean otherMoleculesRemoved, boolean mutation, String method) {
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
		this.otherChainsRemoved = otherChainsRemoved;
		this.otherMoleculesRemoved = otherMoleculesRemoved;
		this.mutation = mutation;
		this.method = method;
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
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append(accesibility).append(sep).append(aa).append(sep).append(pdbID).append(sep).append(positionInPDB)
				.append(sep).append(uniprotACC).append(sep).append(positionInUniprot).append(sep).append(resolution)
				.append(sep).append(atom.toString()).append(sep).append(otherChainsRemoved).append(sep)
				.append(otherMoleculesRemoved).append(sep).append(mutation).append(sep).append(getMethod());
		return sb.toString();
	}

	public static String getToStringHeaders() {
		StringBuilder sb = new StringBuilder();
		sb.append("Surface_accessibility").append(sep).append("AA").append(sep).append("pdb_ID").append(sep)
				.append("position_in_PDB").append(sep).append("uniprot_ACC").append(sep).append("position_in_uniprot")
				.append(sep).append("resolution").append(sep).append(Atom.getToStringHeaders()).append(sep)
				.append("OtherChainsRemoved").append(sep).append("OtherMoleculesRemoved").append(sep).append("Mutated")
				.append(sep).append("Method");
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

			String pdbID = split[2];
			int positionInPDB = Integer.valueOf(split[3]);
			String uniprotACC = split[4];
			int positionInUniprot = Integer.valueOf(split[5]);
			Float resolution = null;
			try {
				resolution = Float.valueOf(split[6]);
			} catch (Exception e) {
			}
			AtomType atomType = AtomType.getByName(split[8]);
			// // get the rest of the splitted items to construct the Atom
			// object
			StringBuilder sb = new StringBuilder();
			String chainID = split[9];
			String atomNumber = split[7];
			sb.append(atomNumber + sep + atomType + sep + chainID);
			Boolean removeOtherChains = Boolean.valueOf(split[10]);
			Boolean removeOtherMolecules = Boolean.valueOf(split[11]);
			Boolean containsMutation = Boolean.valueOf(split[12]);
			String method = split[13];
			// for (int i = 8; i < split.length; i++) {
			// sb.append(split[i]).append(sep);
			// }
			// Atom atom = Atom.getFromString(sb.toString(), aa, positionInPDB);
			Atom atom = Atom.getFromString(sb.toString(), aa, positionInPDB);
			SiteSurfaceAccessibilityReport report = new SiteSurfaceAccessibilityReport(pdbID, accesibility, aa, atom,
					atomType, positionInPDB, uniprotACC, positionInUniprot, resolution, removeOtherChains,
					removeOtherMolecules, containsMutation, method);
			return report;
		} catch (

		Exception e) {
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
	public Float getResolution() {
		return resolution;
	}

	public JsonObject toJson() {
		String aa = getAa();
		double accesibility = getAccesibility();
		int positionInUniprot = getPositionInUniprot();
		int positionInPDB = getPositionInPDB();
		float resolution = getResolution();
		String pdbID = getPdbID();
		JsonObjectBuilder objectBuilder = Json.createObjectBuilder().add("aa", aa).add("accesibility", accesibility)
				.add("positionInUniprot", positionInUniprot).add("positionInPDB", positionInPDB)
				.add("resolution", resolution).add("pdbID", pdbID);
		return objectBuilder.build();
	}

	public boolean isOtherChainsRemoved() {
		return otherChainsRemoved;
	}

	public boolean isOtherMoleculesRemoved() {
		return otherMoleculesRemoved;
	}

	public boolean isMutation() {
		return mutation;
	}

	public String getMethod() {
		return method;
	}

	public String getReportKey() {
		String string = uniprotACC + "-" + getAa() + "-" + getPdbID() + "-" + getAtom().getChainID() + "-"
				+ getAtomType() + "-" + isOtherChainsRemoved() + "-" + isOtherMoleculesRemoved();
		return string;
	}

	public boolean isStored() {
		return stored;
	}

	public void setStored(boolean stored) {
		this.stored = stored;
	}
}
