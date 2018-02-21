package edu.scripps.yates.pdb;

import javax.json.Json;
import javax.json.JsonObject;
import javax.json.JsonObjectBuilder;

import edu.scripps.yates.pdb.model.Atom3D;

public abstract class JMolAtomReport {
	private final String pdbID;
	private final String uniprotACC;
	private final int positionInUniprot;
	private final Atom3D atom;
	private final Float resolution;
	public final static String sep = "\t";
	private final boolean otherChainsRemoved;
	private final boolean otherMoleculesRemoved;
	private final boolean mutation;
	private final String method;
	private boolean stored = false;

	public JMolAtomReport(String pdbID, Atom3D atom, String uniprotACC, int positionInUniprot, Float resolution,
			boolean otherChainsRemoved, boolean otherMoleculesRemoved, boolean mutation, String method) {
		super();
		this.pdbID = pdbID;
		this.atom = atom;
		this.uniprotACC = uniprotACC;
		this.positionInUniprot = positionInUniprot;
		this.resolution = resolution;
		this.otherChainsRemoved = otherChainsRemoved;
		this.otherMoleculesRemoved = otherMoleculesRemoved;
		this.mutation = mutation;
		this.method = method;
	}

	/**
	 * @return the pdbID
	 */
	public String getPdbID() {
		return pdbID;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(pdbID).append(sep).append(uniprotACC).append(sep).append(positionInUniprot).append(sep)
				.append(resolution).append(sep).append(atom.toString()).append(sep).append(otherChainsRemoved)
				.append(sep).append(otherMoleculesRemoved).append(sep).append(mutation).append(sep).append(method);
		return sb.toString();
	}

	public static String getStaticHeaders() {
		StringBuilder sb = new StringBuilder();
		sb.append("pdb_ID").append(sep).append("uniprot_ACC").append(sep).append("position_in_uniprot").append(sep)
				.append("resolution").append(sep).append(Atom3D.getToStringHeaders()).append(sep)
				.append("OtherChainsRemoved").append(sep).append("OtherMoleculesRemoved").append(sep).append("Mutated")
				.append(sep).append("Method");
		return sb.toString();
	}

	public abstract String getToStringHeaders();

	public String getUniprotACC() {
		return uniprotACC;
	}

	public int getPositionInUniprot() {
		return positionInUniprot;
	}

	public Atom3D getAtom() {
		return atom;
	}

	/**
	 * @return the resolution
	 */
	public Float getResolution() {
		return resolution;
	}

	public JsonObject toJson() {
		String aa = getAtom().getAa();

		int positionInUniprot = getPositionInUniprot();
		int positionInPDB = getAtom().getPositionInPDB();
		float resolution = getResolution();
		String pdbID = getPdbID();
		JsonObjectBuilder objectBuilder = Json.createObjectBuilder().add("aa", aa)
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
		String string = uniprotACC + "-" + getAtom().getAa() + "-" + getPdbID() + "-" + getAtom().getChainID() + "-"
				+ getAtom().getAtomType() + "-" + isOtherChainsRemoved() + "-" + isOtherMoleculesRemoved();
		return string;
	}

	public boolean isStored() {
		return stored;
	}

	public void setStored(boolean stored) {
		this.stored = stored;
	}
}
