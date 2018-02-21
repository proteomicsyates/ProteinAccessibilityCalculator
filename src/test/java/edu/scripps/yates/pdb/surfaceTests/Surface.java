package edu.scripps.yates.pdb.surfaceTests;

public class Surface {
	private final String pdbID;
	private final String chainID;
	private final boolean removeOtherChains;
	private final boolean removeOtherMolecules;
	private final double accessibility;
	private final double experimentalAccessibility;
	private final String psmID;
	private final int atomNumber;

	public Surface(String pdbID, String chainID, boolean removeOtherChains, boolean removeOtherMolecules,
			double accessibility, double experimentalAccessibility, String psmID, int atomNumber) {
		this.pdbID = pdbID;
		this.chainID = chainID;
		this.removeOtherChains = removeOtherChains;
		this.removeOtherMolecules = removeOtherMolecules;
		this.accessibility = accessibility;
		this.experimentalAccessibility = experimentalAccessibility;
		this.psmID = psmID;
		this.atomNumber = atomNumber;
	}

	public String getPdbID() {
		return pdbID;
	}

	public String getChainID() {
		return chainID;
	}

	public boolean isRemoveOtherChains() {
		return removeOtherChains;
	}

	public boolean isRemoveOtherMolecules() {
		return removeOtherMolecules;
	}

	public double getAccessibility() {
		return accessibility;
	}

	public double getExperimentalAccessibility() {
		return experimentalAccessibility;
	}

	public String getPSMID() {
		return psmID;
	}

	@Override
	public String toString() {
		return "Surface [pdbID=" + pdbID + ", chainID=" + chainID + ", removeOtherChains=" + removeOtherChains
				+ ", removeOtherMolecules=" + removeOtherMolecules + ", accessibility=" + accessibility
				+ ", experimentalAccessibility=" + experimentalAccessibility;
	}

	public String getExtendedKey() {
		return "Surface [pdbID=" + pdbID + ", chainID=" + chainID + ", removeOtherChains=" + removeOtherChains
				+ ", removeOtherMolecules=" + removeOtherMolecules + ", accessibility=" + accessibility
				+ ", experimentalAccessibility=" + experimentalAccessibility;
	}

	public String getKey() {
		return "Surface [pdbID=" + pdbID + ", chainID=" + chainID + ", atomno=" + atomNumber;
	}

	public String getPsmID() {
		return psmID;
	}

	public int getAtomNumber() {
		return atomNumber;
	}
}
