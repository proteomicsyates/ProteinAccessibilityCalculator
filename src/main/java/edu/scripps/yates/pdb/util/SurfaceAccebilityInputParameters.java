package edu.scripps.yates.pdb.util;

import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.Chain;

public class SurfaceAccebilityInputParameters {
	private final String aa;
	private final AtomType atomType;
	private final int positionInUniprotProtein;
	private final Chain chain;
	private final String uniprotPeptideSeq;
	private final int positionInPeptide;
	private final String uniprotACC;
	private String pdbID;
	private final boolean removeOtherChains;
	private boolean removeOtherMolecules;

	public SurfaceAccebilityInputParameters(String aa, AtomType atomType, int positionInUniprotProtein, Chain chain,
			String uniprotPeptideSeq, int positionInPeptide, String uniprotACC, boolean removeOtherChains,
			boolean removeOtherMolecules) {
		super();
		this.aa = aa;
		this.atomType = atomType;
		this.positionInUniprotProtein = positionInUniprotProtein;
		this.chain = chain;
		this.uniprotPeptideSeq = uniprotPeptideSeq;
		this.positionInPeptide = positionInPeptide;
		this.uniprotACC = uniprotACC;
		if (chain != null) {
			pdbID = chain.getPdbID();
		}
		this.removeOtherChains = removeOtherChains;
		this.removeOtherMolecules = removeOtherMolecules;

	}

	public SurfaceAccebilityInputParameters(String aa, AtomType atomType, String pdbID, String chainID,
			boolean removeOtherChains, boolean removeOtherMolecules) {
		this(aa, atomType, -1, new Chain(pdbID, chainID + "=0-0", -1.0f), null, -1, null, removeOtherChains,
				removeOtherMolecules);
	}

	public SurfaceAccebilityInputParameters(String aa, AtomType atomType, String pdbID, boolean removeOtherChains,
			boolean removeOtherMolecules) {
		this(aa, atomType, -1, null, null, -1, null, removeOtherChains, removeOtherMolecules);
		this.pdbID = pdbID;
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
	 * @return the positionInUniprotProtein
	 */
	public int getPositionInUniprotProtein() {
		return positionInUniprotProtein;
	}

	/**
	 * @return the chain
	 */
	public Chain getChain() {
		return chain;
	}

	public String getPdbID() {
		return pdbID;
	}

	/**
	 * @return the uniprotPeptideSeq
	 */
	public String getUniprotPeptideSeq() {
		return uniprotPeptideSeq;
	}

	/**
	 * @return the positionInPeptide
	 */
	public int getPositionInPeptide() {
		return positionInPeptide;
	}

	/**
	 * @return the uniprotACC
	 */
	public String getUniprotACC() {
		return uniprotACC;
	}

	public boolean isRemoveOtherChains() {
		return removeOtherChains;
	}

	public boolean isRemoveOtherMolecules() {
		return removeOtherMolecules;
	}

	public String getReportKey() {
		return uniprotACC + "-" + getAa() + "-" + getPdbID() + "-" + getChain().getIdentifier() + "-" + getAtomType()
				+ "-" + removeOtherChains + "-" + removeOtherMolecules;
	}
}
