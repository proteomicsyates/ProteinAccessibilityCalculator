package edu.scripps.yates.pdb.util;

import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.Chain;

public class InputParameters {

	private final Map<Character, List<AtomType>> atomTypeMap;
	private final int positionInUniprotProtein;
	private final Chain chain;
	private final String uniprotPeptideSeq;
	private final int positionInPeptide;
	private final String uniprotACC;
	private String pdbID;
	private final boolean removeOtherChains;
	private boolean removeOtherMolecules;

	public InputParameters(Map<Character, List<AtomType>> atomTypeMap, int positionInUniprotProtein, Chain chain,
			String uniprotPeptideSeq, int positionInPeptide, String uniprotACC, boolean removeOtherChains,
			boolean removeOtherMolecules) {
		super();

		this.atomTypeMap = atomTypeMap;
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

	@Override
	public String toString() {
		return "InputParameters [atomTypeMap=" + getAtomTypeMapString() + ", positionInUniprotProtein="
				+ positionInUniprotProtein + ", chain=" + chain + ", uniprotPeptideSeq=" + uniprotPeptideSeq
				+ ", positionInPeptide=" + positionInPeptide + ", uniprotACC=" + uniprotACC + ", pdbID=" + pdbID
				+ ", removeOtherChains=" + removeOtherChains + ", removeOtherMolecules=" + removeOtherMolecules + "]";
	}

	public InputParameters(Map<Character, List<AtomType>> atomTypeMap, String pdbID, String chainID,
			boolean removeOtherChains, boolean removeOtherMolecules) {
		this(atomTypeMap, -1, new Chain(pdbID, chainID + "=0-0", -1.0f), null, -1, null, removeOtherChains,
				removeOtherMolecules);
	}

	public InputParameters(Map<Character, List<AtomType>> atomTypeMap, String pdbID, boolean removeOtherChains,
			boolean removeOtherMolecules) {
		this(atomTypeMap, -1, null, null, -1, null, removeOtherChains, removeOtherMolecules);
		this.pdbID = pdbID;
	}

	/**
	 * @return the atomType
	 */
	public Map<Character, List<AtomType>> getAtomTypeMap() {
		return atomTypeMap;
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
		Chain chain2 = getChain();
		String chainID = null;
		if (chain2 != null) {
			chainID = chain2.getIdentifier();
		}
		return uniprotACC + "-" + getPdbID() + "-" + chainID + "-" + getAtomTypeMapString() + "-" + removeOtherChains
				+ "-" + removeOtherMolecules;
	}

	private String getAtomTypeMapString() {
		StringBuilder sb = new StringBuilder();
		Set<Character> keySet = this.atomTypeMap.keySet();
		for (Character character : keySet) {

			List<AtomType> atomType = atomTypeMap.get(character);
			if (!"".equals(sb.toString())) {
				sb.append("|");
			}
			for (AtomType atomType2 : atomType) {
				sb.append(character + "-" + atomType2);
			}

		}
		return sb.toString();
	}
}
