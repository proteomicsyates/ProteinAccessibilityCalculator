package edu.scripps.yates.pdb.model;

import java.util.Set;

import gnu.trove.set.hash.THashSet;

public class Peptide {
	private String sequence;
	private final Set<Protein> proteins = new THashSet<Protein>();
	private Double ratio;

	public Peptide(String sequence, double ratio) {
		this.sequence = sequence;
		this.ratio = ratio;
	}

	public Peptide(String sequence) {
		this.sequence = sequence;
	}

	/**
	 * @return the sequence
	 */
	public String getSequence() {
		return sequence;
	}

	/**
	 * @param sequence
	 *            the sequence to set
	 */
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public void addProtein(Protein protein) {
		this.proteins.add(protein);
	}

	public Set<Protein> getProteins() {
		return this.proteins;
	}

	/**
	 * @return the ratio
	 */
	public double getRatio() {
		return ratio;
	}

	/**
	 * @param ratio
	 *            the ratio to set
	 */
	public void setRatio(double ratio) {
		this.ratio = ratio;
	}
}
