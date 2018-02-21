package edu.scripps.yates.pdb.model;

import java.util.Set;

import gnu.trove.set.hash.THashSet;

public class Protein {
	private final String acc;
	private final Set<Peptide> peptides = new THashSet<Peptide>();

	public Protein(String acc) {
		this.acc = acc;
	}

	/**
	 * @return the acc
	 */
	public String getAcc() {
		return acc;
	}

	/**
	 * @return the peptides
	 */
	public Set<Peptide> getPeptides() {
		return peptides;
	}

	public void addPeptide(Peptide peptide) {
		this.peptides.add(peptide);
		peptide.addProtein(this);
	}
}
