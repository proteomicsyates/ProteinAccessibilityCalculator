package edu.scripps.yates.pdb.model;

import java.util.Set;

import gnu.trove.set.hash.THashSet;

public class SurfaceProtein {
	private final String acc;
	private final Set<SurfacePeptide> peptides = new THashSet<SurfacePeptide>();

	public SurfaceProtein(String acc) {
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
	public Set<SurfacePeptide> getPeptides() {
		return peptides;
	}

	public void addPeptide(SurfacePeptide peptide) {
		this.peptides.add(peptide);
		peptide.addProtein(this);
	}
}
