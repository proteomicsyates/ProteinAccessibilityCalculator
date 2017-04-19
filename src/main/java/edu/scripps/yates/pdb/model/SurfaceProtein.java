package edu.scripps.yates.pdb.model;

import java.util.HashSet;
import java.util.Set;

public class SurfaceProtein {
	private final String acc;
	private final Set<SurfacePeptide> peptides = new HashSet<SurfacePeptide>();

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
