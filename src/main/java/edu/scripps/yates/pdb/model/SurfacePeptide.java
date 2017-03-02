package edu.scripps.yates.pdb.model;

public class SurfacePeptide {
	private String sequence;
	private double ratio;

	public SurfacePeptide(String sequence, double ratio) {
		super();
		this.sequence = sequence;
		this.ratio = ratio;
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
