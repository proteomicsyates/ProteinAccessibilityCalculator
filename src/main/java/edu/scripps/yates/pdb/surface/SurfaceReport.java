package edu.scripps.yates.pdb.surface;

import edu.scripps.yates.pdb.JMolAtomReport;
import edu.scripps.yates.pdb.model.Atom3D;
import edu.scripps.yates.pdb.util.InputParameters;

public class SurfaceReport extends JMolAtomReport {
	private final Double accessibility;

	public SurfaceReport(Double accesibility, String pdbID, Atom3D atom, String uniprotACC, int positionInUniprot,
			Float resolution, boolean otherChainsRemoved, boolean otherMoleculesRemoved, boolean mutation,
			String method) {
		super(pdbID, atom, uniprotACC, positionInUniprot, resolution, otherChainsRemoved, otherMoleculesRemoved,
				mutation, method);
		this.accessibility = accesibility;
	}

	public SurfaceReport(Double accesibility, InputParameters inputParameters, Atom3D atom, Float resolution,
			boolean mutation, String experimentalMethod) {
		this(accesibility, inputParameters.getPdbID(), atom, inputParameters.getUniprotACC(),
				inputParameters.getPositionInUniprotProtein(), resolution, inputParameters.isRemoveOtherChains(),
				inputParameters.isRemoveOtherMolecules(), mutation, experimentalMethod);
	}

	/**
	 * @return the accessibility
	 */
	public Double getAccessibility() {
		return accessibility;
	}

	/**
	 * Parse a string as opposite function of toString()
	 * sb.append("pdb_ID").append(sep).append("uniprot_ACC").append(sep).append("position_in_uniprot").append(sep)
	 * .append("resolution").append(sep).append(Atom3D.getToStringHeaders()).append(sep)
	 * .append("OtherChainsRemoved").append(sep).append("OtherMoleculesRemoved").append(sep).append("Mutated")
	 * .append(sep).append("Method");
	 * 
	 * @param string
	 * @return
	 */
	public static SurfaceReport getFromString(String string) {
		try {
			final String[] split = string.split(sep);
			int index = 2;

			final String pdbID = split[index++];
			index++;
			final String uniprotACC = split[index++];
			final int positionInUniprot = Integer.valueOf(split[index++]);
			Float resolution = null;
			try {
				resolution = Float.valueOf(split[index++]);
			} catch (final Exception e) {
			}
			// // get the rest of the splitted items to construct the
			// Atom object
			final StringBuilder sb = new StringBuilder();
			for (int i = 0; i < Atom3D.numElementsInPrint(); i++) {
				sb.append(split[index++]).append(sep);
			}

			final Atom3D atom = Atom3D.getFromString(sb.toString());
			final Boolean removeOtherChains = Boolean.valueOf(split[index++]);
			final Boolean removeOtherMolecules = Boolean.valueOf(split[index++]);
			final Boolean containsMutation = Boolean.valueOf(split[index++]);
			final String method = split[index++];
			Double accessibility = null;
			try {
				accessibility = Double.valueOf(split[index++]);
			} catch (final Exception e) {
			}
			final SurfaceReport report = new SurfaceReport(accessibility, pdbID, atom, uniprotACC, positionInUniprot,
					resolution, removeOtherChains, removeOtherMolecules, containsMutation, method);
			return report;
		} catch (

		final Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		String accessibilityString = "-";
		if (accessibility != null) {
			accessibilityString = String.valueOf(accessibility);
		}
		sb.append(super.toString()).append(sep).append(accessibilityString);
		return sb.toString();
	}

	@Override
	public String getToStringHeaders() {
		final StringBuilder sb = new StringBuilder();
		sb.append(JMolAtomReport.getStaticHeaders()).append(sep).append("Surface_accessibility");
		return sb.toString();
	}
}
