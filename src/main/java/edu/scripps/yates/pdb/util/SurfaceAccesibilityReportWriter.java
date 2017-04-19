package edu.scripps.yates.pdb.util;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.pdb.model.SurfacePeptide;
import edu.scripps.yates.pdb.surface.SiteSurfaceAccessibilityReport;
import edu.scripps.yates.pdb.surface.SurfaceAccessibilityProteinReport;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.strings.StringUtils;

public class SurfaceAccesibilityReportWriter {
	/**
	 * Print in a file the information of a surface accessibility sites present
	 * in a particular peptide of a protein
	 *
	 * @param fw
	 * @param accession
	 * @param ratios
	 * @param ratioScoreNames
	 * @param peptideSequence
	 * @param surfaceAccesibilityProteinReport
	 * @param aa
	 * @param b
	 * @throws IOException
	 */
	public static void printReportForPsm(FileWriter fw, Ratio ratio, String psmID, String peptideSequence,
			SurfaceAccessibilityProteinReport surfaceAccesibilityProteinReport, String aa,
			boolean printOnlyTheMostAccessibleSite) throws IOException {
		final String uniprotProteinSequence = surfaceAccesibilityProteinReport.getUniprotProteinSequence();
		final int positionOfPeptideInProtein = uniprotProteinSequence.indexOf(peptideSequence);

		List<Integer> positionsOfAAInPeptide = StringUtils.allPositionsOf(peptideSequence, aa);
		for (Integer positionOfAAInPeptide : positionsOfAAInPeptide) {
			int positionOfAAInProtein = positionOfAAInPeptide + positionOfPeptideInProtein;
			final Set<SiteSurfaceAccessibilityReport> surfaceAccesibilityReports = surfaceAccesibilityProteinReport
					.getSurfaceAccessibilityReportsBySite(positionOfAAInProtein);
			if (surfaceAccesibilityReports != null) {
				if (!printOnlyTheMostAccessibleSite) {
					for (SiteSurfaceAccessibilityReport surfaceAccessibilityReport : surfaceAccesibilityReports) {
						printReportForSite(fw, surfaceAccesibilityProteinReport.getUniprotACC(), ratio, psmID,
								peptideSequence, positionOfAAInPeptide, aa, surfaceAccessibilityReport);
					}
				} else {
					SiteSurfaceAccessibilityReport mostAccessibleSite = null;
					double accessibility = -1;
					for (SiteSurfaceAccessibilityReport surfaceAccessibilityReport : surfaceAccesibilityReports) {
						if (surfaceAccessibilityReport.getAccesibility() > accessibility) {
							accessibility = surfaceAccessibilityReport.getAccesibility();
							mostAccessibleSite = surfaceAccessibilityReport;
						}
					}
					printReportForSite(fw, surfaceAccesibilityProteinReport.getUniprotACC(), ratio, psmID,
							peptideSequence, positionOfAAInPeptide, aa, mostAccessibleSite);
				}
			}
		}

	}

	/**
	 * Print in a file the information of a surface accessibility sites present
	 * in a particular peptide of a protein
	 *
	 * @param fw
	 * @param accession
	 * @param ratios
	 * @param ratioScoreNames
	 * @param peptideSequence
	 * @param surfaceAccesibilityProteinReport
	 * @param aa
	 * @param b
	 * @throws IOException
	 */
	public static void printReportForPeptide(Writer fw, SurfacePeptide peptide,
			SurfaceAccessibilityProteinReport surfaceAccesibilityProteinReport, String aa,
			boolean printOnlyTheMostAccessibleSite) throws IOException {
		final String uniprotProteinSequence = surfaceAccesibilityProteinReport.getUniprotProteinSequence();
		final int positionOfPeptideInProtein = uniprotProteinSequence.indexOf(peptide.getSequence());

		List<Integer> positionsOfAAInPeptide = StringUtils.allPositionsOf(peptide.getSequence(), aa);
		for (Integer positionOfAAInPeptide : positionsOfAAInPeptide) {
			int positionOfAAInProtein = positionOfAAInPeptide + positionOfPeptideInProtein;
			final Set<SiteSurfaceAccessibilityReport> surfaceAccesibilityReports = surfaceAccesibilityProteinReport
					.getSurfaceAccessibilityReportsBySite(positionOfAAInProtein);
			if (surfaceAccesibilityReports != null) {
				if (!printOnlyTheMostAccessibleSite) {
					for (SiteSurfaceAccessibilityReport surfaceAccessibilityReport : surfaceAccesibilityReports) {
						printReportForSite(fw, surfaceAccesibilityProteinReport.getUniprotACC(), peptide,
								positionOfAAInPeptide, aa, surfaceAccessibilityReport);
					}
				} else {
					SiteSurfaceAccessibilityReport mostAccessibleSite = null;
					double accessibility = -1;
					for (SiteSurfaceAccessibilityReport surfaceAccessibilityReport : surfaceAccesibilityReports) {
						if (surfaceAccessibilityReport.getAccesibility() > accessibility) {
							accessibility = surfaceAccessibilityReport.getAccesibility();
							mostAccessibleSite = surfaceAccessibilityReport;
						}
					}
					printReportForSite(fw, surfaceAccesibilityProteinReport.getUniprotACC(), peptide,
							positionOfAAInPeptide, aa, mostAccessibleSite);
				}
			}
		}

	}

	private static void printReportForSite(Writer fw, String accession, SurfacePeptide peptide,
			int positionOfAAInPeptide, String aa, SiteSurfaceAccessibilityReport surfaceAccesibilityReport)
					throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append(peptide.getSequence()).append("\t").append(positionOfAAInPeptide).append("\t");

		String ratioString = printRatio(peptide.getRatio());
		sb.append(ratioString).append("\t");
		sb.append(surfaceAccesibilityReport);
		sb.append("\n");
		fw.write(sb.toString());

	}

	public static String getReportHeader() {
		return "Peptide_seq\tPosition_in_peptide\tRatio\t" + SiteSurfaceAccessibilityReport.getToStringHeaders();
	}

	private static void printReportForSite(FileWriter fw, String accession, Ratio ratio, String psmID,
			String peptideSequence, int positionOfAAInPeptide, String aa,
			SiteSurfaceAccessibilityReport surfaceAccesibilityReport) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append(accession).append("\t").append(psmID).append("\t").append(peptideSequence).append("\t")
				.append(positionOfAAInPeptide).append("\t").append(surfaceAccesibilityReport).append("\t");

		final double value = ratio.getValue();
		String ratioString = printRatio(value);
		sb.append(ratioString);
		sb.append("\n");
		fw.write(sb.toString());

	}

	private static String printRatio(double value) {
		StringBuilder sb = new StringBuilder();
		if (Double.isInfinite(value)) {
			if (Double.compare(value, Double.NEGATIVE_INFINITY) == 0) {
				sb.append(-1000);
			} else {
				sb.append(+1000);
			}
		} else {
			sb.append(value);
		}
		return sb.toString();
	}

}
