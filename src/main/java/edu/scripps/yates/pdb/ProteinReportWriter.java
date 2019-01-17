package edu.scripps.yates.pdb;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.pdb.model.Peptide;
import edu.scripps.yates.pdb.surface.SurfaceReport;
import edu.scripps.yates.utilities.proteomicsmodel.Ratio;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.list.array.TIntArrayList;

public class ProteinReportWriter {
	/**
	 * Print in a file the information of a surface accessibility sites present in a
	 * particular peptide of a protein
	 *
	 * @param fw
	 * @param accession
	 * @param ratios
	 * @param ratioScoreNames
	 * @param peptideSequence
	 * @param proteinReport
	 * @param aa
	 * @param b
	 * @throws IOException
	 */
	public static void printReportForPsm(FileWriter fw, Ratio ratio, String psmID, String peptideSequence,
			ProteinReport proteinReport, char aa, boolean printOnlyTheMostAccessibleSite) throws IOException {
		final String uniprotProteinSequence = proteinReport.getUniprotProteinSequence();
		final int positionOfPeptideInProtein = uniprotProteinSequence.indexOf(peptideSequence);

		final TIntArrayList positionsOfAAInPeptide = StringUtils.allPositionsOf(peptideSequence, aa);
		for (final int positionOfAAInPeptide : positionsOfAAInPeptide.toArray()) {
			final int positionOfAAInProtein = positionOfAAInPeptide + positionOfPeptideInProtein;
			final Set<JMolAtomReport> reports = proteinReport.getReportsBySite(positionOfAAInProtein);
			if (reports != null) {
				if (!printOnlyTheMostAccessibleSite) {
					for (final JMolAtomReport report : reports) {
						printReportForSite(fw, proteinReport.getUniprotACC(), ratio, psmID, peptideSequence,
								positionOfAAInPeptide, report);
					}
				} else {
					JMolAtomReport bestReport = null;
					double accessibility = -1;
					for (final JMolAtomReport report : reports) {
						if (report instanceof SurfaceReport) {
							final SurfaceReport surfaceReport = (SurfaceReport) report;
							if (surfaceReport.getAccessibility() != null
									&& surfaceReport.getAccessibility() > accessibility) {
								accessibility = surfaceReport.getAccessibility();
								bestReport = surfaceReport;
							}
						} else {
							throw new UnsupportedOperationException("This has to be implemented");
						}
					}
					printReportForSite(fw, proteinReport.getUniprotACC(), ratio, psmID, peptideSequence,
							positionOfAAInPeptide, bestReport);
				}
			}
		}

	}

	/**
	 * Print in a file the information of a surface accessibility sites present in a
	 * particular peptide of a protein
	 *
	 * @param fw
	 * @param proteinReport
	 * @throws IOException
	 */
	public static void printReportForPDB(Writer fw, ProteinReport proteinReport) throws IOException {

		final List<JMolAtomReport> reports = proteinReport.getReports();
		if (reports != null) {
			for (final JMolAtomReport surfaceAccessibilityReport : reports) {
				printReport(fw, surfaceAccessibilityReport);
			}
		}
	}

	/**
	 * Print in a file the information of a surface accessibility sites present in a
	 * particular peptide of a protein
	 *
	 * @param fw
	 * @param accession
	 * @param ratios
	 * @param ratioScoreNames
	 * @param peptideSequence
	 * @param proteinReport
	 * @param aa
	 * @param b
	 * @throws IOException
	 */
	public static void printReportForPeptide(Writer fw, Peptide peptide, ProteinReport proteinReport, char aa,
			boolean printOnlyTheMostAccessibleSite) throws IOException {
		final String uniprotProteinSequence = proteinReport.getUniprotProteinSequence();
		final int positionOfPeptideInProtein = uniprotProteinSequence.indexOf(peptide.getSequence());

		final TIntArrayList positionsOfAAInPeptide = StringUtils.allPositionsOf(peptide.getSequence(), aa);
		for (final int positionOfAAInPeptide : positionsOfAAInPeptide.toArray()) {
			final int positionOfAAInProtein = positionOfAAInPeptide + positionOfPeptideInProtein;
			final Set<JMolAtomReport> reports = proteinReport.getReportsBySite(positionOfAAInProtein);
			if (reports != null) {
				if (!printOnlyTheMostAccessibleSite) {
					for (final JMolAtomReport surfaceAccessibilityReport : reports) {
						printReportForSite(fw, peptide, positionOfAAInPeptide, surfaceAccessibilityReport);
					}
				} else {
					JMolAtomReport bestReport = null;
					double accessibility = -1;
					for (final JMolAtomReport report : reports) {
						if (report instanceof SurfaceReport) {
							final SurfaceReport vsurfaceReport = (SurfaceReport) report;
							if (vsurfaceReport.getAccessibility() != null
									&& vsurfaceReport.getAccessibility() > accessibility) {
								accessibility = vsurfaceReport.getAccessibility();
								bestReport = vsurfaceReport;
							}
						} else {
							throw new UnsupportedOperationException("not implemented yet");
						}
					}
					printReportForSite(fw, peptide, positionOfAAInPeptide, bestReport);
				}
			}
		}

	}

	private static void printReportForSite(Writer fw, Peptide peptide, int positionOfAAInPeptide, JMolAtomReport report)
			throws IOException {
		final StringBuilder sb = new StringBuilder();
		sb.append(peptide.getSequence()).append("\t").append(positionOfAAInPeptide).append("\t");

		final String ratioString = printRatio(peptide.getRatio());
		sb.append(ratioString).append("\t");
		sb.append(report);
		sb.append("\n");
		fw.write(sb.toString());
		fw.flush();
	}

	private static void printReport(Writer fw, JMolAtomReport report) throws IOException {
		final StringBuilder sb = new StringBuilder();
		sb.append(report);
		sb.append("\n");
		fw.write(sb.toString());
		fw.flush();
	}

	public static String getPeptideReportHeader(String suffix) {
		return "Peptide_seq\tPosition_in_peptide\tRatio\t" + suffix;
	}

	private static void printReportForSite(FileWriter fw, String accession, Ratio ratio, String psmID,
			String peptideSequence, int positionOfAAInPeptide, JMolAtomReport report) throws IOException {
		final StringBuilder sb = new StringBuilder();
		sb.append(accession).append("\t").append(psmID).append("\t").append(peptideSequence).append("\t")
				.append(positionOfAAInPeptide).append("\t").append(report).append("\t");

		final double value = ratio.getValue();
		final String ratioString = printRatio(value);
		sb.append(ratioString);
		sb.append("\n");
		fw.write(sb.toString());

	}

	private static String printRatio(Double value) {
		if (value == null) {
			return "";
		}
		final StringBuilder sb = new StringBuilder();
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
