package edu.scripps.yates.pdb.surface;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.xml.Entry;
import edu.scripps.yates.pdb.model.SurfaceProtein;

/**
 * Class that manages the surface accessibility reports in a parsistence layer,
 * in this case, a text file. Any query over the surface accesibility will be
 * looked up in the persistence layer. If that is not present, it will call to
 * the {@link SurfaceAccesibilityCalculator}
 *
 * @author Salva
 *
 */
public class SurfaceAccessibilityManager {
	private final static Logger log = Logger.getLogger(SurfaceAccessibilityManager.class);
	private static final String FILE_NAME = "surfaces_accessibilities.csv";
	private static final String SEPARATOR = "\t";
	private final SurfaceAccessibilityCalculator calculator;
	private final File surfaceAccessibilityFile;
	private boolean loaded;
	private final Map<String, SurfaceAccessibilityProteinReport> reports = new HashMap<String, SurfaceAccessibilityProteinReport>();
	private boolean calculateIfNotPresent;

	public SurfaceAccessibilityManager(SurfaceAccessibilityCalculator calc) {
		calculator = calc;
		surfaceAccessibilityFile = new File(
				calc.getPdbParserManager().getPdbFileManager().getParentPath().getAbsolutePath() + File.separator
						+ FILE_NAME);
	}

	public Map<String, SurfaceAccessibilityProteinReport> getProteinAccessibilityReportByProtein(
			Collection<SurfaceProtein> proteins) {
		Map<String, SurfaceAccessibilityProteinReport> ret = new HashMap<String, SurfaceAccessibilityProteinReport>();
		for (SurfaceProtein protein : proteins) {
			final SurfaceAccessibilityProteinReport proteinAccessibilityReportByProtein = getProteinAccessibilityReportByProtein(
					protein);
			if (proteinAccessibilityReportByProtein != null) {
				ret.put(protein.getAcc(), proteinAccessibilityReportByProtein);
			}
		}
		return ret;
	}

	public SurfaceAccessibilityProteinReport getProteinAccessibilityReportByProtein(SurfaceProtein protein) {
		// load data from file if not yet
		loadReportsFromFile();
		// look into the map
		if (!reports.containsKey(protein.getAcc()) && calculateIfNotPresent) {
			log.debug("Getting surface accessibilities for protein: " + protein.getAcc());
			final SurfaceAccessibilityProteinReport surfaceAccessibilityProteinReport = calculator
					.getSurfaceAccesibilityFromProtein(protein);
			if (surfaceAccessibilityProteinReport != null) {

				for (Set<SiteSurfaceAccessibilityReport> reports : surfaceAccessibilityProteinReport
						.getAccessibilitiesByPositionInUniprotSeq().values()) {
					for (SiteSurfaceAccessibilityReport surfaceAccessibilityReport : reports) {
						addReport(surfaceAccessibilityReport);
					}
				}
				if (!surfaceAccessibilityProteinReport.getAccessibilitiesByPositionInUniprotSeq().isEmpty()) {
					appendReportToFile(surfaceAccessibilityProteinReport);
				}
			}
		}
		return reports.get(protein.getAcc());
	}

	private void appendReportToFile(SurfaceAccessibilityProteinReport surfaceAccessibilityProteinReport) {
		PrintWriter out = null;
		try {
			boolean writeHeader = false;
			if (!surfaceAccessibilityFile.exists() || surfaceAccessibilityFile.length() == 0l) {
				// write the header = true
				writeHeader = true;
			}
			log.info("Appending report of protein " + surfaceAccessibilityProteinReport.getUniprotACC() + " to file");
			FileWriter fw = new FileWriter(surfaceAccessibilityFile, true);
			BufferedWriter bw = new BufferedWriter(fw);
			out = new PrintWriter(bw);
			if (writeHeader) {
				out.println(SiteSurfaceAccessibilityReport.getToStringHeaders());
			}
			final Map<Integer, Set<SiteSurfaceAccessibilityReport>> positions = surfaceAccessibilityProteinReport
					.getAccessibilitiesByPositionInUniprotSeq();
			List<Integer> sortedPositions = new ArrayList<Integer>();
			sortedPositions.addAll(positions.keySet());
			Collections.sort(sortedPositions);
			for (Integer position : sortedPositions) {
				final Set<SiteSurfaceAccessibilityReport> reports = positions.get(position);
				for (SiteSurfaceAccessibilityReport surfaceAccessibilityReport : reports) {
					out.println(surfaceAccessibilityReport.toString());
				}
			}

			log.info("Report appended to file");
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			out.close();
		}

	}

	/**
	 * Write all the information in memory to the file, overriding the content
	 * of the file
	 */
	public void dumpToFile() {
		PrintWriter out = null;
		try {
			log.info("Writting into file reports for " + reports.size() + " proteins");
			FileWriter fw = new FileWriter(surfaceAccessibilityFile, false);
			BufferedWriter bw = new BufferedWriter(fw);
			out = new PrintWriter(bw);
			// write headers
			out.println(SiteSurfaceAccessibilityReport.getToStringHeaders());
			List<String> sortedAccessions = new ArrayList<String>();
			sortedAccessions.addAll(reports.keySet());
			Collections.sort(sortedAccessions);
			for (String proteinAcc : sortedAccessions) {
				final SurfaceAccessibilityProteinReport surfaceAccessibilityProteinReport = reports.get(proteinAcc);
				final Map<Integer, Set<SiteSurfaceAccessibilityReport>> positions = surfaceAccessibilityProteinReport
						.getAccessibilitiesByPositionInUniprotSeq();
				List<Integer> sortedPositions = new ArrayList<Integer>();
				sortedPositions.addAll(positions.keySet());
				Collections.sort(sortedPositions);
				for (Integer position : sortedPositions) {
					final Set<SiteSurfaceAccessibilityReport> reports = positions.get(position);
					for (SiteSurfaceAccessibilityReport surfaceAccessibilityReport : reports) {
						out.println(surfaceAccessibilityReport.toString());
					}
				}
			}
			log.info("Reports writed at: " + surfaceAccessibilityFile.getAbsolutePath());
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			out.close();
		}
	}

	private void loadReportsFromFile() {
		if (!loaded) {
			BufferedReader br = null;
			try {
				if (surfaceAccessibilityFile.exists()) {
					FileInputStream fstream = new FileInputStream(surfaceAccessibilityFile);
					DataInputStream in = new DataInputStream(fstream);
					br = new BufferedReader(new InputStreamReader(in));
					String strLine;
					// Read File Line By Line
					boolean firstLine = true;
					while ((strLine = br.readLine()) != null) {
						if (firstLine) {
							firstLine = false;
							continue;
						}
						SiteSurfaceAccessibilityReport report = parseFromLine(strLine);
						if (report != null) {
							addReport(report);
						}
					}
					loaded = true;
					log.info("Surface accesibilities numbers loaded from local file for " + reports.size()
							+ " proteins");
				}
			} catch (IOException e) {
				e.printStackTrace();
			} finally {
				if (br != null) {
					try {
						br.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}

		}

	}

	private void addReport(SiteSurfaceAccessibilityReport report) {
		if (reports.containsKey(report.getUniprotACC())) {
			reports.get(report.getUniprotACC()).addSurfaceAccesibilityReport(report);
		} else {
			final Map<String, Entry> annotatedProtein = calculator.getUplr().getAnnotatedProtein(null,
					report.getUniprotACC());
			if (annotatedProtein.containsKey(report.getUniprotACC())) {
				final Entry entry = annotatedProtein.get(report.getUniprotACC());
				final String proteinSequence = SurfaceAccessibilityCalculator.getUniprotProteinSequence(entry);
				if (proteinSequence != null) {
					SurfaceAccessibilityProteinReport proteinReport = new SurfaceAccessibilityProteinReport(
							report.getUniprotACC(), proteinSequence);
					proteinReport.addSurfaceAccesibilityReport(report);
					reports.put(report.getUniprotACC(), proteinReport);
				}
			}
		}
	}

	private SiteSurfaceAccessibilityReport parseFromLine(String strLine) {
		return SiteSurfaceAccessibilityReport.getFromString(strLine);
	}

	public Map<String, SurfaceAccessibilityProteinReport> getSurfaceAccesibilityFromProteins(
			Set<SurfaceProtein> proteins) {
		Map<String, SurfaceAccessibilityProteinReport> ret = new HashMap<String, SurfaceAccessibilityProteinReport>();
		int counter = 0;
		int percentage = 0;
		for (SurfaceProtein protein : proteins) {
			int newPercentage = counter * 100 / proteins.size();
			if (newPercentage != percentage) {
				log.info(newPercentage + "% of proteins analyzed");
				percentage = newPercentage;
			}
			final SurfaceAccessibilityProteinReport report = getProteinAccessibilityReportByProtein(protein);
			if (report != null) {
				ret.put(protein.getAcc(), report);
				// dumpToFile();
			}
			counter++;
		}
		return ret;
	}

	/**
	 * @return the calculateIfNotPresent
	 */
	public boolean isCalculateIfNotPresent() {
		return calculateIfNotPresent;
	}

	/**
	 * @param calculateIfNotPresent
	 *            the calculateIfNotPresent to set
	 */
	public void setCalculateIfNotPresent(boolean calculateIfNotPresent) {
		this.calculateIfNotPresent = calculateIfNotPresent;
	}

}
