package edu.scripps.yates.pdb.surface;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.channels.FileLock;
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
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;

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
		calculator.setManager(this);
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
		if (protein.getAcc().equals("P09211")) {
			System.out.println(protein);
		}
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
			return surfaceAccessibilityProteinReport;
		}
		return null;

	}

	public SurfaceAccessibilityProteinReport getProteinAccessibilityReportByProtein(String reportKey) {
		// load data from file if not yet
		loadReportsFromFile();
		// look into the map

		return reports.get(reportKey);

	}

	private void appendReportToFile(SurfaceAccessibilityProteinReport surfaceAccessibilityProteinReport) {
		PrintWriter out = null;
		FileLock fileLocker = null;
		try {
			boolean writeHeader = false;
			if (!surfaceAccessibilityFile.exists() || surfaceAccessibilityFile.length() == 0l) {
				// write the header = true
				writeHeader = true;
			}

			FileOutputStream fw = new FileOutputStream(surfaceAccessibilityFile, true);
			fileLocker = fw.getChannel().lock();
			out = new PrintWriter(fw);
			if (writeHeader) {
				out.println(SiteSurfaceAccessibilityReport.getToStringHeaders());
			}
			boolean firstOne = false;
			final Map<Integer, Set<SiteSurfaceAccessibilityReport>> positions = surfaceAccessibilityProteinReport
					.getAccessibilitiesByPositionInUniprotSeq();
			List<Integer> sortedPositions = new ArrayList<Integer>();
			sortedPositions.addAll(positions.keySet());
			Collections.sort(sortedPositions);
			for (Integer position : sortedPositions) {
				final Set<SiteSurfaceAccessibilityReport> reports = positions.get(position);
				for (SiteSurfaceAccessibilityReport surfaceAccessibilityReport : reports) {
					if (!surfaceAccessibilityReport.isStored()) {
						if (!firstOne) {
							log.info("Appending report of protein " + surfaceAccessibilityProteinReport.getUniprotACC()
									+ " to file");
							firstOne = true;
						}
						out.println(surfaceAccessibilityReport.toString());
					}
				}
			}
			if (firstOne) {
				log.info("Report appended to file");
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (out != null) {
				out.close();
			}
			if (fileLocker != null) {
				try {
					fileLocker.release();
				} catch (IOException e) {
					//
				}
			}

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
			FileWriter fw = new FileWriter(surfaceAccessibilityFile, true);
			BufferedWriter bw = new BufferedWriter(fw);
			out = new PrintWriter(bw);
			// write headers
			out.println(SiteSurfaceAccessibilityReport.getToStringHeaders());
			List<String> keys = new ArrayList<String>();
			keys.addAll(reports.keySet());
			Collections.sort(keys);
			for (String key : keys) {
				final SurfaceAccessibilityProteinReport surfaceAccessibilityProteinReport = reports.get(key);
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
							report.setStored(true);
							addReport(report);
						}
					}
					loaded = true;
					log.info("Surface accesibilities numbers loaded from local file for " + reports.size() + " sites");
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
		String key = report.getReportKey();
		if (reports.containsKey(key)) {
			reports.get(key).addSurfaceAccesibilityReport(report);
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
					reports.put(key, proteinReport);
				}
			}
		}
	}

	private SiteSurfaceAccessibilityReport parseFromLine(String strLine) {
		return SiteSurfaceAccessibilityReport.getFromString(strLine);
	}

	public Map<String, SurfaceAccessibilityProteinReport> getSurfaceAccesibilityFromProteins(
			Collection<SurfaceProtein> proteins) {
		Map<String, SurfaceAccessibilityProteinReport> ret = new HashMap<String, SurfaceAccessibilityProteinReport>();

		ProgressCounter counter = new ProgressCounter(proteins.size(), ProgressPrintingType.PERCENTAGE_STEPS, 1);
		for (SurfaceProtein protein : proteins) {
			counter.increment();
			String printIfNecessary = counter.printIfNecessary();
			if (!"".equals(printIfNecessary)) {
				log.info(printIfNecessary);
			}

			final SurfaceAccessibilityProteinReport report = getProteinAccessibilityReportByProtein(protein);
			if (report != null) {
				ret.put(protein.getAcc(), report);
				// dumpToFile();
			}

		}
		return ret;
	}

	/**
	 * @param calculateIfNotPresent
	 *            the calculateIfNotPresent to set
	 */
	public void setCalculateIfNotPresent(boolean calculateIfNotPresent) {
		this.calculateIfNotPresent = calculateIfNotPresent;
	}

}
