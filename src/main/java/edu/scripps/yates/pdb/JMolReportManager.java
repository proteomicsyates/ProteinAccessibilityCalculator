package edu.scripps.yates.pdb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.pdb.model.Protein;
import edu.scripps.yates.pdb.surface.SurfaceCalculator;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;

public abstract class JMolReportManager<R extends ProteinReport<T>, T extends JMolAtomReport> {
	protected final static Logger log = Logger.getLogger(JMolReportManager.class);
	protected final Calculator<R, T> calculator;
	private final File file;
	private boolean loaded;
	private final Map<String, R> reportsByKey = new THashMap<String, R>();

	public JMolReportManager(Calculator<R, T> calculator) {
		this.calculator = calculator;
		this.calculator.setManager(this);
		file = new File(calculator.getPdbParserManager().getPdbFileManager().getParentPath().getAbsolutePath()
				+ File.separator + getFileName());
		log.info("Reports file located at: " + file.getAbsolutePath());
	}

	/**
	 * Gets the name of the file
	 * 
	 * @return
	 */
	public abstract String getFileName();

	public Map<String, R> getProteinReportsByProteins(Collection<Protein> proteins, Entry entry) {
		Map<String, R> ret = new THashMap<String, R>();
		for (Protein protein : proteins) {
			final R proteinReportByProtein = getProteinReportByProtein(protein.getAcc(), null, entry);
			if (proteinReportByProtein != null) {
				ret.put(protein.getAcc(), proteinReportByProtein);
			}
		}
		return ret;
	}

	public R getProteinReportByProtein(Protein protein, String pdbID, Entry entry) {
		return getProteinReportByProtein(protein.getAcc(), pdbID, entry);
	}

	public R getProteinReportByProtein(Protein protein, String pdbID) {
		return getProteinReportByProtein(protein.getAcc(), pdbID, null);
	}

	public R getProteinReportByProtein(String proteinAcc, String pdbID) {
		return getProteinReportByProtein(proteinAcc, pdbID, null);
	}

	public R getProteinReportByProtein(String proteinAcc, String pdbID, Entry entry) {
		// load data from file if not yet
		loadReportsFromFile();
		// look into the map
		log.debug("Getting reports for protein: " + proteinAcc);

		final R proteinReport = calculator.getReportFromProtein(proteinAcc, pdbID, entry);
		if (proteinReport != null) {

			TIntObjectHashMap<Set<T>> reportsByPositionInUniprotSeq = proteinReport.getReportsByPositionInUniprotSeq();
			for (Set<T> reports : reportsByPositionInUniprotSeq.valueCollection()) {
				for (T report : reports) {
					addReport(report);
				}
			}
			if (!proteinReport.getReportsByPositionInUniprotSeq().isEmpty()) {
				appendReportToFile(proteinReport);
			}
			return proteinReport;
		}
		return null;

	}

	public R getReportByKey(String reportKey) {
		// load data from file if not yet
		loadReportsFromFile();
		// look into the map

		return reportsByKey.get(reportKey);

	}

	private void appendReportToFile(R proteinReport) {
		PrintWriter out = null;
		try {
			boolean writeHeader = false;
			if (!file.exists() || file.length() == 0l) {
				// write the header = true
				writeHeader = true;
			}
			FileWriter fw = new FileWriter(file, true);
			BufferedWriter bw = new BufferedWriter(fw);
			out = new PrintWriter(bw);

			boolean firstOne = true;
			final TIntObjectHashMap<Set<T>> positions = proteinReport.getReportsByPositionInUniprotSeq();
			List<Integer> sortedPositions = new ArrayList<Integer>();
			for (int position : positions.keys()) {
				sortedPositions.add(position);
			}
			Collections.sort(sortedPositions);
			for (Integer position : sortedPositions) {
				final Set<T> reports = positions.get(position);
				for (T report : reports) {
					if (!report.isStored()) {
						if (firstOne) {
							log.info("Appending report of protein " + proteinReport.getUniprotACC() + " to file");
							firstOne = false;
							if (writeHeader) {
								out.println(report.getToStringHeaders());
							}
						}
						out.println(report.toString());
					}
				}
			}
			if (!firstOne) {
				log.info("Report appended to file");
			}
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
			log.info("Writting into file reports for " + reportsByKey.size() + " proteins");
			FileWriter fw = new FileWriter(file, true);
			BufferedWriter bw = new BufferedWriter(fw);
			out = new PrintWriter(bw);

			List<String> keys = new ArrayList<String>();
			keys.addAll(reportsByKey.keySet());
			Collections.sort(keys);
			boolean firstOne = true;
			for (String key : keys) {
				final R proteinReport = reportsByKey.get(key);
				final TIntObjectHashMap<Set<T>> positions = proteinReport.getReportsByPositionInUniprotSeq();
				List<Integer> sortedPositions = new ArrayList<Integer>();
				for (int position : positions.keys()) {
					sortedPositions.add(position);
				}

				Collections.sort(sortedPositions);
				for (Integer position : sortedPositions) {
					final Set<T> reports = positions.get(position);
					for (JMolAtomReport surfaceAccessibilityReport : reports) {
						if (firstOne) {
							log.info("Overriding report of protein " + proteinReport.getUniprotACC() + " to file");
							// write header for the first one
							out.println(surfaceAccessibilityReport.getToStringHeaders());
							firstOne = false;
						}
						out.println(surfaceAccessibilityReport.toString());
					}
				}
			}
			log.info("Reports writed at: " + file.getAbsolutePath());
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
				if (file.exists()) {
					// first read the file to get only the uniprot accessions
					// then call to annotate those proteins in uniprot
					// this will avoid the annotation of proteins one by one
					// when a new version of uniprot is released
					annotateAllProteinsInReports();
					FileInputStream fstream = new FileInputStream(file);
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
						T report = parseFromLine(strLine);
						if (report != null) {
							report.setStored(true);
							addReport(report);
						}
					}
					loaded = true;
					log.info(" accesibilities numbers loaded from local file for " + reportsByKey.size() + " sites");
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

	private void annotateAllProteinsInReports() {
		try {
			UniprotProteinLocalRetriever uplr = calculator.getUplr();
			if (uplr != null) {
				log.info(
						"Getting Uniprot annotations for all the proteins in the stored reports, in order to get the protein sequences");
				if (file.exists()) {
					Set<String> accs = Files.lines(Paths.get(file.getAbsolutePath())).map(l -> l.split("\t")[4])
							.collect(Collectors.toSet());

					uplr.getAnnotatedProteins(null, accs);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void addReport(T report) {
		String key = report.getReportKey();
		if (reportsByKey.containsKey(key)) {
			reportsByKey.get(key).addReport(report);
		} else {
			final Map<String, Entry> annotatedProtein = calculator.getUplr().getAnnotatedProtein(null,
					report.getUniprotACC());
			if (annotatedProtein.containsKey(report.getUniprotACC())) {
				final Entry entry = annotatedProtein.get(report.getUniprotACC());
				final String proteinSequence = SurfaceCalculator.getUniprotProteinSequence(entry);
				if (proteinSequence != null) {
					R proteinReport = createProteinReport(report.getUniprotACC(), proteinSequence);
					proteinReport.addReport(report);
					reportsByKey.put(key, proteinReport);
				}
			}
		}
	}

	public abstract R createProteinReport(String proteinacc, String proteinSequence);

	public abstract T parseFromLine(String strLine);

	public Map<String, R> getReportsFromProteins(Collection<Protein> proteins) {
		List<Protein> list = new ArrayList<Protein>();
		list.addAll(proteins);
		return getReportsFromProteins(list, null);
	}

	public Map<String, R> getReportsFromProteins(List<Protein> proteins, List<Entry> entries) {
		Map<String, R> ret = new THashMap<String, R>();

		ProgressCounter counter = new ProgressCounter(proteins.size(), ProgressPrintingType.PERCENTAGE_STEPS, 1);
		int index = 0;
		for (Protein protein : proteins) {
			counter.increment();
			String printIfNecessary = counter.printIfNecessary();
			if (!"".equals(printIfNecessary)) {
				log.info(printIfNecessary);
			}

			final R report = getProteinReportByProtein(protein, null, entries.get(index));
			if (report != null) {
				ret.put(protein.getAcc(), report);
				// dumpToFile();
			}
			index++;
		}
		return ret;
	}

	/**
	 * @param calculateIfNotPresent
	 *            the calculateIfNotPresent to set
	 */
	public void setCalculateIfNotPresent(boolean calculateIfNotPresent) {
	}
}
