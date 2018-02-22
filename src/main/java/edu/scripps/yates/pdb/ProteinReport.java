package edu.scripps.yates.pdb;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.json.Json;
import javax.json.JsonArrayBuilder;
import javax.json.JsonObject;
import javax.json.JsonObjectBuilder;

import org.apache.log4j.Logger;

import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.THashSet;

/**
 * Stores all SurfaceAccesibilityReport from a protein, stored by position in
 * its sequence
 */
public abstract class ProteinReport<T extends JMolAtomReport> {
	private final String uniprotACC;
	private final String uniprotProteinSequence;
	private final TIntObjectHashMap<Set<T>> reportsByPositionInUniprotSeq = new TIntObjectHashMap<Set<T>>();
	private final List<T> reports = new ArrayList<T>();

	private static final Logger log = Logger.getLogger(ProteinReport.class);

	public ProteinReport(String uniprotACC, String uniprotProteinSeq) {
		this.uniprotACC = uniprotACC;
		uniprotProteinSequence = uniprotProteinSeq;
	}

	/**
	 * @return the uniprotACC
	 */
	public String getUniprotACC() {
		return uniprotACC;
	}

	public Set<T> getReportsBySite(int positionInUniprotSeq) {
		return reportsByPositionInUniprotSeq.get(positionInUniprotSeq);
	}

	public void addReport(T report) {

		if (report.getUniprotACC() != null && !report.getUniprotACC().equals(uniprotACC)) {
			log.warn("Reports from different proteins cannot be merged");
			return;
		}

		int positionInUniprotSeq = report.getPositionInUniprot();
		if (positionInUniprotSeq > -1) {
			if (reportsByPositionInUniprotSeq.containsKey(positionInUniprotSeq)) {
				reportsByPositionInUniprotSeq.get(positionInUniprotSeq).add(report);
			} else {
				Set<T> set = new THashSet<T>();
				set.add(report);
				reportsByPositionInUniprotSeq.put(positionInUniprotSeq, set);
			}
		}
		reports.add(report);
	}

	public boolean isEmpty() {
		return reports.isEmpty();
	}

	/**
	 * @return the accesibilitiesByPositionInUniprotSeq
	 */
	public TIntObjectHashMap<Set<T>> getReportsByPositionInUniprotSeq() {
		return reportsByPositionInUniprotSeq;
	}

	/**
	 * @return the uniprotProteinSequence
	 */
	public String getUniprotProteinSequence() {
		return uniprotProteinSequence;
	}

	/**
	 * Get a list of keys formed by individual 'ACC-position'
	 *
	 * @return
	 */
	public List<String> getUniquePositionsInProteinKeys() {
		List<String> list = new ArrayList<String>();
		final int[] keySet = reportsByPositionInUniprotSeq.keys();
		for (int position : keySet) {
			final String key = uniprotACC + "-" + position;
			if (!list.contains(key)) {
				list.add(key);
			}
		}
		return list;
	}

	public List<T> getReports() {
		return reports;
	}

	public Map<String, Set<T>> getReportsByPDBID() {
		final List<T> reports = getReports();
		Map<String, Set<T>> ret = new THashMap<String, Set<T>>();
		for (T report : reports) {
			if (ret.containsKey(report.getPdbID())) {
				ret.get(report.getPdbID()).add(report);
			} else {
				Set<T> set = new THashSet<T>();
				set.add(report);
				ret.put(report.getPdbID(), set);
			}
		}
		return ret;
	}

	public boolean containsReportsForPosition(int positionInUniprotProtein) {
		return reportsByPositionInUniprotSeq.containsKey(positionInUniprotProtein);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		List<Integer> positions = new ArrayList<Integer>();

		for (int position : reportsByPositionInUniprotSeq.keys()) {
			positions.add(position);
		}

		Collections.sort(positions);
		for (int uniprotPosition : positions) {
			final Set<T> accesibilitySet = reportsByPositionInUniprotSeq.get(uniprotPosition);
			for (T report : accesibilitySet) {
				sb.append(report.toString());
				sb.append("\n");
			}
		}

		return sb.toString();
	}

	public JsonObject toJson() {
		JsonObjectBuilder builder = Json.createObjectBuilder().add("uniprotACC", uniprotACC)
				.add("uniprotProteinSequence", uniprotProteinSequence);
		JsonArrayBuilder arrayBuilder = Json.createArrayBuilder();
		builder.add("reports", arrayBuilder);
		for (T report : reports) {
			arrayBuilder.add(report.toJson());
		}
		return builder.build();
	}
}
