package edu.scripps.yates.pdb.surface;

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
public class SurfaceAccessibilityProteinReport {
	private final String uniprotACC;
	private final String uniprotProteinSequence;
	private final TIntObjectHashMap<Set<SiteSurfaceAccessibilityReport>> accessibilitiesByPositionInUniprotSeq = new TIntObjectHashMap<Set<SiteSurfaceAccessibilityReport>>();
	private final List<SiteSurfaceAccessibilityReport> reports = new ArrayList<SiteSurfaceAccessibilityReport>();

	private static final Logger log = Logger.getLogger(SurfaceAccessibilityProteinReport.class);

	public SurfaceAccessibilityProteinReport(String uniprotACC, String uniprotProteinSeq) {
		this.uniprotACC = uniprotACC;
		uniprotProteinSequence = uniprotProteinSeq;
	}

	/**
	 * @return the uniprotACC
	 */
	public String getUniprotACC() {
		return uniprotACC;
	}

	public Set<SiteSurfaceAccessibilityReport> getSurfaceAccessibilityReportsBySite(int positionInUniprotSeq) {
		return accessibilitiesByPositionInUniprotSeq.get(positionInUniprotSeq);
	}

	public void addSurfaceAccesibilityReport(SiteSurfaceAccessibilityReport report) {

		if (report.getUniprotACC() != null && !report.getUniprotACC().equals(uniprotACC)) {
			log.warn("Reports from different proteins cannot be merged");
			return;
		}

		int positionInUniprotSeq = report.getPositionInUniprot();
		if (positionInUniprotSeq > -1) {
			if (accessibilitiesByPositionInUniprotSeq.containsKey(positionInUniprotSeq)) {
				accessibilitiesByPositionInUniprotSeq.get(positionInUniprotSeq).add(report);
			} else {
				Set<SiteSurfaceAccessibilityReport> set = new THashSet<SiteSurfaceAccessibilityReport>();
				set.add(report);
				accessibilitiesByPositionInUniprotSeq.put(positionInUniprotSeq, set);
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
	public TIntObjectHashMap<Set<SiteSurfaceAccessibilityReport>> getAccessibilitiesByPositionInUniprotSeq() {
		return accessibilitiesByPositionInUniprotSeq;
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
		final int[] keySet = accessibilitiesByPositionInUniprotSeq.keys();
		for (int position : keySet) {
			final String key = uniprotACC + "-" + position;
			if (!list.contains(key)) {
				list.add(key);
			}
		}
		return list;
	}

	public List<SiteSurfaceAccessibilityReport> getReports() {
		return reports;
	}

	public Map<String, Set<SiteSurfaceAccessibilityReport>> getReportsByPDBID() {
		final List<SiteSurfaceAccessibilityReport> reports = getReports();
		Map<String, Set<SiteSurfaceAccessibilityReport>> ret = new THashMap<String, Set<SiteSurfaceAccessibilityReport>>();
		for (SiteSurfaceAccessibilityReport report : reports) {
			if (ret.containsKey(report.getPdbID())) {
				ret.get(report.getPdbID()).add(report);
			} else {
				Set<SiteSurfaceAccessibilityReport> set = new THashSet<SiteSurfaceAccessibilityReport>();
				set.add(report);
				ret.put(report.getPdbID(), set);
			}
		}
		return ret;
	}

	public boolean containsReportsForPosition(int positionInUniprotProtein) {
		return accessibilitiesByPositionInUniprotSeq.containsKey(positionInUniprotProtein);
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

		for (int position : accessibilitiesByPositionInUniprotSeq.keys()) {
			positions.add(position);
		}

		Collections.sort(positions);
		for (int uniprotPosition : positions) {
			final Set<SiteSurfaceAccessibilityReport> accesibilitySet = accessibilitiesByPositionInUniprotSeq
					.get(uniprotPosition);
			for (SiteSurfaceAccessibilityReport siteSurfaceAccessibilityReport : accesibilitySet) {
				sb.append(siteSurfaceAccessibilityReport.toString());
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
		for (SiteSurfaceAccessibilityReport report : reports) {
			arrayBuilder.add(report.toJson());
		}
		return builder.build();
	}
}
