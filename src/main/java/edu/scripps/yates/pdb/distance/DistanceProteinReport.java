package edu.scripps.yates.pdb.distance;

import edu.scripps.yates.pdb.ProteinReport;

/**
 * Stores all SurfaceAccesibilityReport from a protein, stored by position in
 * its sequence
 */
public class DistanceProteinReport extends ProteinReport<DistanceReport> {

	public DistanceProteinReport(String uniprotACC, String uniprotProteinSeq) {
		super(uniprotACC, uniprotProteinSeq);
	}

}
