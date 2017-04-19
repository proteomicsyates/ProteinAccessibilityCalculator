package edu.scripps.yates.pdb.util;

import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.Set;

import edu.scripps.yates.pdb.model.Atom;
import edu.scripps.yates.pdb.surface.JMolScript;
import edu.scripps.yates.pdb.surface.SiteSurfaceAccessibilityReport;
import edu.scripps.yates.pdb.surface.SurfaceAccessibilityProteinReport;

public class JMolCommandsUtil {
	private final static DecimalFormat df = new DecimalFormat("0.00");

	public static JMolScript getSelectAndAddLabels(SurfaceAccessibilityProteinReport proteinReport, String pdbID) {
		// delete all molecules not equals to
		JMolScript ret = new JMolScript();
		ret.addCommand("delete not protein");

		ret.addCommand("Color Background White");
		// color to light blue
		ret.addCommand("select all");
		ret.addCommand("color [217,255,255]");

		// add the labels for each site
		Set<SiteSurfaceAccessibilityReport> reports = proteinReport.getReportsByPDBID().get(pdbID);

		Set<String> chainIDs = getChainIDs(reports);
		ret.addCommand(getCommandForDeletingEverythingElse(chainIDs));
		for (SiteSurfaceAccessibilityReport report : reports) {

			// ret.addCommand("save state my_state");
			// ret.addCommand("delete not chain = " +
			// report.getAtom().getChainID());
			ret.addCommand("select atomno = " + report.getAtom().getAtomNumber());
			ret.addCommand("color black");
			ret.addCommand("color label black");
			ret.addCommand("label " + report.getAa() + report.getPositionInUniprot() + ": "
					+ df.format(report.getAccesibility()));
			ret.addCommand("font label 15 serif bold");
			// ret.addCommand("restore state my_state");
		}

		ret.addCommand("select ");

		int i = 1;
		for (SiteSurfaceAccessibilityReport surfaceAccessibilityReport : reports) {
			ret.appendToLastCommand("atomno = " + surfaceAccessibilityReport.getAtom().getAtomNumber());
			if (i < reports.size()) {
				ret.appendToLastCommand(",");
			}
			i++;
		}
		// show isosurface of the selected atoms
		ret.addCommand("isoSurface saSurface");
		return ret;

	}

	private static String getCommandForDeletingEverythingElse(Set<String> chainIDs) {
		StringBuilder sb = new StringBuilder();
		sb.append("delete not (");

		int numChain = 1;
		for (String chainID : chainIDs) {
			if (numChain > 1 && numChain < chainIDs.size()) {
				sb.append(" OR ");
			}
			sb.append("chain = " + chainID);
			numChain++;
		}
		sb.append(")");
		return sb.toString();
	}

	private static Set<String> getChainIDs(Set<SiteSurfaceAccessibilityReport> reports) {
		Set<String> ret = new HashSet<String>();
		for (SiteSurfaceAccessibilityReport surfaceAccessibilityReport : reports) {
			ret.add(surfaceAccessibilityReport.getAtom().getChainID());
		}
		return ret;
	}

	public static JMolScript getDrawSurfaceJMolScriptByAtom(Atom atom) {
		// delete all molecules not equals to
		JMolScript sb = new JMolScript();
		// delete all molecules different than a protein
		sb.addCommand("delete not protein");
		// delete all molecules different than our chain
		sb.addCommand("delete not chain = " + atom.getChainID());
		// select atom number
		sb.addCommand("select atomno = " + atom.getAtomNumber());
		// show isosurfaces
		sb.addCommand("isoSurface saSurface");

		return sb;
	}

	public static JMolScript getSelectAtomScript(Atom atom) {
		JMolScript ret = new JMolScript();

		// select atom number
		ret.addCommand("select atomno = " + atom.getAtomNumber() + "");

		return ret;
	}

	public static JMolScript getCalculateSurfaceScript(Atom atom) {
		JMolScript ret = new JMolScript();
		// calculate isosurface of the selected atom
		ret.addCommand("select atomno=" + atom.getAtomNumber() + ";isoSurface saSurface;isosurface area set 0");

		return ret;
	}

	public static JMolScript getSelectChainJMolScriptByAtom(Atom atom, boolean removeOtherChains,
			boolean removeOtherMolecules) {
		// delete all molecules not equals to
		JMolScript ret = new JMolScript();
		// delete all molecules different than a protein
		if (removeOtherMolecules) {
			ret.addCommand("delete not protein");
		}
		if (removeOtherChains) {
			// delete all molecules different than our chain
			ret.addCommand("delete not chain = " + atom.getChainID() + "");
		}
		return ret;
	}

}
