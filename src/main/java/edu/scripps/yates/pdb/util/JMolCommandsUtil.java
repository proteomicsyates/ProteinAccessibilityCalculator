package edu.scripps.yates.pdb.util;

import java.text.DecimalFormat;
import java.util.Map;
import java.util.Set;

import com.compomics.util.experiment.biology.AminoAcid;

import edu.scripps.yates.pdb.JMolAtomReport;
import edu.scripps.yates.pdb.JMolScript;
import edu.scripps.yates.pdb.ProteinReport;
import edu.scripps.yates.pdb.model.Atom3D;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.surface.SurfaceReport;
import gnu.trove.set.hash.THashSet;

public class JMolCommandsUtil {
	private final static DecimalFormat df = new DecimalFormat("0.00");

	public static JMolScript getSelectAndAddLabels(ProteinReport proteinReport, String pdbID) {
		// delete all molecules not equals to
		JMolScript ret = new JMolScript();
		ret.addCommand("delete not protein");

		ret.addCommand("Color Background White");
		// color to light blue
		ret.addCommand("select all");
		ret.addCommand("color [217,255,255]");

		// add the labels for each site
		Map<String, Set<JMolAtomReport>> reportsByPDBID = proteinReport.getReportsByPDBID();
		Set<JMolAtomReport> reports = reportsByPDBID.get(pdbID);

		Set<String> chainIDs = getChainIDs(reports);
		ret.addCommand(getCommandForDeletingEverythingElse(chainIDs));
		for (JMolAtomReport report : reports) {

			// ret.addCommand("save state my_state");
			// ret.addCommand("delete not chain = " +
			// report.getAtom().getChainID());
			ret.addCommand("select atomno = " + report.getAtom().getAtomNumber());
			ret.addCommand("color black");
			ret.addCommand("color label black");
			if (report instanceof SurfaceReport) {
				SurfaceReport surfaceReport = (SurfaceReport) report;
				Double accessibility = surfaceReport.getAccessibility();
				String formated = "N/A";
				if (accessibility != null) {
					formated = df.format(accessibility);
				}
				ret.addCommand("label " + report.getAtom().getAa() + report.getPositionInUniprot() + ": " + formated);
			} else {
				throw new UnsupportedOperationException("report is not surface report. This has to be Implemented");
			}
			ret.addCommand("font label 15 serif bold");
			// ret.addCommand("restore state my_state");
		}

		ret.addCommand("select ");

		int i = 1;
		for (JMolAtomReport surfaceAccessibilityReport : reports) {
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

	private static Set<String> getChainIDs(Set<JMolAtomReport> reports) {
		Set<String> ret = new THashSet<String>();
		for (JMolAtomReport surfaceAccessibilityReport : reports) {
			ret.add(surfaceAccessibilityReport.getAtom().getChainID());
		}
		return ret;
	}

	public static JMolScript getDrawSurfaceJMolScriptByAtom(Atom3D atom) {
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

	public static JMolScript getSelectAtomScript(Atom3D atom) {
		JMolScript ret = new JMolScript();

		// select atom number
		ret.addCommand("select atomno = " + atom.getAtomNumber() + "");

		return ret;
	}

	public static JMolScript getCalculateSurfaceScript(Atom3D atom) {
		JMolScript ret = new JMolScript();
		// calculate isosurface of the selected atom
		ret.addCommand("select atomno=" + atom.getAtomNumber() + ";isoSurface saSurface;isosurface area set 0");

		return ret;
	}

	public static JMolScript getSelectChainJMolScriptByAtom(Atom3D atom, boolean removeOtherChains,
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

	public static JMolScript getCalculateDistanceScript(String chain1ID, String aa1, AtomType atomType1,
			String chain2ID, String aa2, AtomType atomType2) {
		JMolScript ret = new JMolScript();
		ret.addCommand("measure all ([" + convertToThreeLetterAA(aa1) + "]*" + chain1ID + "." + atomType1.name()
				+ ") ([" + convertToThreeLetterAA(aa2) + "]*" + chain2ID + "." + atomType2.name() + ")");
		return ret;
	}

	public static JMolScript getMeasureListScript() {
		JMolScript ret = new JMolScript();
		ret.addCommand("measure list");
		return ret;
	}

	public static String convertToThreeLetterAA(String aa) {
		return AminoAcid.getAminoAcid(aa.charAt(0)).threeLetterCode;
	}

	public static String convertToOneLetterAA(String threeLetterCodeAA) {
		for (String aminoacid : AminoAcid.getAminoAcids()) {
			if (AminoAcid.getAminoAcid(aminoacid.charAt(0)).threeLetterCode.equalsIgnoreCase(threeLetterCodeAA)) {
				return aminoacid;
			}
		}
		throw new IllegalArgumentException(
				threeLetterCodeAA + " is not recognized as a valid aminoacid three letter code.");
	}

	public static JMolScript getMeasureOffScript() {
		JMolScript ret = new JMolScript();
		ret.addCommand("measure off");
		return ret;
	}

}
