package edu.scripps.yates.pdb.read;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class PDBParserManager {
	// static variables
	private final static Map<String, PDBParser> staticParsersMap = new HashMap<String, PDBParser>();
	private static final Set<String> pdbsNotRetrieved = new HashSet<String>();;

	// non static variables
	private final PDBFileManager pdbFileManager;
	private boolean keepManagersInMemory = false;
	private final Map<String, PDBParser> parsersMap = new HashMap<String, PDBParser>();

	public PDBParserManager(File parentFolder) {
		pdbFileManager = PDBFileManager.getInstance(parentFolder);
	}

	public void clearParsers() {
		parsersMap.clear();
	}

	public PDBParser getPDBParserByPDBID(String pdbID) {
		// if it is in the notretrieved set, return null right away
		if (pdbsNotRetrieved.contains(pdbID)) {
			return null;
		}
		if (!staticParsersMap.containsKey(pdbID) && !parsersMap.containsKey(pdbID)) {
			try {
				// get from local file system
				File pdbFile = pdbFileManager.getPDBFile(pdbID);
				if (pdbFile == null) {
					// get from server
					pdbFile = PDBFileRetriever.getPDBFile(pdbID);
					if (pdbFile != null) {
						// save to local file system
						pdbFile = pdbFileManager.savePDBFile(pdbFile, pdbID);
					}
					if (pdbFile == null) {
						pdbsNotRetrieved.add(pdbID);
						return null;
					}
				}

				PDBParser parser = new PDBParser(pdbFile.getAbsolutePath(), pdbID);
				parsersMap.put(pdbID, parser);
				if (keepManagersInMemory) {
					staticParsersMap.put(pdbID, parser);
				}
				return parser;

			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		if (parsersMap.containsKey(pdbID)) {
			return parsersMap.get(pdbID);
		}
		return staticParsersMap.get(pdbID);

	}

	/**
	 * @return the pdbFileManager
	 */
	public PDBFileManager getPdbFileManager() {
		return pdbFileManager;
	}

	/**
	 * @return the keepManagersInMemory
	 */
	public boolean isKeepManagersInMemory() {
		return keepManagersInMemory;
	}

	/**
	 * @param keepManagersInMemory
	 *            the keepManagersInMemory to set
	 */
	public void setKeepManagersInMemory(boolean keepManagersInMemory) {
		this.keepManagersInMemory = keepManagersInMemory;
	}
}
