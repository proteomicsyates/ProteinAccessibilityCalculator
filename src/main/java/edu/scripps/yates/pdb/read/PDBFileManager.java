package edu.scripps.yates.pdb.read;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.apache.log4j.Logger;

import gnu.trove.map.hash.THashMap;

public class PDBFileManager {
	private static final Logger log = Logger.getLogger(PDBFileManager.class);
	private static Map<File, PDBFileManager> instances = new THashMap<File, PDBFileManager>();
	private final File parentPath;

	private PDBFileManager(File parentPath) {
		if (!parentPath.isDirectory() && parentPath.isFile()) {
			throw new IllegalArgumentException(parentPath.getAbsolutePath() + " is not a folder");
		}
		if (!parentPath.exists()) {
			if (!parentPath.mkdirs()) {
				throw new IllegalArgumentException("Error creating folder " + parentPath.getAbsolutePath());
			}
		}
		this.parentPath = parentPath;
	}

	public static PDBFileManager getInstance(File parentPath) {
		if (!instances.containsKey(parentPath)) {
			final PDBFileManager instance = new PDBFileManager(parentPath);
			instances.put(parentPath, instance);
		}
		return instances.get(parentPath);
	}

	public File getPDBFile(String pdbID) throws IOException {
		final File finalFile = getFile(pdbID);
		if (finalFile.exists()) {
			return finalFile;
		} else {
			return null;
		}
	}

	public File savePDBFile(File pdbFile, String pdbID) throws IOException {
		if (pdbFile != null && pdbFile.exists() && pdbFile.isFile()) {
			final File finalFile = getFile(pdbID);
			pdbFile.renameTo(finalFile);
			log.info("PDB file saved at: " + finalFile);
			return finalFile;
		} else {
			log.warn("Error during saving PDB file " + pdbID);
		}
		return null;
	}

	private File getFile(String pdbID) {
		return new File(parentPath.getAbsolutePath() + File.separator + pdbID + ".pdb");
	}

	/**
	 * @return the parentPath
	 */
	public File getParentPath() {
		return parentPath;
	}
}
