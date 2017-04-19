package edu.scripps.yates.pdb.read;

import java.io.File;

import org.junit.Assert;
import org.junit.Test;

public class PDBFileRetrieverTest {
	@Test
	public void PDBFileRetrieverTest1() {
		String pdbID = "5EDK";
		final File pdbFile = PDBFileRetriever.getPDBFile(pdbID);
		Assert.assertTrue(pdbFile.exists());
		System.out.println(pdbFile.length());

	}

	@Test
	public void PDBFileRetrieverTest2() {
		String pdbID = "5EDK";
		final File pdbFile = PDBFileRetriever.getPDBGZipFile(pdbID);
		Assert.assertTrue(pdbFile.exists());
		System.out.println(pdbFile.length());
	}
}
