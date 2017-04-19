package edu.scripps.yates.pdb.read;

import java.io.IOException;

import org.junit.Test;
import org.springframework.core.io.ClassPathResource;

public class PDBParserTest {

	@Test
	public void parserTest() {
		String path;
		try {
			path = new ClassPathResource("pdb1z6x.ent").getFile().getAbsolutePath();
			PDBParser pdbParser = new PDBParser(path, "pdb1z6x");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
}
