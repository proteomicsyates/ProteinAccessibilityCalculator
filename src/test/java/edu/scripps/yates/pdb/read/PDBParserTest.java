package edu.scripps.yates.pdb.read;

import static org.junit.Assert.fail;

import java.io.IOException;
import java.util.List;

import org.junit.Test;
import org.springframework.core.io.ClassPathResource;

import edu.scripps.yates.pdb.model.Atom3D;
import edu.scripps.yates.pdb.model.AtomType;
import edu.scripps.yates.pdb.model.DBRef;

public class PDBParserTest {

	@Test
	public void parserTest() {
		String path;
		try {
			path = new ClassPathResource("pdb1z6x.ent").getFile().getAbsolutePath();
			PDBParser pdbParser = new PDBParser(path, "pdb1z6x", false);
		} catch (IOException e) {
			e.printStackTrace();
			fail();
		}

	}

	@Test
	public void parserTest2() {
		String path;
		try {
			path = new ClassPathResource("4wnc.pdb").getFile().getAbsolutePath();
			PDBParser pdbParser = new PDBParser(path, "4wnc", false);

			List<String> chainIDs = pdbParser.getChainIDs();
			int i = 1;
			for (String chainID : chainIDs) {
				List<Atom3D> atoms = pdbParser.getAtoms(chainID, "K", AtomType.NZ);
				for (Atom3D atom : atoms) {
					System.out.println(i++ + "\t" + chainID + "\t" + atom);
				}
			}
			List<DBRef> dbRefs = pdbParser.getDBRefs();
			for (DBRef dbRef : dbRefs) {
				System.out.println(dbRef);
			}
		} catch (IOException e) {
			e.printStackTrace();
			fail();
		} catch (NotValidPDBException e) {
			e.printStackTrace();
			fail();
		}

	}

}
