package tripod.chem.indexer;
import java.io.*;
import java.util.*;

import chemaxon.util.MolHandler;
import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class StructureIndexerTest extends TestCase {
    public StructureIndexerTest () {
        super ("StructureIndexerTest");
    }
    
    public static Test suite () {
        return new TestSuite (StructureIndexerTest.class);
    }
    
    public void test1 () {
        assertTrue (true);
    }
}
