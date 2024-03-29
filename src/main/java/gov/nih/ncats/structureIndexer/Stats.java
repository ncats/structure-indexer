package gov.nih.ncats.structureIndexer;

import java.io.*;
import java.util.logging.Logger;

public class Stats {
    static final Logger logger = Logger.getLogger(Stats.class.getName());
    
    public static void main (String[] argv) throws Exception {
        if (argv.length == 0) {
            System.err.println("Usage: Stats INDEXDIR");
            System.exit(1);
        }
        File dir = new File (argv[0]);
        if (!dir.exists())
            throw new IllegalArgumentException ("Bogus index directory!");
        
        StructureIndexer indexer = StructureIndexer.openReadOnly(dir);
        try {
            indexer.stats(System.out);
        }
        finally {
            indexer.shutdown();
        }
    }
}
