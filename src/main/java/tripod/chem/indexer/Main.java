package tripod.chem.indexer;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.concurrent.*;

import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;

public class Main {
    static final Logger logger = Logger.getLogger(Main.class.getName());
    
    public static void main (String[] argv) throws Exception {
        if (argv.length < 2) {
            System.err.println("Usage: Main INDEXDIR FILES...");
            System.exit(1);
        }

        File dir = new File (argv[0]);
        if (!dir.exists())
            dir.mkdirs();

        StructureIndexer indexer = new StructureIndexer (dir);
        try {
            logger.info("Indexing structures...");
            
            long start = System.currentTimeMillis(), total = 0;
            for (int i = 1; i < argv.length; ++i) {
                File file = new File (argv[i]);
                MolImporter mi = new MolImporter (argv[i]);
                int count = 0;
                for (Molecule m = new Molecule (); mi.read(m); ++count) {
                    indexer.add(file.getName(), m);
                }
                logger.info(file.getName()+": "+count);
                total += count;
            }
            logger.info("Indexing time "+String.format
                        ("%1$.2fs",1e-3*(System.currentTimeMillis()-start))
                        +" for "+total+" structures!");
        }
        finally {
            indexer.shutdown();
        }
    }
}
