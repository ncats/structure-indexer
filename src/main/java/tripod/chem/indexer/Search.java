package tripod.chem.indexer;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.concurrent.*;

import chemaxon.formats.MolImporter;

public class Search {
    static final Logger logger = Logger.getLogger(Search.class.getName());
    
    public static void main (String[] argv) throws Exception {
        if (argv.length < 2) {
            System.err.println("Usage: Search INDEX QUERIES...");
            System.exit(1);
        }

        File dir = new File (argv[0]);
        if (!dir.exists())
            throw new IllegalArgumentException
                ("INDEX directory "+dir+" doesn't exist!");

        StructureIndexer indexer = new StructureIndexer (dir);
        try {
            for (int i = 1; i < argv.length; ++i) {
                long start = System.currentTimeMillis();
                StructureIndexer.ResultEnumeration result
                    //= indexer.substructure(argv[i], -1, 3);
                    = indexer.similarity(argv[i], 0.9);
                
                double ellapsed = (System.currentTimeMillis()-start)*1e-3;
                if (result != null) {
                    int count = 0;
                    while (result.hasMoreElements()) {
                        StructureIndexer.Result r = result.nextElement();
                        System.out.println(r.getMol().getName()+"\t"+
                                           r.getMol().getProperty("TANIMOTO"));
                        /*
                          System.out.println(r.getMol().toFormat("smiles:q")
                          +"\t"+r.getMol().getName());
                        */
                        ++count;
                    }
                    logger.info(argv[i]+": "+count+" matches found in "
                                +String.format("%1$.2fs", ellapsed));
                }
            }
        }
        finally {
            indexer.shutdown();
        }
    }
}
