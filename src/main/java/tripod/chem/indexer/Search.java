package tripod.chem.indexer;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;
import java.util.concurrent.*;

import chemaxon.formats.MolImporter;

public class Search {
    static final Logger logger = Logger.getLogger(Search.class.getName());

    public File index;
    public List<String> queries = new ArrayList<String>();
    public Map<String, String> filters = new HashMap<String, String>();
    public String search = "substructure";
    public double threshold = 0.8;
    public boolean list;
    
    public Search (String[] argv) {
        for (int i = 0; i < argv.length; ++i) {
            if (argv[i].charAt(0) == '-') {
                switch (argv[i].charAt(1)) {
                case 'h':
                    usage (System.err);
                    
                case 'f': {
                    String[] toks = argv[++i].split("=");
                    if (toks.length != 2) {
                        logger.warning("Invalid argument for -f option; "
                                       +"argument must be of the form "
                                       +"KEY=VALUE!");
                    }
                    else {
                        filters.put(toks[0], toks[1]);
                    }
                    break;
                }
                case 's': {
                    String arg = argv[++i];
                    if (arg.startsWith("sim") || arg.startsWith("tan")) {
                        search = "similarity";
                    }
                    else if (arg.startsWith("sub")) {
                        search = "substructure";
                    }
                    else {
                        logger.warning("Unknown value \""+arg
                                       +"\" for -s option; must be one of "
                                       +"{similarity, substructure}");
                    }
                    break;
                }
                case 't': {
                    String arg = argv[++i];
                    try {
                        threshold = Double.parseDouble(arg);
                        if (threshold < 0. || threshold > 1.) {
                            logger.warning("Threshold out of range: [0,1]");
                        }
                    }
                    catch (NumberFormatException ex) {
                        logger.warning("Bogus threshold value: "+arg);
                    }
                    break;
                }
                case 'l':
                    list = true;
                    break;
                    
                default:
                    logger.warning("Uknown option: "+argv[i]);
                }
            }
            else if (index == null) {
                index = new File (argv[i]);
                if (!index.exists()) {
                    logger.warning(argv[i]+": doesn't exist!");
                    usage ();
                }
                else if (!index.isDirectory()) {
                    logger.warning(argv[i]+": not a valid index directory!");
                    usage ();
                }
            }
            else {
                queries.add(argv[i]);
            }
        }

        if (index == null)
            usage ();
        
    }

    int process (StructureIndexer.ResultEnumeration result) throws Exception {
        int count = 0;
        while (result.hasMoreElements()) {
            StructureIndexer.Result r = result.nextElement();
            System.out.println
                ("++++ "+r.getId()+" matched = "+r.getSimilarity());
            ++count;
        }
        return count;
    }
    
    void similarity (StructureIndexer indexer, String q) throws Exception {
        long start = System.currentTimeMillis();
        int count = process (indexer.similarity(q, threshold));
        double ellapsed = (System.currentTimeMillis()-start)*1e-3;      
        //Thread.currentThread().sleep(5000);
        logger.info(q+": "+count+" matches found in "
                    +String.format("%1$.2fs", ellapsed));
    }

    void substructure (StructureIndexer indexer, String q)
        throws Exception {
        long start = System.currentTimeMillis();
        int count = process (indexer.substructure(q, 0, 3));
        double ellapsed = (System.currentTimeMillis()-start)*1e-3;      
        //Thread.currentThread().sleep(5000);
        logger.info(q+": "+count+" matches found in "
                    +String.format("%1$.2fs", ellapsed));
    }

    public void exec () throws Exception {
        StructureIndexer indexer = StructureIndexer.openReadOnly(index);
        try {
            if (list) {
                String[] fields = indexer.getFields();
                for (String f : fields) {
                    System.out.println("  \""+f+"\"");
                }
                System.out.println
                    ("## Index contains "+fields.length+" fields!");
            }
            
            for (String q : queries) {
                if (search.equalsIgnoreCase("similarity"))
                    similarity (indexer, q);
                else
                    substructure (indexer, q);
            }
        }
        finally {
            indexer.shutdown();
        }
    }

    static void usage () {
        usage (System.err);
    }
    
    static void usage (PrintStream ps) {
        ps.println
            ("Usage: Search [OPTIONS] INDEX SMILES/SMARTS...");
        ps.println("where OPTIONS can be one or more of the "
                   +"following:");
        ps.println("-h print this message");
        ps.println("-s {similarity|substructure}  search type "
                   +"(default: substructure)");
        ps.println("-t CUTOFF  specify Tanimoto cutoff (default: 0.8) for "
                   +"similarity search");
        ps.println("-f FIELD=VALUE  filter based on the provided filter-value "
                   +"pair. See -l option");
        ps.println(" for a list of valid FIELDs.");
        ps.println("-l list available filter FIELDs");
        
        System.exit(1);
    }
    
    public static void main (String[] argv) throws Exception {
        Search search = new Search (argv);
        search.exec();
    }
}
