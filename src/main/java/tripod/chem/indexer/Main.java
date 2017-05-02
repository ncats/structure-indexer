package tripod.chem.indexer;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;
import java.util.stream.Stream;
import java.util.logging.Level;
import java.util.concurrent.*;

import gov.nih.ncats.chemkit.api.Chemical;
import gov.nih.ncats.chemkit.api.ChemicalReader;
import gov.nih.ncats.chemkit.api.ChemicalReaderFactory;
import gov.nih.ncats.chemkit.api.util.stream.ThrowingStream;

public class Main {
    static final Logger logger = Logger.getLogger(Main.class.getName());

    File index;
    String idField;
    boolean listSource;
    List<File> files = new ArrayList<File>();
    
    public Main (String[] argv) throws IOException {
        for (int i = 0; i < argv.length; ++i) {
            if (argv[i].charAt(0) == '-') {
                switch (argv[i].charAt(1)) {
                case 'h':
                    usage (System.err);
                    
                case 'i':
                    if (argv[i].length() > 2)
                        idField = argv[i].substring(2);
                    else
                        idField = argv[++i];
                    logger.info("Id Field: \""+idField+"\"");
                    break;
                    
                case 'l':
                    listSource = true;
                    break;
                    
                default:
                    logger.warning("Unknown option: "+argv[i]);
                }
            }
            else {
                File f = new File (argv[i]);
                if (index == null) {
                    if (!f.exists())
                        f.mkdirs();
                    index = f;
                }
                else if (f.exists()) {
                    files.add(f);
                }
                else {
                    logger.warning(argv[i]+": Invalid file!");
                }
            }
        }
        
        if (index == null) {
            logger.warning("No INDEX directory specified!");
            usage ();
        }
    }

    public void exec () throws Exception {
        StructureIndexer indexer = StructureIndexer.open(index);
        try {
            if (listSource) {
                Map<String, Integer> sources = indexer.getSources();
                logger.info(sources.size()+" source(s) indexed:");
                for (Map.Entry<String, Integer> me : sources.entrySet()) {
                    System.out.println(me.getKey()+"\t"+me.getValue());
                }
                System.out.println("** Total: "+indexer.size());
            }
            
            if (!files.isEmpty())
                logger.info("Adding structures to index "+index+"...");
            
            long start = System.currentTimeMillis(), total = 0;
            for (File f : files) {
            	 int count = 0;
            	 String source = getSourceNameFrom(f);
            	 try(ThrowingStream<Chemical> stream = ChemicalReaderFactory.newReader(f).stream()){
            		 stream.throwingForEach(chem ->{
	                    String id = chem.getName();
	                    if (id == null) {
	                    	id = String.format("%1$010d", count+1);
	                    }
	                    
	                    indexer.add(source, id, chem);
	                });
            	 }
                logger.info(f.getName()+": "+count+"/"+indexer.size());
                total += count;
            }
            
            if (total > 0) {
                logger.info("Indexing time "+String.format
                            ("%1$.2fs",1e-3*(System.currentTimeMillis()-start))
                            +" for "+total+" structures!");
            }
        }
        finally {
            indexer.shutdown();
        }
    }

	private String getSourceNameFrom(File f) {
		String source;
		 int pos = f.getName().lastIndexOf('.');
		 if (pos > 0) {
		     source = f.getName().substring(0, pos);
		 }else{
			 source = f.getName();
		 }
		return source;
	}
    
    static void usage () {
        usage (System.err);
    }
    
    static void usage (PrintStream ps) {
        ps.println("Usage: Main [OPTIONS] INDEX FILES...");
        ps.println("where INDEX is the index directory and OPTIONS can");
        ps.println("be one or more of the following:");
        ps.println("-h print this message");
        ps.println("-i FIELD  specify the field name to extract ID; if "
                   +"not specified,");
        ps.println("   an autoincrement value is used");
        ps.println("-l print all the source filenames that have been indexed");
        System.exit(1);
    }
    
    public static void main (String[] argv) throws Exception {
        Main m = new Main (argv);
        m.exec();
    }
}
