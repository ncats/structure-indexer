package tripod.chem.indexer;
import java.io.*;
import java.util.*;
import java.util.logging.Logger;
import java.util.logging.Level;

import chemaxon.util.MolHandler;
import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;
import org.apache.lucene.search.Filter;
import org.apache.lucene.search.NumericRangeFilter;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class StructureIndexerTest extends TestCase {
    static final Logger logger =
        Logger.getLogger(StructureIndexerTest.class.getName());
    
    public StructureIndexerTest () {
        super ("StructureIndexerTest");
    }
    
    public static Test suite () {
        return new TestSuite (StructureIndexerTest.class);
    }

    static String random (int size) {
        byte[] b = new byte[size];
        Random r = new Random ();
        r.nextBytes(b);
        StringBuilder sb = new StringBuilder ();
        for (int i = 0; i < b.length; ++i)
            sb.append(String.format("%1$02x", b[i] & 0xff));
        return sb.toString();
    }

    static StructureIndexer createIndexer () throws Exception {
        File temp = File.createTempFile("zzz", "");
        File dir = new File (temp.getParentFile(), random (3));
        dir.mkdir();
        temp.delete();
        return StructureIndexer.open(dir);
    }

    static StructureIndexer createIndexerWithData () throws Exception {
        StructureIndexer indexer = createIndexer ();
        indexer.add("foo", "one", "c1ccccc1");
        indexer.add("bar", "one", "c1ccncc1");
        indexer.add("bar", "two", "OC1CCCC[C@H]1O");
        indexer.add("bar", "three", "CC(=O)Nc1ccc(cc1O)C(O)=O");
        indexer.add("abc", "one", "c1ncncc1");
        indexer.add("abc", "two", "Cc1cc(Cl)nc2N(C3CC3)c3ncccc3C(=O)Nc12");
        return indexer;
    }
    
    public void test1 () {
        StructureIndexer indexer = null;
        try {
            indexer = createIndexer ();
            indexer.add("xxx", "benzene", "c1ccccc1");
            StructureIndexer.ResultEnumeration result =
                indexer.search("benzene", 1);
            assertTrue (result.hasMoreElements());
        }
        catch (Exception ex) {
            ex.printStackTrace();
            assertTrue (false);
        }
        finally {
            if (indexer != null)
                indexer.shutdown();
        }
    }

    public void test2 () {
        StructureIndexer indexer = null;
        try {
            indexer = createIndexer ();
            indexer.add("foobar", "benzene", "c1ccccc1");
            indexer.stats(System.out);
            indexer.remove("foobar", "benzene");
            indexer.stats(System.out);
            StructureIndexer.ResultEnumeration result =
                indexer.search("benzene", 1);
            assertFalse (result.hasMoreElements());
        }
        catch (Exception ex) {
            ex.printStackTrace();
            assertTrue (false);
        }
        finally {
            if (indexer != null)
                indexer.shutdown();
        }
    }

    public void test3 () {
        StructureIndexer indexer = null;
        try {
            indexer = createIndexer ();
            indexer.add("foobar", "benzene", "c1ccccc1");
            StructureIndexer.ResultEnumeration result =
                indexer.substructure("c1ccccc1", 1, 1);
            assertTrue (result.hasMoreElements());
        }
        catch (Exception ex) {
            ex.printStackTrace();
            assertTrue (false);
        }
        finally {
            if (indexer != null)
                indexer.shutdown();
        }
    }

    public void test4 () {
        StructureIndexer indexer = null;
        try {
            indexer = createIndexer ();
            indexer.add("foo", "one", "c1ccccc1");
            indexer.add("bar", "one", "c1ccncc1");
            indexer.add("bar", "two", "c1ncncc1");
            indexer.remove("bar");
            /*
            System.out.println("Waiting for background thread to start...");
            Thread.currentThread().wait(5000);
            for (StructureIndexer.Codebook cb : indexer.getCodebooks()) {
                for (int i = 1; i < cb.size(); ++i)
                    assertTrue (cb.count(i) <= 1);
            }
            */
            assertTrue (indexer.size() == 1);
        }
        catch (Exception ex) {
            ex.printStackTrace();
            assertTrue (false);
        }
        finally {
            if (indexer != null) {
                indexer.shutdown();
            }
        }
    }

    public void test5 () {
        StructureIndexer indexer = null;
        try {
            indexer = createIndexerWithData ();
            logger.info("Index size: "+indexer.size());
            StructureIndexer.ResultEnumeration result =
                indexer.similarity("c1ccnc2Nc3ncccc3C(=O)Nc12", 0.5);
            int c = 0;
            while (result.hasMoreElements()) {
                StructureIndexer.Result r = result.nextElement();
                logger.info(r.getMol().getProperty("TANIMOTO"));
                ++c;
            }
            
            assertTrue (c == 1);
        }
        catch (Exception ex) {
            ex.printStackTrace();
            assertTrue (false);
        }
        finally {
            if (indexer != null)
                indexer.shutdown();
        }
    }

    public void test6 () {
        StructureIndexer indexer = null;
        try {
            indexer = createIndexerWithData ();
            logger.info("Index size: "+indexer.size());
            // filter for molwt >= 110da
            Filter filter = NumericRangeFilter.newDoubleRange
                (StructureIndexer.FIELD_MOLWT, 110., null, true, false);
            StructureIndexer.ResultEnumeration result =
                indexer.similarity(filter, "c1ccnc2Nc3ncccc3C(=O)Nc12", 0.);
            int c = 0;
            while (result.hasMoreElements()) {
                StructureIndexer.Result r = result.nextElement();
                logger.info
                    (String.format("%1$7.3f", r.getMol().getMass())
                     +" "+r.getMol().getProperty("TANIMOTO"));
                ++c;
            }
            
            assertTrue (c == 3);
        }
        catch (Exception ex) {
            ex.printStackTrace();
            assertTrue (false);
        }
        finally {
            if (indexer != null)
                indexer.shutdown();
        }
    }
}
