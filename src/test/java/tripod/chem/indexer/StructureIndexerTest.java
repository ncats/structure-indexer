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
import org.apache.lucene.search.NumericRangeQuery;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import static tripod.chem.indexer.StructureIndexer.*;

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
            ResultEnumeration result = indexer.search("benzene", 1);
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
            ResultEnumeration result = indexer.search("benzene", 1);
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
            indexer.add("benzene", "c1ccccc1");
            indexer.add("benzothiazole", "c1nc2ccccc2s1");
            indexer.add("benzodiazole", "c1nc2ccccc2n1");
            ResultEnumeration result =
                indexer.substructure
                ("c1ccccc1", NumericRangeFilter.newIntRange
                 (FIELD_NATOMS, 6, null, false, false), // filter out benzene
                 NumericRangeFilter.newDoubleRange // filter out benzodiazole
                 (FIELD_MOLWT, 118., null, true, false));
            int count = 0;
            Result r = null;
            while (result.hasMoreElements()) {
                r = result.nextElement();
                ++count;
            }
            assertTrue (count == 1 && r.getId().equals("benzothiazole"));
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
            ResultEnumeration result =
                indexer.similarity("c1ccnc2Nc3ncccc3C(=O)Nc12", 0.5);
            int c = 0;
            while (result.hasMoreElements()) {
                Result r = result.nextElement();
                logger.info(r.getId()+" "+r.getSimilarity());
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
                (FIELD_MOLWT, 110., null, true, false);
            ResultEnumeration result =
                indexer.similarity("c1ccnc2Nc3ncccc3C(=O)Nc12", 0., filter);
            int c = 0;
            while (result.hasMoreElements()) {
                Result r = result.nextElement();
                logger.info
                    (String.format("%1$7.3f", r.getMol().getMass())
                     +" "+r.getSimilarity()+" "+r.getId());
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

    public void test7 () {
        StructureIndexer indexer = null;
        try {
            indexer = createIndexer ();
            MolHandler mh = new MolHandler ("c1ccccc1");
            Molecule mol = mh.getMolecule();
            mol.setProperty("prop1", "foo");
            mol.setProperty("prop2", "123");
            mol.setProperty("prop3", "3.1415926535");
            indexer.add("zzz", "one", mol);
            mh.setMolecule("c1ccncc1");
            mol = mh.getMolecule();
            mol.setProperty("prop1", "456");
            mol.setProperty("prop2", "bar");
            mol.setProperty("prop3", "999");
            indexer.add("zzz", "two", mol);
            
            String[] fields = indexer.getFields();
            for (String f : fields) {
                logger.info("Field \""+f+"\"");
            }
            
            ResultEnumeration result1 =
                indexer.search(NumericRangeQuery.newDoubleRange
                               ("prop3", 3.0, 4.0, true, false));
            ResultEnumeration result2 =
                indexer.search(NumericRangeQuery.newIntRange
                               ("prop3", 0, 1000, false, false));
            ResultEnumeration result3 = indexer.search("bar");
            if (result1.hasMoreElements()
                && result2.hasMoreElements() && result3.hasMoreElements()) {
                Result r1 = result1.nextElement();
                Result r2 = result2.nextElement();
                assertTrue (r1.getId().equals("one")
                            && r2.getId().equals("two")
                            && result3.nextElement().getId().equals("two"));
            }
            else
                assertFalse (true);
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
