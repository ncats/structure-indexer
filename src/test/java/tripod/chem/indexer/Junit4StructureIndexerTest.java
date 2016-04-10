package tripod.chem.indexer;

import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;
import org.apache.lucene.search.Filter;
import org.apache.lucene.search.NumericRangeFilter;
import org.apache.lucene.search.NumericRangeQuery;
import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.IOException;
import java.util.logging.Logger;

import static tripod.chem.indexer.StructureIndexer.*;
import static org.junit.Assert.*;

public class Junit4StructureIndexerTest {

    @Rule
    public TemporaryFolder tmpDir = new TemporaryFolder();

    private StructureIndexer indexer;

    @Before
    public void createIndexer() throws IOException {
        indexer = StructureIndexer.open(tmpDir.getRoot());
    }

    @After
    public void shutdownIndexer() {
        if (indexer != null) {
            indexer.shutdown();
        }
    }


    void createIndexerWithData() throws Exception {

        indexer.add("foo", "one", "c1ccccc1");
        indexer.add("bar", "one", "c1ccncc1");
        indexer.add("bar", "two", "OC1CCCC[C@H]1O");
        indexer.add("bar", "three", "CC(=O)Nc1ccc(cc1O)C(O)=O");
        indexer.add("abc", "one", "c1ncncc1");
        indexer.add("abc", "two", "Cc1cc(Cl)nc2N(C3CC3)c3ncccc3C(=O)Nc12");
    }

    @Test
    public void searchReturnsIndexedRecords() throws Exception {
        indexer.add("xxx", "benzene", "c1ccccc1");
        ResultEnumeration result = indexer.search("benzene", 1);
        assertTrue(result.hasMoreElements());
        assertEquals("benzene", result.nextElement().getId());

    }

    @Test
    public void removeRecordShouldNotBeReturnedFromSearch() throws Exception {

        indexer.add("foobar", "benzene", "c1ccccc1");
        indexer.remove("foobar", "benzene");
        ResultEnumeration result = indexer.search("benzene", 1);
        assertFalse(result.hasMoreElements());

    }

    @Test
    public void filteredSearchShouldOnlyReturnResultsNotFilteredOut() throws Exception {

        indexer.add("benzene", "c1ccccc1");
        indexer.add("benzothiazole", "c1nc2ccccc2s1");
        indexer.add("benzodiazole", "c1nc2ccccc2n1");
        ResultEnumeration result =
                indexer.substructure
                        ("c1ccccc1", NumericRangeFilter.newIntRange
                                        (FIELD_NATOMS, 6, null, false, false), // filter out benzene
                                NumericRangeFilter.newDoubleRange // filter out benzodiazole
                                        (FIELD_MOLWT, 118., null, true, false));
        //only expect 1 element
        assertTrue(result.hasMoreElements());
        assertEquals("benzothiazole", result.nextElement().getId());
        assertFalse("only 1 record should be returned", result.hasMoreElements());

    }

    @Test
    public void removeSourceWithMultipleRecords() throws Exception {

        indexer.add("foo", "one", "c1ccccc1");
        indexer.add("bar", "one", "c1ccncc1");
        indexer.add("bar", "two", "c1ncncc1");
        indexer.remove("bar");

        assertTrue(indexer.size() == 1);


    }

    @Test
    public void similaritySearchWithOnly1Result() throws Exception {

        createIndexerWithData();

        ResultEnumeration result =
                indexer.similarity("c1ccnc2Nc3ncccc3C(=O)Nc12", 0.5);

        //should only get 1 result
        assertTrue(result.hasMoreElements());
        Result result1 = result.nextElement();
        assertEquals("two", result1.getId());
        //dkatzel 2016-01-11
        //bug fix for tanimoto computation lowers similarity from 85 to 80%
        assertEquals(0.8063492063492064, result1.getSimilarity(), 0.00001D);

        assertFalse(result.hasMoreElements());

    }

    @Test
    public void similaritySearchWithMultipleResultsShouldReturnThemAll() throws Exception {
        createIndexerWithData();

        // filter for molwt >= 110da
        Filter filter = NumericRangeFilter.newDoubleRange
                (FIELD_MOLWT, 110D, null, true, false);
        ResultEnumeration result =
                indexer.similarity("c1ccnc2Nc3ncccc3C(=O)Nc12", 0, filter);
        int count = 0;
        while (result.hasMoreElements()) {
            Result r = result.nextElement();
            assertTrue(r.getMol().getMass() > 110D);
            ++count;
        }

        assertEquals(3, count);

    }

    @Test
    public void multipleSearches() throws Exception {

        MolHandler mh = new MolHandler("c1ccccc1");
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

        ResultEnumeration result1 =
                indexer.search(NumericRangeQuery.newDoubleRange
                        ("prop3", 3.0, 4.0, true, false));

        assertTrue(result1.hasMoreElements());
        assertEquals("one", result1.nextElement().getId());
        assertFalse(result1.hasMoreElements());


        ResultEnumeration result2 =
                indexer.search(NumericRangeQuery.newIntRange
                        ("prop3", 0, 1000, false, false));

        assertTrue(result2.hasMoreElements());
        assertEquals("two", result2.nextElement().getId());
        assertFalse(result2.hasMoreElements());


        ResultEnumeration result3 = indexer.search("bar");
        assertTrue(result3.hasMoreElements());
        assertEquals("two", result3.nextElement().getId());
        assertFalse(result3.hasMoreElements());

    }
}
