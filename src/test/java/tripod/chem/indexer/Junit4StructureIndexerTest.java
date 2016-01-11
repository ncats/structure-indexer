package tripod.chem.indexer;

import gov.nih.ncats.chemkit.api.Chemical;

import org.apache.lucene.search.Filter;
import org.apache.lucene.search.NumericRangeFilter;
import org.apache.lucene.search.NumericRangeQuery;
import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.IOException;

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
      //  indexer.add("benzodiazole", "c1nc2ccccc2n1");
     //   indexer.add("benzodiazole", "[nH]1cnc2ccccc12");
        indexer.add("benzodiazole", "c1[nH]c2ccccc2n1");
        
     
        ResultEnumeration result =
                indexer.substructure
                        ("c1ccccc1", NumericRangeFilter.newIntRange
                                        (FIELD_NATOMS, 6, null, false, false), // filter out benzene
                                NumericRangeFilter.newDoubleRange // filter out benzodiazole
                                        (FIELD_MOLWT, 119D, null, true, false));
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
        //different fingerprinting algorithms will have different similarity
        //but it should be close for example 85% vs 81%
       // assertEquals(0.8580392156862745D, result1.getSimilarity(), 0.00001D);
        //dkatzel 2016-01-11
        //noticed structure-indexer tanimoto calculation was wrong should be ~ 76% for cdk, 80% for jchem
        assertTrue(Double.toString(result1.getSimilarity()), result1.getSimilarity() > .75D);
        assertFalse(result.hasMoreElements());

    }
    
    @Test
    public void correctlyFormattedSmilesAromatic() throws Exception {
    	indexer.add("foo", "one", "c1cccn1-c2ccc[nH]2");
       // createIndexerWithData();

        ResultEnumeration result =
                indexer.similarity("c1cccn1-c2ccc[nH]2", 0.5);

        //should only get 1 result
        assertTrue(result.hasMoreElements());
        Result result1 = result.nextElement();
        assertEquals("one", result1.getId());
        assertEquals(1.0D, result1.getSimilarity(), 0.00001D);

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

    	
		Chemical mol =Chemical.createFromSmiles("c1ccccc1");
       
        mol.setProperty("prop1", "foo");
        mol.setProperty("prop2", "123");
        mol.setProperty("prop3", "3.1415926535");
        indexer.add("zzz", "one", mol);
        
        Chemical mol2 =Chemical.createFromSmiles("c1ccncc1");
        
        mol2.setProperty("prop1", "456");
        mol2.setProperty("prop2", "bar");
        mol2.setProperty("prop3", "999");
        indexer.add("zzz", "two", mol2);

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
