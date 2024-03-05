package gov.nih.ncats.structureIndexer;

import static gov.nih.ncats.structureIndexer.StructureIndexer.FIELD_MOLWT;
import static gov.nih.ncats.structureIndexer.StructureIndexer.FIELD_NATOMS;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.apache.lucene.document.DoublePoint;
import org.apache.lucene.document.IntPoint;
import org.apache.lucene.search.Query;
import org.junit.Test;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.structureIndexer.StructureIndexer.Result;
import gov.nih.ncats.structureIndexer.StructureIndexer.ResultEnumeration;

public class Junit4StructureIndexerTest extends AbstractStructureIndexerTest {


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

        addAromatized("benzene", Chemical.createFromSmiles("c1ccccc1"));
        addAromatized("benzothiazole", Chemical.createFromSmiles("c1nc2ccccc2s1"));
        addAromatized("benzodiazole", Chemical.createFromSmiles("N1C=NC2=C1COC1CCCC[C@H]1O=CC=C2"));
        addAromatized("benzodiazole", Chemical.createFromSmiles("[nH]1cnc2ccccc12"));
        addAromatized("benzodiazole", Chemical.createFromSmiles("c1[nH]c2ccccc2n1"));


        ResultEnumeration result =
                indexer.substructure
                        ("c1ccccc1", IntPoint.newRangeQuery
                                        (FIELD_NATOMS, Math.addExact(6, 1), Integer.MAX_VALUE), // filter out benzene
                                DoublePoint.newRangeQuery // filter out benzodiazole
                                        (FIELD_MOLWT, 119D, Double.POSITIVE_INFINITY));
        //only expect 1 element
        assertTrue(result.hasMoreElements());
        Result rr = result.nextElement();

        assertEquals("benzothiazole", rr.getId());
        assertFalse("only 1 record should be returned", result.hasMoreElements());
        assertTrue(rr.similarity > 0.02 && rr.similarity < 1.0);


    }

    @Test
    public void substructureSearchForIdentityShouldMatch100Percent() throws Exception {

        addAromatized("benzene", Chemical.createFromSmiles("c1ccccc1"));
        addAromatized("benzothiazole", Chemical.createFromSmiles("c1nc2ccccc2s1"));
        addAromatized("benzodiazole", Chemical.createFromSmiles("N1C=NC2=C1COC1CCCC[C@H]1O=CC=C2"));
        addAromatized("benzodiazole", Chemical.createFromSmiles("[nH]1cnc2ccccc12"));
        addAromatized("benzodiazole", Chemical.createFromSmiles("c1[nH]c2ccccc2n1"));


        ResultEnumeration result =
                indexer.substructure
                        ("c1nc2ccccc2s1");
        //only expect 1 element
        assertTrue(result.hasMoreElements());
        Result rr = result.nextElement();

        assertEquals("benzothiazole", rr.getId());
        assertFalse("only 1 record should be returned", result.hasMoreElements());
        System.out.println("rr.similarity:" + rr.similarity);
        assertTrue(rr.similarity == 1.0); 


    }

    @Test
    public void substructureSearchFornearIdentityShouldMatchWithHighPercent() throws Exception {

        addAromatized("benzene", Chemical.createFromSmiles("c1ccccc1"));
        addAromatized("benzothiazoleEthyl", Chemical.createFromSmiles("CCC1=CC2=C(SC=N2)C=C1"));
        addAromatized("benzodiazole", Chemical.createFromSmiles("N1C=NC2=C1COC1CCCC[C@H]1O=CC=C2"));
        addAromatized("benzodiazole", Chemical.createFromSmiles("[nH]1cnc2ccccc12"));
        addAromatized("benzodiazole", Chemical.createFromSmiles("c1[nH]c2ccccc2n1"));


        ResultEnumeration result =
                indexer.substructure
                        ("CC1=CC2=C(SC=N2)C=C1");
        //only expect 1 element
        assertTrue(result.hasMoreElements());
        Result rr = result.nextElement();

        assertEquals("benzothiazoleEthyl", rr.getId());
        assertFalse("only 1 record should be returned", result.hasMoreElements());
        assertTrue(rr.similarity > 0.5);

    }

    private void addAromatized(String name, Chemical chem) throws IOException {
        chem.aromatize();
        indexer.add(name, chem);
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
                indexer.similarity("O=C1NC2=C(NC3=NC=CC=C13)N=CC=C2", 0.55);

        
        //should only get 1 result
        assertTrue(result.hasMoreElements());
        Result result1 = result.nextElement();
        assertEquals("two", result1.getId());
        //different fingerprinting algorithms will have different similarity
        assertTrue(Double.toString(result1.getSimilarity()), result1.getSimilarity() > .55D);
        assertFalse(result.hasMoreElements());

    }

    //
    @Test
    public void similaritySearchForVeryCloseThingsShouldStillBeDifferent() throws Exception {
        //C1=CC=C2N=CC=CC2=C1
        //C3=CC4=C(C=C3)C=NC=C4

        indexer.add("foo", "one", "C1=CC=C2C=NC=CC2=C1");
        ResultEnumeration result =
                indexer.similarity("C3=CC=C4N=CC=CC4=C3", 0.1);
        assertTrue(result.hasMoreElements());
        Result result1 = result.nextElement();
        assertEquals("one", result1.getId());
        assertTrue(Double.toString(result1.getSimilarity()), result1.getSimilarity() < 1.00D);
    }
    
    @Test
    public void similaritySearchForVeryCloseLongChainThingsShouldStillBeDifferent() throws Exception {
        indexer.add("foo", "one", "CCCCCCCCCCCCCCCCCCCOCCCCCCCCCCCCCCCCCC");
        ResultEnumeration result =
                indexer.similarity("CCCCCCCCCCCCOCCCCCCCCCCCCCCCCCCCCCCCC", 0.9);
        assertTrue(result.hasMoreElements());
        Result result1 = result.nextElement();
        assertEquals("one", result1.getId());
        assertTrue(Double.toString(result1.getSimilarity()), result1.getSimilarity() < 1.00D);
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
        Query filter = DoublePoint.newRangeQuery
                (FIELD_MOLWT, 110D, Double.POSITIVE_INFINITY);
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


        Chemical mol = Chemical.createFromSmilesAndComputeCoordinates("c1ccccc1");

        mol.setProperty("prop1", "foo");
        mol.setProperty("prop2", "123");
        mol.setProperty("prop3", "3.1415926535");
        indexer.add("zzz", "one", mol);

        Chemical mol2 = Chemical.createFromSmilesAndComputeCoordinates("c1ccncc1");

        mol2.setProperty("prop1", "456");
        mol2.setProperty("prop2", "bar");
        mol2.setProperty("prop3", "999");
        indexer.add("zzz", "two", mol2);

        ResultEnumeration result1 =
                indexer.search(DoublePoint.newRangeQuery
                        ("prop3", 3.0, DoublePoint.nextDown(4.0)));

        assertTrue(result1.hasMoreElements());
        assertEquals("one", result1.nextElement().getId());
        assertFalse(result1.hasMoreElements());


        ResultEnumeration result2 =
                indexer.search(IntPoint.newRangeQuery
                        ("prop3", Math.addExact(0, 1), Math.addExact(1000, -1)));

        assertTrue(result2.hasMoreElements());
        assertEquals("two", result2.nextElement().getId());
        assertFalse(result2.hasMoreElements());


        ResultEnumeration result3 = indexer.search("bar");
        assertTrue(result3.hasMoreElements());
        assertEquals("two", result3.nextElement().getId());
        assertFalse(result3.hasMoreElements());

    }

    @Test
    public void searchExplicitHsAsSmiles() throws Exception {
        indexer.add("one", "C1=CC2=CC=CC=C2C=C1");
        indexer.add("two", "CC1=C2C=CC=CC2=CC=C1");

        ResultEnumeration result =
                indexer.substructure("C(=CC=C1[H])(C(=C1C(=C2)[H])C(=C2)[H])[H]");

        assertTrue(result.hasMoreElements());
        assertEquals("one", result.nextElement().getId());
        assertFalse(result.hasMoreElements());
    }
    @Test
    public void searchExplicitHsAsMol() throws Exception {
        indexer.add("one", "C1=CC2=CC=CC=C2C=C1");
        indexer.add("two", "CC1=C2C=CC=CC2=CC=C1");

        String mol = "\n" +
                "   JSDraw210311918442D\n" +
                "\n" +
                " 14 15  0  0  0  0            999 V2000\n" +
                "   23.6080   -8.6840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   22.2570   -7.9040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   22.2570   -6.3440    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   24.9590   -7.9040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   24.9590   -6.3440    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   23.6080   -5.5640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   26.3100   -5.5640    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   26.3100   -8.6840    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   27.6610   -7.9040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   27.6610   -6.3440    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   23.6080   -4.0040    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   23.6080  -10.2440    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   26.3100   -4.0040    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   26.3100  -10.2440    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  2  0  0  0  0\n" +
                "  2  3  1  0  0  0  0\n" +
                "  1  4  1  0  0  0  0\n" +
                "  4  5  2  0  0  0  0\n" +
                "  5  6  1  0  0  0  0\n" +
                "  6  3  2  0  0  0  0\n" +
                "  5  7  1  0  0  0  0\n" +
                "  4  8  1  0  0  0  0\n" +
                "  8  9  2  0  0  0  0\n" +
                "  9 10  1  0  0  0  0\n" +
                " 10  7  2  0  0  0  0\n" +
                "  6 11  1  0  0  0  0\n" +
                "  1 12  1  0  0  0  0\n" +
                "  7 13  1  0  0  0  0\n" +
                "  8 14  1  0  0  0  0\n" +
                "M  END";

//        System.out.println(mol);
        ResultEnumeration result =
                indexer.substructure(mol);

        assertTrue(result.hasMoreElements());
        assertEquals("one", result.nextElement().getId());
        assertFalse(result.hasMoreElements());
    }

    @Test
    public void substructureDoubleOrAromatic() throws Exception {
        indexer.add("one", "C1=CC=CC=C1");
        indexer.add("two", "C=CC=C");
        indexer.add("two", "CCCCC");

        String mol = "\n" +
                "   JSDraw210311918512D\n" +
                "\n" +
                "  2  1  0  0  0  0            999 V2000\n" +
                "   28.2360   -9.3600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   29.5870   -8.5800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  7  0  0  0  0\n" +
                "M  END";

        ResultEnumeration result =
                indexer.substructure(mol);

        Set<String> actual = new HashSet<>();

        while(result.hasMoreElements()){
            actual.add(result.nextElement().getId());
        }
        Set<String> expected = new HashSet<>();
        expected.add("one");
        expected.add("two");

        assertEquals(expected, actual);
    }
    


    @Test
    public void gsrs1095() throws Exception {
        indexer.add("one", "C=CC1PONC2NOPCC12");
        indexer.add("two",  "C1PONC2NOPC(C12)C3=CC=CC=C3");

        String mol = "\n" +
                "   JSDraw212121918502D\n" +
                "\n" +
                " 12 13  0  0  0  0            999 V2000\n" +
                "   31.5619  -10.4259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   30.2111   -9.6459    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   30.2111   -8.0860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   31.5619   -7.3060    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   31.5619   -5.7461    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   30.2111   -4.9661    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   28.8601   -5.7461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   27.5091   -4.9661    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   26.1581   -5.7461    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   26.1581   -7.3060    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   27.5091   -8.0860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   28.8601   -7.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  4  0  0  0  0\n" +
                "  2  3  1  0  0  0  0\n" +
                "  3  4  1  0  0  0  0\n" +
                "  4  5  1  0  0  0  0\n" +
                "  5  6  1  0  0  0  0\n" +
                "  6  7  1  0  0  0  0\n" +
                "  7  8  1  0  0  0  0\n" +
                "  8  9  1  0  0  0  0\n" +
                "  9 10  1  0  0  0  0\n" +
                " 10 11  1  0  0  0  0\n" +
                " 11 12  1  0  0  0  0\n" +
                "  3 12  1  0  0  0  0\n" +
                "  7 12  1  0  0  0  0\n" +
                "M  END";

        ResultEnumeration result =
                indexer.substructure(mol);

        Set<String> actual = new HashSet<>();

        while(result.hasMoreElements()){
            actual.add(result.nextElement().getId());
        }
        Set<String> expected = new HashSet<>();
        expected.add("two");

        assertEquals(expected, actual);

        String mol2 = "\n" +
                "   JSDraw212121918502D\n" +
                "\n" +
                " 12 13  0  0  0  0            999 V2000\n" +
                "   31.5619  -10.4259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   30.2111   -9.6459    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   30.2111   -8.0860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   31.5619   -7.3060    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   31.5619   -5.7461    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   30.2111   -4.9661    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   28.8601   -5.7461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   27.5091   -4.9661    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   26.1581   -5.7461    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   26.1581   -7.3060    0.0000 P   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   27.5091   -8.0860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   28.8601   -7.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  4  0  0  0  0\n" +
                "  2  3  1  0  0  0  0\n" +
                "  3  4  1  0  0  0  0\n" +
                "  4  5  1  0  0  0  0\n" +
                "  5  6  1  0  0  0  0\n" +
                "  6  7  1  0  0  0  0\n" +
                "  7  8  1  0  0  0  0\n" +
                "  8  9  1  0  0  0  0\n" +
                "  9 10  1  0  0  0  0\n" +
                " 10 11  1  0  0  0  0\n" +
                " 11 12  1  0  0  0  0\n" +
                "  3 12  1  0  0  0  0\n" +
                "  7 12  1  0  0  0  0\n" +
                "M  ALS   2  1 F C   \n" +
                "M  END";

        ResultEnumeration result2 =
                indexer.substructure(mol2);

        Set<String> actual2 = new HashSet<>();

        while(result2.hasMoreElements()){
            actual2.add(result2.nextElement().getId());
        }



        assertEquals(expected, actual2);
    }

}
