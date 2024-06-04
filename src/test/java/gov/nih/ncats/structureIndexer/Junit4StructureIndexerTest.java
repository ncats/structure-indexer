package gov.nih.ncats.structureIndexer;

import static gov.nih.ncats.structureIndexer.StructureIndexer.FIELD_MOLWT;
import static gov.nih.ncats.structureIndexer.StructureIndexer.FIELD_NATOMS;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import gov.nih.ncats.common.io.IOUtil;
import org.apache.commons.io.IOUtils;
import org.apache.lucene.search.NumericRangeQuery;
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
                        ("c1ccccc1", NumericRangeQuery.newIntRange
                                        (FIELD_NATOMS, 6, null, false, false), // filter out benzene
                                NumericRangeQuery.newDoubleRange // filter out benzodiazole
                                        (FIELD_MOLWT, 119D, null, true, false));
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
        Query filter = NumericRangeQuery.newDoubleRange
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

    @Test
    public void addLargeMol() throws Exception {
        String id = "VFB8SJP4RO";
        //String molfileText = IOUtils.toString(this.getClass().getResourceAsStream("mols/VFB8SJP4RO.mol"));
        Chemical testMol = Chemical.parseMol(VFB8SJP4ROMolfileText);
        indexer.add("GSRS", id, testMol);
        ResultEnumeration result = indexer.search(id, 1);
        assertTrue(result.hasMoreElements());
    }

    String VFB8SJP4ROMolfileText = "DALBAVANCIN B1\n" +
            "  Marvin  01132101192D          \n" +
            "\n" +
            "653705  0  0  1  0            999 V2000\n" +
            "   27.7485   -7.5615    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.7021   -2.2680    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   32.0476   -2.2537    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   28.4555   -7.9718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.4529   -2.9703    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   34.1688   -9.1266    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6024   -9.1502    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.1074   -8.4760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.2826   -2.9703    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1688   -7.9575    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   34.8853   -7.5426    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -7.9623    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.6004   -7.5521    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3236  -12.9449    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.4570   -7.5475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6071  -10.8942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1722   -7.5569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9000  -10.4747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3169   -7.9669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.5977  -13.3502    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.2967   -1.5514    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   29.1532  -10.0175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.4671   -1.5514    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   32.0382   -4.5023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9000   -9.7865    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   37.7419   -7.5333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.5977  -14.1752    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   32.0334   -7.5475    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   37.0302   -7.9482    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3263   -4.9124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -7.9529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8838   -7.9718    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -4.9077    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0354  -13.3644    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6071  -11.7099    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0382   -3.6820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0334   -6.1428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3141  -14.5900    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   30.6098   -4.5023    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4475   -9.5321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0354  -14.1846    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   36.3137   -7.5380    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   27.7485   -6.1522    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   29.1532   -9.1926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3236  -12.1200    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   36.3042  -10.3427    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8838   -8.7872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.5316   -2.2680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4618   -4.4929    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7277  -10.3427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.2132   -2.2537    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.4649   -5.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4811   -9.9657    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4475  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8759  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4537   -9.1172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4537   -7.9434    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   31.3263   -5.7326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8979   -4.9171    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -5.7373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1546  -10.7576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8954  -12.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7277   -9.5225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8759   -9.5367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -4.5070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.4649   -4.9218    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8698  -10.4324    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5014   -1.8342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1787  -10.8801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1506   -8.3724    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4570   -6.7225    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1722   -6.7366    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3169   -8.7872    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4395  -10.7576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7419   -6.7131    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1787  -11.7003    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -8.7730    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   36.3137   -6.1380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5910   -9.1973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -6.1522    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1654   -9.5273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.9420   -1.5561    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.5565  -10.7764    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8979   -5.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1735   -4.9030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5014   -1.0094    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5910  -10.0222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1654  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8811  -12.9354    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.4538   -9.6438    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0570   -0.8255    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.7162   -0.8396    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8765  -14.5854    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7310   -9.1172    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1702   -7.5286    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3141  -15.4149    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -5.7279    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -3.6868    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.0272   -5.7420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.9420   -2.9797    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7216  -10.7529    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.8658   -9.4629    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8698  -11.2525    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7472  -14.5995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4395  -11.5826    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4575  -12.1106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8853   -4.4929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1735   -5.7231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8806   -6.1333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -4.9030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.9626   -9.6574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7472  -15.4197    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7848   -2.2444    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   23.6063  -10.2461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.2693  -10.1060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.8867   -7.9386    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.8070   -8.5101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.0330   -9.7909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.0669   -2.2633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0732   -1.8390    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.7883   -1.8484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4999   -2.2539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3565   -2.2537    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.2165   -1.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9282   -2.2539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6448   -1.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   23.3505   -1.8484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.4583   -1.6645    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.5471   -2.4395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   46.8149   -8.6215    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.7485   -7.5615    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.7021   -2.2680    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   32.0476   -2.2537    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   28.4555   -7.9718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.4529   -2.9703    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   34.1688   -9.1266    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6024   -9.1502    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.1074   -8.4760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.2826   -2.9703    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1688   -7.9575    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   34.8853   -7.5426    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -7.9623    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.6004   -7.5521    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3236  -12.9449    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.4570   -7.5475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6071  -10.8942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1722   -7.5569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9000  -10.4747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3169   -7.9669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.5977  -13.3502    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.2967   -1.5514    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   29.1532  -10.0175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.4671   -1.5514    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   32.0382   -4.5023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9000   -9.7865    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   37.7419   -7.5333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.5977  -14.1752    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   32.0334   -7.5475    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   37.0302   -7.9482    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3263   -4.9124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -7.9529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8838   -7.9718    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -4.9077    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0354  -13.3644    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6071  -11.7099    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0382   -3.6820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0334   -6.1428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3141  -14.5900    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   30.6098   -4.5023    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4475   -9.5321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0354  -14.1846    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   36.3137   -7.5380    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   27.7485   -6.1522    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   29.1532   -9.1926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3236  -12.1200    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   36.3042  -10.3427    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8838   -8.7872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.5316   -2.2680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4618   -4.4929    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7277  -10.3427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.2132   -2.2537    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.4649   -5.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4811   -9.9657    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4475  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8759  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4537   -9.1172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4537   -7.9434    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   31.3263   -5.7326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8979   -4.9171    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -5.7373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1546  -10.7576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8954  -12.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7277   -9.5225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8759   -9.5367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -4.5070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.4649   -4.9218    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8698  -10.4324    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5014   -1.8342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1787  -10.8801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1506   -8.3724    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4570   -6.7225    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1722   -6.7366    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3169   -8.7872    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4395  -10.7576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7419   -6.7131    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1787  -11.7003    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -8.7730    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   36.3137   -6.1380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5910   -9.1973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -6.1522    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1654   -9.5273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.9420   -1.5561    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.5565  -10.7764    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8979   -5.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1735   -4.9030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5014   -1.0094    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5910  -10.0222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1654  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8811  -12.9354    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.4538   -9.6438    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0570   -0.8255    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.7162   -0.8396    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8765  -14.5854    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7310   -9.1172    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1702   -7.5286    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3141  -15.4149    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -5.7279    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -3.6868    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.0272   -5.7420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.9420   -2.9797    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7216  -10.7529    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.8658   -9.4629    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8698  -11.2525    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7472  -14.5995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4395  -11.5826    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4575  -12.1106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8853   -4.4929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1735   -5.7231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8806   -6.1333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -4.9030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.9626   -9.6574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7472  -15.4197    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7848   -2.2444    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   23.6063  -10.2461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.2693  -10.1060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.8867   -7.9386    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.8070   -8.5101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.0330   -9.7909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.0669   -2.2633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0732   -1.8390    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.7883   -1.8484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4999   -2.2539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3565   -2.2537    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.2165   -1.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9282   -2.2539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6448   -1.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   23.3505   -1.8484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.4583   -1.6645    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.5471   -2.4395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.7485   -7.5615    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.7021   -2.2680    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   32.0476   -2.2537    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   28.4555   -7.9718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.4529   -2.9703    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   34.1688   -9.1266    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6024   -9.1502    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.1074   -8.4760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.2826   -2.9703    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1688   -7.9575    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   34.8853   -7.5426    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -7.9623    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.6004   -7.5521    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3236  -12.9449    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.4570   -7.5475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6071  -10.8942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1722   -7.5569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9000  -10.4747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3169   -7.9669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.5977  -13.3502    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.2967   -1.5514    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   29.1532  -10.0175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.4671   -1.5514    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   32.0382   -4.5023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9000   -9.7865    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   37.7419   -7.5333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.5977  -14.1752    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   32.0334   -7.5475    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   37.0302   -7.9482    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3263   -4.9124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -7.9529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8838   -7.9718    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -4.9077    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0354  -13.3644    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6071  -11.7099    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0382   -3.6820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0334   -6.1428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3141  -14.5900    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   30.6098   -4.5023    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4475   -9.5321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0354  -14.1846    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   36.3137   -7.5380    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   27.7485   -6.1522    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   29.1532   -9.1926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3236  -12.1200    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   36.3042  -10.3427    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8838   -8.7872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.5316   -2.2680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4618   -4.4929    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7277  -10.3427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.2132   -2.2537    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.4649   -5.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4811   -9.9657    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4475  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8759  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4537   -9.1172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4537   -7.9434    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   31.3263   -5.7326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8979   -4.9171    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -5.7373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1546  -10.7576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8954  -12.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7277   -9.5225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8759   -9.5367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -4.5070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.4649   -4.9218    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8698  -10.4324    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5014   -1.8342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1787  -10.8801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1506   -8.3724    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4570   -6.7225    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1722   -6.7366    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3169   -8.7872    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4395  -10.7576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7419   -6.7131    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1787  -11.7003    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -8.7730    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   36.3137   -6.1380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5910   -9.1973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -6.1522    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1654   -9.5273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.9420   -1.5561    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.5565  -10.7764    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8979   -5.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1735   -4.9030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5014   -1.0094    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5910  -10.0222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1654  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8811  -12.9354    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.4538   -9.6438    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0570   -0.8255    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.7162   -0.8396    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8765  -14.5854    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7310   -9.1172    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1702   -7.5286    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3141  -15.4149    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -5.7279    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -3.6868    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.0272   -5.7420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.9420   -2.9797    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7216  -10.7529    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.8658   -9.4629    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8698  -11.2525    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7472  -14.5995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4395  -11.5826    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4575  -12.1106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8853   -4.4929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1735   -5.7231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8806   -6.1333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -4.9030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.9626   -9.6574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7472  -15.4197    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7848   -2.2444    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   23.6063  -10.2461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.2693  -10.1060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.8867   -7.9386    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.8070   -8.5101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.0330   -9.7909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.0669   -2.2633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0732   -1.8390    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.7883   -1.8484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4999   -2.2539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3565   -2.2537    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.2165   -1.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9282   -2.2539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6448   -1.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   23.3505   -1.8484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.4583   -1.6645    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.5471   -2.4395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.7485   -7.5615    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.7021   -2.2680    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   32.0476   -2.2537    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   28.4555   -7.9718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.4529   -2.9703    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   34.1688   -9.1266    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6024   -9.1502    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.1074   -8.4760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.2826   -2.9703    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1688   -7.9575    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   34.8853   -7.5426    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -7.9623    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.6004   -7.5521    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3236  -12.9449    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.4570   -7.5475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6071  -10.8942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1722   -7.5569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9000  -10.4747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3169   -7.9669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.5977  -13.3502    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.2967   -1.5514    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   29.1532  -10.0175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.4671   -1.5514    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   32.0382   -4.5023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9000   -9.7865    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   37.7419   -7.5333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.5977  -14.1752    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   32.0334   -7.5475    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   37.0302   -7.9482    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3263   -4.9124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -7.9529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8838   -7.9718    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -4.9077    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0354  -13.3644    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6071  -11.7099    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0382   -3.6820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0334   -6.1428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3141  -14.5900    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   30.6098   -4.5023    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4475   -9.5321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0354  -14.1846    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   36.3137   -7.5380    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   27.7485   -6.1522    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   29.1532   -9.1926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3236  -12.1200    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   36.3042  -10.3427    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8838   -8.7872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.5316   -2.2680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4618   -4.4929    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7277  -10.3427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.2132   -2.2537    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.4649   -5.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4811   -9.9657    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4475  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8759  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4537   -9.1172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4537   -7.9434    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   31.3263   -5.7326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8979   -4.9171    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -5.7373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1546  -10.7576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8954  -12.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7277   -9.5225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8759   -9.5367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -4.5070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.4649   -4.9218    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8698  -10.4324    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5014   -1.8342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1787  -10.8801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1506   -8.3724    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4570   -6.7225    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1722   -6.7366    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3169   -8.7872    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4395  -10.7576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7419   -6.7131    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1787  -11.7003    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -8.7730    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   36.3137   -6.1380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5910   -9.1973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -6.1522    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1654   -9.5273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.9420   -1.5561    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.5565  -10.7764    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8979   -5.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1735   -4.9030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5014   -1.0094    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5910  -10.0222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1654  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8811  -12.9354    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.4538   -9.6438    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0570   -0.8255    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.7162   -0.8396    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8765  -14.5854    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7310   -9.1172    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1702   -7.5286    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3141  -15.4149    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -5.7279    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -3.6868    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.0272   -5.7420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.9420   -2.9797    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7216  -10.7529    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.8658   -9.4629    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8698  -11.2525    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7472  -14.5995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4395  -11.5826    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4575  -12.1106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8853   -4.4929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1735   -5.7231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8806   -6.1333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -4.9030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.9626   -9.6574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7472  -15.4197    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7848   -2.2444    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   23.6063  -10.2461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.2693  -10.1060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.8867   -7.9386    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.8070   -8.5101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.0330   -9.7909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.0669   -2.2633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0732   -1.8390    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.7883   -1.8484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4999   -2.2539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3565   -2.2537    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.2165   -1.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9282   -2.2539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6448   -1.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   23.3505   -1.8484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.4583   -1.6645    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.5471   -2.4395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.7485   -7.5615    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.7021   -2.2680    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   32.0476   -2.2537    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   28.4555   -7.9718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.4529   -2.9703    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   34.1688   -9.1266    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6024   -9.1502    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.1074   -8.4760    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.2826   -2.9703    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1688   -7.9575    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   34.8853   -7.5426    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -7.9623    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.6004   -7.5521    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3236  -12.9449    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.4570   -7.5475    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6071  -10.8942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1722   -7.5569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9000  -10.4747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3169   -7.9669    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.5977  -13.3502    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   33.2967   -1.5514    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   29.1532  -10.0175    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.4671   -1.5514    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   32.0382   -4.5023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9000   -9.7865    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   37.7419   -7.5333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.5977  -14.1752    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0\n" +
            "   32.0334   -7.5475    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   37.0302   -7.9482    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3263   -4.9124    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -7.9529    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8838   -7.9718    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -4.9077    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0354  -13.3644    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6071  -11.7099    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0382   -3.6820    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0334   -6.1428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3141  -14.5900    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   30.6098   -4.5023    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4475   -9.5321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0354  -14.1846    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   36.3137   -7.5380    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   27.7485   -6.1522    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   29.1532   -9.1926    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3236  -12.1200    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   36.3042  -10.3427    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8838   -8.7872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.5316   -2.2680    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4618   -4.4929    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7277  -10.3427    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.2132   -2.2537    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.4649   -5.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4811   -9.9657    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4475  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8759  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4537   -9.1172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4537   -7.9434    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
            "   31.3263   -5.7326    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8979   -4.9171    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7452   -5.7373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1546  -10.7576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8954  -12.1152    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7277   -9.5225    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8759   -9.5367    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -4.5070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.4649   -4.9218    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8698  -10.4324    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5014   -1.8342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1787  -10.8801    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1506   -8.3724    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.4570   -6.7225    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1722   -6.7366    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.3169   -8.7872    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4395  -10.7576    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   37.7419   -6.7131    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.1787  -11.7003    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -8.7730    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   36.3137   -6.1380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5910   -9.1973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -6.1522    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1654   -9.5273    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.9420   -1.5561    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.5565  -10.7764    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8979   -5.7420    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1735   -4.9030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5014   -1.0094    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   30.5910  -10.0222    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1654  -10.3474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8811  -12.9354    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.4538   -9.6438    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.0570   -0.8255    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   33.7162   -0.8396    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.8765  -14.5854    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7310   -9.1172    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.1702   -7.5286    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3141  -15.4149    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -5.7279    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.1815   -3.6868    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.0272   -5.7420    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.9420   -2.9797    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   32.7216  -10.7529    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.8658   -9.4629    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.8698  -11.2525    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7472  -14.5995    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   38.4395  -11.5826    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4575  -12.1106    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8853   -4.4929    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.1735   -5.7231    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   34.8806   -6.1333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   35.6018   -4.9030    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.9626   -9.6574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7472  -15.4197    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.7848   -2.2444    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   23.6063  -10.2461    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.2693  -10.1060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   39.8867   -7.9386    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.8070   -8.5101    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   21.0330   -9.7909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.0669   -2.2633    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   29.0732   -1.8390    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   24.7883   -1.8484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   25.4999   -2.2539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   28.3565   -2.2537    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.2165   -1.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   26.9282   -2.2539    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   27.6448   -1.8437    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   23.3505   -1.8484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   31.4583   -1.6645    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   22.5471   -2.4395    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   46.8149   -8.6215    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   46.8149   -8.6215    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   46.8149   -8.6215    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   46.8149   -8.6215    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   46.8149   -8.6215    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   46.8149   -8.6215    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "   46.8149   -8.6215    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "  2  9  1  0  0  0  0\n" +
            "  3  5  1  0  0  0  0\n" +
            "  1  4  1  6  0  0  0\n" +
            "  5 36  1  1  0  0  0\n" +
            "  6 10  1  0  0  0  0\n" +
            "  7  8  1  0  0  0  0\n" +
            "  8  1  1  0  0  0  0\n" +
            "  9  5  1  0  0  0  0\n" +
            " 10 15  1  0  0  0  0\n" +
            " 10 11  1  1  0  0  0\n" +
            " 28 12  1  1  0  0  0\n" +
            " 13 32  1  0  0  0  0\n" +
            " 14 45  1  1  0  0  0\n" +
            " 15 12  1  0  0  0  0\n" +
            " 16 18  2  0  0  0  0\n" +
            " 17  4  1  0  0  0  0\n" +
            " 18 25  1  0  0  0  0\n" +
            " 19 13  1  0  0  0  0\n" +
            " 20 14  1  0  0  0  0\n" +
            " 21 23  1  0  0  0  0\n" +
            " 22 16  1  0  0  0  0\n" +
            " 23  3  1  0  0  0  0\n" +
            " 24 30  1  0  0  0  0\n" +
            " 25  7  1  0  0  0  0\n" +
            " 26 29  1  0  0  0  0\n" +
            " 27 20  1  0  0  0  0\n" +
            " 28 19  1  0  0  0  0\n" +
            " 42 29  1  1  0  0  0\n" +
            " 30 39  1  0  0  0  0\n" +
            " 31 11  1  0  0  0  0\n" +
            " 32 17  1  6  0  0  0\n" +
            " 33 60  1  0  0  0  0\n" +
            " 34 14  1  0  0  0  0\n" +
            " 35 16  1  0  0  0  0\n" +
            " 36 24  1  0  0  0  0\n" +
            " 37 28  1  0  0  0  0\n" +
            " 38 27  1  0  0  0  0\n" +
            " 39 59  1  0  0  0  0\n" +
            " 40  6  2  0  0  0  0\n" +
            " 41 34  1  0  0  0  0\n" +
            " 42 31  1  0  0  0  0\n" +
            " 43  1  1  0  0  0  0\n" +
            " 44 47  2  0  0  0  0\n" +
            " 45 35  1  0  0  0  0\n" +
            " 46 55  1  0  0  0  0\n" +
            " 47 32  1  0  0  0  0\n" +
            "  2 48  1  1  0  0  0\n" +
            " 49 33  1  0  0  0  0\n" +
            " 50 46  1  0  0  0  0\n" +
            " 51  3  1  0  0  0  0\n" +
            " 52 43  1  0  0  0  0\n" +
            " 25 53  1  1  0  0  0\n" +
            " 54 40  1  0  0  0  0\n" +
            " 55 64  2  0  0  0  0\n" +
            " 56 57  1  0  0  0  0\n" +
            " 57 26  1  0  0  0  0\n" +
            " 58 37  1  0  0  0  0\n" +
            " 59 84  1  0  0  0  0\n" +
            " 60 37  2  0  0  0  0\n" +
            " 61 55  1  0  0  0  0\n" +
            " 62 76  1  0  0  0  0\n" +
            " 63 50  1  0  0  0  0\n" +
            " 64  6  1  0  0  0  0\n" +
            " 65 66  1  0  0  0  0\n" +
            " 66 52  2  0  0  0  0\n" +
            " 67 22  2  0  0  0  0\n" +
            " 68 51  1  0  0  0  0\n" +
            " 69 18  1  0  0  0  0\n" +
            " 70  8  2  0  0  0  0\n" +
            " 71 15  2  0  0  0  0\n" +
            " 72 17  2  0  0  0  0\n" +
            " 73 19  2  0  0  0  0\n" +
            " 74 50  2  0  0  0  0\n" +
            " 75 26  2  0  0  0  0\n" +
            " 76 69  2  0  0  0  0\n" +
            " 77 31  2  0  0  0  0\n" +
            " 78 42  1  0  0  0  0\n" +
            " 79 47  1  0  0  0  0\n" +
            " 80 52  1  0  0  0  0\n" +
            " 81 56  1  0  0  0  0\n" +
            " 82 48  2  0  0  0  0\n" +
            " 83 53  2  0  0  0  0\n" +
            " 84 80  2  0  0  0  0\n" +
            " 85 49  1  0  0  0  0\n" +
            " 86 68  2  0  0  0  0\n" +
            " 87 79  2  0  0  0  0\n" +
            " 88 74  1  0  0  0  0\n" +
            " 20 89  1  6  0  0  0\n" +
            " 90 53  1  0  0  0  0\n" +
            " 23 91  1  1  0  0  0\n" +
            " 21 92  1  6  0  0  0\n" +
            " 27 93  1  6  0  0  0\n" +
            " 94 40  1  0  0  0  0\n" +
            " 57 95  1  6  0  0  0\n" +
            " 38 96  1  1  0  0  0\n" +
            " 97 78  1  0  0  0  0\n" +
            " 98 65  1  0  0  0  0\n" +
            " 43 99  1  6  0  0  0\n" +
            "100 48  1  0  0  0  0\n" +
            "101 54  1  0  0  0  0\n" +
            "102115  1  0  0  0  0\n" +
            "103 67  1  0  0  0  0\n" +
            " 41104  1  6  0  0  0\n" +
            "105 74  1  0  0  0  0\n" +
            "106 76  1  0  0  0  0\n" +
            "107 85  1  0  0  0  0\n" +
            "108 85  2  0  0  0  0\n" +
            "109108  1  0  0  0  0\n" +
            "110107  2  0  0  0  0\n" +
            "111114  1  0  0  0  0\n" +
            "112104  1  0  0  0  0\n" +
            "113 68  1  0  0  0  0\n" +
            "114 90  1  0  0  0  0\n" +
            "115111  1  0  0  0  0\n" +
            "116 95  1  0  0  0  0\n" +
            "117102  1  0  0  0  0\n" +
            "118102  1  0  0  0  0\n" +
            "119121  1  0  0  0  0\n" +
            "120113  1  0  0  0  0\n" +
            "121122  1  0  0  0  0\n" +
            "122124  1  0  0  0  0\n" +
            "123120  1  0  0  0  0\n" +
            "124125  1  0  0  0  0\n" +
            "125126  1  0  0  0  0\n" +
            "126123  1  0  0  0  0\n" +
            "127119  1  0  0  0  0\n" +
            "  3128  1  1  0  0  0\n" +
            " 65 59  2  0  0  0  0\n" +
            " 22 44  1  0  0  0  0\n" +
            " 87 67  1  0  0  0  0\n" +
            " 35 62  2  0  0  0  0\n" +
            " 58 30  2  0  0  0  0\n" +
            " 24 33  2  0  0  0  0\n" +
            " 41 38  1  0  0  0  0\n" +
            " 61 54  2  0  0  0  0\n" +
            "110 97  1  0  0  0  0\n" +
            " 97109  2  0  0  0  0\n" +
            "  2 21  1  0  0  0  0\n" +
            " 63 56  2  0  0  0  0\n" +
            " 81 88  2  0  0  0  0\n" +
            "127129  1  0  0  0  0\n" +
            "518521  1  6  0  0  0\n" +
            "525518  1  0  0  0  0\n" +
            "560518  1  0  0  0  0\n" +
            "519526  1  0  0  0  0\n" +
            "519565  1  1  0  0  0\n" +
            "519538  1  0  0  0  0\n" +
            "520522  1  0  0  0  0\n" +
            "540520  1  0  0  0  0\n" +
            "568520  1  0  0  0  0\n" +
            "520645  1  1  0  0  0\n" +
            "534521  1  0  0  0  0\n" +
            "522553  1  1  0  0  0\n" +
            "526522  1  0  0  0  0\n" +
            "523527  1  0  0  0  0\n" +
            "557523  2  0  0  0  0\n" +
            "581523  1  0  0  0  0\n" +
            "524525  1  0  0  0  0\n" +
            "542524  1  0  0  0  0\n" +
            "587525  2  0  0  0  0\n" +
            "527532  1  0  0  0  0\n" +
            "527528  1  1  0  0  0\n" +
            "548528  1  0  0  0  0\n" +
            "545529  1  1  0  0  0\n" +
            "532529  1  0  0  0  0\n" +
            "530549  1  0  0  0  0\n" +
            "536530  1  0  0  0  0\n" +
            "531562  1  1  0  0  0\n" +
            "537531  1  0  0  0  0\n" +
            "551531  1  0  0  0  0\n" +
            "588532  2  0  0  0  0\n" +
            "533535  2  0  0  0  0\n" +
            "539533  1  0  0  0  0\n" +
            "552533  1  0  0  0  0\n" +
            "549534  1  6  0  0  0\n" +
            "589534  2  0  0  0  0\n" +
            "535542  1  0  0  0  0\n" +
            "586535  1  0  0  0  0\n" +
            "545536  1  0  0  0  0\n" +
            "590536  2  0  0  0  0\n" +
            "544537  1  0  0  0  0\n" +
            "537606  1  6  0  0  0\n" +
            "538540  1  0  0  0  0\n" +
            "538609  1  6  0  0  0\n" +
            "584539  2  0  0  0  0\n" +
            "539561  1  0  0  0  0\n" +
            "540608  1  1  0  0  0\n" +
            "541547  1  0  0  0  0\n" +
            "553541  1  0  0  0  0\n" +
            "541550  2  0  0  0  0\n" +
            "542570  1  1  0  0  0\n" +
            "543546  1  0  0  0  0\n" +
            "574543  1  0  0  0  0\n" +
            "592543  2  0  0  0  0\n" +
            "555544  1  0  0  0  0\n" +
            "544610  1  6  0  0  0\n" +
            "554545  1  0  0  0  0\n" +
            "559546  1  1  0  0  0\n" +
            "547556  1  0  0  0  0\n" +
            "575547  2  0  0  0  0\n" +
            "559548  1  0  0  0  0\n" +
            "594548  2  0  0  0  0\n" +
            "564549  1  0  0  0  0\n" +
            "550577  1  0  0  0  0\n" +
            "566550  1  0  0  0  0\n" +
            "558551  1  0  0  0  0\n" +
            "562552  1  0  0  0  0\n" +
            "552579  2  0  0  0  0\n" +
            "575554  1  0  0  0  0\n" +
            "577554  2  0  0  0  0\n" +
            "555613  1  1  0  0  0\n" +
            "558555  1  0  0  0  0\n" +
            "556576  1  0  0  0  0\n" +
            "571557  1  0  0  0  0\n" +
            "611557  1  0  0  0  0\n" +
            "558621  1  6  0  0  0\n" +
            "595559  1  0  0  0  0\n" +
            "569560  1  0  0  0  0\n" +
            "560616  1  6  0  0  0\n" +
            "561564  2  0  0  0  0\n" +
            "563572  1  0  0  0  0\n" +
            "567563  1  0  0  0  0\n" +
            "596564  1  0  0  0  0\n" +
            "599565  2  0  0  0  0\n" +
            "617565  1  0  0  0  0\n" +
            "602566  1  0  0  0  0\n" +
            "580567  1  0  0  0  0\n" +
            "591567  2  0  0  0  0\n" +
            "585568  1  0  0  0  0\n" +
            "583569  2  0  0  0  0\n" +
            "597569  1  0  0  0  0\n" +
            "600570  2  0  0  0  0\n" +
            "607570  1  0  0  0  0\n" +
            "618571  1  0  0  0  0\n" +
            "578571  2  0  0  0  0\n" +
            "572581  2  0  0  0  0\n" +
            "578572  1  0  0  0  0\n" +
            "573574  1  0  0  0  0\n" +
            "598573  1  0  0  0  0\n" +
            "580573  2  0  0  0  0\n" +
            "574612  1  6  0  0  0\n" +
            "576601  1  0  0  0  0\n" +
            "582576  2  0  0  0  0\n" +
            "579593  1  0  0  0  0\n" +
            "582583  1  0  0  0  0\n" +
            "615582  1  0  0  0  0\n" +
            "620584  1  0  0  0  0\n" +
            "604584  1  0  0  0  0\n" +
            "603585  2  0  0  0  0\n" +
            "630585  1  0  0  0  0\n" +
            "593586  2  0  0  0  0\n" +
            "605591  1  0  0  0  0\n" +
            "622591  1  0  0  0  0\n" +
            "623593  1  0  0  0  0\n" +
            "614595  1  0  0  0  0\n" +
            "604596  2  0  0  0  0\n" +
            "601597  2  0  0  0  0\n" +
            "598605  2  0  0  0  0\n" +
            "624602  1  0  0  0  0\n" +
            "625602  2  0  0  0  0\n" +
            "631607  1  0  0  0  0\n" +
            "633612  1  0  0  0  0\n" +
            "627614  1  0  0  0  0\n" +
            "614626  2  0  0  0  0\n" +
            "619632  1  0  0  0  0\n" +
            "634619  1  0  0  0  0\n" +
            "635619  1  0  0  0  0\n" +
            "629621  1  0  0  0  0\n" +
            "627624  2  0  0  0  0\n" +
            "626625  1  0  0  0  0\n" +
            "628631  1  0  0  0  0\n" +
            "632628  1  0  0  0  0\n" +
            "637630  1  0  0  0  0\n" +
            "636638  1  0  0  0  0\n" +
            "644636  1  0  0  0  0\n" +
            "640637  1  0  0  0  0\n" +
            "638639  1  0  0  0  0\n" +
            "639641  1  0  0  0  0\n" +
            "643640  1  0  0  0  0\n" +
            "641642  1  0  0  0  0\n" +
            "642643  1  0  0  0  0\n" +
            "644646  1  0  0  0  0\n" +
            "389392  1  6  0  0  0\n" +
            "396389  1  0  0  0  0\n" +
            "431389  1  0  0  0  0\n" +
            "390397  1  0  0  0  0\n" +
            "390436  1  1  0  0  0\n" +
            "390409  1  0  0  0  0\n" +
            "391393  1  0  0  0  0\n" +
            "411391  1  0  0  0  0\n" +
            "439391  1  0  0  0  0\n" +
            "391516  1  1  0  0  0\n" +
            "405392  1  0  0  0  0\n" +
            "393424  1  1  0  0  0\n" +
            "397393  1  0  0  0  0\n" +
            "394398  1  0  0  0  0\n" +
            "428394  2  0  0  0  0\n" +
            "452394  1  0  0  0  0\n" +
            "395396  1  0  0  0  0\n" +
            "413395  1  0  0  0  0\n" +
            "458396  2  0  0  0  0\n" +
            "398403  1  0  0  0  0\n" +
            "398399  1  1  0  0  0\n" +
            "419399  1  0  0  0  0\n" +
            "416400  1  1  0  0  0\n" +
            "403400  1  0  0  0  0\n" +
            "401420  1  0  0  0  0\n" +
            "407401  1  0  0  0  0\n" +
            "402433  1  1  0  0  0\n" +
            "408402  1  0  0  0  0\n" +
            "422402  1  0  0  0  0\n" +
            "459403  2  0  0  0  0\n" +
            "404406  2  0  0  0  0\n" +
            "410404  1  0  0  0  0\n" +
            "423404  1  0  0  0  0\n" +
            "420405  1  6  0  0  0\n" +
            "460405  2  0  0  0  0\n" +
            "406413  1  0  0  0  0\n" +
            "457406  1  0  0  0  0\n" +
            "416407  1  0  0  0  0\n" +
            "461407  2  0  0  0  0\n" +
            "415408  1  0  0  0  0\n" +
            "408477  1  6  0  0  0\n" +
            "409411  1  0  0  0  0\n" +
            "409480  1  6  0  0  0\n" +
            "455410  2  0  0  0  0\n" +
            "410432  1  0  0  0  0\n" +
            "411479  1  1  0  0  0\n" +
            "412418  1  0  0  0  0\n" +
            "424412  1  0  0  0  0\n" +
            "412421  2  0  0  0  0\n" +
            "413441  1  1  0  0  0\n" +
            "414417  1  0  0  0  0\n" +
            "445414  1  0  0  0  0\n" +
            "463414  2  0  0  0  0\n" +
            "426415  1  0  0  0  0\n" +
            "415481  1  6  0  0  0\n" +
            "425416  1  0  0  0  0\n" +
            "430417  1  1  0  0  0\n" +
            "418427  1  0  0  0  0\n" +
            "446418  2  0  0  0  0\n" +
            "430419  1  0  0  0  0\n" +
            "465419  2  0  0  0  0\n" +
            "435420  1  0  0  0  0\n" +
            "421448  1  0  0  0  0\n" +
            "437421  1  0  0  0  0\n" +
            "429422  1  0  0  0  0\n" +
            "433423  1  0  0  0  0\n" +
            "423450  2  0  0  0  0\n" +
            "446425  1  0  0  0  0\n" +
            "448425  2  0  0  0  0\n" +
            "426484  1  1  0  0  0\n" +
            "429426  1  0  0  0  0\n" +
            "427447  1  0  0  0  0\n" +
            "442428  1  0  0  0  0\n" +
            "482428  1  0  0  0  0\n" +
            "429492  1  6  0  0  0\n" +
            "466430  1  0  0  0  0\n" +
            "440431  1  0  0  0  0\n" +
            "431487  1  6  0  0  0\n" +
            "432435  2  0  0  0  0\n" +
            "434443  1  0  0  0  0\n" +
            "438434  1  0  0  0  0\n" +
            "467435  1  0  0  0  0\n" +
            "470436  2  0  0  0  0\n" +
            "488436  1  0  0  0  0\n" +
            "473437  1  0  0  0  0\n" +
            "451438  1  0  0  0  0\n" +
            "462438  2  0  0  0  0\n" +
            "456439  1  0  0  0  0\n" +
            "454440  2  0  0  0  0\n" +
            "468440  1  0  0  0  0\n" +
            "471441  2  0  0  0  0\n" +
            "478441  1  0  0  0  0\n" +
            "489442  1  0  0  0  0\n" +
            "449442  2  0  0  0  0\n" +
            "443452  2  0  0  0  0\n" +
            "449443  1  0  0  0  0\n" +
            "444445  1  0  0  0  0\n" +
            "469444  1  0  0  0  0\n" +
            "451444  2  0  0  0  0\n" +
            "445483  1  6  0  0  0\n" +
            "447472  1  0  0  0  0\n" +
            "453447  2  0  0  0  0\n" +
            "450464  1  0  0  0  0\n" +
            "453454  1  0  0  0  0\n" +
            "486453  1  0  0  0  0\n" +
            "491455  1  0  0  0  0\n" +
            "475455  1  0  0  0  0\n" +
            "474456  2  0  0  0  0\n" +
            "501456  1  0  0  0  0\n" +
            "464457  2  0  0  0  0\n" +
            "476462  1  0  0  0  0\n" +
            "493462  1  0  0  0  0\n" +
            "494464  1  0  0  0  0\n" +
            "485466  1  0  0  0  0\n" +
            "475467  2  0  0  0  0\n" +
            "472468  2  0  0  0  0\n" +
            "469476  2  0  0  0  0\n" +
            "495473  1  0  0  0  0\n" +
            "496473  2  0  0  0  0\n" +
            "502478  1  0  0  0  0\n" +
            "504483  1  0  0  0  0\n" +
            "498485  1  0  0  0  0\n" +
            "485497  2  0  0  0  0\n" +
            "490503  1  0  0  0  0\n" +
            "505490  1  0  0  0  0\n" +
            "506490  1  0  0  0  0\n" +
            "500492  1  0  0  0  0\n" +
            "498495  2  0  0  0  0\n" +
            "497496  1  0  0  0  0\n" +
            "499502  1  0  0  0  0\n" +
            "503499  1  0  0  0  0\n" +
            "508501  1  0  0  0  0\n" +
            "507509  1  0  0  0  0\n" +
            "515507  1  0  0  0  0\n" +
            "511508  1  0  0  0  0\n" +
            "509510  1  0  0  0  0\n" +
            "510512  1  0  0  0  0\n" +
            "514511  1  0  0  0  0\n" +
            "512513  1  0  0  0  0\n" +
            "513514  1  0  0  0  0\n" +
            "515517  1  0  0  0  0\n" +
            "260263  1  6  0  0  0\n" +
            "267260  1  0  0  0  0\n" +
            "302260  1  0  0  0  0\n" +
            "261268  1  0  0  0  0\n" +
            "261307  1  1  0  0  0\n" +
            "261280  1  0  0  0  0\n" +
            "262264  1  0  0  0  0\n" +
            "282262  1  0  0  0  0\n" +
            "310262  1  0  0  0  0\n" +
            "262387  1  1  0  0  0\n" +
            "276263  1  0  0  0  0\n" +
            "264295  1  1  0  0  0\n" +
            "268264  1  0  0  0  0\n" +
            "265269  1  0  0  0  0\n" +
            "299265  2  0  0  0  0\n" +
            "323265  1  0  0  0  0\n" +
            "266267  1  0  0  0  0\n" +
            "284266  1  0  0  0  0\n" +
            "329267  2  0  0  0  0\n" +
            "269274  1  0  0  0  0\n" +
            "269270  1  1  0  0  0\n" +
            "290270  1  0  0  0  0\n" +
            "287271  1  1  0  0  0\n" +
            "274271  1  0  0  0  0\n" +
            "272291  1  0  0  0  0\n" +
            "278272  1  0  0  0  0\n" +
            "273304  1  1  0  0  0\n" +
            "279273  1  0  0  0  0\n" +
            "293273  1  0  0  0  0\n" +
            "330274  2  0  0  0  0\n" +
            "275277  2  0  0  0  0\n" +
            "281275  1  0  0  0  0\n" +
            "294275  1  0  0  0  0\n" +
            "291276  1  6  0  0  0\n" +
            "331276  2  0  0  0  0\n" +
            "277284  1  0  0  0  0\n" +
            "328277  1  0  0  0  0\n" +
            "287278  1  0  0  0  0\n" +
            "332278  2  0  0  0  0\n" +
            "286279  1  0  0  0  0\n" +
            "279348  1  6  0  0  0\n" +
            "280282  1  0  0  0  0\n" +
            "280351  1  6  0  0  0\n" +
            "326281  2  0  0  0  0\n" +
            "281303  1  0  0  0  0\n" +
            "282350  1  1  0  0  0\n" +
            "283289  1  0  0  0  0\n" +
            "295283  1  0  0  0  0\n" +
            "283292  2  0  0  0  0\n" +
            "284312  1  1  0  0  0\n" +
            "285288  1  0  0  0  0\n" +
            "316285  1  0  0  0  0\n" +
            "334285  2  0  0  0  0\n" +
            "297286  1  0  0  0  0\n" +
            "286352  1  6  0  0  0\n" +
            "296287  1  0  0  0  0\n" +
            "301288  1  1  0  0  0\n" +
            "289298  1  0  0  0  0\n" +
            "317289  2  0  0  0  0\n" +
            "301290  1  0  0  0  0\n" +
            "336290  2  0  0  0  0\n" +
            "306291  1  0  0  0  0\n" +
            "292319  1  0  0  0  0\n" +
            "308292  1  0  0  0  0\n" +
            "300293  1  0  0  0  0\n" +
            "304294  1  0  0  0  0\n" +
            "294321  2  0  0  0  0\n" +
            "317296  1  0  0  0  0\n" +
            "319296  2  0  0  0  0\n" +
            "297355  1  1  0  0  0\n" +
            "300297  1  0  0  0  0\n" +
            "298318  1  0  0  0  0\n" +
            "313299  1  0  0  0  0\n" +
            "353299  1  0  0  0  0\n" +
            "300363  1  6  0  0  0\n" +
            "337301  1  0  0  0  0\n" +
            "311302  1  0  0  0  0\n" +
            "302358  1  6  0  0  0\n" +
            "303306  2  0  0  0  0\n" +
            "305314  1  0  0  0  0\n" +
            "309305  1  0  0  0  0\n" +
            "338306  1  0  0  0  0\n" +
            "341307  2  0  0  0  0\n" +
            "359307  1  0  0  0  0\n" +
            "344308  1  0  0  0  0\n" +
            "322309  1  0  0  0  0\n" +
            "333309  2  0  0  0  0\n" +
            "327310  1  0  0  0  0\n" +
            "325311  2  0  0  0  0\n" +
            "339311  1  0  0  0  0\n" +
            "342312  2  0  0  0  0\n" +
            "349312  1  0  0  0  0\n" +
            "360313  1  0  0  0  0\n" +
            "320313  2  0  0  0  0\n" +
            "314323  2  0  0  0  0\n" +
            "320314  1  0  0  0  0\n" +
            "315316  1  0  0  0  0\n" +
            "340315  1  0  0  0  0\n" +
            "322315  2  0  0  0  0\n" +
            "316354  1  6  0  0  0\n" +
            "318343  1  0  0  0  0\n" +
            "324318  2  0  0  0  0\n" +
            "321335  1  0  0  0  0\n" +
            "324325  1  0  0  0  0\n" +
            "357324  1  0  0  0  0\n" +
            "362326  1  0  0  0  0\n" +
            "346326  1  0  0  0  0\n" +
            "345327  2  0  0  0  0\n" +
            "372327  1  0  0  0  0\n" +
            "335328  2  0  0  0  0\n" +
            "347333  1  0  0  0  0\n" +
            "364333  1  0  0  0  0\n" +
            "365335  1  0  0  0  0\n" +
            "356337  1  0  0  0  0\n" +
            "346338  2  0  0  0  0\n" +
            "343339  2  0  0  0  0\n" +
            "340347  2  0  0  0  0\n" +
            "366344  1  0  0  0  0\n" +
            "367344  2  0  0  0  0\n" +
            "373349  1  0  0  0  0\n" +
            "375354  1  0  0  0  0\n" +
            "369356  1  0  0  0  0\n" +
            "356368  2  0  0  0  0\n" +
            "361374  1  0  0  0  0\n" +
            "376361  1  0  0  0  0\n" +
            "377361  1  0  0  0  0\n" +
            "371363  1  0  0  0  0\n" +
            "369366  2  0  0  0  0\n" +
            "368367  1  0  0  0  0\n" +
            "370373  1  0  0  0  0\n" +
            "374370  1  0  0  0  0\n" +
            "379372  1  0  0  0  0\n" +
            "378380  1  0  0  0  0\n" +
            "386378  1  0  0  0  0\n" +
            "382379  1  0  0  0  0\n" +
            "380381  1  0  0  0  0\n" +
            "381383  1  0  0  0  0\n" +
            "385382  1  0  0  0  0\n" +
            "383384  1  0  0  0  0\n" +
            "384385  1  0  0  0  0\n" +
            "386388  1  0  0  0  0\n" +
            "131134  1  6  0  0  0\n" +
            "138131  1  0  0  0  0\n" +
            "173131  1  0  0  0  0\n" +
            "132139  1  0  0  0  0\n" +
            "132178  1  1  0  0  0\n" +
            "132151  1  0  0  0  0\n" +
            "133135  1  0  0  0  0\n" +
            "153133  1  0  0  0  0\n" +
            "181133  1  0  0  0  0\n" +
            "133258  1  1  0  0  0\n" +
            "147134  1  0  0  0  0\n" +
            "135166  1  1  0  0  0\n" +
            "139135  1  0  0  0  0\n" +
            "136140  1  0  0  0  0\n" +
            "170136  2  0  0  0  0\n" +
            "194136  1  0  0  0  0\n" +
            "137138  1  0  0  0  0\n" +
            "155137  1  0  0  0  0\n" +
            "200138  2  0  0  0  0\n" +
            "140145  1  0  0  0  0\n" +
            "140141  1  1  0  0  0\n" +
            "161141  1  0  0  0  0\n" +
            "158142  1  1  0  0  0\n" +
            "145142  1  0  0  0  0\n" +
            "143162  1  0  0  0  0\n" +
            "149143  1  0  0  0  0\n" +
            "144175  1  1  0  0  0\n" +
            "150144  1  0  0  0  0\n" +
            "164144  1  0  0  0  0\n" +
            "201145  2  0  0  0  0\n" +
            "146148  2  0  0  0  0\n" +
            "152146  1  0  0  0  0\n" +
            "165146  1  0  0  0  0\n" +
            "162147  1  6  0  0  0\n" +
            "202147  2  0  0  0  0\n" +
            "148155  1  0  0  0  0\n" +
            "199148  1  0  0  0  0\n" +
            "158149  1  0  0  0  0\n" +
            "203149  2  0  0  0  0\n" +
            "157150  1  0  0  0  0\n" +
            "150219  1  6  0  0  0\n" +
            "151153  1  0  0  0  0\n" +
            "151222  1  6  0  0  0\n" +
            "197152  2  0  0  0  0\n" +
            "152174  1  0  0  0  0\n" +
            "153221  1  1  0  0  0\n" +
            "154160  1  0  0  0  0\n" +
            "166154  1  0  0  0  0\n" +
            "154163  2  0  0  0  0\n" +
            "155183  1  1  0  0  0\n" +
            "156159  1  0  0  0  0\n" +
            "187156  1  0  0  0  0\n" +
            "205156  2  0  0  0  0\n" +
            "168157  1  0  0  0  0\n" +
            "157223  1  6  0  0  0\n" +
            "167158  1  0  0  0  0\n" +
            "172159  1  1  0  0  0\n" +
            "160169  1  0  0  0  0\n" +
            "188160  2  0  0  0  0\n" +
            "172161  1  0  0  0  0\n" +
            "207161  2  0  0  0  0\n" +
            "177162  1  0  0  0  0\n" +
            "163190  1  0  0  0  0\n" +
            "179163  1  0  0  0  0\n" +
            "171164  1  0  0  0  0\n" +
            "175165  1  0  0  0  0\n" +
            "165192  2  0  0  0  0\n" +
            "188167  1  0  0  0  0\n" +
            "190167  2  0  0  0  0\n" +
            "168226  1  1  0  0  0\n" +
            "171168  1  0  0  0  0\n" +
            "169189  1  0  0  0  0\n" +
            "184170  1  0  0  0  0\n" +
            "224170  1  0  0  0  0\n" +
            "171234  1  6  0  0  0\n" +
            "208172  1  0  0  0  0\n" +
            "182173  1  0  0  0  0\n" +
            "173229  1  6  0  0  0\n" +
            "174177  2  0  0  0  0\n" +
            "176185  1  0  0  0  0\n" +
            "180176  1  0  0  0  0\n" +
            "209177  1  0  0  0  0\n" +
            "212178  2  0  0  0  0\n" +
            "230178  1  0  0  0  0\n" +
            "215179  1  0  0  0  0\n" +
            "193180  1  0  0  0  0\n" +
            "204180  2  0  0  0  0\n" +
            "198181  1  0  0  0  0\n" +
            "196182  2  0  0  0  0\n" +
            "210182  1  0  0  0  0\n" +
            "213183  2  0  0  0  0\n" +
            "220183  1  0  0  0  0\n" +
            "231184  1  0  0  0  0\n" +
            "191184  2  0  0  0  0\n" +
            "185194  2  0  0  0  0\n" +
            "191185  1  0  0  0  0\n" +
            "186187  1  0  0  0  0\n" +
            "211186  1  0  0  0  0\n" +
            "193186  2  0  0  0  0\n" +
            "187225  1  6  0  0  0\n" +
            "189214  1  0  0  0  0\n" +
            "195189  2  0  0  0  0\n" +
            "192206  1  0  0  0  0\n" +
            "195196  1  0  0  0  0\n" +
            "228195  1  0  0  0  0\n" +
            "233197  1  0  0  0  0\n" +
            "217197  1  0  0  0  0\n" +
            "216198  2  0  0  0  0\n" +
            "243198  1  0  0  0  0\n" +
            "206199  2  0  0  0  0\n" +
            "218204  1  0  0  0  0\n" +
            "235204  1  0  0  0  0\n" +
            "236206  1  0  0  0  0\n" +
            "227208  1  0  0  0  0\n" +
            "217209  2  0  0  0  0\n" +
            "214210  2  0  0  0  0\n" +
            "211218  2  0  0  0  0\n" +
            "237215  1  0  0  0  0\n" +
            "238215  2  0  0  0  0\n" +
            "244220  1  0  0  0  0\n" +
            "246225  1  0  0  0  0\n" +
            "240227  1  0  0  0  0\n" +
            "227239  2  0  0  0  0\n" +
            "232245  1  0  0  0  0\n" +
            "247232  1  0  0  0  0\n" +
            "248232  1  0  0  0  0\n" +
            "242234  1  0  0  0  0\n" +
            "240237  2  0  0  0  0\n" +
            "239238  1  0  0  0  0\n" +
            "241244  1  0  0  0  0\n" +
            "245241  1  0  0  0  0\n" +
            "250243  1  0  0  0  0\n" +
            "249251  1  0  0  0  0\n" +
            "257249  1  0  0  0  0\n" +
            "253250  1  0  0  0  0\n" +
            "251252  1  0  0  0  0\n" +
            "252254  1  0  0  0  0\n" +
            "256253  1  0  0  0  0\n" +
            "254255  1  0  0  0  0\n" +
            "255256  1  0  0  0  0\n" +
            "257259  1  0  0  0  0\n" +
            "M  STY  2   1 MUL   2 MUL\n" +
            "M  SAL   1 15   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15\n" +
            "M  SAL   1 15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30\n" +
            "M  SAL   1 15  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45\n" +
            "M  SAL   1 15  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60\n" +
            "M  SAL   1 15  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75\n" +
            "M  SAL   1 15  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90\n" +
            "M  SAL   1 15  91  92  93  94  95  96  97  98  99 100 101 102 103 104 105\n" +
            "M  SAL   1 15 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120\n" +
            "M  SAL   1 15 121 122 123 124 125 126 127 128 129 131 132 133 134 135 136\n" +
            "M  SAL   1 15 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151\n" +
            "M  SAL   1 15 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166\n" +
            "M  SAL   1 15 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181\n" +
            "M  SAL   1 15 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196\n" +
            "M  SAL   1 15 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211\n" +
            "M  SAL   1 15 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226\n" +
            "M  SAL   1 15 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241\n" +
            "M  SAL   1 15 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256\n" +
            "M  SAL   1 15 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271\n" +
            "M  SAL   1 15 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286\n" +
            "M  SAL   1 15 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301\n" +
            "M  SAL   1 15 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316\n" +
            "M  SAL   1 15 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331\n" +
            "M  SAL   1 15 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346\n" +
            "M  SAL   1 15 347 348 349 350 351 352 353 354 355 356 357 358 359 360 361\n" +
            "M  SAL   1 15 362 363 364 365 366 367 368 369 370 371 372 373 374 375 376\n" +
            "M  SAL   1 15 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391\n" +
            "M  SAL   1 15 392 393 394 395 396 397 398 399 400 401 402 403 404 405 406\n" +
            "M  SAL   1 15 407 408 409 410 411 412 413 414 415 416 417 418 419 420 421\n" +
            "M  SAL   1 15 422 423 424 425 426 427 428 429 430 431 432 433 434 435 436\n" +
            "M  SAL   1 15 437 438 439 440 441 442 443 444 445 446 447 448 449 450 451\n" +
            "M  SAL   1 15 452 453 454 455 456 457 458 459 460 461 462 463 464 465 466\n" +
            "M  SAL   1 15 467 468 469 470 471 472 473 474 475 476 477 478 479 480 481\n" +
            "M  SAL   1 15 482 483 484 485 486 487 488 489 490 491 492 493 494 495 496\n" +
            "M  SAL   1 15 497 498 499 500 501 502 503 504 505 506 507 508 509 510 511\n" +
            "M  SAL   1 15 512 513 514 515 516 517 518 519 520 521 522 523 524 525 526\n" +
            "M  SAL   1 15 527 528 529 530 531 532 533 534 535 536 537 538 539 540 541\n" +
            "M  SAL   1 15 542 543 544 545 546 547 548 549 550 551 552 553 554 555 556\n" +
            "M  SAL   1 15 557 558 559 560 561 562 563 564 565 566 567 568 569 570 571\n" +
            "M  SAL   1 15 572 573 574 575 576 577 578 579 580 581 582 583 584 585 586\n" +
            "M  SAL   1 15 587 588 589 590 591 592 593 594 595 596 597 598 599 600 601\n" +
            "M  SAL   1 15 602 603 604 605 606 607 608 609 610 611 612 613 614 615 616\n" +
            "M  SAL   1 15 617 618 619 620 621 622 623 624 625 626 627 628 629 630 631\n" +
            "M  SAL   1 15 632 633 634 635 636 637 638 639 640 641 642 643 644 645 646\n" +
            "M  SPA   1 15   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15\n" +
            "M  SPA   1 15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30\n" +
            "M  SPA   1 15  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45\n" +
            "M  SPA   1 15  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60\n" +
            "M  SPA   1 15  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75\n" +
            "M  SPA   1 15  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90\n" +
            "M  SPA   1 15  91  92  93  94  95  96  97  98  99 100 101 102 103 104 105\n" +
            "M  SPA   1 15 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120\n" +
            "M  SPA   1  9 121 122 123 124 125 126 127 128 129\n" +
            "M  SDI   1  4   20.6130  -15.8397   20.6130   -0.4055\n" +
            "M  SDI   1  4   40.3067   -0.4055   40.3067  -15.8397\n" +
            "M  SMT   1 5\n" +
            "M  SAL   2  8 130 647 648 649 650 651 652 653\n" +
            "M  SPA   2  1 130\n" +
            "M  SDI   2  4   46.3949   -9.0415   46.3949   -8.2015\n" +
            "M  SDI   2  4   47.2349   -8.2015   47.2349   -9.0415\n" +
            "M  SMT   2 8\n" +
            "M  END";
}
