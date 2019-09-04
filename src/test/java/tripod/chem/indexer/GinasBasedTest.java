package tripod.chem.indexer;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import gov.nih.ncats.molwitch.io.ChemFormat;
import tripod.chem.indexer.StructureIndexer.ResultEnumeration;
public class GinasBasedTest extends AbstractStructureIndexerTest{

	/*
	 *  JsonNode form1 = makeChemicalSubstanceJSON("CCCC");
            JsonNode form2 = makeChemicalSubstanceJSON("COC1=CC=CC=C1");
            
            ensurePass( api.submitSubstance(form1));
            ensurePass( api.submitSubstance(form2));
            
            
            SubstanceSearcher searcher = new RestSubstanceSearcher(session);
            SearchResult result=searcher.substructure("C1=CC=CC=C1");
	 */
	
	@Test
	public void findSubstructure() throws Exception{
		indexer.add("id1", "CCCC");
		indexer.add("id2", "COC1=CC=CC=C1");
		
		ResultEnumeration result = indexer.substructure("C1=CC=CC=C1");
		assertTrue(result.hasMoreElements());
		
		assertEquals("COC1=CC=CC=C1", result.nextElement().mol.toSmiles(
				new ChemFormat.SmilesFormatWriterSpecification()
						.setKekulization(ChemFormat.KekulizationEncoding.KEKULE)));
		
		assertFalse(result.hasMoreElements());
		
		
	}
}
