package tripod.chem.indexer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import gov.nih.ncats.chemkit.api.writer.StandardChemFormats;
import gov.nih.ncats.chemkit.api.writer.WriterOptionsBuilder;
import tripod.chem.indexer.StructureIndexer.ResultEnumeration;

import static org.junit.Assert.*;
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
		
		assertEquals("COC1=CC=CC=C1", result.nextElement().mol.formatToString(StandardChemFormats.SMILES, new WriterOptionsBuilder().kekulize().build()));
		
		assertFalse(result.hasMoreElements());
		
		
	}
}
