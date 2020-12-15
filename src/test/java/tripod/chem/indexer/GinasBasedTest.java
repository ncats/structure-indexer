package tripod.chem.indexer;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import gov.nih.ncats.molwitch.Chemical;
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

		Chemical mol = result.nextElement().mol.get();
		mol.clearAtomMaps();
		assertEquals("COC1=CC=CC=C1", mol.toSmiles(
				new ChemFormat.SmilesFormatWriterSpecification()
						.setKekulization(ChemFormat.KekulizationEncoding.KEKULE)));
		
		assertFalse(result.hasMoreElements());
		
		
	}

	@Test
	public void ensureSubstructureSearchHasBasicSmartsSupportForAnyBond() throws Exception{

		indexer.add("1", "COC1=CC=C(O)C2=C(O)C(C)=C3OC(C)(O)C(=O)C3=C12");
		indexer.add("2", "CC1=C2OC(C)(O)C(=O)C2=C3C4=C(C=C(O)C3=C1O)N5C=CC=CC5=N4");

//		ResultEnumeration result = indexer.substructure("[#7,#8]~C1=c2c3c(OC([#6])(O)C3=O)cc(O)c2=C(O)\\C=C/1");
		ResultEnumeration result = indexer.substructure("O~C1=c2c3c(OC([#6])(O)C3=O)cc(O)c2=C(O)\\C=C/1");

		assertTrue(result.hasMoreElements());
		while(result.hasMoreElements()){
			System.out.println(result.nextElement().mol.get().toSmiles(new ChemFormat.SmilesFormatWriterSpecification()
					.setKekulization(ChemFormat.KekulizationEncoding.KEKULE)
						));
		}
	}
}
