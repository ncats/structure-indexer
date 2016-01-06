package tripod.chem.indexer;

import org.junit.ClassRule;
import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;

import gov.nih.ncgc.v3.spi.jchem.JChemChemicalProvider;

@RunWith(Suite.class)
@SuiteClasses({
	StructureIndexTestSuite.class
})
public class JChemStructureIndexTestSuite {

	@ClassRule
	public static ChemicalProviderRule ChemicalProviderRule = new ChemicalProviderRule(new JChemChemicalProvider());
	
}
