package tripod.chem.indexer;

import java.util.Objects;

import org.junit.rules.ExternalResource;

import gov.nih.ncgc.v3.api.ChemicalProvider;

public class ChemicalProviderRule extends ExternalResource{

	private ChemicalProvider oldProvider;
	
	private final ChemicalProvider providerUnderTest;
	
	ChemicalProviderRule(ChemicalProvider providerUnderTest){
		Objects.requireNonNull(providerUnderTest);
		this.providerUnderTest = providerUnderTest;
	}
	@Override
	protected void after() {
		ChemicalProvider.setDefault(oldProvider);
	}

	@Override
	protected void before() throws Throwable {
		oldProvider = ChemicalProvider.getDefault();
		System.out.println("setting default to " + providerUnderTest);
		ChemicalProvider.setDefault(providerUnderTest);
	}
	public ChemicalProvider getProviderUnderTest() {
		return providerUnderTest;
	}

	
}
