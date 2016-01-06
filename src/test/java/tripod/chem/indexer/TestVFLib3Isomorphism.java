package tripod.chem.indexer;

import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

import chemaxon.formats.MolFormatException;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;
import gov.nih.ncgc.v3.api.Chemical2;
import gov.nih.ncgc.v3.api.Chemical2Factory;
import gov.nih.ncgc.v3.spi.Chemical2FactoryImpl;
import gov.nih.ncgc.v3.spi.cdk.CdkChemical2FactoryImpl;
import gov.nih.ncgc.v3.spi.jchem.JChemChemical2FactoryImpl;
import gov.nih.ncgc.v3.spi.jchem.JChemChemicalFactory;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TestVFLib3Isomorphism {

	
	private Molecule parseSmiles(String smiles) throws MolFormatException {
		MolHandler mh = new MolHandler(smiles);
		mh.aromatize();

		Molecule m = mh.getMolecule();
		for (int i = 0; i < m.getAtomCount(); ++i) {
			MolAtom a = m.getAtom(i);
			a.setAtomMap(i + 1);
		}

		return m;
	}

	
	@Test
	public void automorphism() throws IOException{
		String smiles = "CC(C)OC(=O)C1=C(C)NC(N)=C(C1C2=CC(=CC=C2)[N+]([O-])=O)C(=O)OC3CN(C3)C(C4=CC=CC=C4)C5=CC=CC=C5";
		
		Chemical2 mol = createChemicalFor(smiles);
		
		 VFLib3 vf = VFLib3.automorphism (mol);
		 
		  int[][] hits = vf.findAll ();
		  
		  
		  
		  HitSetBuilder expected = new HitSetBuilder();
		  
		  int[][] expectedArray= expected.newHit()
		  			.range(0, 42)
		  			.build()
	  			.newHit()
	  				.range(0, 37)
	  				.range(42, 38)
		  			.build()
		  		.newHit()
		  			.range(0, 31)
	  				.range(36, 32)
	  				.range(37, 42)
		  			.build()
		  		.newHit()
		  			.range(0, 31)
	  				.range(36, 32)
	  				.single(37)
	  				.range(42, 38)
		  			.build()
		  		.newHit()
		  			.range(0, 30)
	  				.range(37,42)
	  				.range(31, 36)
		  			.build()
		  		.newHit()
		  			.range(0, 30)
	  				.range(37,42)
	  				.single(31)
	  				.range(36, 32)
		  			.build()
		  		.newHit()
		  			.range(0, 30)
		  			.single(37)
	  				.range(42,38)
	  				.range(31, 36)
		  			.build()
		  		.newHit()
		  			.range(0, 30)
		  			.single(37)
	  				.range(42,38)
	  				.single(31)
	  				.range(36, 32)
		  			.build()
		  			
		  		.newHit()
		  			.range(0, 26)
		  			.range(29,27)
	  				.range(30,42)
		  			.build()
		  		.newHit()
		  		.range(0, 26)
	  			.range(29,27)
	  				.range(30,37)
	  				.range(42,38)
		  			.build()
		  		.newHit()
		  		.range(0, 26)
	  			.range(29,27)
	  				.range(30,31)
	  				.range(36,32)
	  				.range(37, 42)
		  			.build()
		  			
		  		.newHit()
		  		.range(0, 26)
	  			.range(29,27)
	  				.range(30,31)
	  				.range(36,32)
	  				.single(37)
	  				.range(42,38)
		  			.build()
		  			
		  		.newHit()
		  		.range(0, 26)
	  			.range(29,27)
		  			.single(30)
	  				.range(37,42)
	  				.range(31,36)
		  			.build()
		  		.newHit()
		  		.range(0, 26)
	  			.range(29,27)
		  			.single(30)
	  				.range(37,42)
	  				.single(31)
	  				.range(36,32)
		  			.build()
		  		.newHit()
		  		.range(0, 26)
	  			.range(29,27)
		  			.single(30)
		  			.single(37)
	  				.range(42,38)
	  				.range(31,36)
		  			.build()
		  			
		  		.newHit()
		  		.range(0, 26)
	  			.range(29,27)
		  			.single(30)
		  			.single(37)
	  				.range(42,38)
	  				.single(31)
	  				.range(36,32)
		  			.build()
		  		//match 16
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 42)
		  			.build()
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 37)
		  			.range(42, 38)
		  			.build()
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 31)
		  			.range(36, 32)
		  			.range(37, 42)
		  			.build()
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 31)
		  			.range(36, 32)
		  			.single(37)
		  			.range(42, 38)
		  			.build()
		  			
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 30)
		  			.range(37, 42)
		  			.range(31, 36)
		  			.build()
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 30)
		  			.range(37, 42)
		  			.single(31)
		  			.range(36, 32)
		  			.build()
		  			
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 30)
		  			.single(37)
		  			.range(42, 38)
		  			.range(31, 36)
		  			.build()
		  			.newHit()
		  			.range(2, 0)
		  			.range(3, 30)
		  			.single(37)
		  			.range(42, 38)
		  			.single(31)
		  			.range(36, 32)
		  			.build()
		  			
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 26)
		  			.range(29, 27)
		  			.range(30, 42)
		  			.build()
		  			
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 26)
		  			.range(29, 27)
		  			.range(30, 37)
		  			.range(42, 38)
		  			.build()
		  			
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 26)
		  			.range(29, 27)
		  			.range(30, 31)
		  			.range(36, 32)
		  			.range(37, 42)
		  			.build()
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 26)
		  			.range(29, 27)
		  			.range(30, 31)
		  			.range(36, 32)
		  			.single(37)
		  			.range(42, 38)
		  			.build()
		  			
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 26)
		  			.range(29, 27)
		  			.single(30)
		  			.range(37, 42)
		  			.range(31, 36)
		  			.build()
		  			
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 26)
		  			.range(29, 27)
		  			.single(30)
		  			.range(37, 42)
		  			.single(31)
		  			.range(36, 32)
		  			.build()
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 26)
		  			.range(29, 27)
		  			.single(30)
		  			.single(37)
		  			.range(42, 38)
		  			.range(31, 36)
		  			.build()
		  			
		  		.newHit()
		  			.range(2, 0)
		  			.range(3, 26)
		  			.range(29, 27)
		  			.single(30)
		  			.single(37)
		  			.range(42, 38)
		  			.single(31)
		  			.range(36, 32)
		  			.build()
		  			
		  .build();
		  
		  assertHitsMatchExpected(expectedArray, hits);
	}


	private Chemical2 createChemicalFor(String smiles) throws IOException {
	//	Chemical2FactoryImpl impl = new JChemChemical2FactoryImpl();
		
		Chemical2FactoryImpl impl = new CdkChemical2FactoryImpl();
		
		Chemical2 mol = new Chemical2Factory(impl).createFromSmiles(smiles);
		mol.aromatize();
		return mol;
	}

	
	
	
	@Test
	public void isomorphism1() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(COCCC(NCC1=CC=CC=C1)C2=CC=CC=C2)CC(CC3COC3)OCC4=CC=CC=C4", 
				"CC(=C)O[C@@]12CO[C@@H]1C[C@H](O)[C@]3(C)C2[C@H](OC(=O)c4ccccc4)[C@]5(O)C[C@H](OC(=O)[C@H](O)C(NC(=O)c6ccccc6)c7ccccc7)C(=C([C@@H](OC(=O)C)C3=O)C5(C)C)C");
	    
	    
	    int[][] expectedArray= new HitSetBuilder()
	    								.newHit()
	    									.range(26, 29)
	    									.single(31)
	    									.range(33, 35)
	    									.range(37,48)
	    									.single(24)
	    									.range(14, 13)
	    									.range(4, 7)
	    									.range(15, 16)
	    									.range(18, 23)
	    									.build()
	    								.build();
	    
	    assertHitsMatchExpected(expectedArray, hits);
	    
	}
	
	@Test
	public void isomorphism2() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(COCCC(NCC1=CC=CC=C1)C2=CC=CC=C2)CC(CC3COC3)OCC4=CC=CC=C4", 
				"CC(=O)O[C@H]1C[C@H]2OC[C@@]2(OC(=O)C)[C@H]3[C@H](OC(=O)c4ccccc4)[C@]5(O)C[C@H](OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)C(=C([C@@H](O)C(=O)[C@]13C)C5(C)C)C");
	
		
		int[][] expectedArray= new HitSetBuilder()
						.newHit()
							.range(27, 30)
							.single(32)
							.range(34, 36)
							.range(38,49)
							.single(25)
							.range(15, 14)
							.single(9)
							.range(6, 8)
							.range(16, 17)
							.range(19, 24)
							.build()
						.build();
	
		assertHitsMatchExpected(expectedArray, hits);
	}
	
	@Test
	public void isomorphism3() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(C1=CC=CC=C1)C2=CC=CC=C2", 
				"OC(=O)c1ccc(cc1)C2c3ccccc3c4ccccc24");
		
		
		int[][] expectedArray= new HitSetBuilder()
						.newHit()
							.single(9)
							.range(6, 3)
							.range(8, 7)
							.range(10, 15)							
							.build()
						.newHit()
							.single(9)
							.range(6, 3)
							.range(8, 7)
							.single(21)
							.range(16, 20)							
						.build()
							.newHit()
							.range(9, 15)
							.single(21)
							.range(16, 20)							
							.build()
							
						.build();
	
		assertHitsMatchExpected(expectedArray, hits);
	}
	
	@Test
	public void isomorphism4() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(Cc1ccccc1)NCc2ccccc2", 
				"Oc1ccc(cc1)C2CN(Cc3cc(O)ccc23)C(=O)c4ccccc4");
		
		
		
		int[][] expectedArray= new HitSetBuilder()
				.newHit()
					.range(8,7)
					.range(4, 1)
					.range(6, 5)
					.range(9, 13)
					.range(15, 17)
					.build()
				.newHit()
					.range(8,7)
					.range(4, 1)
					.range(6, 5)
					.single(9)
					.single(18)
					.range(20, 25)
					.build()
				 .newHit()
					.range(8,7)
					.single(17)
					.range(11, 13)
					.range(15, 16)
					.single(9)
					.single(18)
					.range(20, 25)
					.build()
				.build();
		
		assertHitsMatchExpected(expectedArray, hits);
	}
	
	@Test
	public void isomorphism5() throws Exception{
		
		int[][] hits = computeIsomorphismFor("O=C(NCc1ccccc1)c2ccccc2", 
				"Fc1ccc(cc1)C(=O)N2CCN3C(=O)c4ccccc4C23c5ccccc5");
		
		int[][] expectedArray= new HitSetBuilder()
				.newHit()
					.range(8,7)
					.single(9)
					.range(21, 20)
					.range(15, 19)
					.range(4, 1)
					.range(6, 5)
					.build()
				.newHit()
					.range(14, 12)
					.range(21, 27)
					.range(15, 20)
				.build()
				.newHit()
					.range(8,7)
					.single(9)
					.range(21, 27)
					.range(4, 1)
					.range(6, 5)
					.build()
				.build();
		
		assertHitsMatchExpected(expectedArray, hits);
	}
	
	@Test
	public void isomorphism6() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(Cc1ccccc1)N2CCNCC2", 
				"CN1CCN2CC(c3ccc(cc3)[N+]([O-])=O)c4ccccc4C2C1");
		
		
		int[][] expectedArray= new HitSetBuilder()
							.newHit()
								.range(5, 6)
								.range(16, 21)
								.range(4, 1)
								.range(23, 22)
								.build()
							.newHit()
								.range(5, 12)
								.range(4, 1)
								.range(23, 22)
								.build()
							.build();
		
		assertHitsMatchExpected(expectedArray, hits);
	}
	
	@Test
	public void isomorphism7() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(Oc1ccccc1)c2ccccc2", 
				"CCNc1ccc-2c(c1)C(Oc3cccc(OC)c-23)c4ccccc4");
		
		int[][] expectedArray= new HitSetBuilder()
				.newHit()
					.range(9, 15)
					.single(18)
					.range(7, 3)
					.single(8)
					.build()
				.newHit()
					.range(9, 15)
					.range(18, 24)
					.build()
				.build();
		
		assertHitsMatchExpected(expectedArray, hits);
		
	}
	
	@Test
	public void isomorphism8() throws Exception{
		
		int[][] hits = computeIsomorphismFor("c1cn(cn1)C(c2ccccc2)c3ccccc3", 
				"c1cn(cn1)C2(c3ccccc3-c4ccccc24)c5ccccc5");
		
		int[][] expectedArray= new HitSetBuilder()
				.newHit()
					.range(0, 11)
					.range(18, 23)
					.build()
				.newHit()
					.range(0, 11)
					.single(17)
					.range(12, 16)
					.build()
				.newHit()
					.range(0, 5)
					.single(17)
					.range(12, 16)
					.range(18, 23)
					.build()
				.build();
		
		assertHitsMatchExpected(expectedArray, hits);
	
	}
	
	@Test
	public void isomorphism9() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(CNCCc1ccccc1)Cc2ccccc2", 
				"CC(Cc1ccccc1)NC2CC3c4c2cccc4CCc5ccccc35");
		
		int[][] expectedArray= new HitSetBuilder()
				.newHit()
					.range(11, 9)
					.range(1, 8)
					.range(12, 18)
					.build()
				.newHit()
					.range(11, 9)
					.range(1, 8)
					.single(12)
					.single(26)
					.range(21, 25)
					.build()
				.build();
		
		assertHitsMatchExpected(expectedArray, hits);
	}
	
	@Test
	public void isomorphism10() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(NSNCc1ccccc1)OCc2ccccc2", 
				"O=C1N(Cc2ccccc2)S(=O)(=O)N(COCc3ccccc3)c4ccccc14");
		
		int[][] expectedArray= new HitSetBuilder()
				.newHit()
					.range(14, 13)
					.single(10)
					.range(2, 9)
					.range(15, 22)
					.build()
				.newHit()
					.range(14, 13)
					.single(10)
					.range(2, 1)
					.single(28)
					.range(23, 27)
					.range(15, 22)
					.build()
				.build();
		
		assertHitsMatchExpected(expectedArray, hits);
	
	}
	
	@Test
	public void isomorphism11() throws Exception{
		
		int[][] hits = computeIsomorphismFor("O=C(CCCCC1CCCCC1)NCc2ccccc2", 
				"CC12CC(O)C3C(CCC4=CC(=O)C=CC34C)C1CC[C@]2(O)C(O)C(=O)NCc5ccccc5");
		
		
		int[][] expectedArray= new HitSetBuilder()
						.newHit()
							.range(25, 24)
							.single(22)
							.single(20)
							.single(1)
							.single(17)
							.range(6, 5)
							.single(15)
							.range(9, 7)
							.range(26, 33)
							.build()
						.newHit()
							.range(25, 24)
							.single(22)
							.single(20)
							.range(19, 17)
							.range(1, 3)
							.range(5, 6)
							.range(26, 33)
							.build()
				.build();
		
		assertHitsMatchExpected(expectedArray, hits);
	
	}
	
	@Test
	public void isomorphism12() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(Cc1ccccc1)Nc2ccccc2", 
				"O=C(Nc1ccccc1)C2C(=O)N3c4c2cccc4Sc5ccccc35");
		
		int[][] expectedArray= new HitSetBuilder()
					.newHit()
						.range(10, 9)
						.range(14, 13)
						.range(18, 15)
						.single(12)
						.single(25)
						.range(20, 24)
						.build()
					.newHit()
						.single(1)
						.single(9)
						.range(14, 13)
						.range(18, 15)
						.range(2, 8)
						.build()
				.build();
		
		assertHitsMatchExpected(expectedArray, hits);
	}

	@Test
	public void isomorphism13() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(c1ccccc1)[n+]2ccn3CCC(Cc23)c4ccoc4", 
				"[Cl-].C(c1ccccc1)n2cc[n+]3C4CC(C(c23)c5ccccc45)(c6ccoc6)c7ccoc7");
		
		int[][] expectedArray= new HitSetBuilder()
					.newHit()
						.range(1, 16)
						.range(23, 27)
						.build()
					.newHit()
						.range(1, 16)
						.range(28, 32)
						.build()
				.build();
		assertHitsMatchExpected(expectedArray, hits);	
	}
	
	@Test
	public void isomorphism14() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(CC1CCCC1)OC2COCC3C(CCCC23)C4CCC=CC4", 
				"C[C@H]1[C@@H](O)[C@@]2(O)OC[C@@]34[C@@H](C[C@H]5C(=CC(=O)[C@@H](O)[C@]5(C)[C@@H]23)C)OC(=O)[C@H](OC(=O)CC6(O)C(C)(C)CCC6(C)C)[C@@H]14");
		
		int[][] expectedArray= new HitSetBuilder()
					.newHit()
					.single(27)
						.range(29, 30)
						.single(32)
						.range(35, 37)
						.range(26, 25)
						.range(23, 22)
						.range(9, 8)
						.single(20)
						.single(4)
						.range(2, 1)
						.single(40)
						.single(18)
						.single(16)
						.range(14, 11)
						.build()
					
				.build();
		
		assertHitsMatchExpected(expectedArray, hits);
		
		
	
	}
	
	@Test
	public void isomorphism15() throws Exception{
		
		int[][] hits = computeIsomorphismFor("C(CC1CCCC1)OC2COCC3C(CCCC23)C4CCC=CC4", 
				"C[C@H]1[C@@H](O)[C@@]2(O)OC[C@@]34[C@@H](C[C@H]5C(=CC(=O)[C@@H](O)[C@]5(C)[C@@H]23)C)OC(=O)[C@H](OC(=O)CC6(O)CCCC6)[C@@H]14"
				);
		
		int[][] expectedArray= new HitSetBuilder()
				.newHit()
				.single(27)
					.range(29, 30)
					.range(32, 35)
					.range(26, 25)
					.range(23, 22)
					.range(9, 8)
					.single(20)
					.single(4)
					.range(2, 1)
					.single(36)
					.single(18)
					.single(16)
					.range(14, 11)
					.build()
				
			.build();
	
	assertHitsMatchExpected(expectedArray, hits);
	}
	
	private int[][] computeIsomorphismFor(String queryStr, String targetStr) throws IOException {
		
		VFLib3 vf = VFLib3.subgraphIsomorphism (createChemicalFor(queryStr), createChemicalFor(targetStr));
		
	    int[][] hits = vf.findAll (true);
		return hits;
	}
	
	
	private void assertHitsMatchExpected(int[][] expectedArray, int[][] hits) throws ArrayComparisonFailure {
		assertEquals(expectedArray.length, hits.length);
		 
		  for(int i=0; i<expectedArray.length; i++){
			  assertArrayEquals(toString(i, hits[i]), expectedArray[i], hits[i]);
		  }
	}
	private static String toString(int[][] hits){
		StringBuilder builder = new StringBuilder();
		for(int i=0; i< hits.length; i++){
			builder.append(toString(i, hits[i]));
			builder.append("\n\n");
		}
		return builder.toString();
	}
	private static String toString(int j, int[] hit){
		StringBuilder builder = new StringBuilder();
		builder.append ("Matched " + j + ":");
         for (int i = 0; i < hit.length; ++i) {
            builder.append(" ").append(i + 1).append(":").append(hit[i] + 1);
            
         }
         builder.append("\n");
        builder.append(Arrays.toString(hit));
        return builder.toString();
	}
	
	private static class HitSetBuilder{
		List<HitBuilder> builders = new ArrayList<>();
		
		public HitBuilder newHit(){
			HitBuilder builder = new HitBuilder(this);
			builders.add(builder);
			return builder;
		}
		
		public HitSetBuilder add(HitBuilder duplicate){
			builders.add(duplicate);
			return this;
		}
		
		public int[][] build(){
			int[][] ret = new int[builders.size()][];
			for(int j=0; j< ret.length; j++){
				HitBuilder builder = builders.get(j);
				int[] array = new int[builder.hits.size()];
				int i=0;
				for(Integer hit : builder.hits){
					array[i++] = hit.intValue();
				}
				ret[j] = array;
			}
			return ret;
		}
	}
	
	private static class HitBuilder{
		private final List<Integer> hits = new ArrayList<>();
		private HitSetBuilder parent;
		
		HitBuilder(HitSetBuilder parent){
			this.parent = parent;
		}
		public HitBuilder range(int from, int to){
			if(to > from){
				for(int i=from; i<=to; i++){
					hits.add(i);
				}
			}else{
				for(int i=from; i>=to; i--){
					hits.add(i);
				}
			}
			return this;
		}
		public HitBuilder single(int i){			
			hits.add(i);
			return this;
		}
		public HitSetBuilder build(){
			return parent;
			/*
			int[] ret = new int[hits.size()];
			int i=0;
			for(Integer hit : hits){
				ret[i++] = hit.intValue();
			}
			return ret;
			*/
		}
		
	}
	
	
	
}
