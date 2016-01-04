package tripod.chem.indexer;

import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

import chemaxon.formats.MolFormatException;
import chemaxon.sss.search.MolSearch;
import chemaxon.sss.search.SearchException;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class TestVFLib2Isomorphism {

	/*
	 *  static void testIsomorphism () throws Exception {
        automorphism ("CC(C)OC(=O)C1=C(C)NC(N)=C(C1C2=CC(=CC=C2)[N+]([O-])=O)C(=O)OC3CN(C3)C(C4=CC=CC=C4)C5=CC=CC=C5");
    }
	 */
	
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
	public void automorphismSearch() throws SearchException, MolFormatException{
		Molecule m = parseSmiles("CC(C)OC(=O)C1=C(C)NC(N)=C(C1C2=CC(=CC=C2)[N+]([O-])=O)C(=O)OC3CN(C3)C(C4=CC=CC=C4)C5=CC=CC=C5");
		
		 System.out.println ("---- automorphism (MolSearch)  ----");
	        MolSearch ms = new MolSearch ();
	        ms.setQuery(m);
	        ms.setTarget(m);
	        
	        int[][] hits = ms.findAll();
	        
	        
	        HitSetBuilder expected = new HitSetBuilder();
			  
			  int[][] expectedArray= expected.newHit()
								  			.range(0, 42)
								  			.build()
								  		.build();
			  assertHitsMatchExpected(expectedArray, hits);
	      
	       
	}
	
	@Test
	public void automorphism() throws MolFormatException{
		Molecule mol = parseSmiles("CC(C)OC(=O)C1=C(C)NC(N)=C(C1C2=CC(=CC=C2)[N+]([O-])=O)C(=O)OC3CN(C3)C(C4=CC=CC=C4)C5=CC=CC=C5");
		 VFLib2 vf = VFLib2.automorphism (mol);
		 
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

	private int[][] computeIsomorphismFor(String queryStr, String targetStr) throws MolFormatException {
		MolHandler mh = new MolHandler();
		
		mh.setMolecule(queryStr);
		mh.aromatize();
		
		Molecule query = mh.getMolecule();
		
		mh.setMolecule(targetStr);
		mh.aromatize();
		
		Molecule target = mh.getMolecule();
		
		VFLib2 vf = VFLib2.subgraphIsomorphism (query, target);
		
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
