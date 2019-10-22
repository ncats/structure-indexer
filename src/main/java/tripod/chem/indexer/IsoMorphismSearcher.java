package tripod.chem.indexer;

import java.util.Objects;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.search.MolSearcher;
import gov.nih.ncats.molwitch.search.MolSearcherFactory;

public class IsoMorphismSearcher {

	private final Chemical query;
	private MolSearcher molSearcher;

	public IsoMorphismSearcher(Chemical query) {
		Objects.requireNonNull(query);
		this.query = query;
		molSearcher = MolSearcherFactory.create(query);
	}
	
	public int[][] search(Chemical target){
		return VFLib3.subgraphIsomorphism(query, target)
				.findAll(true);
	}
	
	public int[] findMax(Chemical target){

		return molSearcher.search(target).orElse(new int[0]);
		/*
		target.aromatize();
		int[][] hits = search(target);
//		System.out.println("find max hits # hits =" + hits.length);
		if(hits.length ==0){
			return new int[0];
		}
		
		int bestOffset =0;
		int bestLength= hits[0].length;
		for(int i=1; i< hits.length; i++){
			if(hits[i].length > bestLength){
				bestOffset=i;
			}
		}
		return hits[bestOffset];

		 */
	}
	
}
