package gov.nih.ncats.structureIndexer;

import java.io.IOException;
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
		molSearcher = MolSearcherFactory.create(query).orElse(null);
	}
	
	public int[][] search(Chemical target){
		return VFLib3.subgraphIsomorphism(query, target)

				.findAll(true);
	}
	
	public int[] findMax(Chemical target){
	  
		if(molSearcher == null){
			int[][] search = search(target);
			if(search.length >0) {
				return search[0];
			}
			return new int[0];
		}
		return molSearcher.search(target).orElse(new int[0]);

	}
	
}
