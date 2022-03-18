package gov.nih.ncats.structureIndexer;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.stream.Collectors;

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Chemical;



public class ECFingerprint{
	private int nBits;
	private int MAX_LENGTH=8;
	private int BITS_PER_STRING=1;
	
	private boolean wholeBit=true;
	
	private Map<Integer,HashMap<Integer,Boolean>> FPbit = new HashMap<Integer,HashMap<Integer,Boolean>>();
	private ArrayList<Integer> onBits = new ArrayList<Integer>();
	private Stack<Atom> atomStack= new Stack<Atom>();
	private Stack<Integer> bondStack= new Stack<Integer>();
	private Stack<String> descriptor= new Stack<String>();
	
	public ECFingerprint(int nBits, int MAX_LENGTH, int BITS_PER_STRING){
		this.nBits=nBits;
		this.MAX_LENGTH=MAX_LENGTH;
		this.BITS_PER_STRING=BITS_PER_STRING;
	}
	
	public long[] getFingerprint(Chemical c){
		return myFingerprint(c,nBits);
	}
	
	private static void assignMap(Chemical c){
		int i=1;
		
		for(Atom a:c.getAtoms()){
			a.setAtomToAtomMap(i);
			i++;
		}
	}
	public static void unAssignMap(Chemical c){
		for(Atom a:c.getAtoms()){
			a.setAtomToAtomMap(0);
		}
	}
	//MY STUFF ... kinda dumb
	private long[] myFingerprint(Chemical c,int nBits){
		c.makeHydrogensImplicit();
		FPbit = new HashMap<Integer,HashMap<Integer,Boolean>>();
		onBits = new ArrayList<Integer>();
		List<Atom> atomList = new ArrayList<Atom>();
		List<Atom> visitedAtoms = new ArrayList<Atom>();
		List<Atom> cRing = new ArrayList<Atom>();
				
		assignMap(c);
		for(Atom atom : c.getAtoms()){
			atomList = new ArrayList<Atom>();
			visitedAtoms = new ArrayList<Atom>();
			atomList.add(atom);
			String desc=atom.getSymbol();
			int depth=0;
			this.addDescriptor(desc);
			while(depth<MAX_LENGTH){
				visitedAtoms.addAll(atomList);	
				for(Atom gAtom : atomList){
					for(Atom nAtom : gAtom.getNeighbors()){
						if(!visitedAtoms.contains(nAtom))
							cRing.add(nAtom);
					}
				}
				if(cRing.size()<=0)break;
				Collections.sort(cRing, new Comparator<Atom>(){
					@Override
					public int compare(Atom arg0, Atom arg1) {
						return asString(arg0).compareTo(asString(arg1));
					}});
				desc+="(" +makeSTR(cRing) + ")";
				atomList=cRing;
				cRing = new ArrayList<Atom>();	
				depth++;
				this.addDescriptor(desc);
			}	
		}
	
		long[] fprints= new long[nBits/64];
		long sum=0;
		for(int i:onBits){
			sum+=i;
			fprints[i/64]=fprints[i/64]|(1<<i%64);
		}
		if(wholeBit){
			int wBit=fold(sum);
			fprints[wBit/64]=fprints[wBit/64]|(1<<wBit%64);
		}
		
		unAssignMap(c);
		return fprints;
	}
	private static String asString(Atom ca){
		return ca.getSymbol() + ca.getImplicitHCount();
	}
	private static String makeSTR(List<Atom> cList){
		String ret="";
		for(Atom ca: cList){
			ret+=asString(ca);
		}
		return ret;
	}
	private void addDescriptor(String s){
		for(int i=0;i<BITS_PER_STRING;i++){
			int pos=hashPos(s + "?" + i);
			this.addReferences(pos);
		}
	}
	private void addReferences(int pos){
		onBits.add(pos);
		HashMap<Integer,Boolean> cAtoms=FPbit.get(pos);
		if(cAtoms==null){
			cAtoms= new HashMap<Integer, Boolean>();
			FPbit.put(pos,cAtoms);
		}
		for(Atom a2:atomStack){
			cAtoms.put(a2.getAtomToAtomMap().orElse(0),true);
		}
	}
	private int hashPos(String s){
		return fold(s.hashCode());
	}
	private int fold(long n){
		return Math.abs(
			     (int)(n%(nBits))
			       );
	}
	
	public int getNumberBits() {
		return this.nBits;
	}
	
	public void setNumberBits(int nBits) {
		this.nBits=nBits;
	}
}
