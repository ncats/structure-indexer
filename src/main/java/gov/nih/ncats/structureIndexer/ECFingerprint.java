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
    private boolean encodeWhole=false;
    private Map<Integer,HashMap<Integer,Boolean>> FPbit = new HashMap<Integer,HashMap<Integer,Boolean>>();
    private ArrayList<Integer> onBits = new ArrayList<Integer>();
    private Stack<Atom> atomStack= new Stack<Atom>();

    public ECFingerprint(int nBits, int MAX_LENGTH, int BITS_PER_STRING, boolean encodeWhole){
        this.nBits=nBits;
        this.MAX_LENGTH=MAX_LENGTH;
        this.BITS_PER_STRING=BITS_PER_STRING;
        this.encodeWhole=encodeWhole;
    }

    public long[] getFingerprint(Chemical c){
        return myFingerprint(c,nBits);
    }
    private void assignMap(Chemical c){
        int i=1;
        for(Atom a:c.getAtoms()){
            a.setAtomToAtomMap(i);
            i++;
        }
    }
    private void unAssignMap(Chemical c){
        for(Atom a:c.getAtoms()){
            a.setAtomToAtomMap(0);
        }
    }
    private long[] myFingerprint(Chemical c,int nBits){
        c.makeHydrogensImplicit();
        FPbit = new HashMap<Integer,HashMap<Integer,Boolean>>();
        onBits = new ArrayList<Integer>();
        ArrayList<Atom> atomList = new ArrayList<Atom>();
        ArrayList<Atom> visitedAtoms = new ArrayList<Atom>();
        ArrayList<Atom> cRing = new ArrayList<Atom>();

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
                Collections.sort(cRing, Comparator.comparing(a->asString(a)));
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
        if(encodeWhole) {
            int b=fold(sum);
            fprints[b/64]=fprints[b/64]|(1<<b%64);
        }
        unAssignMap(c);
        return fprints;
    }
    private String asString(Atom ca){
        return ca.getSymbol() + ca.getImplicitHCount();
    }
    private String makeSTR(List<Atom> cList){
        StringBuilder ret=new StringBuilder();
        for(Atom ca: cList){
            ret.append(asString(ca));
        }
        return ret.toString();
    }
    private void addDescriptor(String s){
        for(int i=0;i<BITS_PER_STRING;i++){
            int pos=hashPos(s+"?" +i);
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
    private int fold(long l) {
        return (int)Math.abs(l%(nBits));
    }
}
