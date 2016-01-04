// $Id: Base26.java 3478 2009-10-27 21:53:37Z nguyenda $

package tripod.chem.indexer.util;

import java.util.logging.Logger;
import java.util.logging.Level;

import java.util.UUID;
import java.security.SecureRandom;

public class Base26 {
    private static final Logger logger = 
	Logger.getLogger(Base26.class.getName());

    private final static String[] BIGRAM = new String[676]; // 26^2
    private static final int BASE = Character.getNumericValue('A');
    private static final String ALPHA = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    private static final int[][] LUT = new int[26][26];
    static {
	int len = ALPHA.length(), k = 0;

	for (int i = 0; i < len; ++i) {
	    char a = ALPHA.charAt(i);
	    int m = Character.getNumericValue(a) - BASE;
	    for (int j = 0; j < len; ++j) {
		char b = ALPHA.charAt(j);
		int n = Character.getNumericValue(b) - BASE;
		BIGRAM[k] = String.valueOf(a) + b;
		LUT[m][n] = k;
		++k;
	    }
	}
    }

    private static int pos (char c) {
	int p = Character.getNumericValue(c) - BASE;
	if (p < 0 || p >= 26) {
	    throw new IllegalArgumentException 
		("Character " + c + " is not part of the Base26 alphabet!");
	}
	return p;
    }

    public static String encode64 (byte[] b) {
	return encode64 (b, 0);
    }

    public static String encode32 (byte[] b) {
	return encode32 (b, 0);
    }

    public static String encode16 (byte[] b) {
	return encode16 (b, 0);
    }

    public static String encode16 (byte[] b, int o) {
	if (b.length - o < 2) {
	    throw new IllegalArgumentException ("Not enough data");
	}
	int c0 = ((b[o  ] & 0xff) | ((b[o+1] & 0x01) << 8)); // 8 + 1 = 9
	int c1 = (b[o+1] & 0xfe) >> 1; // 7 
	//System.out.println(c0 + " " + c1);
	return BIGRAM[c1] + BIGRAM[c0];
    }

    public static String encode (short x) {
	byte[] b = new byte[2];
	b[1] = (byte)((x & 0xff00) >> 8);
	b[0] = (byte)(x & 0x00ff);
	return encode16 (b);
    }

    public static String encode64 (byte[] b, int o) {
	if (b.length - o < 8) {
	    throw new IllegalArgumentException ("Not enough data");
	}
	int c0 = ((b[o  ] & 0xff) | ((b[o+1] & 0x01) << 8)); // 8 + 1 = 9
	int c1 = ((b[o+1] & 0xfe) | ((b[o+2] & 0x03) << 8)) >> 1; // 7 + 2 = 18
	int c2 = ((b[o+2] & 0xfc) | ((b[o+3] & 0x07) << 8)) >> 2; // 6 + 3 = 27
	int c3 = ((b[o+3] & 0xf8) | ((b[o+4] & 0x0f) << 8)) >> 3; // 5 + 4 = 36
	int c4 = ((b[o+4] & 0xf0) | ((b[o+5] & 0x1f) << 8)) >> 4; // 4 + 5 = 45
	int c5 = ((b[o+5] & 0xe0) | ((b[o+6] & 0x3f) << 8)) >> 5; // 3 + 6 = 54
	int c6 = ((b[o+6] & 0xc0) | ((b[o+7] & 0x7f) << 8)) >> 6; // 2 + 7 = 63
	int c7 = (b[o+7] & 0x80) >> 7; // 64
	/*
	logger.info(c0 + " " + c1 + " " + c2 + " " + c3 + " " 
			 + c4 + " " + c5 + " " + c6 + " " +c7);
	*/
	return  BIGRAM[c7] + BIGRAM[c6] + BIGRAM[c5] 
	    + BIGRAM[c4] + BIGRAM[c3] + BIGRAM[c2] 
			 + BIGRAM[c1] + BIGRAM[c0];
    }

    public static String encode32 (byte[] b, int o) {
	if (b.length - 0 < 4) {
	    throw new IllegalArgumentException ("Not enough data");
	}
	int c0 = ((b[o  ] & 0xff) | ((b[o+1] & 0x01) << 8)); // 8 + 1 = 9
	int c1 = ((b[o+1] & 0xfe) | ((b[o+2] & 0x03) << 8)) >> 1; // 7 + 2 = 18
	int c2 = ((b[o+2] & 0xfc) | ((b[o+3] & 0x07) << 8)) >> 2; // 6 + 3 = 27
	int c3 = ((b[o+3] & 0xf8) >> 3); 
	//logger.info(c0 + " " + c1 + " " + c2 + " " + c3);
	return BIGRAM[c3] + BIGRAM[c2] + BIGRAM[c1] + BIGRAM[c0];
    }

    public static String encode (long x) {
	byte[] b = new byte[8];
	b[7] = (byte)((x & 0xff00000000000000l) >> 56);
	b[6] = (byte)((x & 0x00ff000000000000l) >> 48);
	b[5] = (byte)((x & 0x0000ff0000000000l) >> 40);
	b[4] = (byte)((x & 0x000000ff00000000l) >> 32);
	b[3] = (byte)((x & 0x00000000ff000000l) >> 24);
	b[2] = (byte)((x & 0x0000000000ff0000l) >> 16);
	b[1] = (byte)((x & 0x000000000000ff00l) >>  8);
	b[0] = (byte)( x & 0x00000000000000ffl);
	/*
	logger.info(b[0]+" "+b[1]+" "+b[2]+" "+b[3]+" "+b[4]
		    +" "+b[5]+" "+b[6]+" "+b[7]);
	*/
	return encode64 (b);
    }

    public static String encode32 (long x) {
	if ((x & 0xffffffff00000000l) != 0) {
	    logger.log(Level.WARNING, "Not all bits for "+x 
		       + " is 32-bit encodable!");
	}
	byte[] b = new byte[4];
	b[3] = (byte)((x & 0x00000000ff000000l) >> 24);
	b[2] = (byte)((x & 0x0000000000ff0000l) >> 16);
	b[1] = (byte)((x & 0x000000000000ff00l) >>  8);
	b[0] = (byte)( x & 0x00000000000000ffl);
	return encode32 (b);
    }

    public static String encode (int x) {
	byte[] b = new byte[4];
	b[3] = (byte)((x & 0xff000000) >> 24);
	b[2] = (byte)((x & 0x00ff0000) >> 16);
	b[1] = (byte)((x & 0x0000ff00) >>  8);
	b[0] = (byte)(x & 0xff);
	return encode32 (b);
    }

    public static String encode16 (int x) {
	if ((x & 0xffff0000) != 0) {
	    logger.log(Level.WARNING, "Not all bits for "+x 
		       + " is 16-bit encodable!");
	}
	byte[] b = new byte[2];
	b[1] = (byte)((x & 0x0000ff00) >> 8);
	b[0] = (byte)( x & 0x000000ff);
	//System.out.print(x + " " + b[1] + " " + b[0] + " " );
	return encode16 (b);
    }

    public static int decode16 (String x) {
	if (x.length() < 4) {
	    throw new IllegalArgumentException ("Invalid input for decode16");
	}
	char[] s = x.toCharArray();
	
	int c1 = LUT[pos (s[0])][pos (s[1])] << 1;
	int c0 = LUT[pos (s[2])][pos (s[3])];
	int b0 = c0 & 0xff;
	int b1 = ((c0 & 0xff00) >> 8) | (c1 & 0xff);
	//System.out.println(x + " " + c0 + " " + c1 + " " + (b0 | b1<<8));
	return b0 | b1<<8;
    }

    public static int decode32 (String x) {
	if (x.length() < 8) {
	    throw new IllegalArgumentException ("Invalid input for decode32");
	}
	char[] s = x.toCharArray();
	
	int c3 = LUT[pos (s[0])][pos (s[1])] << 3;
	int c2 = LUT[pos (s[2])][pos (s[3])] << 2;
	int c1 = LUT[pos (s[4])][pos (s[5])] << 1;
	int c0 = LUT[pos (s[6])][pos (s[7])];
	//logger.info(c0 + " " + c1+ " " + c2 + " " + c3);
	int a = (c0 & 0xff);
	int b = ((c0 & 0xff00) >> 8) | (c1 & 0xff);
	int c = ((c1 & 0xff00) >> 8) | (c2 & 0xff);
	int d = ((c2 & 0xff00) >> 8) | (c3 & 0xff);
	//logger.info(x + " " + a + " " + b + " " + c + " " + d);
	return a | (b << 8) | (c << 16) | (d << 24);
    }

    public static long decode64 (String x) {
	if (x.length() < 16) {
	    throw new IllegalArgumentException ("Invalid input for decode64");
	}
	char[] s = x.toCharArray();
	
	int c7 = LUT[pos (s[0])][pos (s[1])] << 7;
	int c6 = LUT[pos (s[2])][pos (s[3])] << 6;
	int c5 = LUT[pos (s[4])][pos (s[5])] << 5;
	int c4 = LUT[pos (s[6])][pos (s[7])] << 4;
	int c3 = LUT[pos (s[8])][pos (s[9])] << 3;
	int c2 = LUT[pos (s[10])][pos (s[11])] << 2;
	int c1 = LUT[pos (s[12])][pos (s[13])] << 1;
	int c0 = LUT[pos (s[14])][pos (s[15])];
	/*
	logger.info(c0+" ["+((c0&0xff00)>>8)+","+(c0&0xff)+"] "+
		    c1+" ["+((c1&0xff00)>>8)+","+(c1&0xff)+"] "+
		    c2+" ["+((c2&0xff00)>>8)+","+(c2&0xff)+"] "+
		    c3+" ["+((c3&0xff00)>>8)+","+(c3&0xff)+"] "+
		    c4+" ["+((c4&0xff00)>>8)+","+(c4&0xff)+"] "+
		    c5+" ["+((c5&0xff00)>>8)+","+(c5&0xff)+"] "+
		    c6+" ["+((c6&0xff00)>>8)+","+(c6&0xff)+"] "+
		    c7+" ["+((c7&0xff00)>>8)+","+(c7&0xff)+"]");
	*/
	long a = (c0 & 0xff);
	long b = ((c0 & 0xff00) >> 8) | (c1 & 0xff);
	long c = ((c1 & 0xff00) >> 8) | (c2 & 0xff);
	long d = ((c2 & 0xff00) >> 8) | (c3 & 0xff);
	long e = ((c3 & 0xff00) >> 8) | (c4 & 0xff);
	long f = ((c4 & 0xff00) >> 8) | (c5 & 0xff);
	long g = ((c5 & 0xff00) >> 8) | (c6 & 0xff);
	long h = ((c6 & 0xff00) >> 8) | (c7 & 0xff);
	/*
	logger.info((byte)a+" "+(byte)b+" "+(byte)c+" "+(byte)d+" "
		    +(byte)e+" "+(byte)f+" "+(byte)g+" "+(byte)h);
	*/
	return a | (b<<8) | (c<<16) | (d<<24) 
	    | (e<<32) | (f<<40) | (g<<48) | (h<<56);
    }

    public static void _main (String[] argv) throws Exception {
	long x = -7908098135001315402l;
	String y = encode (x);
	System.out.println(x + " " + y + " " + decode64 (y));
    }

    public static void main (String[] argv) throws Exception {
	for (int x = -100; x <= 100; ++x) {
	    String y = encode (x);
	    System.out.printf("%1$s %2$4d %3$d\n", y, x, decode32 (y));
	}
	System.out.println
	    (Integer.MAX_VALUE + " " + encode (Integer.MAX_VALUE) + 
	     " " + decode32 (encode (Integer.MAX_VALUE)));
	System.out.println
	    (Integer.MIN_VALUE + " " + encode (Integer.MIN_VALUE) + 
	     " " + decode32 (encode (Integer.MIN_VALUE)));
	System.out.println(Long.MAX_VALUE + " " + encode (Long.MAX_VALUE)
			   + " " + decode64 (encode (Long.MAX_VALUE)));
	System.out.println(Long.MIN_VALUE + " " + encode (Long.MIN_VALUE)
			   + " " + decode64 (encode (Long.MIN_VALUE)));
	SecureRandom rand = new SecureRandom ();
	byte[] b = new byte[8];
	for (int x = 0; x < 100; ++x) {
	    rand.nextBytes(b);
	    UUID uuid = UUID.randomUUID();
	    long z = uuid.getMostSignificantBits();
	    long y = uuid.getLeastSignificantBits();
	    long w = z ^ y;
	    System.out.println(encode64 (b) + " " 
			       + encode (z) 
			       + (z != decode64 (encode (z))?"*":"") 
			       + " "
			       + encode (y) 
			       + (y != decode64 (encode (y)) ? "*":"")
			       + " " 
			       + encode (w)
			       + (w != decode64 (encode (w)) ? "*":""));
	}

	java.util.Set<String> u = new java.util.HashSet<String>();
	for (int x = 0; x <= (1<<16); ++x) {
	    String v = encode16 (x);
	    int y = decode16 (v);
	    if (y != x) {
		System.out.println("** bogus codec for " + x);
	    }
	    /*
	    System.out.println(v + " " + x + " " + y + " " + (y==x)
			       + " " + decode32 (encode (x)));
	    */
	    u.add(v);
	}
	System.out.println(u.size() + " unique codes!");
	long id = 336472450229454154l;
	System.out.println(id + " => " + Base26.encode(id) + " => " 
			   + Base26.decode64(Base26.encode(id)));
    }    

    public static class Decode {
	public static void main (String[] argv) throws Exception {
	    for (String arg : argv) {
		System.out.println(arg + ": " + decode64 (arg));
	    }
	}
    }
}
