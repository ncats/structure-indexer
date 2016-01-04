// $Id: Base26T.java 3983 2010-02-01 19:38:52Z nguyenda $

package tripod.chem.indexer.util;

import java.util.logging.Logger;
import java.util.logging.Level;

import java.util.UUID;
import java.security.SecureRandom;

public class Base26T {
    private static final Logger logger = 
	Logger.getLogger(Base26T.class.getName());

    // 2^14 = 16384  ~ 26^3 = 17576
    private final static String[] TRIGRAM = new String[17576]; // 26^3
    private static final int BASE = Character.getNumericValue('A');
    private static final String ALPHA = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    private static final int[][][] LUT = new int[26][26][26];
    static {
	int len = ALPHA.length(), p = 0;

	for (int i = 0; i < len; ++i) {
	    char a = ALPHA.charAt(i);
	    int m = Character.getNumericValue(a) - BASE;
	    for (int j = 0; j < len; ++j) {
		char b = ALPHA.charAt(j);
		int n = Character.getNumericValue(b) - BASE;
		for (int l = 0; l < len; ++l) {
		    char c = ALPHA.charAt(l);
		    int k = Character.getNumericValue(c) - BASE;
		    TRIGRAM[p] = String.valueOf(a) + b + c;
		    LUT[m][n][k] = p++;
		}
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
	return encode64 (b, false);
    }

    public static String encode64 (byte[] b, boolean padding) {
	return encode64 (b, 0, padding);
    }

    public static String encode32 (byte[] b) {
	return encode32 (b, false);
    }

    public static String encode32 (byte[] b, boolean padding) {
	return encode32 (b, 0, padding);
    }

    public static String encode16 (byte[] b) {
	return encode16 (b, false);
    }

    public static String encode16 (byte[] b, boolean padding) {
	return encode16 (b, 0, padding);
    }

    public static String encode (short x) {
	return encode (x, false);
    }

    public static String encode (short x, boolean padding) {
	byte[] b = new byte[2];
	b[1] = (byte)((x & 0xff00) >> 8);
	b[0] = (byte)(x & 0x00ff);
	return encode16 (b, padding);
    }


    public static String encode (long x) {
	return encode (x, false);
    }

    public static String encode (long x, boolean padding) {
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
	return encode64 (b, padding);
    }

    public static String encode32 (long x) {
	return encode32 (x, false);
    }

    public static String encode32 (long x, boolean padding) {
	/*
	if ((x & 0xffffffff00000000l) != 0) {
	    logger.log(Level.WARNING, "Not all bits for "+x 
		       + " is 32-bit encodable!");
	}
	*/
	byte[] b = new byte[4];
	b[3] = (byte)((x & 0x00000000ff000000l) >> 24);
	b[2] = (byte)((x & 0x0000000000ff0000l) >> 16);
	b[1] = (byte)((x & 0x000000000000ff00l) >>  8);
	b[0] = (byte)( x & 0x00000000000000ffl);
	return encode32 (b, padding);
    }

    public static String encode (int x) {
	return encode (x, false);
    }
    public static String encode (int x, boolean padding) {
	byte[] b = new byte[4];
	b[3] = (byte)((x & 0xff000000) >> 24);
	b[2] = (byte)((x & 0x00ff0000) >> 16);
	b[1] = (byte)((x & 0x0000ff00) >>  8);
	b[0] = (byte)(x & 0xff);
	return encode32 (b, padding);
    }

    public static String encode16 (int x) {
	return encode16 (x, false);
    }
    public static String encode16 (int x, boolean padding) {
	/*
	if ((x & 0xffff0000) != 0) {
	    logger.log(Level.WARNING, "Not all bits for "+x 
		       + " is 16-bit encodable!");
	}
	*/
	byte[] b = new byte[2];
	b[1] = (byte)((x & 0x0000ff00) >> 8);
	b[0] = (byte)( x & 0x000000ff);
	//System.out.print(x + " " + b[1] + " " + b[0] + " " );
	return encode16 (b, padding);
    }

    public static String encode16 (byte[] b, int o, boolean padding) {
	if (b.length - o < 2) {
	    throw new IllegalArgumentException ("Not enough data");
	}
	int c0 = ((b[o  ] & 0xff) | ((b[o+1] & 0x3f) << 8)); // 8 + 6 = 14
	int c1 = (b[o+1] & 0xc0) >> 6; // 2
	//System.out.println(c0 + " " + c1);
	if (c1 != 0 || padding) {
	    return TRIGRAM[c1] + TRIGRAM[c0];
	}
	return TRIGRAM[c0];
    }

    public static short decode16 (String x) {
	int len = x.length();
	if (len == 0 || (len % 3) != 0) {
	    throw new IllegalArgumentException 
		("Invalid input length for decode16: "+len);
	}
	char[] s = x.toCharArray();

	int c0 = LUT[pos (s[len-3])][pos (s[len-2])][pos (s[len-1])];	
	len -= 3;

	int c1 = 0;
	if (len > 0) {
	    c1 = LUT[pos (s[len-3])][pos (s[len-2])][pos (s[len-1])] << 6;
	}

	int b0 = (short)(c0 & 0xff);
	int b1 = (short)(((c0 & 0xff00) >> 8) | (c1 & 0xff));
	//System.out.println(x + " " + c0 + " " + c1 + " " + (b0 | b1<<8));
	return (short)(b0 | b1<<8);
    }

    public static String encode32 (byte[] b, int o, boolean padding) {
	if (b.length - o < 4) {
	    throw new IllegalArgumentException ("Not enough data");
	}
	int c0 = ((b[o  ] & 0xff) | ((b[o+1] & 0x3f) << 8)); // 8 + 6 = 14
	int c1 = ((b[o+1] & 0xc0) | ((b[o+2] & 0xff) << 8)
		  | ((b[o+3] & 0x0f) << 16)) >> 6; // 2 + 8 + 4 = 28
	int c2 = (b[o+3] & 0xf0); // 4 = 32

	//logger.info(c0 + " " + c1 + " " + c2 + " " + c3);
	StringBuilder sb = new StringBuilder ();
	if (c2 != 0 || padding) {
	    sb.append(TRIGRAM[c2]);
	    sb.append(TRIGRAM[c1]);
	}
	else if (c1 != 0) {
	    sb.append(TRIGRAM[c1]);
	}
	sb.append(TRIGRAM[c0]);

	return sb.toString();
    }

    public static int decode32 (String x) {
	int len = x.length();
	if (len == 0 || (len % 3) != 0) {
	    throw new IllegalArgumentException 
		("Invalid input length for decode32: "+len);
	}
	char[] s = x.toCharArray();

	int c0 = LUT[pos (s[len-3])][pos (s[len-2])][pos (s[len-1])];
	len -= 3;

	int c1 = 0;
	if (len > 0) {
	    c1 = LUT[pos (s[len-3])][pos (s[len-2])][pos (s[len-1])] << 6;
	    len -= 3;
	}
	
	int c2 = 0;
	if (len > 0) {
	    c2 = LUT[pos (s[len-3])][pos (s[len-2])][pos (s[len-1])];
	}

	//logger.info(c0 + " " + c1+ " " + c2 + " " + c3);
	int a = (c0 & 0xff);
	int b = ((c0 & 0xff00) >> 8) | (c1 & 0xff);
	int c = (c1 & 0xff00) >> 8;
	int d = ((c1 & 0xffff0000) >> 16) | (c2 & 0xff);

	//logger.info(x + " " + a + " " + b + " " + c + " " + d);
	return a | (b << 8) | (c << 16) | (d << 24);
    }

    public static String encode64 (byte[] b, int o, boolean padding) {
	if (b.length - o < 8) {
	    throw new IllegalArgumentException ("Not enough data");
	}
	int c0 = ((b[o  ] & 0xff) | ((b[o+1] & 0x3f) << 8)); // 8 + 6 = 14
	int c1 = ((b[o+1] & 0xc0) | ((b[o+2] & 0xff) << 8) 
		  | ((b[o+3] & 0x0f) << 16)) >> 6; // 2 + 8 + 4 = 28
	int c2 = ((b[o+3] & 0xf0) | ((b[o+4] & 0xff) << 8)
		  | ((b[o+5] & 0x03) << 16)) >> 4; // 4 + 8 + 2 = 42
	int c3 = ((b[o+5] & 0xfc) | ((b[o+6] & 0xff) << 8)) >> 2; // 6 + 8 = 56
	int c4 = (b[o+7] & 0xff); // 8 = 64

	//logger.info(c0 + " " + c1 + " " + c2 + " " + c3 + " " + c4);
	StringBuilder sb = new StringBuilder ();
	if (c4 != 0 || padding) {
	    sb.append(TRIGRAM[c4]);
	    sb.append(TRIGRAM[c3]);
	    sb.append(TRIGRAM[c2]);
	    sb.append(TRIGRAM[c1]);
	}
	else if (c3 != 0) {
	    sb.append(TRIGRAM[c3]);
	    sb.append(TRIGRAM[c2]);
	    sb.append(TRIGRAM[c1]);
	}
	else if (c2 != 0) {
	    sb.append(TRIGRAM[c2]);
	    sb.append(TRIGRAM[c1]);
	}
	else if (c1 != 0) {
	    sb.append(TRIGRAM[c1]);
	}
	sb.append(TRIGRAM[c0]);

	return sb.toString();
    }

    public static long decode64 (String x) {
	int len = x.length();
	if (len == 0 || (len % 3) != 0) {
	    throw new IllegalArgumentException 
		("Invalid input length for decode64: "+len);
	}
	char[] s = x.toCharArray();

	int c0 = LUT[pos (s[len-3])][pos (s[len-2])][pos (s[len-1])];
	len -= 3;

	int c1 = 0;
	if (len > 0) {
	    c1 = LUT[pos (s[len-3])][pos (s[len-2])][pos (s[len-1])] << 6;
	    len -= 3;
	}

	int c2 = 0;
	if (len > 0) {
	    c2 = LUT[pos (s[len-3])][pos (s[len-2])][pos (s[len-1])] << 4;
	    len -= 3;
	}
	
	int c3 = 0;
	if (len > 0) {
	    c3 = LUT[pos (s[len-3])][pos (s[len-2])][pos (s[len-1])] << 2;
	    len -= 3;
	}

	int c4 = 0;
	if (len > 0) {
	     c4 = LUT[pos (s[len-3])][pos (s[len-2])][pos (s[len-1])];
	     len -= 3;
	}

 	long b0 = (c0 & 0xff);
	long b1 = ((c0 & 0xff00) >> 8) | (c1 & 0xff);
	long b2 = (c1 & 0xff00) >> 8;
	long b3 = ((c1 & 0xffff0000) >> 16) | (c2 & 0xff);
	long b4 = (c2 & 0xff00) >> 8;
	long b5 = ((c2 & 0xffff0000) >> 16) | (c3 & 0xff);
	long b6 = (c3 & 0xff00) >> 8;
	long b7 = c4 & 0xff;

	return b0 | (b1<<8) | (b2<<16) | (b3<<24) 
	    | (b4<<32)  | (b5<<40) |(b6<<48) | (b7<<56);
    }

    public static String random64 () {
	UUID uuid = UUID.randomUUID();
	long x = uuid.getMostSignificantBits();
	long y = uuid.getLeastSignificantBits();
	return encode (x^y);
    }

    public static void _main (String[] argv) throws Exception {
	long x = -7908098135001315402l;
	String y = encode (x);
	System.out.println(x + " " + y + " " + decode64 (y));
    }

    public static void main (String[] argv) throws Exception {
	for (short x = Short.MIN_VALUE; x < Short.MAX_VALUE; ++x) {
	    String y = encode (x);
	    short z = decode16 (y);
	    if (z != x || Math.abs(x) <= 100) {
		System.out.printf("%1$s %2$4d %3$d\n", y, x, z);
	    }
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
	boolean pad = true;
	for (int x = 0; x < 100; ++x) {
	    rand.nextBytes(b);
	    UUID uuid = UUID.randomUUID();
	    long z = uuid.getMostSignificantBits();
	    long y = uuid.getLeastSignificantBits();
	    long w = z ^ y;
	    System.out.println(encode64 (b, pad) + " " 
			       + encode (z, pad) 
			       + (z != decode64 (encode (z,pad))?"*":"") 
			       + " "
			       + encode (y, pad) 
			       + (y != decode64 (encode (y,pad)) ? "*":"")
			       + " " 
			       + encode (w, pad)
			       + (w != decode64 (encode (w,pad)) ? "*":""));
	}

	java.util.Set<String> u = new java.util.HashSet<String>();
	for (int x = Short.MIN_VALUE; x <= Short.MAX_VALUE; ++x) {
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
	System.out.println(id + " => " + encode (id) + " => " 
			   + decode64 (encode (id)));
    }    

    public static class Decode {
	public static void main (String[] argv) throws Exception {
	    for (String arg : argv) {
		System.out.println(arg + ": " + decode64 (arg));
	    }
	}
    }

    public static class Encode {
	public static void main (String[] argv) throws Exception {
	    for (String arg : argv) {
		System.out.println
		    (arg + ": " + encode (Long.parseLong(arg)));
	    }
	}
    }

    public static class InChICollision {
	public static void main (String[] argv) throws Exception {
	    String[] algo = new String[]{"SHA-1","MD5","SHA-256"};
	    for (String a :algo) {
		System.out.println(a);
	    java.security.MessageDigest md = 
		java.security.MessageDigest.getInstance(a);

	    String major1 = "C63H95ClO21/c1-33(19-42(67)18-17-35(3)64)20-53-55(72)57-39(7)58(79-53)59(73)63(75)31-51(70)37(5)52(85-63)16-14-12-13-15-44-22-43(68)27-61(81-44)29-47(76-11)23-45(82-61)25-50(69)38(6)56(78-41(9)66)36(4)34(2)21-49-28-60(10,74)32-62(84-49)30-48(77-40(8)65)24-46(83-62)26-54(71)80-57/h13,15,17-18,36-39,42-49,51-53,55-59,67-68,70,72-75H,1-3,12,14,16,19-32H2,4-11H3";
	    String minor1 = "/b15-13-,18-17+/t36-,37+,38+,39-,42-,43-,44+,45+,46-,47+,48+,49+,51+,52+,53-,55-,56-,57+,58+,59+,60+,61+,62-,63-/m1/s1";

	    String major2 = "C63H95ClO21/c1-33(19-42(67)18-17-35(3)64)20-53-55(72)57-39(7)58(79-53)59(73)63(75)31-51(70)37(5)52(85-63)16-14-12-13-15-44-22-43(68)27-61(81-44)29-47(76-11)23-45(82-61)25-50(69)38(6)56(78-41(9)66)36(4)34(2)21-49-28-60(10,74)32-62(84-49)30-48(77-40(8)65)24-46(83-62)26-54(71)80-57/h13,15,17-18,36-39,42-49,51-53,55-59,67-68,70,72-75H,1-3,12,14,16,19-32H2,4-11H3";
	    String minor2 = "/b15-13-,18-17+/t36-,37+,38-,39+,42-,43+,44+,45-,46-,47-,48+,49+,51-,52-,53+,55+,56+,57+,58-,59-,60+,61-,62-,63-/m0/s1";

	    // ICXJVZHDZFXYQC-IVINTITMSA-N

	    System.out.println(toHex (md.digest(major1.getBytes())));
	    System.out.println(toHex (md.digest((minor1+minor1).getBytes())));

	    System.out.println("--");
	    System.out.println(toHex (md.digest(major2.getBytes())));
	    System.out.println(toHex (md.digest((minor2+minor2).getBytes())));
	    System.out.println();
	    }
    	}

	public static String toHex (byte[] h) {
	    StringBuilder sb = new StringBuilder ();
	    for (byte b : h) {
		sb.append(String.format("%1$02x", b));
	    }
	    return sb.toString();
	}
    }

    public static class EncodeDouble {
	public static void main (String[] argv) throws Exception {
	    for (String s : argv) {
		double dv = Double.parseDouble(s);
		String enc = encode (Double.doubleToLongBits(dv));
		System.out.println(dv + " <=> " + enc + " <=> " 
				   + Double.longBitsToDouble(decode64 (enc)));
	    }
	}
    }
}
