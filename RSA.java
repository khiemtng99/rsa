package mm1;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.Random;
import java.util.StringJoiner;


public class RSA {	
	private BigInteger p;
	private BigInteger q;
	private BigInteger n; 
	private BigInteger e;
	private BigInteger d;  
	private int bitLength;
	final static int RADIX = 16;
	
	//--------------------constructor---------------------//
	private RSA(int bitLength) {
		this.bitLength = bitLength;
		System.out.println("Generating prime numbers...");
		final long startTime1 = System.nanoTime();
		generatePrimeNumbers();
		final long endTime1 = System.nanoTime();
		System.out.println("Completed! Time for generating 2 prime numbers with " + bitLength + "-bit long is: " + (endTime1 - startTime1) + "ns (~" + (float)(endTime1 - startTime1)/1000000 + "ms)");
		System.out.println("Creating keys...");
		final long startTime2 = System.nanoTime();
		createKeys();
		final long endTime2 = System.nanoTime();
		System.out.println("Completed! Time for creating keys is: " + (endTime2 - startTime2) + "ns (~" + (float)(endTime2 - startTime2)/1000000 + "ms)");
	}
	
	//--------------------get value of p, q, n, d, e---------------------//
	private BigInteger getp() {
		return (p) ;
	}

	private BigInteger getq() {
		return (q) ;
	}

	private BigInteger getn() {
		return (n) ;
	}
	
	private BigInteger getd() {
		return (d) ;
	}
	
	private BigInteger gete() {
		return (e) ;
	}
	
	
	//---------------compare 2 big integers-----------//
	private int compare(BigInteger so1, BigInteger so2) {	
		String str1 = so1.toString();
		String str2 = so2.toString();
		int mark = 1;
		if (str1.contains("-") && str2.contains("-")) {
		    mark = -1;
		}
		if (str1.contains("-") && !str2.contains("-"))
		    return -1;
		else if (!str1.contains("-") && str2.contains("-"))
		    return 1;
        // Calculate lengths of both string
        int n1 = str1.length(), n2 = str2.length();
        if(n1 < n2)
            return -1;
        else if(n1 > n2)
            return 1;
        else {
        	for(int i = 0; i < n1; i++) {
                if (str1.charAt(i) < str2.charAt(i))
                    return -1*mark;
                else if (str1.charAt(i) > str2.charAt(i))
                    return 1*mark;
        	}
        }
        return 0;
    }
	
	// ------------------------Addition of 2 Big Integers ------------------------//
	private BigInteger add(BigInteger op1, BigInteger op2) {
		// Before proceeding further, make sure length
        // of str2 is larger.
        
        String sOp1 = op1.toString();
        String sOp2 = op2.toString();
        
        int signum = 0;
        
        if (sOp2.equals("0"))                           //op2 = 0
            return op1;
        if (sOp1.equals("0"))                           //op2 = 0
            return op2;
        
        if (sOp1.contains("-") && !sOp2.contains("-")) {  // op1 <0, op2 > 0 
            BigInteger num1 = new BigInteger(sOp1.substring(1));
            BigInteger num3 = subtract(op2,num1);
            return new BigInteger(num3.toString());
        }
        else if (!sOp1.contains("-") && sOp2.contains("-")) {   //op1 > 0, op2 < 0
            BigInteger num2 = new BigInteger(sOp2.substring(1));
            BigInteger num3 = subtract(op1,num2);
            return num3;
        }
        else if (sOp1.contains("-") && sOp2.contains("-")) {        // op1, op2 both <0, add as normal
            String sTemp = sOp1.substring(1);
            sOp1 = sTemp;
            sTemp = sOp2.substring(1);
            sOp2 = sTemp;
            signum = -1;
        } 
        else if (!sOp1.contains("-") && !sOp2.contains("-")) {      // op1, op2 both >0, subtract as normal
            signum = 1;
        }
        
        //////////////////////////////
        
        if (sOp1.length() > sOp2.length()) {
            String sTemp = sOp1;
            sOp1 = sOp2;
            sOp2 = sTemp;
        }
     
        // Take an empty String for storing result
        String sResult = "";
        
        // Calculate length of both String
        int opLength1 = sOp1.length(), opLength2 = sOp2.length(); // opLength1 < opLength2
     
        // Reverse both of Strings
        sOp1 = new StringBuilder(sOp1).reverse().toString();
        sOp2 = new StringBuilder(sOp2).reverse().toString();
     
        int carry = 0;
        for (int i = 0; i < opLength1; i++) {
            // Do school mathematics, compute sum of
            // current digits and carry
            int sum = ((int)(sOp1.charAt(i) - '0') +
                        (int)(sOp2.charAt(i) - '0') + carry);
            sResult += (char)(sum % 10 + '0');
     
            // Calculate carry for next step
            carry = sum / 10;
        }
     
        // Add remaining digits of larger number
        for (int i = opLength1; i < opLength2; i++) {
            int sum = ((int)(sOp2.charAt(i) - '0') + carry);
            sResult += (char)(sum % 10 + '0');
            carry = sum / 10;
        }
     
        // Add remaining carry
        if (carry > 0)
            sResult += (char)(carry + '0');
     
        // reverse resultant String
        String sResultReversed = new StringBuilder(sResult).reverse().toString();
        if (signum == -1)
            sResultReversed = "-"+sResultReversed;
        BigInteger result = new BigInteger(sResultReversed);
        return result;
	}
 
	// ------------------------Subtraction of 2 Big Integers ------------------------//
	private BigInteger subtract(BigInteger op1, BigInteger op2) {
		String sOp1 = op1.toString();
        String sOp2 = op2.toString();
        String big ="", little ="";
        int signum = 0;
        
        if (sOp2.equals("0"))                           //op2 = 0
            return op1;
        if (sOp1.equals("0") && sOp2.contains("-"))     //op1 = 0, op2 < 0
            return new BigInteger(sOp2.substring(1));
        if (sOp1.equals("0") && !sOp2.contains("-"))    //op1 = 0, op2 > 0
            return new BigInteger("-"+sOp2);
        if (sOp1.equals(sOp2))                          //op1 = op2
            return BigInteger.ZERO;
        
        if (sOp1.contains("-") && !sOp2.contains("-")) {  // op1 <0, op2 > 0 
            BigInteger num1 = new BigInteger(sOp1.substring(1));
            BigInteger num3 = add(num1,op2);
            return new BigInteger("-"+num3.toString());
        }
        else if (!sOp1.contains("-") && sOp2.contains("-")) {   //op1 > 0, op2 < 0
            BigInteger num2 = new BigInteger(sOp2.substring(1));
            BigInteger num3 = add(op1,num2);
            return num3;
        }
        else if (sOp1.contains("-") && sOp2.contains("-")) {        // op1, op2 both <0, subtract as normal
            
            if (compare(op1,op2) == -1)     { //op1 < op2 => hieu = - ( abs(op1) - abs(op2) )
                big = sOp1.substring(1);
                little = sOp2.substring(1);
                signum = -1;
            } 
            else { //op1 > op2 => hieu = ( abs(op1) - abs(op2) )
                big = sOp2.substring(1);
                little = sOp1.substring(1);
                signum = 1;
            }
        } 
        else if (!sOp1.contains("-") && !sOp2.contains("-")) {      // op1, op2 both >0, subtract as normal
            if (compare(op1,op2) == -1) {       // op1 < op2 => hieu = - ( abs(op2) - abs(op1) )
                big = sOp2;
                little = sOp1;
                signum = -1;
            }
            else {  //op1 > op2 => hieu =  abs(op1) - abs(op2) 
                big = sOp1;
                little = sOp2;
                signum = 1;
            }
        }
  
 
        // Take an empty string for storing result
        String sResult = "";
 
        // Calculate length of both string
        int numLength1 = big.length(), numLength2 = little.length();
 
        // Reverse both of strings
        big = new StringBuilder(big).reverse().toString();
        little = new StringBuilder(little).reverse().toString();
 
        int carry = 0;
 
        // Run loop till small string length
        // and subtract digit of str1 to str2
        for (int i = 0; i < numLength2; i++) {
            // Do school mathematics, compute difference of
            // current digits
            int sub
                = ((int)(big.charAt(i) - '0') - (int)(little.charAt(i) - '0') - carry);
 
            // If subtraction is less then zero
            // we add then we add 10 into sub and
            // take carry as 1 for calculating next step
            if (sub < 0) {
                sub = sub + 10;
                carry = 1;
            }
            else
                carry = 0;
 
            sResult += (char)(sub + '0');
        }
        
        
        // subtract remaining digits of larger number
        for (int i = numLength2; i < numLength1; i++) {
            int sub = ((int)(big.charAt(i) - '0') - carry);
 
            if (sub < 0) {
                sub = sub + 10;
                carry = 1;
            }
            else
                carry = 0;
 
            sResult += (char)(sub + '0');
        }
        
        
        // reverse resultant string
        String sResultReversed = new StringBuilder(sResult).reverse().toString();
        
        BigInteger resultNum;
        if (signum == 1)  {
            resultNum = new BigInteger(sResultReversed);
            return resultNum; 
        }
        else if (signum == -1 ) {
            resultNum =  new BigInteger("-"+sResultReversed);
            return resultNum; 
        }
        return BigInteger.ZERO;
    }
    
    
    // ------------------------Multiplying 2 Big Integers ------------------------//
	private BigInteger multiply(BigInteger op1, BigInteger op2) {
		String sOp1 = op1.toString();
		String sOp2 = op2.toString();
		String mark = "+";
		// check if negative
		if ((sOp1.charAt(0) == '-' || sOp2.charAt(0) == '-') &&
        (sOp1.charAt(0) != '-' || sOp2.charAt(0) != '-'))
            mark = "-";
 
        if (sOp1.charAt(0) == '-')
            sOp1 = sOp1.substring(1);
       
        if (sOp2.charAt(0) == '-')
            sOp2 = sOp2.substring(1);
        int len1 = sOp1.length();
        int len2 = sOp2.length();
     
        // will keep the result number in vector
        // in reverse order
        int result[] = new int[len1 + len2];
     
        // Below two indexes are used to
        // find positions in result.
        int i_n1 = 0;
        int i_n2 = 0;
         
        // Go from right to left in num1
        for (int i = len1 - 1; i >= 0; i--)
        {
            int carry = 0;
            int n1 = sOp1.charAt(i) - '0';
     
            // To shift position to left after every
            // multipliccharAtion of a digit in num2
            i_n2 = 0;
             
            // Go from right to left in num2            
            for (int j = len2 - 1; j >= 0; j--)
            {
                // Take current digit of second number
                int n2 = sOp2.charAt(j) - '0';
     
                // Multiply with current digit of first number
                // and add result to previously stored result
                // charAt current position.
                int sum = n1 * n2 + result[i_n1 + i_n2] + carry;
     
                // Carry for next itercharAtion
                carry = sum / 10;
     
                // Store result
                result[i_n1 + i_n2] = sum % 10;
     
                i_n2++;
            }
     
            // store carry in next cell
            if (carry > 0)
                result[i_n1 + i_n2] += carry;
     
            // To shift position to left after every
            // multipliccharAtion of a digit in num1.
            i_n1++;
        }
     
        // ignore '0's from the right
        int i = result.length - 1;
        while (i >= 0 && result[i] == 0)
        i--;

     
        // genercharAte the result String
        String str3 = "";
         
        while (i >= 0)
            str3 += (result[i--]);
        // when 1 * 0 = 0
        if (str3.length() == 0)
            return BigInteger.ZERO;
     
        BigInteger ketqua = new BigInteger(mark + str3);
	    return ketqua;
    }
    
    // ------------------------Dividing 2 Big Integers ------------------------//
//	private BigInteger divide(BigInteger num1, BigInteger num2) {
//    	/*-----------------------------------------------------------*/
//    	
//    	return num1.divide(num2);
//    }
    
    
    // ------------------------generate a prime number given the bitlength------------------------//
	private BigInteger genRandBigPrimeNum(int bitLength) {
    	//number of iteration rounds to verify Miller Rabin
    	int rounds = 0;
    	if (bitLength < 512) {
            rounds = 15;
        } else if (bitLength < 768) {
            rounds = 8;
        } else if (bitLength < 1024) {
            rounds = 4;
        } else {
            rounds = 2;
        }
    	
    	int[] prePrimeList = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                 31, 37, 41, 43, 47, 53, 59, 61, 67, 
                 71, 73, 79, 83, 89, 97, 101, 103, 
                 107, 109, 113, 127, 131, 137, 139, 
                 149, 151, 157, 163, 167, 173, 179, 
                 181, 191, 193, 197, 199, 211, 223,
                 227, 229, 233, 239, 241, 251, 257,
                 263, 269, 271, 277, 281, 283, 293,
                 307, 311, 313, 317, 331, 337, 347, 349};
    	
    	
    	//---------generate prime candidate and verify primary property------------//
    	//int byteLength = (int)(((long)bitLength+7)/8); // calculte number of bytes to generate random bytes, avoid overflow
    	//byte[] randomBits = new byte[byteLength]; //array to store random bytes
    	
    	boolean preTest;
    	BigInteger divisor = BigInteger.ZERO;
    	while (true) {
    		//generate a SecureRandom
    		Random random = new Random();
    		BigInteger p = new BigInteger(bitLength, random);
    		
            //-------------1st/ low-level pre-test with prime number list-----------//
            preTest = true;
           
            for (int i = 0; i < prePrimeList.length; i++) {
            	divisor = BigInteger.valueOf(prePrimeList[i]);
            	if ((p.mod(divisor)==BigInteger.valueOf(0)) 
            			&& (compare(multiply(divisor,divisor),p)<=0)) {
            		preTest = false; 
            		break;
            	}
            }
            if (!preTest) {
            	continue;
            }
            else {
            	//-----------------2nd/ pass pre-test, continue to test Miller Rabin-----------------//
                if (!testMillerRabin(p,rounds,random)) {
                	continue;
                }
                else {
                	return p;
                }
            }
    	}
    }
    
	private boolean testMillerRabin(BigInteger p, int iterations, Random rnd) {
        /*return true if BigInteger p pass test Miller Rabin, return false if not*/
    	
    	// Find a and m such that m is odd and p == 1 + 2^a * m
    	BigInteger ONE = BigInteger.valueOf(1);
    	BigInteger TWO = BigInteger.valueOf(2);
        BigInteger evenComponent = subtract(p,ONE);
        BigInteger m = evenComponent;
        int a = m.getLowestSetBit();
        m = m.shiftRight(a);

        for (int i=0; i < iterations; i++) {
            // Pick a uniform random number b in range of [1,p]
            BigInteger b;
            do {
                b = new BigInteger(p.bitLength(), rnd);
            } while (compare(b, ONE) <= 0 || compare(b, p) >= 0);

            int j = 0;
            BigInteger z = powAndMod(b,m,p); //compute  z = (b^m) % p
            // Below loop mainly runs 'm-1' times.
            while (!((j == 0 && z.equals(ONE)) || z.equals(evenComponent))) { //Do following while x doesn't become p-1.
                if (j > 0 && z.equals(ONE) || ++j == a)
                    return false;
                z = powAndMod(z,TWO,p);   //z = (z*z) % p, if z = 1 return false, if z = p-1 return true
            }
        }
        return true;
    }
    
    
    /*----------------------tinh gcd------------------------------------*/
	private BigInteger gcd(BigInteger num1, BigInteger num2) {
    	BigInteger r;
    	while(compare(num2, BigInteger.ZERO) == 1) {
    		r = num1.mod(num2);
    		num1 = num2;
    		num2 = r;
    	}
    	return num1;
    }
    
    /*---------------------------------sinh 2 so nguyen to lon p, q---------------------------------------------------*/
	private void generatePrimeNumbers() {
        this.p = genRandBigPrimeNum(bitLength);
        do {
        	this.q = genRandBigPrimeNum(bitLength);
    	} while( compare(this.q, this.p ) == 0 ) ;
    }

	private void createKeys() {   
        this.n = multiply(this.p, this.q);
        // part of public key
        BigInteger phi = multiply(subtract(this.p, BigInteger.ONE), subtract(this.q, BigInteger.ONE));
    	//Find e
        BigInteger i = phi.divide(BigInteger.TWO);
        while(compare(i, phi) < 0) {
        	if(compare(gcd(i, phi), BigInteger.ONE) == 0) {
        		this.e = i;
        		if(compare(multiply(this.e, this.e).mod(phi), BigInteger.ONE) == 0) {
        			i = add(i, BigInteger.ONE);
        			continue;
        		}
        		break;
        	}
        	i = add(i, BigInteger.ONE);
        }
        // private key d
        this.d = findModInverse();
    }
    
    // find d with e, q, p
	private BigInteger findModInverse() {
    	BigInteger e = this.e;
        BigInteger phi = multiply(subtract(this.p, BigInteger.ONE), subtract(this.q, BigInteger.ONE));
        BigInteger phi0 = phi;
        BigInteger y = BigInteger.ZERO; 
        BigInteger x = BigInteger.ONE;
        if (compare(phi, BigInteger.ONE) == 0)
                return BigInteger.ZERO;
 
        while (compare(BigInteger.ONE, e) == -1) {
            // q1 is quotient
            BigInteger q1 = e.divide(phi);
 
            BigInteger t = phi;
 
            // m is remainder now, process
            // same as Euclid's algo
            phi = e.mod(phi);
            e = t;
            t = y;

            // Update x and y
            y = x.subtract(q1.multiply(y));//y = subtract(x, multiply(q1, y));
            x = t;
        }
 
        // Make x positive
        if (compare(x, BigInteger.ZERO) == -1)
            x = add(x, phi0);
        return x ;
    }

	//-------------------------- chuyen doi 1 ky tu sang chuoi binary theo unicode---------------------------//
	private String charToBinary(char ch) { 
    	String bin = Integer.toBinaryString(ch);
    	while(bin.length() < 16) {
    		bin = '0' + bin;
    	}
    	return bin;
    }
    
    //-------------------------- chuyen doi 1 chuoi ky tu sang chuoi binary theo unicode -----------------------//
	private String textToBinary(String str) {
    	String s = "";
    	for(int i = 0; i < str.length(); i++) {
    		s = s + charToBinary(str.charAt(i));
    	}
    	return s;
    }
    
    //--------------------------  chuyen doi 1 chuoi binary sang so nguyen lon -----------------------//
	private BigInteger binaryToDecimal(String binaryString) {
    	return new BigInteger(binaryString, 2);
    }
    
    //-------------------------- chuyen doi 1 so nguyen lon sang 1 chuoi binary ----------------------//
	private String decimalToBinary(BigInteger bigInt) {
    	String s = bigInt.toString(2);
    	while(s.length() % 16 != 0) {
    		s = '0' + s;
    	}
    	return s;
    }
    
    //-------------------------- chuyen 1 chuoi binary ve lai chuoi unicode ----------------------//
	private String binaryToText(String binString) {
        String text = "";
        for (int i = 0; i < binString.length() / 16; i++) {
            int a = Integer.parseInt(binString.substring(16 * i, (i + 1) * 16), 2);
            text += (char)(a);
        }
        return text;
    }
  
    //--------------------------  a^b mod c----------------------//
	private BigInteger powAndMod(BigInteger base, BigInteger exponent, BigInteger modulus) {
    	/*using repeated squaring - source: http://cs.brown.edu/courses/cs007/comp/node2.html */
    	BigInteger result = BigInteger.ONE;
		String binaryString = exponent.toString(2);
		  
		if(compare(base, BigInteger.ZERO) == 0)
			return BigInteger.ZERO;
		if(compare(exponent, BigInteger.ZERO) == 0 && compare(modulus, BigInteger.ONE) == 0)
			return BigInteger.ZERO;
		for(int i = 0; i < exponent.bitLength(); i++) {
			result = multiply(result, result).mod(modulus);
			if(binaryString.charAt(i) == '1') {
				result = multiply(result, base).mod(modulus);
			}
		}
		return result;
    }
    
	//--------------------------  encryption ----------------------//
	private String encryptAlgorithm(String message) throws Exception {
        BigInteger bigIntCiphertext = BigInteger.ZERO;
        BigInteger bigIntMessage = binaryToDecimal(textToBinary(message));
        // dieu kien: 0 < M < n
        if(compare(BigInteger.ZERO, bigIntMessage) <= 0 && compare(bigIntMessage, this.n) < 0) { 
            bigIntCiphertext = powAndMod(bigIntMessage, this.e, this.n);
        }
        else {
            throw new Exception("The message must be smaller than n.");
        }
        return bigIntCiphertext.toString(RADIX);
    }
    
    private String[] encrypt(String message) throws Exception {
        String[] ciphertextArray;
        // * ----------------------------- if message is a number -----------------*//
        if (message.matches("^-?[0-9]+$")) {       
        	int mlen = message.length();
        	int midIndex;
        	if ((mlen*16) >= this.bitLength) { // string M converted to number m > n then break M into 2 submessages to encrypt
        		if (mlen %2 == 0) midIndex = mlen / 2;
        		else midIndex = (mlen-1)/2;
        		String[] array1 = encrypt(message.substring(0,midIndex));   // array1 stores encryption of submessage1
        		String[] array2 = encrypt(message.substring(midIndex,mlen));  // array2 stores encryption of submessage2
        		ciphertextArray = Arrays.copyOf(array1, array1.length + array2.length);
        	    System.arraycopy(array2, 0, ciphertextArray, array1.length, array2.length);  //join 2 arrays and store result in ciphertextArray
        	}
        	else {
        		ciphertextArray = new String[1];
                ciphertextArray[0] = encryptAlgorithm(message);
        	}
        }
        // *---------------else if message contains characters other than numeric characters-----------------*//
        else {                                      
            ciphertextArray = new String[message.length()];
            for (int i = 0; i < message.length(); i++) {
                String temp = String.valueOf(message.charAt(i));
                ciphertextArray[i] = encryptAlgorithm(temp);
            }
        }
        return ciphertextArray;
    }
    
    //--------------------------  decryption-----------------------//
	private String decryptAlgorithm(String ciphertext) {
        BigInteger bigIntPlaintext = BigInteger.ZERO;
        BigInteger bigIntCiphertext = new BigInteger(ciphertext, RADIX);
        bigIntPlaintext = powAndMod(bigIntCiphertext, this.d, this.n);
        return binaryToText(decimalToBinary(bigIntPlaintext));
    }
    
    private String decrypt(String[] ciphertextArray) {
        String[] deciphertextArray = new String[ciphertextArray.length];
        for (int i=0;i<ciphertextArray.length;i++) {
            deciphertextArray[i] = decryptAlgorithm(ciphertextArray[i]);
        }
        StringJoiner joiner = new StringJoiner("");
          for(int i = 0; i < deciphertextArray.length; i++) {
              joiner.add(deciphertextArray[i]);
          }
        return joiner.toString();
    }
    
    //--------------------------  main program -----------------------//
    public static void main(String[] args) throws IOException, Exception, NumberFormatException {		
    	while(true) {
    		try {
        		InputStreamReader read=new InputStreamReader(System.in);  
            	BufferedReader in=new BufferedReader(read);  
            	
            	System.out.print(">>> Enter the number of bits to start generating large primes: ");    	
                int numBits = Integer.parseInt(in.readLine());
                
        		RSA rsa = new RSA(numBits);
        		System.out.println("---------------------------------------------------------------------");
        		System.out.println(">>> Generated prime numbers p and q:");
        		System.out.println("[p] = " + rsa.getp().toString(10).toUpperCase());
        		System.out.println("[q] = " + rsa.getq().toString(10).toUpperCase());
        		System.out.println("---------------------------------------------------------------------");
        		
        		System.out.println(">>> The public key is the pair (n, e) which will be published.");
        		System.out.println("[n] = " + rsa.getn().toString(10).toUpperCase());
        		System.out.println("[e] = " + rsa.gete().toString(10).toUpperCase());
        		System.out.println("---------------------------------------------------------------------");

        		System.out.println(">>> The private key is the pair (n, d) which will be kept private.");
        		System.out.println("[n] = " + rsa.getn().toString(10).toUpperCase());
        		System.out.println("[d] = " + rsa.getd().toString(10).toUpperCase());
        		System.out.println("---------------------------------------------------------------------");
        		
        		while(true) {
        			try {
        				// Get message (plaintext) from user
        				System.out.println(">>> Please enter message:");
        				String plaintext = in.readLine();
        				System.out.println("-------------------------------------------------------------");
        				
        				// Encrypt Message
                        System.out.println("Encrypting...");
                        final long startTime3 = System.nanoTime();  
                        String[] ciphertextArray = rsa.encrypt(plaintext);
                        final long endTime3 = System.nanoTime();
                        System.out.println("Encryption completed! Time for encrypting is: " + (endTime3 - startTime3) + "ns (~" + (float)(endTime3 - startTime3)/1000000 + "ms)");
                        System.out.println(">>> Ciphertext:");
                        for (int i = 0; i < ciphertextArray.length;i++) 
                            System.out.print(ciphertextArray[i]+" ");
                        System.out.println("\n-------------------------------------------------------------");
                        
        				// Decrypt Ciphertext
        				System.out.println("Decrypting...");
        				final long startTime4 = System.nanoTime();
        				String recoveredPlaintext = rsa.decrypt(ciphertextArray);
        				final long endTime4 = System.nanoTime();
        				System.out.println("Decryption completed! Time for decrypting is: " + (endTime4 - startTime4) + "ns (~" + (float)(endTime4 - startTime4)/1000000 + "ms)");
        				System.out.println(">>> Recovered plaintext:");
        				System.out.println(recoveredPlaintext);
        				System.out.println("_______________________________________________________________________________________________________________________________________" );
        			} catch(Exception e) {
        				System.out.println(e);
        			}
        		}
        	} catch(NumberFormatException e) {
        		System.out.println("You must to enter an integer.");
        	}	
    	}
    }
}
