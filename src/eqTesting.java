import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;

public class eqTesting {
	// alphabet
	public static int alphabetSize;
	public static Character[] alphabet;
	// maps letters in alphabet to a corresponding index
	public static HashMap<Character, Integer> letterToIndex;
	
	// array of tests used to see if a hypothesis is correct
	public static String[] tests;
	// used in fillTests
	public static int index;
	public static int maxTestSize;
	
	// γ of target function
	public static int[] fy1;
	// set of nxn μ's of target function
	public static int[][][] fu1;
	// size of target function
	public static int r1;
	
	// γ of target function
	public static int[] fy2;
	// set of nxn μ's of target function
	public static int[][][] fu2;
	// size of target function
	public static int r2;
	
	// length of "long" final check
	public static int checkLength;
		
	public static void main(String[] args) throws IOException {
		initialize();
		if(eq() && finalCheck())
			System.out.println("Equivalent!");
		else
			System.out.println("Not equivalent");

	}
	
	public static boolean eq() {
		// checks to see if the hypothesis and target functions have the same output for all words in tests
		for(int i=0;i<maxTestSize;i++) {
			if(MQH(tests[i], fy1, fu1, r1)!=MQH(tests[i], fy2, fu2, r2))
				return false;
		}
		return true;
	}
	
	public static void initialize() throws IOException {
		// read in from input file
		BufferedReader f = new BufferedReader(new FileReader("input4.txt"));
		
		// read in the alphabet
		alphabetSize = Integer.parseInt(f.readLine());
		StringTokenizer st = new StringTokenizer(f.readLine());
		alphabet = new Character[alphabetSize];
		for(int i=0;i<alphabetSize;i++)
			alphabet[i] = st.nextToken().charAt(0);
		
		// size of the target function
		r1 = Integer.parseInt(f.readLine());
		
		// γ of the target function
		fy1 = new int[r1];
		st = new StringTokenizer(f.readLine());
		for(int i=0;i<r1;i++)
			fy1[i] = Integer.parseInt(st.nextToken());
		
		// set of μ's for the target function
		fu1 = new int[alphabet.length][r1][r1];
		for(int i=0;i<alphabetSize;i++) {
			for(int j=0;j<r1;j++) {
				st = new StringTokenizer(f.readLine());
				for(int k=0;k<r1;k++)
					fu1[i][j][k] = Integer.parseInt(st.nextToken());
			}
		}
		f.close();
		
		
		
		f = new BufferedReader(new FileReader("file1.txt"));
		
		// size of the target function
		r2 = Integer.parseInt(f.readLine());
		
		// γ of the target function
		fy2 = new int[r2];
		st = new StringTokenizer(f.readLine());
		for(int i=0;i<r2;i++)
			fy2[i] = Integer.parseInt(st.nextToken());
		
		// set of μ's for the target function
		fu2 = new int[alphabet.length][r2][r2];
		for(int i=0;i<alphabetSize;i++) {
			for(int j=0;j<r2;j++) {
				st = new StringTokenizer(f.readLine());
				for(int k=0;k<r2;k++)
					fu2[i][j][k] = Integer.parseInt(st.nextToken());
			}
		}
		f.close();
		
		// build the list of tests
		
		// will run EQ through maxTestSize number of tests
		maxTestSize = 100000;
		tests = new String[maxTestSize];
		index = 0;
		
		// each time fill is called, add all tests of length wordLen to tests
		// stops when tests is full
		int wordLen = 0;
		while(fillTests("",0,wordLen++));
		
		// maps each letter in alphabet to an index
		letterToIndex = new HashMap<Character, Integer>();
		for(int i=0;i<alphabetSize;i++)
			letterToIndex.put(alphabet[i], i);
	}
	
	public static boolean fillTests(String word, int curLen, int wordLen) {
		// completely filled up tests
		if(index == maxTestSize)
			return false;
		
		// found a test of length wordLen
		if(curLen==wordLen) {
			tests[index++] = word;
			return true;
		}
		
		// create every word of +1 length
		boolean cont = true;
		for(int i=0;i<alphabetSize;i++) {
			if(!fillTests(word + alphabet[i], curLen+1, wordLen)) {
				cont = false;
				break;
			}
		}
		return cont;
	}
	
	public static int MQH(String w, int[] hy, int[][][] hu, int l) {
		// MQ for the current hypothesis
		int[][] cur = new int[l][l];
		for(int i=0;i<l;i++)
			cur[i][i] = 1;
		
		for(int i=0;i<w.length();i++)
			cur = nxnMatrixMult(cur, hu[letterToIndex.get(w.charAt(i))]);
		
		return dotProduct(cur[0],hy);
	}
	
	public static int[][] nxnMatrixMult(int[][] arr1, int[][] arr2){
		int size = arr1.length;
		int[][] out = new int[size][size];
		for(int i=0;i<size;i++) {
			for(int j=0;j<size;j++) {
				for(int k=0;k<size;k++)
					out[i][j] += arr1[i][k]*arr2[k][j];
				out[i][j] = mod2(out[i][j]);
			}
		}
		return out;
	}
	
	public static int dotProduct(int[] v1, int[] v2) {
		int sum = 0;
		for(int i=0;i<v1.length;i++)
			sum += v1[i]*v2[i];
		return mod2(sum);
	}
	
	public static int mod2(int n) {
		if(n%2==0)
			return 0;
		else
			return 1;
	}
	
	public static String genLongTest() {
		int size = checkLength;
		String longTest = "";
		for(int i=0;i<size;i++)
			longTest += alphabet[(int)(Math.random()*alphabetSize)];
		return longTest;
	}
	
	public static boolean finalCheck() {
		for(int i=0;i<100;i++) {
			// creates a "long" test and sees if the hypothesis and target function have the same output
			String longTest = genLongTest();
			if(MQH(longTest, fy1, fu1, r1)!=MQH(longTest, fy2, fu2, r2))
				return false;
		}
		return true;
	}
}
