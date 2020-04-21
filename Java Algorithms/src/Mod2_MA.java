/*
 * Author: Nevin George
 * Advisor: Dana Angluin
 * Program Description: The algorithm takes in as input a mod-2-MA and prints to stdout the MA obtained after 
 * learning the input function through a series of membership and equivalence queries. The motivation behind this 
 * algorithm originally arose from Angluin's exact learning model described in her paper "Learning regular sets 
 * from queries and counterexamples."
 * 
 * References:
 * 1 Amos Beimel, Francesco Bergadano, Nader H. Bshouty, Eyal Kushilevitz, Stefano Varric- chio. Learning 
 *   functions represented as multiplicity automata. J. ACM, 47(3):506–530, May 2000.
 * 2 Dana Angluin. Learning regular sets from queries and counterexamples. Inf. Comput., 75(2):87–106, 1987.
 * 3 Dana Angluin, Timos Antonopoulos, Dana Fisman. Strongly Unambiguous Büchi Automata Are Polynomially 
 *   Predictable with Membership Queries. 28th International Conference on Computer Science Logic, 8:1–8:17, 2020.
 * 4 Michael Thon and Herbert Jaeger. Links Between Multiplicity Automata, Observable Operator Models and 
 *   Predictive State Representations — a Unified Learning Framework. Journal of Machine Learning Research, 
 *   16(4):103−147, 2015.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.StringTokenizer;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class Mod2_MA {
	
	// if true, displays the observation table as it is being constructed
	public static boolean verbose;
	
	// alphabet
	public static int alphabetSize;
	public static Character[] alphabet;
	// maps letters in alphabet to a corresponding index
	public static HashMap<Character, Integer> letterToIndex;
	
	// γ of target function
	public static double[] fy;
	// set of nxn μ's of target function
	public static double[][][] fu;
	// size of target function
	public static int r;
	
	// Hankel matrix
	public static ArrayList<Integer> F;
	public static int sizeF;
	// maps words to indices of F
	public static HashMap<String, Integer> wordToIndex;
	
	// row indices of the observation table
	public static ArrayList<String> X;
	// column indices of the observation table
	public static ArrayList<String> Y;
	// size of X and Y
	public static int l;
	
	// counter-example
	public static String z;
	// subset of z
	public static String w;
	// additional character of the prefix ω + σ
	public static char sig;
	// experiment
	public static String y;
		
	// γ of the learned function
	public static double[] resulty;
	// set of nxn μ's of the learned function
	public static double[][][] resultu;
	
	// takes in user input
	public static Scanner in;
	
	public static void main(String[] args) throws Exception {
		initialize();
		run();
		
		// performs desired operations with the learned mod-2-MA
		operations(true);
	}
	
	public static void initialize() throws Exception {
		/* The input file must have the following form (no line separation, characters are space separated, 
		 * and lines beginning with // are ignored):
		 * <alphabet size>
		 * <characters in the alphabet>
		 * <size of the target function (r)>
		 * <γ of the target function (fy)>
		 * List of μ's for each character in the alphabet, each μ appears in a rxr grid
		 * 
		 * Example input files can be found in the GitHub repository.
		 */
		
		// reads in file name + optional flag -v from stdin
		System.out.println("Input file name and optional flag -v (e.g. input1.txt or input1.txt -v)");
		in = new Scanner(System.in);
		String[] arrInput = in.nextLine().split(" ");
		verbose = false;
		if(arrInput.length == 2 && arrInput[1].equals("-v"))
			verbose = true;
		BufferedReader f = new BufferedReader(new FileReader(arrInput[0]));
		System.out.println();
		
		// alphabet size
		alphabetSize = Integer.parseInt(readInput(f));
		
		// alphabet
		StringTokenizer st = new StringTokenizer(readInput(f));
		alphabet = new Character[alphabetSize];
		for(int i=0;i<alphabetSize;i++) {
			String letter = st.nextToken();
			if(letter.length()!=1)
				throwException(f,"Invalid input: invalid character in the alphabet");
			alphabet[i] = letter.charAt(0);
		}
		if(st.hasMoreTokens())
			throwException(f,"Invalid input: alphabet size exceeds the specified size");
		
		// maps each letter in alphabet to an index
		letterToIndex = new HashMap<Character, Integer>();
		for(int i=0;i<alphabetSize;i++)
			letterToIndex.put(alphabet[i], i);
		
		// size of the target function
		r = Integer.parseInt(readInput(f));
		
		// γ of the target function
		st = new StringTokenizer(readInput(f));
		fy = new double[r];
		for(int i=0;i<r;i++)
			fy[i] = Integer.parseInt(st.nextToken());
		if(st.hasMoreTokens())
			throwException(f,"Invalid input: γ length exceeds the specified size");
		
		// set of μ's for the target function
		fu = new double[alphabetSize][r][r];
		for(int i=0;i<alphabetSize;i++) {
			for(int j=0;j<r;j++) {
				st = new StringTokenizer(readInput(f));
				for(int k=0;k<r;k++)
					fu[i][j][k] = Integer.parseInt(st.nextToken());
				if(st.hasMoreTokens())
					throwException(f,"Invalid input: μ size exceeds the specified size");
			}
		}
		if(readInput(f) != null)
			throwException(f,"Invalid input: μ size exceeds the specified size");
		
		f.close();
	}
	
	public static String readInput(BufferedReader f) throws IOException {
		String line = f.readLine();
		// ignores lines beginning with "//"
		while(line != null && line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();
		return line;
	}
	
	public static void throwException(BufferedReader f, String message) throws Exception {
		// properly closes input streams
		if(f != null)
			f.close();
		if(in != null)
			in.close();
		throw new Exception(message);
	}
	
	public static void run() throws Exception {
		// initializes the rows and columns of the observation table
		X = new ArrayList<String>();
		Y = new ArrayList<String>();
		X.add("");
		Y.add("");
		l = 1;
		
		// initializes the Hankel matrix
		F = new ArrayList<Integer>();
		sizeF = 0;
		wordToIndex = new HashMap<String, Integer>();
		
		if(verbose)
			System.out.println("Results after individual queries\n--------------------------------");
		
		/* f("") cannot equal 0 (otherwise can't form a linearly independent basis of elements in X).
		 * The algorithm instead begins with a 2x2 matrix of full rank.
		 */
		if(MQ("")==0) {
			double[] hy = createHY();
			double[][][] setOfHu = createHU();
			
			// generates a counter-example z
			if(!EQ(hy, setOfHu)) {
				l++;
				X.add(z);
				Y.add(z);
			}
		}
		
		// runs the algorithm
		learnMA();
		
		// statistical final check of equivalence
		if(finalCheck())
			displayResults();
		else
			throwException(null,"Algorithm failed: failed final check.");
	}
	
	public static void learnMA() throws Exception {
		if(verbose)
			displayQueries();
		
		// creates the γ for the hypothesis
		double[] hy = createHY();
		// creates the set of μ's for the hypothesis
		double[][][] hu = createHU();
		
		// sees if the hypothesis = target function, if so returns the hypothesis
		if(EQ(hy, hu)) {
			resulty = hy;
			resultu = hu;
			return;
		}
		
		w = "";
		sig = 0;
		y = "";
		// attempts to calculate ω, σ, and y, if it cannot find a y that works it throws an exception
		calcWSigY(hu);
		
		if(l==r)
			throwException(null,"Algorithm failed: size of the hypothesis exceeds that of the target function.");
		
		// updates l, X, and Y
		l++;
		X.add(w);
		Y.add(sig+y);
		
		learnMA();
	}
	
	public static double[] createHY() {
		// γ is the set of results obtained after performing membership queries on the indices in X
		double[] y = new double[l];
		for(int i=0;i<l;i++)
			y[i] = MQ(X.get(i));
		return y;
	}
	
	public static double[][][] createHU() {
		/*
		 * For every s, define a matrix μ by letting its i-th row be the coefficients of the vector F_{xi+σ}(y) 
		 * when expressed as a linear combination of the vectors F_x1 to F_xl (such coefficients exist as F_x1 
		 * to F_xl are linearly independent).
		*/
		
		double[][][] setOfU = new double[alphabetSize][l][l];
		for(int c=0;c<alphabetSize;c++) {
			// calculuates μ_c
			char sig = alphabet[c];
			
			// calculates the vectors F_xi and F_{xi+σ}
			double[][] F_xi = new double[l][l];
			double[][] F_xi_sigma = new double[l][l];
			for(int i=0;i<l;i++) {
				for(int j=0;j<l;j++) {
					F_xi[j][i] = MQ(X.get(i) + Y.get(j));
					F_xi_sigma[i][j] = MQ(X.get(i) + sig + Y.get(j));
				}
			}
			
			// solves the matrix equation using LU Decomposition
			RealMatrix coefficients = new Array2DRowRealMatrix(F_xi);
			DecompositionSolver solver = new LUDecomposition(coefficients).getSolver();
			for(int i=0;i<l;i++) {
				RealVector constants = new ArrayRealVector(F_xi_sigma[i]);
				try {
					RealVector solution = solver.solve(constants);
					for(int j=0;j<l;j++)
						setOfU[c][i][j] = mod2(solution.getEntry(j));
				}
				// matrix is not invertible
				catch(Exception e) {
					for(int j=0;j<l;j++)
						setOfU[c][i][j] = 0;
				}
			}
			
		}
		return setOfU;
	}
	
	public static int MQ(String w) {
		// MQ for the target function
		
		// MQ(ω) was previously calculated and is in the Hankel matrix
		if(wordToIndex.get(w) != null)
			return F.get(wordToIndex.get(w));
		
		// initializes cur as the rxr identity matrix
		RealMatrix cur = MatrixUtils.createRealIdentityMatrix(r);
		
		// multiplies cur by the corresponding μ for each letter in ω
		for(int i=0;i<w.length();i++) {
			cur = cur.multiply(MatrixUtils.createRealMatrix(fu[letterToIndex.get(w.charAt(i))]));
			// accounts for rounding issues with larger words
			if(i!=0 && i%25==0) {
				for(int j=0;j<r;j++) {
					for(int k=0;k<r;k++)
						cur.setEntry(j, k, mod2(cur.getEntry(j, k)));
				}
			}
		}
		
		// multiplies the final result with γ
		int out = mod2(cur.getRowVector(0).dotProduct(new ArrayRealVector(fy)));
		
		// adds MQ(ω) to the Hankel matrix
		wordToIndex.put(w, sizeF++);
		F.add(out);
		
		return out;
	}
	
	public static int MQResults(String w) {
		// MQ for the learned function
		
		// initializes cur as the lxl identity matrix
		RealMatrix cur = MatrixUtils.createRealIdentityMatrix(l);
		
		// multiplies cur by the corresponding μ for each letter in ω
		for(int i=0;i<w.length();i++) {
			cur = cur.multiply(MatrixUtils.createRealMatrix(resultu[letterToIndex.get(w.charAt(i))]));
			// accounts for rounding issues with larger words
			if(i!=0 && i%25==0) {
				for(int j=0;j<l;j++) {
					for(int k=0;k<l;k++)
						cur.setEntry(j, k, mod2(cur.getEntry(j, k)));
				}
			}
		}
		
		// multiplies the final result with γ
		return mod2(cur.getRowVector(0).dotProduct(new ArrayRealVector(resulty)));
	}
	
	public static boolean EQ(double[] hy, double[][][] hu) {		 
		/* EQ constructs the MA formed by combining the target function and hypothesis.
		 * 
		 * Each μ(ω) of the combined MA has the following form (the 0's representing block 0 matrices):
		 * |fu(ω)   0  |
		 * |  0   hu(ω)|
		 * 
		 * γ has the form [fy hy] (fy and hy are joined together in the same vector).
		 * 
		 * The initial vector is the initial vectors of the target and hypothesis joined together.
		 * 
		 * This MA represents the XOR of the target function and hypothesis.
		 * We will construct a basis for the set span(μ(ω)γ : ω∈Σ*).
		 * If for every vector v in the basis v[0] + v[r] == 0, then the MA outputs 0 for every possible input 
		 * word and hence the target and hypothesis are equivalent.
		 * If there exists a vector v in the basis such that v[0] + v[r] == 1, the corresponding ω of v will be 
		 * returned as the counter-example.
		 */
		
		// set of μ for the combined MA
		double[][][] setOfU = new double[alphabetSize][r+l][r+l];
		for(int i=0;i<alphabetSize;i++) {
			for(int j=0;j<r+l;j++) {
				for(int k=0;k<r+l;k++) {
					// fu forms the upper left block of μ
					if(j<r && k<r)
						setOfU[i][j][k] = fu[i][j][k];
					// hu forms the lower right block of μ
					else if(j>=r && k>=r)
						setOfU[i][j][k] = hu[i][j-r][k-r];
					// everything else is 0
					else
						setOfU[i][j][k] = 0;
				}
			}
		}
		
		// γ for the combined MA
		double[] y = new double[r+l];
		for(int i=0;i<r+l;i++) {
			// γ has the form [fy hy]
			if(i<r)
				y[i] = fy[i];
			else
				y[i] = hy[i-r];
		}
		
		// To form the basis, we will follow algorithm 1 detailed in the paper by Thon and Jaeger.
		// Basis for the set span(μ(ω)γ : ω∈Σ*)
		double[][] B = new double[r+l][r+l];
		int sizeB = 0;
		// Contains the corresponding ω for every element in B
		ArrayList<String> WB = new ArrayList<String>();
		
		// Set with elements to try to add to B, begin with y
		ArrayList<double[]> C = new ArrayList<double[]>();
		// Contains the corresponding ω for every element in C
		ArrayList<String> WC = new ArrayList<String>();
		C.add(y);
		WC.add("");
		int sizeC = 1;
		
		while(sizeC>0) {
			// element to test
			double[] w = C.remove(0);
			String s = WC.remove(0);
			sizeC--;
			
			// tests if ω is linearly independent of B
			if(linInd(w, B, sizeB)) {
				// found a counter-example
				if(mod2(w[0]+w[r]) == 1) {
					z = s;
					return false;
				}
				
				// extends B
				B[sizeB++] = w;
				WB.add(s);
				
				// adds {μ(σ)ω | σ∈Σ} to C
				for(int i=0;i<alphabetSize;i++) {
					RealMatrix m = MatrixUtils.createRealMatrix(setOfU[i]);
					RealVector p = MatrixUtils.createRealVector(w);
					double[] v = m.operate(p).toArray();
					C.add(v);
					WC.add(alphabet[i]+s);
					sizeC++;
				}
			}
		}
		
		// v[0] + v[r] == 0 for every vector v in the basis, so the target and hypothesis are equivalent
		return true;
	}
	
	public static boolean linInd(double[] w, double[][] B, int sizeB) {
		if(sizeB==0)
			return true;
		
		// forms the augmented matrix B|w
		int numRows = w.length;
		int numCols = sizeB+1;
		
		RealMatrix m = MatrixUtils.createRealMatrix(numRows, numCols);
		for(int i=0;i<sizeB;i++)
			m.setColumn(i,B[i]);
		m.setColumn(numCols-1, w);
		
		// put the augmented matrix in rref
		int r=0;
		for(int c=0;c<numCols && r<numRows;c++) {
			int j = r;
			for(int i=r+1;i<numRows;i++)
				if(mod2(m.getEntry(i, c)) >mod2(m.getEntry(j, c)))
					j = i;
			if(mod2(m.getEntry(j, c)) == 0)
				continue;

			RealMatrix temp = m.getRowMatrix(j);
			m.setRowMatrix(j,m.getRowMatrix(r));
			m.setRowMatrix(r,temp);

			for(int i=0;i<numRows;i++) {
				if(i!=r) {
					int t = mod2(m.getEntry(i, c));
					for(j=0;j<numCols;j++)
						m.setEntry(i, j, mod2(m.getEntry(i,j) - (t * m.getEntry(r, j))));
				}
			}
			r++;
		}
		
		// finds the index of the last 1 in the last column (if exists)
		int index = -1;
		for(int i=numRows-1;i>=0;i--) {
			if(mod2(m.getEntry(i, numCols-1)) == 1) {
				index = i;
				break;
			}
		}
		
		// last vector is the 0 vector, in span(B)
		if(index == -1)
			return false;
		
		// checks whether in span
		for(int j=0;j<numCols-1;j++) {
			if(mod2(m.getEntry(index, j)) == 1)
				return false;
		}
		
		// linearly independent
		return true;
	}
	
	public static void calcWSigY(double[][][] hu) throws Exception {
		// goes through every possible prefix of z starting with ω = "" and σ = (first character of ω)
		// prefix = ω + σ
		for(int i=0;i<z.length();i++) {
			if(i!=0)
				w = z.substring(0,i);
			sig = z.charAt(i);
			
			// calculates μ(ω)
			RealMatrix u = MatrixUtils.createRealIdentityMatrix(l);
			for(int n=0;n<w.length();n++)
				u = u.multiply(MatrixUtils.createRealMatrix(hu[letterToIndex.get(w.charAt(n))]));
			
			// goes through every possible value of y in Y
			// equation is F_{ω+σ}(y) != sum((μ(ω)_1,i) * F_{xi+σ}(y))
			for(int j=0;j<l;j++) {
				y = Y.get(j);
				
				int left = MQ(w+sig+y);
				
				int right = 0;
				for(int k=0;k<l;k++)
					right += mod2(u.getEntry(0, k)) * MQ(X.get(k) + sig + y);
				right = mod2(right);
				
				// found a solution, values we want to return are set using global variables
				if(left!=right)
					return;
			}
		}
		
		throwException(null,"Algorithm failed: didn't find a suitable omega, sigma, and gamma.");
	}

	public static void displayResults() {
		System.out.println("Learned mod-2-MA");
		System.out.println("----------------");
		// prints γ
		System.out.print("y: ");
		String s = "";
		for(int i=0;i<resulty.length;i++)
			s += mod2(resulty[i]) + " ";
		System.out.println(s + "\n");
		
		// prints the μ's
		System.out.println("Set of u:\n");
		for(int i=0;i<resultu.length;i++) {
			System.out.println("Letter " + alphabet[i]);
			for(int j=0;j<resultu[i].length;j++) {
				s = "";
				for(int k=0;k<resultu[i].length;k++)
					s += mod2(resultu[i][j][k]) + " ";
				System.out.println(s);
			}
			System.out.println();
		}
	}
	
	public static void displayQueries() {
		// displays the observation table
		System.out.println("l = " + l);
		System.out.print("Rows: ɛ ");
		for(int i=1;i<X.size();i++)
			System.out.print(X.get(i) + " ");
		System.out.println();
		
		System.out.print("Cols: ɛ ");
		for(int i=1;i<Y.size();i++)
			System.out.print(Y.get(i) + " ");
		
		System.out.println("\nTable:");
		for(int i=0;i<X.size();i++) {
			for(int j=0;j<Y.size();j++)
				System.out.print(MQ(X.get(i) + Y.get(j)) + " ");
			System.out.println();
		}
		System.out.println();
	}

	public static String genTest(int len) {
		// adds len number of random characters in alphabet to test
		String test = "";
		for(int i=0;i<len;i++)
			test += alphabet[(int)(Math.random()*alphabetSize)];
		return test;
	}
	
	public static boolean finalCheck() {
		// creates 40 tests of length 1-50
		// checks whether the learned function and target function have the same output
		for(int i=1;i<=40;i++) {
			String test = genTest((int)(Math.random()*50)+1);
			if(MQ(test)!=MQResults(test))
				return false;
		}
		return true;
	}
	
	public static void operations(boolean inMod2MA) {
		System.out.println("Available operations for the learned Mod-2-MA (enter \"quit\" to terminate):");
		if(inMod2MA)
			System.out.println("- Test whether a word is accepted: \"-a <word>\"");
		else
			System.out.println("- Test whether a word of the form u$v in (L)_$ is accepted: \"-a <word>\"");
		while(true) {
			// reads in cmd
			String line = in.nextLine();
			
			// terminate the program
			if(line.equals("quit")) {
				in.close();
				break;
			}
			
			String[] input = line.split(" ");
			if(input.length == 0)
				System.out.println("Invalid cmd");
			
			// tests whether a word is accepted
			if(input[0].equals("-a")) {
				if(input.length>2)
					System.out.println("Invalid cmd");
				else {
					// if no word is passed, the test is the empty string
					String test = "";
					if(input.length==2)
						test = input[1];
					
					if(!inAlphabet(test, inMod2MA))
						System.out.println("Inputted word is not in the language.");
					else if(MQResults(test) == 1)
						System.out.println("Accepted");
					else
						System.out.println("Not accepted");
				}
			}
			else
				System.out.println("Invalid cmd");
		}
	}
	
	public static boolean inAlphabet(String word, boolean inMod2MA) {
		// a word in (L)_$ must contain exactly one '$'
		boolean containsDollar = false;
		for(int i=0;i<word.length();i++) {
			if(letterToIndex.get(word.charAt(i)) == null)
				return false;
			if(word.charAt(i) == '$') {
				if(containsDollar)
					return false;
				else
					containsDollar = true;
			}
		}
		if(!inMod2MA && !containsDollar)
			return false;
		return true;
	}
	
	public static int mod2(double n) {
		int temp = (int) Math.round(n);
		if(temp%2==0)
			return 0;
		else
			return 1;
	}
}