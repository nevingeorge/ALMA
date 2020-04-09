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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.StringTokenizer;

public class Mod2_MA {
	
	// alphabet
	public static int alphabetSize;
	public static Character[] alphabet;
	// maps letters in alphabet to a corresponding index
	public static HashMap<Character, Integer> letterToIndex;
	
	// γ of target function
	public static int[] fy;
	// set of nxn μ's of target function
	public static int[][][] fu;
	// size of target function
	public static int r;
	// γ the algorithm produces
	public static int[] resulty;
	// set of nxn μ's the algorithm produces
	public static int[][][] resultu;
	
	// counter-example
	public static String z;
	// subset of z
	public static String w;
	// additional character of the prefix ω + σ
	public static char sig;
	// experiment
	public static String y;
	
	// length of current X and Y
	public static int l;
	
	// length of "long" final check
	public static int checkLength;
	// number of "long" strings to check
	public static int numToCheck;
	
	public static void main(String[] args) throws Exception {
		initialize();
		
		// contain the indices xi and yi
		ArrayList<String> X = new ArrayList<String>();
		ArrayList<String> Y = new ArrayList<String>();
		X.add("");
		Y.add("");
		
		/* f("") cannot equal 0 (otherwise can't form a linearly independent basis of elements in X).
		 * The algorithm instead begins with a 2x2 matrix of full rank.
		 */
		if(MQ("")==0) {
			int[] hy = createHY(l, X);
			int[][][] setOfHu = createHU(X, Y, l);
			
			// generates a counter-example z
			if(!EQ(hy, setOfHu, l)) {
				X.add(z);
				Y.add(z);
				l++;
			}
		}
		
		// runs the algorithm
		learnMA(X, Y);
		
		// additional final check, checks equivalence of numToCheck "long" randomly generated strings
		if(finalCheck())
			displayResults();
		else
			throw new Exception("Algorithm failed: failed final check.");
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
		
		// reads from the input file in quotations, ***EDIT THE FILE NAME DEPENDING ON THE INTENDED INPUT FILE***
		BufferedReader f = new BufferedReader(new FileReader("input4.txt"));
		
		// alphabet size
		String line = f.readLine();
		while(line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();
		
		alphabetSize = Integer.parseInt(line);
		
		// alphabet
		line = f.readLine();
		while(line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();
		StringTokenizer st = new StringTokenizer(line);
		
		alphabet = new Character[alphabetSize];
		for(int i=0;i<alphabetSize;i++) {
			String letter = st.nextToken();
			if(letter.length()!=1) {
				f.close();
				throw new Exception("Invalid input: invalid character in the alphabet");
			}
			alphabet[i] = letter.charAt(0);
		}
		
		if(st.hasMoreTokens()) {
			f.close();
			throw new Exception("Invalid input: alphabet size exceeds the specified size");
		}
		
		// size of the target function
		line = f.readLine();
		while(line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();
		
		r = Integer.parseInt(line);
		
		// γ of the target function
		line = f.readLine();
		while(line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();
		st = new StringTokenizer(line);
		
		fy = new int[r];
		for(int i=0;i<r;i++)
			fy[i] = Integer.parseInt(st.nextToken());
		
		if(st.hasMoreTokens()) {
			f.close();
			throw new Exception("Invalid input: γ length exceeds the specified size");
		}
		
		// set of μ's for the target function
		fu = new int[alphabetSize][r][r];
		for(int i=0;i<alphabetSize;i++) {
			for(int j=0;j<r;j++) {
				line = f.readLine();
				while(line.charAt(0) == '/' && line.charAt(1) == '/')
					line = f.readLine();
				st = new StringTokenizer(line);
				
				for(int k=0;k<r;k++)
					fu[i][j][k] = Integer.parseInt(st.nextToken());
				
				if(st.hasMoreTokens()) {
					f.close();
					throw new Exception("Invalid input: μ size exceeds the specified size");
				}
			}
		}
		
		line = f.readLine();
		while(line!=null) {
			if(line.charAt(0) != '/' && line.charAt(1) != '/') {
				f.close();
				throw new Exception("Invalid input: μ size exceeds the specified size");
			}
			line = f.readLine();
		}
		f.close();
		
		// general variables to initialize
		
		checkLength = 100;
		numToCheck = 20;
		
		// initial size of X and Y
		l = 1;
		
		// maps each letter in alphabet to an index
		letterToIndex = new HashMap<Character, Integer>();
		for(int i=0;i<alphabetSize;i++)
			letterToIndex.put(alphabet[i], i);
	}
	
	public static void learnMA(ArrayList<String> X, ArrayList<String> Y) throws Exception {
		// creates the γ for the hypothesis
		int[] hy = createHY(l, X);
		// creates the set of μ's for the hypothesis
		int[][][] hu = createHU(X, Y, l);
		
		// sees if the hypothesis = target function, if so returns the hypothesis
		if(EQ(hy, hu, l)) {
			resulty = hy;
			resultu = hu;
			return;
		}
		
		w = "";
		sig = 0;
		y = "";
		// attempts to calculate ω, σ, and y, if it cannot find a y that works it throws an exception
		calcWSigY(l, hu, X, Y);
		
		if(l==r)
			throw new Exception("Algorithm failed: size of the hypothesis exceeds that of the target function.");
		
		// updates l, X, and Y for next recursive call
		l++;
		X.add(w);
		Y.add(sig+y);
		learnMA(X, Y);
	}
	
	public static int[] createHY(int l, ArrayList<String> X) {
		// γ is the set of results obtained after performing membership queries on the indices in X
		int[] y = new int[l];
		for(int i=0;i<l;i++) 
			y[i] = MQ(X.get(i));
		return y;
	}
	
	public static int[][][] createHU(ArrayList<String> X, ArrayList<String> Y, int l) {
		/*
		 * For every s, define a matrix μ by letting its i-th row be the coefficients of the vector F_{xi+σ}(y) 
		 * when expressed as a linear combination of the vectors F_x1 to F_xl (such coefficients exist as F_x1 
		 * to F_xl are linearly independent).
		*/
		int[][][] setOfU = new int[alphabetSize][l][l];
		for(int c=0;c<alphabetSize;c++) {
			// calculuates μ_c
			char sig = alphabet[c];
			int[][] u = new int[l][l];
			
			// calculates the vectors F_x1 to F_xl
			int[][] rowsOfFxi = new int[l][l];
			for(int j=0;j<l;j++) {
				int[] Fxi = new int[l];
				for(int k=0;k<l;k++)
					Fxi[k] = MQ(X.get(j) + Y.get(k));
				rowsOfFxi[j] = Fxi;
			}
			
			int[][] inverseRoF;
			// if rowsOfFxi is not invertible, inverseRoF is a 0 matrix
			if(determinant(rowsOfFxi)==0) {
				inverseRoF = new int[l][l];
				for(int j=0;j<l;j++) {
					for(int k=0;k<l;k++)
						inverseRoF[j][k] = 0;
				}
			}
			else
				inverseRoF = inverse(rowsOfFxi);
			
			for(int i=0;i<l;i++) {
				// calculates F_{xi+σ}
				int[] Fxisig = new int[l];
				for(int j=0;j<l;j++)
					Fxisig[j] = MQ(X.get(i) + sig + Y.get(j));
				
				// calculates μ
				u[i] = rowxnxnMatrixMult(Fxisig, inverseRoF);
			}
			
			setOfU[c] = u;
		}
		return setOfU;
	}
	
	public static int MQ(String w) {
		// MQ for the target function
		
		// initializes cur as the rxr identity matrix
		int[][] cur = new int[r][r];
		for(int i=0;i<r;i++)
			cur[i][i] = 1;
		
		// multiplies cur by the corresponding μ for each letter in ω
		for(int i=0;i<w.length();i++)
			cur = nxnMatrixMult(cur, fu[letterToIndex.get(w.charAt(i))]);
		
		// multiplies the final result with γ
		return dotProduct(cur[0],fy);
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
	
	public static boolean EQ(int[] hy, int[][][] hu, int l) {		 
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
		int[][][] setOfU = new int[alphabetSize][r+l][r+l];
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
		int[] y = new int[r+l];
		for(int i=0;i<r+l;i++) {
			// γ has the form [fy hy]
			if(i<r)
				y[i] = fy[i];
			else
				y[i] = hy[i-r];
		}
		
		// To form the basis, we will follow algorithm 1 detailed in the paper by Thon and Jaeger.
		// Basis for the set span(μ(ω)γ : ω∈Σ*)
		ArrayList<int[]> B = new ArrayList<int[]>();
		// Contains the corresponding ω for every element in B
		ArrayList<String> WB = new ArrayList<String>();
		
		// Set with elements to try to add to B, begin with y
		ArrayList<int[]> C = new ArrayList<int[]>();
		// Contains the corresponding ω for every element in C
		ArrayList<String> WC = new ArrayList<String>();
		C.add(y);
		WC.add("");
		int sizeC = 1;
		
		while(sizeC>0) {
			// element to test
			int[] w = C.remove(0);
			String s = WC.remove(0);
			sizeC--;
			
			// tests if ω is linearly independent of B
			if(linInd(w, B)) {
				// found a counter-example
				if(mod2(w[0]+w[r]) == 1) {
					z = s;
					return false;
				}
				
				// extends B
				B.add(w);
				WB.add(s);
				
				// adds {μ(σ)ω | σ∈Σ} to C
				for(int i=0;i<alphabetSize;i++) {
					C.add(nxnxnx1MatrixMult(setOfU[i],w));
					WC.add(alphabet[i]+s);
					sizeC++;
				}
			}
		}
		
		// v[0] + v[r] == 0 for every vector v in the basis, so the target and hypothesis are equivalent
		return true;
	}
	
	public static void calcWSigY(int l, int[][][] hu, ArrayList<String> X, ArrayList<String> Y) throws Exception {
		// goes through every possible prefix of z starting with ω = "" and σ = (first character of ω)
		// prefix = ω + σ
		for(int i=0;i<z.length();i++) {
			if(i!=0)
				w = z.substring(0,i);
			sig = z.charAt(i);
			
			// calculates μ(ω)
			int[][] u = new int[l][l];
			for(int m=0;m<l;m++)
				u[m][m] = 1;
			
			for(int n=0;n<w.length();n++)
				u = nxnMatrixMult(u, hu[letterToIndex.get(w.charAt(n))]);
			
			// goes through every possible value of y in Y
			// equation is F_{ω+σ}(y) != sum((μ(ω)_1,i) * F_{xi+σ}(y))
			for(int j=0;j<l;j++) {
				y = Y.get(j);
				int left = MQ(w+sig+y);
				
				int right = 0;
				for(int k=0;k<l;k++)
					right += u[0][k] * MQ(X.get(k) + sig + y);
				right = mod2(right);
				
				// found a solution, values we want to return are set using global variables
				if(left!=right)
					return;
			}
		}
		
		throw new Exception("Algorithm failed: didn't find a suitable omega, sigma, and gamma.");
	}

	public static void displayResults() {
		// prints γ
		System.out.print("y: ");
		String s = "";
		for(int i=0;i<resulty.length;i++)
			s += resulty[i] + " ";
		System.out.println(s + "\n");
		
		// prints the μ's
		System.out.println("Set of u:\n");
		for(int i=0;i<resultu.length;i++) {
			System.out.println("Letter " + alphabet[i]);
			for(int j=0;j<resultu[i].length;j++) {
				s = "";
				for(int k=0;k<resultu[i].length;k++)
					s += resultu[i][j][k] + " ";
				System.out.println(s);
			}
			System.out.println();
		}
	}

	public static String genLongTest() {
		String longTest = "";
		
		// adds checkLength random characters in alphabet to longTest
		for(int i=0;i<checkLength;i++)
			longTest += alphabet[(int)(Math.random()*alphabetSize)];
		return longTest;
	}
	
	public static boolean finalCheck() {
		for(int i=0;i<numToCheck;i++) {
			// creates a "long" test and checks whether the hypothesis and target function have the same output
			String longTest = genLongTest();
			if(MQ(longTest)!=MQH(longTest, resulty, resultu, l))
				return false;
		}
		return true;
	}
	
	public static int mod2(int n) {
		if(n%2==0)
			return 0;
		else
			return 1;
	}
	
	
	
	// functions for linear algebra

	public static int determinant(int[][] arr) {
		int result = 0;
		if (arr.length == 1) {
			result = arr[0][0];
			return result;
		}
		if (arr.length == 2) {
			result = arr[0][0] * arr[1][1] - arr[0][1] * arr[1][0];
			return result;
		}
		for (int i = 0; i < arr[0].length; i++) {
			int temp[][] = new int[arr.length - 1][arr[0].length - 1];

			for (int j = 1; j < arr.length; j++) {
				for (int k = 0; k < arr[0].length; k++) {

					if (k < i) {
						temp[j - 1][k] = arr[j][k];
					} else if (k > i) {
						temp[j - 1][k - 1] = arr[j][k];
					}
				}
			}
			result += arr[0][i] * Math.pow(-1, (int) i) * determinant(temp);
		}
		return result;
	}
	
	public static int[][] inverse(int[][] a){
		int n = a.length;
		int x[][] = new int[n][n];
		int b[][] = new int[n][n];
        int index[] = new int[n];
        for (int i=0; i<n; ++i) 
            b[i][i] = 1;
 
        gaussian(a, index);
 
        for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
                for (int k=0; k<n; ++k)
                    b[index[j]][k]
                    	    -= a[index[j]][i]*b[index[i]][k];
 
        for (int i=0; i<n; ++i) 
        {
            x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
            for (int j=n-2; j>=0; --j) 
            {
                x[j][i] = b[index[j]][i];
                for (int k=j+1; k<n; ++k) 
                {
                    x[j][i] -= a[index[j]][k]*x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
        
        for(int i=0;i<n;i++) {
        	for(int j=0;j<n;j++)
        		x[i][j] = mod2(x[i][j]);
        }
 
        return x;
	}
	
	public static void gaussian(int a[][], int index[]) 
    {
        int n = index.length;
        int c[] = new int[n];
 
        for (int i=0; i<n; ++i) 
            index[i] = i;
 
        for (int i=0; i<n; ++i) 
        {
        	int c1 = 0;
            for (int j=0; j<n; ++j) 
            {
            	int c0 = Math.abs(a[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }
 
        int k = 0;
        for (int j=0; j<n-1; ++j) 
        {
        	int pi1 = 0;
            for (int i=j; i<n; ++i) 
            {
            	int pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1) 
                {
                    pi1 = pi0;
                    k = i;
                }
            }
 
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<n; ++i) 	
            {
            	int pj = a[index[i]][j]/a[index[j]][j];
 
                a[index[i]][j] = pj;
 
                for (int l=j+1; l<n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }
	
	public static int[] rowxnxnMatrixMult(int[] arr1, int[][] arr2) {
		int size = arr1.length;
		int[] out = new int[size];
		for(int i=0;i<size;i++) {
			for(int j=0;j<size;j++)
				out[i] += arr1[j]*arr2[j][i];
			out[i] = mod2(out[i]);
		}
		return out;
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
	
	public static int[] nxnxnx1MatrixMult(int[][] arr1, int[] arr2) {
		int size = arr1.length;
		int[] out = new int[size];
		for(int i=0;i<size;i++) {
			for(int j=0;j<size;j++)
				out[i] += arr1[i][j]*arr2[j];
			out[i] = mod2(out[i]);
		}
		return out;
	}
	
	public static int dotProduct(int[] v1, int[] v2) {
		int sum = 0;
		for(int i=0;i<v1.length;i++)
			sum += v1[i]*v2[i];
		return mod2(sum);
	}
	
	public static int[][] rref(int[][] matrix) {
		int[][] rref = new int[matrix.length][];
		for (int i = 0; i < matrix.length; i++)
			rref[i] = Arrays.copyOf(matrix[i], matrix[i].length);

		int r = 0;
		for (int c = 0; c < rref[0].length && r < rref.length; c++) {
			int j = r;
			for (int i = r + 1; i < rref.length; i++)
				if (rref[i][c] > rref[j][c])
					j = i;
			if (rref[j][c] == 0)
				continue;

			int[] temp = rref[j];
			rref[j] = rref[r];
			rref[r] = temp;

			for (int i = 0; i < rref.length; i++) {
				if (i != r) {
					int t = rref[i][c];
					for (j = 0; j < rref[0].length; j++)
						rref[i][j] = mod2(rref[i][j] - (t * rref[r][j]));
				}
			}
			r++;
		}

		return rref;
	}
	
	public static boolean linInd(int[] w, ArrayList<int[]> B) {
		// forms augmented matrix B|w
		int numRows = w.length;
		int numCols = B.size()+1;
		int[][] arr = new int[numRows][numCols];
		for(int j=0;j<numCols-1;j++) {
			for(int i=0;i<numRows;i++)
				arr[i][j] = B.get(j)[i];
		}
		for(int i=0;i<numRows;i++)
			arr[i][numCols-1] = w[i];
		
		int[][] rrefArr = rref(arr);
		
		// finds the index of the last 1 in the last column (if exists)
		int index = -1;
		for(int i=numRows-1;i>=0;i--) {
			if(rrefArr[i][numCols-1] == 1) {
				index = i;
				break;
			}
		}
		
		// last vector is the 0 vector, in span(B)
		if(index == -1)
			return false;
		
		// checks whether in span
		for(int j=0;j<numCols-1;j++) {
			if(rrefArr[index][j] == 1)
				return false;
		}
		
		// linearly independent
		return true;
	}

}