/*
 * Creator: Nevin George
 * Advisor: Dana Angluin
 * Program Description: The algorithm takes in as input a mod-2-MA and prints to stdout the MA
 * obtained after learning the input function through a series of membership and equivalence queries.
 * The algorithm is explained in more detail in the paper "Learning functions represented as multiplicity 
 * automata" by Beimel et al. The motivation behind this algorithm originally arose from Angluin's 
 * exact learning model described in her paper "Learning regular sets from queries and counterexamples."
 * 
 * References:
 * 1 Amos Beimel, Francesco Bergadano, Nader H. Bshouty, Eyal Kushilevitz, Stefano Varric- chio. Learning 
 *   functions represented as multiplicity automata. J. ACM, 47(3):506–530, May 2000.
 * 2 Dana Angluin. Learning regular sets from queries and counterexamples. Inf. Comput., 75(2):87–106, 1987.
 * 3 Dana Angluin, Timos Antonopoulos, Dana Fisman. Strongly Unambiguous Büchi Automata Are Polynomially 
 *   Predictable with Membership Queries. 28th International Conference on Computer Science Logic, 8:1–8:17, 2020.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

public class Mod2_MA {
	
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
	
	// boolean checking whether the algorithm terminates
	public static boolean fail;
	
	public static void main(String[] args) throws IOException {
		initialize();
		
		// contain the indices xi and yi
		ArrayList<String> X = new ArrayList<String>();
		ArrayList<String> Y = new ArrayList<String>();
		X.add("");
		Y.add("");
		
		/* f("") cannot equal 0 (otherwise can't form a linearly independent basis of elements in X).
		 * Start the algorithm instead with a 2x2 matrix of full rank.
		 */
		if(MQ("")==0) {
			int[] hy = createHY(l, X);
			int[][][] setOfHu = createHU(X, Y, l);
			// generate a counter-example z
			if(!EQ(hy, setOfHu, l)) {
				X.add(z);
				Y.add(z);
				l++;
			}
		}
		
		// run the algorithm
		learnSUBA(X, Y);
		if(fail) {
			System.out.println("The algorithm did not stop.");
			return;
		}
		
		// additional final check, checks equivalence of 100 "long" randomly generated strings
		if(finalCheck())
			displayResults();
		else
			System.out.println("Failed final check.");
	}
	
	public static void initialize() throws IOException{
		/* input file must have the form
		 * <alphabetSize>
		 * <characters in alphabet (space separated)>
		 * <r>
		 * <fy (space separated)>
		 * list of μ's for each character in the alphabet, each μ appears in a rxr grid (space separated)
		 * no line separation between any of the above
		 */
		
		// read in from input file
		BufferedReader f = new BufferedReader(new FileReader("input3.txt"));
		
		// read in the alphabet
		alphabetSize = Integer.parseInt(f.readLine());
		StringTokenizer st = new StringTokenizer(f.readLine());
		alphabet = new Character[alphabetSize];
		for(int i=0;i<alphabetSize;i++)
			alphabet[i] = st.nextToken().charAt(0);
		
		// size of the target function
		r = Integer.parseInt(f.readLine());
		
		// γ of the target function
		fy = new int[r];
		st = new StringTokenizer(f.readLine());
		for(int i=0;i<r;i++)
			fy[i] = Integer.parseInt(st.nextToken());
		
		// set of μ's for the target function
		fu = new int[alphabet.length][r][r];
		for(int i=0;i<alphabetSize;i++) {
			for(int j=0;j<r;j++) {
				st = new StringTokenizer(f.readLine());
				for(int k=0;k<r;k++)
					fu[i][j][k] = Integer.parseInt(st.nextToken());
			}
		}
		f.close();
		
		// general variables to initialize
		
		fail = false;
		checkLength = 100;
		
		// initial size of X and Y
		l = 1;
		
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
	
	public static void learnSUBA(ArrayList<String> X, ArrayList<String> Y) {
		// creates the γ for the hypothesis
		int[] hy = createHY(l, X);
		// creates the set of μ's for the hypothesis
		int[][][] setOfHu = createHU(X, Y, l);
		
		// sees if the hypothesis = target function, if so returns the hypothesis
		if(EQ(hy, setOfHu, l)) {
			resulty = hy;
			resultu = setOfHu;
			return;
		}
		
		w = "";
		sig = 0;
		y = "";
		// tries to calculate ω, σ, and y, if cannot find a y that works throws an exception
		try {
			calcWSigY(l, setOfHu, X, Y);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		// algorithm failed
		if(l==r) {
			fail = true;
			return;
		}
		
		// updates l, X, and Y for next recursive call
		l++;
		X.add(w);
		Y.add(sig+y);
		learnSUBA(X, Y);
	}
	
	public static int[] createHY(int l, ArrayList<String> X) {
		// γ is the set of membership queries of the indices in X
		int[] y = new int[l];
		for(int i=0;i<l;i++) 
			y[i] = MQ(X.get(i));
		return y;
	}
	
	public static int[][][] createHU(ArrayList<String> X, ArrayList<String> Y, int l) {
		/*
		 * For every s, define a matrix μ by letting its i-th row be the coefficients of the vector F_{xi+σ}(y) when expressed 
		 * as a linear combination of the vectors F_x1 to F_xl (such coefficients exist as F_x1 to F_xl are linearly independent).
		*/
		int[][][] setOfU = new int[alphabetSize][l][l];
		for(int c=0;c<alphabetSize;c++) {
			// calculuate the μ for character c in the alphabet
			char sig = alphabet[c];
			int[][] u = new int[l][l];
			
			// calculate the vectors F_x1 to F_xl
			int[][] rowsOfFxi = new int[l][l];
			for(int j=0;j<l;j++) {
				int[] Fxi = new int[l];
				for(int k=0;k<l;k++)
					Fxi[k] = MQ(X.get(j) + Y.get(k));
				rowsOfFxi[j] = Fxi;
			}
			
			int[][] inverseRoF;
			// if rowsOfFxi is not invertible, let inverseRoF be all 0's
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
				// find F_{xi+σ}
				int[] Fxisig = new int[l];
				for(int j=0;j<l;j++)
					Fxisig[j] = MQ(X.get(i) + sig + Y.get(j));
				
				// calculate μ
				u[i] = rowxnxnMatrixMult(Fxisig, inverseRoF);
			}
			
			setOfU[c] = u;
		}
		return setOfU;
	}
	
	public static int MQ(String w) {
		// MQ for the target function
		
		// initialize cur as the rxr identity matrix
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
		// checks to see if the hypothesis and target functions have the same output for all words in tests
		for(int i=0;i<maxTestSize;i++) {
			if(MQ(tests[i])!=MQH(tests[i], hy, hu, l)) {
				z = tests[i];
				return false;
			}
		}
		return true;
	}
	
	public static void calcWSigY(int l, int[][][] setOfHu, ArrayList<String> X, ArrayList<String> Y) throws Exception {
		// goes through every possible prefix of z starting with ω = "" and σ = (first character of ω)
		// prefix = ω + σ
		for(int i=0;i<z.length();i++) {
			if(i!=0)
				w = z.substring(0,i);
			sig = z.charAt(i);
			
			// calculate μ(ω)
			int[][] u = new int[l][l];
			for(int m=0;m<l;m++)
				u[m][m] = 1;
			
			for(int n=0;n<w.length();n++)
				u = nxnMatrixMult(u, setOfHu[letterToIndex.get(w.charAt(n))]);
			
			// go through every possible value of y in Y
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
		throw new Exception("Didn't find a suitable omega, sigma, and gamma.");
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
 
        // Transform the matrix into an upper triangle
        gaussian(a, index);
 
        // Update the matrix b[i][j] with the ratios stored
        for (int i=0; i<n-1; ++i)
            for (int j=i+1; j<n; ++j)
                for (int k=0; k<n; ++k)
                    b[index[j]][k]
                    	    -= a[index[j]][i]*b[index[i]][k];
 
        // Perform backward substitutions
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
 
        // Initialize the index
        for (int i=0; i<n; ++i) 
            index[i] = i;
 
        // Find the rescaling factors, one from each row
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
 
        // Search the pivoting element from each column
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
 
            // Interchange rows according to the pivoting order
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<n; ++i) 	
            {
            	int pj = a[index[i]][j]/a[index[j]][j];
 
                // Record pivoting ratios below the diagonal
                a[index[i]][j] = pj;
 
                // Modify other elements accordingly
                for (int l=j+1; l<n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }
	
	public static int[] rowxnxnMatrixMult(int[] arr1, int[][] arr2) {
		int[] out = new int[arr1.length];
		for(int i=0;i<arr1.length;i++) {
			int sum = 0;
			for(int j=0;j<arr1.length;j++)
				sum += arr1[j]*arr2[j][i];
			out[i] = mod2(sum);
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
	
	public static int dotProduct(int[] v1, int[] v2) {
		int sum = 0;
		for(int i=0;i<v1.length;i++)
			sum += v1[i]*v2[i];
		return mod2(sum);
	}
	
	public static int rank(int[][] target) {
		double EPS = 1E-9;
		int n = target.length;
	    int m = target[0].length;

	    int rank = 0;
	    boolean[] row_selected = new boolean[n];;
	    for (int i = 0; i < m; ++i) {
	        int j;
	        for (j = 0; j < n; ++j) {
	            if (!row_selected[j] && Math.abs(target[j][i]) > EPS)
	                break;
	        }

	        if (j != n) {
	            ++rank;
	            row_selected[j] = true;
	            for (int p = i + 1; p < m; ++p)
	                target[j][p] /= target[j][i];
	            for (int k = 0; k < n; ++k) {
	                if (k != j && Math.abs(target[k][i]) > EPS) {
	                    for (int p = i + 1; p < m; ++p)
	                        target[k][p] -= target[j][p] * target[k][i];
	                }
	            }
	        }
	    }
	    return rank;
	}

}