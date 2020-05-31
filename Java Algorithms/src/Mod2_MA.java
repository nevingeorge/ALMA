/*
 * Author: Nevin George
 * Advisor: Dana Angluin
 * Program Description: The program takes in as input a mod-2-MA and prints to stdout the mod-2-MA obtained after 
 * learning the input function through a series of membership and equivalence queries.
 * 
 * References:
 * 1 Amos Beimel, Francesco Bergadano, Nader H. Bshouty, Eyal Kushilevitz, Stefano Varric- chio. Learning 
 *   functions represented as multiplicity automata. J. ACM, 47(3):506–530, May 2000.
 * 2 Dana Angluin. Learning regular sets from queries and counterexamples. Inf. Comput., 75(2):87–106, 1987.
 * 3 Dana Angluin, Timos Antonopoulos, Dana Fis`man. Strongly Unambiguous Büchi Automata Are Polynomially 
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

import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class Mod2_MA {
	
	// if true, display the observation table as it is constructed
	public static boolean verbose;
	
	// alphabet
	public static Character[] alphabet;
	// map each letter in the alphabet to an index
	public static HashMap<Character, Integer> letterToIndex;
	
	// γ of the target function
	public static double[] inputY;
	// set of nxn μ's of the target function
	public static double[][][] inputU;
	// size of the target function
	public static int inputR;
	
	// γ of the minimized target function
	public static double[] minY;
	// set of nxn μ's of the minimized target function
	public static double[][][] minU;
	// size of the minimized target function
	public static int minR;
	// row indices of the observation table of the minimized target function
	public static ArrayList<String> rowIndices;
	// column indices/experiments of the observation table of the minimized target function
	public static ArrayList<String> colIndices;
	
	// Hankel matrix
	public static HashMap<String, Integer> F;
	
	// row indices of the observation table
	public static ArrayList<String> X;
	// column indices of the observation table
	public static ArrayList<String> Y;
	// size of X and Y
	public static int l;
	
	// counter-example
	public static String z;
		
	// γ of the learned function
	public static double[] resultY;
	// set of nxn μ's of the learned function
	public static double[][][] resultU;
	
	// use in EQ to avoid testing the same word
	public static int rowStartIndex;
	public static int colStartIndex;
	
	// take in user input
	public static Scanner in;
	
	public static void main(String[] args) throws Exception {
		// read in the mod-2-MA
		initialize();
		
		// minimize the mod-2-MA
		minimize();
				
		// run the learning algorithm
		run();

		// statistical final check of equivalence
		if(finalCheck(50,40))
			displayResults();
		else
			throwException(null,"Algorithm failed: failed final check.");
		
		// perform operations on the learned mod-2-MA
		operations();
	}
	
	public static void initialize() throws Exception {
		/* The input file must have the following format (no line separation, entries are space separated, 
		 * and lines beginning with // are ignored):
		 * <characters in the alphabet>
		 * <size of the target function (r)>
		 * <γ of the target function (fy)>
		 * List of μ's for each character in the alphabet, each μ appears in a rxr grid
		 * 
		 * Example input files can be found in the GitHub repository.
		 */
		
		// read in file name + optional flag -v from stdin
		System.out.println("Input file name and optional flag -v (e.g. Mod2_MA_input1.txt or Mod2_MA_input1.txt -v)");
		in = new Scanner(System.in);
		String[] arrInput = in.nextLine().split(" ");
		verbose = false;
		if(arrInput.length > 2)
			throwException(null,"Invalid input: too many inputs passed");
		if(arrInput.length == 2) {
			if(arrInput[1].equals("-v"))
				verbose = true;
			else
				throwException(null,"Invalid input: invalid flag");
		}
		BufferedReader f = new BufferedReader(new FileReader(arrInput[0]));
		System.out.println();
		
		// alphabet
		StringTokenizer st = new StringTokenizer(readInput(f));
		ArrayList<Character> tempAlphabet = new ArrayList<Character>();
		while(st.hasMoreTokens()) {
			String letter = st.nextToken();
			if(letter.length()!=1)
				throwException(f,"Invalid input: invalid character in the alphabet");
			tempAlphabet.add(letter.charAt(0));
		}
		alphabet = new Character[tempAlphabet.size()];
		for(int i=0;i<tempAlphabet.size();i++)
			alphabet[i] = tempAlphabet.get(i);
		
		// map each letter in the alphabet to an index
		letterToIndex = new HashMap<Character, Integer>();
		for(int i=0;i<alphabet.length;i++)
			letterToIndex.put(alphabet[i], i);
		
		// size of the target function
		inputR = Integer.parseInt(readInput(f));
		
		// γ of the target function
		st = new StringTokenizer(readInput(f));
		inputY = new double[inputR];
		for(int i=0;i<inputR;i++)
			inputY[i] = Integer.parseInt(st.nextToken());
		if(st.hasMoreTokens())
			throwException(f,"Invalid input: γ length exceeds the specified size");
		
		// set of μ's for the target function
		inputU = new double[alphabet.length][inputR][inputR];
		for(int i=0;i<alphabet.length;i++) {
			for(int j=0;j<inputR;j++) {
				st = new StringTokenizer(readInput(f));
				for(int k=0;k<inputR;k++)
					inputU[i][j][k] = Integer.parseInt(st.nextToken());
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
		// ignore lines beginning with "//"
		while(line != null && line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();
		return line;
	}
	
	// properly close input streams
	public static void throwException(BufferedReader f, String message) throws Exception {
		if(f != null)
			f.close();
		if(in != null)
			in.close();
		throw new Exception(message);
	}
	
	// follow an adapted version of algorithm 2 in Thon and Jaeger to minimize the target function
	public static void minimize() throws OutOfRangeException, Exception {
		// basis for the state space
		ArrayList<String> stateIndices = new ArrayList<String>();
		RealMatrix B_1 = basis(inputY, inputU, stateIndices, true);
		// basis for the co-state space
		ArrayList<String> co_stateIndices = new ArrayList<String>();
		RealMatrix B_2 = basis(inputY, inputU, co_stateIndices, false);
		
		// observation table B_1 x B_2
		RealMatrix T_1 = MatrixUtils.createRealMatrix(B_1.getRowDimension(), B_2.getRowDimension());
		for(int i=0;i<B_1.getRowDimension();i++) {
			for(int j=0;j<B_2.getRowDimension();j++)
				T_1.setEntry(i, j, B_1.getRowVector(i).dotProduct(B_2.getRowVector(j)));
		}
		
		// obtain the smallest set of linearly independent rows and columns from T
		rowIndices = new ArrayList<String>();
		RealMatrix T_2 = basisT(T_1, stateIndices, rowIndices, true);
		colIndices = new ArrayList<String>();
		RealMatrix T_3 = basisT(T_2, co_stateIndices, colIndices, false);
		
		// case where T_3 = [[0]] (T_3 is singular, must be treated separately)
		if(T_3.getRowDimension()==1 && T_3.getEntry(0, 0)==0) {
			minR = 1;
			minY = new double[1];
			minY[0] = 0;
			minU = new double[alphabet.length][1][1];
			for(int i=0;i<alphabet.length;i++)
				minU[i][0][0] = 1;
			return;
		}
		
		DecompositionSolver solver = new solver(T_3).getSolver();
		RealMatrix Tinverse = solver.getInverse();
		
		// initialize the Hankel matrix
		F = new HashMap<String, Integer>();
		
		// minU = xSigma*Tinverse, where xSigma is the matrix where row_i = row of the observation table indexed by x_i+σ
		// temporarily set the minimized values to the input values because MQ relies on the minimized values
		minR = inputR;
		minY = inputY;
		minU = inputU;
		double[][][] tempMinU = new double[alphabet.length][T_3.getRowDimension()][T_3.getRowDimension()];
		for(int i=0;i<alphabet.length;i++) {
			RealMatrix xSigma = MatrixUtils.createRealMatrix(T_3.getRowDimension(), T_3.getRowDimension());
			for(int j=0;j<T_3.getRowDimension();j++) {
				for(int k=0;k<T_3.getRowDimension();k++)
					xSigma.setEntry(j, k, MQ(rowIndices.get(j)+alphabet[i]+colIndices.get(k)));
			}
			
			tempMinU[i] = xSigma.multiply(Tinverse).getData();
		}
		
		minU = new double[alphabet.length][T_3.getRowDimension()][T_3.getRowDimension()];
		for(int i=0;i<alphabet.length;i++) {
			for(int j=0;j<T_3.getRowDimension();j++) {
				for(int k=0;k<T_3.getRowDimension();k++)
					minU[i][j][k] = mod2(tempMinU[i][j][k]);
			}
		}
		
		minR = T_3.getRowDimension();
		
		// γ = T_3*e_1, where e_1 is a standard unit basis vector
		double[] temp = new double[minR];
		temp[0] = 1;
		RealVector e_1 = MatrixUtils.createRealVector(temp);
		minY = T_3.operate(e_1).toArray();
		for(int i=0;i<minR;i++)
			minY[i] = mod2(minY[i]);
		
		// use in EQ to avoid testing the same word
		rowStartIndex = 0;
		colStartIndex = 0;
	}
	
	// follow algorithm 1 detailed in Thon and Jaeger to form the basis for the state/co-state space
	public static RealMatrix basis(double[] hy, double[][][] hu, ArrayList<String> indices, boolean stateSpace) {
		// basis
		double[][] B = new double[hy.length][hy.length];
		int sizeB = 0;
		
		// set with elements to try to add to B
		ArrayList<double[]> C = new ArrayList<double[]>();
		// begin with ω_i = (1,0,0,...,0)
		if(stateSpace) {
			double[] w_i = new double[hy.length];
			w_i[0] = 1;
			C.add(w_i);
		}
		// begin with γ
		else
			C.add(hy);
		int sizeC = 1;
		
		// contain the corresponding ω for every element in C
		ArrayList<String> WC = new ArrayList<String>();
		WC.add("");
		
		while(sizeC>0) {
			// element to test
			double[] w = C.remove(0);
			String s = WC.remove(0);
			sizeC--;
			
			// test if ω is linearly independent of B
			if(linInd(w, B, sizeB)) {	
				// extend B
				B[sizeB++] = w;
				indices.add(s);
				
				// add {ωμ(σ) | σ∈Σ} to C
				for(int i=0;i<alphabet.length;i++) {
					RealMatrix m = MatrixUtils.createRealMatrix(hu[i]);
					RealVector p = MatrixUtils.createRealVector(w);
					double[] v;
					// basis for the set span(ω_iμ(ω) : ω∈Σ*)
					if(stateSpace) {
						v = m.preMultiply(p).toArray();
						WC.add(s+alphabet[i]);
					}
					// basis for the set span(μ(ω)γ : ω∈Σ*)
					else {
						v = m.operate(p).toArray();
						WC.add(alphabet[i]+s);
					}
					for(int j=0;j<v.length;j++)
						v[j] = mod2(v[j]);
					C.add(v);
					sizeC++;
				}
			}
		}
		
		return MatrixUtils.createRealMatrix(B).getSubMatrix(0, sizeB-1, 0, hy.length-1);
	}
	
	// find a maximal subset of linearly independent rows/columns of T
	public static RealMatrix basisT(RealMatrix T, ArrayList<String> oldIndices, ArrayList<String> newIndices, boolean rows) {	
		if(!rows)
			T = T.transpose();
		
		double[][] newT = new double[T.getRowDimension()][T.getColumnDimension()];
		int sizeT = 0;
		
		for(int i=0;i<T.getRowDimension();i++) {
			// extend the current subset of linearly independent rows/columns
			if(linInd(T.getRow(i), newT, sizeT)) {
				newT[sizeT++] = T.getRow(i);
				newIndices.add(oldIndices.get(i));
			}
		}
		
		if(rows)
			return MatrixUtils.createRealMatrix(newT).getSubMatrix(0, sizeT-1, 0, T.getColumnDimension()-1);
		else
			return MatrixUtils.createRealMatrix(newT).getSubMatrix(0, sizeT-1, 0, T.getColumnDimension()-1).transpose();
	}
	
	// test whether the vector w is in the span of the set B
	public static boolean linInd(double[] w, double[][] B, int sizeB) {
		if(sizeB==0)
			return true;
		
		// form the augmented matrix B|w
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
				if(mod2(m.getEntry(i, c))>mod2(m.getEntry(j, c)))
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
		
		// find the index of the last 1 in the last column (if exists)
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
		
		// check whether in span
		for(int j=0;j<numCols-1;j++) {
			if(mod2(m.getEntry(index, j)) == 1)
				return false;
		}
		
		// linearly independent
		return true;
	}
	
	public static void run() throws Exception {	
		// initialize the rows and columns of the observation table
		X = new ArrayList<String>();
		Y = new ArrayList<String>();
		X.add("");
		Y.add("");
		l = 1;
		
		// initialize the Hankel matrix
		if(F==null)
			F = new HashMap<String, Integer>();
		
		/* f("") cannot equal 0 (otherwise can't form a linearly independent basis of elements in X).
		 * The algorithm instead begins with a 2x2 matrix of full rank.
		 */
		if(MQ("")==0) {
			double[] hy = createHY();
			double[][][] hu = createHU();
			
			// generate a counter-example z
			if(!EQ(hy, hu)) {
				l++;
				X.add(z);
				Y.add(z);
			}
		}
		
		if(verbose) {
			System.out.println("Results after individual queries\n--------------------------------");
			displayTable();
		}
		
		// run the algorithm
		learnMA();
	}
	
	public static void learnMA() throws Exception {
		// create the γ for the hypothesis
		double[] hy = createHY();
		// create the set of μ's for the hypothesis
		double[][][] hu = createHU();
		
		// see if the hypothesis = target function, if so return the hypothesis
		if(EQ(hy, hu)) {
			resultY = hy;
			resultU = hu;
			return;
		}
		
		// attempt to calculate ω, σ, and y
		// if there is no y that works throw an exception
		calcWSigY(hu);
		
		learnMA();
	}
	
	public static double[] createHY() throws Exception {
		// γ is the set of results obtained after performing membership queries on the indices in X
		double[] y = new double[l];
		for(int i=0;i<l;i++)
			y[i] = MQ(X.get(i));
		return y;
	}
	
	public static double[][][] createHU() throws Exception {
		/*
		 * For every s, define a matrix μ by letting its i-th row be the coefficients of the vector F_{xi+σ}(y) 
		 * when expressed as a linear combination of the vectors F_x1 to F_xl (such coefficients exist as F_x1 
		 * to F_xl are linearly independent).
		*/
		
		double[][][] hu = new double[alphabet.length][l][l];
		for(int c=0;c<alphabet.length;c++) {
			// calculuate μ_c
			char sigma = alphabet[c];
			
			// calculate the vectors F_xi and F_{xi+σ}
			double[][] F_xi = new double[l][l];
			double[][] F_xi_sigma = new double[l][l];
			for(int i=0;i<l;i++) {
				for(int j=0;j<l;j++) {
					F_xi[j][i] = MQ(X.get(i) + Y.get(j));
					F_xi_sigma[i][j] = MQ(X.get(i) + sigma + Y.get(j));
				}
			}
			
			// solve the matrix equation using LU Decomposition
			RealMatrix coefficients = new Array2DRowRealMatrix(F_xi);
			DecompositionSolver solver = new solver(coefficients).getSolver();
			for(int i=0;i<l;i++) {
				RealVector constants = new ArrayRealVector(F_xi_sigma[i]);
				try {
					RealVector solution = solver.solve(constants);
					for(int j=0;j<l;j++)
						hu[c][i][j] = solution.getEntry(j);
				}
				// matrix is not invertible
				catch(Exception e) {
					for(int j=0;j<l;j++)
						hu[c][i][j] = 0;
				}
			}
		}
		return hu;
	}
	
	// MQ for the target function
	public static int MQ(String w) throws Exception {		
		// MQ(ω) was previously calculated and is in the Hankel matrix
		if(F.get(w) != null)
			return F.get(w);
		
		int out = 0;
		
		// use a membership query function defined in NBA.java
		if(NBA.F!=null)
			out = NBA.MQ(w);
		
		// use a membership query function defined in MQ.java
		else if(arbitrary.MQarbitrary!=null) {
			try {
				out = (int) arbitrary.MQarbitrary.invoke(null,w);
			} catch (Exception e) {
				throwException(null, "Invalid input: invalid membership query function");
			} 
		}
		
		// perform the membership query on the target mod-2-MA
		else {
			// initialize cur as the rxr identity matrix
			RealMatrix cur = MatrixUtils.createRealIdentityMatrix(minR);
			
			// multiply cur by the corresponding μ for each letter in ω
			for(int i=0;i<w.length();i++) {
				cur = cur.multiply(MatrixUtils.createRealMatrix(minU[letterToIndex.get(w.charAt(i))]));
				// account for rounding issues with larger words
				if(i!=0 && i%25==0) {
					for(int j=0;j<minR;j++) {
						for(int k=0;k<minR;k++)
							cur.setEntry(j, k, mod2(cur.getEntry(j, k)));
					}
				}
			}
			
			// multiply the final result with γ
			out = mod2(cur.getRowVector(0).dotProduct(new ArrayRealVector(minY)));
		}
		
		// add MQ(ω) to the Hankel matrix
		F.put(w, out);
		
		return out;
	}
	
	// MQ for the current hypothesis
	public static int MQH(double[] hy, double[][][] hu, String w) {		
		// initialize cur as the lxl identity matrix
		RealMatrix cur = MatrixUtils.createRealIdentityMatrix(hy.length);
		
		// multiply cur by the corresponding μ for each letter in ω
		for(int i=0;i<w.length();i++) {
			cur = cur.multiply(MatrixUtils.createRealMatrix(hu[letterToIndex.get(w.charAt(i))]));
			// account for rounding issues with larger words
			if(i!=0 && i%25==0) {
				for(int j=0;j<hy.length;j++) {
					for(int k=0;k<hy.length;k++)
						cur.setEntry(j, k, mod2(cur.getEntry(j, k)));
				}
			}
		}
		
		// multiply the final result with γ
		return mod2(cur.getRowVector(0).dotProduct(new ArrayRealVector(hy)));
	}
	
	public static boolean EQ(double[] hy, double[][][] hu) throws Exception {
		// use a statistical equivalence query
		if(NBA.F!=null || arbitrary.MQarbitrary!=null)
			return arbitrary.EQapprox(hy, hu);
		
		// test every element in T_3, the observation table for the minimized target function
		for(int i=rowStartIndex;i<rowIndices.size();i++) {
			for(int j=colStartIndex;j<colIndices.size();j++) {
				String test = rowIndices.get(i)+colIndices.get(j);
				if(MQ(test)!=MQH(hy, hu, test)) {
					// update rowStartIndex and colStartIndex to avoid testing the same words in the next EQ
					if(j==colIndices.size()-1) {
						rowStartIndex = i+1;
						colStartIndex = 0;
					}
					else {
						rowStartIndex = i;
						colStartIndex = j+1;
					}
					
					z = test;
					return false;
				}
			}
		}
		
		return true;
	}
	
	public static void calcWSigY(double[][][] hu) throws Exception {
		// prefix of z = ω + σ
		String w = "";
		char sigma = 0;
		// experiment
		String y = "";
		
		// go through every possible prefix of z starting with ω = "" and σ = (first character of ω)
		for(int i=0;i<z.length();i++) {
			if(i!=0)
				w = z.substring(0,i);
			sigma = z.charAt(i);
			
			// calculate μ(ω)
			RealMatrix mu = MatrixUtils.createRealIdentityMatrix(l);
			for(int n=0;n<w.length();n++)
				mu = mu.multiply(MatrixUtils.createRealMatrix(hu[letterToIndex.get(w.charAt(n))]));
			
			// check if F_ω = sum(μ(ω)_1,i * F_xi)
			boolean failed = false;
			for(int j=0;j<l;j++) {
				int sum = 0;
				for(int k=0;k<l;k++)
					sum = mod2(sum + mu.getEntry(0, k)*MQ(X.get(k)+Y.get(j)));
				
				if(MQ(w+Y.get(j)) != sum) {
					failed = true;
					break;
				}
			}
			if(failed)
				continue;
			
			// go through every possible value of y in Y
			// check if F_{ω+σ}(y) != sum(μ(ω)_1,i * F_{xi+σ}(y))
			for(int j=0;j<l;j++) {
				y = Y.get(j);
			
				int sum = 0;
				for(int k=0;k<l;k++)
					sum = mod2(sum + mu.getEntry(0, k)*MQ(X.get(k) + sigma + y));
				
				// found a solution
				if(MQ(w+sigma+y) != sum) {		
					if(l==minR)
						throwException(null,"Algorithm failed: size of the hypothesis exceeds that of the target function.");
					// update l, X, and Y
					l++;
					X.add(w);
					Y.add(sigma+y);
					
					// display the updated observation table
					if(verbose)
						displayTable();
					
					return;
				}
			}
		}

		throwException(null,"Algorithm failed: didn't find a suitable omega, sigma, and gamma.");
	}

	public static void displayResults() {
		System.out.println("Learned mod-2-MA");
		System.out.println("----------------");
		// print γ
		System.out.print("y: ");
		String s = "";
		for(int i=0;i<resultY.length;i++)
			s += mod2(resultY[i]) + " ";
		System.out.println(s + "\n");
		
		// print the μ's
		System.out.println("Set of u:\n");
		for(int i=0;i<resultU.length;i++) {
			System.out.println("Letter " + alphabet[i]);
			for(int j=0;j<resultU[i].length;j++) {
				s = "";
				for(int k=0;k<resultU[i].length;k++)
					s += mod2(resultU[i][j][k]) + " ";
				System.out.println(s);
			}
			System.out.println();
		}
	}
	
	public static void displayTable() throws Exception {
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
		// add len number of random characters in alphabet to test
		String test = "";
		for(int i=0;i<len;i++)
			test += alphabet[(int)(Math.random()*alphabet.length)];
		return test;
	}

	// perform a statistical equivalence query between the target and learned mod-2-MA's as a final safety check
	public static boolean finalCheck(int maxTestLen, int numTests) throws Exception {
		// create numTests tests of length at most maxTestLen
		for(int i=1;i<=numTests;i++) {
			String test = genTest((int)(Math.random()*(maxTestLen+1)));
			if(MQH(inputY,inputU,test)!=MQH(resultY, resultU, test))
				return false;
		}
		return true;
	}
	
	// perform operations on the learned mod-2-MA
	public static void operations() {
		System.out.println("Available operations for the learned Mod-2-MA (enter \"quit\" to terminate):");
		System.out.println("- Test whether a word is accepted: \"-a <word>\"\n  If the language is (L)_$, words must be of the form u$v.");
		while(true) {
			// read in cmd
			String line = in.nextLine();
			
			// terminate the program
			if(line.equals("quit")) {
				in.close();
				break;
			}
			
			String[] input = line.split(" ");
			if(input.length == 0)
				System.out.println("Invalid cmd");
			
			// test whether a word is accepted
			if(input[0].equals("-a")) {
				if(input.length>2)
					System.out.println("Invalid cmd");
				else {
					// if no word is passed, the test is the empty string
					String test = "";
					if(input.length==2)
						test = input[1];
					
					if(!inAlphabet(test))
						System.out.println("Inputted word is not in the language.");
					else if(MQH(resultY, resultU, test) == 1)
						System.out.println("Accepted");
					else
						System.out.println("Not accepted");
				}
			}
			else
				System.out.println("Invalid cmd");
		}
	}
	
	public static boolean inAlphabet(String word) {
		for(int i=0;i<word.length();i++) {
			if(letterToIndex.get(word.charAt(i)) == null)
				return false;
		}
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