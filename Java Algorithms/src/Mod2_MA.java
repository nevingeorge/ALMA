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

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class Mod2_MA {
	
	// if true, displays the observation table as it is being constructed
	public static boolean verbose;
	
	// alphabet
	public static Character[] alphabet;
	// maps each letter in the alphabet to an index
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
	
	// takes in user input
	public static Scanner in;
	
	public static void main(String[] args) throws Exception {
		// reads in the mod-2-MA
		initialize();
		
		// minimizes the mod-2-MA
		minimize();
				
		// runs the learning algorithm
		run();
		
		// statistical final check of equivalence
		if(finalCheck(50,40))
			displayResults();
		else
			throwException(null,"Algorithm failed: failed final check.");
		
		// performs desired operations with the learned mod-2-MA
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
		
		// reads in file name + optional flag -v from stdin
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
		
		// maps each letter in the alphabet to an index
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
	
	public static void minimize() {
		minR = inputR;
		minY = inputY;
		minU = inputU;
		
		// TODO
		
		/*
		// follows algorithm 2 detailed in Thon and Jaeger
		
		// basis for the state space
		RealMatrix phi = basis1(inputY, inputU);
		// pseudo-inverse
		DecompositionSolver solver1 = new mod2LUDecomposition(phi.multiply(phi.transpose())).getSolver();
		RealMatrix phiInverse = phi.transpose().multiply(solver1.getInverse());
		
		// (phi-inverse)*M*(phi)
		RealVector sig2 = phi.getRowVector(0);
		RealMatrix[] hu2 = new RealMatrix[inputU.length];
		for(int i=0;i<inputU.length;i++) {
			RealMatrix u = MatrixUtils.createRealMatrix(inputU[i]);
			hu2[i] = phiInverse.multiply(u).multiply(phi);
		}
		RealVector y2 = phiInverse.operate(MatrixUtils.createRealVector(inputY));
		
		// basis for the co-state space
		RealMatrix pi = basis2(sig2, hu2);
		// pseudo-inverse
		DecompositionSolver solver2 = new mod2LUDecomposition(pi.transpose().multiply(pi)).getSolver();
		RealMatrix piInverse = solver2.getInverse().multiply(pi.transpose());
		
		
		// (pi)*M'*(pi-inverse)
		RealVector sig3 = piInverse.preMultiply(sig2);
		RealMatrix[] hu3 = new RealMatrix[inputU.length];
		for(int i=0;i<inputU.length;i++)
			hu3[i] = pi.multiply(hu2[i]).multiply(piInverse);
		RealVector y3 = pi.operate(y2);
		*/
		
		/*
		// basis for the state space
		RealMatrix phi = basis1(inputY, inputU);
		// pseudo-inverse
		RealMatrix phiInverse = (new solver(phi).getSolver()).getInverse();
		
		// (phi-inverse)*M*(phi)
		RealVector sig2 = phi.getRowVector(0);
		RealMatrix[] hu2 = new RealMatrix[inputU.length];
		for(int i=0;i<inputU.length;i++) {
			RealMatrix u = MatrixUtils.createRealMatrix(inputU[i]);
			hu2[i] = phiInverse.multiply(u).multiply(phi);
		}
		RealVector y2 = phiInverse.operate(MatrixUtils.createRealVector(inputY));
		
		// basis for the co-state space
		RealMatrix pi = basis2(sig2, hu2);
		// pseudo-inverse
		RealMatrix piInverse = (new solver(pi.transpose()).getSolver()).getInverse().transpose();
		
		// (pi)*M'*(pi-inverse)
		RealVector sig3 = piInverse.preMultiply(sig2);
		RealMatrix[] hu3 = new RealMatrix[inputU.length];
		for(int i=0;i<inputU.length;i++)
			hu3[i] = pi.multiply(hu2[i]).multiply(piInverse);
		RealVector y3 = pi.operate(y2);
		
		// size of the minimized mod-2-MA
		minR = sig3.getDimension();
		
		// checks if σ = (1,0,0,...,0)
		boolean correct = true;
		for(int i=0;i<minR;i++) {
			if((i==0 && sig3.getEntry(0)==0) || (i!=0 && sig3.getEntry(0)==1)) {
				correct = false;
				break;
			}
		}
		// need to change σ to (1,0,0,...,0)
		if(!correct) {
			// basis for (F_2)^minR satisfying σ*col_0(rho)=1 and σ*col_i(rho)=0 for 1<=i<minR
			RealMatrix rho = MatrixUtils.createRealMatrix(minR, minR);
			
			boolean first = true;
			int curCol = 1;
			int firstOnePos = -1;
			for(int i=0;i<minR;i++) {
				if(sig3.getEntry(i) == 0) {
					// adds the standard basis vector e_i to rho
					for(int j=0;j<minR;j++) {
						if(j==i)
							rho.setEntry(j, curCol, 1);
						else
							rho.setEntry(j, curCol, 0);
					}
					curCol++;
				}
				else if(first) {
					// sets the first column of rho to the standard basis vector e_i
					for(int j=0;j<minR;j++) {
						if(j==i)
							rho.setEntry(j, 0, 1);
						else
							rho.setEntry(j, 0, 0);
					}
					first = false;
					firstOnePos = i;
				}
				else {
					// adds to rho the vector with all 0's except for a 1 at indices firstOnePos and i
					for(int j=0;j<minR;j++) {
						if(j==firstOnePos || j==i)
							rho.setEntry(j, curCol, 1);
						else
							rho.setEntry(j, curCol, 0);
					}
					curCol++;
				}
			}
			
			// pseudo-inverse
			RealMatrix rhoInverse = (new solver(rho).getSolver()).getInverse();
			
			// (rho)*M''*(rho-inverse)
			for(int i=0;i<inputU.length;i++) {
				RealMatrix u = rhoInverse.multiply(hu3[i]).multiply(rho);
				minU[i] = new double[u.getColumnDimension()][u.getColumnDimension()];
				for(int j=0;j<u.getColumnDimension();j++) {
					for(int k=0;k<u.getColumnDimension();k++)
						minU[i][j][k] = u.getEntry(j, k);
				}
			}
			minY = rhoInverse.operate(y3).toArray();
		}
		else {
			for(int i=0;i<inputU.length;i++) {
				int size = hu3[0].getColumnDimension();
				minU[i] = new double[size][size];
				for(int j=0;j<size;j++) {
					for(int k=0;k<size;k++)
						minU[i][j][k] = hu3[i].getEntry(j, k);
				}
			}
			minY = y3.toArray();
		}
		*/
	}
	
	public static RealMatrix basis1(double[] hy, double[][][] hu) {
		// To form the basis, we will follow algorithm 1 detailed in the paper by Thon and Jaeger.
		// basis for the set span(μ(ω)γ : ω∈Σ*)
		double[][] B = new double[hy.length][hy.length];
		int sizeB = 0;
		
		// Set with elements to try to add to B, begin with y
		ArrayList<double[]> C = new ArrayList<double[]>();
		C.add(hy);
		int sizeC = 1;
		
		while(sizeC>0) {
			// element to test
			double[] w = C.remove(0);
			sizeC--;
			
			// tests if ω is linearly independent of B
			if(linInd(w, B, sizeB)) {	
				// extends B
				B[sizeB++] = w;
				
				// adds {μ(σ)ω | σ∈Σ} to C
				for(int i=0;i<alphabet.length;i++) {
					RealMatrix m = MatrixUtils.createRealMatrix(hu[i]);
					RealVector p = MatrixUtils.createRealVector(w);
					double[] v = m.operate(p).toArray();
					for(int j=0;j<v.length;j++)
						v[j] = mod2(v[j]);
					C.add(v);
				}
			}
		}
		
		return MatrixUtils.createRealMatrix(B).transpose().getSubMatrix(0, hy.length-1, 0, sizeB-1);
	}
	
	public static RealMatrix basis2(RealVector sig2, RealMatrix[] hu2) {
		// To form the basis, we will follow algorithm 1 detailed in the paper by Thon and Jaeger.
		// basis for the set span(σ'μ'(ω) : ω∈Σ*)
		double[][] B = new double[sig2.getDimension()][sig2.getDimension()];
		int sizeB = 0;
		
		// Set with elements to try to add to B, begin with σ'
		ArrayList<double[]> C = new ArrayList<double[]>();
		C.add(sig2.toArray());
		int sizeC = 1;
		
		while(sizeC>0) {
			// element to test
			double[] w = C.remove(0);
			sizeC--;
			
			// tests if ω is linearly independent of B
			if(linInd(w, B, sizeB)) {	
				// extends B
				B[sizeB++] = w;
				
				// adds {ωμ'(σ) | σ∈Σ} to C
				for(int i=0;i<alphabet.length;i++) {
					RealVector p = MatrixUtils.createRealVector(w);
					double[] v = hu2[i].preMultiply(p).toArray();
					for(int j=0;j<v.length;j++)
						v[j] = mod2(v[j]);
					C.add(v);
				}
			}
		}
		
		return MatrixUtils.createRealMatrix(B).getSubMatrix(0, sizeB-1, 0, sig2.getDimension()-1);
	}
	
	public static void run() throws Exception {	
		// initializes the rows and columns of the observation table
		X = new ArrayList<String>();
		Y = new ArrayList<String>();
		X.add("");
		Y.add("");
		l = 1;
		
		// initializes the Hankel matrix
		F = new HashMap<String, Integer>();
		
		/* f("") cannot equal 0 (otherwise can't form a linearly independent basis of elements in X).
		 * The algorithm instead begins with a 2x2 matrix of full rank.
		 */
		if(MQ("")==0) {
			double[] hy = createHY();
			double[][][] hu = createHU();
			
			// generates a counter-example z
			if(!EQ(hy, hu)) {
				l++;
				X.add(z);
				Y.add(z);
			}
		}
		
		if(verbose) {
			System.out.println("Results after individual queries\n--------------------------------");
			displayQueries();
		}
		
		// runs the algorithm
		learnMA();
	}
	
	public static void learnMA() throws Exception {
		// creates the γ for the hypothesis
		double[] hy = createHY();
		// creates the set of μ's for the hypothesis
		double[][][] hu = createHU();
		
		// sees if the hypothesis = target function, if so returns the hypothesis
		if(EQ(hy, hu)) {
			resultY = hy;
			resultU = hu;
			return;
		}
		
		// attempts to calculate ω, σ, and y
		// if it cannot find a y that works it throws an exception
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
			// calculuates μ_c
			char sigma = alphabet[c];
			
			// calculates the vectors F_xi and F_{xi+σ}
			double[][] F_xi = new double[l][l];
			double[][] F_xi_sigma = new double[l][l];
			for(int i=0;i<l;i++) {
				for(int j=0;j<l;j++) {
					F_xi[j][i] = MQ(X.get(i) + Y.get(j));
					F_xi_sigma[i][j] = MQ(X.get(i) + sigma + Y.get(j));
				}
			}
			
			// solves the matrix equation using LU Decomposition
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
	
	public static int MQ(String w) throws Exception {
		// MQ for the target function
		
		// MQ(ω) was previously calculated and is in the Hankel matrix
		if(F.get(w) != null)
			return F.get(w);
		
		int out = 0;
		
		// uses the membership query function for NBA's defined in NBA.java
		if(NBA.F!=null)
			out = NBA.MQ(w);
		
		// uses a membership query function defined in MQ.java
		else if(arbitrary.MQarbitrary!=null) {
			try {
				out = (int) arbitrary.MQarbitrary.invoke(null,w);
			} catch (Exception e) {
				throwException(null, "Invalid input: invalid membership query function");
			} 
		}
		
		// performs the membership query on the target mod-2-MA
		else {
			// initializes cur as the rxr identity matrix
			RealMatrix cur = MatrixUtils.createRealIdentityMatrix(minR);
			
			// multiplies cur by the corresponding μ for each letter in ω
			for(int i=0;i<w.length();i++) {
				cur = cur.multiply(MatrixUtils.createRealMatrix(minU[letterToIndex.get(w.charAt(i))]));
				// accounts for rounding issues with larger words
				if(i!=0 && i%25==0) {
					for(int j=0;j<minR;j++) {
						for(int k=0;k<minR;k++)
							cur.setEntry(j, k, mod2(cur.getEntry(j, k)));
					}
				}
			}
			
			// multiplies the final result with γ
			out = mod2(cur.getRowVector(0).dotProduct(new ArrayRealVector(minY)));
		}
		
		// adds MQ(ω) to the Hankel matrix
		F.put(w, out);
		
		return out;
	}
	
	public static int MQH(double[] hy, double[][][] hu, String w) {
		// MQ for the current hypothesis
		
		// initializes cur as the lxl identity matrix
		RealMatrix cur = MatrixUtils.createRealIdentityMatrix(hy.length);
		
		// multiplies cur by the corresponding μ for each letter in ω
		for(int i=0;i<w.length();i++) {
			cur = cur.multiply(MatrixUtils.createRealMatrix(hu[letterToIndex.get(w.charAt(i))]));
			// accounts for rounding issues with larger words
			if(i!=0 && i%25==0) {
				for(int j=0;j<hy.length;j++) {
					for(int k=0;k<hy.length;k++)
						cur.setEntry(j, k, mod2(cur.getEntry(j, k)));
				}
			}
		}
		
		// multiplies the final result with γ
		return mod2(cur.getRowVector(0).dotProduct(new ArrayRealVector(hy)));
	}
	
	public static boolean EQ(double[] hy, double[][][] hu) throws Exception {
		// uses a statistical equivalence query
		if(NBA.F!=null || arbitrary.MQarbitrary!=null)
			return arbitrary.EQapprox(hy, hu);
		
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
		double[][][] mu = new double[alphabet.length][minR+l][minR+l];
		for(int i=0;i<alphabet.length;i++) {
			for(int j=0;j<minR+l;j++) {
				for(int k=0;k<minR+l;k++) {
					// fu forms the upper left block of μ
					if(j<minR && k<minR)
						mu[i][j][k] = minU[i][j][k];
					// hu forms the lower right block of μ
					else if(j>=minR && k>=minR)
						mu[i][j][k] = hu[i][j-minR][k-minR];
					// everything else is 0
					else
						mu[i][j][k] = 0;
				}
			}
		}
		
		// γ for the combined MA
		double[] y = new double[minR+l];
		for(int i=0;i<minR+l;i++) {
			// γ has the form [fy hy]
			if(i<minR)
				y[i] = minY[i];
			else
				y[i] = hy[i-minR];
		}
		
		// To form the basis, we will follow algorithm 1 detailed in the paper by Thon and Jaeger.
		// Basis for the set span(μ(ω)γ : ω∈Σ*)
		double[][] B = new double[minR+l][minR+l];
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
				if(mod2(w[0]+w[minR]) == 1) {
					z = s;
					return false;
				}
				
				// extends B
				B[sizeB++] = w;
				WB.add(s);
				
				// adds {μ(σ)ω | σ∈Σ} to C
				for(int i=0;i<alphabet.length;i++) {
					RealMatrix m = MatrixUtils.createRealMatrix(mu[i]);
					RealVector p = MatrixUtils.createRealVector(w);
					double[] v = m.operate(p).toArray();
					for(int j=0;j<v.length;j++)
						v[j] = mod2(v[j]);
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
		
		// version that calculates linear independence using orthogonal projections
		/*
		if(sizeB==0)
			return true;

		// forms the matrix where the columns are the vectors in the basis B
		RealMatrix m = MatrixUtils.createRealMatrix(w.length, sizeB);
		for(int i=0;i<sizeB;i++)
			m.setColumn(i,B[i]);
		
		// calculates the orthogonal projection matrix P of the basis B
		RealMatrix P = m.multiply(MatrixUtils.inverse(m.transpose().multiply(m))).multiply(m.transpose());
		
		// if Pω = ω, ω is linearly dependent with B
		RealVector Pw = P.operate(MatrixUtils.createRealVector(w)); 
		for(int i=0;i<w.length;i++) {
			if(Math.round(Pw.getEntry(i)) != Math.round(w[i]))
				return true;
		}
		
		return false;
		*/
	}
	
	public static void calcWSigY(double[][][] hu) throws Exception {
		// prefix of z = ω + σ
		String w = "";
		char sigma = 0;
		// experiment
		String y = "";
		
		// goes through every possible prefix of z starting with ω = "" and σ = (first character of ω)
		for(int i=0;i<z.length();i++) {
			if(i!=0)
				w = z.substring(0,i);
			sigma = z.charAt(i);
			
			// calculates μ(ω)
			RealMatrix mu = MatrixUtils.createRealIdentityMatrix(l);
			for(int n=0;n<w.length();n++)
				mu = mu.multiply(MatrixUtils.createRealMatrix(hu[letterToIndex.get(w.charAt(n))]));
			
			// checks if F_ω = sum(μ(ω)_1,i * F_xi)
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
			
			// goes through every possible value of y in Y
			// checks if F_{ω+σ}(y) != sum(μ(ω)_1,i * F_{xi+σ}(y))
			for(int j=0;j<l;j++) {
				y = Y.get(j);
			
				int sum = 0;
				for(int k=0;k<l;k++)
					sum = mod2(sum + mu.getEntry(0, k)*MQ(X.get(k) + sigma + y));
				
				// found a solution
				if(MQ(w+sigma+y) != sum) {		
					if(l==minR)
						throwException(null,"Algorithm failed: size of the hypothesis exceeds that of the target function.");
					// updates l, X, and Y
					l++;
					X.add(w);
					Y.add(sigma+y);
					
					// displays the updated observation table
					if(verbose)
						displayQueries();
					
					return;
				}
			}
		}

		throwException(null,"Algorithm failed: didn't find a suitable omega, sigma, and gamma.");
		
		// version that traverses through all prefixes of z before calling a new equivalence query
		/*
		// prefix of z = ω + σ
		String w = "";
		char sig = 0;
		// experiment
		String y = "";
		boolean noSoln = true;
		
		// goes through every possible prefix of z starting with ω = "" and σ = (first character of ω)
		for(int i=0;i<z.length();i++) {
			if(i!=0)
				w = z.substring(0,i);
			sig = z.charAt(i);
			
			// calculates μ(ω)
			RealMatrix u = MatrixUtils.createRealIdentityMatrix(l);
			for(int n=0;n<w.length();n++)
				u = u.multiply(MatrixUtils.createRealMatrix(hu[letterToIndex.get(w.charAt(n))]));
			
			// checks if F_ω = sum(μ(ω)_1,i * F_xi)
			boolean failed = false;
			for(int j=0;j<l;j++) {
				int sum = 0;
				for(int k=0;k<l;k++)
					sum = mod2(sum + u.getEntry(0, k)*MQ(X.get(k)+Y.get(j)));
				
				if(MQ(w+Y.get(j)) != sum) {
					failed = true;
					break;
				}
			}
			if(failed)
				continue;
			
			// goes through every possible value of y in Y
			// checks if F_{ω+σ}(y) != sum(μ(ω)_1,i * F_{xi+σ}(y))
			for(int j=0;j<l;j++) {
				y = Y.get(j);
			
				int sum = 0;
				for(int k=0;k<l;k++)
					sum = mod2(sum + u.getEntry(0, k)*MQ(X.get(k) + sig + y));
				
				// found a solution
				if(MQ(w+sig+y) != sum) {		
					if(l==r)
						throwException(null,"Algorithm failed: size of the hypothesis exceeds that of the target function.");
					// updates l, X, and Y
					l++;
					X.add(w);
					Y.add(sig+y);
					
					// displays the updated observation table
					if(verbose)
						displayQueries();
					
					noSoln = false;
					hu = createHU();
					break;
				}
			}
		}
		
		if(noSoln)
			throwException(null,"Algorithm failed: didn't find a suitable omega, sigma, and gamma.");
		*/
	}

	public static void displayResults() {
		System.out.println("Learned mod-2-MA");
		System.out.println("----------------");
		// prints γ
		System.out.print("y: ");
		String s = "";
		for(int i=0;i<resultY.length;i++)
			s += mod2(resultY[i]) + " ";
		System.out.println(s + "\n");
		
		// prints the μ's
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
	
	public static void displayQueries() throws Exception {
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
			test += alphabet[(int)(Math.random()*alphabet.length)];
		return test;
	}
	
	public static boolean finalCheck(int maxTestLen, int numTests) throws Exception {
		// creates numTests tests of length at most maxTestLen
		// checks whether the learned function and target function have the same output
		for(int i=1;i<=numTests;i++) {
			String test = genTest((int)(Math.random()*(maxTestLen+1)));
			if(MQH(inputY,inputU,test)!=MQH(resultY, resultU, test))
				return false;
		}
		return true;
	}

	public static void operations() {
		System.out.println("Available operations for the learned Mod-2-MA (enter \"quit\" to terminate):");
		System.out.println("- Test whether a word is accepted: \"-a <word>\"\n  If the language is (L)_$, words must be of the form u$v.");
		while(true) {
			// reads in cmd
			String line = in.nextLine();
			
			// terminates the program
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