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
 * 3 Dana Angluin, Timos Antonopoulos, Dana Fisman. Strongly Unambiguous Büchi Automata Are Polynomially 
 *   Predictable with Membership Queries. 28th International Conference on Computer Science Logic, 8:1–8:17, 2020.
 * 4 Michael Thon and Herbert Jaeger. Links Between Multiplicity Automata, Observable Operator Models and 
 *   Predictive State Representations — a Unified Learning Framework. Journal of Machine Learning Research, 
 *   16(4):103−147, 2015.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.StringTokenizer;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class M2MA {
	
	/*
	 * Matrices are represented as sparse matrices using a HashMap<Integer, ArrayList<Integer>>.
	 * Each element in the hashmap consists of a (Integer, ArrayList<Integer>) pair.
	 * 
	 * map.get(0) returns an arraylist of size 2, where the first element is the row dimension, and
	 * the second element is the column dimension.
	 * I.e. map.get(0).get(0) = row dimension, map.get(0).get(1) = column dimension
	 * 
	 * Every row/column in the matrix is represented as an arraylist, where the arraylist contains
	 * the indices (starting from 1 instead of 0) corresponding to locations of the 1's.
	 * E.g. [1, 0, 0, 1, 0, 0, 0, 1] -> [1, 4, 8]
	 * 
	 * For any positive integer r, map.get(r) returns the arraylist corresponding to row r.
	 * For any negative integer c, map.get(c) returns the arraylist corresponding to row -c.
	 * 
	 * Useful sparse matrix functions are implemented such as initializing matrices, adding elements, 
	 * and multiplying matrices.
	 */
	
	// if true, displays the observation table as it is constructed
	public static boolean observationTableFlag;
	// if true, displays information on the progress of the minimization algorithm
	public static boolean minProgressFlag;
	// if true, only displays the minimized dimension and terminates
	public static boolean minDimensionFlag;
	// if true, displays the dimension of an equivalent minimal DFA
	public static boolean dfaFlag;
	
	public static String[] alphabet;
	// maps each letter in the alphabet to an index
	public static HashMap<String, Integer> letterToIndex;
	
	// input mod-2-MA
	public static int inputSize;
	public static HashMap<Integer, ArrayList<Integer>> inputFinalVector;
	public static HashMap<Integer, ArrayList<Integer>>[] inputTransitionMatrices;
	public static HashMap<String, Integer> Hankel;
	
	// minimized mod-2-MA
	public static HashMap<Integer, ArrayList<Integer>> minFinalVector;
	public static HashMap<Integer, ArrayList<Integer>>[] minTransitionMatrices;
	public static int minSize;
	// row and column indices of the minimized mod-2-MA's observation table
	public static ArrayList<String> minRowIndices;
	public static ArrayList<String> minColIndices;
	
	// mod-2-MA being learned
	public static int learnedSize;
	// row and column indices of the observation table being constructed
	public static ArrayList<String> learnedRowIndices;
	public static ArrayList<String> learnedColIndices;
	public static String counterExample;

	// learned mod-2-MA
	public static HashMap<Integer, ArrayList<Integer>> resultFinalVector;
	public static HashMap<Integer, ArrayList<Integer>>[] resultTransitionMatrices;
	
	// used in EQ to avoid testing the same word
	public static boolean[][] tested;
	public static boolean[][][][] testedExtension;
	
	public static Scanner in;
	public static long startTime;
	
	// true if running minimize.java
	public static boolean inMinimize = false;
	
	public static void main(String[] args) throws Exception {
		readInput();
		
		minimize();
				
		learn();
		
		if (minSize != learnedSize) {
			throwException(null, "Algorithm failed: the learned mod-2-MA has a different dimension "
					+ "(" + learnedSize + ") than the minimized mod-2-MA (" + minSize + ").");
		}

		if (finalCheck(25, 1000, false)) {
			displayResults();
		} else {
			throwException(null, "Algorithm failed: failed final check.");
		}
		
		displayRuntime();
		
		operationsOnLearnedMA();
	}
	
	@SuppressWarnings("unchecked")
	public static void readInput() throws Exception {	
		BufferedReader f;
		if (inMinimize) {
			System.out.println("Input file name and optional flags -ma (e.g. M2MA_input1.txt, M2MA_input1.txt -m)");
			f = getFile(false, true, false, true);
		} else {
			System.out.println("Input file name and optional flags -vma (e.g. M2MA_input1.txt -v, M2MA_input1.txt -ma, M2MA_input1.txt -vma)");
			f = getFile(true, true, false, true);
		}
		
		readAlphabet(f, false);
		
		inputSize = Integer.parseInt(readFile(f));
		
		StringTokenizer st = new StringTokenizer(readFile(f));
		inputFinalVector = initialize(1, inputSize);
		for (int i=1; i<=inputSize; i++) {
			if (Integer.parseInt(st.nextToken()) == 1) {
				addElement(inputFinalVector, 1, i);
			}
		}
		if (st.hasMoreTokens()) {
			throwException(f, "Invalid input: final vector length exceeds the specified size.");
		}
		
		inputTransitionMatrices = new HashMap[alphabet.length];
		for (int i=0; i<alphabet.length; i++) {
			HashMap<Integer, ArrayList<Integer>> transitionMatrix = initialize(inputSize, inputSize);
			
			for (int j=1; j<=inputSize; j++) {
				st = new StringTokenizer(readFile(f));
				for (int k=1; k<=inputSize; k++) {
					if (Integer.parseInt(st.nextToken()) == 1) {
						addElement(transitionMatrix, j, k);
					}
				}
				if (st.hasMoreTokens()) {
					throwException(f, "Invalid input: transition matrix size exceeds the specified size.");
				}
			}
			
			inputTransitionMatrices[i] = transitionMatrix;
		}
		if (readFile(f) != null) {
			throwException(f, "Invalid input: transition matrix size exceeds the specified size.");
		}
		
		f.close();
	}
	
	public static BufferedReader getFile(boolean vFlag, boolean mFlag, boolean dFlag, boolean aFlag) throws Exception {
		in = new Scanner(System.in);
		String[] arrInput = in.nextLine().split(" ");
		startTime = System.nanoTime();
		
		if (arrInput.length > 2) {
			throwException(null, "Invalid input: too many inputs passed.");
		}
		
		// handle flags
		observationTableFlag = false;
		minProgressFlag = false;
		minDimensionFlag = false;
		dfaFlag = false;
		if (arrInput.length == 2) {
			if (vFlag && arrInput[1].contains("v")) {
				observationTableFlag = true;
			}
			if (mFlag && arrInput[1].contains("m")) {
				minProgressFlag = true;
			}
			if (dFlag && arrInput[1].contains("d")) {
				minDimensionFlag = true;
			}
			if (aFlag && arrInput[1].contains("a")) {
				dfaFlag = true;
			}
		}
		
		if (minProgressFlag && minDimensionFlag) {
			throwException(null, "Invalid input: cannot have both the -m and -d flags.");
		}
		
		System.out.println();
		
		return new BufferedReader(new FileReader(arrInput[0]));
	}
	
	public static String readFile(BufferedReader f) throws IOException {
		String line = f.readLine();
		
		// ignore empty lines and lines beginning with "//"
		while (line != null && (line.equals("") || (line.charAt(0) == '/' && line.charAt(1) == '/' ))) {
			line = f.readLine();
		}
		return line;
	}
	
	public static void readAlphabet(BufferedReader f, boolean inLDollar) throws Exception {
		ArrayList<String> tempAlphabet = new ArrayList<String>();
		
		StringTokenizer st = new StringTokenizer(readFile(f));
		while (st.hasMoreTokens()) {
			String letter = st.nextToken();
			if(inLDollar && letter.equals("$")) {
				throwException(f, "Invalid input: invalid character in the alphabet.");
			}
			tempAlphabet.add(letter);
		}
		
		if (inLDollar) {
			tempAlphabet.add("$");
		}
		
		alphabet = new String[tempAlphabet.size()];
		for (int i=0; i<tempAlphabet.size(); i++) {
			alphabet[i] = tempAlphabet.get(i);
		}
		
		// map each letter in the alphabet to an index
		letterToIndex = new HashMap<String, Integer>();
		for (int i=0; i<alphabet.length; i++) {
			letterToIndex.put(alphabet[i], i);
		}
	}
	
	// properly closes input streams
	public static void throwException(BufferedReader f, String message) throws Exception {
		if (f != null) {
			f.close();
		}
		if (in != null) {
			in.close();
		}
		throw new Exception(message);
	}
	
	@SuppressWarnings("unchecked")
	public static void minimize() throws Exception {
		int inConvert = 0;
		if (convert.results != null) {
			if (convert.inM2MA) {
				inConvert = 1;
			} else {
				inConvert = 2;
			}
		}
		
		if ((inConvert == 0) && minProgressFlag) {
			System.out.println("Mod-2-MA to minimize:");
			System.out.println("---------------------");
			
			System.out.println("Dimension: " + inputSize + '\n');
			
			System.out.print("Final Vector: ");
			displayMatrix(inputFinalVector);
			
			System.out.println("Transition Matrices:\n");
			for (int i=0; i<inputTransitionMatrices.length; i++) {
				System.out.println("Letter " + alphabet[i]);
				displayMatrix(inputTransitionMatrices[i]);
			}
			
			System.out.println("Minimization in progress...");
			System.out.println("-------------------------");
		} else if ((inConvert == 0) && minDimensionFlag) {
			System.out.println("Minimization in progress...");
		}
		
		ArrayList<String> stateSpaceBasisIndices = new ArrayList<String>();
		HashMap<String, HashMap<Integer, ArrayList<Integer>>> stateSpaceIndexToVector = new HashMap<String, HashMap<Integer, ArrayList<Integer>>>();
		HashMap<Integer, ArrayList<Integer>> stateSpaceBasis = basis(inputFinalVector, inputTransitionMatrices, stateSpaceIndexToVector, stateSpaceBasisIndices, true);
		
		if ((inConvert == 0) && (minProgressFlag || minDimensionFlag)) {
			System.out.println("Created the state space.");
		}
		
		ArrayList<String> coStateSpaceBasisIndices = new ArrayList<String>();
		HashMap<String, HashMap<Integer, ArrayList<Integer>>> coStateSpaceIndexToVector = new HashMap<String, HashMap<Integer, ArrayList<Integer>>>();
		HashMap<Integer, ArrayList<Integer>> coStateSpaceBasis = basis(inputFinalVector, inputTransitionMatrices, coStateSpaceIndexToVector, coStateSpaceBasisIndices, false);
		
		if ((inConvert == 0) && (minProgressFlag || minDimensionFlag)) {
			System.out.println("Created the co-state space.");
			if (minProgressFlag) {
				System.out.println();
			}
		}
		
		// (state space x co-state space) observation table
		HashMap<Integer, ArrayList<Integer>> observationTable = multiply(stateSpaceBasis, coStateSpaceBasis);
		
		if ((inConvert == 0) && minDimensionFlag) {
			System.out.println("Created the observation table.");
		}
		
		if ((inConvert == 0) && minProgressFlag) {
			System.out.println("Observation table:" );
			System.out.println("Dimension: " + stateSpaceBasis.get(0).get(0) + " x " + coStateSpaceBasis.get(0).get(1));
			System.out.println("Rows: " + displayIndices(stateSpaceBasisIndices));
			System.out.println("Cols: " + displayIndices(coStateSpaceBasisIndices));
			displayMatrix(observationTable);
		}
		
		// obtain the smallest set of linearly independent rows and columns from observationTable
		minRowIndices = new ArrayList<String>();
		HashMap<Integer, ArrayList<Integer>> linIndRowsObservationTable = linIndSubMatrixRows(observationTable, stateSpaceBasisIndices, minRowIndices);
		
		minSize = linIndRowsObservationTable.get(0).get(0);
		
		if (inConvert == 1) {
			return;
		}
		
		if ((inConvert == 0) && minDimensionFlag) {
			System.out.println("Minimized dimension: " + minSize);
			
			if (!dfaFlag) {
				if (in != null) {
					in.close();
				}
				System.exit(0);
			}
		}
		
		minColIndices = new ArrayList<String>();
		HashMap<Integer, ArrayList<Integer>> minObservationTable = linIndSubMatrixCols(linIndRowsObservationTable, coStateSpaceBasisIndices, minColIndices);
		
		if ((inConvert == 0) && minProgressFlag) {
			System.out.println("Minimized observation table:");
			System.out.println("Dimension: " + minSize);
			System.out.println("Rows: " + displayIndices(minRowIndices));
			System.out.println("Cols: " + displayIndices(minColIndices));
			System.out.println("Table:");
			displayMatrix(minObservationTable);
		}
		
		// case where minObservationTable = [[0]] (singular, must be treated separately)
		if ((int) minObservationTable.get(0).get(0) == 1 && minObservationTable.get(1) == null) {
			minFinalVector = initialize(1, 1);
			minTransitionMatrices = new HashMap[alphabet.length];		
			for (int i=0; i<alphabet.length; i++) {
				minTransitionMatrices[i] = initialize(1, 1);
			}
			
			if (inConvert == 0) {
				if (minProgressFlag) {
					System.out.println("Minimization completed.\n");
				}
				
				tested = new boolean[1][1];	
			}
			
			return;
		}
		
		DecompositionSolver solver = new solver(sparseToReal(minObservationTable)).getSolver();
		HashMap<Integer, ArrayList<Integer>> tableInverse = realToSparse(solver.getInverse());
		
		Hankel = new HashMap<String, Integer>();
		
		// minTransitionMatrices = xSigma*tableInverse, where xSigma is the matrix where row_i = row_(x_i+σ) of the observation table
		minTransitionMatrices = new HashMap[alphabet.length];
		int dim = minObservationTable.get(0).get(0);
		for (int i=0; i<alphabet.length; i++) {	
			minTransitionMatrices[i] = initialize(dim, dim);
			HashMap<Integer, ArrayList<Integer>> xSigma = initialize(dim, dim);
			
			for (int j=0; j<dim; j++) {			
				for (int k=0; k<dim; k++) {
					HashMap<Integer, ArrayList<Integer>> stateVector = stateSpaceIndexToVector.get(minRowIndices.get(j));
					HashMap<Integer, ArrayList<Integer>> coStateVector = coStateSpaceIndexToVector.get(minColIndices.get(k));
					
					if (dotProduct(multiply(stateVector, inputTransitionMatrices[i]).get(1), coStateVector.get(-1)) == 1) {
						addElement(xSigma, j+1, k+1);
					}
				}
			}

			minTransitionMatrices[i] = multiply(xSigma, tableInverse);
		}
		
		// minFinalVector is the first column of minObservationTable
		minFinalVector = initialize(1, minObservationTable.get(0).get(0));
		if (minObservationTable.get(-1) != null) {
			for (int num : minObservationTable.get(-1)) {
				addElement(minFinalVector, 1, num);
			}
		}
		
		if (inConvert == 0) {
			// used in EQ to avoid testing the same word
			tested = new boolean[minRowIndices.size()][minColIndices.size()];
			
			if (minProgressFlag) {
				System.out.println("Minimization completed.\n");
			}
			
			if (dfaFlag) {
				System.out.println("Number of states of equivalent minimal DFA: " + dimensionMinDFA(true) + "\n");
				
				if (minDimensionFlag) {
					if (in != null) {
						in.close();
					}
					System.exit(0);
				}
				
			}
		}
	}
	
	public static String displayIndices(ArrayList<String> indices) {
		String out = "";
		for (String index : indices) {
			if (index.length() == 0) {
				out += "ɛ ";
			} else {
				out += removeSpaces(index) + " ";
			}
		}	
		return out;
	}
	
	// follows algorithm 1 detailed in Thon and Jaeger to form the basis for the state/co-state space
	public static HashMap<Integer, ArrayList<Integer>> basis(HashMap<Integer, ArrayList<Integer>> hypothesisFinalVector, HashMap<Integer, ArrayList<Integer>>[] hypothesisTransitionMatrices, HashMap<String, HashMap<Integer, ArrayList<Integer>>> indexToVector, ArrayList<String> indices, boolean stateSpace) throws Exception {
		// set with elements to try to add to the basis
		ArrayList<HashMap<Integer, ArrayList<Integer>>> tests = new ArrayList<HashMap<Integer, ArrayList<Integer>>>();
		
		HashMap<Integer, ArrayList<Integer>> basis;
		int sizeBasis = 0;
		HashMap<Integer, ArrayList<Integer>> hypothesisFinalVectorTranpose = null;
		
		if (stateSpace) {
			basis = initialize(1, hypothesisFinalVector.get(0).get(1));
			
			// begin with ω_i = (1,0,0,...,0)
			HashMap<Integer, ArrayList<Integer>> w_i = initialize(1, hypothesisFinalVector.get(0).get(1));
			addElement(w_i, 1, 1);
			tests.add(w_i);
		} else {
			basis = initialize(hypothesisFinalVector.get(0).get(1), 1);		
			
			// begin with the transpose of hypothesisFinalVector
			hypothesisFinalVectorTranpose = initialize(hypothesisFinalVector.get(0).get(1), 1);
			if (hypothesisFinalVector.get(1) != null) {
				for (int col : hypothesisFinalVector.get(1)) {
					addElement(hypothesisFinalVectorTranpose, col, 1);
				}
			}
			
			tests.add(hypothesisFinalVectorTranpose);
		}
		
		int sizeTests = 1;
		
		// contains the corresponding string for every element in tests
		ArrayList<String> testStrings = new ArrayList<String>();
		testStrings.add("");
		
		ArrayList<int[]> operations = new ArrayList<int[]>();
		
		while (sizeTests > 0) {
			HashMap<Integer, ArrayList<Integer>> test = tests.remove(0);
			String testString = testStrings.remove(0);
			sizeTests--;
			
			if (linInd(stateSpace ? test.get(1) : test.get(-1), operations, sizeBasis, hypothesisFinalVector.get(0).get(1))) {	
				// extend the basis
				sizeBasis++;
				if (stateSpace) {
					basis.get(0).set(0, sizeBasis);
					for (int num : test.get(1)) {
						addElement(basis, sizeBasis, num);
					}
				} else {
					basis.get(0).set(1, sizeBasis);
					for (int num : test.get(-1)) {
						addElement(basis, num, sizeBasis);
					}
				}
				
				indices.add(testString);
				indexToVector.put(testString, test);
				
				// add to tests the one-letter extensions of test
				for (int i=0; i<alphabet.length; i++) {
					HashMap<Integer, ArrayList<Integer>> newTest;
					
					if (stateSpace) {
						// basis for the set span((initial vector) * (transitionMatrix_ω) : ω∈Σ*)
						newTest = multiply(test, hypothesisTransitionMatrices[i]);
						testStrings.add(addStrings(testString, alphabet[i]));
					} else {
						// basis for the set span((transitionMatrix_ω) * (final vector) : ω∈Σ*)
						newTest = multiply(hypothesisTransitionMatrices[i], test);
						testStrings.add(addStrings(alphabet[i], testString));
					}
					
					tests.add(newTest);
					sizeTests++;
				}
			}
		}
		
		if (sizeBasis == 0) {
			indices.add("");
			indexToVector.put("", hypothesisFinalVectorTranpose);
		}
		
		return basis;
	}
	
	// finds a maximal submatrix of linearly independent rows/columns of the observation table
	public static HashMap<Integer, ArrayList<Integer>> linIndSubMatrixRows(HashMap<Integer, ArrayList<Integer>> observationTable, ArrayList<String> oldIndices, ArrayList<String> newIndices) throws Exception {
		HashMap<Integer, ArrayList<Integer>> newObservationTable = initialize(0, observationTable.get(0).get(1));
		int sizeT = 0;
		
		ArrayList<int[]> operations = new ArrayList<int[]>();
		
		for (int row : observationTable.keySet()) {
			if (row > 0) {
				// extend the current subset of linearly independent rows/columns
				if (linInd(observationTable.get(row), operations, sizeT, newObservationTable.get(0).get(1))) {
					sizeT++;
					newObservationTable.get(0).set(0, sizeT);
					
					for (int num : observationTable.get(row)) {
						addElement(newObservationTable, sizeT, num);
					}
					
					newIndices.add(oldIndices.get(row - 1));
				}
			}
		}
		
		if ((int) newObservationTable.get(0).get(0) == 0 || (int) newObservationTable.get(0).get(1) == 0) {
			newObservationTable.get(0).set(0, 1);
			newObservationTable.get(0).set(1, 1);
			newIndices.add("");
		}
		
		return newObservationTable;
	}
	
	// finds a maximal submatrix of linearly independent rows/columns of the observation table
	public static HashMap<Integer, ArrayList<Integer>> linIndSubMatrixCols(HashMap<Integer, ArrayList<Integer>> observationTable, ArrayList<String> oldIndices, ArrayList<String> newIndices) throws Exception {
		HashMap<Integer, ArrayList<Integer>> newObservationTable = initialize(observationTable.get(0).get(0), 0);
		int sizeT = 0;
		
		ArrayList<int[]> operations = new ArrayList<int[]>();
		
		for (int col : observationTable.keySet()) {
			if (col < 0) {
				// extend the current subset of linearly independent rows/columns
				if (linInd(observationTable.get(col), operations, sizeT, observationTable.get(0).get(0))) {
					sizeT++;
					newObservationTable.get(0).set(1, sizeT);
					
					for (int num : observationTable.get(col)) {
						addElement(newObservationTable, num, sizeT);
					}
					
					newIndices.add(oldIndices.get((col * -1) - 1));
				}
			}
		}
		
		if ((int) newObservationTable.get(0).get(0) == 0 || (int) newObservationTable.get(0).get(1) == 0) {
			newObservationTable.get(0).set(0, 1);
			newObservationTable.get(0).set(1, 1);
			newIndices.add("");
		}
		
		return newObservationTable;
	}
	
	/*
	 * Returns false if vector is in the span of the basis, true otherwise.
	 * 
	 * Operations consists of all the elementary row operations needed to convert the basis into rref.
	 * These operations are applied to vector, and the additional elementary row operations needed to
	 * convert the augmented matrix basis+vector into rref are added to operations.
	 * 
	 * Each int[] in operations has size 3, where int[0] = 0 for a swap operation and
	 * 1 for a subtract operation, and int[1] = row_1 and int[2] = row_2, meaning either
	 * you swap row_1 and row_2 or subtract row_1 from row_2. 
	 */
	public static boolean linInd(ArrayList<Integer> vector, ArrayList<int[]> operations, int sizeBasis, int numRows) {	
		if (vector == null || vector.size() == 0) {
			return false;
		}
		
		boolean[] booleanVector = new boolean[numRows + 1];
		for (int n : vector) {
			booleanVector[n] = true;
		}
		
		// apply all of the previously done elementary row operations to vector
		for (int[] operation : operations) {
			if ((int) operation[0] == 0) {
				// swap operation
				if (booleanVector[operation[1]] && !booleanVector[operation[2]]) {
					booleanVector[operation[1]] = false;
					booleanVector[operation[2]] = true;
				} else if (!booleanVector[operation[1]] && booleanVector[operation[2]]) {
					booleanVector[operation[1]] = true;
					booleanVector[operation[2]] = false;
				}
			} else {
				// subtract operation 
				if (booleanVector[operation[1]]) {
					booleanVector[operation[2]] = !booleanVector[operation[2]];
				}
			}
		}
		
		ArrayList<Integer> rows = new ArrayList<Integer>();
		int index = -1;
		
		// convert booleanVector to a sparse vector
		// find the first 1 in a row after sizeBasis
		for (int row = 1; row <= numRows; row++) {
			if (booleanVector[row]) {
				if (row > sizeBasis && index == -1) {
					index = row;
				}
				rows.add(row);
			}
		}
		
		// no 1 was found, so vector is linearly dependent with the basis
		if (index == -1) {
			return false;
		}
		
		// add the operations needed to convert the augmented matrix basis+vector into rref
		
		// swap operations
		if (index != sizeBasis + 1) {
			// need to perform a swap
			int[] swap = new int[3];
			swap[0] = 0;
			swap[1] = index;
			swap[2] = sizeBasis + 1;
			operations.add(swap);
		}
		
		// subtract operations
		for (int row : rows) {
			if (row != index) {
				// need to subtract two rows
				int[] subtract = new int[3];
				subtract[0] = 1;
				subtract[1] = sizeBasis + 1;
				subtract[2] = row;
				operations.add(subtract);
			}
		}
		
		return true;
	}
	
	public static void learn() throws Exception {	
		learnedRowIndices = new ArrayList<String>();
		learnedColIndices = new ArrayList<String>();
		learnedRowIndices.add("");
		learnedColIndices.add("");
		learnedSize = 1;
		
		if (Hankel == null) {
			Hankel = new HashMap<String, Integer>();
		}
		
		/* 
		 * F("") cannot equal 0 (otherwise can't form a linearly independent basis of elements in learnedRowIndices).
		 * The algorithm instead begins with a 2x2 matrix of full rank.
		 */
		if (MQ("") == 0) {
			HashMap<Integer, ArrayList<Integer>> hypothesisFinalVector = createHypothesisFinalVector();
			HashMap<Integer, ArrayList<Integer>>[] hypothesisTransitionMatrices = createHypothesisTransitionMatrices();
			
			if (!EQ(hypothesisFinalVector, hypothesisTransitionMatrices)) {
				learnedSize++;
				learnedRowIndices.add(counterExample);
				learnedColIndices.add(counterExample);
			}
		}
		
		if (observationTableFlag) {
			System.out.println("Observation table after individual queries\n------------------------------------------");
			displayTable();
		}

		learnMain();
	}
	
	public static void learnMain() throws Exception {
		HashMap<Integer, ArrayList<Integer>> hypothesisFinalVector = createHypothesisFinalVector();
		HashMap<Integer, ArrayList<Integer>>[] hypothesisTransitionMatrices = createHypothesisTransitionMatrices();
		
		if (EQ(hypothesisFinalVector, hypothesisTransitionMatrices)) {
			resultFinalVector = hypothesisFinalVector;
			resultTransitionMatrices = hypothesisTransitionMatrices;
			return;
		}
		
		growObservationTable(hypothesisTransitionMatrices);
		
		learnMain();
	}
	
	public static HashMap<Integer, ArrayList<Integer>> createHypothesisFinalVector() throws Exception {
		HashMap<Integer, ArrayList<Integer>> hypothesisFinalVector = initialize(1, learnedSize);
		for (int i=0; i<learnedSize; i++) {
			if (MQ(learnedRowIndices.get(i)) == 1) {
				addElement(hypothesisFinalVector, 1, i+1);
			}
		}
		return hypothesisFinalVector;
	}
	
	@SuppressWarnings("unchecked")
	public static HashMap<Integer, ArrayList<Integer>>[] createHypothesisTransitionMatrices() throws Exception {
		/*
		 * For every letter in alphabet, define a transition matrix by letting its i-th row be the coefficients 
		 * of the vector F_{xi+letter}(y) when expressed as a linear combination of the row vectors of F (such 
		 * coefficients exist as the row vectors are linearly independent).
		 */
		HashMap<Integer, ArrayList<Integer>>[] hypothesisTransitionMatrices = new HashMap[alphabet.length];
		for (int c=0; c<alphabet.length; c++) {
			hypothesisTransitionMatrices[c] = initialize(learnedSize, learnedSize);
			String letter = alphabet[c];
			
			double[][] F_xi = new double[learnedSize][learnedSize];
			double[][] F_xi_letter = new double[learnedSize][learnedSize];
			for (int i=0; i<learnedSize; i++) {
				for (int j=0; j<learnedSize; j++) {
					F_xi[j][i] = MQ(addStrings(learnedRowIndices.get(i), learnedColIndices.get(j)));
					F_xi_letter[i][j] = MQ(addStrings(addStrings(learnedRowIndices.get(i), letter), learnedColIndices.get(j)));
				}
			}
			
			// solve the matrix equation using LU Decomposition
			RealMatrix coefficients = new Array2DRowRealMatrix(F_xi);
			DecompositionSolver solver = new solver(coefficients).getSolver();
			for (int i=0; i<learnedSize; i++) {
				RealVector constants = new ArrayRealVector(F_xi_letter[i]);
				try {
					RealVector solution = solver.solve(constants);
					for(int j=0; j<learnedSize; j++) {
						if (mod2(solution.getEntry(j)) == 1) {
							addElement(hypothesisTransitionMatrices[c], i + 1, j + 1);
						}
					}
				} catch(Exception e) {
					// matrix is not invertible
					hypothesisTransitionMatrices[c] = initialize(learnedSize, learnedSize);
				}
			}
		}
		return hypothesisTransitionMatrices;
	}
	
	// MQ for the target function
	public static int MQ(String word) throws Exception {	
		// MQ(ω) was previously calculated and is in the Hankel matrix
		if (Hankel.get(word) != null) {
			return Hankel.get(word);
		}
		
		int out = 0;
		
		// NBA.java and arbitrary.java use their own MQ functions
		if(NBA.NBAFinalStates != null) {
			out = NBA.MQ(word);
		} else if(arbitrary.MQMethod != null) {
			try {
				out = (int) arbitrary.MQMethod.invoke(null, word);
			} catch (Exception e) {
				throwException(null, "Invalid input: invalid membership query function.");
			} 
		} else {
			HashMap<Integer, ArrayList<Integer>> current = identity(minSize);
			
			String[] wordArr = word.split(" ");
			if (word.length() == 0) {
				wordArr = new String[0];
			}
			
			for (int i=0; i<wordArr.length; i++) {
				current = multiply(current, minTransitionMatrices[letterToIndex.get(wordArr[i])]);
			}
			
			out = dotProduct(current.get(1), minFinalVector.get(1));
		}
		
		Hankel.put(word, out);
		
		return out;
	}
	
	// MQ for any given final vector and set of transition matrices
	public static int MQArbitrary(HashMap<Integer, ArrayList<Integer>> finalVector, HashMap<Integer, ArrayList<Integer>>[] transitionMatrices, String word) throws Exception {	
		HashMap<Integer, ArrayList<Integer>> current = identity(finalVector.get(0).get(1));
		
		String[] wordArr = word.split(" ");
		if (word.length() == 0) {
			wordArr = new String[0];
		}
		
		for (int i=0; i<wordArr.length; i++) {
			current = multiply(current, transitionMatrices[letterToIndex.get(wordArr[i])]);
		}
		
		return dotProduct(current.get(1), finalVector.get(1));
	}
	
	public static boolean EQ(HashMap<Integer, ArrayList<Integer>> hypothesisFinalVector, HashMap<Integer, ArrayList<Integer>>[] hypothesisTransitionMatrices) throws Exception {
		// NBA.java and arbitrary.java use statistical EQ's
		if (NBA.NBAFinalStates != null || arbitrary.MQMethod != null) {
			return arbitrary.EQstatistical(hypothesisFinalVector, hypothesisTransitionMatrices);
		}
		
		// test every element in the observation table of the minimized mod-2-MA
		for (int i=0; i<minRowIndices.size(); i++) {
			for (int j=0; j<minColIndices.size(); j++) {
				if (!tested[i][j]) {
					// update tested to avoid testing the same words in the next EQ
					tested[i][j] = true;
					
					String test = addStrings(minRowIndices.get(i), minColIndices.get(j));
	
					if (MQ(test) != MQArbitrary(hypothesisFinalVector, hypothesisTransitionMatrices, test)) {
						counterExample = test;
						return false;
					}
				}
			}
		}
		
		if (learnedSize == minSize) {
			return true;
		}
		
		System.out.println("Testing one-letter extensions.");
		
		testedExtension = new boolean[minRowIndices.size()][minColIndices.size()][alphabet.length+1][alphabet.length+1];
		
		// test the one-letter extensions of the row and column indices of the observation table
		for (int i=0; i<minRowIndices.size(); i++) {
			for (int j=0; j<minColIndices.size(); j++) {
				for (int a1=0; a1<=alphabet.length; a1++) {
					String letter1;
					if (a1 == alphabet.length) {
						letter1 = ""; 
					} else {
						letter1 = alphabet[a1];
					}
					
					for (int a2=0; a2<=alphabet.length; a2++) {
						if (!testedExtension[i][j][a1][a2] && (a1 != alphabet.length || a2 != alphabet.length)) {
							testedExtension[i][j][a1][a2] = true;
							
							String letter2;
							if (a2 == alphabet.length) {
								letter2 = ""; 
							} else {
								letter2 = alphabet[a2];
							}
							
							String test = addStrings(addStrings(addStrings(minRowIndices.get(i), letter1), minColIndices.get(j)), letter2);
							if (MQ(test) != MQArbitrary(hypothesisFinalVector, hypothesisTransitionMatrices, test)) {
								counterExample = test;
								return false;
							}
						}
					}
				}
			}
		}
		
		return true;
	}
	
	public static void growObservationTable(HashMap<Integer, ArrayList<Integer>>[] hypothesisTransitionMatrices) throws Exception {
		// prefix of the counter-example = ω + σ
		String w = "";
		String sigma = "";
		// experiment
		String y = "";
		
		String[] counterExampleArr = counterExample.split(" ");
		if (counterExample.length() == 0) {
			counterExampleArr = new String[0];
		}
		
		// go through every possible prefix of the counter-example starting with ω = "" and σ = (first character of ω)
		for (int i=0; i<counterExampleArr.length; i++) {
			if (i != 0) {
				w = "";
				for (int j=0; j<i-1; j++) {
					w += counterExampleArr[j] + " ";
				}
				if (i >= 1) {
					w += counterExampleArr[i - 1];
				}
			}
			sigma = counterExampleArr[i];
			
			HashMap<Integer, ArrayList<Integer>> transitionMatrix_w = identity(learnedSize);
			for (int n=0; n<i; n++) {
				transitionMatrix_w = multiply(transitionMatrix_w, hypothesisTransitionMatrices[letterToIndex.get(counterExampleArr[n])]);
			}
			
			// if F is the Hankel matrix, check if F_ω = sum(μ(ω)_1,i * F_xi)
			for (int j=0; j<learnedSize; j++) {
				int sum = 0;
				for (int k=0; k<learnedSize; k++) {
					sum = mod2(sum + getEntry(transitionMatrix_w, 1, k+1) * MQ(addStrings(learnedRowIndices.get(k), learnedColIndices.get(j))));
				}
				if (MQ(addStrings(w, learnedColIndices.get(j))) != sum) {
					break;
				}
			}
			
			// go through every possible value of y in learnedColIndices
			// check if F_{ω+σ}(y) != sum(μ(ω)_1,i * F_{xi+σ}(y))
			for (int j=0; j<learnedSize; j++) {
				y = learnedColIndices.get(j);
			
				int sum = 0;
				for (int k=0; k<learnedSize; k++) {
					sum = mod2(sum + getEntry(transitionMatrix_w, 1, k+1) * MQ(addStrings(addStrings(learnedRowIndices.get(k), sigma), y)));
				}
				
				// found a solution
				if (MQ(addStrings(addStrings(w, sigma), y)) != sum) {		
					if (learnedSize == minSize) {
						throwException(null, "Algorithm failed: size of the hypothesis exceeds that of the target function.");
					}

					learnedSize++;
					learnedRowIndices.add(w);
					learnedColIndices.add(addStrings(sigma, y));
					
					if (observationTableFlag) {
						displayTable();
					}
					
					return;
				}
			}
		}

		throwException(null, "Algorithm failed: didn't find a suitable omega, sigma, and gamma.");
	}

	public static void displayResults() {
		System.out.println("Learned mod-2-MA");
		System.out.println("----------------");
		
		System.out.println("Dimension: " + learnedSize + '\n');
		
		System.out.print("Final Vector: ");
		displayMatrix(resultFinalVector);
		
		System.out.println("Transition Matrices:\n");
		for (int i=0; i<resultTransitionMatrices.length; i++) {
			System.out.println("Letter " + alphabet[i]);
			displayMatrix(resultTransitionMatrices[i]);
		}
	}
	
	public static void displayTable() throws Exception {
		System.out.println("Size: " + learnedSize);
		System.out.print("Rows: ɛ ");
		for (int i=1; i<learnedRowIndices.size(); i++) {
			System.out.print(removeSpaces(learnedRowIndices.get(i)) + " ");
		}
		System.out.println();
		
		System.out.print("Cols: ɛ ");
		for (int i=1; i<learnedColIndices.size(); i++) {
			System.out.print(removeSpaces(learnedColIndices.get(i)) + " ");
		}
		
		System.out.println("\nTable:");
		for (int i=0; i<learnedRowIndices.size(); i++) {
			for (int j=0; j<learnedColIndices.size(); j++) {
				System.out.print(MQ(addStrings(learnedRowIndices.get(i), learnedColIndices.get(j))) + " ");
			}
			System.out.println();
		}
		System.out.println();
	}

	public static String genTest(int len, boolean smallerAlphabet) {
		int length = alphabet.length;
		if (smallerAlphabet) {
			length--;
		}
		
		String test = "";
		for (int i=0; i<len-1; i++) {
			test += alphabet[(int) (Math.random() * length)] + " ";
		}
		if (len >= 1) {
			test +=  alphabet[(int) (Math.random() * length)];
		}
		return test;
	}

	// performs a statistical EQ between the target and final mod-2-MA
	public static boolean finalCheck(int maxTestLen, int numTests, boolean inMinimize) throws Exception {
		for (int i=1; i<=numTests; i++) {
			String test = genTest((int) (Math.random() * (maxTestLen + 1)), false);
			
			if (!inMinimize && MQArbitrary(inputFinalVector, inputTransitionMatrices, test) != MQArbitrary(resultFinalVector, resultTransitionMatrices, test)) {
				return false;
			} else if (inMinimize && M2MA.MQArbitrary(M2MA.inputFinalVector, M2MA.inputTransitionMatrices, test) != M2MA.MQArbitrary(M2MA.minFinalVector, M2MA.minTransitionMatrices, test)) {
				return false;
			}
		}
		return true;
	}
	
	// returns the number of states of a minimal DFA equivalent to the minimized M2MA
	// The number of states of the minimal DFA is the number of reachable states of the minimized M2MA.
	// Uses sparse ArrayLists instead of the HashMap representation for the tests to avoid out of memory issues.
	public static int dimensionMinDFA(boolean minAutomata) throws Exception {
		HashMap<Integer, ArrayList<Integer>> orderedColSets = obtainColumnSets(minAutomata);
		
		HashSet<ArrayList<Integer>> reachable = new HashSet<ArrayList<Integer>>();
		int sizeReachable = 0;
		
		ArrayList<ArrayList<Integer>> tests = new ArrayList<ArrayList<Integer>>();
		
		// begin with ω_i = (1,0,0,...,0)
		ArrayList<Integer> w_i = new ArrayList<Integer>();
		w_i.add(1);
		tests.add(w_i);
		int numTests = 1;
		
		while (numTests > 0) {
			ArrayList<Integer> test = tests.remove(0);
			numTests--;
			
			if (!reachable.contains(test)) {
				reachable.add(test);
				sizeReachable++;
				
				// add to tests the one-letter extensions of test
				for (int i=0; i<alphabet.length; i++) {
					if (minAutomata) {
						tests.add(multiplyDimensionMinDFA(test, minTransitionMatrices[i], orderedColSets.get(i)));
					} else {
						tests.add(multiplyDimensionMinDFA(test, resultTransitionMatrices[i], orderedColSets.get(i)));
					}
					numTests++;
				}
			}
		}
		
		return sizeReachable;
	}
	
	public static HashMap<Integer, ArrayList<Integer>> obtainColumnSets(boolean minAutomata) {
		HashMap<Integer, ArrayList<Integer>> orderedColSets = new HashMap<Integer, ArrayList<Integer>>();
		
		for (int i = 0; i < alphabet.length; i++) {
			ArrayList<Integer> colSet = new ArrayList<Integer>();
			int prevMax = 0;
			boolean needsSorting = false;
			
			HashMap<Integer, ArrayList<Integer>> transitionMatrix;
			
			if (minAutomata) {
				transitionMatrix = minTransitionMatrices[i];
			} else {
				transitionMatrix = resultTransitionMatrices[i];
			}
			
			for (int c : transitionMatrix.keySet()) {
				if (c < 0) {
					colSet.add(c);
					
					if (!needsSorting && c > prevMax) {
						needsSorting = true;
					}
					
					prevMax = c;
				}
			}
			
			if (needsSorting) {
				Collections.sort(colSet);
			}
			
			orderedColSets.put(i, colSet);
		}
		
		return orderedColSets;
	}
	
	public static ArrayList<Integer> multiplyDimensionMinDFA(ArrayList<Integer> test, HashMap<Integer, ArrayList<Integer>> transitionMatrix, ArrayList<Integer> colSet) {
		ArrayList<Integer> out = new ArrayList<Integer>();
		
		for (int c : colSet) {
			if (dotProduct(test, transitionMatrix.get(c)) == 1) {
				out.add(-1*c);
			}
		}
		
		return out;
	}
	
	public static void displayRuntime() {
		long endTime = System.nanoTime();
		double totalTime = (endTime - startTime) / Math.pow(10, 9);
		DecimalFormat df = new DecimalFormat("0.00");
		System.out.println("Ran in " + df.format(totalTime) + "s.\n");
	}
	
	// performs operations on the learned mod-2-MA
	public static void operationsOnLearnedMA() throws Exception {
		System.out.println("Input a word to test if it is accepted.");
		System.out.println("Words are space-separated strings of letters (e.g. for the alphabet {a, b0}, a word is \"a b0 b0\").");
		System.out.println("If the language is (L)_$, words must be of the form u$v.");
		System.out.println("Enter \"quit\" to terminate.");
		
		while (true) {
			// read in cmd
			String test = in.nextLine().trim();
			
			// terminate the program
			if (test.equals("quit")) {
				in.close();
				break;
			}
			
			// test whether a word is accepted
			// if no word is passed, the test is the empty string		
			if (!inAlphabet(test)) {
				System.out.println("Inputted word is not in the language.");
			} else if(MQArbitrary(resultFinalVector, resultTransitionMatrices, test) == 1) {
				System.out.println("Accepted");
			} else {
				System.out.println("Not accepted");
			}
		}
	}
	
	public static boolean inAlphabet(String word) {
		String[] wordArr = word.split(" ");
		if (word.length() == 0) {
			return true;
		}
		for (int i=0; i<wordArr.length; i++) {
			if (letterToIndex.get(wordArr[i]) == null) {
				return false;
			}
		}
		return true;
	}
	
	public static int mod2(double n) {
		int temp = (int) Math.round(n);
		if (temp%2 == 0) {
			return 0;
		} else {
			return 1;
		}
	}
	
	// adds two strings with a white space in between (trimmed)
	public static String addStrings(String s1, String s2) {
		String out = s1;

		if (s2.length() > 0) {
			if (s1.length() > 0) {
				out += " ";
			}
			out += s2;
		}
		return out;
	}

	public static String removeSpaces(String s) {
		String out = "";
		for (int i=0; i<s.length(); i++) {
			char c = s.charAt(i);
			
			if (c != ' ') {
				out += c;
			}
		}
		return out;
	}
	
	/* Sparse matrix operations. */
	
	public static HashMap<Integer, ArrayList<Integer>> initialize(int numRows, int numCols) {
		HashMap<Integer, ArrayList<Integer>> arr = new HashMap<Integer, ArrayList<Integer>>();
		ArrayList<Integer> dim = new ArrayList<Integer>();
		dim.add(numRows);
		dim.add(numCols);
		arr.put(0, dim);
		return arr;
	}
	
	// multiplies a sparse mxn matrix by a sparse nxp matrix
	public static HashMap<Integer, ArrayList<Integer>> multiply(HashMap<Integer, ArrayList<Integer>> arr1, HashMap<Integer, ArrayList<Integer>> arr2) throws Exception {
		if ((int) arr1.get(0).get(1) != (int) arr2.get(0).get(0)) {
			throwException(null, "Multiplied matrices of invalid dimension.");
		}

		HashMap<Integer, ArrayList<Integer>> out = initialize(arr1.get(0).get(0), arr2.get(0).get(1));
		ArrayList<Integer> colSet = null;
		
		int rowsTraversed = 0;
		for (int r : arr1.keySet()) {
			if (rowsTraversed >= arr1.get(0).get(0)) {
				break;
			}
			
			if (r > 0) {
				rowsTraversed++;
				
				if (colSet == null) {
					colSet = new ArrayList<Integer>();
					
					int colsTraversed = 0;
					for (int c : arr2.keySet()) {
						if (colsTraversed == arr2.get(0).get(1)) {
							break;
						}
						if (c < 0) {
							colsTraversed++;
							colSet.add(c);
							
							if (dotProduct(arr1.get(r), arr2.get(c)) == 1) {
								addElement(out, r, c * -1);
							}
						}
					}
				} else {
					for (int c : colSet) {
						if (dotProduct(arr1.get(r), arr2.get(c)) == 1) {
							addElement(out, r, c * -1);
						}
					}
				}
			}
		}
		
		return out;
	}
	
	// returns the dot product of two sparse vectors
	public static int dotProduct(ArrayList<Integer> v1, ArrayList<Integer> v2) {
		if (v1 == null || v2 == null) {
			return 0;
		}
		
		HashSet<Integer> set = new HashSet<Integer>();
		for (int n : v1) {
			set.add(n);
		}
		
		int count = 0;
		
		for (int n : v2) {
			if (set.contains(n)) {
				count++;
			}
		}
		
		return count % 2;
	}
	
	// creates the identity matrix
	public static HashMap<Integer, ArrayList<Integer>> identity(int size) throws Exception {
		HashMap<Integer, ArrayList<Integer>> identity = initialize(size, size);
		for (int i=1; i<=size; i++) {
			addElement(identity, i, i);
		}
		return identity;
	}
	
	// adds a 1 to a sparse array at the position (row, col)
	public static void addElement(HashMap<Integer, ArrayList<Integer>> arr, int row, int col) throws Exception {
		if (row < 1 || row > arr.get(0).get(0) || col < 1 || col > arr.get(0).get(1)) {
			throwException(null, "Added an invalid element to a matrix.");
		}
		if (arr.get(row) == null) {
			ArrayList<Integer> newRow = new ArrayList<Integer>();
			newRow.add(col);
			arr.put(row, newRow);
		} else {
			arr.get(row).add(col);
		}
		
		if (arr.get(col * -1) == null) {
			ArrayList<Integer> newCol = new ArrayList<Integer>();
			newCol.add(row);
			arr.put(col * -1, newCol);
		} else {
			arr.get(col * -1).add(row);
		} 
	}
	
	// removes a 1 from a sparse array at the position (row, col)
	public static void removeElement(HashMap<Integer, ArrayList<Integer>> arr, int row, int col) throws Exception {
		if (row < 1 || row > arr.get(0).get(0) || col < 1 || col > arr.get(0).get(1) || arr.get(row) == null || arr.get(col * -1) == null) {
			throwException(null, "Attempted to remove an invalid element from a matrix.");
		}
		
		arr.get(row).remove(Integer.valueOf(col));
		arr.get(col * -1).remove(Integer.valueOf(row));
	}
	
	// returns the element at position (row, col)
	public static int getEntry(HashMap<Integer, ArrayList<Integer>> arr, int row, int col) {
		if (arr.get(row) == null) {
			return 0;
		}
		
		for (int num : arr.get(row)) {
			if (num == col) {
				return 1;
			}
		}
		
		return 0;
	}
	
	// displays a sparse matrix
	public static void displayMatrix(HashMap<Integer, ArrayList<Integer>> arr) {
		ArrayList<Integer> dim = arr.get(0);
		
		if (dim.get(0) == 0 || dim.get(1) == 0) {
			System.out.println();
			return;
		}
		
		for (int row = 1; row <= dim.get(0); row++) {
			String rowStr = "";
			
			if (arr.get(row) == null) {
				for (int col = 1; col <= dim.get(1); col++) {
					rowStr += 0 + " ";
				}
			} else {
				ArrayList<Integer> rowVector = arr.get(row);
				int pos = 0;
				for (int col = 1; col <= dim.get(1); col++) {
					if (pos < rowVector.size() && rowVector.get(pos) == col) {
						rowStr += 1 + " ";
						pos++;
					} else {
						rowStr += 0 + " ";
					}
				}
			}
			
			System.out.println(rowStr);
		}
		System.out.println();
	}

	// converts a sparse matrix into a real matrix
	public static RealMatrix sparseToReal(HashMap<Integer, ArrayList<Integer>> arr) {
		RealMatrix out = MatrixUtils.createRealMatrix(arr.get(0).get(0), arr.get(0).get(1));
		
		for (int row : arr.keySet()) {
			if (row > 0) {
				for (int col : arr.get(row)) {
					out.setEntry(row - 1, col - 1, 1);
				}
			}
		}
		
		return out;
	}
	
	// converts a real matrix into a sparse matrix
	public static HashMap<Integer, ArrayList<Integer>> realToSparse(RealMatrix arr) throws Exception {
		HashMap<Integer, ArrayList<Integer>> out = initialize(arr.getRowDimension(), arr.getColumnDimension());
		
		for (int r = 0; r < arr.getRowDimension(); r++) {
			for (int c = 0; c < arr.getColumnDimension(); c++) {
				if ((int) arr.getEntry(r, c) == 1) {
					addElement(out, r + 1, c + 1);
				}
			}
		}
		
		return out;
	}
	
	// converts a int[][] into a sparse matrix
	public static HashMap<Integer, ArrayList<Integer>> intArrayToSparse(int[][] arr) throws Exception {
		HashMap<Integer, ArrayList<Integer>> out = initialize(arr.length, arr[0].length);
		
		for (int r = 0; r < arr.length; r++) {
			for (int c = 0; c < arr[0].length; c++) {
				if ((int) arr[r][c] == 1) {
					addElement(out, r + 1, c + 1);
				}
			}
		}
		
		return out;
	}
}