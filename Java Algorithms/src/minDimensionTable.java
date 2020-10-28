import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.StringTokenizer;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class minDimensionTable {

	public static void main(String[] args) throws Exception {
		System.out.println("Input file name.");;
		Scanner in = new Scanner(System.in);
		BufferedReader f = new BufferedReader(new FileReader(in.nextLine()));
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter("dimensionTable.txt")));
		out.println("SUBA Dimension, Minimized Dimension");
		System.out.println("SUBA Dimension, Minimized Dimension");
		
		int numSUBA = Integer.parseInt(Mod2_MA.readFile(f));
		
		for (int i = 0; i < numSUBA; i++) {
			readInput(f);
			SUBA.UFAtoMod2MA();
			int minimizedDimension = minimization();
			out.println(SUBA.SUBAStates + " " + minimizedDimension);
			System.out.println(SUBA.SUBAStates + " " + minimizedDimension);
		}
		
		out.close();
		in.close();
	}
	
	@SuppressWarnings("unchecked")
	public static void readInput(BufferedReader f) throws Exception {
		// UFAStates = SUBAStates U (SUBAStates x SUBAStates x {0,1})
		SUBA.SUBAStates = Integer.parseInt(Mod2_MA.readFile(f));
		SUBA.UFAStates = SUBA.SUBAStates + SUBA.SUBAStates * SUBA.SUBAStates * 2;
		
		// alphabet ΣU{$}
		StringTokenizer st = new StringTokenizer(Mod2_MA.readFile(f));
		ArrayList<String> tempAlphabet = new ArrayList<String>();
		while (st.hasMoreTokens()) {
			String letter = st.nextToken();
			if(letter.equals("$")) {
				Mod2_MA.throwException(f, "Invalid input: invalid character in the alphabet.");
			}
			tempAlphabet.add(letter);
		}
		tempAlphabet.add("$");
		Mod2_MA.alphabet = new String[tempAlphabet.size()];
		for (int i=0; i<tempAlphabet.size(); i++) {
			Mod2_MA.alphabet[i] = tempAlphabet.get(i);
		}
		
		// map each letter in alphabet to an index
		Mod2_MA.letterToIndex = new HashMap<String, Integer>();
		for (int i=0; i<Mod2_MA.alphabet.length; i++) {
			Mod2_MA.letterToIndex.put(Mod2_MA.alphabet[i], i);
		}
		
		st = new StringTokenizer(Mod2_MA.readFile(f));
		SUBA.SUBAFinalStates = new boolean[SUBA.SUBAStates+1];
		while (st.hasMoreTokens()) {
			int state = Integer.parseInt(st.nextToken());
			if (1 <= state && state <= SUBA.SUBAStates && !SUBA.SUBAFinalStates[state]) {
				SUBA.SUBAFinalStates[state] = true;
			} else {
				Mod2_MA.throwException(f,"Invalid input: invalid or duplicate final state.");
			}
		}
		
		/* 
		 * Following the paper by Bousquet and Löding, UFATransitions contains (where q,p,p'∈SUBAStates)
		 * - SUBATransitions
		 * - all transitions of the form (q,$,(q,q,0))
		 * - all transitions of the form ((q,p,i),a,(q,p',i')), where (p,a,p')∈SUBATransitions, and
		 * 	 i' = 1 if p'∈SUBAFinalStates and i if p'∉SUBAFinalStates
		 * 
		 * Transitions will be stored in a (UFAStates x alphabetSize x UFAStates) adjacency matrix.
		 * The first states of UFATransitions will be SUBAStates.
		 * The remaining states will be of the form (q_j,q_k,i), where q_j,q_k∈SUBAStates and i∈{0,1}.
		 * State (q_j,q_k,i) will be found at index (2*SUBAStates*j)+(2*k)-(SUBAStates)+(i-1) of UFATransitions.
		*/
		int numTransitions = Integer.parseInt(Mod2_MA.readFile(f));
		if (numTransitions > ((Mod2_MA.alphabet.length - 1) * SUBA.SUBAStates * SUBA.SUBAStates)) {
			Mod2_MA.throwException(f, "Invalid input: invalid number of transitions.");
		}
		
		/*
		 *  For each index (q,a) where q∈SUBAStates and a∈alphabet, transition_SUBA[q][a] is an ArrayList containing 
		 *  all of the reachable states from (q,a).
		 *  The alphabet for the SUBA does not include $.
		 */
		SUBA.SUBATransitions = new ArrayList[SUBA.SUBAStates + 1][Mod2_MA.alphabet.length - 1];
		for (int i=1; i<=SUBA.SUBAStates; i++) {
			for (int j=0; j<Mod2_MA.alphabet.length-1; j++) {
				SUBA.SUBATransitions[i][j] = new ArrayList<Integer>();
			}
		}
		
		// (start state, letter, end state)
		SUBA.UFATransitions = new boolean[SUBA.UFAStates + 1][Mod2_MA.alphabet.length][SUBA.UFAStates + 1];
		
		// lines of the form q_j a q_k, where q_j,q_k∈SUBAStates and a∈alphabet
		for (int i=0; i<numTransitions; i++) {
			st = new StringTokenizer(Mod2_MA.readFile(f));
			int p_start = Integer.parseInt(st.nextToken());
			
			String letter = st.nextToken();
			
			int a = Mod2_MA.letterToIndex.get(letter);
			int p_end = Integer.parseInt(st.nextToken());
			if (p_start < 1 || p_start > SUBA.SUBAStates || p_end < 1 || p_end > SUBA.SUBAStates) {
				Mod2_MA.throwException(f, "Invalid input: invalid transition.");
			}
			
			// SUBATransitions ⊆ UFATransitions 
			SUBA.SUBATransitions[p_start][a].add(p_end);
			SUBA.UFATransitions[p_start][a][p_end] = true;
			
			// transitions of the form ((q,p,i),a,(q,p',i'))
			// p'∈SUBAFinalStates so i'=1
			if (SUBA.SUBAFinalStates[p_end]) {
				for (int q=1; q<=SUBA.SUBAStates; q++) {
					// ((q,p,0),a,(q,p',1))
					SUBA.UFATransitions[SUBA.getIndex(q, p_start, 0)][a][SUBA.getIndex(q, p_end, 1)] = true;
					// ((q,p,1),a,(q,p',1))
					SUBA.UFATransitions[SUBA.getIndex(q, p_start, 1)][a][SUBA.getIndex(q, p_end, 1)] = true;
				}
			}
			// p'∉SUBAFinalStates so i'=i
			else {
				for (int q=1; q<=SUBA.SUBAStates; q++) {
					// ((q,p,0),a,(q,p',0))
					SUBA.UFATransitions[SUBA.getIndex(q, p_start, 0)][a][SUBA.getIndex(q, p_end, 0)] = true;
					// ((q,p,1),a,(q,p',1))
					SUBA.UFATransitions[SUBA.getIndex(q, p_start, 1)][a][SUBA.getIndex(q, p_end, 1)] = true;
				}
			}
		}
		
		// transitions for the UFA of the form (q,$,(q,q,0)), where q∈SUBAStates
		// final states for the UFA of the form (q,q,1), where q∈SUBAStates
		SUBA.UFAFinalStates = new boolean[SUBA.UFAStates+1];
		for (int q=1; q<=SUBA.SUBAStates; q++) {
			SUBA.UFATransitions[q][Mod2_MA.letterToIndex.get("$")][SUBA.getIndex(q, q, 0)] = true;
			SUBA.UFAFinalStates[SUBA.getIndex(q, q, 1)] = true;
		}
	}

	public static int minimization() {
		ArrayList<String> stateSpaceBasisIndices = new ArrayList<String>();
		HashMap<String, double[]> stateSpaceIndexToVector = new HashMap<String, double[]>();
		RealMatrix stateSpaceBasis = Mod2_MA.sparseBasis(Mod2_MA.inputFinalVector, Mod2_MA.inputTransitionMatrices, stateSpaceIndexToVector, stateSpaceBasisIndices, true);
		
		ArrayList<String> coStateSpaceBasisIndices = new ArrayList<String>();
		HashMap<String, double[]> coStateSpaceIndexToVector = new HashMap<String, double[]>();
		RealMatrix coStateSpaceBasis = Mod2_MA.sparseBasis(Mod2_MA.inputFinalVector, Mod2_MA.inputTransitionMatrices, coStateSpaceIndexToVector, coStateSpaceBasisIndices, false);
		
		RealMatrix observationTable = MatrixUtils.createRealMatrix(stateSpaceBasis.getRowDimension(), coStateSpaceBasis.getRowDimension());
		for (int row=0; row<stateSpaceBasis.getRowDimension(); row++) {
			for (int col=0; col<coStateSpaceBasis.getRowDimension(); col++) {
				observationTable.setEntry(row, col, Mod2_MA.mod2(stateSpaceBasis.getRowVector(row).dotProduct(coStateSpaceBasis.getRowVector(col))));
			}
		}
		
		Mod2_MA.minRowIndices = new ArrayList<String>();
		RealMatrix linIndRowsObservationTable = Mod2_MA.sparseLinIndSubMatrix(observationTable, stateSpaceBasisIndices, Mod2_MA.minRowIndices, true);
		
		return linIndRowsObservationTable.getRowDimension();
	}
}
