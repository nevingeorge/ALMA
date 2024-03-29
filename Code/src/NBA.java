/*
 * Author: Nevin George
 * Advisor: Dana Angluin
 * Program Description: The program takes in as input an NBA and prints to stdout the M2MA obtained after 
 * learning the NBA through a series of membership and statistical equivalence queries.
 */

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class NBA {
	
	// NBA
	public static int NBAStates;
	public static ArrayList<Integer>[][] NBATransitions;
	public static boolean[] NBAFinalStates;

	public static void main(String[] args) throws Exception {
		System.out.println("Program Description:");
		System.out.println("The program takes in as input an NBA and prints to stdout the M2MA obtained after"
				+ " learning the NBA through a series of membership and statistical equivalence queries.");
		System.out.println("\n----------------------------\n");
		
		readInput();

		M2MA.learn();
		
		M2MA.displayResults();
		
		M2MA.displayRuntime();
		
		M2MA.operationsOnLearnedMA();
	}
	
	@SuppressWarnings("unchecked")
	public static void readInput() throws Exception {
		System.out.println("Input file name and optional flag -v (e.g. NBA_input1.txt or NBA_input1.txt -v)");
		
		BufferedReader f = M2MA.getFile(true, false, false, false);
		
		arbitrary.EQMaxTestLen = Integer.parseInt(M2MA.readFile(f));	
		arbitrary.EQNumTests = Integer.parseInt(M2MA.readFile(f));
		arbitrary.EQLimit = Integer.parseInt(M2MA.readFile(f));
		arbitrary.EQNumPerformed = 0;

		NBAStates = Integer.parseInt(M2MA.readFile(f));
		
		// alphabet ΣU{$}
		M2MA.readAlphabet(f, true);
		
		StringTokenizer st = new StringTokenizer(M2MA.readFile(f));
		NBAFinalStates = new boolean[NBAStates+1];
		while (st.hasMoreTokens()) {
			int state = Integer.parseInt(st.nextToken());
			if(1<=state && state<=NBAStates && !NBAFinalStates[state]) {
				NBAFinalStates[state] = true;
			} else {
				M2MA.throwException(f, "Invalid input: invalid or duplicate final state.");
			}
		}
		
		int numTransitions = Integer.parseInt(M2MA.readFile(f));
		if(numTransitions > ((M2MA.alphabet.length - 1) * NBAStates * NBAStates)) {
			M2MA.throwException(f, "Invalid input: invalid number of transitions.");
		}
		
		// (start state, letter, end state)
		NBATransitions = new ArrayList[NBAStates+1][M2MA.alphabet.length-1];
		for (int i=1; i<=NBAStates; i++) {
			for (int j=0; j<M2MA.alphabet.length-1; j++) {
				NBATransitions[i][j] = new ArrayList<Integer>();
			}
		}
		
		// lines of the form q_i a q_j, where q_i,q_j∈NBAStates and a∈alphabet
		for (int i=0; i<numTransitions; i++) {
			st = new StringTokenizer(M2MA.readFile(f));
			int p_start = Integer.parseInt(st.nextToken());
			
			String letter = st.nextToken();
			
			int a = M2MA.letterToIndex.get(letter);
			int p_end = Integer.parseInt(st.nextToken());
			if (p_start < 1 || p_start > NBAStates || p_end < 1 || p_end > NBAStates) {
				M2MA.throwException(f, "Invalid input: invalid transition.");
			}

			NBATransitions[p_start][a].add(p_end);
		}
		
		if (M2MA.readFile(f) != null) {
			M2MA.throwException(f, "Invalid input: more transitions inputted than specified.");
		}
		
		f.close();
	}
	
	public static int MQ(String w) {
		// ω must contain exactly one $
		int dollarIndex = -1;
		for (int i=0; i<w.length(); i++) {
			if (w.charAt(i) == '$' && dollarIndex == -1) {
				dollarIndex = i;
			} else if(w.charAt(i) == '$' && dollarIndex != -1) {
				return 0;
			}
		}
		// $ must appear in ω and cannot be at the final index (otherwise the periodic string v is empty)
		if (dollarIndex == -1 || dollarIndex == w.length()-1) {
			return 0;
		}
		
		String u = w.substring(0, dollarIndex).trim();
		String v = w.substring(dollarIndex+1).trim();
		
		/*
		 * States are represented as length-2 integer arrays.
		 * The first element is the state number.
		 * The second element is either 0 or 1, with 1 indicating that a final state was passed to reach the state (or if the state itself is final).
		 */
		
		// initial state by default is q_1
		ArrayList<int[]> initialState = new ArrayList<int[]>();
		int[] state = {1, 0};
		initialState.add(state);
		
		// states reachable from initialState on one u
		ArrayList<int[]> readU = readStr(u, initialState);
		
		// states reachable from initialState on one u and a non-negative number of v's
		ArrayList<int[]> reachable = readV(v, readU);
		reachable.addAll(readU);
		refine(reachable);
		
		// see if any of the states in reachable have an accepting loop
		for (int i=0; i<reachable.size(); i++) {
			if (acceptingLoop(reachable.get(i), v)) {
				return 1;
			}
		}
		return 0;
	}
	
	// returns true if, starting at state, after reading some positive number of v's it returns to state while having passed a final state
	public static boolean acceptingLoop(int[] state, String v) {
		if(NBAFinalStates[state[0]]) {
			state[1] = 1;
		} else {
			state[1] = 0;
		}
		ArrayList<int[]> initialState = new ArrayList<int[]>();
		initialState.add(state);
		
		// find all of the states reachable from state on a positive number of v's
		ArrayList<int[]> reachable = readV(v, initialState);
		
		// see if reachable contains (state, 1), which indicates an accepting loop
		for (int i=0; i<reachable.size(); i++) {
			if((state[0] == reachable.get(i)[0]) && (reachable.get(i)[1] == 1)) {
				return true;
			}
		}
		return false;
	}
	
	// returns an ArrayList with all of the states reachable from a state in states on a positive number of v's
	public static ArrayList<int[]> readV(String v, ArrayList<int[]> states) {
		// read one v (must read a positive number of v's)
		ArrayList<int[]> reachable = readStr(v, states);
		
		// keep reading v's until we find no new states upon reading another v
		while (true)	{
			ArrayList<int[]> readNextV = readStr(v, reachable);
			
			// check if readNextV has states not found in reachable
			if (subset(readNextV, reachable)) {
				return reachable;
			}
			refine(reachable);
		}
	}
	
	// returns an ArrayList with all of the states reachable from startStates on str
	public static ArrayList<int[]> readStr(String str, ArrayList<int[]> startStates) {
		if (str.length() == 0) {
			return refine(startStates);
		}
		
		ArrayList<int[]> nextStates = new ArrayList<int[]>();
		boolean[][] visited = new boolean[NBAStates+1][2];
		
		String firstLetter = "";
		int pos = 0;
		while (pos < str.length() && str.charAt(pos) != ' ') {
			firstLetter += str.charAt(pos++);
		}
		
		// update nextStates with all of the states reachable from startStates on the first letter of str
		for (int i=0; i<startStates.size(); i++) {
			ArrayList<Integer> curTransition = NBATransitions[startStates.get(i)[0]][M2MA.letterToIndex.get(firstLetter)];
			for (int j=0; j<curTransition.size(); j++) {
				int nextState = curTransition.get(j);
				int mark = 0;
				
				// either came from a state that has passed a final state or is a final state itself
				if (startStates.get(i)[1] == 1 || NBAFinalStates[nextState]) {
					mark = 1;
				}
				
				if (!visited[nextState][1] && (!visited[nextState][0] || mark == 1)) {
					visited[nextState][mark] = true;
					int[] markedState = {nextState, mark};
					nextStates.add(markedState);
				}
			}
		}
		
		return readStr(str.substring(Math.min(firstLetter.length() + 1, str.length())), nextStates);
	}
	
	// refines states so that it contains unique values 
	// if both (state, 0) and (state, 1) are in states, only (state, 1) is kept
	public static ArrayList<int[]> refine(ArrayList<int[]> states) {
		// nothing to refine, makes reading u="" more efficient
		if (states.size() <= 1) {
			return states;
		}
		
		// place the elements of states in a boolean array
		boolean[][] seen = new boolean[NBAStates+1][2];
		for (int i=0; i<states.size(); i++) {
			if (states.get(i)[1] == 0) {
				seen[states.get(i)[0]][0] = true; 
			} else {
				seen[states.get(i)[0]][1] = true; 
			}
		}
		
		// add unique, "maximal" states to out
		ArrayList<int[]> out = new ArrayList<int[]>();
		for (int i=1; i<=NBAStates; i++) {
			if (seen[i][1]) {
				int[] markedState = {i,1};
				out.add(markedState);
			} else if (seen[i][0]) {
				int[] markedState = {i,0};
				out.add(markedState);
			}
		}
		return out;
	}

	// returns true if arr1 is a subset of arr2
	// adds the elements in arr1 that are not in arr2 to arr2
	public static boolean subset(ArrayList<int[]> arr1, ArrayList<int[]> arr2) {
		boolean out = true;
		
		// place the elements of arr1 in a boolean array
		boolean[][] inArr2 = new boolean[NBAStates+1][2];
		for(int i=0; i<arr2.size(); i++) {
			if (arr2.get(i)[1] == 0) {
				inArr2[arr2.get(i)[0]][0] = true;
			} else {
				inArr2[arr2.get(i)[0]][1] = true;
			}
		}
		
		// arr2 = arr1 U arr2
		for (int i=0; i<arr1.size(); i++) {
			if (!inArr2[arr1.get(i)[0]][arr1.get(i)[1]]) {
				arr2.add(arr1.get(i));
				out = false;
			}
		}
		return out;
	}
}
