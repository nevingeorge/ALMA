import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class convert {
	
	public static boolean inSUBA;
	public static boolean inNBA;
	public static boolean inM2MA;
	
	public static int[][] results;
	
	public static BufferedReader f;
	
	public static void main(String[] args) throws Exception {
		readInput();
		
		if (inSUBA) {
			convertSUBA();
		} else {
			convertNBA();
		}
		
		f.close();
		
		displayResults();
	}
	
	public static void readInput() throws Exception {
		inSUBA = false;
		inNBA = false;
		inM2MA = false;
		
		M2MA.in = new Scanner(System.in);
		
		System.out.println("Enter 1 to convert SUBA, 2 to convert NBA, and 3 to convert DBA.");
		int automataToConvert = Integer.parseInt(M2MA.in.nextLine());
		
		if (automataToConvert == 1) {
			inSUBA = true;
		} else if (automataToConvert == 2) {
			inNBA = true;
		} else if (automataToConvert != 3) {
			M2MA.throwException(null, "Invalid number entered.");
		}
		
		System.out.println("Enter 1 to convert to a M2MA, and 2 to convert to a DFA.");
		int convertToAutomata = Integer.parseInt(M2MA.in.nextLine());
		
		if (convertToAutomata == 1) {
			inM2MA = true;
		} else if (convertToAutomata != 2) {
			M2MA.throwException(null, "Invalid number entered.");
		}
		
		System.out.println("Input file name.");
		f = new BufferedReader(new FileReader(M2MA.in.nextLine()));
		M2MA.in.close();
		
		// automata of size 1 to 30, int[][0] = count, int[][1] = sum of converted automata dimension
		results = new int[31][2];
	}
	
	public static void convertSUBA() throws Exception {
		int numSUBA = Integer.parseInt(M2MA.readFile(f));
		
		int prevPercentComplete = 0;
		System.out.println("0% complete.");
		
		for (int i = 0; i < numSUBA; i++) {
			int percentComplete = (int) ((((double) i) / numSUBA) * 100);
			if (percentComplete >= prevPercentComplete + 5) {
				System.out.println(percentComplete + "% complete.");
				prevPercentComplete = percentComplete;
			}
			
			SUBA.SUBAtoUFA();
			
			SUBA.UFAtoMod2MA();
			
			M2MA.minimize();

			results[SUBA.SUBAStates][0]++;
			
			if (inM2MA) {
				results[SUBA.SUBAStates][1] += M2MA.minSize;
			} else {
				int dim = M2MA.dimensionMinDFA(true);
				results[SUBA.SUBAStates][1] += dim;
				System.out.println(dim);
			}
		}
		System.out.println("100% complete.\n");
	}

	public static void convertNBA() throws Exception {
		arbitrary.EQMaxTestLen = Integer.parseInt(M2MA.readFile(f));	
		arbitrary.EQNumTests = Integer.parseInt(M2MA.readFile(f));
		arbitrary.EQLimit = Integer.parseInt(M2MA.readFile(f));
		
		// alphabet ΣU{$}
		M2MA.readAlphabet(f, true);		
		
		int numLines = Integer.parseInt(M2MA.readFile(f));	
		
		int prevPercentComplete = 0;
		System.out.println("0% complete.");
		
		for (int i = 0; i < numLines; i++) {			
			int percentComplete = (int) ((((double) i) / numLines) * 100);
			if (percentComplete >= prevPercentComplete + 5) {
				System.out.println(percentComplete + "% complete.");
				prevPercentComplete = percentComplete;
			}
			
			String[] line = M2MA.readFile(f).split(" ");
			int numNBA = Integer.parseInt(line[0]);
			int numStates = Integer.parseInt(line[1]);
			int numTransitionsToRemove = Integer.parseInt(line[2]);
			int numFinal = Integer.parseInt(line[3]);
			
			runNBA(numNBA, numStates, numTransitionsToRemove, numFinal);
		}
		System.out.println("100% complete.\n");
	}
	
	@SuppressWarnings("unchecked")
	public static void runNBA(int numTimes, int numStates, int numTransitionsToRemove, int numFinalStates) throws Exception {
		// for every reachable state add it to an array (also boolean array checkoff). Length of the array ends up being NBAStates.
		for (int i = 0; i < numTimes; i++) {
			arbitrary.EQNumPerformed = 0;
			
			ArrayList<int[]>[] tempTransitions = new ArrayList[numStates + 1];
			ArrayList<int[]>[] reverseTempTransitions = new ArrayList[numStates + 1];
			for (int j = 1; j <= numStates; j++) {
				tempTransitions[j] = new ArrayList<int[]>();
				reverseTempTransitions[j] = new ArrayList<int[]>();
			}
			
			if (inNBA) {
				NBAtransitions(numStates, numTransitionsToRemove, tempTransitions, reverseTempTransitions);
			} else {
				DBAtransitions(numStates, numTransitionsToRemove, tempTransitions, reverseTempTransitions);
			}
			
			// find all reachable states from state 1
			HashSet<Integer> reachableStates = new HashSet<Integer>();
			reachable(1, tempTransitions, reachableStates);
			
			HashSet<Integer> tempFinalStates = new HashSet<Integer>();
			finalStates(reachableStates, numFinalStates, tempFinalStates);
			
			// Find all the reachable states from every final state, reverse-deterministically.
			// If it can't reach itself, remove its final property.
			// If it reaches itself, add all the reachable states to finalReachableStates.
			HashSet<Integer> finalReachableStates = new HashSet<Integer>();
			HashSet<Integer> finalFinalStates = new HashSet<Integer>();
			
			for (int finalState : tempFinalStates) {
				HashSet<Integer> tempFinalReachableStates = new HashSet<Integer>();
				
				finalReachable(true, false, finalState, finalState, finalReachableStates, tempFinalReachableStates, reachableStates, reverseTempTransitions);
				
				if (tempFinalReachableStates.contains(finalState)) {
					finalFinalStates.add(finalState);
				}
			}			
			
			if (finalReachableStates.size() != 0) {
				NBA.NBAStates = finalReachableStates.size();
				
				HashMap<Integer, Integer> tempToRealState = new HashMap<Integer, Integer>();
				boolean[] tempStatesFound = new boolean[numStates + 1];
				
				for (int state : finalReachableStates) {
					tempStatesFound[state] = true;
				}
				
				int count = 1;
				for (int j = 1; j <= numStates; j++) {
					if (tempStatesFound[j]) {
						tempToRealState.put(j, count);
						count++;
					}
				}
				
				NBA.NBAFinalStates = new boolean[NBA.NBAStates + 1];
				for (int state : finalFinalStates) {
					NBA.NBAFinalStates[tempToRealState.get(state)] = true;
				}
				
				NBA.NBATransitions = new ArrayList[NBA.NBAStates + 1][M2MA.alphabet.length - 1];
				for (int j = 1; j <= NBA.NBAStates; j++) {
					for (int k = 0; k < M2MA.alphabet.length - 1; k++) {
						NBA.NBATransitions[j][k] = new ArrayList<Integer>();
					}
				}
				
				for (int state1 = 1; state1 <= numStates; state1++) {
					if (finalReachableStates.contains(state1)) {
						for (int[] transition : tempTransitions[state1]) {
							if (finalReachableStates.contains(transition[1])) {
								NBA.NBATransitions[tempToRealState.get(state1)][transition[0]].add(tempToRealState.get(transition[1]));
							}
						}
					}
				}
				
				M2MA.Hankel = null;
				M2MA.learn();
				
				results[NBA.NBAStates][0]++;
				if (inM2MA) {
					results[NBA.NBAStates][1] += M2MA.learnedSize;
				} else {
					results[NBA.NBAStates][1] += M2MA.dimensionMinDFA(false);
				}
			}
		}
	}
	
	public static void finalReachable(boolean firstState, boolean passedFinalState, int startState, int finalState, HashSet<Integer> finalReachableStates, HashSet<Integer> tempFinalReachableStates, HashSet<Integer> reachableStates, ArrayList<int[]>[] reverseTempTransitions) {
		if (!firstState && startState == finalState && !passedFinalState) {
			passedFinalState = true;
			
			for (int state : tempFinalReachableStates) {
				finalReachableStates.add(state);
			}
		}
		
		if (passedFinalState) {
			finalReachableStates.add(startState);
		}
		
		if (!firstState) {
			tempFinalReachableStates.add(startState);
		}
		
		for (int[] reverseTransition : reverseTempTransitions[startState]) {
			if (reachableStates.contains(reverseTransition[1]) && !tempFinalReachableStates.contains(reverseTransition[1])) {
				finalReachable(false, passedFinalState, reverseTransition[1], finalState, finalReachableStates, tempFinalReachableStates, reachableStates, reverseTempTransitions);
			}
		}
	}
	
	public static void NBAtransitions(int numStates, int numTransitionsToRemove, ArrayList<int[]>[] tempTransitions, ArrayList<int[]>[] reverseTempTransitions) {
		boolean[][][] unusedTransitions = new boolean[numStates + 1][M2MA.alphabet.length - 1][numStates + 1];
		
		for (int i = 0; i < numTransitionsToRemove; i++) {
			int state1 = (int) (Math.random() * numStates) + 1;
			int letter = (int) (Math.random() * (M2MA.alphabet.length - 1));
			int state2 = (int) (Math.random() * numStates) + 1;
			
			if (!unusedTransitions[state1][letter][state2]) {	
				unusedTransitions[state1][letter][state2] = true;
			}
		}
		
		for (int state1 = 1; state1 <= numStates; state1++) {
			for (int letter = 0; letter < M2MA.alphabet.length - 1; letter++) {
				for (int state2 = 1; state2 <= numStates; state2++) {
					if (!unusedTransitions[state1][letter][state2]) {
						int[] transition = new int[2];
						transition[0] = letter;
						transition[1] = state2;
						tempTransitions[state1].add(transition);
						
						int[] reverseTransition = new int[2];
						reverseTransition[0] = letter;
						reverseTransition[1] = state1;
						reverseTempTransitions[state2].add(reverseTransition);
					}
				}
			}
		}
	}
	
	public static void DBAtransitions(int numStates, int numTransitionsToRemove, ArrayList<int[]>[] tempTransitions, ArrayList<int[]>[] reverseTempTransitions) {
		boolean[][] unusedTransitions = new boolean[numStates + 1][M2MA.alphabet.length - 1];
		
		for (int i = 0; i < numTransitionsToRemove; i++) {
			int state1 = (int) (Math.random() * numStates) + 1;
			int letter = (int) (Math.random() * (M2MA.alphabet.length - 1));
			
			if (!unusedTransitions[state1][letter]) {
				unusedTransitions[state1][letter] = true;
			}
		}
		
		for (int state1 = 1; state1 <= numStates; state1++) {
			for (int letter = 0; letter < M2MA.alphabet.length - 1; letter++) {
				if (!unusedTransitions[state1][letter]) {
					int state2 = (int) (Math.random() * numStates) + 1;
					
					int[] transition = new int[2];
					transition[0] = letter;
					transition[1] = state2;
					tempTransitions[state1].add(transition);
					
					int[] reverseTransition = new int[2];
					reverseTransition[0] = letter;
					reverseTransition[1] = state1;
					reverseTempTransitions[state2].add(reverseTransition);
				}
			}
		}
	}
	
	public static void reachable(int startState, ArrayList<int[]>[] tempTransitions, HashSet<Integer> reachableStates) {
		reachableStates.add(startState);
		
		for (int[] transition : tempTransitions[startState]) {
			if (!reachableStates.contains(transition[1])) {
				reachable(transition[1], tempTransitions, reachableStates);
			}
		}
	}
	
	public static void finalStates(HashSet<Integer> reachableStates, int numFinalStates, HashSet<Integer> tempFinalStates) {	
		int size = reachableStates.size();
		numFinalStates = Math.min(size, numFinalStates);
		HashSet<Integer> randomSet = new HashSet<Integer>();
		
		for (int i = 0; i < numFinalStates; i++) {
			int randomNum = (int) (Math.random() * size);
			
			int count = 0;
			while (count < size) {
				if (!randomSet.contains((randomNum + count) % size)) {
					randomSet.add((randomNum + count) % size);
					break;
				}
				count++;
			}
		}
		
		int count = 0;
		int numAdded = 0;
		for (int state : reachableStates) {
			if (numAdded >= numFinalStates) {
				break;
			}
			
			if (randomSet.contains(count)) {
				tempFinalStates.add(state);
				numAdded++;
			}
			
			count++;
		}
	}

	public static void displayResults() {
		System.out.println("Results");
		System.out.println("-------");
		System.out.println("Number of states: number of converted automata, sum of converted automata sizes, average converted automata size\n");
		for (int i = 1; i <= 30; i++) {
			if (results[i][0] != 0) {
				System.out.println(i + ": " + results[i][0] + ", " + results[i][1] + ", " + ((double) results[i][1])/results[i][0]);
			}
		}
	}
}
