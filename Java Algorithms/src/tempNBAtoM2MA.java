import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class tempNBAtoM2MA {
	
	public static int[][] results;
	
	public static void main(String[] args) throws Exception {
		arbitrary.EQMaxTestLen = 30;	
		arbitrary.EQNumTests = 1000;
		arbitrary.EQLimit = 100000;
		
		Mod2_MA.alphabet = new String[4];
		Mod2_MA.alphabet[0] = "a";
		Mod2_MA.alphabet[1] = "b";
		Mod2_MA.alphabet[2] = "c";
		Mod2_MA.alphabet[3]	= "$";			
		
		Mod2_MA.letterToIndex = new HashMap<String, Integer>();
		for (int i=0; i<Mod2_MA.alphabet.length; i++) {
			Mod2_MA.letterToIndex.put(Mod2_MA.alphabet[i], i);
		}
		
		// NBAs of size 1 to 30, int[][0] = count, int[][1] = sum of M2MA dimension
		results = new int[31][2];
		
		System.out.println("0% complete");
		
		run(10, 5, 10, 2, true);
		run(10, 5, 15, 2, true);
		run(10, 5, 20, 2, true);
		run(10, 5, 25, 2, true);
		run(10, 5, 30, 2, true);
		
		System.out.println("10% complete");
		
		run(10, 5, 10, 3, true);
		run(10, 5, 15, 3, true);
		run(10, 5, 20, 3, true);
		run(10, 5, 25, 3, true);
		run(10, 5, 30, 3, true);
		
		System.out.println("20% complete");
		
		run(10, 10, 20, 3, true);
		run(10, 10, 25, 3, true);
		run(10, 10, 30, 3, true);
		run(10, 10, 35, 3, true);
		run(10, 10, 40, 3, true);
		run(10, 10, 45, 3, true);
		run(10, 10, 50, 3, true);
		
		System.out.println("34% complete");
		
		run(10, 10, 20, 4, true);
		run(10, 10, 25, 4, true);
		run(10, 10, 30, 4, true);
		run(10, 10, 35, 4, true);
		run(10, 10, 40, 4, true);
		run(10, 10, 45, 4, true);
		run(10, 10, 50, 4, true);
		
		System.out.println("48% complete");
		
		run(10, 15, 30, 4, true);
		run(10, 15, 40, 4, true);
		run(10, 15, 50, 4, true);
		run(10, 15, 60, 4, true);
		run(10, 15, 70, 4, true);
		run(10, 15, 80, 4, true);
		
		System.out.println("60% complete");
		
		run(10, 15, 30, 5, true);
		run(10, 15, 40, 5, true);
		run(10, 15, 50, 5, true);
		run(10, 15, 60, 5, true);
		run(10, 15, 70, 5, true);
		run(10, 15, 80, 5, true);
		
		System.out.println("72% complete");
		
		run(10, 20, 40, 5, true);
		run(10, 20, 50, 5, true);
		run(10, 20, 60, 5, true);
		run(10, 20, 70, 5, true);
		run(10, 20, 80, 5, true);
		run(10, 20, 90, 5, true);
		run(10, 20, 100, 5, true);
		
		System.out.println("86% complete");
		
		run(10, 20, 40, 6, true);
		run(10, 20, 50, 6, true);
		run(10, 20, 60, 6, true);
		run(10, 20, 70, 6, true);
		run(10, 20, 80, 6, true);
		run(10, 20, 90, 6, true);
		run(10, 20, 100, 6, true);
		
		System.out.println("100% complete\n");
		
		System.out.println("NBA to M2MA results.");
		for (int i = 1; i <= 30; i++) {
			if (results[i][0] != 0) {
				System.out.println(i + ": " + results[i][0] + " " + results[i][1] + " " + ((double) results[i][1])/results[i][0]);
			}
		}
	}
	
	@SuppressWarnings("unchecked")
	public static void run(int numTimes, int numStates, int numTransitions, int numFinalStates, boolean inNBA) throws Exception {
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
				NBAtransitions(numStates, numTransitions, tempTransitions, reverseTempTransitions);
			} else {
				tempDBAtoM2MA.DBAtransitions(numStates, numTransitions, tempTransitions, reverseTempTransitions);
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
				
				NBA.NBATransitions = new ArrayList[NBA.NBAStates + 1][Mod2_MA.alphabet.length - 1];
				for (int j = 1; j <= NBA.NBAStates; j++) {
					for (int k = 0; k < Mod2_MA.alphabet.length - 1; k++) {
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
				
				Mod2_MA.Hankel = null;
				Mod2_MA.learn();
				
				results[NBA.NBAStates][0]++;
				results[NBA.NBAStates][1] += Mod2_MA.learnedSize;
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
	
	public static void NBAtransitions(int numStates, int numTransitions, ArrayList<int[]>[] tempTransitions, ArrayList<int[]>[] reverseTempTransitions) {
		boolean[][][] usedTransitions = new boolean[numStates + 1][Mod2_MA.alphabet.length - 1][numStates + 1];
		
		for (int i = 0; i < numTransitions; i++) {
			int state1 = (int) (Math.random() * numStates) + 1;
			int letter = (int) (Math.random() * (Mod2_MA.alphabet.length - 1));
			int state2 = (int) (Math.random() * numStates) + 1;
			
			if (!usedTransitions[state1][letter][state2]) {
				int[] transition = new int[2];
				transition[0] = letter;
				transition[1] = state2;
				tempTransitions[state1].add(transition);
				
				int[] reverseTransition = new int[2];
				reverseTransition[0] = letter;
				reverseTransition[1] = state1;
				reverseTempTransitions[state2].add(reverseTransition);
				
				usedTransitions[state1][letter][state2] = true;
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
}
