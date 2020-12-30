import java.util.ArrayList;
import java.util.HashMap;

public class tempDBAtoM2MA {

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
		tempNBAtoM2MA.results = new int[31][2];

		System.out.println("0% complete");

		tempNBAtoM2MA.run(10, 5, 10, 2, false);
		tempNBAtoM2MA.run(10, 5, 11, 2, false);
		tempNBAtoM2MA.run(10, 5, 12, 2, false);
		tempNBAtoM2MA.run(10, 5, 13, 2, false);
		tempNBAtoM2MA.run(10, 5, 14, 2, false);
		tempNBAtoM2MA.run(10, 5, 15, 2, false);
		
		System.out.println("12% complete");
		
		tempNBAtoM2MA.run(10, 5, 10, 3, false);
		tempNBAtoM2MA.run(10, 5, 11, 3, false);
		tempNBAtoM2MA.run(10, 5, 12, 3, false);
		tempNBAtoM2MA.run(10, 5, 13, 3, false);
		tempNBAtoM2MA.run(10, 5, 14, 3, false);
		tempNBAtoM2MA.run(10, 5, 15, 3, false);
		
		System.out.println("24% complete");
		
		tempNBAtoM2MA.run(10, 10, 20, 3, false);
		tempNBAtoM2MA.run(10, 10, 22, 3, false);
		tempNBAtoM2MA.run(10, 10, 24, 3, false);
		tempNBAtoM2MA.run(10, 10, 26, 3, false);
		tempNBAtoM2MA.run(10, 10, 28, 3, false);
		tempNBAtoM2MA.run(10, 10, 30, 3, false);
		
		System.out.println("36% complete");
		
		tempNBAtoM2MA.run(10, 10, 20, 4, false);
		tempNBAtoM2MA.run(10, 10, 22, 4, false);
		tempNBAtoM2MA.run(10, 10, 24, 4, false);
		tempNBAtoM2MA.run(10, 10, 26, 4, false);
		tempNBAtoM2MA.run(10, 10, 28, 4, false);
		tempNBAtoM2MA.run(10, 10, 30, 4, false);
		
		System.out.println("48% complete");
		
		tempNBAtoM2MA.run(10, 15, 30, 4, false);
		tempNBAtoM2MA.run(10, 15, 33, 4, false);
		tempNBAtoM2MA.run(10, 15, 36, 4, false);
		tempNBAtoM2MA.run(10, 15, 39, 4, false);
		tempNBAtoM2MA.run(10, 15, 42, 4, false);
		tempNBAtoM2MA.run(10, 15, 45, 4, false);
		
		System.out.println("60% complete");
		
		tempNBAtoM2MA.run(10, 15, 30, 5, false);
		tempNBAtoM2MA.run(10, 15, 33, 5, false);
		tempNBAtoM2MA.run(10, 15, 36, 5, false);
		tempNBAtoM2MA.run(10, 15, 39, 5, false);
		tempNBAtoM2MA.run(10, 15, 42, 5, false);
		tempNBAtoM2MA.run(10, 15, 45, 5, false);
		
		System.out.println("72% complete");
	
		tempNBAtoM2MA.run(10, 20, 42, 5, false);
		tempNBAtoM2MA.run(10, 20, 45, 5, false);
		tempNBAtoM2MA.run(10, 20, 48, 5, false);
		tempNBAtoM2MA.run(10, 20, 51, 5, false);
		tempNBAtoM2MA.run(10, 20, 54, 5, false);
		tempNBAtoM2MA.run(10, 20, 57, 5, false);
		tempNBAtoM2MA.run(10, 20, 60, 5, false);
		
		System.out.println("86% complete");
		
		tempNBAtoM2MA.run(10, 20, 42, 6, false);
		tempNBAtoM2MA.run(10, 20, 45, 6, false);
		tempNBAtoM2MA.run(10, 20, 48, 6, false);
		tempNBAtoM2MA.run(10, 20, 51, 6, false);
		tempNBAtoM2MA.run(10, 20, 54, 6, false);
		tempNBAtoM2MA.run(10, 20, 57, 6, false);
		tempNBAtoM2MA.run(10, 20, 60, 6, false);
		
		System.out.println("100% complete\n");
		
		System.out.println("DBA to M2MA results.");
		for (int i = 1; i <= 30; i++) {
			if (tempNBAtoM2MA.results[i][0] != 0) {
				System.out.println(i + ": " + tempNBAtoM2MA.results[i][0] + " " + tempNBAtoM2MA.results[i][1] + " " + ((double) tempNBAtoM2MA.results[i][1])/tempNBAtoM2MA.results[i][0]);
			}
		}
	}
	
	public static void DBAtransitions(int numStates, int numTransitions, ArrayList<int[]>[] tempTransitions, ArrayList<int[]>[] reverseTempTransitions) {
		boolean[][] usedTransitions = new boolean[numStates + 1][Mod2_MA.alphabet.length - 1];
		
		for (int i = 0; i < numTransitions; i++) {
			int state1 = (int) (Math.random() * numStates) + 1;
			int letter = (int) (Math.random() * (Mod2_MA.alphabet.length - 1));
			int state2 = (int) (Math.random() * numStates) + 1;
			
			if (!usedTransitions[state1][letter]) {
				int[] transition = new int[2];
				transition[0] = letter;
				transition[1] = state2;
				tempTransitions[state1].add(transition);
				
				int[] reverseTransition = new int[2];
				reverseTransition[0] = letter;
				reverseTransition[1] = state1;
				reverseTempTransitions[state2].add(reverseTransition);
				
				usedTransitions[state1][letter] = true;
			}
		}
	}
}