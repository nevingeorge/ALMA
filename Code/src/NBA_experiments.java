import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class NBA_experiments {

	public static void main(String[] args) throws Exception {
		System.out.println("Running NBA.jar on the following input files:");
		System.out.println("NBA_exp1.txt: input NBA size = 2, number of tests/EQ = 10000, expected run time = ~1.02s");
		System.out.println("NBA_exp2.txt: input NBA size = 4, number of tests/EQ = 10000, expected run time = ~.41s");
		System.out.println("NBA_exp3.txt: input NBA size = 6, number of tests/EQ = 10000, expected run time = ~1.80s");
		System.out.println("NBA_exp4.txt: input NBA size = 8, number of tests/EQ = 100000, expected run time = ~14.00s");
		System.out.println("NBA_exp5.txt: input NBA size = 10, number of tests/EQ = 100000, expected run time = ~1833.65s");
		System.out.println();
		
		M2MA.observationTableFlag = false;
		M2MA.minProgressFlag = false;
		M2MA.minDimensionFlag = false;
		M2MA.dfaFlag = false;
		
		for (int i = 1; i <= 5; i++) {
			String fileName = "NBA_exp" + Integer.toString(i) + ".txt";
			System.out.println("Running " + fileName + ": ");
			
			M2MA.startTime = System.nanoTime();
			
			NBAReadInputWithFile(new BufferedReader(new FileReader(fileName)));
			
			M2MA.Hankel = null;
			
			M2MA.learn();
			
			System.out.println("Learned M2MA dimension: " + M2MA.learnedSize);
			
			M2MA.displayRuntime();
		}
		
		System.out.println("\nProgram complete.");
	}
	
	@SuppressWarnings("unchecked")
	public static void NBAReadInputWithFile(BufferedReader f) throws Exception {
		arbitrary.EQMaxTestLen = Integer.parseInt(M2MA.readFile(f));	
		arbitrary.EQNumTests = Integer.parseInt(M2MA.readFile(f));
		arbitrary.EQLimit = Integer.parseInt(M2MA.readFile(f));
		arbitrary.EQNumPerformed = 0;

		NBA.NBAStates = Integer.parseInt(M2MA.readFile(f));
		
		// alphabet ΣU{$}
		M2MA.readAlphabet(f, true);
		
		StringTokenizer st = new StringTokenizer(M2MA.readFile(f));
		NBA.NBAFinalStates = new boolean[NBA.NBAStates+1];
		while (st.hasMoreTokens()) {
			int state = Integer.parseInt(st.nextToken());
			if(1<=state && state<=NBA.NBAStates && !NBA.NBAFinalStates[state]) {
				NBA.NBAFinalStates[state] = true;
			} else {
				M2MA.throwException(f, "Invalid input: invalid or duplicate final state.");
			}
		}
		
		int numTransitions = Integer.parseInt(M2MA.readFile(f));
		if(numTransitions > ((M2MA.alphabet.length - 1) * NBA.NBAStates * NBA.NBAStates)) {
			M2MA.throwException(f, "Invalid input: invalid number of transitions.");
		}
		
		// (start state, letter, end state)
		NBA.NBATransitions = new ArrayList[NBA.NBAStates+1][M2MA.alphabet.length-1];
		for (int i=1; i<=NBA.NBAStates; i++) {
			for (int j=0; j<M2MA.alphabet.length-1; j++) {
				NBA.NBATransitions[i][j] = new ArrayList<Integer>();
			}
		}
		
		// lines of the form q_i a q_j, where q_i,q_j∈NBAStates and a∈alphabet
		for (int i=0; i<numTransitions; i++) {
			st = new StringTokenizer(M2MA.readFile(f));
			int p_start = Integer.parseInt(st.nextToken());
			
			String letter = st.nextToken();
			
			int a = M2MA.letterToIndex.get(letter);
			int p_end = Integer.parseInt(st.nextToken());
			if (p_start < 1 || p_start > NBA.NBAStates || p_end < 1 || p_end > NBA.NBAStates) {
				M2MA.throwException(f, "Invalid input: invalid transition.");
			}

			NBA.NBATransitions[p_start][a].add(p_end);
		}
		
		if (M2MA.readFile(f) != null) {
			M2MA.throwException(f, "Invalid input: more transitions inputted than specified.");
		}
		
		f.close();
	}
}
