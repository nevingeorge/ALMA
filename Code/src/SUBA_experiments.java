import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class SUBA_experiments {

	public static void main(String[] args) throws FileNotFoundException, IOException, Exception {
		System.out.println("Running SUBA.jar on the following input files:");
		System.out.println("SUBA_exp1.txt: language = a Σ* (Σ* b Σ*)^ω, input SUBA size = 2, expected run time = ~.22s");
		System.out.println("SUBA_exp2.txt: language = Σ* a Σ^5 a b^ω, input SUBA size = 8, expected run time = ~.51s");
		System.out.println("SUBA_exp3.txt: language = ((a+b)* (a (a+b) a (a+b) c + b (a+b) b (a+b) d))^ω, input SUBA size = 9, expected run time = ~57.87s");
		System.out.println("SUBA_exp4.txt: language = (a^* a^4 b)^w, input SUBA size = 5, expected run time = ~1.29s");
		System.out.println("SUBA_exp5.txt: language = (a^* a^5 b)^w, input SUBA size = 6, expected run time = ~3.17s");
		System.out.println("SUBA_exp6.txt: language = (a^* a^6 b)^w, input SUBA size = 7, expected run time = ~9.47s");
		System.out.println("SUBA_exp7.txt: language = a^ω, input SUBA size = 1, expected run time = ~0.16s");
		System.out.println("SUBA_exp8.txt: language = (ab^5)^ω, input SUBA size = 6, expected run time = ~1.21s");
		System.out.println("SUBA_exp9.txt: language = (ab^10)^ω, input SUBA size = 11, expected run time = ~31.59s");
		System.out.println("SUBA_exp10.txt: language = (ab^15)^ω, input SUBA size = 16, expected run time = ~470.74s");
		System.out.println("SUBA_exp11.txt: language = (ab^20)^ω, input SUBA size = 21, expected run time = ~3064.32s");
		System.out.println();
		
		M2MA.observationTableFlag = false;
		M2MA.minProgressFlag = false;
		M2MA.minDimensionFlag = false;
		M2MA.dfaFlag = false;
		
		for (int i = 1; i <= 11; i++) {
			String fileName = "SUBA_exp" + Integer.toString(i) + ".txt";
			System.out.println("Running " + fileName + ": ");
			
			M2MA.startTime = System.nanoTime();
			
			SUBAtoUFAwithFile(new BufferedReader(new FileReader(fileName)));
			
			SUBA.UFAtoMod2MA();

			M2MA.minimize();
			
			M2MA.learn();
			
			if (M2MA.minSize != M2MA.learnedSize) {
				M2MA.throwException(null, "Algorithm failed: the learned mod-2-MA has a different dimension "
						+ "(" + M2MA.learnedSize + ") than the minimized mod-2-MA (" + M2MA.minSize + ").");
			}
			
			if (SUBA.finalCheck(25,1000)) {
				System.out.println("Learned M2MA dimension: " + M2MA.learnedSize);
				
				M2MA.displayRuntime();
			} else {
				M2MA.throwException(null, "Failed final check.");
			}
		}
		
		System.out.println("\nProgram complete.");
	}
	
	@SuppressWarnings("unchecked")
	public static void SUBAtoUFAwithFile(BufferedReader f) throws IOException, Exception {
		// UFAStates = SUBAStates U (SUBAStates x SUBAStates x {0,1})
		SUBA.SUBAStates = Integer.parseInt(M2MA.readFile(f));
		SUBA.UFAStates = SUBA.SUBAStates + SUBA.SUBAStates * SUBA.SUBAStates * 2;
		
		// alphabet ΣU{$}
		M2MA.readAlphabet(f, true);
		
		StringTokenizer st = new StringTokenizer(M2MA.readFile(f));
		SUBA.SUBAFinalStates = new boolean[SUBA.SUBAStates+1];
		while (st.hasMoreTokens()) {
			int state = Integer.parseInt(st.nextToken());
			if (1 <= state && state <= SUBA.SUBAStates && !SUBA.SUBAFinalStates[state]) {
				SUBA.SUBAFinalStates[state] = true;
			} else {
				M2MA.throwException(f,"Invalid input: invalid or duplicate final state.");
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
		int numTransitions = Integer.parseInt(M2MA.readFile(f));
		if (numTransitions > ((M2MA.alphabet.length - 1) * SUBA.SUBAStates * SUBA.SUBAStates)) {
			M2MA.throwException(f, "Invalid input: invalid number of transitions.");
		}
		
		/*
		 *  For each index (q,a) where q∈SUBAStates and a∈alphabet, transition_SUBA[q][a] is an ArrayList containing 
		 *  all of the reachable states from (q,a).
		 *  The alphabet for the SUBA does not include $.
		 */
		SUBA.SUBATransitions = new ArrayList[SUBA.SUBAStates + 1][M2MA.alphabet.length - 1];
		for (int i=1; i<=SUBA.SUBAStates; i++) {
			for (int j=0; j<M2MA.alphabet.length-1; j++) {
				SUBA.SUBATransitions[i][j] = new ArrayList<Integer>();
			}
		}
		
		// (start state, letter, end state)
		SUBA.UFATransitions = new boolean[SUBA.UFAStates + 1][M2MA.alphabet.length][SUBA.UFAStates + 1];
		
		// lines of the form q_j a q_k, where q_j,q_k∈SUBAStates and a∈alphabet
		for (int i=0; i<numTransitions; i++) {
			st = new StringTokenizer(M2MA.readFile(f));
			int p_start = Integer.parseInt(st.nextToken());
			
			String letter = st.nextToken();
			
			int a = M2MA.letterToIndex.get(letter);
			int p_end = Integer.parseInt(st.nextToken());
			if (p_start < 1 || p_start > SUBA.SUBAStates || p_end < 1 || p_end > SUBA.SUBAStates) {
				M2MA.throwException(f, "Invalid input: invalid transition.");
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
			SUBA.UFATransitions[q][M2MA.letterToIndex.get("$")][SUBA.getIndex(q, q, 0)] = true;
			SUBA.UFAFinalStates[SUBA.getIndex(q, q, 1)] = true;
		}
		
		if (convert.f == null) {
			if (M2MA.readFile(f) != null) {
				M2MA.throwException(f, "Invalid input: more transitions inputted than specified.");
			}

			f.close();
		}
	}

}
