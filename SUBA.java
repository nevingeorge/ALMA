/*
 * Author: Nevin George
 * Advisor: Dana Angluin
 * Program Description: The algorithm takes in as input a SUBA of the form (Q, Σ, Q_in, ∆, F) and converts it
 * into an equivalent UFA of the form (Q',ΣU{$},Q_in,∆',F'). The resulting UFA is then converted into an
 * equivalent mod-2-MA and learned using the Mod2_MA.java program.
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
 * 5 N. Bousquet and C. Löding. Equivalence and inclusion problem for strongly unambiguous büchi automata. In 
 *   Language and Automata Theory and Applications, 4th International Conference, LATA 2010, Trier, Germany, May 
 *   24-28, 2010. Proceedings, pages 118–129, 2010. 
 *   URL: https: //doi.org/10.1007/978-3-642-13089-2_10, doi:10.1007/978-3-642-13089-2\_10.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.StringTokenizer;

public class SUBA {
	// UFA = (Q',ΣU{$},Q_in,∆',F')
	
	// Q' has states q_1 to q_numStates
	public static int numStates;
	
	// alphabet ΣU{$}
	public static int alphabetSize;
	public static Character[] alphabet;
	
	// initial states Q_in
	public static boolean[] Q_in;
	
	// transition matrix Δ'
	public static boolean[][][] transition;
	
	// final states F'
	public static boolean[] F;

	public static void main(String[] args) throws Exception {
		initializeUFA();
		System.out.println("done");
	}
	
	public static void initializeUFA() throws Exception {
		/* The input file containing the SUBA (Q, Σ, Q_in, ∆, F) must have the following format (no line 
		 * separation, characters are space separated, and lines beginning with // are ignored):
		 * <number of states (Q)>
		 * <alphabet size>
		 * <characters in the alphabet>
		 * <initial states (Q_in)>
		 * <final states (F)>
		 * <number of transitions>
		 * The remaining lines are the transitions, with each line having the form q_j a q_k, where q_j,q_k∈Q
		 * and a is a letter in the alphabet.
		 * 
		 * Example input files can be found in the GitHub repository.
		 */
		
		// converts the input SUBA into the UFA described in Bousquet and Löding of the form (Q',ΣU{$},Q_in,∆',F')
		
		// reads from the input file in quotations, ***EDIT THE FILE NAME DEPENDING ON THE INTENDED INPUT FILE***
		BufferedReader f = new BufferedReader(new FileReader("SUBA_input1.txt"));
		
		// Q' = Q U (Q x Q x {0,1})
		String line = f.readLine();
		while(line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();
		
		int Q = Integer.parseInt(line);
		numStates = Q + Q*Q*2;
		
		
		
		// alphabet size
		line = f.readLine();
		while(line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();
		
		// add {$} to the language
		alphabetSize = Integer.parseInt(line) + 1;
		
		// alphabet ΣU{$}
		line = f.readLine();
		while(line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();
		StringTokenizer st = new StringTokenizer(line);
		
		alphabet = new Character[alphabetSize];
		for(int i=0;i<alphabetSize-1;i++) {
			String letter = st.nextToken();
			if(letter.length()!=1 || letter.charAt(0) == '$') {
				f.close();
				throw new Exception("Invalid input: invalid character in the alphabet");
			}
			alphabet[i] = letter.charAt(0);
		}
		alphabet[alphabetSize-1] = '$';
		
		if(st.hasMoreTokens()) {
			f.close();
			throw new Exception("Invalid input: alphabet size exceeds the specified size");
		}
		
		// maps each letter in alphabet to an index
		HashMap<Character, Integer> letterToIndex = new HashMap<Character, Integer>();
		for(int i=0;i<alphabetSize;i++)
			letterToIndex.put(alphabet[i], i);
		
		
		
		// initial states (same initial states as input SUBA)
		line = f.readLine();
		while(line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();
		st = new StringTokenizer(line);
		
		Q_in = new boolean[numStates+1];
		while(st.hasMoreTokens()) {
			int state = Integer.parseInt(st.nextToken());
			if(1<=state && state<=Q && !Q_in[state])
				Q_in[state] = true;
			else {
				f.close();
				throw new Exception("Invalid input: invalid or duplicate input state");
			}
		}
		

		
		// final states for input SUBA
		line = f.readLine();
		while(line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();
		st = new StringTokenizer(line);
		
		boolean[] F_SUBA = new boolean[Q+1];
		while(st.hasMoreTokens()) {
			int state = Integer.parseInt(st.nextToken());
		if(1<=state && state<=Q && !F_SUBA[state])
				F_SUBA[state] = true;
			else {
				f.close();
				throw new Exception("Invalid input: invalid or duplicate final state");
			}
		}
		
		
		
		/* Following the paper by Bousquet and Löding, Δ' contains (where q,p,p'∈Q)
		 * - all transitions from Δ
		 * - all transitions of the form (q,$,(q,q,0))
		 * - all transitions of the form ((q,p,i),a,(q,p',i')), where (p,a,p')∈Δ, and
		 * 	 i' = 1 if p'∈F and i if p'∉F
		 * 
		 * Transitions will be stored in a numStates x alphabetSize x numStates adjacency matrix.
		 * Let the number of states in the SUBA be Q.
		 * The first Q states of Δ' will be the states of the SUBA.
		 * The remaining states will be of the form (q_j,q_k,i), where q_j,q_k∈Q and i∈{0,1}.
		 * State (q_j,q_k,i) will be found at index (2*Q*j)+(2*k)-(Q)+(i-1) of Δ'.
		*/
		
		// number of transitions
		line = f.readLine();
		while(line.charAt(0) == '/' && line.charAt(1) == '/')
			line = f.readLine();

		int numTransitions = Integer.parseInt(line);
		if(numTransitions<1 || numTransitions>((alphabetSize-1)*Q*Q)) {
			f.close();
			throw new Exception("Invalid input: invalid number of transitions");
		}
		
		// transition matrix Δ' = (start state, letter, end state) 
		transition = new boolean[numStates+1][alphabetSize][numStates+1];
		
		// lines of the form q_j a q_k, where q_j,q_k∈Q
		for(int i=0;i<numTransitions;i++) {
			line = f.readLine();
			while(line.charAt(0) == '/' && line.charAt(1) == '/')
				line = f.readLine();
			
			if(line.length() != 5) {
				f.close();
				throw new Exception("Invalid input: invalid transition");
			}
				
			int p_start = Character.getNumericValue(line.charAt(0));
			int a = letterToIndex.get(line.charAt(2));
			int p_end = Character.getNumericValue(line.charAt(4));
			
			if(p_start<1 || p_start>Q || p_end<1 || p_end>Q) {
				f.close();
				throw new Exception("Invalid input: invalid transition");
			}
			
			// transition from Δ
			transition[p_start][a][p_end] = true;
			
			// transitions of the form ((q,p,i),a,(q,p',i'))
			// p'∈F so i'=1
			if(F_SUBA[p_end]) {
				for(int q=1;q<=Q;q++) {
					// ((q,p,0),a,(q,p',1))
					transition[2*Q*q+2*p_start-Q-1][a][2*Q*q+2*p_end-Q] = true;
					// ((q,p,1),a,(q,p',1))
					transition[2*Q*q+2*p_start-Q][a][2*Q*q+2*p_end-Q] = true;
				}
			}
			// p'∉F so i'=i
			else {
				for(int q=1;q<=Q;q++) {
					// ((q,p,0),a,(q,p',0))
					transition[2*Q*q+2*p_start-Q-1][a][2*Q*q+2*p_end-Q-1] = true;
					// ((q,p,1),a,(q,p',1))
					transition[2*Q*q+2*p_start-Q][a][2*Q*q+2*p_end-Q] = true;
				}
			}
		}

		
		
		// transitions of the form (q,$,(q,q,0)), where q∈Q
		// final states for UFA of the form (q,q,1), where q∈Q
		F = new boolean[numStates+1];
		for(int q=1;q<=Q;q++) {
			transition[q][letterToIndex.get('$')][2*Q*q+2*q-Q-1] = true;
			F[2*Q*q+2*q-Q] = true;
		}
		
		line = f.readLine();
		while(line!=null) {
			if(line.charAt(0) != '/' && line.charAt(1) != '/') {
				f.close();
				throw new Exception("Invalid input: more transitions inputted than specified");
			}
			line = f.readLine();
		}
		f.close();
	}
}
