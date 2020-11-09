/*
 * Author: Nevin George
 * Advisor: Dana Angluin
 * Program Description: The program displays to stdout the mod-2-MA learned using a membership query method
 * specified in MQ.java and statistical equivalence queries. The program can be used to approximately learn any
 * type of automata, provided that MQ.java contains the desired automata's membership query function.
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
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;

public class arbitrary {
	
	// membership query function to call in MQ.java
	public static Method MQMethod;
	
	// EQ settings
	public static int EQMaxTestLen;
	public static int EQNumTests;
	public static int EQLimit;
	public static int EQNumPerformed;

	public static void main(String[] args) throws Exception {
		readInput();
		
		Mod2_MA.learn();
		
		Mod2_MA.displayResults();
		
		Mod2_MA.displayRuntime();
		
		Mod2_MA.operationsOnLearnedMA();
	}
	
	public static void readInput() throws Exception {
		System.out.println("Input file name and optional flag -v (e.g. arb_input1.txt or arb_input1.txt -v)");

		BufferedReader f = Mod2_MA.getFile(true, false, false);
		
		// membership query function to call in MQ.java
		try {
			MQMethod = (new MQ()).getClass().getMethod(Mod2_MA.readFile(f), String.class);
		} catch (Exception e) {
			Mod2_MA.throwException(f, "Invalid input: invalid membership query function name.");
		}
		
		EQMaxTestLen = Integer.parseInt(Mod2_MA.readFile(f));
		EQNumTests = Integer.parseInt(Mod2_MA.readFile(f));
		EQLimit = Integer.parseInt(Mod2_MA.readFile(f));
		EQNumPerformed = 0;
		Mod2_MA.readAlphabet(f, false);
				
		f.close();
	}
	
	public static boolean EQstatistical(HashMap<Integer, ArrayList<Integer>> hypothesisFinalVector, HashMap<Integer, ArrayList<Integer>>[] hypothesisTransitionMatrices) throws Exception {
		int numFail = 0;
		for (int i=0; i<EQNumTests; i++) {
			String test = Mod2_MA.genTest((int) (Math.random() * (EQMaxTestLen + 1)), false);
			
			if (Mod2_MA.MQ(test) != Mod2_MA.MQArbitrary(hypothesisFinalVector, hypothesisTransitionMatrices, test)) {
				// found a counter-example
				// count the number of counter-examples
				if (EQNumPerformed == EQLimit-1) {
					numFail++;
				} else {
					EQNumPerformed++;
					Mod2_MA.counterExample = test;
					return false;
				}
			}
		}
		
		// performs EQlimit equivalence queries
		if (EQNumPerformed == EQLimit-1 && numFail != 0) {
			Mod2_MA.resultFinalVector = hypothesisFinalVector;
			Mod2_MA.resultTransitionMatrices = hypothesisTransitionMatrices;
			Mod2_MA.displayResults();
			
			System.out.println("Reached equivalence query limit.\nFinal equivalence query failed on " + numFail + " out of " + EQNumTests + " tests.");
			System.exit(0);
		}

		return true;
	}
}
