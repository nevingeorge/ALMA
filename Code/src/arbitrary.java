/*
 * Author: Nevin George
 * Advisor: Dana Angluin
 * Program Description: The program displays to stdout the M2MA learned using a membership query function
 * specified in MQ.java and statistical equivalence queries. The program can be used to approximately learn any
 * type of automata, provided that MQ.java contains the desired automata's membership query function.
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
		System.out.println("Program Description:");
		System.out.println("The program displays to stdout the M2MA learned using a membership query function"
				+ " specified in MQ.java and statistical equivalence queries.");
		System.out.println("The program can be used to approximately learn any"
				+ " type of automata, provided that MQ.java contains the desired automata's membership query function.");
		System.out.println("\n----------------------------\n");
		
		readInput();
		
		M2MA.learn();
		
		M2MA.displayResults();
		
		M2MA.displayRuntime();
		
		M2MA.operationsOnLearnedMA();
	}
	
	public static void readInput() throws Exception {
		System.out.println("Input file name and optional flag -v (e.g. arb_input1.txt or arb_input1.txt -v)");

		BufferedReader f = M2MA.getFile(true, false, false, false);
		
		// membership query function to call in MQ.java
		try {
			MQMethod = (new MQ()).getClass().getMethod(M2MA.readFile(f), String.class);
		} catch (Exception e) {
			M2MA.throwException(f, "Invalid input: invalid membership query function name.");
		}
		
		EQMaxTestLen = Integer.parseInt(M2MA.readFile(f));
		EQNumTests = Integer.parseInt(M2MA.readFile(f));
		EQLimit = Integer.parseInt(M2MA.readFile(f));
		EQNumPerformed = 0;
		M2MA.readAlphabet(f, false);
				
		f.close();
	}
	
	// approximate equivalence query using a given number of randomly generated tests of some max length
	public static boolean EQstatistical(HashMap<Integer, ArrayList<Integer>> hypothesisFinalVector, HashMap<Integer, ArrayList<Integer>>[] hypothesisTransitionMatrices) throws Exception {
		int numFail = 0;
		for (int i=0; i<EQNumTests; i++) {
			String test = M2MA.genTest((int) (Math.random() * (EQMaxTestLen + 1)), false);
			
			if (M2MA.MQ(test) != M2MA.MQArbitrary(hypothesisFinalVector, hypothesisTransitionMatrices, test)) {
				// found a counter-example
				// count the number of counter-examples
				if (EQNumPerformed == EQLimit-1) {
					numFail++;
				} else {
					EQNumPerformed++;
					M2MA.counterExample = test;
					return false;
				}
			}
		}
		
		// performs EQlimit equivalence queries
		if (EQNumPerformed == EQLimit-1 && numFail != 0) {
			M2MA.resultFinalVector = hypothesisFinalVector;
			M2MA.resultTransitionMatrices = hypothesisTransitionMatrices;
			M2MA.displayResults();
			
			System.out.println("Reached equivalence query limit.\nFinal equivalence query failed on " + numFail + " out of " + EQNumTests + " tests.");
			System.exit(0);
		}

		return true;
	}
}
