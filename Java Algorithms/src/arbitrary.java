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
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.StringTokenizer;

public class arbitrary {
	
	// specifies which membership query method in MQ.java to call
	public static int inputType;
	
	// maximum length of a test in the statistical equivalence query
	public static int maxTestLen;
	// number of tests the statistical equivalence query will check
	public static int numTests;

	public static void main(String[] args) throws Exception {
		// reads in the input file
		initialize();
		
		// runs the learning algorithm using the specified membership query function in MQ.java
		Mod2_MA.run();
		
		// displays the learned mod-2-MA
		Mod2_MA.displayResults();
		
		// performs desired operations with the learned mod-2-MA
		Mod2_MA.operations(false);
	}
	
	public static void initialize() throws Exception {
		/* The input file must have the following format (no line separation, entries are space separated, 
		 * and lines beginning with // are ignored):
		 * <characters in the alphabet>
		 * <membership query function number (found in MQ.java)>
		 * <max word length and number of tests>
		 * 
		 * Example input files can be found in the GitHub repository.
		 */
		
		// reads in file name + optional flag -v from stdin
		System.out.println("Input file name and optional flag -v (e.g. arb_input1.txt or arb_input1.txt -v)");
		Mod2_MA.in = new Scanner(System.in);
		String[] arrInput = Mod2_MA.in.nextLine().split(" ");
		Mod2_MA.verbose = false;
		if(arrInput.length > 2)
			Mod2_MA.throwException(null,"Invalid input: too many inputs passed");
		if(arrInput.length == 2) {
			if(arrInput[1].equals("-v"))
				Mod2_MA.verbose = true;
			else
				Mod2_MA.throwException(null,"Invalid input: invalid flag");
		}
		BufferedReader f = new BufferedReader(new FileReader(arrInput[0]));
		System.out.println();
		
		// alphabet
		StringTokenizer st = new StringTokenizer(Mod2_MA.readInput(f));
		ArrayList<Character> tempAlphabet = new ArrayList<Character>();
		while(st.hasMoreTokens()) {
			String letter = st.nextToken();
			if(letter.length()!=1)
				Mod2_MA.throwException(f,"Invalid input: invalid character in the alphabet");
			tempAlphabet.add(letter.charAt(0));
		}
		Mod2_MA.alphabet = new Character[tempAlphabet.size()];
		for(int i=0;i<tempAlphabet.size();i++)
			Mod2_MA.alphabet[i] = tempAlphabet.get(i);
		
		// maps each letter in the alphabet to an index
		Mod2_MA.letterToIndex = new HashMap<Character, Integer>();
		for(int i=0;i<Mod2_MA.alphabet.length;i++)
			Mod2_MA.letterToIndex.put(Mod2_MA.alphabet[i], i);
		
		// specifies which membership query function in MQ.java to call
		inputType = Integer.parseInt(Mod2_MA.readInput(f));
		// currently supports non-negative method numbers <= 0
		if(inputType<0 || inputType>0)
			Mod2_MA.throwException(f, "Invalid input: invalid membership query number");
		
		st = new StringTokenizer(Mod2_MA.readInput(f));
		// maximum length of a test in the statistical equivalence query
		maxTestLen = Integer.parseInt(st.nextToken());
		// number of tests the statistical equivalence query will check
		numTests = Integer.parseInt(st.nextToken());
		if(st.hasMoreTokens())
			Mod2_MA.throwException(f, "Invalid input: too many inputs passed");
		
		f.close();
	}
	
	public static boolean EQapprox(double[] hy, double[][][] hu) {
		// creates numTests tests of length at most maxWordLen
		// checks if the hypothesis and target function have the same output
		for(int i=0;i<numTests;i++) {
			// generates a random length for the test from 0-maxWordLen
			int len = (int)(Math.random()*(maxTestLen+1));
			
			// adds len number of random characters in alphabet to test
			String test = "";
			for(int j=0;j<len;j++)
				test += Mod2_MA.alphabet[(int)(Math.random()*Mod2_MA.alphabet.length)];
			
			if(Mod2_MA.MQ(test)!=Mod2_MA.MQH(hy,hu,test)) {
				// found a counter-example
				Mod2_MA.z = test;
				return false;
			}
		}
		return true;
	}
}
