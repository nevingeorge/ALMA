import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class M2MA_experiments {
	
	public static int alphabetSize;
	public static int minDim;
	public static int maxDim;
	public static int dimIncrement;
	public static int numM2MAs;

	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws Exception {
		System.out.println("Program Description:");
		System.out.println("The program will run M2MA.java on a set of randomly generated M2MAs of dimensions within a certain range.");
		System.out.println("For example, if the range given is [10, 40] and the increment is 10, the program will run M2MA.java on M2MAs of dimension 10, 20, 30, and 40.");
		System.out.println("\n----------------------------\n");
		
		readInput();
		
		// initialize M2MA parameters
		M2MA.observationTableFlag = false;
		M2MA.minProgressFlag = false;
		M2MA.minDimensionFlag = false;
		M2MA.dfaFlag = false;
		M2MA.displayFlag = false;
		
		M2MA.alphabet = new String[alphabetSize];
		for (int i = 0; i < alphabetSize; i++) {
			M2MA.alphabet[i] = Integer.toString(i);
		}
		M2MA.letterToIndex = new HashMap<String, Integer>();
		for (int i = 0; i < alphabetSize; i++) {
			M2MA.letterToIndex.put(M2MA.alphabet[i], i);
		}
		
		for (int dim = minDim; dim <= maxDim; dim += dimIncrement) {
			M2MA.inputSize = dim;
			double sumRuntimes = 0;
			
			for (int i = 0; i < numM2MAs; i++) {
				// initialize M2MA parameters
				M2MA.inputFinalVector = M2MA.initialize(1, M2MA.inputSize);
				for (int j = 1; j <= M2MA.inputSize; j++) {
					if (Math.random() < .5) {
						M2MA.addElement(M2MA.inputFinalVector, 1, j);
					}
				}
				
				M2MA.inputTransitionMatrices = new HashMap[alphabetSize];
				for (int j = 0; j < alphabetSize; j++) {
					HashMap<Integer, ArrayList<Integer>> transitionMatrix = M2MA.initialize(M2MA.inputSize, M2MA.inputSize);
					
					for (int k = 1; k <= M2MA.inputSize; k++) {
						for (int l = 1; l <= M2MA.inputSize; l++) {
							if (Math.random() < .5) {
								M2MA.addElement(transitionMatrix, k, l);
							}
						}
					}
					
					M2MA.inputTransitionMatrices[j] = transitionMatrix;
				}
				
				long startTime = System.nanoTime();
				
				// run the M2MA.java algorithms
				M2MA.minimize();
				
				M2MA.learn();
				
				if (M2MA.minSize != M2MA.learnedSize) {
					M2MA.throwException(null, "Algorithm failed: the learned mod-2-MA has a different dimension "
							+ "(" + M2MA.learnedSize + ") than the minimized mod-2-MA (" + M2MA.minSize + ").");
				}

				if (M2MA.finalCheck(25, 1000, false)) {
					long endTime = System.nanoTime();
					sumRuntimes += endTime - startTime;
				} else {
					M2MA.throwException(null, "Algorithm failed: failed final check.");
				}
			}
			
			double averageTime = sumRuntimes / numM2MAs;
			double convertedTime = averageTime / Math.pow(10,  9);
			double roundedTime = ((int) (convertedTime * 100))/100.0;
			System.out.println("Average run time of M2MAs of dimension " + dim + ": " + roundedTime);
		}
	}
	
	public static void readInput() {
		Scanner scan = new Scanner(System.in);
		
		System.out.println("Enter the size of the alphabet.");
		alphabetSize = Integer.parseInt(scan.nextLine());
		
		System.out.println("Enter the minimum dimension in the range.");
		minDim = Integer.parseInt(scan.nextLine());
		
		System.out.println("Enter the maximum dimension in the range.");
		maxDim = Integer.parseInt(scan.nextLine());
		
		System.out.println("Enter the increment.");
		dimIncrement = Integer.parseInt(scan.nextLine());
		
		System.out.println("Enter the number of M2MAs to learn for each dimension.");
		numM2MAs = Integer.parseInt(scan.nextLine());
		scan.close();
		
		System.out.println("\n----------------------------\n");
		System.out.println("Running experiments:");
	}

}
