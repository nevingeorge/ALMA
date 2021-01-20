/*
 * Author: Nevin George
 * Advisor: Dana Angluin
 * Program Description: The program generates statistics such as mean and median for the results of convert.java.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

public class statistics {
	
	public static boolean inDFA;
	
	public static int[][] meanResultsM2MA;
	public static int[][] meanResultsDFA;
	
	public static ArrayList<Integer>[] medianResultsM2MA;
	public static ArrayList<Integer>[] medianResultsDFA;

	public static void main(String[] args) throws Exception {
		Scanner in = new Scanner(System.in);
		
		System.out.println("Input number of files to read.");
		int numFiles = Integer.parseInt(in.nextLine());
		
		String[] fileNames = new String[numFiles];
		for (int i = 0; i < numFiles; i++) {
			System.out.println("Input file name.");
			fileNames[i] = in.nextLine();
		}
		
		inDFA = false;
		
		BufferedReader firstFile = new BufferedReader(new FileReader(fileNames[0]));
		if (firstFile.readLine().contains("DFA")) {
			inDFA = true;
		}
		firstFile.close();
		
		System.out.println("Enter 1 to give each file equal weighting or 2 to combine "
				+ "the input from all the text files into one larger text file before calculating the statistics.");
		int equalWeighting = Integer.parseInt(in.nextLine());
	
		in.close();
		
		if (equalWeighting == 1) {
			// int[][][0] = mean results, int[][][1] = median results
			int[][][] fileResultsM2MA = new int[numFiles][31][2];
			int[][][] fileResultsDFA = new int[numFiles][31][2];
			
			for (int i = 0; i < numFiles; i++) {
				initialize();
				
				readFile(fileNames[i]);
				
				for (int j = 1; j <= 30; j++) {
					if (meanResultsM2MA[j][0] > 0) {
						fileResultsM2MA[i][j][0] = calculateMean(meanResultsM2MA[j]);
						fileResultsM2MA[i][j][1] = calculateMedian(medianResultsM2MA[j]);
					}
					
					if (inDFA && meanResultsDFA[j][0] > 0) {
						fileResultsDFA[i][j][0] = calculateMean(meanResultsDFA[j]);
						fileResultsDFA[i][j][1] = calculateMedian(medianResultsDFA[j]);
					}
				}
			}
			
			System.out.println("\nExperiment Results");
			System.out.println("------------------");
			
			displayResultsEqualWeighting(numFiles, fileResultsM2MA, "M2MA");
			
			if (inDFA) {
				displayResultsEqualWeighting(numFiles, fileResultsDFA, "DFA");
			}
		} else {
			initialize();
			
			for (int i = 0; i < numFiles; i++) {
				readFile(fileNames[i]);
			}
			
			System.out.println("\nExperiment Results");
			System.out.println("------------------");
			
			displayResultsCombined(meanResultsM2MA, medianResultsM2MA, "M2MA");
			
			if (inDFA) {
				displayResultsCombined(meanResultsDFA, medianResultsDFA, "DFA");
			}
		}
	}
	
	@SuppressWarnings("unchecked")
	public static void initialize() {
		meanResultsM2MA = new int[31][2];
		meanResultsDFA = new int[31][2];
		
		medianResultsM2MA = new ArrayList[31];
		for (int i = 1; i <= 30; i++) {
			medianResultsM2MA[i] = new ArrayList<Integer>();
		}
		
		medianResultsDFA = new ArrayList[31];
		for (int i = 1; i <= 30; i++) {
			medianResultsDFA[i] = new ArrayList<Integer>();
		}
	}
	
	public static void readFile(String fileName) throws IOException {
		BufferedReader f = new BufferedReader(new FileReader(fileName));
		f.readLine();
		String line = f.readLine();
		
		if (line != null) {		
			String[] lineArr = line.split(" ");
			
			while (true) {
				int initialSize = Integer.parseInt(lineArr[0]);
						
				int convertedM2MASize = Integer.parseInt(lineArr[1]);
				meanResultsM2MA[initialSize][0]++;
				meanResultsM2MA[initialSize][1] += convertedM2MASize;
				medianResultsM2MA[initialSize].add(convertedM2MASize);
				
				if (inDFA) {
					int convertedDFASize = Integer.parseInt(lineArr[2]);
					meanResultsDFA[initialSize][0]++;
					meanResultsDFA[initialSize][1] += convertedDFASize;
					medianResultsDFA[initialSize].add(convertedDFASize);
				}
				
				line = f.readLine();
				if (line == null) {
					break;
				}
				lineArr = line.split(" ");
			}
		}
		
		f.close();
	}
	
	public static void displayResultsEqualWeighting(int numFiles, int[][][] fileResults, String convertedAutomataName) {
		System.out.println("\nInitial automata size: mean size of converted " + convertedAutomataName + ", median size of converted " + convertedAutomataName);
		for (int j = 1; j <= 30; j++) {
			int countMean = 0;
			int sumMean = 0;
			int countMedian = 0;
			int sumMedian = 0;
			
			for (int i = 0; i < numFiles; i++) {
				if (fileResults[i][j][0] > 0) {
					countMean++;
					sumMean += fileResults[i][j][0];
				}
				
				if (fileResults[i][j][1] > 0) {
					countMedian++;
					sumMedian += fileResults[i][j][1];
				}
			}
			
			if (countMean > 0) {
				int mean = (int) Math.round(sumMean / (double) countMean);
				int median = (int) Math.round(sumMedian / (double) countMedian);
				System.out.println(j + ": " + mean + ", " + median);
			}
		}
	}
	
	public static void displayResultsCombined(int[][] meanResults, ArrayList<Integer>[] medianResults, String convertedAutomataName) {
		System.out.println("\nInitial automata size: mean size of converted " + convertedAutomataName + ", median size of converted " + convertedAutomataName);
		for (int i = 1; i <= 30; i++) {
			if (meanResults[i][0] > 0) {
				int mean = calculateMean(meanResults[i]);
				int median = calculateMedian(medianResults[i]);
				
				System.out.println(i + ": " + mean + ", " + median);
			}
		}
	}
	
	
	public static int calculateMean(int[] arr) {
		return (int) Math.round(arr[1] / (double) arr[0]);
	}

	public static int calculateMedian(ArrayList<Integer> arr) {
		Collections.sort(arr);
		
		int len = arr.size();
		
		if (len % 2 == 0) {
			return (int) Math.round((arr.get((len / 2) - 1) + arr.get(len / 2)) / 2.0);
		} else {
			return arr.get(len / 2);
		}
	}

}
