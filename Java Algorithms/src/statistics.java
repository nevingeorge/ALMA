/*
 * Author: Nevin George
 * Advisor: Dana Angluin
 * Program Description: Used in the paper by Angluin et. al. in FoSSaCS 2022. The program generates statistics such as mean and median for the results of convert.java.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

public class statistics {
	
	public static int numFiles;
	public static int weighting;
	public static boolean inDFA;
	
	public static ArrayList<Integer>[][] resultsM2MA;
	public static ArrayList<Integer>[][] resultsDFA;
	public static int[][] counts;
	public static int[][] weights;

	public static void main(String[] args) throws Exception {
		Scanner in = new Scanner(System.in);
		
		System.out.println("Input number of files to read.");
		numFiles = Integer.parseInt(in.nextLine());
		
		String[] fileNames = new String[numFiles];
		for (int i = 0; i < numFiles; i++) {
			System.out.println("Input file name.");
			fileNames[i] = in.nextLine();
		}
		
		System.out.println("Enter 1 to give each file equal weight, or 2 to give each input automata equal weight.");
		weighting = Integer.parseInt(in.nextLine());
		in.close();
		
		inDFA = containsDFAResults(fileNames[0]);
		
		initialize();
		
		for (int i = 0; i < numFiles; i++) {
			readFile(fileNames[i], i);
		}
		
		displayResults(resultsM2MA, "M2MA");
		
		if (inDFA) {
			displayResults(resultsDFA, "DFA");
		}
	}
	
	public static boolean containsDFAResults(String firstFileName) throws IOException {
		BufferedReader firstFile = new BufferedReader(new FileReader(firstFileName));
		String line = firstFile.readLine();
		firstFile.close();
		
		if (line.contains("DFA")) {
			return true;
		}
		return false;
	}
	
	@SuppressWarnings("unchecked")
	public static void initialize() {
		resultsM2MA = new ArrayList[31][numFiles];
		for (int i = 1; i <= 30; i++) {
			for (int j = 0; j < numFiles; j++) {
				resultsM2MA[i][j] = new ArrayList<Integer>();
			}
		}
		
		resultsDFA = new ArrayList[31][numFiles];
		for (int i = 1; i <= 30; i++) {
			for (int j = 0; j < numFiles; j++) {
				resultsDFA[i][j] = new ArrayList<Integer>();
			}
		}
		
		counts = new int[31][numFiles];
	}
	
	public static void readFile(String fileName, int fileNumber) throws IOException {
		BufferedReader f = new BufferedReader(new FileReader(fileName));
		f.readLine();
		String line = f.readLine();
		
		if (line != null) {		
			String[] lineArr = line.split(" ");
			
			while (true) {
				int initialSize = Integer.parseInt(lineArr[0]);
						
				int convertedM2MASize = Integer.parseInt(lineArr[1]);
				resultsM2MA[initialSize][fileNumber].add(convertedM2MASize);
				
				if (inDFA) {
					int convertedDFASize = Integer.parseInt(lineArr[2]);
					resultsDFA[initialSize][fileNumber].add(convertedDFASize);
				}
				
				counts[initialSize][fileNumber]++;
				
				line = f.readLine();
				if (line == null) {
					break;
				}
				lineArr = line.split(" ");
			}
		}
		
		f.close();
	}
	
	public static void displayResults(ArrayList<Integer>[][] results, String convertedAutomataName) {
		System.out.println("\n" + convertedAutomataName + " Results");
		System.out.println("------------");
		System.out.println("Initial automata size: mean, median, standard deviation, min, max, number of automata");
		
		for (int i = 1; i <= 30; i++) {
			if (containsData(i)) {
				double[] weights = calculateWeights(i);
						
				double mean = calculateMean(results[i], weights, i);
				double median = calculateMedian(results[i], weights, i);
				double stdDev = calculateStandardDeviation(results[i], mean, weights, i);
				int min = calculateMin(results[i]);
				int max = calculateMax(results[i]);
				int count = calculateCount(i);
				
				System.out.println(i + ": " + round(mean, 2) + ", " + round(median, 2) + ", " 
						+ round(stdDev, 2) + ", " + min + ", " + max + ", " + count);
			}
		}
	}
	
	public static boolean containsData(int size) {
		for (int i = 0; i < counts[size].length; i++) {
			if (counts[size][i] > 0) {
				return true;
			}
		}
		return false;
	}
	
	public static double[] calculateWeights(int size) {
		double[] weights = new double[numFiles];
		
		if (weighting == 2) {
			for (int i = 0; i < numFiles; i++) {
				weights[i] = 1;
			}
			return weights;
		}
		
		double max = 0;
		for (int i = 0; i < numFiles; i++) {
			if (counts[size][i] > max) {
				max = counts[size][i];
			}
		}
		
		for (int i = 0; i < numFiles; i++) {
			if (counts[size][i] > 0) {
				weights[i] = max / counts[size][i];
			}
		}
		
		return weights;
	}
	
	public static double calculateMean(ArrayList<Integer>[] arr, double[] weights, int size) {
		double weightsSum = 0;
		double sum = 0;
		
		for (int i = 0; i < numFiles; i++) {
			if (counts[size][i] > 0) {
				for (int n : arr[i]) {
					weightsSum += weights[i];
					sum += n * weights[i];
				}
			}
		}
		
		return sum / weightsSum;
	}

	public static double calculateMedian(ArrayList<Integer>[] arr, double[] weights, int size) {
		if (weighting == 1) {
			int count = 0; 
			double sum = 0;
			
			for (int i = 0; i < numFiles; i++) {
				if (counts[size][i] > 0) {
					count++;
					sum += medianFormula(arr[i]);
				}
			}
			
			return sum / count;
		} else {
			ArrayList<Integer> aggregateArr = new ArrayList<Integer>();
			
			for (int i = 0; i < numFiles; i++) {
				if (counts[size][i] > 0) {
					for (int n : arr[i]) {
						aggregateArr.add(n);
					}
				}
			}
			
			return medianFormula(aggregateArr);
		}
	}
	
	public static double medianFormula(ArrayList<Integer> arr) {
		Collections.sort(arr);
		
		int len = arr.size();
		
		if (len % 2 == 0) {
			return (arr.get((len / 2) - 1) + arr.get(len / 2)) / 2.0;
		} else {
			return arr.get(len / 2);
		}
	}
	
	public static double calculateStandardDeviation(ArrayList<Integer>[] arr, double mean, double[] weights, int size) {
		double count = 0;
		double weightsSum = 0;
		double sum = 0;
		
		for (int i = 0; i < numFiles; i++) {
			if (counts[size][i] > 0) {
				for (int n : arr[i]) {
					count++;
					weightsSum += weights[i];
					sum += weights[i] * Math.pow(n - mean, 2);
				}
			}
		}
		
		return Math.sqrt(sum / (((count - 1) / count) * weightsSum));
	}
	
	public static int calculateMin(ArrayList<Integer>[] arr) {
		int min = Integer.MAX_VALUE;
		
		for (int i = 0; i < numFiles; i++) {
			for (int n : arr[i]) {
				if (n < min) {
					min = n;
				}
			}
		}
		
		return min;
	}
	
	public static int calculateMax(ArrayList<Integer>[] arr) {
		int max = 0;
		
		for (int i = 0; i < numFiles; i++) {
			for (int n : arr[i]) {
				if (n > max) {
					max = n;
				}
			}
		}
		
		return max;
	}

	public static double round(double num, int numDecimalPlaces) {
		double pow10 = Math.pow(10, numDecimalPlaces);
		return Math.round(pow10 * num) / pow10;
	}
	
	public static int calculateCount(int size) {
		int count = 0;
		for (int i = 0; i < numFiles; i++) {
			count += counts[size][i];
		}
		return count;
	}
}
