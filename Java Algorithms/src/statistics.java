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
	
	public static ArrayList<Integer>[] resultsM2MA;
	public static ArrayList<Integer>[] resultsDFA;

	public static void main(String[] args) throws Exception {
		Scanner in = new Scanner(System.in);
		
		System.out.println("Input number of files to read.");
		int numFiles = Integer.parseInt(in.nextLine());
		
		String[] fileNames = new String[numFiles];
		for (int i = 0; i < numFiles; i++) {
			System.out.println("Input file name.");
			fileNames[i] = in.nextLine();
		}
		in.close();
		
		inDFA = containsDFAResults(fileNames[0]);
		
		initialize();
		
		for (int i = 0; i < numFiles; i++) {
			readFile(fileNames[i]);
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
		resultsM2MA = new ArrayList[31];
		for (int i = 1; i <= 30; i++) {
			resultsM2MA[i] = new ArrayList<Integer>();
		}
		
		resultsDFA = new ArrayList[31];
		for (int i = 1; i <= 30; i++) {
			resultsDFA[i] = new ArrayList<Integer>();
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
				resultsM2MA[initialSize].add(convertedM2MASize);
				
				if (inDFA) {
					int convertedDFASize = Integer.parseInt(lineArr[2]);
					resultsDFA[initialSize].add(convertedDFASize);
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
	
	public static void displayResults(ArrayList<Integer>[] results, String convertedAutomataName) {
		System.out.println("\n" + convertedAutomataName + " Results");
		System.out.println("------------");
		System.out.println("Initial automata size: mean, median, standard deviation");
		for (int i = 1; i <= 30; i++) {
			if (results[i].size() > 0) {
				double mean = calculateMean(results[i]);
				double median = calculateMedian(results[i]);
				double stdDev = calculateStandardDeviation(results[i], mean);
				
				System.out.println(i + ": " + mean + ", " + median + ", " + stdDev);
			}
		}
	}
	
	
	public static double calculateMean(ArrayList<Integer> arr) {
		int count = 0;
		double sum = 0;
		
		for (int n : arr) {
			count++;
			sum += n;
		}
		
		double mean = sum / count;
		
		return Math.round(mean * 100) / 100.0;
	}

	public static double calculateMedian(ArrayList<Integer> arr) {
		Collections.sort(arr);
		
		int len = arr.size();
		
		if (len % 2 == 0) {
			return (arr.get((len / 2) - 1) + arr.get(len / 2)) / 2.0;
		} else {
			return arr.get(len / 2);
		}
	}
	
	public static double calculateStandardDeviation(ArrayList<Integer> arr, double mean) {
		int count = 0;
		double sum = 0;
		
		for (int n : arr) {
			count++;
			sum += Math.pow(n - mean, 2);
		}
		
		double stdDev = Math.sqrt(sum / (count - 1));
		
		return Math.round(stdDev * 100) / 100.0;
	}
}
