import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

public class statistics {

	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws Exception {
		Scanner in = new Scanner(System.in);
		
		System.out.println("Input number of files to read.");
		int numFiles = Integer.parseInt(in.nextLine());
		
		String[] fileNames = new String[numFiles];
		for (int i = 0; i < numFiles; i++) {
			System.out.println("Input file name.");
			fileNames[i] = in.nextLine();
		}
		
		int[][] meanResultsM2MA = new int[31][2];
		int[][] meanResultsDFA = new int[31][2];
		
		ArrayList<Integer>[] medianResultsM2MA = new ArrayList[31];
		for (int i = 1; i <= 30; i++) {
			medianResultsM2MA[i] = new ArrayList<Integer>();
		}
		
		ArrayList<Integer>[] medianResultsDFA = new ArrayList[31];
		for (int i = 1; i <= 30; i++) {
			medianResultsDFA[i] = new ArrayList<Integer>();
		}
		
		boolean DFAresults = false;
		
		for (int i = 0; i < numFiles; i++) {
			BufferedReader f = new BufferedReader(new FileReader(fileNames[i]));
			f.readLine();
			String line = f.readLine();
			
			if (line != null) {		
				boolean inDFA = false;
				String[] lineArr = line.split(" ");
				if (lineArr.length == 3) {
					inDFA = true;
					DFAresults = true;
				}
				
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
		in.close();
		
		System.out.println("\nExperiment Results");
		System.out.println("------------------");
		
		
		System.out.println("Initial automata size: mean size of converted M2MA, median size of converted M2MA");
		for (int i = 1; i <= 30; i++) {
			if (meanResultsM2MA[i][0] > 0) {
				int mean = (int) Math.round(meanResultsM2MA[i][1] / (double) meanResultsM2MA[i][0]);
				
				Collections.sort(medianResultsM2MA[i]);
				
				int len = medianResultsM2MA[i].size();
				int median;
				
				if (len % 2 == 0) {
					median = (int) Math.round((medianResultsM2MA[i].get((len / 2) - 1) + medianResultsM2MA[i].get(len / 2)) / 2.0);
				} else {
					median = medianResultsM2MA[i].get(len / 2);
				}
				
				System.out.println(i + ": " + mean + ", " + median);
			}
		}
		
		if (DFAresults) {
			System.out.println("\nInitial automata size: mean size of converted DFA, median size of converted DFA");
			for (int i = 1; i <= 30; i++) {
				if (meanResultsDFA[i][0] > 0) {
					int mean = (int) Math.round(meanResultsDFA[i][1] / (double) meanResultsDFA[i][0]);
					
					Collections.sort(medianResultsDFA[i]);
					
					int len = medianResultsDFA[i].size();
					int median;
					
					if (len % 2 == 0) {
						median = (int) Math.round((medianResultsDFA[i].get((len / 2) - 1) + medianResultsDFA[i].get(len / 2)) / 2.0);
					} else {
						median = medianResultsDFA[i].get(len / 2);
					}
					
					System.out.println(i + ": " + mean + ", " + median);
				}
			}
		}
	}

}
