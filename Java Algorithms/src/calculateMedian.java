import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

public class calculateMedian {

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
		
		int[][] meanResults = new int[31][2];
		
		ArrayList<Integer>[] medianResults = new ArrayList[31];
		for (int i = 1; i <= 30; i++) {
			medianResults[i] = new ArrayList<Integer>();
		}
		
		for (int i = 0; i < numFiles; i++) {
			BufferedReader f = new BufferedReader(new FileReader(fileNames[i]));
			f.readLine();
			String line = f.readLine();
			
			if (line != null) {		
				boolean inM2MA = false;
				
				String[] lineArr = line.split(" ");
				if (lineArr.length == 2) {
					inM2MA = true;
				}
				
				while (true) {
					int initialSize = Integer.parseInt(lineArr[0]);
					meanResults[initialSize][0]++;
					
					int convertedSize;
					if (inM2MA) {
						convertedSize = Integer.parseInt(lineArr[1]);
					} else {
						convertedSize = Integer.parseInt(lineArr[2]);
					}
					
					meanResults[initialSize][1] += convertedSize;
					medianResults[initialSize].add(convertedSize);
					
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
		
		System.out.println("Initial automata size: mean size of converted automata, median of converted automata size");
		for (int i = 1; i <= 30; i++) {
			if (meanResults[i][0] > 0) {
				int mean = (int) Math.round(meanResults[i][1] / (double) meanResults[i][0]);
				
				Collections.sort(medianResults[i]);
				
				int len = medianResults[i].size();
				int median;
				
				if (len % 2 == 0) {
					median = (int) Math.round((medianResults[i].get((len / 2) - 1) + medianResults[i].get(len / 2)) / 2.0);
				} else {
					median = medianResults[i].get(len / 2);
				}
				
				System.out.println(i + ": " + mean + ", " + median);
			}
		}
	}

}
