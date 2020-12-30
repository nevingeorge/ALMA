import java.util.HashMap;

public class DBAtoDFA {
	
	public static void main(String[] args) throws Exception {
		arbitrary.EQMaxTestLen = 30;	
		arbitrary.EQNumTests = 1000;
		arbitrary.EQLimit = 100000;
		
		M2MA.alphabet = new String[4];
		M2MA.alphabet[0] = "a";
		M2MA.alphabet[1] = "b";
		M2MA.alphabet[2] = "c";
		M2MA.alphabet[3]	= "$";			
		
		M2MA.letterToIndex = new HashMap<String, Integer>();
		for (int i=0; i<M2MA.alphabet.length; i++) {
			M2MA.letterToIndex.put(M2MA.alphabet[i], i);
		}
		
		// NBAs of size 1 to 30, int[][0] = count, int[][1] = sum of M2MA dimension
		NBAtoM2MA.results = new int[31][2];

		System.out.println("0% complete");

		NBAtoM2MA.run(20, 5, 10, 2, false, false);
		NBAtoM2MA.run(20, 5, 11, 2, false, false);
		NBAtoM2MA.run(20, 5, 12, 2, false, false);
		NBAtoM2MA.run(20, 5, 13, 2, false, false);
		NBAtoM2MA.run(20, 5, 14, 2, false, false);
		NBAtoM2MA.run(20, 5, 15, 2, false, false);
		
		System.out.println("12% complete");
		
		NBAtoM2MA.run(20, 5, 10, 3, false, false);
		NBAtoM2MA.run(20, 5, 11, 3, false, false);
		NBAtoM2MA.run(20, 5, 12, 3, false, false);
		NBAtoM2MA.run(20, 5, 13, 3, false, false);
		NBAtoM2MA.run(20, 5, 14, 3, false, false);
		NBAtoM2MA.run(20, 5, 15, 3, false, false);
		
		System.out.println("24% complete");
		
		NBAtoM2MA.run(20, 10, 20, 3, false, false);
		NBAtoM2MA.run(20, 10, 22, 3, false, false);
		NBAtoM2MA.run(20, 10, 24, 3, false, false);
		NBAtoM2MA.run(20, 10, 26, 3, false, false);
		NBAtoM2MA.run(20, 10, 28, 3, false, false);
		NBAtoM2MA.run(20, 10, 30, 3, false, false);
		
		System.out.println("36% complete");
		
		NBAtoM2MA.run(20, 10, 20, 4, false, false);
		NBAtoM2MA.run(20, 10, 22, 4, false, false);
		NBAtoM2MA.run(20, 10, 24, 4, false, false);
		NBAtoM2MA.run(20, 10, 26, 4, false, false);
		NBAtoM2MA.run(20, 10, 28, 4, false, false);
		NBAtoM2MA.run(20, 10, 30, 4, false, false);
		
		System.out.println("48% complete");
		
		NBAtoM2MA.run(20, 15, 30, 4, false, false);
		NBAtoM2MA.run(20, 15, 33, 4, false, false);
		NBAtoM2MA.run(20, 15, 36, 4, false, false);
		NBAtoM2MA.run(20, 15, 39, 4, false, false);
		NBAtoM2MA.run(20, 15, 42, 4, false, false);
		NBAtoM2MA.run(20, 15, 45, 4, false, false);
		
		System.out.println("60% complete");
		
		NBAtoM2MA.run(20, 15, 30, 5, false, false);
		NBAtoM2MA.run(20, 15, 33, 5, false, false);
		NBAtoM2MA.run(20, 15, 36, 5, false, false);
		NBAtoM2MA.run(20, 15, 39, 5, false, false);
		NBAtoM2MA.run(20, 15, 42, 5, false, false);
		NBAtoM2MA.run(20, 15, 45, 5, false, false);
		
		System.out.println("72% complete");
	
		NBAtoM2MA.run(20, 20, 42, 5, false, false);
		NBAtoM2MA.run(20, 20, 45, 5, false, false);
		NBAtoM2MA.run(20, 20, 48, 5, false, false);
		NBAtoM2MA.run(20, 20, 51, 5, false, false);
		NBAtoM2MA.run(20, 20, 54, 5, false, false);
		NBAtoM2MA.run(20, 20, 57, 5, false, false);
		NBAtoM2MA.run(20, 20, 60, 5, false, false);
		
		System.out.println("86% complete");
		
		NBAtoM2MA.run(20, 20, 42, 6, false, false);
		NBAtoM2MA.run(20, 20, 45, 6, false, false);
		NBAtoM2MA.run(20, 20, 48, 6, false, false);
		NBAtoM2MA.run(20, 20, 51, 6, false, false);
		NBAtoM2MA.run(20, 20, 54, 6, false, false);
		NBAtoM2MA.run(20, 20, 57, 6, false, false);
		NBAtoM2MA.run(20, 20, 60, 6, false, false);

		System.out.println("100% complete\n");
		
		System.out.println("DBA to DFA results.");
		for (int i = 1; i <= 30; i++) {
			if (NBAtoM2MA.results[i][0] != 0) {
				System.out.println(i + ": " + NBAtoM2MA.results[i][0] + " " + NBAtoM2MA.results[i][1] + " " + ((double) NBAtoM2MA.results[i][1])/NBAtoM2MA.results[i][0]);
			}
		}
	}
}
