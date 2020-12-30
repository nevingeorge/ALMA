import java.util.HashMap;

public class NBAtoDFA {

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

		NBAtoM2MA.run(5, 5, 10, 2, true, false);
		NBAtoM2MA.run(5, 5, 15, 2, true, false);
		NBAtoM2MA.run(5, 5, 20, 2, true, false);
		NBAtoM2MA.run(5, 5, 25, 2, true, false);
		NBAtoM2MA.run(5, 5, 30, 2, true, false);
		
		System.out.println("10% complete");
		
		NBAtoM2MA.run(20, 5, 10, 3, true, false);
		NBAtoM2MA.run(20, 5, 15, 3, true, false);
		NBAtoM2MA.run(20, 5, 20, 3, true, false);
		NBAtoM2MA.run(20, 5, 25, 3, true, false);
		NBAtoM2MA.run(20, 5, 30, 3, true, false);
		
		System.out.println("20% complete");
		
		NBAtoM2MA.run(20, 10, 20, 3, true, false);
		NBAtoM2MA.run(20, 10, 25, 3, true, false);
		NBAtoM2MA.run(20, 10, 30, 3, true, false);
		NBAtoM2MA.run(20, 10, 35, 3, true, false);
		NBAtoM2MA.run(20, 10, 40, 3, true, false);
		NBAtoM2MA.run(20, 10, 45, 3, true, false);
		NBAtoM2MA.run(20, 10, 50, 3, true, false);
		
		System.out.println("34% complete");
		
		NBAtoM2MA.run(20, 10, 20, 4, true, false);
		NBAtoM2MA.run(20, 10, 25, 4, true, false);
		NBAtoM2MA.run(20, 10, 30, 4, true, false);
		NBAtoM2MA.run(20, 10, 35, 4, true, false);
		NBAtoM2MA.run(20, 10, 40, 4, true, false);
		NBAtoM2MA.run(20, 10, 45, 4, true, false);
		NBAtoM2MA.run(20, 10, 50, 4, true, false);
		
		System.out.println("48% complete");
		
		NBAtoM2MA.run(20, 15, 30, 4, true, false);
		NBAtoM2MA.run(20, 15, 40, 4, true, false);
		NBAtoM2MA.run(20, 15, 50, 4, true, false);
		NBAtoM2MA.run(20, 15, 60, 4, true, false);
		NBAtoM2MA.run(20, 15, 70, 4, true, false);
		NBAtoM2MA.run(20, 15, 80, 4, true, false);
		
		System.out.println("60% complete");
		
		NBAtoM2MA.run(20, 15, 30, 5, true, false);
		NBAtoM2MA.run(20, 15, 40, 5, true, false);
		NBAtoM2MA.run(20, 15, 50, 5, true, false);
		NBAtoM2MA.run(20, 15, 60, 5, true, false);
		NBAtoM2MA.run(20, 15, 70, 5, true, false);
		NBAtoM2MA.run(20, 15, 80, 5, true, false);
		
		System.out.println("72% complete");
		
		NBAtoM2MA.run(20, 20, 40, 5, true, false);
		NBAtoM2MA.run(20, 20, 50, 5, true, false);
		NBAtoM2MA.run(20, 20, 60, 5, true, false);
		NBAtoM2MA.run(20, 20, 70, 5, true, false);
		NBAtoM2MA.run(20, 20, 80, 5, true, false);
		NBAtoM2MA.run(20, 20, 90, 5, true, false);
		NBAtoM2MA.run(20, 20, 100, 5, true, false);
		
		System.out.println("86% complete");
		
		NBAtoM2MA.run(20, 20, 40, 6, true, false);
		NBAtoM2MA.run(20, 20, 50, 6, true, false);
		NBAtoM2MA.run(20, 20, 60, 6, true, false);
		NBAtoM2MA.run(20, 20, 70, 6, true, false);
		NBAtoM2MA.run(20, 20, 80, 6, true, false);
		NBAtoM2MA.run(20, 20, 90, 6, true, false);
		NBAtoM2MA.run(20, 20, 100, 6, true, false);
	
		System.out.println("100% complete\n");
		
		System.out.println("NBA to DFA results.");
		for (int i = 1; i <= 30; i++) {
			if (NBAtoM2MA.results[i][0] != 0) {
				System.out.println(i + ": " + NBAtoM2MA.results[i][0] + " " + NBAtoM2MA.results[i][1] + " " + ((double) NBAtoM2MA.results[i][1])/NBAtoM2MA.results[i][0]);
			}
		}
	}
}
