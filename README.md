To run the programs, in Terminal/Command Prompt change your current directory to the Java Algorithms folder and input the command "java -jar" followed by "Mod2_MA.jar", "SUBA.jar", or "arbitrary.jar". The program will then ask for the input file name and an optional flag -v. Make sure the desired input file is in the same folder as the jar file. If the optional flag is entered, the program will display the observation table as it is being constructed.

## Algorithm for learning mod-2-multiplicity automata
Mod2_MA.java takes in as input a mod-2-MA and prints to stdout the MA obtained after learning the input function through a series of membership and equivalence queries.

### Format of Input File
The input file is a text document containing the specifications of the target function. It must have the following format (no line separation, entries are space separated, and lines beginning with // are ignored):

Line 1: characters in the alphabet

Line 2: size of the target function (r)

Line 3: γ of the target function (fy)

Lines 4-end: List of μ's for each character in the alphabet, with each μ appearing in a rxr grid

Example input files can be found in the repository.

## Algorithm for learning strongly unambiguous Büchi automata (SUBA)
SUBA.java takes in as input a SUBA of n states and converts it into an equivalent UFA of 2n<sup>2</sup>+n states. The UFA is then converted into an equivalent mod-2-MA of the same size and learned using Mod2_MA.java.

### Format of Input File
The input file is a text document containing the specifications of a SUBA of the form (Q, Σ, ∆, F). The file must have the following format (no line separation, entries are space separated, and lines beginning with // are ignored):

Line 1: number of states (Q)

Line 2: characters in the alphabet

Line 3: final states (F)

Line 4: number of transitions

Lines 5-end: transitions - each line has the form q_i a q_j, where q_i,q_j∈Q and a∈Σ.

By default the only initial state of the SUBA (and therefore also the UFA) is q_1.

Example input files can be found in the repository.

## Algorithm for learning arbitrary automata
arbitrary.java displays to stdout the mod-2-MA learned using a membership query method specified in MQ.java and statistical equivalence queries. The program can be used to approximately learn any type of automata, provided that MQ.java contains the desired automata's membership query function.

### Format of Input File
The input file is a text document that specifies the desired membership query function in MQ.java and level of approximation for the statistical equivalence queries. The file must have the following format (no line separation, entries are space separated, and lines beginning with // are ignored):

Line 1: characters in the alphabet

Line 2: name of the desired membership query function in MQ.java

Line 3: maximum length of a test in the statistical equivalence query

Line 4: number of tests the statistical equivalence query will check

Line 5: limit on the number of equivalence queries to run

Example input files can be found in the repository.

## Author: Nevin George

## Advisor: Dana Angluin

## References
Amos Beimel, Francesco Bergadano, Nader H. Bshouty, Eyal Kushilevitz, Stefano Varric- chio. Learning functions represented    as multiplicity automata. *J. ACM*, 47(3):506–530, May 2000.

Dana Angluin. Learning regular sets from queries and counterexamples. *Inf. Comput.*, 75(2):87–106, 1987.

Dana Angluin, Timos Antonopoulos, Dana Fisman. Strongly Unambiguous Büchi Automata Are Polynomially Predictable with Membership Queries. *28th International Conference on Computer Science Logic*, 8:1–8:17, 2020.

Michael Thon and Herbert Jaeger. Links Between Multiplicity Automata, Observable Operator Models and Predictive State Representations — a Unified Learning Framework. *Journal of Machine Learning Research*, 16(4):103−147, 2015.

N. Bousquet and C. Löding. Equivalence and inclusion problem for strongly unambiguous büchi automata. In *Language and Automata Theory and Applications, 4th International Conference, LATA 2010, Trier, Germany, May 24-28, 2010. Proceedings,* pages 118–129, 2010. URL: https: //doi.org/10.1007/978-3-642-13089-2_10, doi:10.1007/978-3-642-13089-2\_10.
