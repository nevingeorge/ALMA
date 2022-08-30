# ALMA: Automata Learner using Modulo 2 Multiplicity Automata
To run the programs, in Terminal/Command Prompt change your current directory to the Java Algorithms folder and input the command "java -jar" followed by the name of the desired jar file. The program will then ask for the input file name and optional flags. Make sure the desired input file is in the same folder as the jar file. Each input file must be a text document following the format required of its intended program (specific details described below). Each file must have no line separation, entries must be space separated, and lines beginning with // are ignored. Example input files for all of the programs can be found in the repository.

Optional flags:\
-v - display more verbose information regarding the procedures and outputs of the algorithms\
-m - display the progress of the minimization algorithm\
-d - only display the dimension of the minimized M2MA\
-a - display the number of states of a minimal DFA equivalent to the minimized M2MA

## Learning modulo 2 multiplicity automata (M2MA)
M2MA.java takes in as input an M2MA and prints to stdout the M2MA obtained after learning the input function through a series of membership and equivalence queries.

### Input File Format
Contains the specifications of an M2MA.

Line 1: alphabet

Line 2: dimension

Line 3: final vector

Lines 4-end: transition matrices for each character in the alphabet

By default, the initial vector is (1,0,0,...,0).

## Learning strongly unambiguous Büchi automata (SUBA)
SUBA.java takes in as input a SUBA of n states and converts it into an equivalent UFA of 2n<sup>2</sup>+n states. The UFA is then converted into an equivalent M2MA of the same size and learned using M2MA.java.

### Input File Format
Contains the specifications of a SUBA of the form (Q, Σ, ∆, F).

Line 1: number of states

Line 2: alphabet

Line 3: final states

Line 4: number of transitions

Lines 5-end: transitions - each line has the form q_i a q_j, where q_i,q_j∈Q and a∈Σ.

By default, the only initial state of the SUBA (and therefore also the UFA) is q_1.

## Learning non-deterministic Büchi automata (NBA)
NBA.java takes in as input an NBA and prints to stdout the M2MA obtained after learning the NBA through a series of membership and statistical equivalence queries.

### Input File Format
Contains the specifications of an NBA of the form (Q, Σ, ∆, F) and the desired level of approximation for the statistical equivalence queries.

Line 1: maximum length of a test in the statistical equivalence query

Line 2: number of tests the statistical equivalence query will check

Line 3: limit on the number of equivalence queries to run

Line 4: number of states

Line 5: alphabet

Line 6: final states

Line 7: number of transitions

Lines 8-end: transitions - each line has the form q_i a q_j, where q_i,q_j∈Q and a∈Σ.

By default the only initial state of the NBA is q_1.

## Learning arbitrary automata
arbitrary.java displays to stdout the M2MA learned using a membership query method specified in MQ.java and statistical equivalence queries. The program can be used to approximately learn any type of automata, provided that MQ.java contains the desired automata's membership query function.

### Input File Format
Contains the name of the desired membership query function in MQ.java and level of approximation for the statistical equivalence queries.

Line 1: name of the desired membership query function in MQ.java

Line 2: maximum length of a test in the statistical equivalence query

Line 3: number of tests the statistical equivalence query will check

Line 4: limit on the number of equivalence queries to run

Line 5: alphabet

## Minimizing automata
minimize.java takes in as input an M2MA or SUBA and prints to stdout the M2MA obtained after minimizing the input function (in the SUBA case, it first converts the function into an equivalent UFA then M2MA). The format for the M2MA/SUBA inputs were described earlier.

## Converting SUBA, NBA, and DBA to M2MA and DFA
convert.jar takes in as input a series of SUBA, NBA, or DBA. For each input automaton, the program adds to an output file the size of a minimal M2MA or DFA that accepts the same language. The program displays to stdout the average converted M2MA/DFA size for each input omega automaton size. Also, statistics.jar can be used to obtain more detailed statistics on the results of multiple output files representing the same conversion (e.g. SUBA->DFA or NBA->M2MA).

### SUBA Input File Format
Line 1: number of SUBA in the input file

Lines 2-end: input SUBA - follows the same format as SUBA.jar

### NBA/DBA Input File Format
Line 1: maximum length of a test in the statistical equivalence query

Line 2: number of tests the statistical equivalence query will check

Line 3: limit on the number of equivalence queries to run

Line 4: alphabet

Line 5: number of lines to follow

Lines 6-end: lines of the form (number of NBA/DBA to generate, max number of states, max number of transitions to remove, number of final states)

## Author: Nevin George

## Advisor: Dana Angluin

## References
Angluin, D. Learning regular sets from queries and counterexamples. Inf. Comput., 75(2):87–106 (1987)

Angluin, D. Queries and concept learning. \textit{Mach. Learn. 2}, 4, 319–342 (1988)

Angluin, D., Antonopoulos, T., Fisman, D.:Strongly unambiguous Büchi automata are polynomially predictable with membership queries. In: 28th EACSL Annual Conference on Computer Science Logic, CSL. pp. 8:1–8:17 (2020)

Angluin, D., Antonopoulos, T., Fisman, D., George, N.:Representing Regular Languages of Infinite Words Using Mod 2 Multiplicity Automata. FoSSaCS (2022)

Beimel, A., Bergadano, F., Bshouty, N.H., Kushilevitz, E., Varricchio, S.: Learning
functions represented as multiplicity automata. J. ACM 47(3), 506–530 (May 2000)

Bousquet, N., Löding, C.: Equivalence and inclusion problem for strongly unambiguous Büchi automata. In: Language and Automata Theory and Applications, 4th International Conference, LATA. Proceedings. pp. 118–129 (2010).

Büchi, J.R.: On a decision method in restricted second order arithmetic. In: International Congress on Logic, Methodology and Philosophy of Science, Stanford University Press (1962) 1–11

Calbrix, H., Nivat, M., Podelski, A.: Ultimately periodic words of rational $\omega$-languages. In: Proceedings of the 9th International Conference on Mathematical Foundations of Programming Semantics. pp. 554–566. Springer-Verlag (1994)

Carton, O., Michel, M.: Unambiguous Büchi automata. Theor. Comput. Sci. 297(1–3) (2003) 37–81

Michael R. Thon and Herbert Jaeger. Links between multiplicity automata, observable operator models and predictive state representations: a unified learning framework. \textit{Journal of Machine Learning Research}, 16:103–147 (2015)

Sakarovitch, J.: Elements of Automata Theory. Cambridge University Press, USA (2009)
Amos Beimel, Francesco Bergadano, Nader H. Bshouty, Eyal Kushilevitz, Stefano Varric- chio. Learning functions represented    as multiplicity automata. *J. ACM*, 47(3):506–530, May 2000.