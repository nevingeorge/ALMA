## Algorithm for Learning Mod-2-Multiplicity Automata
Mod2_MA.java takes as input a mod-2-MA and prints to stdout the MA obtained after learning the input function through a series of membership and equivalence queries. The motivation behind this algorithm originally arose from Angluin's exact learning model described in her paper "Learning regular sets from queries and counterexamples."

### Format of Input File
The input file is a text document containing the specifications of the target function. It must have the following format (no line separation, characters are space separated, and lines beginning with // are ignored):

Line 1: alphabet size

Line 2: characters in the alphabet

Line 3: size of the target function (r)

Line 4: γ of the target function (fy)

Lines 5-end: List of μ's for each character in the alphabet, with each μ appearing in a rxr grid

Example input files can be found in the respository.

In Mod2_MA.java, change the file name in the initialize function to the name of the intended input file.

## Algorithm for Learning Strongly Unambiguous Büchi Automata (SUBA)
SUBA.java takes as input a SUBA and converts it into an equivalent UFA. The resulting UFA is then converted into an equivalent mod-2-MA and learned using Mod2_MA.java.

### Format of Input File
The input file is a text document containing the specifications of a SUBA of the form (Q, Σ, ∆, F). The file must have the following format (no line separation, characters are space separated, and lines beginning with // are ignored):

Line 1: number of states (Q)

Line 2: alphabet size

Line 3: characters in the alphabet

Line 4: final states (F)

Line 5: number of transitions

Lines 6-end: transitions - each line has the form q_i a q_j, where q_i,q_j∈Q and a∈Σ.

By default the only initial state of the SUBA (and therefore also the UFA) is q_1.

Example input files can be found in the GitHub repository.

In SUBA.java, change the file name in the SUBAtoUFA function to the name of the intended input file.
  
## Author: Nevin George

## Advisor: Dana Angluin

## References
Amos Beimel, Francesco Bergadano, Nader H. Bshouty, Eyal Kushilevitz, Stefano Varric- chio. Learning functions represented    as multiplicity automata. *J. ACM*, 47(3):506–530, May 2000.

Dana Angluin. Learning regular sets from queries and counterexamples. *Inf. Comput.*, 75(2):87–106, 1987.

Dana Angluin, Timos Antonopoulos, Dana Fisman. Strongly Unambiguous Büchi Automata Are Polynomially Predictable with Membership Queries. *28th International Conference on Computer Science Logic*, 8:1–8:17, 2020.

Michael Thon and Herbert Jaeger. Links Between Multiplicity Automata, Observable Operator Models and Predictive State Representations — a Unified Learning Framework. *Journal of Machine Learning Research*, 16(4):103−147, 2015.

N. Bousquet and C. Löding. Equivalence and inclusion problem for strongly unambiguous büchi automata. In *Language and Automata Theory and Applications, 4th International Conference, LATA 2010, Trier, Germany, May 24-28, 2010. Proceedings,* pages 118–129, 2010. URL: https: //doi.org/10.1007/978-3-642-13089-2_10, doi:10.1007/978-3-642-13089-2\_10.
