# Algorithm for Learning Mod-2-MA

The algorithm takes in as input a mod-2-MA and prints to stdout the MA obtained after learning the input function through a series of membership and equivalence queries. The algorithm is explained in more detail in the paper "Learning functions represented as multiplicity automata" by Beimel et al. The motivation behind this algorithm originally arose from Angluin's exact learning model described in her paper "Learning regular sets from queries and counterexamples."

### Format of Input File
The input file is a text document containing the specifications of the target function. It must have the following format (there is no separation between the lines and all the characters are space separated):

<size of the alphabet>
  
<characters in the alphabet>
  
<size of the target function (r)>

<γ of the target function (fy)>

Rest of the file:
List of μ's for each character in the alphabet, with each μ appears in a rxr grid

### Author: Nevin George

### Advisor: Dana Angluin

### References
Amos Beimel, Francesco Bergadano, Nader H. Bshouty, Eyal Kushilevitz, Stefano Varric- chio. Learning functions represented    as multiplicity automata. *J. ACM*, 47(3):506–530, May 2000.

Dana Angluin. Learning regular sets from queries and counterexamples. *Inf. Comput.*, 75(2):87–106, 1987.

Dana Angluin, Timos Antonopoulos, Dana Fisman. Strongly Unambiguous Büchi Automata Are Polynomially Predictable with Membership Queries. *28th International Conference on Computer Science Logic*, 8:1–8:17, 2020.
