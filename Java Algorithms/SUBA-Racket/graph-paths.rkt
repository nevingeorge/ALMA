#lang racket
; Procedures that take an edge-labeled directed graph, two vertices v1 and v2 and a length,
; and return whether there exists a directed path from v1 to v2 of that length, and if so,
; return a randomly selected path of that length from v1 to v2 or a sequence of edge labels
; for such a randomly selected path.

; Dana Angluin
; October 2020

(provide
 aut->digraph
 random-path random-path-labels)

(require "general-utilities.rkt")
(require "aut-definitions.rkt")
(require "aut-procedures.rkt")

; example aut with 3 states and 8 transitions, nondeterministic
(define ex-aut 
(aut '(a b c)
     '(1 2 3)
     '(1)
     (list
      (entry '(2 c) 1)
      (entry '(3 b) 1)
      (entry '(3 c) 1)
      (entry '(2 a) 2)
      (entry '(2 b) 3)
      (entry '(1 a) 2)
      (entry '(1 b) 3)
      (entry '(3 b) 2))))

; An edge-labeled directed graph is represented as an object (Racket-style: a procedure)
; with a set of vertices V, and a set of labeled directed edges E.
; V is a list, and E is represented by two hash tables:
; a hash table E-labels mapping the list (v1 v2) to the list of labels
; (which may have duplicates) for the edge from v1 to v2, and
; a hash table E-counts mapping the list (v1 v2) to the number of labels
; of edges from v1 to v2.
; Edges with no labels are not recorded.
; There is a auxiliary list of powers of the relation E-counts, which give the
; numbers of paths of a given length from v1 to v2, to avoid recomputing it.
; More powers are added to this list as necessary.

; Create a labeled directed graph object from an automaton.
; V is the list of automaton states.
; E-labels is the hash table for entries (state1 state2) giving the list of symbols
; (possibly with duplicates) in transitions ((state1 symbol) state2).
; E-counts is the hash table for entries (state1 state2), where the value
; is the number of alphabet symbols with transitions ((state1 symbol) state2).

(define (aut->digraph aut1)
  (let* [(V (aut-states aut1))
         (trans (aut-trans aut1))
         (E-labels (make-hash))
         (E-counts (make-hash))
         (path-table-list (list E-counts))]
    (for-each
     (lambda (tuple)
       (let* [(source-state (first (entry-key tuple)))
              (symbol (second (entry-key tuple)))
              (dest-state (entry-value tuple))]
         (hash-set! E-labels
                    (list source-state dest-state)
                    (cons symbol (hash-ref E-labels (list source-state dest-state) '())))
         (hash-set! E-counts
                    (list source-state dest-state)
                    (+ 1 (hash-ref E-counts (list source-state dest-state) 0)))))
     trans)
    (lambda (cmd . args)
      (case cmd
        [(show)
         (display "vertices: ")(displayln V)
         (displayln "hash table of edge labels:")
         (displayln E-labels)
         (displayln "hash table of edge counts:")
         (displayln E-counts)]
        [(get-V)
         V]
        [(get-labels)
         (hash-ref E-labels (first args) '())]
        [(get-tables)
         path-table-list]
        [(n-paths)
         (let* [(v1 (first args))
                (v2 (second args))
                (len (third args))]
           (if (= len 0)
               (if (equal? v1 v2) 1 0)
               (begin
                 (set! path-table-list
                       (extend-powers-of V path-table-list len))
                 (hash-ref (list-ref path-table-list (- len 1)) (list v1 v2) 0))))]
        [else
         (display "digraph unknown command: ")
         (displayln (cons cmd args))]))))
    
;> (define g1 (aut->digraph ex-aut))
;> (g1 'show)
;vertices: (1 2 3)
;hash table of edge labels:
;#hash(((3 1) . (c b)) ((3 2) . (b)) ((2 3) . (b)) ((2 2) . (a)) ((2 1) . (c)) ((1 2) . (a)) ((1 3) . (b)))
;hash table of edge counts:
;#hash(((3 1) . 2) ((3 2) . 1) ((2 3) . 1) ((2 2) . 1) ((2 1) . 1) ((1 2) . 1) ((1 3) . 1))


; Given a digraph, two vertices v1 and v2 and a  positive length len,
; create (if necessary) the tables for paths of lengths up to len,
; and return either 'none, if there are no paths of that length from
; v1 to v2 in the digraph, or a randomly chosen such path as a list
; of edges ((v1 w2) (w2 w3) ... (wk v2)).

(define (random-path graph v1 v2 len)
  (let* [(result (random-path-rec graph v1 v2 len))]
    (if (equal? result 'none)
        'none
        (reverse result))))
 
(define (random-path-rec graph v1 v2 len)
  (let* [(V (graph 'get-V))]
;    (display "random-path-rec ")
;    (displayln (list v1 v2 len))
    (cond
      [(zero? (graph 'n-paths v1 v2 len))
        'none]
      [(= len 0)
       '()]
      [(= len 1)
       (list (list v1 v2))]
      [else
       (let* [(products
               (map (lambda (int-v)
                      (* (graph 'n-paths v1 int-v (- len 1))
                         (graph 'n-paths int-v v2 1)))
                    V))
              (v-choice (random-select-counts V products))]
         (cons (list v-choice v2)
               (random-path-rec
                graph
                v1
                v-choice
                (- len 1))))])))

; Given a digraph, two vertices v1 and v2 and a  positive length len,
; call (random-path graph v1 v2 len) to get either 'none or
; a random path, and, if it is a path, return a list of labels
; randomly drawn for the edges on the path.

(define (random-path-labels graph v1 v2 len)
  (let* [(result (random-path graph v1 v2 len))]
    (if (equal? result 'none)
        'none
        (map (lambda (edge)
               (pick (graph 'get-labels edge)))
             result))))

; Given two (hash) tables that have the same sets rows and cols
; multiply them and return the result as a hash table.
; If table1 gives the number of paths of length 1 and
; table2 gives the number of paths of length k, then multiplying table1 times
; table2 gives the number of paths of length k+1.

(define (multiply-tables V table1 table2)
  (let* [(new-table (make-hash))]
    (for-each
     (lambda (row)
       (for-each
        (lambda (col)
          (for-each
           (lambda (k)
             (let* [(v1 (hash-ref table1 (list row k) 0))
                    (v2 (hash-ref table2 (list k col) 0))
                    (product (* v1 v2))]
               (if (not (zero? product))
                   (hash-set! new-table
                              (list row col)
                              (+ product (hash-ref new-table (list row col) 0)))
                   #t)))
           V))
        V))
     V)
    new-table))

; Extend a given nonempty list of consecutive powers of a transition relation
; (A^1 A^2 ... A^k) to length len, that is, (A^1 A^2 ... A^{len})

(define (extend-powers-of V powers-lst len)
  (let *[(n-powers (length powers-lst))]
    (if (<= len n-powers)
        powers-lst
        (powers-of-loop
         (first powers-lst) V (reverse powers-lst) (- len n-powers)))))

(define (powers-of-loop table1 V r-powers-lst remaining-length)
  (if (<= remaining-length 0)
      (reverse r-powers-lst)
      (powers-of-loop table1
                      V
                      (cons (multiply-tables
                             V table1 (first r-powers-lst))
                            r-powers-lst)
                      (- remaining-length 1))))