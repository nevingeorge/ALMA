#lang racket
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Automaton procedures: DFAs, NFAs, NBAs

(provide
 display-aut display-aut-f
 random-aut random-reverse-det-aut
 random-aut-f random-aut-f-k
 aut-next-states aut-string-states
 aut-reachable aut-reaches
 reverse-transition deterministic-trans? complete-trans?
 aut-symbols->letters)

(require "general-utilities.rkt")
(require "aut-definitions.rkt")
(require "sequence-utilities.rkt")

; display an automaton (aut)

(define (display-aut a)
  (display "alphabet: ")(display (aut-alphabet a))(newline)
  (display "states: ")(display (aut-states a))(newline)
  (display "initial states: ")(display (aut-init-states a))(newline)
  (display "transitions:")(newline)
  (for-each (lambda (item)(display item)(newline))(aut-trans a)))

; display an aut-f

(define (display-aut-f mach)
  (display "AUT-F")(newline)
  (display-aut (aut-f-aut mach))
  (display "final set of states: ")
  (display (aut-f-final-states mach))(newline))


; Create a random (nondeterministic) automaton, given
; alphabet, number of states, number of transitions.
; The states are the positive integers from 1 to the number of states.
; The initial states are just '(1).
; The transitions are chosen randomly as ((state1 symbol) state2)
; and are not necessarily distinct.

; Example:
;> (random-aut '(a b) 4 10)
;(aut
; '(a b)
; '(1 2 3 4)
; '(1)
; (list
;  (entry '(1 a) 4)
;  (entry '(4 b) 2)
;  (entry '(4 b) 3)
;  (entry '(1 a) 4)
;  (entry '(2 b) 4)
;  (entry '(3 b) 2)
;  (entry '(3 a) 2)
;  (entry '(2 b) 1)
;  (entry '(2 b) 3)
;  (entry '(1 b) 1)))

(define (random-aut alphabet number-of-states number-of-transitions)
  (let* ((states (range 1 (+ 1 number-of-states)))
         (transitions
          (map (lambda (transition)
                 (entry (list (pick states) (pick alphabet)) 
		 	(pick states)))
               (range 1 (+ 1 number-of-transitions)))))
    (aut alphabet states '(1) transitions)))

; Create a random reverse-deterministic automaton
; with the specified number of missing transitions.
; This is used as part of generating random SUBAs.

(define (random-reverse-det-aut alphabet number-of-states number-missing)
;  (displayln "random-reverse-det-aut")
  (let* ((states (range 1 (+ 1 number-of-states)))
         (transitions
          (apply append
                 (map (lambda (state)
                        (map (lambda (symbol)
                               (let ((second-state (pick states)))
                                 (entry (list second-state symbol) state)))
                             alphabet))
                      states)))
         (removed-transitions (pick-k number-missing transitions))
         (remaining-transitions (set-diff transitions removed-transitions)))
    
    (aut alphabet states '(1) remaining-transitions)))


; Create a random aut-f given alphabet, number of states and number of transitions.
; First ceate a random automaton (with random-aut),
; and then make the final states a randomly chosen subset of the states.
; Example:
;> (random-aut-f '(a b) 4 10)
;(aut-f
; (aut
;  '(a b)
;  '(1 2 3 4)
;  '(1)
;  (list
;   (entry '(3 b) 2)
;   (entry '(3 b) 3)
;   (entry '(2 a) 3)
;   (entry '(3 a) 3)
;   (entry '(3 a) 3)
;   (entry '(4 a) 1)
;   (entry '(2 b) 4)
;   (entry '(1 a) 3)
;   (entry '(4 b) 4)
;   (entry '(3 a) 4)))
; '(1 2))

(define (random-aut-f alphabet number-of-states number-of-transitions)
  (let* ((aut (random-aut alphabet number-of-states number-of-transitions))
         (final-states (random-subset (aut-states aut))))
    (aut-f aut final-states)))

; Create a random aut-f given alphabet, number of states and number of transitions.
; The number of final states is specified -- a random subset of the states with
; the given cardinality is chosen for the set of final states.

(define (random-aut-f-k
         alphabet number-of-states number-of-transitions number-of-final-states)
  (let* ((aut1 (random-aut alphabet number-of-states number-of-transitions)))
    (aut-f
     aut1
     (pick-k number-of-final-states (aut-states aut1)))))

;; Example aut-f (which is deterministic)
;; As an NBA, it accepts strings of a's and b's with infinitely many a's

(define nba-inf-as
  (aut-f
   (aut
    '(a b)
    '(0 1)
    '(0)
    (list
     (entry '(0 a) 0)
     (entry '(0 b) 1)
     (entry '(1 a) 0)
     (entry '(1 b) 1)))
   '(0)))

; Another example aut-f, which is also deterministic.
; As an NBA, it accepts strings of a's and b's with infinitely many a's and infinitely many b's

(define nba-inf-as-inf-bs
  (aut-f
   (aut
    '(a b)
    '(0 1 2)
    '(0)
    (list
     (entry '(0 a) 1)
     (entry '(0 b) 2)
     (entry '(1 a) 1)
     (entry '(1 b) 0)
     (entry '(2 a) 0)
     (entry '(2 b) 2)))
   '(0)))
    
; Find automaton next states
; from a given state with a given input symbol.
;> (aut-next-states 'b 3 aut-ex1)
;'(4 3)

(define (aut-next-states symbol state aut1)
  (let ((trans (aut-trans aut1)))
    (remove-duplicates
     (lookup-all (list state symbol) trans))))
        
; Find automaton states reached from a given state
; with a given input string of symbols.
; Example:
;> (aut-string-states '(b b a) 3 aut-ex1)
;'(3 2)

(define (aut-string-states str state aut1)
  (if (null? str)
      (list state)
      (let* [(next-states (aut-next-states (first str) state aut1))]
        (remove-duplicates
         (apply append
                (map (lambda(next-state)
                       (aut-string-states (rest str) next-state aut1))
                     next-states))))))

; Create a hash table for the transition relation of an AUT
; in which the key is a source state and the value is the list of
; all the destination states it can reach on any symbol.

(define (make-trans-hash all-states trans)
;  (displayln "make-trans-hash")
  (let* [(ht (make-hash))]
    (for-each
     (lambda (state)
       (hash-set! ht state '()))
     all-states)
    (for-each
     (lambda (tuple)
       (let [(source-state (first (entry-key tuple)))
             (dest-state (entry-value tuple))]
       (hash-set! ht source-state
                  (cons dest-state (hash-ref ht source-state)))))
     trans)
    ht))

; Given an automaton and a set of states, determine all states
; reachable from the given states.  Uses hash tables for speed.

(define (aut-reachable aut1 states)
;  (displayln "aut-reachable")
  (let* [(trans-hash (make-trans-hash (aut-states aut1) (aut-trans aut1)))
         (seen-hash (make-hash))]
    (for-each
     (lambda (state)
       (hash-set! seen-hash state #t))
     states)
    (aut-reachable-loop trans-hash seen-hash states)))

(define (aut-reachable-loop trans-hash seen-hash to-do)
  (if (empty? to-do)
      (hash-keys seen-hash)
      (let* [(current-state (first to-do))
             (state-successors (hash-ref trans-hash current-state '()))
             (new-states
              (filter (lambda (state)
                        (if (hash-ref seen-hash state #f)
                            #f
                            (begin
                              (hash-set! seen-hash state #t)
                              #t)))
                      state-successors))]
        (aut-reachable-loop
         trans-hash seen-hash (append new-states (rest to-do))))))

; Given an automaton and a set of states, determine all states
; that reach any of the given states.  Reverses the transition relation
; and calls aut-reachable-loop for speed.

(define (aut-reaches aut1 states)
;  (displayln "aut-reaches")
;  (check-aut aut1)
  (let* [(r-trans (reverse-transition (aut-trans aut1)))
         (trans-hash (make-trans-hash (aut-states aut1) r-trans))
         (seen-hash (make-hash))]
    (for-each
     (lambda (state)
       (hash-set! seen-hash state #t))
     states)
    (aut-reachable-loop trans-hash seen-hash states)))

; error check that all the states in a transition occur in states
(define (check-aut aut1)
  (let* [(states (aut-states aut1))
         (trans (aut-trans aut1))]
    (for-each
     (lambda (tuple)
       (if (not (member? (entry-value tuple) states))
           (begin (displayln "check-aut")
                  (displayln (entry-value tuple))
                  (error "check-aut entry-value not in states"))
           #t))
     trans)))

; Reverse a transition relation.
; ((q1 a) q2) becomes ((q2 a) q1)

(define (reverse-transition trans)
  (map (lambda (item)
         (let ((key (entry-key item))
               (value (entry-value item)))
           (entry
            (list value (second key))
            (first key))))
       trans))

; Test a transition relation for being deterministic.
(define (deterministic-trans? trans)
  (let* ((dedup-entries (remove-duplicates trans))
         (keys (map entry-key dedup-entries)))
    (equal? keys (remove-duplicates keys))))
  
; Test a transition relation for being complete.
; (Not used?)
(define (complete-trans? alphabet states trans)
  (let ((result (map (lambda (item) (lookup item trans))
                     (c-product states alphabet))))
    (if (member #f result)
        #f
        #t)))

; Convert an automaton over an arbitrary alphabet to use symbols a, b, c, ...,
; if possible.  Error returned if too many alphabet symbols (more than 26.)

(define (aut-symbols->letters aut1)
  (let* [(alphabet (aut-alphabet aut1))
         (len-alphabet (length alphabet))]
    (if (> len-alphabet 26)
        (begin
          (display "Too many symbols to convert to letters: ")
          (displayln (length alphabet))
          (error "symbols->letters"))
        (let* [(letters '(a b c d e f g h i j k l m n o p q r s t u v w x y z))
               (used-letters (take letters len-alphabet))
               (table (map (lambda (sym letter)
                             (entry sym letter))
                           alphabet
                           used-letters))]
          (aut
           used-letters
           (aut-states aut1)
           (aut-init-states aut1)
           (map (lambda (tuple)
                  (entry (list (first (entry-key tuple))
                               (lookup (second (entry-key tuple)) table))
                         (entry-value tuple)))
                (aut-trans aut1)))))))
