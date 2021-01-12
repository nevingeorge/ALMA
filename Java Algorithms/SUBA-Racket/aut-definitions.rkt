#lang racket
; Data structures for automata (aut) and automata with final states (aut-f)
; Dana Angluin
; 9/27/2020

(provide aut aut? aut-alphabet aut-states aut-init-states aut-trans
         aut-f aut-f? aut-f-aut aut-f-final-states
         aut-f-alphabet aut-f-states aut-f-init-states aut-f-trans)

(require "general-utilities.rkt")
(require "sequence-utilities.rkt")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; An automaton is represented as a structure of
; alphabet, states, initial states, and transition relation.
; The transition relation may be nondeterministic.

; alphabet: a list of symbols
; states: a list of distinct values
; initial states: a list of states
; transition relation: a list of entries, 
;     each entry a list ((state1 symbol) state2)

(struct aut (alphabet states init-states trans) #:transparent)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Automaton with final states

; This consists of an automaton and a set of final states.
; It can be interpreted as either
; a nondeterministic finite acceptor (NFA), or
; a nondeterministic Buchi acceptor (NBA).

(struct aut-f (aut final-states) #:transparent)

; Definitions of selectors for components of an aut-f
; accessed by reaching into the automaton component.

(define aut-f-alphabet (lambda (mach) (aut-alphabet (aut-f-aut mach))))
(define aut-f-states (lambda (mach) (aut-states (aut-f-aut mach))))
(define aut-f-init-states (lambda (mach) (aut-init-states (aut-f-aut mach))))
(define aut-f-trans (lambda (mach) (aut-trans (aut-f-aut mach))))

; Example automaton (nondeterministic):

(define aut-ex1
  (aut '(a b)
       '(1 2 3 4)
       '(1 3)
        (list
         (entry '(1 a) 2)
         (entry '(1 a) 3)
         (entry '(1 b) 1)
         (entry '(2 a) 1)
         (entry '(2 b) 4)
         (entry '(3 b) 4)
         (entry '(3 b) 3)
         (entry '(4 a) 3)
         (entry '(4 a) 2))))


